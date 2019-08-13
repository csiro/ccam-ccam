! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module cloudmod

! prognostic cloud fraction scheme based on Tiedtke from GFDL-CM3.

! This module works with the LDR cloud microphysics scheme with
! prognostic cloud liquid, frozen/snow and rain condensate.
    
implicit none
    
private
public progcloud
public convectivecloudfrac, convectivecloudarea

real, save :: u00ramp = 0.01

contains
    
subroutine progcloud(qc,qtot,press,rho,fice,qs,t,rhcrit, &
                     dpsldt,nettend,stratcloud)

use const_phys           ! Physical constants
use parm_m               ! Model configuration

implicit none

include 'kuocom.h'       ! Convection parameters

integer k, kl
real, dimension(:,:), intent(inout) :: qc ! condensate = qf + ql
real, dimension(:,:), intent(in) :: qtot, rho, fice, qs, t, rhcrit, press
real, dimension(:,:), intent(in) :: dpsldt
real, dimension(:,:), intent(inout) :: nettend
real, dimension(:,:), intent(inout) :: stratcloud
real, dimension(size(qc,1)) :: aa, bb, cc, at, a_dt, b_dt, cf1, cfeq, cfbar
real, dimension(size(qc,1)) :: qv, omega, hlrvap, dqsdT, gamma, xf, dqs
real erosion_scale

kl = size(qc,2)

stratcloud = max( min( stratcloud, 1. ), 0. )

! background erosion scale in 1/secs
erosion_scale = 1.E-6

!if ( ncloud>=5 ) then
!  ! convert convective mass flux from half levels to full levels
!  do k = 1,kl-1
!    cmflx(:,k) = rathb(k)*fluxtot(:,k)+ratha(k)*fluxtot(:,k+1)
!  end do
!  cmflx(:,kl) = rathb(kl)*fluxtot(:,kl)
!else ! ncloud==4
!  ! use convective area fraction in leoncld.f, instead of convective mass flux
!  cmflx = 0.
!end if

! calculate dqs = ((omega + grav*Mc)/(cp*rho)+nettend)*dqsdT*dt
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma
! MJT notes - GFDL AM adds (stratcloud+0.5*at*da)*gamma term

! Change in saturated volume fraction
! da = -0.5*(1.-cf)^2*dqs/(qs-qv)
! MJT notes - Tiedtke 93 does not use 0.5

! gamma = L/cp*dqsdT

! Follow GFDL AM approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! dqs*dqs*AA + dqs*BB + CC = 0
! AA = 0.25*gamma*(1-cf)^2/(qs-qv)
! BB = -(1+gamma*cf)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

do k = 1,kl
  qv = qtot(:,k) - qc(:,k)  
  ! calculate vertical velocity, dqs/dT and gamma 
  omega = press(:,k)*dpsldt(:,k)
  hlrvap = (hl+fice(:,k)*hlf)/rvap
  dqsdT = qs(:,k)*hlrvap/(t(:,k)**2)
  gamma = (hlcp+fice(:,k)*hlfcp)*dqsdT
  
  xf = max(min( (qv/qs(:,k) - rhcrit(:,k) - u00ramp ) / ( 2.*u00ramp ), 1. ), 0. ) ! MJT suggestion
  
  !cc = ((omega + grav*cmflx(:,k))/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT
  cc = (omega/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT ! neglect cmflx
  at = 1.-stratcloud(:,k)
  aa = 0.5*at*at/max( qs(:,k)-qv, 1.e-20 )
  bb = 1.+gamma*stratcloud(:,k)
  where ( cc<=0. .and. xf>0. )
    !dqs = ( bb - sqrt( bb*bb - 2.*gamma*xf*aa*cc ) ) / ( gamma*xf*aa ) ! GFDL style
    !dqs = min( dqs, cc/(1. + 0.5*bb) )                                 ! GFDL style
    dqs = 2.*cc/( bb + sqrt( bb*bb - 2.*gamma*xf*aa*cc ) ) ! alternative form of quadratic equation
                                                           ! note that aa and bb have been multipled by 2 and -1, respectively.
    ! Large scale cloud formation via condensation (A)
    a_dt = -xf*aa*dqs
  elsewhere
    ! da = 0, so dqs can be solved from a linear equation
    dqs = cc/bb
    ! Large scale cloud formation via condensation (A)
    a_dt = 0.
  end where

  ! Large scale cloud destruction via erosion (B)
  b_dt = stratcloud(:,k)*erosion_scale*dt*max(qs(:,k)-qv, 1.e-20)/max(qc(:,k), 1.e-20)

  ! Integrate
  !   dcf/dt = (1-cf)*A - cf*B
  ! to give (use cf' = A-cf*(A+B))
  !   cf(t=1) = cfeq + (cf(t=0) - cfeq)*exp(-(A+B)*dt)
  !   cfeq = A/(A+B)
  ! Average cloud fraction over the interval t=tau to t=tau+1
  !   cfbar = cfeq - (cf(t=1) - cf(t=0))/((A+B)*dt)
  ! cfeq is the equilibrum cloud fraction that is approached with
  ! a time scale of 1/(A+B)
  where ( a_dt>1.e-20 .or. b_dt>1.e-20 )
    cfeq  = a_dt/(a_dt+b_dt)
    cf1   = cfeq + (stratcloud(:,k) - cfeq)*exp(-a_dt-b_dt)
    cfbar = cfeq + (stratcloud(:,k) - cf1 )/(a_dt+b_dt)
  elsewhere
    cfeq  = stratcloud(:,k)
    cf1   = stratcloud(:,k)
    cfbar = stratcloud(:,k)
  end where

  ! Change in condensate
  ! dqc = -dqs*(stratcloud+0.5*da) = -dqs*cfbar
  ! MJT notes - missing erosion term -cfbar*erosion_scale*dt*(qs-qv)
  qc(:,k) = max(min( qc(:,k) - max(cfbar,1.e-20)*dqs, qtot(:,k)-qgmin ), 0. )

  ! Change in cloud fraction
  where ( qc(:,k)>1.e-20 )
    stratcloud(:,k) = max(min( cf1, 1.), 1.e-20 )
  elsewhere
    ! MJT notes - cloud fraction is maintained (da=0.) while condesate evaporates (dqc<0.) until
    ! the condesate dissipates
    stratcloud(:,k) = 0.
    qc(:,k) = 0.
  end where

  ! Reset tendency and mass flux for next time-step
  nettend(:,k) = 0.
  
end do

return
end subroutine progcloud

subroutine convectivecloudfrac(clcon,kbsav,ktsav,condc,cldcon)

use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer k, kl
real, dimension(:,:), intent(out) :: clcon
real, dimension(:), intent(out), optional :: cldcon
real, dimension(size(clcon,1)) :: cldcon_temp
real, dimension(size(clcon,1)) :: n, cldcon_local
integer, dimension(:), intent(in) :: kbsav
integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc

! MJT notes - This is an old parameterisation from NCAR.  acon and
! bcon represent shallow and deep convection, respectively.  It can
! be argued that acon should be zero in CCAM to avoid discontinuous
! evolution of the model.  Furthermore, a much better fit can be
! obtained from the mass flux, rather than rainfall.  It also
! should be noted that acon and bcon are likely to depend on the
! spatial resolution.

kl = size(clcon,2)

cldcon_temp = 0. ! for cray compiler
call convectivecloudarea(cldcon_temp,ktsav,condc)
if ( present(cldcon) ) then
  cldcon = cldcon_temp
end if

! Impose cloud overlap assumption
if ( nmr>=1 ) then
  cldcon_local = cldcon_temp               ! maximum overlap
else
  n = 1./real(max(ktsav-kbsav,1))
  cldcon_local = 1. - (1.-cldcon_temp)**n  ! random overlap
end if

do k = 1,kl
  where( k<kbsav+1 .or. k>ktsav )
    clcon(:,k) = 0.
  elsewhere
    clcon(:,k) = cldcon_local
  end where
end do

return
end subroutine convectivecloudfrac

pure subroutine convectivecloudarea(cldcon,ktsav,condc)

use newmpar_m        ! Grid parameters
use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc
real, dimension(:), intent(out) :: cldcon

where ( ktsav<kl-1 )
  cldcon = min( acon+bcon*log(1.+condc*86400./dt), 0.8 ) !NCAR
elsewhere
  cldcon = 0.
end where

return
end subroutine convectivecloudarea

end module cloudmod
