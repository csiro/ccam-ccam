module cloudmod

! prognostic cloud scheme based on Tiedtke from GFDL-CM3.
    
implicit none
    
private
public progcloud, cloudmod_init, combinecloudfrac
public stratcloud, nettend

real, dimension(:,:), allocatable, save :: stratcloud  ! prognostic cloud fraction
real, dimension(:,:), allocatable, save :: nettend     ! change in temperature from radiation and vertical mixing

contains

subroutine cloudmod_init(ifull,iextra,kl,ncloud)

implicit none

integer, intent(in) :: ifull, iextra, kl, ncloud

if (ncloud>=3) then
  allocate(stratcloud(ifull+iextra,kl))
  allocate(nettend(ifull,kl))
  stratcloud=0.
  nettend=0.
end if

return
end subroutine cloudmod_init
    
subroutine progcloud(cloudfrac,qc,qtot,ps,rho,fice,qs,t)

use kuocomb_m            ! JLM convection
use sigs_m               ! Atmosphere sigma levels
use vvel_m               ! Additional vertical velocity

implicit none

include 'newmpar.h'      ! Grid parameters
include 'const_phys.h'   ! Physical constants
include 'kuocom.h'       ! Convection parameters
include 'parm.h'         ! Model configuration

real, dimension(ifull,kl), intent(out) :: cloudfrac
real, dimension(ifull,kl), intent(inout) :: qc ! condensate = qf + ql
real, dimension(ifull,kl), intent(in) :: qtot, rho, fice, qs, t
real, dimension(ifull), intent(in) :: ps
real, dimension(ifull,kl) :: u00p, erosion_scale
real, dimension(ifull,kl) :: da, dqs, cfbar
real, dimension(ifull,kl) :: cf1, cfeq, a_dt, b_dt
real, dimension(ifull,kl) :: dqsdT, mflx, gamma
real, dimension(ifull,kl) :: aa, bb, cc, omega
real, dimension(ifull,kl) :: cmflx
integer k

stratcloud(1:ifull,:)=max( min( stratcloud(1:ifull,:), 1. ), 0. )

! Critical relative humidity (neglected profile option)
u00p(:,:) = 0.8

! background erosion scale
erosion_scale(:,:) = 1.E-6

! calculaate vertical velocity, dqs/dT and gamma
do k=1,kl
  omega(:,k) = ps(1:ifull)*dpsldt(:,k)
end do
gamma = (hl+fice*hlf)/rvap
dqsdT = qs*gamma/(t*t)
if ( ncloud>=4 ) then
  do k=1,kl-1
    cmflx(:,k)=rathb(k)*fluxtot(:,k)+ratha(k)*fluxtot(:,k+1)
  end do
  cmflx(:,kl)=rathb(kl)*fluxtot(:,kl)
else ! ncloud==3
  cmflx=0.
end if

! calculate dqs = (((omega + grav*Mc)/(cp*rho)+nettend)*dqsdT*dt)
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma

! Follow GFDL CM3 approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! MJT notes - if we use tliq, then nettend and omega should also be modified
! to account for changes in ql and qf

! dqs = (-BB + sqrt( BB*BB - 4*AA*CC ))/(2*AA)
! AA = 0.25*gamma*(1-cf)^2/(qs-qtot)
! BB = -(1+gamma*cf)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

cc = ((omega + grav*cmflx)/(cp*rho)+nettend)*dt*dqsdT
aa = 0.5*gamma*(1.-stratcloud(1:ifull,:))*(1.-stratcloud(1:ifull,:))/max(qs-qtot,1.E-20)
bb = 1.+gamma*stratcloud(1:ifull,:)
where ( cc<=0. .and. qtot>u00p*qs )
  dqs = 2.*cc/( bb + sqrt( bb*bb - 2.*aa*cc ) ) ! alternative form of quadratic equation
                                                ! note that aa has been multipled by 2.
elsewhere
  ! da = 0, so dqs can be solved from a linear equation
  dqs = cc/bb
end where

! Change in saturated volume fraction
! da = - 0.5*(1.-cf)^2*dqs/(qs-qtot)
where( dqs<=0. .and. qtot>u00p*qs )
  da = -aa*dqs/gamma
elsewhere
  da = 0.
end where

! Large scale cloud formation (A)
a_dt = da/max( 1.-stratcloud(1:ifull,:), 1.e-20 )

! Large scale cloud destruction (B)
b_dt = stratcloud(1:ifull,:)*erosion_scale*dt*max(qs-qtot,0.)/max(qc,1.e-8 )

! Integrate
!   dcf/dt = (1-cf)*A - cf*B
! to give (use cf' = A-cf*(A+B))
!   cf(t=1) = cfeq + (cf(t=0) - cfeq)*exp(-(A+B)*dt)
!   cfeq = A/(A+B)
! Average cloud fraction over the interval t=tau to t=tau+1
!   cfbar = cfeq - (cf(t=1) - cf(t=0))/((A+B)*dt)
! cfeq is the equilibrum cloud fraction that is approached with
! a time scale of 1/(A+B)
where ( a_dt>1.E-20 .or. b_dt>1.E-20 )
  cfeq  = a_dt/(a_dt+b_dt)
  cf1   = min(max( cfeq + (stratcloud(1:ifull,:) - cfeq)*exp(-a_dt-b_dt), 0.), 1.)
  cfbar = min(max( cfeq + (stratcloud(1:ifull,:) - cf1 )/(a_dt+b_dt),     0.), 1.)
elsewhere
  cfeq  = stratcloud(1:ifull,:)
  cf1   = stratcloud(1:ifull,:)
  cfbar = stratcloud(1:ifull,:)
end where

! Change in condensate
! dqc = -dqs*(stratcloud+0.5*da) = -dqs*cfbar
qc = qc - cfbar*dqs
qc = min( max( qc, 0. ), qtot - qgmin )

! Change in cloud fraction
where( qc>0. )
  stratcloud(1:ifull,:) = max( min( cf1, 1. ), 1.E-8 )
elsewhere
  stratcloud(1:ifull,:) = 0.
end where
cloudfrac = stratcloud(1:ifull,:)

! Reset tendency and mass flux for next time-step
nettend=0.

return
end subroutine progcloud


! This subroutine combines large scale and subgrid scale cloud fractions

subroutine combinecloudfrac

use cfrac_m          ! Cloud fraction
use kuocomb_m        ! JLM convection
use morepbl_m        ! Additional boundary layer diagnostics

implicit none

include 'newmpar.h'  ! Grid parameters
include 'kuocom.h'   ! Convection parameters
include 'parm.h'     ! Model configuration

integer k
real, dimension(ifull) :: cldcon
real, dimension(ifull,kl) :: clcon
real, dimension(ifull) :: n, crnd

if ( ncloud>=4 ) then
  cfrac(:,:)=stratcloud(1:ifull,:)
else

  ! estimate convective cloud fraction from leoncld.f

  ! MJT notes - This is an old parameterisation from NCAR.  acon and
  ! bcon represent shallow and deep convection, respectively.  It can
  ! be argued that acon should be zero in CCAM to avoid discontinuous
  ! evolution of the model.  Furthermore, a much better fit can be
  ! obtained from the mass flux, rather than rainfall.  It also
  ! should be noted that acon and bcon are likely to depend on the
  ! spatial resolution.
  where ( ktsav<kl-1 )
    cldcon=min(acon+bcon*log(1.+condc*86400./dt),0.8) !NCAR
  elsewhere
    cldcon=0.
  end where

  ! Impose cloud overlap assumption
  if ( nmr>=1 ) then
    do k=1,kl
      where( k<kbsav(:) .or. k>ktsav(:) )
        clcon(:,k)=0.
      elsewhere
        clcon(:,k)=cldcon ! maximum overlap
      end where
    end do
  else
    n=1./real(ktsav-kbsav+1)
    crnd=1.-(1.-cldcon)**n
    do k=1,kl
      where( k<kbsav .or. k>ktsav )
        clcon(:,k)=0.
      elsewhere
        clcon(:,k)=cldcon  ! random overlap
      end where
    end do
  end if

  cfrac(:,:)=stratcloud(1:ifull,:)*(1.-clcon(:,:))+clcon(:,:)
  
end if

return
end subroutine combinecloudfrac

end module cloudmod