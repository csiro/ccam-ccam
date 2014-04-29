module cloudmod

! prognostic cloud scheme based on Tiedtke from GFDL-CM3.
    
implicit none
    
private
public progcloud, cloudmod_init, combinecloudfrac
public stratcloud, nettend, cmflx

real, dimension(:,:), allocatable, save :: stratcloud  ! prognostic cloud fraction
real, dimension(:,:), allocatable, save :: nettend     ! change in temperature from radiation and vertical mixing
real, dimension(:,:), allocatable, save :: cmflx       ! mass flux from convection

contains

subroutine cloudmod_init(ifull,iextra,kl,ncloud)

implicit none

integer, intent(in) :: ifull, iextra, kl, ncloud

if (ncloud>=3) then
  allocate(stratcloud(ifull+iextra,kl))
  allocate(nettend(ifull,kl),cmflx(ifull,kl))
  stratcloud=0.
  nettend=0.
  cmflx=0.
end if

return
end subroutine cloudmod_init
    
subroutine progcloud(cloudfrac,qc,qv,t,ps,rho,fice)

use estab                ! Liquid saturation function
use sigs_m               ! Atmosphere sigma levels
use vvel_m               ! Additional vertical velocity

implicit none

include 'newmpar.h'      ! Grid parameters
include 'const_phys.h'   ! Physical constants
include 'parm.h'         ! Model configuration

real, dimension(ifull,kl), intent(out) :: cloudfrac
real, dimension(ifull,kl), intent(inout) :: qc ! condensate = qf + ql
real, dimension(ifull,kl), intent(in) :: qv, t, rho, fice
real, dimension(ifull), intent(in) :: ps
real, dimension(ifull,kl) :: u00p, erosion_scale
real, dimension(ifull,kl) :: da, dqs, cfbar
real, dimension(ifull,kl) :: cf1, cfeq, a_dt, b_dt
real, dimension(ifull,kl) :: qs, qtot, dqsdT, mflx, gamma
real, dimension(ifull,kl) :: aa, bb, cc, omega, hlrvap
real es, pk
integer iq, k

stratcloud(1:ifull,:)=max( min( stratcloud(1:ifull,:), 1. ), 0. )
qtot = qc + qv

! saturated water mixing ratio
do k=1,kl
  do iq=1,ifull
    es=establ(t(iq,k))
    pk=ps(iq)*sig(k)
    qs(iq,k)=0.622*es/(pk-es)
  end do
end do

! Critical relative humidity (neglected profile option)
u00p(:,:) = 0.8

! background erosion scale
erosion_scale(:,:) = 1.E-6

! calculate the condensation without supersaturation

do k=1,kl
  omega(:,k) = ps(1:ifull)*dpsldt(:,k)
end do
hlrvap = (hl+fice*hlf)/rvap
dqsdT  = qs*hlrvap/(t*t)
gamma  = (hl+fice*hlf)/rvap

! calculate dqs = (((omega + grav*Mc)/(cp*rho)+nettend)*dqsdT*dt)
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma

! Follow GFDL CM3 approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! dqs = (-BB + sqrt( BB*BB - 4*AA*CC ))/(2*AA)
! AA = 0.25*gamma*(1-cf)^2/(qs-qv)
! BB = -(1+gamma*cf)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

cc = ((omega + grav*cmflx)/(cp*rho)+nettend)*dt*dqsdT
where (cc<=0. .and. qv>u00p*qs)
  aa = 0.25*gamma*(1.-stratcloud(1:ifull,:))*(1.-stratcloud(1:ifull,:))/max(qs-qv,1.E-20)
  bb = -(1.+gamma*stratcloud(1:ifull,:))
  dqs = ( -bb - sqrt( bb*bb - 4.*aa*cc ) ) / (2.*aa)
  !dqs = min( tmp, cc/(1.+0.5*gamm*(1.+stratcloud(1:ifull,:))) )
elsewhere
  ! da = 0
  dqs = cc/(1.+gamma*stratcloud(1:ifull,:))
end where

! Change in saturated volume fraction
where( dqs<=0. .and. qv>u00p*qs )
  da = -0.5*(1.-stratcloud(1:ifull,:))*(1.-stratcloud(1:ifull,:))*dqs/max(qs-qv,1.E-20)
elsewhere
  da = 0.
end where

! Large scale cloud formation (A)
a_dt = da/max(1.-stratcloud(1:ifull,:),1.e-10)

! Large scale cloud destruction (B)
b_dt = stratcloud(1:ifull,:)*erosion_scale*dt*max(qs-qv,0.)/max(qc,1.e-8 )

! Integrate
!   dcf/dt = (1-cf)*A - cf*B
! to give
!   cf(t=1) = cfeq - (cfeq - cf(t=0))*exp(-(A+B)*dt)
!   cfeq = A/(A+B)
! Average cloud fraction over the interval t=tau to t=tau+1
!   cfbar = cfeq - (cf(t=1) - cf(t=0))/(dt*(A+B))
! cfeq is the equilibrum cloud fraction that is approached with
! a time scale of 1/(A+B)
where (a_dt>1.e-8.or.b_dt>1.e-8)
  cfeq = a_dt/(a_dt + b_dt)
  cf1 = min(max( cfeq - (cfeq - stratcloud(1:ifull,:))*exp(-(a_dt + b_dt)), 0.), 1.)
  cfbar = min(max( cfeq - (cf1 - stratcloud(1:ifull,:))/(a_dt + b_dt), 0.), 1.)
elsewhere
  cf1 = stratcloud(1:ifull,:)
  cfbar = stratcloud(1:ifull,:)
end where

! Change in condensate
! dqc/dt = -dqs*(stratcloud+0.5*da) = -dqs*cfbar
qc = qc - cfbar*dqs
qc = min( max( qc, 0. ), qtot - qgmin )

! Change in cloud fraction
stratcloud(1:ifull,:) = min( max( cf1, 0. ), 1. )
where(qc>0.)
  stratcloud(1:ifull,:) = 1.E-8
end where
cloudfrac = stratcloud(1:ifull,:)

! Reset tendency and mass flux for next time-step
nettend=0.
cmflx=0.

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

integer iq, k
real, dimension(ifull) :: cldcon
real, dimension(ifull,kl) :: clcon

! estimate convective cloud fraction from leoncld.f
where (ktsav<kl-1)
  cldcon=min(acon+bcon*log(1.+condc*86400./dt),0.8) !NCAR
elsewhere
  cldcon=0.
end where
if (nmr>=1) then
  do iq=1,ifull
    do k=1,kbsav(iq)
      clcon(iq,k)=0.
    end do
    do k=kbsav(iq)+1,ktsav(iq)
      clcon(iq,k)=cldcon(iq) ! maximum overlap
    end do
    do k=ktsav(iq)+1,kl
      clcon(iq,k)=0.
    end do
  end do  
else
  do iq=1,ifull
    do k=1,kbsav(iq)
      clcon(iq,k)=0.
    end do
    do k=kbsav(iq)+1,ktsav(iq)
      clcon(iq,k)=1.-(1.-cldcon(iq))**(1./real(ktsav(iq)-kbsav(iq)+2)) !Random overlap
    end do
    do k=ktsav(iq)+1,kl
      clcon(iq,k)=0.
    end do
  end do  
end if

cfrac(:,:)=stratcloud(1:ifull,:)*(1.-clcon(:,:))+clcon(:,:)

return
end subroutine combinecloudfrac

end module cloudmod