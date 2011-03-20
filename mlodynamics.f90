! These subroutines handle dynamics for the Mixed-Layer-Ocean model

module mlodynamics

implicit none

private
public mlodiffusion,mlorouter,mlohadv,watbdy

real, dimension(:), allocatable, save :: watbdy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion
subroutine mlodiffusion

use cc_mpi
use indices_m
use map_m
use mlo
use soil_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer k,i
integer, parameter :: salfilt = 0 ! Additional salinity filter (0=off, 1=Katzfey)
real, parameter :: k_smag = 2. ! 0.4 in mom3, 2. in Griffies (2000)
real hdif
real, dimension(ifull+iextra) :: uc,vc,wc,ee,dz
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: u,v,base
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc
logical, dimension(ifull+iextra) :: wtr

hdif=dt*(k_smag/pi)**2
ee=0.
where(land(1:ifull))
  ee(1:ifull)=1.
end where
call bounds(ee)
wtr=ee.lt.0.5

if (.not.any(wtr(1:ifull))) return

emi=1./em(1:ifull)**2
do k=1,wlev
  call mloexpdep(1,dz(1:ifull),k,0)
  call bounds(dz)
  dz=max(dz,1.E-3)

  ! neglect bathymetry following component for now
  ! u=ax*uc+ay*vc+az*wc
  ! dudx*dx=u(ie)-u(iw)=ax*(uc(ie)-uc(iw))+ay*(vc(ie)-vc(iw))+az*(wc(ie)-wc(iw))
  ! dudy*dy=u(in)-u(is)=ax*(uc(in)-uc(is))+ay*(vc(in)-vc(is))+az*(wc(in)-wc(is))
  u=0.
  call mloexport(2,u,k,0)
  uc(1:ifull)=ax(1:ifull)*u
  vc(1:ifull)=ay(1:ifull)*u
  wc(1:ifull)=az(1:ifull)*u
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)
  dudx=0.
  dudy=0.
  where (wtr(1:ifull).and.wtr(ie).and.wtr(iw))
    dudx=(ax(1:ifull)*(uc(ie)-uc(iw))  &
         +ay(1:ifull)*(vc(ie)-vc(iw))  &
         +az(1:ifull)*(wc(ie)-wc(iw))) &
         *0.5*em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(ie))
    dudx=(ax(1:ifull)*(uc(ie)-uc(1:ifull))  &
         +ay(1:ifull)*(vc(ie)-vc(1:ifull))  &
         +az(1:ifull)*(wc(ie)-wc(1:ifull))) &
         *em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(iw))
    dudx=(ax(1:ifull)*(uc(1:ifull)-uc(iw))  &
         +ay(1:ifull)*(vc(1:ifull)-vc(iw))  &
         +az(1:ifull)*(wc(1:ifull)-wc(iw))) &
         *em(1:ifull)/ds  
  end where
  where (wtr(1:ifull).and.wtr(in).and.wtr(is))       
    dudy=(ax(1:ifull)*(uc(in)-uc(is))  &
         +ay(1:ifull)*(vc(in)-vc(is))  &
         +az(1:ifull)*(wc(in)-wc(is))) &
         *0.5*em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(in))
    dudy=(ax(1:ifull)*(uc(in)-uc(1:ifull))  &
         +ay(1:ifull)*(vc(in)-vc(1:ifull))  &
         +az(1:ifull)*(wc(in)-wc(1:ifull))) &
         *em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(is))
    dudy=(ax(1:ifull)*(uc(1:ifull)-uc(is))  &
         +ay(1:ifull)*(vc(1:ifull)-vc(is))  &
         +az(1:ifull)*(wc(1:ifull)-wc(is))) &
         *em(1:ifull)/ds
  end where

  ! v=bx*uc+by*vc+bz*wc
  ! dvdx*dx=v(ie)-v(iw)=bx*(uc(ie)-uc(iw))+by*(vc(ie)-vc(iw))+bz*(wc(ie)-wc(iw))
  ! dvdy*dy=v(in)-v(is)=bx*(uc(in)-uc(is))+by*(vc(in)-vc(is))+bz*(wc(in)-wc(is))
  v=0.
  call mloexport(3,v,k,0)
  uc(1:ifull)=bx(1:ifull)*v
  vc(1:ifull)=by(1:ifull)*v
  wc(1:ifull)=bz(1:ifull)*v
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)
  dvdx=0.
  dvdy=0.
  where (wtr(1:ifull).and.wtr(ie).and.wtr(iw))
    dvdx=(bx(1:ifull)*(uc(ie)-uc(iw))  &
         +by(1:ifull)*(vc(ie)-vc(iw))  &
         +bz(1:ifull)*(wc(ie)-wc(iw))) &
         *0.5*em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(ie))
    dvdx=(bx(1:ifull)*(uc(ie)-uc(1:ifull))  &
         +by(1:ifull)*(vc(ie)-vc(1:ifull))  &
         +bz(1:ifull)*(wc(ie)-wc(1:ifull))) &
         *em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(iw))
    dvdx=(bx(1:ifull)*(uc(1:ifull)-uc(iw))  &
         +by(1:ifull)*(vc(1:ifull)-vc(iw))  &
         +bz(1:ifull)*(wc(1:ifull)-wc(iw))) &
         *em(1:ifull)/ds
  end where
  where (wtr(1:ifull).and.wtr(in).and.wtr(is))
    dvdy=(bx(1:ifull)*(uc(in)-uc(is))  &
         +by(1:ifull)*(vc(in)-vc(is))  &
         +bz(1:ifull)*(wc(in)-wc(is))) &
         *0.5*em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(in))
    dvdy=(bx(1:ifull)*(uc(in)-uc(1:ifull))  &
         +by(1:ifull)*(vc(in)-vc(1:ifull))  &
         +bz(1:ifull)*(wc(in)-wc(1:ifull))) &
         *em(1:ifull)/ds
  elsewhere (wtr(1:ifull).and.wtr(is))
    dvdy=(bx(1:ifull)*(uc(1:ifull)-uc(is))  &
         +by(1:ifull)*(vc(1:ifull)-vc(is))  &
         +bz(1:ifull)*(wc(1:ifull)-wc(is))) &
         *em(1:ifull)/ds
  end where

  uc(1:ifull) = ax(1:ifull)*u + bx(1:ifull)*v
  vc(1:ifull) = ay(1:ifull)*u + by(1:ifull)*v
  wc(1:ifull) = az(1:ifull)*u + bz(1:ifull)*v
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)

  ! Smagorinsky
  cc=(dudx-dvdy)**2+(dudy+dvdx)**2
  t_kh(1:ifull)=sqrt(cc)*hdif*emi  ! this one with em in D terms
  call bounds(t_kh)
  xfact=0.
  yfact=0.
  where (wtr(1:ifull).and.wtr(ie))
    xfact(1:ifull) = (t_kh(ie)+t_kh)*.5
  end where
  where (wtr(1:ifull).and.wtr(in))
    yfact(1:ifull) = (t_kh(in)+t_kh)*.5
  end where
  call boundsuv(xfact,yfact)
  
  base= ( emi*dz(1:ifull) +                    &
          xfact(1:ifull)*dz(ie) +              & 
          xfact(iwu)*dz(iw) +                  &
          yfact(1:ifull)*dz(in) +              &
          yfact(isv)*dz(is) )

  ucc = ( uc(1:ifull)*emi*dz(1:ifull) +        &
          xfact(1:ifull)*uc(ie)*dz(ie) +       &
          xfact(iwu)*uc(iw)*dz(iw) +           &
          yfact(1:ifull)*uc(in)*dz(in) +       &
          yfact(isv)*uc(is)*dz(is) ) / base
  vcc = ( vc(1:ifull)*emi*dz(1:ifull) +        &
          xfact(1:ifull)*vc(ie)*dz(ie) +       &
          xfact(iwu)*vc(iw)*dz(iw) +           &
          yfact(1:ifull)*vc(in)*dz(in) +       &
          yfact(isv)*vc(is)*dz(is) ) / base
  wcc = ( wc(1:ifull)*emi*dz(1:ifull) +        &
          xfact(1:ifull)*wc(ie)*dz(ie) +       &
          xfact(iwu)*wc(iw)*dz(iw) +           &
          yfact(1:ifull)*wc(in)*dz(in) +       &
          yfact(isv)*wc(is)*dz(is) ) / base
  u = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u,k,0)
  call mloimport(3,v,k,0)
   
  do i=0,1
    ee=0.
    call mloexport(i,ee(1:ifull),k,0)
    call bounds(ee)
    ff = ( ee(1:ifull)*emi*dz(1:ifull) +        &
           xfact(1:ifull)*ee(ie)*dz(ie) +       &
           xfact(iwu)*ee(iw)*dz(iw) +           &
           yfact(1:ifull)*ee(in)*dz(in) +       &
           yfact(isv)*ee(is)*dz(is) ) / base
    call mloimport(i,ff,k,0)
  end do
  
  ! Jack Katzfey salinity filter
  if (salfilt.eq.1) then
    ee(1:ifull)=ff
    call bounds(ee)
    xfact=0.
    yfact=0.
    where (wtr(1:ifull).and.wtr(ie))
      xfact(1:ifull)=1.
    end where
    where (wtr(1:ifull).and.wtr(in))
      yfact(1:ifull)=1.
    end where
    call boundsuv(xfact,yfact)
    ff = ( ee(1:ifull)*emi*dz(1:ifull) +        &
           xfact(1:ifull)*ee(ie)*dz(ie) +       &
           xfact(iwu)*ee(iw)*dz(iw) +           &
           yfact(1:ifull)*ee(in)*dz(in) +       &
           yfact(isv)*ee(is)*dz(is) ) / base
    call mloimport(1,ff,k,0)
  end if
  
end do

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing
! This version assumes that the orography is sufficently resolved
! to determine the river flow
subroutine mlorouter

use arrays_m
use cable_ccam 
use cc_mpi
use indices_m
use map_m
use mlo
use nsibd_m
use soil_m
use soilsnow_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer i,iq
integer, dimension(ifull,4) :: xp
real, dimension(ifull) :: newwat
real, dimension(ifull,4) :: dp,slope,mslope,vel,flow
real :: xx,yy

if (.not.allocated(watbdy)) then
  allocate(watbdy(ifull+iextra))
  watbdy=0.
end if

call bounds(watbdy)

newwat=watbdy(1:ifull)
xp(:,1)=in
xp(:,2)=ie
xp(:,3)=is
xp(:,4)=iw
do i=1,4
  dp(:,i)=0.5*(ds/em(1:ifull)+ds/em(xp(:,i)))
  slope(:,i)=(zs(1:ifull)-zs(xp(:,i)))/dp(:,i)
end do

! outflow
mslope=max(slope,0.)
vel=0.35*sqrt(mslope/0.00005) ! from Hal's Mk3.5 scheme
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=watbdy(1:ifull)/(dp(:,i)/(-vel(:,i)*dt)-1.) ! (mm)
end do
newwat=newwat+sum(flow,2)
  
! inflow
mslope=max(-slope,0.)
vel=0.35*sqrt(mslope/0.00005) ! from Hal's Mk3.5 scheme
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=watbdy(xp(:,i))/(dp(:,i)/(vel(:,i)*dt)-1.) ! (mm)
  flow(:,i)=flow(:,i)*em(xp(:,i))*em(xp(:,i))/(em(1:ifull)*em(1:ifull)) ! correct for changing grid size
end do
newwat=newwat+sum(flow,2)
  
! basin 
do iq=1,ifull
  if (all(slope(iq,:).lt.-1.E-10).and.land(iq)) then
  
    ! runoff is inserted into soil moisture since CCAM has discrete land and sea points
    xx=watbdy(iq)
    if (nsib.eq.4.or.nsib.eq.6.or.nsib.eq.7) then
      call cableinflow(iq,xx)
    else
      yy=min(xx,(ssat(isoilm(iq))-wb(iq,ms))*1000.*zse(ms))
      wb(iq,ms)=wb(iq,ms)+yy/(1000.*zse(ms))
      xx=max(xx-yy,0.)
    end if
    newwat(iq)=newwat(iq)-watbdy(iq)+xx

    ! runoff is removed (assume some sub-grid scale lake is present)
    !newwat(iq)=newwat(iq)-watbdy(iq)
  end if
end do

watbdy(1:ifull)=max(newwat,0.)

return
end subroutine mlorouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements some basic hydrostatic dynamics
! (note that this subroutine is a work in progress)

subroutine mlohadv

use arrays_m
use cc_mpi
use indices_m
use map_m
use mlo
use soil_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer l,ii
real, dimension(ifull+iextra) :: ee
real, dimension(ifull) :: xp,xm
real, dimension(ifull) :: dpsdx,dpsdy
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
real, dimension(ifull+iextra,wlev) :: dep,rho
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s
real, dimension(ifull,wlev) :: nuh,nvh
real, dimension(ifull,wlev) :: dzdx,dzdy
real, dimension(ifull,wlev) :: drdx,drdy,drdz
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
real, dimension(:,:), allocatable, save :: oldu1,oldv1
logical, dimension(ifull+iextra) :: wtr
integer, parameter :: lmax=0 ! 1=predictor-only, 2=predictor-corrector

!Define land/sea mask
ee=0.
where(land(1:ifull))
  ee(1:ifull)=1.
end where
call bounds(ee)
wtr=ee.lt.0.5

! Initialise arrays
do ii=1,wlev
  call mloexpdep(0,dep(1:ifull,ii),ii,0)
  call bounds(dep(:,ii))
  call mloexport(0,w_t(:,ii),ii,0)
  call mloexport(1,w_s(:,ii),ii,0)
  call mloexport(2,w_u(:,ii),ii,0)
  call mloexport(3,w_v(:,ii),ii,0)
end do
call bounds(ps)

dep=max(dep,1.E-3)

! Calculate depth gradient
dzdx=0.
dzdy=0.
do ii=1,wlev
  where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
    dzdx(:,ii)=(dep(ie,ii)-dep(iw,ii))*0.5*em(1:ifull)/ds
  elsewhere(wtr(1:ifull).and.wtr(ie))
    dzdx(:,ii)=(dep(ie,ii)-dep(1:ifull,ii))*em(1:ifull)/ds
  elsewhere(wtr(1:ifull).and.wtr(iw))
    dzdx(:,ii)=(dep(1:ifull,ii)-dep(iw,ii))*em(1:ifull)/ds
  end where
  where(wtr(1:ifull).and.wtr(in).and.wtr(is))
    dzdy(:,ii)=(dep(in,ii)-dep(is,ii))*0.5*em(1:ifull)/ds
  elsewhere(wtr(1:ifull).and.wtr(in))
    dzdy(:,ii)=(dep(in,ii)-dep(1:ifull,ii))*em(1:ifull)/ds
  elsewhere(wtr(1:ifull).and.wtr(is))
    dzdy(:,ii)=(dep(1:ifull,ii)-dep(is,ii))*em(1:ifull)/ds
  end where
end do

! Calculate surface pressure gradients
! (free surface gradients are neglected for now)
dpsdx=0.
dpsdy=0.
where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
  dpsdx=(ps(ie)-ps(iw))*0.5*em(1:ifull)/ds
elsewhere(wtr(1:ifull).and.wtr(ie))
  dpsdx=(ps(ie)-ps(1:ifull))*em(1:ifull)/ds
elsewhere(wtr(1:ifull).and.wtr(iw))
  dpsdx=(ps(1:ifull)-ps(iw))*em(1:ifull)/ds
end where
where(wtr(1:ifull).and.wtr(in).and.wtr(is))
  dpsdy=(ps(in)-ps(is))*0.5*em(1:ifull)/ds
elsewhere(wtr(1:ifull).and.wtr(in))
  dpsdy=(ps(in)-ps(1:ifull))*em(1:ifull)/ds
elsewhere(wtr(1:ifull).and.wtr(is))
  dpsdy=(ps(1:ifull)-ps(is))*em(1:ifull)/ds
end where

! Calculate pressure and gradient for slow mode
call mloexpdensity(rho(1:ifull,:),w_t,w_s,dep(1:ifull,:),0)
do ii=1,wlev
  call bounds(rho(:,ii))
end do
drdz(:,1)=(rho(1:ifull,2)-rho(1:ifull,1))/max(dep(1:ifull,2)-dep(1:ifull,1),1.E-20)
do ii=2,wlev-1
  drdz(:,ii)=(rho(1:ifull,ii+1)-rho(1:ifull,ii-1))/max(dep(1:ifull,ii+1)-dep(1:ifull,ii-1),1.E-20)
end do
drdz(:,wlev)=(rho(1:ifull,wlev)-rho(1:ifull,wlev-1))/max(dep(1:ifull,wlev)-dep(1:ifull,wlev-1),1.E-20)

drdx=0.
drdy=0.
do ii=1,wlev
  where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
    drdx(:,ii)=(rho(ie,ii)-rho(iw,ii))*0.5*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(ie))
    drdx(:,ii)=(rho(ie,ii)-rho(1:ifull,ii))*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(iw))
    drdx(:,ii)=(rho(1:ifull,ii)-rho(iw,ii))*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  end where
  where(wtr(1:ifull).and.wtr(in).and.wtr(is))
    drdy(:,ii)=(rho(in,ii)-rho(is,ii))*0.5*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(in))
    drdy(:,ii)=(rho(in,ii)-rho(1:ifull,ii))*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(is))
    drdy(:,ii)=(rho(1:ifull,ii)-rho(is,ii))*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  end where
end do

xp=1.
xm=0.
where (wtr(1:ifull))
  xp=1.+(dt*f(1:ifull))**2
  xm=dt*f(1:ifull)
end where

! Split U and V coriolis and pressure gradient terms
! (Here we differentiate for horizontal X and Y.  The correction for levels is built into rho, but the varying bottom
!  depth is then removed)
do ii=1,wlev
  nu(1:ifull,ii)=(w_u(:,ii)+xm*w_v(:,ii)-dt*grav*dep(1:ifull,ii)*(drdx(:,ii)+xm*drdy(:,ii))/rho(1:ifull,ii) &
                                              -dt*(dpsdx+xm*dpsdy)/rho(1:ifull,ii))/xp
  nv(1:ifull,ii)=(w_v(:,ii)-xm*w_u(:,ii)-dt*grav*dep(1:ifull,ii)*(drdy(:,ii)-xm*drdx(:,ii))/rho(1:ifull,ii) &
                                              -dt*(dpsdy-xm*dpsdx)/rho(1:ifull,ii))/xp
end do
w_u=nu(1:ifull,:)
w_v=nv(1:ifull,:)
nt(1:ifull,:)=w_t
ns(1:ifull,:)=w_s

if (.not.allocated(oldu1)) then
  allocate(oldu1(ifull,wlev))
  allocate(oldv1(ifull,wlev))
  oldu1=w_u
  oldv1=w_v
end if

! prep for conservation fix
! ...

do l=1,lmax ! predictor-corrector loop

  ! update bounds of prognostic variables
  do ii=1,wlev
    call bounds(nu(:,ii))
    call bounds(nv(:,ii))
    call bounds(nt(:,ii))
    call bounds(ns(:,ii))
  end do

  if (l.eq.1) then
    nuh(1:ifull,:)=1.5*w_u-0.5*oldu1 ! U at t+1/2
    nvh(1:ifull,:)=1.5*w_v-0.5*oldv1 ! V at t+1/2
  else
    nuh(1:ifull,:)=0.5*nu(1:ifull,:)+0.5*w_u ! U at t+1/2
    nvh(1:ifull,:)=0.5*nv(1:ifull,:)+0.5*w_v ! V at t+1/2
  end if

  ! Calculate depature points
  call mlodeps(nuh,nvh,x3d,y3d,z3d)

!  Vertical advection
!  call mlovadv

!  Convert (u,v) to cartesian coordinates (U,V,W)
!  do ii=1,wlev
!    cou(:,ii)=ax*ou(:,ii)+bx*ov(:,ii)
!    cov(:,ii)=ay*ou(:,ii)+by*ov(:,ii)
!    cow(:,ii)=az*ou(:,ii)+bz*ov(:,ii)
!  end do

!  Horizontal advection for U,V,W
!  call mloints(cou,nface,xg,yg)
!  call mloints(cov,nface,xg,yg)
!  call mloints(cow,nface,xg,yg)

!  Rotate vector to arrival point
!  call mlorot

!  Convert (U,V,W) back to conformal cubic coordinates
!  do ii=1,wlev
!    nu(:,ii)=ax*cou(:,ii)+ay*cov(:,ii)+az*cow(:,ii)
!    nv(:,ii)=bx*cou(:,ii)+by*cov(:,ii)+bz*cow(:,ii)
!  end do

!  Horizontal advector for T,S
!  call mloints(nt,nface,xg,yg)
!  call mloints(ns,nface,xg,yg)

end do

! fix conservation
! ...

oldu1=w_u
oldv1=w_v

do ii=1,wlev
  call mloimport(0,nt(1:ifull,ii),ii,0)
  call mloimport(1,ns(1:ifull,ii),ii,0)
  call mloimport(2,nu(1:ifull,ii),ii,0)
  call mloimport(3,nv(1:ifull,ii),ii,0)
end do

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(ubar,vbar,x3d,y3d,z3d)

use cc_mpi
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer ii
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real*8, dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do ii=1,wlev
  uc(:,ii)=(ax*ubar(:,ii)+bx*vbar(:,ii))*dt/rearth ! unit sphere 
  vc(:,ii)=(ay*ubar(:,ii)+by*vbar(:,ii))*dt/rearth ! unit sphere 
  wc(:,ii)=(az*ubar(:,ii)+bz*vbar(:,ii))*dt/rearth ! unit sphere 
  x3d(:,ii)=x-uc(:,ii) ! 1st guess
  y3d(:,ii)=y-vc(:,ii)
  z3d(:,ii)=z-wc(:,ii)
end do

! convert to grid point numbering
!do ii=1,wlev
!  call mlotoij5(ii,x3d(:,ii),y3d(:,ii),z3d(:,ii))
!end do
! Share off processor departure points.
!call deptsync(nface,xg,yg)

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices

subroutine mlotoij5(ii,x3d,y3d,z3d)

implicit none

include 'newmpar.h'

integer, intent(in) :: ii
real*8, dimension(ifull), intent(inout) :: x3d,y3d,z3d

return
end subroutine mlotoij5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tide

!subroutine tide(eta,lon,lat,mtimer)
!
!implicit none
!
!include 'newmpar.h'
!
!integer i
!real, dimension(ifull), intent(out) :: eta
!real, dimension(ifull), intent(in) :: lon,lat,mtimer
!real, dimension(ifull) :: slon,slat,stime
!real, dimension(4) :: aa,ab,ba,bb,wa,wb
!real, parameter :: pi = 3.1415927
!
!aa(1)=0.141565
!aa(2)=0.100661
!aa(3)=0.046848
!aa(4)=0.019273
!ab(1)=0.242334
!ab(2)=0.112743
!ab(3)=0.046397
!ab(4)=0.030684
!ba(1)=0.736
!ba(2)=0.695
!ba(3)=0.706
!ba(4)=0.695
!bb(1)=0.693
!bb(2)=0.693
!bb(3)=0.693
!bb(4)=0.693
!wa(1)=0.7292117
!wa(2)=0.6750774
!wa(3)=0.7252295
!wa(4)=0.6495854
!wb(1)=1.405189
!wb(2)=1.454441
!wb(3)=1.378797
!wb(4)=1.458423
!
!slon=pi*lon/180.
!slat=pi*lat/180.
!stime=mtimer/1440.
!
!! Could pre-compute trig terms
!eta=0.
!do i=1,4
!  eta=eta+ba(i)*aa(i)*cos(slat)**2*(cos(wa(i)*stime)*cos(2.*slon)-sin(wa(i)*stime)*sin(2.*slon))
!  eta=eta+bb(i)*ab(i)*sin(2.*slat)*(cos(wb(i)*stime)*cos(2.*slon)-sin(wb(i)*stime)*sin(2.*slon))
!end do
!
!return
!end subroutine tide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine fills ocean data over land points
!subroutine mlofill(x)
!
!use cc_mpi
!use indices_m
!use soil_m
!
!implicit none
!
!include 'newmpar.h'
!include 'mpif.h'
!
!integer globalc,localc,ierr,iq,lnum
!real miss,lsum
!real, dimension(ifull), intent(inout) :: x
!real, dimension(ifull+iextra) :: xx,yy
!logical, dimension(ifull+iextra) :: smap
!
!miss=999999.
!
!xx=miss
!where (.not.land)
!  xx(1:ifull)=x
!end where
!globalc=1
!
!! technically we only need 1 pass of this fill to ensure all water points have a non-trival neighbour
!do while (globalc.gt.0)
!  call bounds(xx)
!  smap=abs(xx-miss).lt.0.1
!  yy=xx
!  do iq=1,ifull
!    if (smap(iq)) then
!      lsum=0.
!      lnum=0
!      if (.not.smap(in(iq))) then
!        lsum=lsum+xx(in(iq))
!        lnum=lnum+1
!      end if
!      if (.not.smap(ie(iq))) then
!        lsum=lsum+xx(ie(iq))
!        lnum=lnum+1
!      end if
!      if (.not.smap(is(iq))) then
!        lsum=lsum+xx(is(iq))
!        lnum=lnum+1
!      end if
!      if (.not.smap(iw(iq))) then
!        lsum=lsum+xx(iw(iq))
!        lnum=lnum+1
!      end if
!      if (lnum.gt.0) then
!        yy(iq)=lsum/real(lnum)
!      end if
!    end if
!  end do
!  xx=yy
!  smap=abs(xx-miss).lt.0.1
!  localc=count(smap(1:ifull))
!  call MPI_AllReduce(localc,globalc,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!end do
!
!return
!end subroutine mlofill

end module mlodynamics