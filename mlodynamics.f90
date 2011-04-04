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
  
  base= ( emi +                         &
          xfact(1:ifull) +              & 
          xfact(iwu) +                  &
          yfact(1:ifull) +              &
          yfact(isv) )

  ucc = ( uc(1:ifull)*emi +             &
          xfact(1:ifull)*uc(ie) +       &
          xfact(iwu)*uc(iw) +           &
          yfact(1:ifull)*uc(in) +       &
          yfact(isv)*uc(is) ) / base
  vcc = ( vc(1:ifull)*emi +             &
          xfact(1:ifull)*vc(ie) +       &
          xfact(iwu)*vc(iw) +           &
          yfact(1:ifull)*vc(in) +       &
          yfact(isv)*vc(is) ) / base
  wcc = ( wc(1:ifull)*emi +             &
          xfact(1:ifull)*wc(ie) +       &
          xfact(iwu)*wc(iw) +           &
          yfact(1:ifull)*wc(in) +       &
          yfact(isv)*wc(is) ) / base
  u = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u,k,0)
  call mloimport(3,v,k,0)
   
  do i=0,1
    ee=0.
    call mloexport(i,ee(1:ifull),k,0)
    call bounds(ee)
    ff = ( ee(1:ifull)*emi +             &
           xfact(1:ifull)*ee(ie) +       &
           xfact(iwu)*ee(iw) +           &
           yfact(1:ifull)*ee(in) +       &
           yfact(isv)*ee(is) ) / base
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
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer l,ll,ii,intsch
integer, dimension(ifull,wlev) :: nface
real, dimension(ifull+iextra) :: ee,neta,suu,svv
real, dimension(ifull) :: xp,xm,dpsdx,dpsdy,w_e
real, dimension(ifull) :: div,dedx,dedy
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
real, dimension(ifull+iextra,wlev) :: cou,cov,cow
real, dimension(ifull+iextra,wlev) :: dep,rhobar
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,rho,dz,dzbar
real, dimension(ifull,wlev) :: nuh,nvh
real, dimension(ifull,wlev) :: dzdx,dzdy
real, dimension(ifull,wlev) :: drdx,drdy,drdz
real, dimension(ifull,wlev) :: xg,yg
real, dimension(:,:), allocatable, save :: oldu1,oldv1
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
integer, parameter :: lmax=1 ! 0=no-advection, 1=predictor-only, 2=predictor-corrector
integer, parameter :: llmax=60 ! iterations for surface height
real, parameter :: alpha = 0.02

! new z levels for including free surface eta (effectively sigma levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth

intsch=mod(ktau,2)

!Define land/sea mask
ee=0.
where(land(1:ifull))
  ee(1:ifull)=1.
end where
call bounds(ee)
wtr=ee.lt.0.5

cou=0.
cov=0.
cow=0.
rhobar=1030.
dep=0.
dz=0.
w_t=293.
w_s=0.
w_u=0.
w_v=0.
w_e=0.
nt=293.
ns=30.
nu=0.
nv=0.
neta=0.

! ADVECT WATER ------------------------------------------------------
do ii=1,wlev
  call mloexpdep(0,dep(1:ifull,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
  call mloexport(0,w_t(:,ii),ii,0)
  call mlofill(w_t(:,ii))
  call mloexport(1,w_s(:,ii),ii,0)
  call mlofill(w_s(:,ii))
  call mloexport(2,w_u(:,ii),ii,0)
  call mloexport(3,w_v(:,ii),ii,0)
end do
call mloexport(4,w_e,0,0)
call bounds(dep)
call bounds(ps)

dep=max(dep,1.E-3)
dz=max(dz,1.E-3/real(wlev))

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
call mloexpdensity(rho,w_t,w_s,dep(1:ifull,:),0)
rhobar(1:ifull,1)=rho(:,1)*dz(:,1)
dzbar(:,1)=dz(:,1)
do ii=2,wlev
  rhobar(1:ifull,ii)=rhobar(1:ifull,ii-1)+rho(:,ii)*dz(:,ii)
  dzbar(:,ii)=dzbar(:,ii-1)+dz(:,ii)
end do
rhobar(1:ifull,:)=rhobar(1:ifull,:)/dzbar
call bounds(rhobar)

drdz(:,1)=(rhobar(1:ifull,2)-rhobar(1:ifull,1))/max(dep(1:ifull,2)-dep(1:ifull,1),1.E-20)
do ii=2,wlev-1
  drdz(:,ii)=(rhobar(1:ifull,ii+1)-rhobar(1:ifull,ii-1))/max(dep(1:ifull,ii+1)-dep(1:ifull,ii-1),1.E-20)
end do
drdz(:,wlev)=(rhobar(1:ifull,wlev)-rhobar(1:ifull,wlev-1))/max(dep(1:ifull,wlev)-dep(1:ifull,wlev-1),1.E-20)

! simple treatment of land-sea mask and bathmetry with masking gradients
drdx=0.
drdy=0.
do ii=1,wlev-1
  where(wtr(1:ifull).and.wtr(ie).and.wtr(iw).and.(dep(ie,wlev).le.dep(1:ifull,ii)).and.(dep(iw,wlev).le.dep(1:ifull,ii)))
    drdx(:,ii)=(rhobar(ie,ii)-rhobar(iw,ii))*0.5*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(ie).and.(dep(ie,wlev).le.dep(1:ifull,ii)))
    drdx(:,ii)=(rhobar(ie,ii)-rhobar(1:ifull,ii))*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(iw).and.(dep(iw,wlev).le.dep(1:ifull,ii)))
    drdx(:,ii)=(rhobar(1:ifull,ii)-rhobar(iw,ii))*em(1:ifull)/ds-dzdx(:,ii)*drdz(:,ii)
  end where
  where(wtr(1:ifull).and.wtr(in).and.wtr(is).and.(dep(in,wlev).le.dep(1:ifull,ii)).and.(dep(is,wlev).le.dep(1:ifull,ii)))
    drdy(:,ii)=(rhobar(in,ii)-rhobar(is,ii))*0.5*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(in).and.(dep(in,wlev).le.dep(1:ifull,ii)))
    drdy(:,ii)=(rhobar(in,ii)-rhobar(1:ifull,ii))*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  elsewhere(wtr(1:ifull).and.wtr(is).and.(dep(is,wlev).le.dep(1:ifull,ii)))
    drdy(:,ii)=(rhobar(1:ifull,ii)-rhobar(is,ii))*em(1:ifull)/ds-dzdy(:,ii)*drdz(:,ii)
  end where
end do

! coriolis terms
xp=1.
xm=0.
where (wtr(1:ifull))
  xp=1.+(dt*f(1:ifull))**2
  xm=dt*f(1:ifull)
end where

! Split U and V coriolis and pressure gradient terms (slow terms without free surface)
do ii=1,wlev
  nu(1:ifull,ii)=(w_u(:,ii)+xm*w_v(:,ii)-dt*grav*dep(1:ifull,ii)*(drdx(:,ii)+xm*drdy(:,ii))/rho(:,ii) &
                                              -dt*(dpsdx+xm*dpsdy)/rho(:,ii))/xp
  nv(1:ifull,ii)=(w_v(:,ii)-xm*w_u(:,ii)-dt*grav*dep(1:ifull,ii)*(drdy(:,ii)-xm*drdx(:,ii))/rho(:,ii) &
                                              -dt*(dpsdy-xm*dpsdx)/rho(:,ii))/xp
end do
w_u=nu(1:ifull,:)
w_v=nv(1:ifull,:)

if (.not.allocated(oldu1)) then
  allocate(oldu1(ifull,wlev))
  allocate(oldv1(ifull,wlev))
  oldu1=w_u
  oldv1=w_v
end if

! prep for conservation fix
! ...

do l=1,lmax ! predictor-corrector loop

  if (l.eq.1) then
    nuh(1:ifull,:)=1.5*w_u-0.5*oldu1 ! U at t+1/2
    nvh(1:ifull,:)=1.5*w_v-0.5*oldv1 ! V at t+1/2
  else
    nuh(1:ifull,:)=0.5*nu(1:ifull,:)+0.5*w_u ! U at t+1/2
    nvh(1:ifull,:)=0.5*nv(1:ifull,:)+0.5*w_v ! V at t+1/2
  end if

  ! Calculate depature points
  call mlodeps(nuh,nvh,nface,xg,yg,x3d,y3d,z3d)

! Vertical advection
! call mlovadv

! Convert (u,v) to cartesian coordinates (U,V,W)
  do ii=1,wlev
    cou(1:ifull,ii)=ax(1:ifull)*w_u(:,ii)+bx(1:ifull)*w_v(:,ii)
    cov(1:ifull,ii)=ay(1:ifull)*w_u(:,ii)+by(1:ifull)*w_v(:,ii)
    cow(1:ifull,ii)=az(1:ifull)*w_u(:,ii)+bz(1:ifull)*w_v(:,ii)
  end do

! Horizontal advection for U,V,W
  call bounds(cou)
  call bounds(cov)
  call bounds(cow)
  call mloints(cou,intsch,nface,xg,yg,2)
  call mloints(cov,intsch,nface,xg,yg,2)
  call mloints(cow,intsch,nface,xg,yg,2)

! Rotate vector to arrival point
  call mlorot(cou,cov,cow,x3d,y3d,z3d)

! Convert (U,V,W) back to conformal cubic coordinates
  do ii=1,wlev
    where (wtr(1:ifull))
      nu(1:ifull,ii)=ax(1:ifull)*cou(1:ifull,ii)+ay(1:ifull)*cov(1:ifull,ii)+az(1:ifull)*cow(1:ifull,ii)
      nv(1:ifull,ii)=bx(1:ifull)*cou(1:ifull,ii)+by(1:ifull)*cov(1:ifull,ii)+bz(1:ifull)*cow(1:ifull,ii)
    elsewhere
      nu(1:ifull,ii)=0.
      nv(1:ifull,ii)=0.
    end where
  end do
  
  ! Horizontal advection for T,S
  nt(1:ifull,:)=w_t
  ns(1:ifull,:)=w_s
  call bounds(nt)
  call bounds(ns)
  call mloints(nt,intsch,nface,xg,yg,2)
  call mloints(ns,intsch,nface,xg,yg,5)
  ns=max(ns,0.)
  
  ! fix for sigma levels (need to include map factor...)
  ! ...

  ! stagger nu and nv
  call mlostaguv(nu(1:ifull,:),nv(1:ifull,:),nuh(1:ifull,:),nvh(1:ifull,:))
  nu(1:ifull,:)=nuh(1:ifull,:)
  nv(1:ifull,:)=nvh(1:ifull,:)

  cou(1:ifull,:)=nu(1:ifull,:)
  cov(1:ifull,:)=nv(1:ifull,:)
  neta(1:ifull)=w_e
  
  ! Now solve for the coupled eta, nu and nv equations.
  ! This should be replaced with AB approach and SOR
  do ll=1,llmax

    ! update free surface height
    suu=0.
    svv=0.
    do ii=1,wlev
        suu(1:ifull)=suu(1:ifull)+dz(:,ii)*nu(1:ifull,ii)
        svv(1:ifull)=svv(1:ifull)+dz(:,ii)*nv(1:ifull,ii)
    end do
    where (.not.wtr(1:ifull))
      suu=0.
      svv=0.
    end where
    call boundsuv(suu,svv)
    div(:)=(suu(1:ifull)-suu(iwu)+svv(1:ifull)-svv(isv))*em(1:ifull)/ds
    neta(1:ifull)=alpha*(w_e-dt*div)/(1.+dt*div/dzbar(:,wlev))+(1.-alpha)*neta(1:ifull)
    neta(1:ifull)=max(neta(1:ifull),-dzbar(:,wlev))
    where (.not.wtr(1:ifull))
      neta(1:ifull)=0.
    end where
    call bounds(neta)

    dedx=0.
    dedy=0.
    where(wtr(1:ifull).and.wtr(ie))
      dedx(:)=(neta(ie)-neta(1:ifull))*em(1:ifull)/ds
    end where
    where(wtr(1:ifull).and.wtr(in))
      dedy(:)=(neta(in)-neta(1:ifull))*em(1:ifull)/ds
    end where

    ! update fast pressure gradient terms (predictor-corrector)
    do ii=1,wlev
      nu(1:ifull,ii)=cou(1:ifull,ii)-dt*grav*dep(1:ifull,ii)*rhobar(:,ii)*dedx/(rho(:,ii)*dzbar(:,wlev))
      nv(1:ifull,ii)=cov(1:ifull,ii)-dt*grav*dep(1:ifull,ii)*rhobar(:,ii)*dedy/(rho(:,ii)*dzbar(:,wlev))
      where (.not.wtr(1:ifull))
        nu(1:ifull,ii)=0.
        nv(1:ifull,ii)=0.
      end where
    end do
    
  end do

  ! unstagger nu and nv
  call mlounstaguv(nu(1:ifull,:),nv(1:ifull,:),nuh(1:ifull,:),nvh(1:ifull,:))

  ! update remaining pressure term
  do ii=1,wlev
    nu(1:ifull,ii)=nuh(1:ifull,ii)-dt*grav*dep(1:ifull,ii)*drdx(:,ii)*neta(1:ifull)/(rho(:,ii)*dzbar(:,wlev))
    nv(1:ifull,ii)=nvh(1:ifull,ii)-dt*grav*dep(1:ifull,ii)*drdy(:,ii)*neta(1:ifull)/(rho(:,ii)*dzbar(:,wlev))
  end do
  
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
call mloimport(4,neta(1:ifull),0,0)

! ADVECT ICE --------------------------------------------------------
!do ii=1,wlev
!  call mloexpice(w_u(:,1),8,0)
!  call mloexpice(w_v(:,1),9,0)
!end do

!if (.not.allocated(oldiceu1)) then
!  allocate(oldiceu1(ifull))
!  allocate(oldicev1(ifull))
!  oldiceu1=w_u(:,1)
!  oldicev1=w_v(:,1)
!end if

!do l=1,lmax
!
!  if (l.eq.1) then
!    nuh(1:ifull,1)=1.5*w_u(:,1)-0.5*oldiceu1 ! U at t+1/2
!    nvh(1:ifull,1)=1.5*w_v(:,1)-0.5*oldicev1 ! V at t+1/2
!  else
!    nuh(1:ifull,1)=0.5*nu(1:ifull,1)+0.5*w_u(:,1) ! U at t+1/2
!    nvh(1:ifull,1)=0.5*nv(1:ifull,1)+0.5*w_v(:,1) ! V at t+1/2
!  end if
!
!  ! Calculate depature points
!  call mlodeps(nuh(:,1),nvh(:,1),nface(:,1),xg(:,1),yg(:,1))
!
!! Account for ice deformation
!!...
!
!! Horizontal advector for ice fraction
!  nt(:,1)=i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newfracice=nt(:,1)
!
!! Horizontal advector for ice volume
!  nt(:,1)=i_dic*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newdic=nt(:,1)/newfracice
!
!! Horizontal advector for snow volume
!  nt(:,1)=i_dsn*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newdsn=nt(:,1)/newfracice
!
!! Horizontal advector for energy store
!  nt(:,1)=i_sto*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newsto=nt(:,1)/newfracice
!
!! Horizontal advector for temperatures 
!  nt(:,1)=i_tsurf*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newtsurf=nt(:,1)/newfracice
!  nt(:,1)=i_t0*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt0=nt(:,1)/newfracice
!  nt(:,1)=i_t1*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt1=nt(:,1)/newfracice
!  nt(:,1)=i_t2*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt2=nt(:,1)/newfracice
!  nt(:,1)=i_t3*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt3=nt(:,1)/newfracice

!
!end do

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(ubar,vbar,nface,xg,yg,x3d,y3d,z3d)

use cc_mpi
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer ii,intsch,n
integer, dimension(ifull,wlev), intent(out) :: nface
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real, dimension(ifull,wlev), intent(out) :: xg,yg
real*8, dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc
real, dimension(ifull+iextra,wlev) :: temp
integer, parameter :: nguess = 2

intsch=mod(ktau,2)

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
do ii=1,wlev
  call mlotoij5(x3d(:,ii),y3d(:,ii),z3d(:,ii),nface(:,ii),xg(:,ii),yg(:,ii))
end do
! Share off processor departure points.
call deptsync(nface,xg,yg)

do n=1,nguess
  temp(1:ifull,:) = uc(1:ifull,:)
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    x3d(1:ifull,ii) = x(1:ifull) - 0.5*(uc(1:ifull,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = vc(1:ifull,:)
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    y3d(1:ifull,ii) = y(1:ifull) - 0.5*(vc(1:ifull,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = wc(1:ifull,:)
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    z3d(1:ifull,ii) = z(1:ifull) - 0.5*(wc(1:ifull,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do

  do ii=1,wlev
    call mlotoij5(x3d(:,ii),y3d(:,ii),z3d(:,ii),nface(:,ii),xg(:,ii),yg(:,ii))
  end do
  !     Share off processor departure points.
  call deptsync(nface,xg,yg)
end do

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices
! This code is from depts.f

subroutine mlotoij5(x3d,y3d,z3d,nface,xg,yg)

use bigxy4_m
use cc_mpi
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmgeom.h'

integer loop,iq,i,j,is,js
integer, dimension(ifull), intent(out) :: nface
real, dimension(ifull), intent(out) :: xg,yg
real, dimension(ifull) :: xstr,ystr,zstr
real, dimension(ifull) :: denxyz,xd,yd,zd
real, dimension(ifull) :: ri,rj
real dxx,dxy,dyx,dyy
real*8, dimension(ifull), intent(inout) :: x3d,y3d,z3d
real*8, dimension(ifull) :: den
real*8 alf,alfonsch
real*8, parameter :: one = 1.
integer, parameter :: nmaploop = 3

!     if necessary, transform (x3d, y3d, z3d) to equivalent
!     coordinates (xstr, ystr, zstr) on regular gnomonic panels
if(schmidt.eq.1.)then
   xstr=x3d
   ystr=y3d
   zstr=z3d
else      ! (schmidt.ne.1.)
   alf=(one-schmidt**2)/(one+schmidt**2)
   alfonsch=2.*schmidt/(one+schmidt**2)  ! same but bit more accurate
   den=one-alf*z3d ! to force real*8
   xstr=x3d*(alfonsch/den)
   ystr=y3d*(alfonsch/den)
   zstr=    (z3d-alf)/den
endif     ! (schmidt.ne.1.)

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1
denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
xd=xstr/denxyz
yd=ystr/denxyz
zd=zstr/denxyz

where (abs(xstr-denxyz).lt.1.E-6)
  nface(:)    =0
  xg(:) =      yd
  yg(:) =      zd
elsewhere (abs(xstr+denxyz).lt.1.E-6)
  nface(:)    =3
  xg(:) =     -zd
  yg(:) =     -yd
elsewhere (abs(zstr-denxyz).lt.1.E-6)
  nface(:)    =1
  xg(:) =      yd
  yg(:) =     -xd
elsewhere (abs(zstr+denxyz).lt.1.E-6)
  nface(:)    =4
  xg(:) =      xd
  yg(:) =     -yd
elsewhere (abs(ystr-denxyz).lt.1.E-6)
  nface(:)    =2
  xg(:) =     -zd
  yg(:) =     -xd
elsewhere
  nface(:)    =5
  xg(:) =      xd
  yg(:) =      zd
end where

!     use 4* resolution grid il --> 4*il
xg=min(max(-.99999,xg),.99999)
yg=min(max(-.99999,yg),.99999)
!      first guess for ri, rj and nearest i,j
ri=1.+(1.+xg)*2.*real(il_g)
rj=1.+(1.+yg)*2.*real(il_g)
do loop=1,nmaploop
  do iq=1,ifull
    i=nint(ri(iq))
    j=nint(rj(iq))
    is=nint(sign(1.,ri(iq)-real(i)))
    js=nint(sign(1.,rj(iq)-real(j)))
!       predict new value for ri, rj
    dxx=xx4(i+is,j)-xx4(i,j)
    dyx=xx4(i,j+js)-xx4(i,j)
    dxy=yy4(i+is,j)-yy4(i,j)
    dyy=yy4(i,j+js)-yy4(i,j)       
    den(iq)=dxx*dyy-dyx*dxy
    ri(iq)=real(i)+real(is)*((xg(iq)-xx4(i,j))*dyy-(yg(iq)-yy4(i,j))*dyx)/den(iq)
    rj(iq)=real(j)+real(js)*((yg(iq)-yy4(i,j))*dxx-(xg(iq)-xx4(i,j))*dxy)/den(iq)
  end do
enddo  ! loop loop
!      expect xg, yg to range between .5 and il+.5
xg=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
yg=.25*(rj+3.) -.5  ! -.5 for stag

return
end subroutine mlotoij5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate depature points for semi-Lagrangian advection
! This code is from ints.f

subroutine mloints(s,intsch,nface,xg,yg,nfield)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'    ! has mh_bs
include 'mpif.h'

integer, intent(in) :: nfield,intsch
integer idel,iq,jdel,nn
integer i,j,k,n,ind,ip,jp,iproc,ierr
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(ifull+iextra,wlev), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev) :: sx
real, dimension(4) :: r
real a3,a4,c1,c2,c3,c4,cmax,cmin,sss,xxg,yyg
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do k=1,wlev

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,k) = s(ind(i,j,n),k)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop
            
!   this is intsb           EW interps done first
!   first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do n=1,npan
      do j=1,jpan
        sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
        sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
        sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
        sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,k) = s(is(ind(i,1,n)),k)
        sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
        sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
        sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
      enddo            ! i loop
!        for ew interpolation, sometimes need (different from ns):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,k) = s(lwws(n),k)
      sx(0,0,n,k) = s(lws(n),k)
      sx(0,-1,n,k) = s(lwss(n),k)
      sx(ipan+1,0,n,k) = s(les(n),k)
      sx(ipan+2,0,n,k) = s(lees(n),k)
      sx(ipan+1,-1,n,k) = s(less(n),k)
      sx(-1,jpan+1,n,k) = s(lwwn(n),k)
      sx(0,jpan+2,n,k) = s(lwnn(n),k)
      sx(ipan+2,jpan+1,n,k) = s(leen(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lenn(n),k)
      sx(0,jpan+1,n,k)    = s(iwn(ind(1,jpan,n)),k)
      sx(ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k)
    enddo               ! n loop

    if(nfield<mh_bs)then
      do iq=1,ifull    ! non Berm-Stan option
!                 Convert face index from 0:npanels to array indices
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          s(iq,k) = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          s(iq,k) = ((1.-yyg)*((2.-yyg)*        &
                    ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                    -yyg*(1.+yyg)*r(4)/3.)      &
                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    else                ! (nfield<mh_bs)
      do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sss = ((1.-yyg)*((2.-yyg)*        &
                ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                 -yyg*(1.+yyg)*r(4)/3.)     &
                 +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
        s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
      enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
            
  end do                 ! k

! Loop over points that need to be calculated for other processes

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2) + 0.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)* &        ! MJT memory
            ((1.+yyg)*r(2)-yyg*r(1)/3.)              &
            -yyg*(1.+yyg)*r(4)/3.)                   &
            +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    end do              ! iproc loop
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sss = ((1.-yyg)*((2.-yyg)*        &
                ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                -yyg*(1.+yyg)*r(4)/3.)      &
                +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
        sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
      enddo            ! iq loop
    end do              ! iproc loop
  endif                  ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do k=1,wlev              ! A single main k loop uses cache better
    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,k)=s(ind(i,j,n),k)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop

    do n=1,npan
      do j=1,jpan
        sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
        sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
        sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
        sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,k) = s(is(ind(i,1,n)),k)
        sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
        sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
        sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
      enddo            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,k)=s(lsww(n),k)
      sx(0,0,n,k) = s(lsw(n),k)
      sx(0,-1,n,k) = s(lssw(n),k)
      sx(ipan+2,0,n,k) = s(lsee(n),k)
      sx(ipan+1,-1,n,k) = s(lsse(n),k)
      sx(-1,jpan+1,n,k) = s(lnww(n),k)
      sx(0,jpan+1,n,k) = s(lnw(n),k)
      sx(0,jpan+2,n,k) = s(lnnw(n),k)
      sx(ipan+2,jpan+1,n,k) = s(lnee(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lnne(n),k)
      sx(ipan+1,0,n,k)    = s(ise(ind(ipan,1,n)),k)
      sx(ipan+1,jpan+1,n,k) = s(ine(ind(ipan,jpan,n)),k)
    enddo               ! n loop

    if(nfield<mh_bs)then
      do iq=1,ifull    ! non Berm-Stan option
!       Convert face index from 0:npanels to array indices
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          s(iq,k) = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          s(iq,k) = ((1.-xxg)*((2.-xxg)*        &
                    ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                   -xxg*(1.+xxg)*r(4)/3.)       &
                   +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    else                ! (nfield<mh_bs)
      do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sss = ((1.-xxg)*((2.-xxg)*        &
                ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                -xxg*(1.+xxg)*r(4)/3.)      &
                +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
        s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
      enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
  end do                 ! k

! For other processes
  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2)+ 0.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*  & ! MJT memory
            ((1.+xxg)*r(2)-xxg*r(1)/3.)               &
            -xxg*(1.+xxg)*r(4)/3.)                    &
            +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    end do              ! iproc
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sss = ((1.-xxg)*((2.-xxg)*        &
                ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                -xxg*(1.+xxg)*r(4)/3.)      &
                +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
        sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
      enddo            ! iq loop
    end do              ! iproc
  endif                  ! (nfield<mh_bs)  .. else ..

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync(s)

return
end subroutine mloints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use mlo
use xyzinfo_m

implicit none

include 'newmpar.h'

integer k
real, dimension(ifull,wlev), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real*8, dimension(ifull,wlev), intent(in) :: x3d,y3d,z3d

do k=1,wlev
!         cross product n1xn2 into vec1
  vec1x = y3d(:,k)*z - y*z3d(:,k)
  vec1y = z3d(:,k)*x - z*x3d(:,k)
  vec1z = x3d(:,k)*y - x*y3d(:,k)
  denb = vec1x**2 + vec1y**2 + vec1z**2
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
  where (denb>1.e-4)
    vecdot = x3d(:,k)*x + y3d(:,k)*y + z3d(:,k)*z
    vec2x = x3d(:,k)*vecdot - x
    vec2y = y3d(:,k)*vecdot - y
    vec2z = z3d(:,k)*vecdot - z
    vec3x = x3d(:,k) - vecdot*x
    vec3y = y3d(:,k) - vecdot*y
    vec3z = z3d(:,k) - vecdot*z
    vdot1 = (vec1x*cou(:,k) + vec1y*cov(:,k) + vec1z*cow(:,k))/denb
    vdot2 = (vec2x*cou(:,k) + vec2y*cov(:,k) + vec2z*cow(:,k))/denb
    cou(:,k) = vdot1*vec1x + vdot2*vec3x
    cov(:,k) = vdot1*vec1y + vdot2*vec3y
    cow(:,k) = vdot1*vec1z + vdot2*vec3z
  end where
end do ! k

return
end subroutine mlorot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stagger u and v

subroutine mlostaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'

integer k,itn
real, dimension(ifull,wlev), intent(in) :: u,v
real, dimension(ifull,wlev), intent(out) :: uout,vout
real, dimension(ifull+iextra,wlev) :: uin,vin
real, dimension(ifull+iextra,wlev) :: ua,va,ud,vd
integer, parameter :: itnmax=3

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,wlev
  ud(1:ifull,k)= uin(iwu,k)/10.+uin(1:ifull,k)+uin(ieu,k)/2.
  vd(1:ifull,k)= vin(isv,k)/10.+vin(1:ifull,k)+vin(inv,k)/2.
enddo
call boundsuv(ud,vd)

do k=1,wlev
  ua(1:ifull,k)=ud(1:ifull,k)-ud(iwu,k)/2. ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(isv,k)/2. ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)

  do k=1,wlev
    uin(1:ifull,k)=(ud(1:ifull,k)-.5*ud(iwu,k)-ua(ieu,k)/10. +ua(iwwu,k)/4.)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-.5*vd(isv,k)-va(inv,k)/10. +va(issv,k)/4.)/.95
  enddo

  call boundsuv(uin,vin,nrows=2)
  do k=1,wlev
    ua(1:ifull,k)=(ud(1:ifull,k)-.5*ud(iwu,k)-uin(ieu,k)/10. +uin(iwwu,k)/4.)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-.5*vd(isv,k)-vin(inv,k)/10. +vin(issv,k)/4.)/.95
  end do
end do                 ! itn=1,itnmax

uout=ua(1:ifull,:)
vout=va(1:ifull,:)

return
end subroutine mlostaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unstagger u and v

subroutine mlounstaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'

integer k,itn
real, dimension(ifull,wlev), intent(in) :: u,v
real, dimension(ifull,wlev), intent(out) :: uout,vout
real, dimension(ifull+iextra,wlev) :: uin,vin
real, dimension(ifull+iextra,wlev) :: ua,va,ud,vd
integer, parameter :: itnmax=3

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,wlev
  ud(1:ifull,k)= uin(ieu,k)/10.+uin(1:ifull,k)+uin(iwu,k)/2.
  vd(1:ifull,k)= vin(inv,k)/10.+vin(1:ifull,k)+vin(isv,k)/2.
enddo
call boundsuv(ud,vd)

do k=1,wlev
  ua(1:ifull,k)=ud(1:ifull,k)-ud(ieu,k)/2. ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(inv,k)/2. ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)

  do k=1,wlev
    uin(1:ifull,k)=(ud(1:ifull,k)-.5*ud(ieu,k)-ua(iwu,k)/10. +ua(ieeu,k)/4.)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-.5*vd(inv,k)-va(isv,k)/10. +va(innv,k)/4.)/.95
  enddo
  call boundsuv(uin,vin,nrows=2)
  do k=1,wlev
    ua(1:ifull,k)=(ud(1:ifull,k)-.5*ud(ieu,k)-uin(iwu,k)/10. +uin(ieeu,k)/4.)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-.5*vd(inv,k)-vin(isv,k)/10. +vin(innv,k)/4.)/.95
  enddo
enddo                  ! itn=1,itnmax
      
uout=ua(1:ifull,:)
vout=va(1:ifull,:)

return
end subroutine mlounstaguv

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
subroutine mlofill(x)

use cc_mpi
use indices_m
use soil_m

implicit none

include 'newmpar.h'
include 'mpif.h'

integer globalc,localc,ierr,iq,lnum
real miss,lsum
real, dimension(ifull), intent(inout) :: x
real, dimension(ifull) :: yy
real, dimension(ifull+iextra) :: xx
logical, dimension(ifull+iextra) :: smap

miss=999999.

xx=miss
where (.not.land)
  xx(1:ifull)=x
end where
globalc=1

! technically we only need 1 pass of this fill to ensure all water points have a non-trival neighbour
do while (globalc.gt.0)
  call bounds(xx)
  smap=abs(xx-miss).lt.0.1
  yy=xx(1:ifull)
  do iq=1,ifull
    if (smap(iq)) then
      lsum=0.
      lnum=0
      if (.not.smap(in(iq))) then
        lsum=lsum+xx(in(iq))
        lnum=lnum+1
      end if
      if (.not.smap(ie(iq))) then
        lsum=lsum+xx(ie(iq))
        lnum=lnum+1
      end if
      if (.not.smap(is(iq))) then
        lsum=lsum+xx(is(iq))
        lnum=lnum+1
      end if
      if (.not.smap(iw(iq))) then
        lsum=lsum+xx(iw(iq))
        lnum=lnum+1
      end if
      if (lnum.gt.0) then
        yy(iq)=lsum/real(lnum)
      end if
    end if
  end do
  xx(1:ifull)=yy
  smap(1:ifull)=abs(xx(1:ifull)-miss).lt.0.1
  localc=count(smap(1:ifull))
  call MPI_AllReduce(localc,globalc,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
end do

x=xx(1:ifull)

return
end subroutine mlofill

end module mlodynamics