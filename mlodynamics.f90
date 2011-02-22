! These subroutines handle dynamics for the Mixed-Layer-Ocean model

module mlodynamics

implicit none

private
public mlodiffusion,mlorouter,watbdy

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
real, parameter :: k_smag = 2. ! 0.4 in mom3
real hdif
real, dimension(ifull+iextra) :: uc,vc,wc,ee
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: u,v
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc
logical, dimension(ifull+iextra) :: wtr

hdif=dt*(k_smag/pi)**2
ee=0.
where(land)
  ee=1.
end where
call bounds(ee)
wtr=ee.lt.0.5

if (.not.any(wtr(1:ifull))) return

emi=1./em(1:ifull)**2
do k=1,wlev
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
  cc=max(cc,1.E-10)
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

  ucc = ( uc(1:ifull)*emi +                    &
          xfact(1:ifull)*uc(ie) +              &
          xfact(iwu)*uc(iw) +                  &
          yfact(1:ifull)*uc(in) +              &
          yfact(isv)*uc(is) ) /                &
        ( emi + xfact(1:ifull) + xfact(iwu) +  &
          yfact(1:ifull)+yfact(isv) )
  vcc = ( vc(1:ifull)*emi +                    &
          xfact(1:ifull)*vc(ie) +              &
          xfact(iwu)*vc(iw) +                  &
          yfact(1:ifull)*vc(in) +              &
          yfact(isv)*vc(is) ) /                &
        ( emi + xfact(1:ifull) + xfact(iwu) +  &
          yfact(1:ifull)+yfact(isv) )
  wcc = ( wc(1:ifull)*emi +                    &
          xfact(1:ifull)*wc(ie) +              &
          xfact(iwu)*wc(iw) +                  &
          yfact(1:ifull)*wc(in) +              &
          yfact(isv)*wc(is) ) /                &
        ( emi + xfact(1:ifull) + xfact(iwu) +  &
          yfact(1:ifull) + yfact(isv) )
  u = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u,k,0)
  call mloimport(3,v,k,0)
   
  do i=0,1
    ee=0.
    call mloexport(i,ee(1:ifull),k,0)
    call bounds(ee)
    ff = ( ee(1:ifull)*emi +                    &
           xfact(1:ifull)*ee(ie) +              &
           xfact(iwu)*ee(iw) +                  &
           yfact(1:ifull)*ee(in) +              &
           yfact(isv)*ee(is) ) /                &
         ( emi + xfact(1:ifull) + xfact(iwu) +  &
           yfact(1:ifull)+yfact(isv))
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
    ff = ( ee(1:ifull)*emi +                    &
           xfact(1:ifull)*ee(ie) +              &
           xfact(iwu)*ee(iw) +                  &
           yfact(1:ifull)*ee(in) +              &
           yfact(isv)*ee(is) ) /                &
         ( emi + xfact(1:ifull) + xfact(iwu) +  &
           yfact(1:ifull)+yfact(isv))
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
!subroutine mlohadv
!
!use cc_mpi
!use map_m
!use mlo
!
!implicit none
!
!include 'newmpar.h'
!
!integer l,ii
!real, dimension(ifull+iextra) :: ee
!real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
!real, dimension(ifull+iextra,wlev) :: p,dep
!real, dimension(ifull,wlev) :: dzdx,dzdy
!logical, dimension(ifull) :: wtr
!
!Define land/sea mask
!ee=0.
!where(land)
!  ee=1.
!end where
!call bounds(ee)
!wtr=ee.lt.0.5
!
!dep(1:ifull,:)=depth
!
! Initialise arrays
!do ii=1,wlev
!  call mloexport(0,nt(1:ifull,ii),ii,0)
!  call mloexport(1,ns(1:ifull,ii),ii,0)
!  call mloexport(2,nu(1:ifull,ii),ii,0)
!  call mloexport(3,nv(1:ifull,ii),ii,0)
!end do
!
!do l=1,2 ! predictor-corrector loop
!
!  ! update bounds
!  do ii=1,wlev
!    call bounds(nu(:,ii))
!    call bounds(nv(:,ii))
!    call bounds(nt(:,ii))
!    call bounds(ns(:,ii))
!  end do
!
!  ! Calculate depth gradient
!  dzdx=0.
!  dzdy=0.
!  do ii=1,wlev
!    call bounds(dep(:,ii))
!    where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
!      dzdx(:,ii)=(dep(ie,ii)-dep(iw,ii))*0.5*em(:)/ds
!    elsewhere(wtr(1:ifull).and.wtr(ie))
!      dzdx(:,ii)=(dep(ie,ii)-dep(1:ifull,ii))*em(:)/ds
!    elsewhere(wtr(1:ifull).and.wtr(iw))
!      dzdx(:,ii)=(dep(1:ifull,ii)-dep(iw,ii)*em(:)/ds
!    end where
!    where(wtr(1:ifull).and.wtr(in).and.wtr(is))
!      dzdy(:,ii)=(dep(in,ii)-dep(is,ii))*0.5*em(:)/ds
!    elsewhere(wtr(1:ifull).and.wtr(in))
!      dzdy(:,ii)=(dep(in,ii)-dep(1:ifull,ii))*em(:)/ds
!    elsewhere(wtr(1:ifull).and.wtr(is))
!      dzdy(:,ii)=(dep(1:ifull,ii)-dep(is,ii)*em(:)/ds
!    end where
!  end do
!
!  ! Calculate pressure and gradient
!  do ii=1,wlev
!    rho(:,ii)=
!    p(:,ii)=grav*depth(1:ifull,ii)*rho(:,ii)
!    call bounds(p(:,ii))
!  end do
!  dpdz(:,1)=
!  do ii=2,wlev-1
!    dpdz(:,ii)=
!  end do
!  dpdz(:,wlev)=
!
!  dpdx=0.
!  dpdy=0.
!  do ii=1,wlev
!    where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
!      dpdx(:,ii)=(p(ie,ii)-p(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dpdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(ie))
!      dpdx(:,ii)=(p(ie,ii)-p(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dpdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(iw))
!      dpdx(:,ii)=(p(1:ifull,ii)-p(iw,ii)*em(:)/ds+dzdx(:,ii)*dpdz(:,ii)
!    end where
!    where(wtr(1:ifull).and.wtr(in).and.wtr(is))
!      dpdy(:,ii)=(p(in,ii)-p(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dpdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(in))
!      dpdy(:,ii)=(p(in,ii)-p(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dpdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(is))
!      dpdy(:,ii)=(p(1:ifull,ii)-p(is,ii)*em(:)/ds+dzdy(:,ii)*dpdz(:,ii)
!    end where
!  end do
!
!
!  ! Split U and V coriolis and pressure gradient terms
!  xp=1.+(dt*a_f)**2
!  xm=dt*a_f
!  do ii=1,wlev
!    ou=(w_u(:,ii)+xm*w_v(:,ii)-dt*(dpdx(:,ii)+xm*dpdy(:,ii))/rho(:,ii))/xp
!    ov=(w_v(:,ii)-xm*w_u(:,ii)-dt*(dpdy(:,ii)-xm*dpdx(:,ii))/rho(:,ii))/xp
!  end do
!
!  if (l.eq.1) then
!    nu=ou
!    nv=ov
!  end if
!
!  ! Convert (u,v) to cartesian coordinates (U,V,W) for advection terms
!  cou=ax*ou+bx*ov
!  cov=ay*ou+by*ov
!  cow=az*ou+bz*ov
!  cnu=ax*nu+bx*nv
!  cnv=ay*nu+by*nv
!  cnw=az*nu+bz*nv
!
!  ! Calculate vertical velocity
!  do ii=1,wlev
!    dudz(:,ii)=
!    dvdz(:,ii)=
!    dwdz(:,ii)=
!  end do
!
!  ! Calculate velocity gradients
!  dudx=0.
!  dvdx=0.
!  dwdx=0.
!  dudy=0.
!  dvdy=0.
!  dwdy=0.
!  do ii=1,wlev
!    where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
!      dudx(:,ii)=(cnu(ie,ii)-cnu(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dudz(:,ii)
!      dvdx(:,ii)=(cnv(ie,ii)-cnv(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dvdz(:,ii)
!      dwdx(:,ii)=(cnw(ie,ii)-cnw(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dwdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(ie))
!      dudx(:,ii)=(cnu(ie,ii)-cnu(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dudz(:,ii)
!      dvdx(:,ii)=(cnv(ie,ii)-cnv(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dvdz(:,ii)
!      dwdx(:,ii)=(cnw(ie,ii)-cnw(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dwdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(iw))
!      dudx(:,ii)=(cnu(1:ifull,ii)-cnu(iw,ii)*em(:)/ds+dzdx(:,ii)*dudz(:,ii)
!      dvdx(:,ii)=(cnv(1:ifull,ii)-cnv(iw,ii)*em(:)/ds+dzdx(:,ii)*dvdz(:,ii)
!      dwdx(:,ii)=(cnw(1:ifull,ii)-cnw(iw,ii)*em(:)/ds+dzdx(:,ii)*dwdz(:,ii)
!    end where
!    where(wtr(1:ifull).and.wtr(in).and.wtr(is))
!      dudy(:,ii)=(cnu(in,ii)-cnu(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dudz(:,ii)
!      dvdy(:,ii)=(cnv(in,ii)-cnv(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dvdz(:,ii)
!      dwdy(:,ii)=(cnw(in,ii)-cnw(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dwdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(in))
!      dudy(:,ii)=(cnu(in,ii)-cnu(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dudz(:,ii)
!      dvdy(:,ii)=(cnv(in,ii)-cnv(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dvdz(:,ii)
!      dwdy(:,ii)=(cnw(in,ii)-cnw(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dwdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(is))
!      dudy(:,ii)=(cnu(1:ifull,ii)-cnu(is,ii)*em(:)/ds+dzdy(:,ii)*dudz(:,ii)
!      dvdy(:,ii)=(cnv(1:ifull,ii)-cnv(is,ii)*em(:)/ds+dzdy(:,ii)*dvdz(:,ii)
!      dwdy(:,ii)=(cnw(1:ifull,ii)-cnw(is,ii)*em(:)/ds+dzdy(:,ii)*dwdz(:,ii)
!    end where
!  end do
!
!  ! Stagger u and v
!  do i=1,wlev
!    snu(:,ii)=
!    snv(:,ii)=
!  end do
!
!  ! Horizontal terms
!  do ii=1,wlev
!    cuu(:)=cou(:,ii)+dt*(snu(:,ii)*dudx(:,ii)+snv(:,ii)*dudy(:,ii))    ! need cartesian version
!    cvv(:)=cov(:,ii)+dt*(snu(:,ii)*dvdx(:,ii)+snv(:,ii)*dvdy(:,ii))
!    cww(:)=cow(:,ii)+dt*(snu(:,ii)*dwdx(:,ii)+snv(:,ii)*dwdy(:,ii))
!    uu(:,ii)=ax*cuu+ay*cvv+az*cww
!    vv(:,ii)=bx*cuu+by*cvv+bz*cww
!  end do
!
!
!  ! Vertical terms
!  w_w(:,ii)=
!  w_u(:,ii)=w_u(:,ii)+dt*w_w(:,ii)*dudz(:,ii)
!  w_v(:,ii)=w_v(:,ii)+dt*w_w(:,ii)*dvdz(:,ii)
!
!  ! Calculate vertical gradients for temperature and salinity
!  dtdz(:,1)=
!  dsdz(:,1)=
!  do ii=1,wlev
!    dtdz(:,ii)=
!    dsdz(:,ii)=
!  end do
!  dtdz(:,wlev)=
!  dsdx(:,wlev)=
!
!  ! Calculate temperature and salinity gradients
!  dtdx=0.
!  dsdx=0.
!  dtdy=0.
!  dsdy=0.
!  do ii=1,wlev
!    where(wtr(1:ifull).and.wtr(ie).and.wtr(iw))
!      dtdx(:,ii)=(nt(ie,ii)-nt(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dtdz(:,ii)
!      dsdx(:,ii)=(ns(ie,ii)-ns(iw,ii))*0.5*em(:)/ds+dzdx(:,ii)*dsdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(ie))
!      dtdx(:,ii)=(nt(ie,ii)-nt(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dtdz(:,ii)
!      dsdx(:,ii)=(ns(ie,ii)-ns(1:ifull,ii))*em(:)/ds+dzdx(:,ii)*dsdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(iw))
!      dtdx(:,ii)=(nt(1:ifull,ii)-nt(iw,ii)*em(:)/ds+dzdx(:,ii)*dtdz(:,ii)
!      dsdx(:,ii)=(ns(1:ifull,ii)-ns(iw,ii)*em(:)/ds+dzdx(:,ii)*dsdz(:,ii)
!    end where
!    where(wtr(1:ifull).and.wtr(in).and.wtr(is))
!      dtdy(:,ii)=(nt(in,ii)-nt(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dtdz(:,ii)
!      dsdy(:,ii)=(ns(in,ii)-ns(is,ii))*0.5*em(:)/ds+dzdy(:,ii)*dsdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(in))
!      dtdy(:,ii)=(nt(in,ii)-nt(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dtdz(:,ii)
!      dsdy(:,ii)=(ns(in,ii)-ns(1:ifull,ii))*em(:)/ds+dzdy(:,ii)*dsdz(:,ii)
!    elsewhere(wtr(1:ifull).and.wtr(is))
!      dtdy(:,ii)=(nt(1:ifull,ii)-nt(is,ii)*em(:)/ds+dzdy(:,ii)*dtdz(:,ii)
!      dsdy(:,ii)=(ns(1:ifull,ii)-ns(is,ii)*em(:)/ds+dzdy(:,ii)*dsdz(:,ii)
!    end where
!  end do
!
!  ! Advect temperature and salinity
!  do ii=1,wlev
!    tt(:,ii)=w_t(:,ii)+dt*(su(:,ii)*dtdx(:,ii)+sv(:,ii)*dtdy(:,ii))
!    ss(:,ii)=w_s(:,ii)+dt*(su(:,ii)*dsdx(:,ii)+sv(:,ii)*dsdy(:,ii))
!  end do
!
!  w_t(:,ii)=w_t(:,ii)+dt*w_w(:,ii)*dtdz(:,ii)
!  w_s(:,ii)=w_s(:,ii)+dt*w_w(:,ii)*dsdz(:,ii)
!
!  ! Update rho for corrector step
!  nu=uu
!  nv=vv
!  nt=tt
!  ns=ss
!end do
!
!do k=1,wlev
!  call mloimport(0,tt(:,ii),k,0)
!  call mloimport(1,ss(:,ii),k,0)
!  call mloimport(2,uu(:,ii),k,0)
!  call mloimport(3,vv(:,ii),k,0)
!end do
!
!return
!end subroutine mlohadv

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