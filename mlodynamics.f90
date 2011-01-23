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
use cc_mpi
use indices_m
use map_m
use mlo

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer i,iq
integer, dimension(ifull,4) :: xp
real, dimension(ifull) :: newwat
real, dimension(ifull,4) :: dp,slope,mslope,vel,flow

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
  if (all(slope(iq,:).lt.-1.E-10)) then
    newwat(iq)=newwat(iq)-watbdy(iq) ! sink - assume sub grid scale lake is present
  end if
end do

watbdy(1:ifull)=max(newwat,0.)

return
end subroutine mlorouter

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