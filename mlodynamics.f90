! These subroutines handle dynamics for the MLO

! This subroutine processes horizontal diffusion
subroutine mlodiffusion

use cc_mpi
use indices_m
use map_m
use mlo
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer k,i
real, parameter :: k_smag = 0.4
real hdif
real, dimension(ifull+iextra) :: uc,vc,wc,ee
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: u,v
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc

hdif=dt*(k_smag/pi)**2

emi=1./em(1:ifull)**2
do k=1,wlev
  ! neglect terrain following component for now
  ! u=ax*uc+ay*vc+az*wc
  ! dudx=u(ie)-u(iw)=ax*(uc(ie)-uc(iw))+ay*(vc(ie)-vc(iw))+az*(wc(ie)-wc(iw))
  ! dudy=u(in)-u(is)=ax*(uc(in)-uc(is))+ay*(vc(in)-vc(is))+az*(wc(in)-wc(is))
  u=0.
  call mloexport(2,u,k,0)
  call mlofill(u)
  uc(1:ifull)=ax(1:ifull)*u
  vc(1:ifull)=ay(1:ifull)*u
  wc(1:ifull)=az(1:ifull)*u
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)
  dudx=(ax(1:ifull)*(uc(ie)-uc(iw))  &
       +ay(1:ifull)*(vc(ie)-vc(iw))  &
       +az(1:ifull)*(wc(ie)-wc(iw))) &
       *0.5*em(1:ifull)/ds
  dudy=(ax(1:ifull)*(uc(in)-uc(is))  &
       +ay(1:ifull)*(vc(in)-vc(is))  &
       +az(1:ifull)*(wc(in)-wc(is))) &
       *0.5*em(1:ifull)/ds

  ! v=bx*uc+by*vc+bz*wc
  ! dvdx=v(ie)-v(iw)=bx*(uc(ie)-uc(iw))+by*(vc(ie)-vc(iw))+bz*(wc(ie)-wc(iw))
  ! dvdy=v(in)-v(is)=bx*(uc(in)-uc(is))+by*(vc(in)-vc(is))+bz*(wc(in)-wc(is))
  v=0.
  call mloexport(3,v,k,0)
  call mlofill(v)
  uc(1:ifull)=bx(1:ifull)*v
  vc(1:ifull)=by(1:ifull)*v
  wc(1:ifull)=bz(1:ifull)*v
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)
  dvdx=(bx(1:ifull)*(uc(ie)-uc(iw))  &
       +by(1:ifull)*(vc(ie)-vc(iw))  &
       +bz(1:ifull)*(wc(ie)-wc(iw))) &
       *0.5*em(1:ifull)/ds
  dvdy=(bx(1:ifull)*(uc(in)-uc(is))  &
       +by(1:ifull)*(vc(in)-vc(is))  &
       +bz(1:ifull)*(wc(in)-wc(is))) &
       *0.5*em(1:ifull)/ds

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
  xfact(1:ifull) = (t_kh(ie)+t_kh)*.5
  yfact(1:ifull) = (t_kh(in)+t_kh)*.5
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
    call mloexport(i,ee,k,0)
    call mlofill(ee)
    call bounds(ee)
    ff = ( ee(1:ifull)*emi +                    &
           xfact(1:ifull)*ee(ie) +              &
           xfact(iwu)*ee(iw) +                  &
           yfact(1:ifull)*ee(in) +              &
           yfact(isv)*ee(is) ) /                &
         ( emi + xfact(1:ifull) + xfact(iwu) +  &
           yfact(1:ifull)+yfact(isv))
    call mloimport(i,ee,k,0)
  end do
  
end do

return
end subroutine mlodiffusion

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
real, dimension(ifull+iextra) :: xx,yy
logical, dimension(ifull+iextra) :: smap

miss=999999.

xx=miss
where (.not.land)
  xx(1:ifull)=x
end where
globalc=1

! technically we only need 1 pass of this fill to ensure all water points have a non-trival neighbour
!do while (globalc.gt.0)
  call bounds(xx)
  smap=abs(xx-miss).lt.0.1
  yy=xx
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
  xx=yy
  !smap=abs(xx-miss).lt.0.1
  !localc=count(smap(1:ifull))
  !call MPI_AllReduce(localc,globalc,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!end do

return
end subroutine mlofill