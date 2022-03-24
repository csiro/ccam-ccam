! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module mlovadvtvd

implicit none

private
public mlovadv, mlontvd

integer, save      :: mlontvd     = 1       ! Vertical advection limiter (0=MC, 1=Superbee)

contains
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs vertical advection based on JLMs TVD scheme

subroutine mlovadv(dtin,ww,uu,vv,ss,tt,mm,depdum,idzdum,wtr,cnum)

use cc_mpi
use mlo
use newmpar_m

implicit none

integer, intent(in) :: cnum
integer ii,iq,its_g
integer, dimension(ifull) :: its
real, intent(in) :: dtin
real, dimension(ifull) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depdum,idzdum
real, dimension(:,:), intent(inout) :: uu,vv,ss,tt,mm
real, dimension(ifull,wlev) :: dzdum
logical, dimension(ifull+iextra,wlev), intent(in) :: wtr

call START_LOG(watervadv_begin)

dzdum = max(idzdum(:,:),1.E-10)
  
! reduce time step to ensure stability
dtnew(:) = dtin
do ii = 1,wlev-1
  do iq = 1,ifull    
    if ( wtr(iq,ii) ) then
      ! this trick works if dzdum(iq,ii)<dzdum(iq,ii+1)
      dtnew(iq)=min(dtnew(iq),0.3*dzdum(iq,ii)/max(abs(ww(iq,ii)),1.E-12))
    end if
  end do
end do
its(:) = int(dtin/(dtnew(:)+0.01))+1
dtnew(:) = dtin/real(its(:))

its_g=maxval(its(:))
if (its_g>500) then
  write(6,*) "MLOVERT myid,cnum,its_g",myid,cnum,its_g
end if

#ifdef GPU
!$omp target data map(to:its,dtnew,ww,depdum,dzdum)
#else
!$omp parallel sections
#endif
!$acc data create(its,dtnew,ww,depdum,dzdum)
!$acc update device(its,dtnew,ww,depdum,dzdum)

#ifndef GPU
!$omp section
#endif
call mlotvd(its,dtnew,ww,uu,depdum,dzdum)
#ifndef GPU
!$omp section
#endif
call mlotvd(its,dtnew,ww,vv,depdum,dzdum)
#ifndef GPU
!$omp section
#endif
call mlotvd(its,dtnew,ww,ss,depdum,dzdum)
#ifndef GPU
!$omp section
#endif
call mlotvd(its,dtnew,ww,tt,depdum,dzdum)
#ifndef GPU
!$omp section
#endif
call mlotvd(its,dtnew,ww,mm,depdum,dzdum)

#ifdef GPU
!$omp end target data
#else
!$omp end parallel sections
#endif
!$acc wait
!$acc end data
  
ss(1:ifull,:)=max(ss(1:ifull,:),0.)
tt(1:ifull,:)=max(tt(1:ifull,:),-wrtemp)

call END_LOG(watervadv_end)

return
end subroutine mlovadv

subroutine mlotvd(its,dtnew,ww,uu,depdum,dzdum)

use cc_acc, only : async_length
use mlo
use newmpar_m

implicit none

integer ii,i,iq,kp,kx
integer, save :: async_counter = -1
integer, dimension(ifull), intent(in) :: its
real, dimension(ifull), intent(in) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depdum,dzdum
real, dimension(:,:), intent(inout) :: uu
real, dimension(ifull,0:wlev) :: ff
real, dimension(ifull,0:wlev) :: delu
real fl,fh,cc,rr

! f=(w*u) at half levels
! du/dt = u*dw/dz-df/dz = -w*du/dz

async_counter = mod( async_counter+1, async_length )

if ( mlontvd==0 ) then ! MC

#ifdef GPU
  !$omp target enter data map(to:uu) map(alloc:delu,ff)
#endif
  !$acc enter data create(uu,delu,ff) async(async_counter)
  !$acc update device(uu,dzdum) async(async_counter)

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq)
#endif
#else
  !$acc parallel loop collapse(2) present(delu,uu,ff,dzdum) async(async_counter)
#endif
  do ii = 1,wlev-1
    do iq = 1,ifull
      ff(iq,ii) = 0.  
      if ( dzdum(iq,ii+1)>1.e-4 ) then
        delu(iq,ii) = uu(iq,ii+1) - uu(iq,ii)
      else
        delu(iq,ii) = 0.  
      end if    
    end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
  !$omp target teams distribute parallel do schedule(static) private(iq)
#endif
#else
  !$acc end parallel loop
  !$acc parallel loop present(ff,delu) async(async_counter)
#endif
  do iq = 1,ifull
    ff(iq,0) = 0.
    ff(iq,wlev) = 0.
    delu(iq,0) = 0.
    delu(iq,wlev) = 0.
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

! TVD part
#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq,kp,kx,rr,fl,cc,fh)
#endif
#else
  !$acc parallel loop collapse(2) present(ff,ww,delu,uu,dtnew,depdum,dzdum) async(async_counter)
#endif
  do ii = 1,wlev-1
    do iq = 1,ifull
      ! +ve ww is downwards to the ocean floor
      kp = nint(sign(1.,ww(iq,ii)))
      kx = ii+(1-kp)/2 !  k for ww +ve,  k+1 for ww -ve
      rr=delu(iq,ii-kp)/(delu(iq,ii)+sign(1.E-20,delu(iq,ii)))
      fl=ww(iq,ii)*uu(iq,kx)
      cc = max(0.,min(2.*rr, 0.5+0.5*rr,2.)) ! MC
      fh = ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1))             &
        - 0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq) &
        /max(depdum(iq,ii+1)-depdum(iq,ii),1.E-10)
      if ( dzdum(iq,ii+1)>1.e-4 ) then
        ff(iq,ii) = fl + cc*(fh-fl)
      end if
     !ff(iq,ii)=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1)) ! explicit
   end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq)
#endif
#else
  !$acc end parallel loop
  !$acc parallel loop collapse(2) present(ff,uu,ww,dtnew,dzdum) async(async_counter)
#endif
  do ii = 1,wlev
    do iq = 1,ifull
      if ( dzdum(iq,ii)>1.e-4 ) then  
        uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1))   &
                           -ff(iq,ii)+ff(iq,ii-1))/dzdum(iq,ii)
      end if   
    end do  
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do schedule(static) private(iq,i,ii,kp,kx,rr,fl,cc,fh)
#endif
#else
  !$acc parallel loop present(its,delu,uu,ww,ff,dtnew,depdum,dzdum) async(async_counter)
#endif
  do iq = 1,ifull
    do i = 2,its(iq)
      do ii=1,wlev-1
        if ( dzdum(iq,ii+1)>1.e-4 ) then
          delu(iq,ii) = uu(iq,ii+1) - uu(iq,ii)
        else
          delu(iq,ii) = 0.  
        end if 
      end do
      ! TVD part
      do ii=1,wlev-1
        ! +ve ww is downwards to the ocean floor
        kp = nint(sign(1.,ww(iq,ii)))
        kx = ii+(1-kp)/2 !  k for ww +ve,  k+1 for ww -ve
        rr=delu(iq,ii-kp)/(delu(iq,ii)+sign(1.E-20,delu(iq,ii)))
        fl=ww(iq,ii)*uu(iq,kx)
        cc = max(0.,min(2.*rr, 0.5+0.5*rr,2.)) ! MC
        fh=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1))          &
          -0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq) &
          /max(depdum(iq,ii+1)-depdum(iq,ii),1.E-10)
        if ( dzdum(iq,ii+1)>1.e-4 ) then 
          ff(iq,ii)=fl+cc*(fh-fl)
        end if 
      end do
      do ii=1,wlev
        if ( dzdum(iq,ii)>1.e-4 ) then  
          uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1)) &
                                        -ff(iq,ii)+ff(iq,ii-1))/dzdum(iq,ii)
        end if  
      end do
    end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

#ifdef GPU
  !$omp target exit data map(from:uu)
#endif
  !$acc update self(uu) async(async_counter)
  !$acc exit data delete(uu,delu,ff) async(async_counter)

else if ( mlontvd==1 ) then ! Superbee

#ifdef GPU
  !$omp target enter data map(to:uu) map(alloc:delu,ff)
#endif
  !$acc enter data create(uu,delu,ff) async(async_counter)
  !$acc update device(uu) async(async_counter)

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq)
#endif
#else
  !$acc parallel loop collapse(2) present(delu,uu,ff,dzdum) async(async_counter)
#endif
  do ii = 1,wlev-1
    do iq = 1,ifull
      ff(iq,ii) = 0.  
      if ( dzdum(iq,ii+1)>1.e-4 ) then
        delu(iq,ii) = uu(iq,ii+1) - uu(iq,ii)
      else
        delu(iq,ii) = 0.  
      end if    
    end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
  !$omp target teams distribute parallel do schedule(static) private(iq)
#endif
#else
  !$acc end parallel loop
  !$acc parallel loop present(ff,delu) async(async_counter)
#endif
  do iq = 1,ifull
    ff(iq,0) = 0.
    ff(iq,wlev) = 0.
    delu(iq,0) = 0.
    delu(iq,wlev) = 0.
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

! TVD part
#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq,kp,kx,rr,fl,cc,fh)
#endif
#else
  !$acc parallel loop collapse(2) present(ff,ww,delu,uu,dtnew,depdum,dzdum) async(async_counter)
#endif
  do ii = 1,wlev-1
    do iq = 1,ifull
      ! +ve ww is downwards to the ocean floor
      kp = nint(sign(1.,ww(iq,ii)))
      kx = ii+(1-kp)/2 !  k for ww +ve,  k+1 for ww -ve
      rr = delu(iq,ii-kp)/(delu(iq,ii)+sign(1.E-20,delu(iq,ii)))
      fl = ww(iq,ii)*uu(iq,kx)
      cc = max(0.,min(1.,2.*rr),min(2.,rr)) ! superbee
      fh = ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1))             &
        - 0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq) &
        /max(depdum(iq,ii+1)-depdum(iq,ii),1.E-10)
      if ( dzdum(iq,ii+1)>1.e-4 ) then
        ff(iq,ii) = fl + cc*(fh-fl)
      end if  
     !ff(iq,ii)=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1)) ! explicit
   end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
  !$omp target teams distribute parallel do collapse(2) schedule(static) private(ii,iq)
#endif
#else
  !$acc end parallel loop
  !$acc parallel loop collapse(2) present(ff,uu,ww,dtnew,dzdum) async(async_counter)
#endif
  do ii = 1,wlev
    do iq = 1,ifull
      if ( dzdum(iq,ii)>1.e-4 ) then  
        uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1))   &
                           -ff(iq,ii)+ff(iq,ii-1))/dzdum(iq,ii)
      end if   
    end do  
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

#ifdef _OPENMP
#ifdef GPU
  !$omp target teams distribute parallel do schedule(static) private(iq,i,ii,kp,kx,rr,fl,cc,fh)
#endif
#else
  !$acc parallel loop present(its,delu,uu,ww,ff,dtnew,depdum,dzdum) async(async_counter)
#endif
  do iq = 1,ifull
    do i = 2,its(iq)
      do ii = 1,wlev-1
        if ( dzdum(iq,ii+1)>1.e-4 ) then
          delu(iq,ii) = uu(iq,ii+1) - uu(iq,ii)
        else
          delu(iq,ii) = 0.  
        end if 
      end do
      ! TVD part
      do ii = 1,wlev-1
        ! +ve ww is downwards to the ocean floor
        kp = nint(sign(1.,ww(iq,ii)))
        kx = ii+(1-kp)/2 !  k for ww +ve,  k+1 for ww -ve
        rr=delu(iq,ii-kp)/(delu(iq,ii)+sign(1.E-20,delu(iq,ii)))
        fl=ww(iq,ii)*uu(iq,kx)
        cc=max(0.,min(1.,2.*rr),min(2.,rr)) ! superbee
        fh=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1))          &
          -0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq) &
          /max(depdum(iq,ii+1)-depdum(iq,ii),1.E-10)
        if ( dzdum(iq,ii+1)>1.e-4 ) then 
          ff(iq,ii)=fl+cc*(fh-fl)
        end if  
      end do
      do ii = 1,wlev
        if ( dzdum(iq,ii)>1.e-4 ) then  
          uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1)) &
                                        -ff(iq,ii)+ff(iq,ii-1))/dzdum(iq,ii)
        end if  
      end do
    end do
  end do
#ifdef _OPENMP
#ifdef GPU
  !$omp end target teams distribute parallel do
#endif
#else
  !$acc end parallel loop
#endif

#ifdef GPU
  !$omp target exit data map(from:uu)
#endif
  !$acc update self(uu) async(async_counter)
  !$acc exit data delete(uu,delu,ff) async(async_counter)
    
else
  write(6,*) "ERROR: Unknown option mlontvd ",mlontvd
  stop
endif

return
end subroutine mlotvd

end module mlovadvtvd
