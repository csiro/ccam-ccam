! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module mlodiffg

implicit none

private
public mlodiffusion,mlodiffusion_work
public mlodiff,ocndelphi,ocnsmag

integer, parameter :: nf          = 2       ! power for horizontal diffusion reduction factor
integer, save      :: mlodiff     = 0       ! diffusion (0=all, 1=scalars only)
real, save      :: ocndelphi      = 150.    ! horizontal diffusion reduction factor gradient (1.e6 = disabled)
real, save      :: ocnsmag        = 1.      ! horizontal diffusion (2. in Griffies (2000), 1.-1.4 in POM (Mellor 2004), 1. in SHOC)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion

use mlo
use newmpar_m

implicit none

real, dimension(ifull,wlev) :: u,v,tt,ss

! extract data from MLO
u=0.
v=0.
tt=0.
ss=0.
call mloexport3d(0,tt,0)
call mloexport3d(1,ss,0)
call mloexport3d(2,u,0)
call mloexport3d(3,v,0)

call mlodiffusion_work(u,v,tt,ss)

call mloimport3d(0,tt,0)
call mloimport3d(1,ss,0)
call mloimport3d(2,u,0)
call mloimport3d(3,v,0)

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion_work(u,v,tt,ss)

use cc_mpi
use const_phys
use indices_m
use map_m
use mlo
use mlodynamicsarrays_m
use nharrs_m, only : lrestart
use newmpar_m
use parm_m
use soil_m
use vecsuv_m

implicit none

integer iq, k
real hdif
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull+iextra,wlev) :: work_tt, work_ss
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact,dep
real, dimension(ifull+iextra,wlev) :: w_ema
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull,wlev), intent(inout) :: u,v,tt,ss
real, dimension(ifull,wlev) :: workdata2
real, dimension(ifull) :: dwdx, dwdy
real base
real dudx,dvdx,dudy,dvdy
real nu,nv,nw
real, dimension(ifull) :: emi
real tx_fact, ty_fact

call START_LOG(waterdiff_begin)

if ( abs(nmlo)>=3 ) then
  do k = 1,wlev  
    uau(1:ifull,k) = (av_vmod*u(1:ifull,k)+(1.-av_vmod)*oldu1(1:ifull,k))*ee(1:ifull,k)
    uav(1:ifull,k) = (av_vmod*v(1:ifull,k)+(1.-av_vmod)*oldv1(1:ifull,k))*ee(1:ifull,k)
  end do  
else
  do k = 1,wlev  
    uau(1:ifull,k) = u(1:ifull,k)*ee(1:ifull,k)
    uav(1:ifull,k) = v(1:ifull,k)*ee(1:ifull,k)
  end do
end if

! Define diffusion scale and grid spacing
hdif = dt*(ocnsmag/pi)**2
emi = dd(1:ifull)/em(1:ifull)

! calculate shear from EMA
call mlo_ema(dt,"uvw")
w_ema = 0.
do k = 1,wlev
  call mloexport("w_ema",w_ema(:,k),k,0)
end do
call bounds(w_ema)
do k = 2,wlev-1
  do iq = 1,ifull  
    dwdx(iq) = 0.5*((w_ema(ie(iq),k)-w_ema(iq,k))*emu(iq)*eeu(iq,k) &
                   +(w_ema(iq,k)-w_ema(iw(iq),k))*emu(iwu(iq))*eeu(iwu(iq),k))/ds
    dwdy(iq) = 0.5*((w_ema(in(iq),k)-w_ema(iq,k))*emv(iq)*eev(iq,k) &
                   +(w_ema(iq,k)-w_ema(is(iq),k))*emv(isv(iq))*eev(isv(iq),k))/ds
  end do  
  call mloimport("dwdx",dwdx,k,0)  
  call mloimport("dwdy",dwdy,k,0)  
end do

! calculate diffusion following Smagorinsky
call boundsuv(uau,uav,allvec=.true.)
do k = 1,wlev
  do iq = 1,ifull
    dudx = 0.5*((uau(ieu(iq),k)-uau(iq,k))*emu(iq)        &
               +(uau(iq,k)-uau(iwu(iq),k))*emu(iwu(iq)))/ds
    dudy = 0.5*((uau(inu(iq),k)-uau(iq,k))*emv(iq)        &
               +(uau(iq,k)-uau(isu(iq),k))*emv(isv(iq)))/ds
    dvdx = 0.5*((uav(iev(iq),k)-uav(iq,k))*emu(iq)        &
               +(uav(iq,k)-uav(iwv(iq),k))*emu(iwu(iq)))/ds
    dvdy = 0.5*((uav(inv(iq),k)-uav(iq,k))*emv(iq)        &
               +(uav(iq,k)-uav(isv(iq),k))*emv(isv(iq)))/ds
    t_kh(iq,k) = sqrt(dudx**2+dvdy**2+0.5*(dudy+dvdx)**2)*hdif*emi(iq)
  end do
end do
call bounds(t_kh(:,1:wlev),nehalf=.true.)


! reduce diffusion errors where bathymetry gradients are steep
if ( mlosigma>=0 .and. mlosigma<=3 ) then
  ! sigma levels  
  do k = 1,wlev
    dep(1:ifull,k) = gosig(1:ifull,k)*dd(1:ifull) !+ gosig(1:ifull,k)*eta(1:ifull)
  end do  
  call bounds(dep,nehalf=.true.)
  do k = 1,wlev
    do iq = 1,ifull
      tx_fact = 1./(1.+(abs(dep(ie(iq),k)-dep(iq,k))/ocndelphi)**nf)
      ty_fact = 1./(1.+(abs(dep(in(iq),k)-dep(iq,k))/ocndelphi)**nf)
      xfact(iq,k) = 0.5*(t_kh(iq,k)+t_kh(ie(iq),k))*tx_fact*eeu(iq,k) ! reduction factor
      yfact(iq,k) = 0.5*(t_kh(iq,k)+t_kh(in(iq),k))*ty_fact*eev(iq,k) ! reduction factor
    end do
  end do
else
  ! z* levels
  do k = 1,wlev
    xfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh(ie,k))*eeu(1:ifull,k)
    yfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh(in,k))*eev(1:ifull,k)
  end do
end if
call boundsuv(xfact,yfact,stag=-9)


! pre-process boundaries

if ( mlodiff==0 ) then
  ! Laplacian diffusion terms (closure #1)
  do k = 1,wlev
    duma(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(duma(:,:,1:3))
else if ( mlodiff==1 ) then
  ! no diffusion applied to momentum
else
  write(6,*) "ERROR: Unknown option for mlodiff = ",mlodiff
  call ccmpi_abort(-1)
end if

! Potential temperature and salinity
! MJT notes - only apply salinity diffusion to salt water
work_tt(1:ifull,:) = tt(1:ifull,:)
workdata2(1:ifull,:) = ss(1:ifull,:)
work_ss(1:ifull,:) = ss(1:ifull,:) - 34.72
where( workdata2(1:ifull,:)<2. )
  work_ss(1:ifull,:) = 0.
end where
call bounds(work_tt)
call bounds(work_ss)


! perform diffusion

#ifdef GPU
!$omp target data map(to:xfact,yfact,emi,iwu,isv,in,is,ie,iw)
#else
!$omp parallel
!$omp sections
#endif
!$acc data create(xfact,yfact,emi,iwu,isv,in,is,ie,iw)
!$acc update device(xfact,yfact,emi,iwu,isv,in,is,ie,iw)

#ifndef GPU    
!$omp section
#endif
if ( mlodiff==0 ) then
  ! UX  
  call mlodiffcalc(duma(:,:,1),xfact,yfact,emi)
endif

#ifndef GPU    
!$omp section
#endif
if ( mlodiff==0 ) then
  ! VY  
  call mlodiffcalc(duma(:,:,2),xfact,yfact,emi)
endif

#ifndef GPU    
!$omp section
#endif
if ( mlodiff==0 ) then
  ! WZ  
  call mlodiffcalc(duma(:,:,3),xfact,yfact,emi)
endif

#ifndef GPU    
!$omp section
#endif
! potential temperature
call mlodiffcalc(work_tt,xfact,yfact,emi)

#ifndef GPU    
!$omp section
#endif
! salinity
call mlodiffcalc(work_ss,xfact,yfact,emi)

#ifdef GPU
!$omp end target data
#else
!$omp end sections
!$omp end parallel
#endif
!$acc wait
!$acc end data


! post-processing

if ( mlodiff==0 ) then
  do k = 1,wlev
    do iq = 1,ifull
      u(iq,k) = ax(iq)*duma(iq,k,1) + ay(iq)*duma(iq,k,2) + az(iq)*duma(iq,k,3)
      v(iq,k) = bx(iq)*duma(iq,k,1) + by(iq)*duma(iq,k,2) + bz(iq)*duma(iq,k,3)
    end do
  end do
end if
do k = 1,wlev
  do iq = 1,ifull
    tt(iq,k) = max(work_tt(iq,k), -wrtemp)
    ss(iq,k) = work_ss(iq,k) + 34.72
    if ( workdata2(iq,k)<2. ) then
      ss(iq,k) = workdata2(iq,k)
    end if  
  end do
end do

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion_work

subroutine mlodiffcalc(work,xfact,yfact,emi)    

use indices_m
use mlo
use newmpar_m

implicit none

integer k, iq
integer, parameter :: async_length = 2
integer, save :: async_counter = -1
real, dimension(ifull+iextra,wlev), intent(in) :: xfact, yfact
real, dimension(ifull), intent(in) :: emi
real, dimension(ifull+iextra,wlev), intent(inout) :: work
real, dimension(ifull+iextra,wlev) :: ans
real base

async_counter = mod(async_counter+1, async_length)

#ifdef _OPENMP
#ifdef GPU
!$omp target enter data map(to:work) map(alloc:ans)
!$omp target teams distribute parallel do collapse(2) schedule(static) private(k,iq,base)
#endif
#else
!$acc enter data create(work,ans) async(async_counter)
!$acc update device(work) async(async_counter)
!$acc parallel loop collapse(2) present(work,ans,xfact,yfact,emi,iwu,isv,in,is,ie,iw) async(async_counter)
#endif
do k = 1,wlev
   do iq = 1,ifull
     base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                   +yfact(iq,k)+yfact(isv(iq),k)  
     ans(iq,k) = ( emi(iq)*work(iq,k) +               &
                   xfact(iq,k)*work(ie(iq),k) +       &
                   xfact(iwu(iq),k)*work(iw(iq),k) +  &
                   yfact(iq,k)*work(in(iq),k) +       &
                   yfact(isv(iq),k)*work(is(iq),k) )  &
                / base
  end do
end do
#ifdef _OPENMP
#ifdef GPU
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do collapse(2) schedule(static) private(k,iq)
#endif
#else
!$acc end parallel loop
!$acc parallel loop collapse(2) present(work,ans) async(async_counter)
#endif
do k = 1,wlev
   do iq = 1,ifull
     work(iq,k) = ans(iq,k)
  end do
end do
#ifdef _OPENMP
#ifdef GPU
!$omp end target teams distribute parallel do
!$omp target exit data map(from:work)
#endif
#else
!$acc end parallel loop
!$acc update self(work) async(async_counter)
!$acc exit data delete(work,ans) async(async_counter)
#endif

return
end subroutine mlodiffcalc

end module mlodiffg
