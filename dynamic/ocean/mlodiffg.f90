! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public mlodiff,ocnsmag,mlodiff_numits

integer, save   :: mlodiff        = 0       ! diffusion (0=all laplican, 1=scalars laplcian, 2=momentum laplacian,
                                            !            10=all biharmonic, 11=scalars biharmonic, 12=momentum biharmonic )
integer, save   :: mlodiff_numits = 10      ! Number of iterations required for bi-harmonic diffusion
real, save      :: ocnsmag        = 1.      ! horizontal diffusion for Bilinear (2. in Griffies (2000))

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion

use cc_mpi, only : ccmpi_abort
use mlo_ctrl
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
use mlo_ctrl
use mlodynamicsarrays_m
use nharrs_m, only : lrestart
use newmpar_m
use parm_m
use soil_m
use vecsuv_m

implicit none

integer iq, k
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull+iextra,wlev) :: ttl, ssl
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact,dep
real, dimension(ifull+iextra,wlev) :: w_ema
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull,wlev), intent(inout) :: u,v,tt,ss
real, dimension(ifull,wlev) :: workdata2
real, dimension(ifull) :: dwdx, dwdy
real hdif, base
real dudx,dvdx,dudy,dvdy
real nu,nv,nw
real, dimension(ifull) :: emi
real tx_fact, ty_fact

call START_LOG(waterdiff_begin)

duma = 0.

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
if ( mlodiff>=0 .and. mlodiff<=9 ) then
  ! Laplacian only (redefine ocnsmag)
  hdif = dt*(ocnsmag/pi)**2
else if ( mlodiff>=10 .and. mlodiff<=19 ) then
  ! Biharmonic and Laplacian
  hdif = sqrt(0.125)*ocnsmag/pi
else
  write(6,*) "ERROR: Unknown option mlodiff = ",mlodiff
  call ccmpi_abort(-1)
end if

emi(1:ifull) = 1./em(1:ifull)

! calculate shear for EMA
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
    dudx = 0.5*((uau(ieu(iq),k)-uau(iq,k))*emu(iq)*eeu(iq,k)                &
               +(uau(iq,k)-uau(iwu(iq),k))*emu(iwu(iq))*eeu(iwu(iq),k))/ds
    dudy = 0.5*((uau(inu(iq),k)-uau(iq,k))*emv(iq)*eev(iq,k)                &
               +(uau(iq,k)-uau(isu(iq),k))*emv(isv(iq))*eev(isv(iq),k))/ds
    dvdx = 0.5*((uav(iev(iq),k)-uav(iq,k))*emu(iq)*eeu(iq,k)                &
               +(uav(iq,k)-uav(iwv(iq),k))*emu(iwu(iq))*eeu(iwu(iq),k))/ds
    dvdy = 0.5*((uav(inv(iq),k)-uav(iq,k))*emv(iq)*eev(iq,k)                &
               +(uav(iq,k)-uav(isv(iq),k))*emv(isv(iq))*eev(isv(iq),k))/ds
    t_kh(iq,k) = sqrt(dudx**2+dvdy**2+0.5*(dudy+dvdx)**2)*emi(iq)*ee(iq,k)
  end do
end do
call bounds(t_kh(:,1:wlev),nehalf=.true.)


! reduce diffusion errors where bathymetry gradients are steep
if ( mlosigma>=0 .and. mlosigma<=3 ) then
  write(6,*) "ERROR: Unsupported option for mlosigma = ",mlodiff
  call ccmpi_abort(-1)
else
  ! z* levels
  do k = 1,wlev
    xfact(1:ifull,k) = sqrt(0.5*(t_kh(1:ifull,k)+t_kh(ie,k)))*eeu(1:ifull,k)
    yfact(1:ifull,k) = sqrt(0.5*(t_kh(1:ifull,k)+t_kh(in,k)))*eev(1:ifull,k)
  end do
end if
call boundsuv(xfact,yfact,stag=-9)


! perform diffusion
if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
  ! UX, VX, WX
  ! Laplacian diffusion terms (closure #1)
  do k = 1,wlev
    do iq = 1,ifull
      duma(iq,k,1) = ( ax(iq)*u(iq,k) + bx(iq)*v(iq,k) )*ee(iq,k)
      duma(iq,k,2) = ( ay(iq)*u(iq,k) + by(iq)*v(iq,k) )*ee(iq,k)
      duma(iq,k,3) = ( az(iq)*u(iq,k) + bz(iq)*v(iq,k) )*ee(iq,k)
    end do
  end do
  call bounds(duma(:,:,1:3))
end if
if ( mlodiff==0 .or. mlodiff==1 .or. mlodiff==10 .or. mlodiff==11 ) then
  ! potential temperature and salinity
  ! MJT notes - only apply salinity diffusion to salt water
  ttl(1:ifull,:) = tt(1:ifull,:)
  workdata2(1:ifull,:) = ss(1:ifull,:)
  ssl(1:ifull,:) = ss(1:ifull,:) - 34.72
  where( workdata2(1:ifull,:)<2. )
    ssl(1:ifull,:) = 0.
  end where  
  call bounds(ttl)
  call bounds(ssl)
end if


if ( mlodiff>=10 .and. mlodiff<20 ) then
  !Biharmonic version

  if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
    call mlodiffcalc(duma(:,:,1),xfact,yfact,emi,ee,hdif)
    call mlodiffcalc(duma(:,:,2),xfact,yfact,emi,ee,hdif)
    call mlodiffcalc(duma(:,:,3),xfact,yfact,emi,ee,hdif)
  end if
  if ( mlodiff==0 .or. mlodiff==1 .or. mlodiff==10 .or. mlodiff==11 ) then
    call mlodiffcalc(ttl,xfact,yfact,emi,ee,hdif)
    call mlodiffcalc(ssl,xfact,yfact,emi,ee,hdif)
  end if

else if ( mlodiff>=0 .and. mlodiff<10 ) then
  ! Laplacian version
    
  !$omp parallel
  !$omp sections

  !$omp section
  if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
    call mlodifflap(duma(:,:,1),xfact,yfact,emi,ee,hdif)
  end if
  !$omp section
  if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
    call mlodifflap(duma(:,:,2),xfact,yfact,emi,ee,hdif)
  end if
  !$omp section
  if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
    call mlodifflap(duma(:,:,3),xfact,yfact,emi,ee,hdif)
  end if
  !$omp section
  if ( mlodiff==0 .or. mlodiff==1 .or. mlodiff==10 .or. mlodiff==11 ) then
    call mlodifflap(ttl,xfact,yfact,emi,ee,hdif)
  end if
  !$omp section
  if ( mlodiff==0 .or. mlodiff==1 .or. mlodiff==10 .or. mlodiff==11 ) then
    call mlodifflap(ssl,xfact,yfact,emi,ee,hdif)
  end if

  !$omp end sections
  !$omp end parallel

else  
  write(6,*) "ERROR: Unknown mlodiff option mlodiff=",mlodiff
  call ccmpi_abort(-1)
end if  


if ( mlodiff==0 .or. mlodiff==2 .or. mlodiff==10 .or. mlodiff==12 ) then
  do k = 1,wlev
    do iq = 1,ifull
      u(iq,k) = ax(iq)*duma(iq,k,1) + ay(iq)*duma(iq,k,2) + az(iq)*duma(iq,k,3)
      v(iq,k) = bx(iq)*duma(iq,k,1) + by(iq)*duma(iq,k,2) + bz(iq)*duma(iq,k,3)
    end do
  end do
end if
if ( mlodiff==0 .or. mlodiff==1 .or. mlodiff==10 .or. mlodiff==11 ) then
  do k = 1,wlev
    do iq = 1,ifull
      tt(iq,k) = max(ttl(iq,k), -wrtemp)
      ss(iq,k) = ssl(iq,k) + 34.72
      if ( workdata2(iq,k)<2. ) then
        ss(iq,k) = workdata2(iq,k)
      end if  
    end do
  end do
end if  

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion_work

!------------------------------------------------------------------------------
subroutine mlodiffcalc(work,xfact,yfact,emi,ee,hdif)

use cc_mpi, only : bounds, ccmpi_abort
use indices_m
use mlo_ctrl
use newmpar_m
use parm_m, only : dt

implicit none

integer k, iq, its
real, dimension(ifull+iextra,wlev), intent(in) :: xfact, yfact
real, dimension(ifull), intent(in) :: emi
real, dimension(:,:), intent(inout) :: work
real, dimension(:,:), intent(in) :: ee
real, dimension(ifull+iextra,wlev) :: ans
real, dimension(ifull) :: ansl
real, dimension(ifull,wlev) :: work_save
real, intent(in) :: hdif
real base, xfact_iwu, yfact_isv

! Bi-Harmonic diffusion.  iterative version.

ans = 0. 
work_save(1:ifull,1:wlev) = work(1:ifull,1:wlev)

do its = 1,mlodiff_numits
  
  if ( its>1 ) then
    call bounds(work)
  end if

  ! Estimate Laplacian at t+1
  !$omp parallel do schedule(static) private(k,iq,xfact_iwu,yfact_isv,base)
  do k = 1,wlev
    do iq = 1,ifull  
      xfact_iwu = xfact(iwu(iq),k)
      yfact_isv = yfact(isv(iq),k)
      base = hdif*xfact(iq,k) + hdif*xfact_iwu &
           + hdif*yfact(iq,k) + hdif*yfact_isv
      ans(iq,k) = ( -base*work(iq,k) +                    &
                     hdif*xfact(iq,k)*work(ie(iq),k) +    &
                     hdif*xfact_iwu*work(iw(iq),k) +      &
                     hdif*yfact(iq,k)*work(in(iq),k) +    &
                     hdif*yfact_isv*work(is(iq),k) )      &
                  / sqrt(emi(iq))
    end do   
  end do
  !$omp end parallel do

  call bounds(ans)

  ! Estimate Laplacian^2 (= Grad^4) at t+1
  !$omp parallel do schedule(static) private(k,iq,xfact_iwu,yfact_isv,base)
  do k = 1,wlev
    do iq = 1,ifull  
      xfact_iwu = xfact(iwu(iq),k)
      yfact_isv = yfact(isv(iq),k)
      base = hdif*xfact(iq,k) + hdif*xfact_iwu &
           + hdif*yfact(iq,k) + hdif*yfact_isv
      work(iq,k) = work_save(iq,k) - dt*(                   &
                    -base*ans(iq,k) +                       &
                    hdif*xfact(iq,k)*ans(ie(iq),k) +        &
                    hdif*xfact_iwu*ans(iw(iq),k) +          &
                    hdif*yfact(iq,k)*ans(in(iq),k) +        &
                    hdif*yfact_isv*ans(is(iq),k) ) / sqrt(emi(iq))
    end do    
    do iq = 1,ifull
      if ( ee(iq,k)>0.5 ) then
        work(iq,k) = work_save(iq,k)
      end if
    end do
  end do  
  !$omp end parallel do
    
end do ! its

return
end subroutine mlodiffcalc

subroutine mlodifflap(work,xfact,yfact,emi,ee,ldif)

use indices_m
use mlo_ctrl
use newmpar_m

implicit none

integer k, iq
real, dimension(ifull+iextra,wlev), intent(in) :: xfact, yfact
real, dimension(ifull), intent(in) :: emi
real, dimension(:,:), intent(inout) :: work
real, dimension(:,:), intent(in) :: ee
real, dimension(ifull) :: ansl
real, intent(in) :: ldif
real base, xfact_iwu, yfact_isv

! Laplacian diffusion

do k = 1,wlev
  do iq = 1,ifull
    xfact_iwu = xfact(iwu(iq),k)
    yfact_isv = yfact(isv(iq),k)
    base = emi(iq)+ldif*xfact(iq,k)**2+ldif*xfact_iwu**2  &
                  +ldif*yfact(iq,k)**2+ldif*yfact_isv**2
    ansl(iq) = ( emi(iq)*work(iq,k) +                     &
                 ldif*xfact(iq,k)**2*work(ie(iq),k) +     &
                 ldif*xfact_iwu**2*work(iw(iq),k) +       &
                 ldif*yfact(iq,k)**2*work(in(iq),k) +     &
                 ldif*yfact_isv**2*work(is(iq),k) )       &
               / base
  end do
  do iq = 1,ifull
    if ( ee(iq,k)>0.5 ) then
      work(iq,k) = ansl(iq)
    end if
  end do
end do  
    
return
end subroutine mlodifflap

end module mlodiffg