! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
integer, save      :: mlodiff     = 0       ! diffusion (0=all laplican, 1=scalars laplican, 10=all biharmonic, 11=scalars biharmonic )
real, save      :: ocndelphi      = 150.    ! horizontal diffusion reduction factor gradient (1.e6 = disabled)
real, save      :: ocnsmag        = 1.      ! horizontal diffusion (2. in Griffies (2000), 1.-1.4 in POM (Mellor 2004), 1. in SHOC)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion

use cc_mpi, only : ccmpi_abort
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
real, dimension(ifull+iextra,wlev,2) :: work
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
if ( mlodiff>=0 .and. mlodiff<=9 ) then
  hdif = dt*(ocnsmag/pi)**2
else if ( mlodiff>=10 .and. mlodiff<=19 ) then
  hdif = 0.125*(ocnsmag/pi)**2  
end if
if ( mlosigma>=0 .and. mlosigma<=3 ) then
  write(6,*) "ERROR: Unsupported option for mlosigma = ",mlodiff
  call ccmpi_abort(-1)
else
  ! z* levels  
  emi = 1./em(1:ifull)
end if

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
    t_kh(iq,k) = sqrt(dudx**2+dvdy**2+0.5*(dudy+dvdx)**2)*hdif*emi(iq)*ee(iq,k)
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
    xfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh(ie,k))*eeu(1:ifull,k)
    yfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh(in,k))*eev(1:ifull,k)
  end do
end if
if ( mlodiff>=10 .and. mlodiff<=19 ) then
  ! Take sqrt as we apply the Laplican twice for Grad^4  
  xfact(1:ifull,1:wlev) = sqrt(xfact(1:ifull,1:wlev))
  yfact(1:ifull,1:wlev) = sqrt(yfact(1:ifull,1:wlev))
end if
call boundsuv(xfact,yfact,stag=-9)

! pre-process boundaries

if ( mlodiff==0 .or. mlodiff==10 ) then
  ! Laplacian diffusion terms (closure #1)
  do k = 1,wlev
    duma(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
else if ( mlodiff==1 .or. mlodiff==11 ) then
  ! no diffusion applied to momentum
else
  write(6,*) "ERROR: Unknown option for mlodiff = ",mlodiff
  call ccmpi_abort(-1)
end if

! Potential temperature and salinity
! MJT notes - only apply salinity diffusion to salt water
work(1:ifull,:,1) = tt(1:ifull,:)
workdata2(1:ifull,:) = ss(1:ifull,:)
work(1:ifull,:,2) = ss(1:ifull,:) - 34.72
where( workdata2(1:ifull,:)<2. )
  work(1:ifull,:,2) = 0.
end where


! perform diffusion

#ifdef GPU
!$acc data create(xfact,yfact,emi)
!$acc update device(xfact,yfact,emi)
#endif

if ( mlodiff==0 .or. mlodiff==10 ) then
  ! UX, VX, WX
  call mlodiffcalc(duma(:,:,1:3),xfact,yfact,emi,3)
end if
! potential temperature and salinity
call mlodiffcalc(work(:,:,1:2),xfact,yfact,emi,2)

#ifdef GPU
!$acc end data
#endif


! post-processing

if ( mlodiff==0 .or. mlodiff==10 ) then
  do k = 1,wlev
    do iq = 1,ifull
      u(iq,k) = ax(iq)*duma(iq,k,1) + ay(iq)*duma(iq,k,2) + az(iq)*duma(iq,k,3)
      v(iq,k) = bx(iq)*duma(iq,k,1) + by(iq)*duma(iq,k,2) + bz(iq)*duma(iq,k,3)
    end do
  end do
end if
do k = 1,wlev
  do iq = 1,ifull
    tt(iq,k) = max(work(iq,k,1), -wrtemp)
    ss(iq,k) = work(iq,k,2) + 34.72
    if ( workdata2(iq,k)<2. ) then
      ss(iq,k) = workdata2(iq,k)
    end if  
  end do
end do

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion_work

subroutine mlodiffcalc(work,xfact,yfact,emi,ntr)    

use cc_acc, only : async_length
use cc_mpi, only : bounds_send, bounds_recv, maxcolour, iqx, &
                   ifull_colour, ifull_colour_border,        &
                   ccmpi_abort, nagg
use indices_m
use map_m
use mlo
use mlodynamicsarrays_m, only : ee
use newmpar_m
use parm_m, only : dt, ds

implicit none

integer, intent(in) :: ntr
integer k, iq, iqc, nc, its
integer nstart, nend, nlen, nn, np
integer, parameter :: num_its = 4
integer async_counter
real, dimension(ifull+iextra,wlev), intent(in) :: xfact, yfact
real, dimension(ifull), intent(in) :: emi
real, dimension(ifull+iextra,wlev,ntr), intent(inout) :: work
real, dimension(ifull+iextra,wlev,nagg) :: ans
real, dimension(ifull,wlev,nagg) :: work_save
real base, xfact_iwu, yfact_isv

if ( mlodiff>=0 .and. mlodiff<=9 ) then
    
  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr )
    nlen = nend - nstart + 1  
    
    call bounds_send(work(:,:,nstart:nend))

    ! update non-boundary grid points
    ! here we use the coloured indices to identify interior and boundary points  
    do nn = 1,nlen
      np = nn - 1 + nstart
#ifndef GPU
      !$omp parallel do schedule(static) private(k,iq,xfact_iwu,yfact_isv,base)
#else
      async_counter = mod(nn-1,async_length)
      !$acc parallel loop collapse(2) copyin(work(:,:,np)) copyout(ans(1:ifull,:,nn)) &
      !$acc   present(xfact,yfact,emi,in,is,ie,iw,iwu,isv)                            &  
      !$acc   private(xfact_iwu,yfact_isv,base) async(async_counter)
#endif
      do k = 1,wlev
        do iq = 1,ifull
          xfact_iwu = xfact(iwu(iq),k)
          yfact_isv = yfact(isv(iq),k)
          base = emi(iq)+xfact(iq,k)+xfact_iwu  &
                        +yfact(iq,k)+yfact_isv
          ans(iq,k,nn) = ( emi(iq)*work(iq,k,np) +               &
                           xfact(iq,k)*work(ie(iq),k,np) +       &
                           xfact_iwu*work(iw(iq),k,np) +         &
                           yfact(iq,k)*work(in(iq),k,np) +       &
                           yfact_isv*work(is(iq),k,np) )         &
                        / base
        end do  
      end do
#ifndef GPU
      !$omp end parallel do
#else
      !$acc end parallel loop
#endif
    end do
#ifdef GPU
    !$acc wait
#endif

    call bounds_recv(work(:,:,nstart:nend))

    ! update boundary grid points
    ! here we use the coloured indices to identify interior and boundary points  
    do nn = 1,nlen
      np = nn - 1 + nstart
      do k = 1,wlev
        do nc = 1,maxcolour  
          do iqc = 1,ifull_colour_border(nc)
            iq = iqx(iqc,nc)  
            xfact_iwu = xfact(iwu(iq),k)
            yfact_isv = yfact(isv(iq),k)
            base = emi(iq)+xfact(iq,k)+xfact_iwu  &
                          +yfact(iq,k)+yfact_isv
            ans(iq,k,nn) = ( emi(iq)*work(iq,k,np) +               &
                             xfact(iq,k)*work(ie(iq),k,np) +       &
                             xfact_iwu*work(iw(iq),k,np) +         &
                             yfact(iq,k)*work(in(iq),k,np) +       &
                             yfact_isv*work(is(iq),k,np) )         &
                          / base
          end do  
        end do
        do iq = 1,ifull
          work(iq,k,np) = ans(iq,k,nn)
        end do
      end do  
    end do

  end do ! nstart
    
else if ( mlodiff>=10 .and. mlodiff<=19 ) then

  ! Bi-Harmonic diffusion.  iterative version.

  ans = 0.
    
#ifdef GPU
  !$acc enter data create(work_save)   
#endif
  
  do nstart = 1,ntr,nagg
    nend = min(nstart + nagg - 1, ntr )
    nlen = nend - nstart + 1

    work_save(1:ifull,1:wlev,1:nlen) = work(1:ifull,1:wlev,nstart:nend)
#ifdef GPU
    !$acc update device(work_save)
#endif
    
    do its = 1,num_its  
      
      call bounds_send(work(:,:,nstart:nend))

      do nn = 1,nlen
        np = nn - 1 + nstart  
#ifndef GPU
        !$omp parallel do schedule(static) private(iq,k,xfact_iwu,yfact_isv,base)
#else
        async_counter = mod(nn-1,async_length)
        !$acc parallel loop collapse(2) copyin(work(:,:,np)) copyout(ans(1:ifull,:,nn)) &
        !$acc   present(xfact,yfact,emi,in,is,ie,iw,iwu,isv)                            &
        !$acc   private(xfact_iwu,yfact_isv,base) async(async_counter)
#endif
        do k = 1,wlev
          do iq = 1,ifull  
            xfact_iwu = xfact(iwu(iq),k)
            yfact_isv = yfact(isv(iq),k)
            base = xfact(iq,k) + xfact_iwu + yfact(iq,k) + yfact_isv
            ans(iq,k,nn) = ( -base*work(iq,k,np) +              &
                             xfact(iq,k)*work(ie(iq),k,np) +    &
                             xfact_iwu*work(iw(iq),k,np) +      &
                             yfact(iq,k)*work(in(iq),k,np) +    &
                             yfact_isv*work(is(iq),k,np) ) / sqrt(emi(iq))
          end do   
        end do
#ifndef GPU
        !$omp end parallel do
#else
        !$acc end parallel loop
#endif
      end do
#ifdef GPU
      !$acc wait
#endif

      call bounds_recv(work(:,:,nstart:nend))

      do nn = 1,nlen
        np = nn - 1 + nstart  
        do k = 1,wlev
          do nc = 1,maxcolour  
            do iqc = 1,ifull_colour_border(nc)
              iq = iqx(iqc,nc)
              xfact_iwu = xfact(iwu(iq),k)
              yfact_isv = yfact(isv(iq),k)
              base = xfact(iq,k) + xfact_iwu + yfact(iq,k) + yfact_isv
              ans(iq,k,nn) = ( -base*work(iq,k,np) +              &
                               xfact(iq,k)*work(ie(iq),k,np) +    &
                               xfact_iwu*work(iw(iq),k,np) +      &
                               yfact(iq,k)*work(in(iq),k,np) +    &
                               yfact_isv*work(is(iq),k,np) ) / sqrt(emi(iq))
            end do  
          end do
        end do  
      end do

      call bounds_send(ans(:,:,1:nlen))

      do nn = 1,nlen
        np = nn - 1 + nstart
#ifndef GPU
        !$omp parallel do schedule(static) private(iq,k,xfact_iwu,yfact_isv,base)
#else
        async_counter = mod(nn-1,async_length)  
        !$acc parallel loop collapse(2) copyin(ans(:,:,nn)) copyout(work(1:ifull,:,np))  &
        !$acc   present(work_save,xfact,yfact,emi,in,is,ie,iw,iwu,isv)                   &
        !$acc   private(xfact_iwu,yfact_isv,base) async(async_counter)
#endif        
        do k = 1,wlev
          do iq = 1,ifull  
            xfact_iwu = xfact(iwu(iq),k)
            yfact_isv = yfact(isv(iq),k)
            base = xfact(iq,k) + xfact_iwu + yfact(iq,k) + yfact_isv
            work(iq,k,np) = work_save(iq,k,nn) - dt*(               &
                              -base*ans(iq,k,nn) +                  &
                              xfact(iq,k)*ans(ie(iq),k,nn) +        &
                              xfact_iwu*ans(iw(iq),k,nn) +          &
                              yfact(iq,k)*ans(in(iq),k,nn) +        &
                              yfact_isv*ans(is(iq),k,nn) ) / sqrt(emi(iq))
          end do    
        end do
#ifndef GPU
        !$omp end parallel do
#else
        !$acc end parallel loop
#endif
      end do
#ifdef GPU
      !$acc wait
#endif

      call bounds_recv(ans(:,:,1:nlen))

      do nn = 1,nlen
        np = nn - 1 + nstart
        do k = 1,wlev
          do nc = 1,maxcolour  
            do iqc = 1,ifull_colour_border(nc)
              iq = iqx(iqc,nc)  
              xfact_iwu = xfact(iwu(iq),k)
              yfact_isv = yfact(isv(iq),k)
              base = xfact(iq,k) + xfact_iwu + yfact(iq,k) + yfact_isv
              work(iq,k,np) = work_save(iq,k,nn) - dt*(            &
                             -base*ans(iq,k,nn) +                  &
                             xfact(iq,k)*ans(ie(iq),k,nn) +        &
                             xfact_iwu*ans(iw(iq),k,nn) +          &
                             yfact(iq,k)*ans(in(iq),k,nn) +        &
                             yfact_isv*ans(is(iq),k,nn) ) / sqrt(emi(iq))
            end do  
          end do  
          do iq = 1,ifull
            if ( ee(iq,k)<0.5 ) then
              work(iq,k,np) = work_save(iq,k,nn)
            end if
          end do
        end do  
      end do
    
    end do ! its
    
  end do   ! nstart  
  
#ifdef GPU
  !$acc exit data delete(work_save)
#endif

else
  write(6,*) "ERROR: Unknown mlodiff option ",mlodiff
  call ccmpi_abort(-1)
end if

return
end subroutine mlodiffcalc

end module mlodiffg
