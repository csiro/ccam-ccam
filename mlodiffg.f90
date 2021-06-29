! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
use newmpar_m
use parm_m
use soil_m
use vecsuv_m

implicit none

integer iq, k
real hdif
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact,dep
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull,wlev), intent(inout) :: u,v,tt,ss
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
call boundsuv(uau,uav,allvec=.true.)

! Define diffusion scale and grid spacing
hdif = dt*(ocnsmag/pi)**2
emi = dd(1:ifull)/em(1:ifull)

! calculate diffusion following Smagorinsky
!$acc parallel loop collapse(2) copyin(emu,emv,emi,uau,uav) copyout(t_kh(1:ifull,1:wlev)) &
!$acc   present(ieu,iwu,inu,isu,iev,iwv,inv,isv)
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
!$acc end parallel loop
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


if ( mlodiff==0 ) then
  ! Laplacian diffusion terms (closure #1)
  do k = 1,wlev
    duma(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(duma(:,:,1:3))
  do k = 1,wlev
    do iq = 1,ifull
      base = emi(iq) + xfact(iq,k) + xfact(iwu(iq),k) &
                     + yfact(iq,k) + yfact(isv(iq),k)
      nu = ( duma(iq,k,1)*emi(iq) +                       &
             xfact(iq,k)*duma(ie(iq),k,1) +               &
             xfact(iwu(iq),k)*duma(iw(iq),k,1) +          &
             yfact(iq,k)*duma(in(iq),k,1) +               &
             yfact(isv(iq),k)*duma(is(iq),k,1) )          &
           / base
      nv = ( duma(iq,k,2)*emi(iq) +                       &
             xfact(iq,k)*duma(ie(iq),k,2) +               &
             xfact(iwu(iq),k)*duma(iw(iq),k,2) +          &
             yfact(iq,k)*duma(in(iq),k,2) +               &
             yfact(isv(iq),k)*duma(is(iq),k,2) )          &
           / base
      nw = ( duma(iq,k,3)*emi(iq) +                       &
             xfact(iq,k)*duma(ie(iq),k,3) +               &
             xfact(iwu(iq),k)*duma(iw(iq),k,3) +          &
             yfact(iq,k)*duma(in(iq),k,3) +               &
             yfact(isv(iq),k)*duma(is(iq),k,3) )          &
           / base
      u(iq,k) = ax(iq)*nu + ay(iq)*nv + az(iq)*nw
      v(iq,k) = bx(iq)*nu + by(iq)*nv + bz(iq)*nw
    end do
  end do

else if ( mlodiff==1 ) then
  ! no diffusion applied to momentum
else
  write(6,*) "ERROR: Unknown option for mlodiff = ",mlodiff
  call ccmpi_abort(-1)
end if

! Potential temperature and salinity
duma(1:ifull,:,1) = tt(1:ifull,:)
duma(1:ifull,:,2) = ss(1:ifull,:) - 34.72
call bounds(duma(:,:,1:2))
!$acc parallel loop collapse(2) copyin(emi,xfact,yfact,duma(:,:,1:2)) copyout(tt,ss)
do k = 1,wlev
  do iq = 1,ifull
    base = emi(iq) + xfact(iq,k) + xfact(iwu(iq),k) &
                   + yfact(iq,k) + yfact(isv(iq),k)
    tt(iq,k) = ( duma(iq,k,1)*emi(iq) +               &
                 xfact(iq,k)*duma(ie(iq),k,1) +       &
                 xfact(iwu(iq),k)*duma(iw(iq),k,1) +  &
                 yfact(iq,k)*duma(in(iq),k,1) +       &
                 yfact(isv(iq),k)*duma(is(iq),k,1) )  &
               / base
    tt(iq,k) = max(tt(iq,k), -wrtemp)
    ss(iq,k) = ( duma(iq,k,2)*emi(iq) +               &
                 xfact(iq,k)*duma(ie(iq),k,2) +       &
                 xfact(iwu(iq),k)*duma(iw(iq),k,2) +  &
                 yfact(iq,k)*duma(in(iq),k,2) +       &
                 yfact(isv(iq),k)*duma(is(iq),k,2) )  &
               / base
    ss(iq,k) = max(ss(iq,k)+34.72, 0.)
  end do
end do
!$acc end parallel loop

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion_work

end module mlodiffg
