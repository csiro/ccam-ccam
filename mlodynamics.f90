! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! These subroutines handle dynamics for the Mixed-Layer-Ocean model

! - Horizontal diffusion
! - Ocean dynamics
! - Ice dynamics

! This module links to mlo.f90 which solves for 1D physics of ocean
! and ice processes.  Currently the code assumes the hydrostatic
! approximation which is reasonably valid to 1km resolution.

! Ocean and sea-ice dynamics are based on the R-grid used by CCAM.
! sigma-s coordinates are used by the ocean to improve coastal
! regions.  Flexible nudging options are used for error correction
! (see nesting.f90).
    
! Several versions of the pressure gradient terms are avaliable and
! are specified using mlojacobi.  mlo_rtest can also be used to
! check that the grid is sufficently smooth for sigma coordinates.

module mlodynamics

implicit none

private
public mlodiffusion,mlohadv,mlodyninit
public gosig,gosigh,godsig,ocnsmag,ocneps
public mlodiff,usetide,mlojacobi,mlo_rtest
public usepice
public dd
public nstagoffmlo,mstagf,koff

real, dimension(:), allocatable, save :: ee,eeu,eev,dd,ddu,ddv
real, dimension(:), allocatable, save :: gosig,gosigh,godsig
real, dimension(:,:), allocatable, save :: stwgt
integer, save :: nstagoffmlo    = 0       ! staggering offset
integer, save :: usetide        = 1       ! tidal forcing (0=off, 1=on)
integer, parameter :: icemode   = 2       ! ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, parameter :: mstagf    = 30      ! alternating staggering (0=off left, -1=off right, >0 alternating)
integer, parameter :: koff      = 1       ! time split stagger relative to A-grid (koff=0) or C-grid (koff=1)
integer, parameter :: nf        = 2       ! power for horizontal diffusion reduction factor
integer, parameter :: itnmax    = 6       ! number of interations for reversible staggering
integer, parameter :: nxtrrho   = 1       ! Estimate rho at t+1 (0=off, 1=on)
integer, save      :: usepice   = 0       ! include ice in surface pressure (0=without ice, 1=with ice)
integer, save      :: mlodiff   = 0       ! diffusion (0=all, 1=scalars only)
integer, save      :: mlojacobi = 1       ! density gradient method (0=off, 1=non-local spline, 2=non-local linear, 3=Song)
real, parameter :: rhosn      = 330.      ! density snow (kg m^-3)
real, parameter :: rhoic      = 900.      ! density ice  (kg m^-3)
real, parameter :: grav       = 9.80616   ! gravitational constant (m s^-2)
real, parameter :: delphi     = 150.      ! horizontal diffusion reduction factor gradient
real, save      :: ocnsmag    = 1.        ! horizontal diffusion (2. in Griffies (2000), 1.-1.4 in POM (Mellor 2004), 1. in SHOC)
real, save      :: ocneps     = 0.1       ! semi-implicit off-centring term
real, parameter :: maxicefrac = 0.999     ! maximum ice fraction
real, parameter :: tol        = 5.E-3     ! Tolerance for SOR solver (water)
real, parameter :: itol       = 1.E1      ! Tolerance for SOR solver (ice)
real, save      :: mlo_rtest  = 0.        ! check bathymetry (0.=off, 0.2=recommended)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises mlo dynamical arrays
!
subroutine mlodyninit

use cc_mpi
use indices_m
use map_m
use mlo, only : wlev,mloexpdep
use mlodynamicsarrays_m
use newmpar_m
use parm_m
use parmdyn_m
use soil_m

implicit none

integer ii,iq
real rmax_l, r0max_g, r1max_g
real, dimension(ifull,0:wlev) :: dephl
real, dimension(ifull,wlev) :: dep,dz
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn
real, dimension(wlev) :: sig,sigh,dsig
real, dimension(3*wlev) :: dumz,gdumz
logical, dimension(ifull+iextra) :: wtr

! prep land-sea mask
allocate( ee(ifull+iextra) )
allocate( eeu(ifull+iextra), eev(ifull+iextra) )
ee  = 0.
eeu = 0.
eev = 0. 
where ( .not.land )
  ee(1:ifull) = 1.
end where
call bounds(ee,nrows=2)
where ( ee>1.5 .or. ee<=0.5 )
  ee = 0.
end where
eeu(1:ifull) = ee(1:ifull)*ee(ie)
eev(1:ifull) = ee(1:ifull)*ee(in)
call boundsuv(eeu,eev,nrows=2)

wtr = abs(ee-1.)<0.

! The following is for the in-line MLO model ------------------------

! Calculate depth arrays (free suface term is included later)
allocate( dd(ifull+iextra) )
allocate( ddu(ifull+iextra), ddv(ifull+iextra) )
dd = 0.
ddu = 0.
ddv = 0.
dep = 0.
dz = 0.
do ii = 1,wlev
  call mloexpdep(0,dep(:,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
end do
dephl(:,0) = 0.
dephl(:,1) = dz(:,1)
do ii = 2,wlev
  dephl(:,ii) = dephl(:,ii-1) + dz(:,ii)
end do
dd(1:ifull) = dephl(:,wlev)
dd(1:ifull) = max(dd(1:ifull), 1.E-8)
call bounds(dd,nrows=2)
ddu(1:ifull) = 0.5*(dd(1:ifull)+dd(ie))
ddv(1:ifull) = 0.5*(dd(1:ifull)+dd(in))
where ( abs(eeu(1:ifull))<1.e-20 )
  ddu(1:ifull) = 1.E-8
end where
where ( abs(eev(1:ifull))<1.e-20 )
  ddv(1:ifull) = 1.E-8
end where
call boundsuv(ddu,ddv,nrows=2)

! allocate memory for mlo dynamics arrays
call mlodynamicsarrays_init(ifull,iextra,wlev)

if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then

  ! Precompute weights for calculating staggered gradients
  allocate(stwgt(ifull+iextra,4))
  stwgt=0.
  where (wtr(in).and.wtr(ien).and.wtr(ie).and.wtr(1:ifull))
    stwgt(1:ifull,1)=1.
  end where
  where (wtr(is).and.wtr(ies).and.wtr(ie).and.wtr(1:ifull))
    stwgt(1:ifull,2)=1.
  end where
  where (wtr(in).and.wtr(ine).and.wtr(ie).and.wtr(1:ifull))
    stwgt(1:ifull,3)=1.
  end where
  where (wtr(iw).and.wtr(inw).and.wtr(in).and.wtr(1:ifull))
    stwgt(1:ifull,4)=1.
  end where
  call boundsuv(stwgt(:,1:2),stwgt(:,3:4))
  
end if ! abs(nmlo)>=3.and.abs(nmlo)<=9

! sigma coordinates should be the same for all iq
allocate(gosig(wlev),gosigh(wlev),godsig(wlev))
sig=0.
sigh=0.
dsig=0.
do iq=1,ifull
  if (.not.land(iq)) then
    do ii=1,wlev
      sig(ii)=max(sig(ii),dep(iq,ii)/dd(iq))
      sigh(ii)=max(sigh(ii),dephl(iq,ii)/dd(iq))
      dsig(ii)=max(dsig(ii),dz(iq,ii)/dd(iq))
    end do
  end if
end do
do iq=1,ifull
  if (.not.land(iq)) then
    do ii=1,wlev
      if (abs(sig(ii)*dd(iq)-dep(iq,ii))>0.01    .or. &
          abs(sigh(ii)*dd(iq)-dephl(iq,ii))>0.01 .or. &
          abs(dsig(ii)*dd(iq)-dz(iq,ii))>0.01) then
        write(6,*) "ERROR: Invalid sigma levels for mlodynamics"
        call ccmpi_abort(-1)
      end if
    end do
  end if
end do
dumz(1:wlev)          = sig
dumz(wlev+1:2*wlev)   = sigh
dumz(2*wlev+1:3*wlev) = dsig
call ccmpi_allreduce(dumz(1:3*wlev),gdumz(1:3*wlev),"max",comm_world)
gosig  = gdumz(1:wlev)
gosigh = gdumz(wlev+1:2*wlev)
godsig = gdumz(2*wlev+1:3*wlev)


! r-test
rmax_l = 0.
do iq = 1,ifull
  if ( dd(iq)>0.1 ) then
    if ( dd(in(iq))>0.1 ) then  
      rmax_l = max( rmax_l, abs(dd(in(iq))-dd(iq))/(dd(in(iq))+dd(iq)) )
    end if
    if ( dd(ie(iq))>0.1 ) then
      rmax_l = max( rmax_l, abs(dd(ie(iq))-dd(iq))/(dd(ie(iq))+dd(iq)) )
    end if
  end if  
end do
call ccmpi_reduce(rmax_l,r0max_g,"max",0,comm_world)
if ( myid==0 ) then
  r1max_g = 0.  
  do ii = 1,wlev-1
    r1max_g = max( r1max_g, r0max_g*(gosig(ii)+gosig(ii+1))/abs(gosig(ii)-gosig(ii+1)) )
  end do  
  write(6,*) "MLODYNAMICS rtest: r0max_g, r1max_g = ",r0max_g,r1max_g
  if ( mlo_rtest>0. ) then
    if ( r0max_g>mlo_rtest*1.1 ) then
      write(6,*) "ERROR: mlodynamics rtest failed."
      write(6,*) "Please run ocnbath with bathfilt=t and rtest=",mlo_rtest
      call ccmpi_abort(-1)
    end if
  end if  
end if

return
end subroutine mlodyninit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion

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

integer k, iq
integer kp1, km1
real hdif
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull,wlev) :: u,v,tt,ss
real, dimension(ifull,wlev) :: fs,base,outu,outv
real, dimension(ifull,wlev) :: xfact_iwu, yfact_isv
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: nu,nv,nw
real, dimension(ifull) :: tx_fact,ty_fact
real, dimension(ifull) :: emi, emu_w, eeu_w, emv_s, eev_s
real, dimension(ifull) :: dd_e, dd_n, dd_w, dd_s
real, dimension(ifull) :: duma_n, duma_s, duma_e, duma_w
real, dimension(ifull) :: v_n, v_s, u_e, u_w
real, dimension(ifull) :: v_e, v_w, u_n, u_s
real, dimension(ifull) :: vm1_n, vm1_s, um1_e, um1_w
real, dimension(ifull) :: vm1_e, vm1_w, um1_n, um1_s
real, dimension(ifull) :: ddi, ddi_n, ddi_e, ddim1, ddim1_n, ddim1_e
real, dimension(ifull) :: ddi_s, ddi_w, ddim1_s, ddim1_w
real, dimension(ifull) :: t_kh_n, t_kh_e

call START_LOG(waterdiff_begin)

! extract data from MLO
u=0.
v=0.
tt=0.
ss=0.
call mloexport3d(0,tt,0)
call mloexport3d(1,ss,0)
call mloexport3d(2,u,0)
call mloexport3d(3,v,0)
if ( abs(nmlo)>=3 ) then
  do k = 1,wlev  
    uau(1:ifull,k) = (av_vmod*u(1:ifull,k)+(1.-av_vmod)*oldu1(1:ifull,k))*ee(1:ifull)
    uav(1:ifull,k) = (av_vmod*v(1:ifull,k)+(1.-av_vmod)*oldv1(1:ifull,k))*ee(1:ifull)
  end do  
else
  do k = 1,wlev  
    uau(1:ifull,k) = u(1:ifull,k)*ee(1:ifull)
    uav(1:ifull,k) = v(1:ifull,k)*ee(1:ifull)
  end do
end if
call boundsuv(uau,uav,allvec=.true.)

! Define diffusion scale and grid spacing
hdif = dt*(ocnsmag/pi)**2
emi = dd(1:ifull)/em(1:ifull)

! calculate diffusion following Smagorinsky
call unpack_svwu(emu,emv,emv_s,emu_w)
call unpack_svwu(eeu,eev,eev_s,eeu_w)
do k = 1,wlev
  call unpack_nveu(uau(:,k),uav(:,k),v_n,u_e)  
  call unpack_svwu(uau(:,k),uav(:,k),v_s,u_w)
  dudx = 0.5*((u_e-uau(1:ifull,k))*emu(1:ifull)*eeu(1:ifull)        &
             +(uau(1:ifull,k)-u_w)*emu_w*eeu_w)/ds
  dudy = 0.5*((uau(inu,k)-uau(1:ifull,k))*emv(1:ifull)*eev(1:ifull) &
             +(uau(1:ifull,k)-uau(isu,k))*emv_s*eev_s)/ds
  dvdx = 0.5*((uav(iev,k)-uav(1:ifull,k))*emu(1:ifull)*eeu(1:ifull) &
             +(uav(1:ifull,k)-uav(iwv,k))*emu_w*eeu_w)/ds
  dvdy = 0.5*((v_n-uav(1:ifull,k))*emv(1:ifull)*eev(1:ifull)        &
             +(uav(1:ifull,k)-v_s)*emv_s*eev_s)/ds

  !t_kh(1:ifull,k) = sqrt((dudx-dvdy)**2+(dudy+dvdx)**2)*hdif*emi
  t_kh(1:ifull,k) = sqrt(dudx**2+dvdy**2+0.5*(dudy+dvdx)**2)*hdif*emi
end do
!t_kh(1:ifull,wlev+1) = etain(1:ifull)
call bounds(t_kh(:,1:wlev),nehalf=.true.)
!eta(:) = t_kh(:,wlev+1)

! reduce diffusion errors where bathymetry gradients are strong
call unpack_ne(dd,dd_n,dd_e)
do k = 1,wlev
  !depadj = gosig(k)*max(dd+eta,minwater) ! neglect eta
  tx_fact = 1./(1.+(gosig(k)*abs(dd_e-dd(1:ifull))/delphi)**nf)
  ty_fact = 1./(1.+(gosig(k)*abs(dd_n-dd(1:ifull))/delphi)**nf)

  call unpack_ne(t_kh(:,k),t_kh_n,t_kh_e)
  xfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh_e)*tx_fact*eeu(1:ifull) ! reduction factor
  yfact(1:ifull,k) = 0.5*(t_kh(1:ifull,k)+t_kh_n)*ty_fact*eev(1:ifull) ! reduction factor
end do
call boundsuv(xfact,yfact,stag=-9)

do k = 1,wlev
  call unpack_svwu(xfact(:,k),yfact(:,k),yfact_isv(:,k),xfact_iwu(:,k))  
  base(:,k) = emi + xfact(1:ifull,k) + xfact_iwu(:,k) + yfact(1:ifull,k) + yfact_isv(:,k)
end do

if ( mlodiff==0 ) then
  ! Laplacian diffusion terms (closure #1)
  do k = 1,wlev
    duma(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(duma(:,:,1:3))
  do k = 1,wlev
    call unpack_nsew(duma(:,k,1),duma_n,duma_s,duma_e,duma_w)  
    nu = ( duma(1:ifull,k,1)*emi +                      &
           xfact(1:ifull,k)*duma_e +                    &
           xfact_iwu(1:ifull,k)*duma_w +                &
           yfact(1:ifull,k)*duma_n +                    &
           yfact_isv(1:ifull,k)*duma_s ) / base(:,k)
    call unpack_nsew(duma(:,k,2),duma_n,duma_s,duma_e,duma_w)  
    nv = ( duma(1:ifull,k,2)*emi +                      &
           xfact(1:ifull,k)*duma_e +                    &
           xfact_iwu(1:ifull,k)*duma_w +                &
           yfact(1:ifull,k)*duma_n +                    &
           yfact_isv(1:ifull,k)*duma_s ) / base(:,k)
    call unpack_nsew(duma(:,k,3),duma_n,duma_s,duma_e,duma_w)  
    nw = ( duma(1:ifull,k,3)*emi +                      &
           xfact(1:ifull,k)*duma_e +                    &
           xfact_iwu(1:ifull,k)*duma_w +                &
           yfact(1:ifull,k)*duma_n +                    &
           yfact_isv(1:ifull,k)*duma_s ) / base(:,k)
    outu(1:ifull,k) = ax(1:ifull)*nu + ay(1:ifull)*nv + az(1:ifull)*nw
    outv(1:ifull,k) = bx(1:ifull)*nu + by(1:ifull)*nv + bz(1:ifull)*nw
  end do

!  Laplacian diffusion and viscosity terms (closure #2)
!  duma(1:ifull,:,1)=u(1:ifull,:)
!  duma(1:ifull,:,2)=v(1:ifull,:)
!  call boundsuv(duma(:,:,1),duma(:,:,2),allvec=.true.)
!  do k=1,wlev
!    outu(:,k)=(duma(1:ifull,k,1)*emi+2.*xfact(1:ifull,k)*duma(ieu,k,1)+2.*xfact(iwu,k)*duma(iwu,k,1) &
!      +yfact(1:ifull,k)*duma(inu,k,1)+yfact(isv,k)*duma(isu,k,1)                                     &
!      +(yfact(1:ifull,k)-yfact(isv,k))*0.5*(duma(iev,k,2)-duma(iwv,k,2))                             &
!      +t_kh(1:ifull,k)*0.5*(duma(inv,k,2)+duma(iev,k,2)-duma(isv,k,2)-duma(iwv,k,2)))                &
!      /base(:,k)
!    outv(:,k)=(duma(1:ifull,k,2)*emi+2.*yfact(1:ifull,k)*duma(inv,k,2)+2.*yfact(isv,k)*duma(isv,k,2) &
!      +xfact(1:ifull,k)*duma(iev,k,2)+xfact(iwu,k)*duma(iwv,k,2)                                     &
!      +(xfact(1:ifull,k)-xfact(iwu,k))*0.5*(duma(inu,k,1)-duma(isu,k,1))                             &
!      +t_kh(1:ifull,k)*0.5*(duma(inu,k,1)+duma(ieu,k,1)-duma(isu,k,1)-duma(iwu,k,1)))                &
!      /base(:,k)
!  end do

  call mloimport3d(2,outu,0)
  call mloimport3d(3,outv,0)
  
else if ( mlodiff==1 ) then
  outu(1:ifull,:) = u(1:ifull,:)
  outv(1:ifull,:) = v(1:ifull,:)
  ! no need to call mloimport3d for u and v
else
  write(6,*) "ERROR: Unknown option for mlodiff = ",mlodiff
  call ccmpi_abort(-1)
end if
  
! Potential temperature and salinity
duma(1:ifull,:,1) = tt(1:ifull,:)
duma(1:ifull,:,2) = ss(1:ifull,:) - 34.72
call bounds(duma(:,:,1:2))
do k=1,wlev
  call unpack_nsew(duma(:,k,1),duma_n,duma_s,duma_e,duma_w)  
  fs(:,k) = ( duma(1:ifull,k,1)*emi +                      &
              xfact(1:ifull,k)*duma_e +                    &
              xfact_iwu(1:ifull,k)*duma_w +                &
              yfact(1:ifull,k)*duma_n +                    &
              yfact_isv(1:ifull,k)*duma_s ) / base(:,k)
end do
fs = max(fs, -wrtemp)
call mloimport3d(0,fs,0)
do k=1,wlev
  call unpack_nsew(duma(:,k,2),duma_n,duma_s,duma_e,duma_w)    
  fs(:,k) = ( duma(1:ifull,k,2)*emi +                      &
              xfact(1:ifull,k)*duma_e +                    &
              xfact_iwu(1:ifull,k)*duma_w +                &
              yfact(1:ifull,k)*duma_n +                    &
              yfact_isv(1:ifull,k)*duma_s ) / base(:,k)
end do
fs = max(fs+34.72, 0.)
call mloimport3d(1,fs,0)

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements some basic hydrostatic dynamics for the
! ocean and ice.  The ocean component employs the R-grid design used
! in CCAM semi-Lagragain dynamics, but uses sigma-z vertical
! coordinates.  The staggered/unstaggered pivoting has been modified
! for the land-sea boundary.  Sea-ice is advected using an upwind
! scheme.  Internal sea-ice pressure follows a cavitating fluid
! apporximation.
subroutine mlohadv

use arrays_m
use cc_mpi
use const_phys
use infile
use indices_m
use latlong_m
use map_m
use mlo
use mlodynamicsarrays_m
use newmpar_m
use nharrs_m, only : lrestart
use parm_m
use parmdyn_m
use soil_m
use soilsnow_m
use vecsuv_m

implicit none

integer ii,totits
integer jyear,jmonth,jday,jhour,jmin,mins
integer tyear,jstart, iq, mspec_mlo, mspeca_mlo
integer, dimension(ifull,wlev) :: nface
real maxglobseta,maxglobip,hdt,dtin_mlo
real alph_p
real, dimension(2) :: delpos, delneg
real, dimension(ifull+iextra) :: neta,pice,imass,xodum
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv
real, dimension(ifull+iextra) :: snu,sou,spu,squ,ssu,snv,sov,spv,sqv,ssv
real, dimension(ifull+iextra) :: ibu,ibv,icu,icv,spnet,oeu,oev,tide
real, dimension(ifull+iextra) :: ipmax
real, dimension(ifull) :: i_u,i_v,i_sto,ndum
real, dimension(ifull) :: w_e,xps
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: dipdxu,dipdyu,dipdxv,dipdyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv,odum
real, dimension(ifull) :: imu,imv
real, dimension(ifull) :: dumf,dumg
real, dimension(ifull) :: dd_e, dd_n
real, dimension(ifull) :: ddu_e, ddv_n, ddu_w, ddv_s
real, dimension(ifull) :: pice_n, pice_e, pice_s, pice_w
real, dimension(ifull) :: tide_n, tide_s, tide_e, tide_w
real, dimension(ifull) :: imass_n, imass_e
real, dimension(ifull) :: neta_n, neta_s, neta_e, neta_w
real, dimension(ifull) :: ipice_n, ipice_s, ipice_e, ipice_w
real, dimension(ifull) :: dd_isv, dd_iwu, em_isv, em_iwu
real, dimension(ifull) :: oev_isv, oeu_iwu, cc_isv, cc_iwu
real, dimension(ifull) :: eo_isv, eo_iwu, ni_isv, ni_iwu
real, dimension(ifull) :: sp_isv, sp_iwu, sq_isv, sq_iwu, qu_isv, qu_iwu
real, dimension(ifull) :: so_isv, so_iwu, ss_isv, ss_iwu
real, dimension(ifull) :: dnetadx, dnetady, ddddx, ddddy
real, dimension(ifull) :: sdiv
real, dimension(ifull+iextra,wlev,3) :: cou
real, dimension(ifull+iextra,wlev+1) :: eou,eov
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns,dzdum_rho
real, dimension(ifull+iextra,wlev) :: dalpha,dbeta
real, dimension(ifull+iextra,wlev) :: ccu,ccv
real, dimension(ifull+iextra,9) :: dumc,dumd
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull,wlev,2) :: mfixdum
real, dimension(ifull,wlev+1) :: tau,tav,ttau,ttav
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv
real, dimension(ifull,wlev) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull,wlev) :: depdum,dzdum
real, dimension(ifull,wlev) :: dd_adv,mps
real, dimension(ifull,0:wlev) :: nw
real, dimension(ifull,4) :: i_it
real, dimension(ifull,3) :: gamm
real, dimension(wlev) :: neg_godsig
real(kind=8), dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical lleap

! Vertical coordinates are defined as:
!   sigma = (z+neta)/(D+neta)
! which was optimised for coastal modelling problems.

! Instead of packing ocean points, we keep the original
! grid and define a land/sea mask.  This allows us to
! use the same index and map factor arrays from the
! atmospheric dynamical core.

! save time-step
dtin_mlo = dt
mspeca_mlo = 1

! Define land/sea mask
wtr(:) = ee(:)>0.5

! Default values
w_t   = 293.16-wrtemp ! potential water temperature delta at tau=t
w_s   = 34.72         ! water salinity at tau=t
w_u   = 0.            ! u component of water current at tau=t
w_v   = 0.            ! v component of water current at tau=t
w_e   = 0.            ! free surface height at tau=t
i_it  = 273.16        ! ice temperature
i_sto = 0.            ! ice brine storage
i_u   = 0.            ! u component of ice velocity
i_v   = 0.            ! v component of ice velocity
nw    = 0.            ! water vertical velocity
pice  = 0.            ! ice pressure for cavitating fluid
imass = 0.            ! ice mass
tide  = 0.            ! tidal forcing
nt    = 0.            ! new water temperature
ns    = 0.            ! new water salinity
nu    = 0.            ! new u component of water current
nv    = 0.            ! new v component of water current
neta  = 0.            ! new free surface height
cou   = 0.            ! working array
dumc  = 0.            ! working array
eou   = 0.            ! working array
eov   = 0.            ! working array
niu   = 0.            ! new u component of ice velocity
niv   = 0.            ! new v component of ice velocity
nfracice = 0.
ndic     = 0.
ndsn     = 0.

! EXPORT WATER AND ICE DATA FROM MLO ------------------------------------------
call mloexport3d(0,w_t,0)
call mloexport3d(1,w_s,0)
call mloexport3d(2,w_u,0)
call mloexport3d(3,w_v,0)
call mloexport(4,w_e,0,0)
do ii = 1,4
  call mloexpice(i_it(:,ii),ii,0)
end do
call mloexpice(nfracice,5,0)
call mloexpice(ndic,6,0)
call mloexpice(ndsn,7,0)
call mloexpice(i_sto,8,0)
call mloexpice(i_u,9,0)
call mloexpice(i_v,10,0)
where (wtr(1:ifull))
  fracice(1:ifull) = nfracice(1:ifull)  
  sicedep(1:ifull) = ndic(1:ifull)
  snowd(1:ifull)   = ndsn(1:ifull)*1000.
end where  

! estimate tidal forcing (assumes leap days)
if ( usetide==1 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins,allleap=.true.)
  jstart=0
  if ( jyear>1900 ) then
    do tyear=1900,jyear-1
      call mloleap(tyear,lleap)
      if ( lleap ) jstart=jstart+1
      jstart=jstart+365
    end do
  else if ( jyear<1900 ) then
    do tyear=1899,jyear,-1
      call mloleap(tyear,lleap)
      if ( lleap ) jstart=jstart-1
      jstart=jstart-365
    end do
  end if
  mins=mins+720 ! base time is 12Z 31 Dec 1899
  call mlotide(tide,rlongg,rlatt,mins,jstart)
end if

! initialise t+1 variables with t data
neta(1:ifull)=w_e
nu(1:ifull,:)=w_u
nv(1:ifull,:)=w_v
nt(1:ifull,:)=w_t
ns(1:ifull,:)=w_s
nit(1:ifull,:)=i_it
nsto(1:ifull)=i_sto
niu(1:ifull)=i_u
niv(1:ifull)=i_v

! initialise save arrays
if ( ktau==1 .and. .not.lrestart ) then
  oldu1 = nu(1:ifull,:)
  oldv1 = nv(1:ifull,:)
  oldu2 = nu(1:ifull,:)
  oldv2 = nv(1:ifull,:)
  ipice = 0. ! ice free drift solution
  mspeca_mlo = 2
  dt = 0.5*dtin_mlo
end if

#ifdef mlodebug
if ( any( nt(1:ifull,:)+wrtemp<100. .or. nt(1:ifull,:)+wrtemp>400. ) ) then
  write(6,*) "ERROR: nt is out of range at start of mlodynamics"
  write(6,*) "minval,maxval ",minval(nt(1:ifull,:)),maxval(nt(1:ifull,:))
  write(6,*) "minloc,maxloc ",minloc(nt(1:ifull,:)),maxloc(nt(1:ifull,:))
  call ccmpi_abort(-1)
end if
if ( any(abs(nu(1:ifull,:))>20.) .or. any(abs(nv(1:ifull,:))>20.) ) then
  write(6,*) "ERROR: current out-of-range at start of mlodynamics"
  write(6,*) "u ",minval(nu(1:ifull,:)),maxval(nu(1:ifull,:))
  write(6,*) "v ",minval(nv(1:ifull,:)),maxval(nv(1:ifull,:))
  stop
end if
#endif

do mspec_mlo = mspeca_mlo,1,-1

  hdt = 0.5*dt  
    
  if ( mspeca_mlo==2 .and. mspec_mlo==1 ) then
    nu(1:ifull,:) = oldu1
    nv(1:ifull,:) = oldv1
    w_t = nt(1:ifull,:)
    w_s = ns(1:ifull,:)
    w_e = neta(1:ifull)
  end if

  ! surface pressure and ice mass
  ! (assume ice velocity is 'slow' compared to 'fast' change in neta)
  imass(1:ifull) = ndic(1:ifull)*rhoic + ndsn(1:ifull)*rhosn  ! ice mass per unit area (kg/m^2), unstaggered at time t
  pice(1:ifull) = ps(1:ifull)                                 ! pressure due to atmosphere at top of water column (unstaggered at t)
  if ( usepice==1 ) then
    pice(1:ifull) = pice(1:ifull) + grav*nfracice(1:ifull)*imass(1:ifull) ! include ice in surface pressure
  end if
  ! Limit minimum ice mass for ice velocity calculation.  Hence, we can solve for the ice velocity at
  ! grid points where the ice is not yet present.  
  imass(1:ifull) = max( imass(1:ifull), 100. )
  ! maximum pressure for cavitating fluid
  select case(icemode)
    case(2)
      ! cavitating fluid
      ipmax(1:ifull) = 27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))*ee(1:ifull)
    case(1)
      ! incompressible fluid
      ipmax(1:ifull) = 9.E9*ee(1:ifull)
    case DEFAULT
      ! free drift
      ipmax(1:ifull) = 0.
  end select

  ! update scalar bounds and various gradients (aggregate fields for MPI)
  dumc(1:ifull,1) = neta(1:ifull)
  dumc(1:ifull,2) = pice(1:ifull)
  dumc(1:ifull,3) = tide(1:ifull)
  dumc(1:ifull,4) = imass(1:ifull)
  dumc(1:ifull,5) = ipmax(1:ifull)
  call bounds(dumc(:,1:5),corner=.true.)
  neta(ifull+1:ifull+iextra)  = dumc(ifull+1:ifull+iextra,1)
  pice(ifull+1:ifull+iextra)  = dumc(ifull+1:ifull+iextra,2)
  tide(ifull+1:ifull+iextra)  = dumc(ifull+1:ifull+iextra,3)
  imass(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,4)
  ipmax(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,5)

  call unpack_nsew(pice,pice_n,pice_s,pice_e,pice_w)
  call unpack_nsew(tide,tide_n,tide_s,tide_e,tide_w)
  call unpack_ne(imass,imass_n,imass_e)
  call unpack_ne(neta,neta_n,neta_e)
  call unpack_ne(dd,dd_n,dd_e)
  call unpack_svwu(ddu,ddv,dd_isv,dd_iwu)
  call unpack_svwu(emu,emv,em_isv,em_iwu)

  ! surface pressure
!$omp simd
  do iq = 1,ifull
    tnu(iq) = 0.5*( pice_n(iq) + pice(ien(iq)) )
    tee(iq) = 0.5*( pice(iq)   + pice_e(iq)    )
    tsu(iq) = 0.5*( pice_s(iq) + pice(ies(iq)) )
    tev(iq) = 0.5*( pice_e(iq) + pice(ine(iq)) )
    tnn(iq) = 0.5*( pice(iq)   + pice_n(iq)    )
    twv(iq) = 0.5*( pice_w(iq) + pice(inw(iq)) )
  end do  
  dpsdxu = (pice_e-pice(1:ifull))*emu(1:ifull)/ds ! staggered at time t
  dpsdyu = 0.5*(stwgt(1:ifull,1)*(tnu-tee)+stwgt(1:ifull,2)*(tee-tsu))*emu(1:ifull)/ds
  dpsdxv = 0.5*(stwgt(1:ifull,3)*(tev-tnn)+stwgt(1:ifull,4)*(tnn-twv))*emv(1:ifull)/ds
  dpsdyv = (pice_n-pice(1:ifull))*emv(1:ifull)/ds

  ! tides
!$omp simd
  do iq = 1,ifull
    tnu(iq) = 0.5*( tide_n(iq) + tide(ien(iq)) )
    tee(iq) = 0.5*( tide(iq)   + tide_e(iq)    )
    tsu(iq) = 0.5*( tide_s(iq) + tide(ies(iq)) )
    tev(iq) = 0.5*( tide_e(iq) + tide(ine(iq)) )
    tnn(iq) = 0.5*( tide(iq)   + tide_n(iq)    )
    twv(iq) = 0.5*( tide_w(iq) + tide(inw(iq)) )
  end do  
  dttdxu = (tide_e-tide(1:ifull))*emu(1:ifull)/ds ! staggered
  dttdyu = 0.5*(stwgt(1:ifull,1)*(tnu-tee)+stwgt(1:ifull,2)*(tee-tsu))*emu(1:ifull)/ds
  dttdxv = 0.5*(stwgt(1:ifull,3)*(tev-tnn)+stwgt(1:ifull,4)*(tnn-twv))*emv(1:ifull)/ds
  dttdyv = (tide_n-tide(1:ifull))*emv(1:ifull)/ds


  call START_LOG(watereos_begin)

  ! Calculate adjusted depths and thicknesses
  !xodum = max( dd(:)+neta(:), minwater )
  xodum = dd(:)
  do ii = 1,wlev
    depdum(1:ifull,ii)           = gosig(ii)*xodum(1:ifull)
    dzdum(1:ifull,ii)            = godsig(ii)*xodum(1:ifull)
    dzdum_rho(1:ifull+iextra,ii) = godsig(ii)*dd(1:ifull+iextra)
  end do

  ! Calculate normalised density coeffs dalpha and dbeta (unstaggered at time t)
  ! (Assume free surface correction is small so that changes in the compression 
  ! effect due to neta can be neglected.  Consequently, the neta dependence is 
  ! separable in the iterative loop)
  call bounds(nt,corner=.true.)
  call bounds(ns,corner=.true.)
  ! rho is the pertubation from wrtrho, which we assume is much smaller than wrtrho
  call mloexpdensity(ccu,dalpha,dbeta,nt,ns,dzdum_rho,pice,0,rawrho=.true.) ! rho=ccu (not used)

  ! Calculate normalised density gradients
  ! method 2: Use potential temperature and salinity Jacobians (see Shchepetkin and McWilliams 2003)
  call tsjacobi(nt,ns,dalpha,dbeta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
  drhobardxu(:,1) = drhobardxu(:,1)*godsig(1)
  drhobardxv(:,1) = drhobardxv(:,1)*godsig(1)
  drhobardyu(:,1) = drhobardyu(:,1)*godsig(1)
  drhobardyv(:,1) = drhobardyv(:,1)*godsig(1)
  do ii = 2,wlev
    drhobardxu(:,ii) = drhobardxu(:,ii-1) + drhobardxu(:,ii)*godsig(ii)
    drhobardxv(:,ii) = drhobardxv(:,ii-1) + drhobardxv(:,ii)*godsig(ii)
    drhobardyu(:,ii) = drhobardyu(:,ii-1) + drhobardyu(:,ii)*godsig(ii)
    drhobardyv(:,ii) = drhobardyv(:,ii-1) + drhobardyv(:,ii)*godsig(ii)
  end do
  do ii = 1,wlev
    drhobardxu(:,ii) = drhobardxu(:,ii)/gosigh(ii)
    drhobardxv(:,ii) = drhobardxv(:,ii)/gosigh(ii)
    drhobardyu(:,ii) = drhobardyu(:,ii)/gosigh(ii)
    drhobardyv(:,ii) = drhobardyv(:,ii)/gosigh(ii)
  end do

  call END_LOG(watereos_end)


  ! ADVECT WATER AND ICE ----------------------------------------------
  ! Water currents are advected using semi-Lagrangian advection
  ! based on McGregor's CCAM advection routines.
  ! Velocity is set to zero at ocean boundaries.
  
  ! Estimate currents at t+1/2 for semi-Lagrangian advection
  do ii = 1,wlev
    nuh(:,ii) = (15.*nu(1:ifull,ii)-10.*oldu1(:,ii)+3.*oldu2(:,ii))*ee(1:ifull)/8. ! U at t+1/2
    nvh(:,ii) = (15.*nv(1:ifull,ii)-10.*oldv1(:,ii)+3.*oldv2(:,ii))*ee(1:ifull)/8. ! V at t+1/2
  end do
  nuh = min( max( nuh, -50. ), 50. )
  nvh = min( max( nvh, -50. ), 50. )

  ! save arrays for extrapolating currents at next time-step
  oldu2(1:ifull,1:wlev) = oldu1(1:ifull,1:wlev)
  oldv2(1:ifull,1:wlev) = oldv1(1:ifull,1:wlev)
  oldu1(1:ifull,1:wlev) = nu(1:ifull,1:wlev)
  oldv1(1:ifull,1:wlev) = nv(1:ifull,1:wlev)


  ! Define horizontal transport
  ! dH(phi)/dt = d(phi)/dt + u*d(phi)/dx + v*d(phi)/dy
  
  ! Continuity equation
  ! 1/(neta+D) dH(neta+D)/dt + du/dx + dv/dy + dw/dz = 0
  ! dH(neta+D)/dt + (neta+D)*(du/dx+dv/dy) + dw/dsig = 0
  ! dH(neta)/dt + neta*(du/dx+dv/dy) + d(D*u)/dx + d(D*v)/dy + dw/dsig = 0  ! where D is now included in flux form
  ! dH(neta)/dt + d(D*u)/dx + d(D*v)/dy + dw/dsig = 0 (approximate)
  
  ! Lagrangian version of the contunity equation with D included in flux form
  !   [neta + (1+eps)*0.5*dt*(d(D*u)/dx+d(D*v)/dy+dw/dsig)]^(t+1)
  ! = [neta - (1-eps)*0.5*dt*(d(D*u)/dx+d(D*v)/dy+dw/dsig)]^(t*)
  
  ! Vertical velocity is then (sum_0_sig is the integral from 0 to sig)

  ! w = -(neta+D)*sum_0_sig(du/dx+dv/dy,dsig) + (neta+D)*sig*sum_0_1(du/dx+dv/dy,dsig)
  !     -sum_0_sig(u*d(neta+D)/dx+v*d(neta+D)/dy),dsig) + sum_0_1(u*d(neta+D)/dx+v*d(neta+D)/dy,dsig)
  !   = -sum_0_sig(d(u*(neta+D))/dx+d(v*(neta+D))/dy,dsig) + sig*sum_0_1(d(u*(neta+D))/dx+d(v*(neta+D))/dy,dsig)
  !   = -sum_0_sig(d(u*D)/dx+d(v*D)/dy,dsig) + sig*sum_0_1(d(u*D)/dx+d(v*D)/dy,dsig)  (approximate)

  ! Advect continuity equation to tstar
  ! Calculate velocity on C-grid for consistancy with iterative free surface calculation
  ! nw is at half levels with nw(:,1) at the bottom of level 1
  ! positive nw is moving downwards to the ocean floor
  ! Assume Boussinesq fluid and treat density in continuity equation as constant
  ! true vertical velocity = nw-u*((1-sig)*deta/dx-sig*d(dd)/dx)-v*((1-sig)*deta/dy-sig*d(dd)/dy)-(1-sig)*deta/dt
  call mlostaguv(nu(:,1:wlev),nv(:,1:wlev),eou(:,1:wlev),eov(:,1:wlev))
  ! surface height at staggered coordinate
  eou(1:ifull,wlev+1) = 0.5*(neta(1:ifull)+neta_e)*eeu(1:ifull) ! height at staggered coordinate
  eov(1:ifull,wlev+1) = 0.5*(neta(1:ifull)+neta_n)*eev(1:ifull) ! height at staggered coordinate
  call boundsuv(eou(:,1:wlev+1),eov(:,1:wlev+1),stag=-9)
  oeu(1:ifull+iextra) = eou(1:ifull+iextra,wlev+1)
  oev(1:ifull+iextra) = eov(1:ifull+iextra,wlev+1)
  ! vertically integrate currents
  ccu(:,1) = eou(:,1)*godsig(1)
  ccv(:,1) = eov(:,1)*godsig(1)
  do ii = 2,wlev
    ccu(:,ii) = ccu(:,ii-1) + eou(:,ii)*godsig(ii)
    ccv(:,ii) = ccv(:,ii-1) + eov(:,ii)*godsig(ii)
  end do
  ! calculate vertical velocity (use flux form)
  call unpack_svwu(oeu,oev,oev_isv,oeu_iwu)
  call unpack_svwu(ccu(:,wlev),ccv(:,wlev),cc_isv,cc_iwu)
  sdiv(:) = (ccu(1:ifull,wlev)*ddu(1:ifull)/emu(1:ifull)    &
            -cc_iwu*dd_iwu/em_iwu                           &
            +ccv(1:ifull,wlev)*ddv(1:ifull)/emv(1:ifull)    &
            -cc_isv*dd_isv/em_isv)                          &
            *em(1:ifull)**2/ds
  do ii = 1,wlev-1
    call unpack_svwu(ccu(:,ii),ccv(:,ii),cc_isv,cc_iwu)  
    nw(:,ii) = ee(1:ifull)*(sdiv(:)*gosigh(ii)              &
               - (ccu(1:ifull,ii)*ddu(1:ifull)/emu(1:ifull) &
                 -cc_iwu*dd_iwu/em_iwu                      &
                 +ccv(1:ifull,ii)*ddv(1:ifull)/emv(1:ifull) &
                 -cc_isv*dd_isv/em_isv)                     &
                 *em(1:ifull)**2/ds)
  end do

  ! compute contunity equation horizontal transport terms
  do ii = 1,wlev
    call unpack_svwu(eou(:,ii),eov(:,ii),eo_isv,eo_iwu)  
    mps(1:ifull,ii) = neta(1:ifull) - (1.-ocneps)*0.5*dt*ee(1:ifull) &
                     *((eou(1:ifull,ii)*ddu(1:ifull)/emu(1:ifull)    &
                       -eo_iwu*dd_iwu/em_iwu                         &
                       +eov(1:ifull,ii)*ddv(1:ifull)/emv(1:ifull)    &
                       -eo_isv*dd_isv/em_isv)                        &
                       *em(1:ifull)**2/ds                            &
                      +(nw(:,ii)-nw(:,ii-1))/godsig(ii))
  end do
      
  ! update vertical velocity at full levels
  dnetadx = (oeu(1:ifull)/emu(1:ifull)-oeu_iwu/em_iwu)*em(1:ifull)**2/ds
  dnetady = (oev(1:ifull)/emv(1:ifull)-oev_isv/em_isv)*em(1:ifull)**2/ds
  ddddx = (ddu(1:ifull)/emu(1:ifull)-dd_iwu/em_iwu)*em(1:ifull)**2/ds
  ddddy = (ddv(1:ifull)/emv(1:ifull)-dd_isv/em_isv)*em(1:ifull)**2/ds
  do ii = 1,wlev
    w_ocn(:,ii) = ee(1:ifull)*(0.5*(nw(:,ii-1)+nw(:,ii))                       &
                   - nu(1:ifull,ii)*((1.-gosig(ii))*dnetadx + gosig(ii)*ddddx) &
                   - nv(1:ifull,ii)*((1.-gosig(ii))*dnetady + gosig(ii)*ddddy))
                  !- (1.-gosig(ii))*dnetadt ! neglect for now 
  end do  

  
  ! ocean
  ! Prepare pressure gradient terms at t=t and incorporate into velocity field
  detadxu = (neta_e-neta(1:ifull))*emu(1:ifull)/ds
  detadyv = (neta_n-neta(1:ifull))*emv(1:ifull)/ds
  do ii = 1,wlev
    tau(:,ii) = grav*gosig(ii)*ddu(1:ifull)*drhobardxu(:,ii)/wrtrho  &
               + dpsdxu/wrtrho + grav*dttdxu + grav*detadxu ! staggered
    tav(:,ii) = grav*gosig(ii)*ddv(1:ifull)*drhobardyv(:,ii)/wrtrho  &
               + dpsdyv/wrtrho + grav*dttdyv + grav*detadyv  
  end do
  ! ice
  !tau(:,wlev+1)=grav*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds ! staggered
  !tav(:,wlev+1)=grav*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
  !call mlounstaguv(tau,tav,ttau,ttav,toff=1)
  call mlounstaguv(tau(:,1:wlev),tav(:,1:wlev),ttau(:,1:wlev),ttav(:,1:wlev),toff=1)
  ! ocean
  do ii = 1,wlev
    uau(:,ii) = nu(1:ifull,ii) + (1.-ocneps)*0.5*dt*( f(1:ifull)*nv(1:ifull,ii)-ttau(:,ii)) ! unstaggered
    uav(:,ii) = nv(1:ifull,ii) + (1.-ocneps)*0.5*dt*(-f(1:ifull)*nu(1:ifull,ii)-ttav(:,ii))
    uau(:,ii) = uau(:,ii)*ee(1:ifull)
    uav(:,ii) = uav(:,ii)*ee(1:ifull)
  end do
  ! ice
  snu(1:ifull) = i_u !-dt*ttau(:,wlev+1)
  snv(1:ifull) = i_v !-dt*ttav(:,wlev+1)
  

  ! Vertical advection (first call for 0.5*dt)
  call mlovadv(hdt,nw,uau,uav,ns,nt,mps,depdum,dzdum,wtr(1:ifull),1)

  
#ifdef mlodebug
  if ( any( nt(1:ifull,:)+wrtemp<100. .or. nt(1:ifull,:)+wrtemp>400. ) ) then
    write(6,*) "ERROR: nt is out of range after first vertical advection"
    write(6,*) "minval,maxval ",minval(nt(1:ifull,:)),maxval(nt(1:ifull,:))
    write(6,*) "minloc,maxloc ",minloc(nt(1:ifull,:)),maxloc(nt(1:ifull,:))
    call ccmpi_abort(-1)
  end if
  if ( any(abs(uau(1:ifull,:))>20.) .or. any(abs(uav(1:ifull,:))>20.) ) then
    write(6,*) "ERROR: current out-of-range after vertical advection 1"
    write(6,*) "u ",minval(uau(1:ifull,:)),maxval(uau(1:ifull,:))
    write(6,*) "v ",minval(uav(1:ifull,:)),maxval(uav(1:ifull,:))
    stop
  end if
#endif


  ! Calculate depature points
  call mlodeps(dt,nuh,nvh,nface,xg,yg,x3d,y3d,z3d,wtr)

  ! Convert (u,v) to cartesian coordinates (U,V,W)
  do ii = 1,wlev
    cou(1:ifull,ii,1) = ax(1:ifull)*uau(:,ii) + bx(1:ifull)*uav(:,ii)
    cou(1:ifull,ii,2) = ay(1:ifull)*uau(:,ii) + by(1:ifull)*uav(:,ii)
    cou(1:ifull,ii,3) = az(1:ifull)*uau(:,ii) + bz(1:ifull)*uav(:,ii)
  end do
  ! Horizontal advection for U, V, W
  call mlob2ints_uv(cou(:,:,1:3),nface,xg,yg,wtr)
  ! Rotate vector to arrival point
  call mlorot(cou(:,:,1),cou(:,:,2),cou(:,:,3),x3d,y3d,z3d)
  ! Convert (U,V,W) back to conformal cubic coordinates
  do ii = 1,wlev
    uau(:,ii) = ax(1:ifull)*cou(1:ifull,ii,1) + ay(1:ifull)*cou(1:ifull,ii,2) + az(1:ifull)*cou(1:ifull,ii,3)
    uav(:,ii) = bx(1:ifull)*cou(1:ifull,ii,1) + by(1:ifull)*cou(1:ifull,ii,2) + bz(1:ifull)*cou(1:ifull,ii,3)
    uau(:,ii) = uau(:,ii)*ee(1:ifull)
    uav(:,ii) = uav(:,ii)*ee(1:ifull)
  end do

  ! Horizontal advection for continuity
  cou(1:ifull,1:wlev,1) = mps(1:ifull,1:wlev)
  call mlob2ints(cou(:,:,1:1),nface,xg,yg,wtr)
  mps(1:ifull,1:wlev) = cou(1:ifull,1:wlev,1)

  ! Horizontal advection for T and S
  do ii = 1,wlev
    cou(1:ifull,ii,1) = nt(1:ifull,ii)
    cou(1:ifull,ii,2) = ns(1:ifull,ii) - 34.72
  end do
  call mlob2ints_bs(cou(:,:,1:2),nface,xg,yg,wtr)
  do ii = 1,wlev
    nt(1:ifull,ii) = max( cou(1:ifull,ii,1), -wrtemp )
    ns(1:ifull,ii) = max( cou(1:ifull,ii,2) + 34.72, 0. )
  end do


#ifdef mlodebug
  if ( any( nt(1:ifull,:)+wrtemp<100. .or. nt(1:ifull,:)+wrtemp>400. ) ) then
    write(6,*) "ERROR: nt is out of range after horizontal advection"
    write(6,*) "minval,maxval ",minval(nt(1:ifull,:)),maxval(nt(1:ifull,:))
    write(6,*) "minloc,maxloc ",minloc(nt(1:ifull,:)),maxloc(nt(1:ifull,:))
    call ccmpi_abort(-1)
  end if
  if ( any(abs(uau(1:ifull,:))>20.) .or. any(abs(uav(1:ifull,:))>20.) ) then
    write(6,*) "ERROR: current out-of-range after horizontal advection"
    write(6,*) "u ",minval(uau(1:ifull,:)),maxval(uau(1:ifull,:))
    write(6,*) "v ",minval(uav(1:ifull,:)),maxval(uav(1:ifull,:))
    stop
  end if
#endif


  ! Vertical advection (second call for 0.5*dt)
  ! use explicit nw and depdum,dzdum from t=tau step (i.e., following JLM in CCAM atmospheric dynamics)
  ! Could use nuh and nvh to estimate nw at t+1/2, but also require an estimate of neta at t+1/2
  call mlovadv(hdt,nw,uau,uav,ns,nt,mps,depdum,dzdum,wtr(1:ifull),2)


#ifdef mlodebug
  if ( any( nt(1:ifull,:)+wrtemp<100. .or. nt(1:ifull,:)+wrtemp>400. ) ) then
    write(6,*) "ERROR: nt is out of range after second vertical advection"
    write(6,*) "minval,maxval ",minval(nt(1:ifull,:)),maxval(nt(1:ifull,:))
    write(6,*) "minloc,maxloc ",minloc(nt(1:ifull,:)),maxloc(nt(1:ifull,:))
    call ccmpi_abort(-1)
  end if
  if ( any(abs(uau(1:ifull,:))>20.) .or. any(abs(uav(1:ifull,:))>20.) ) then
    write(6,*) "ERROR: current out-of-range after vertical advection 2"
    write(6,*) "u ",minval(uau(1:ifull,:)),maxval(uau(1:ifull,:))
    write(6,*) "v ",minval(uav(1:ifull,:)),maxval(uav(1:ifull,:))
    stop
  end if
#endif

  xps(:) = mps(1:ifull,1)*godsig(1)
  do ii = 2,wlev
    xps(:) = xps(:) + mps(1:ifull,ii)*godsig(ii)
  end do
  xps(:) = xps(:)*ee(1:ifull)


  call START_LOG(watereos_begin)

  ! Approximate normalised density rhobar at t+1 (unstaggered, using T and S at t+1)
  if ( nxtrrho==1 ) then
    call bounds(nt,corner=.true.)
    call bounds(ns,corner=.true.)
    ! rho is the pertubation from wrtrho
    call mloexpdensity(ccu,dalpha,dbeta,nt,ns,dzdum_rho,pice,0,rawrho=.true.) ! rho=ccu

    ! update normalised density gradients
    call tsjacobi(nt,ns,dalpha,dbeta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
    drhobardxu(:,1) = drhobardxu(:,1)*godsig(1)
    drhobardxv(:,1) = drhobardxv(:,1)*godsig(1)
    drhobardyu(:,1) = drhobardyu(:,1)*godsig(1)
    drhobardyv(:,1) = drhobardyv(:,1)*godsig(1)
    do ii = 2,wlev
      drhobardxu(:,ii) = drhobardxu(:,ii-1) + drhobardxu(:,ii)*godsig(ii)
      drhobardxv(:,ii) = drhobardxv(:,ii-1) + drhobardxv(:,ii)*godsig(ii)
      drhobardyu(:,ii) = drhobardyu(:,ii-1) + drhobardyu(:,ii)*godsig(ii)
      drhobardyv(:,ii) = drhobardyv(:,ii-1) + drhobardyv(:,ii)*godsig(ii)
    end do
    do ii = 1,wlev
      drhobardxu(:,ii) = drhobardxu(:,ii)/gosigh(ii)
      drhobardxv(:,ii) = drhobardxv(:,ii)/gosigh(ii)
      drhobardyu(:,ii) = drhobardyu(:,ii)/gosigh(ii)
      drhobardyv(:,ii) = drhobardyv(:,ii)/gosigh(ii)
    end do
  end if

  call END_LOG(watereos_end)


  ! FREE SURFACE CALCULATION ----------------------------------------

  ! Prepare integral terms
  sou = 0.
  spu = 0.
  squ = 0.
  ssu = 0.
  sov = 0.
  spv = 0.
  sqv = 0.
  ssv = 0.

  ! ocean
  ! Precompute U,V current and integral terms at t+1
  do ii = 1,wlev
    tau(:,ii) = uau(:,ii) + (1.+ocneps)*0.5*dt*f(1:ifull)*uav(:,ii)
    tav(:,ii) = uav(:,ii) - (1.+ocneps)*0.5*dt*f(1:ifull)*uau(:,ii)
  end do
  ! ice
  tau(:,wlev+1) = snu(1:ifull) + dt*f(1:ifull)*snv(1:ifull) ! unstaggered
  tav(:,wlev+1) = snv(1:ifull) - dt*f(1:ifull)*snu(1:ifull)
  call mlostaguv(tau(:,1:wlev+1),tav(:,1:wlev+1),ttau(:,1:wlev+1),ttav(:,1:wlev+1))
  ! ocean
  odum = 1./(1.+((1.+ocneps)*0.5*dt*fu(1:ifull))**2)
  odum = odum*eeu(1:ifull)
  dumf = -(1.+ocneps)*0.5*dt*odum
  do ii = 1,wlev
    ccu(1:ifull,ii) = ttau(:,ii)*odum ! staggered
  end do
  odum = 1./(1.+((1.+ocneps)*0.5*dt*fv(1:ifull))**2)
  odum = odum*eev(1:ifull)
  dumg = -(1.+ocneps)*0.5*dt*odum
  do ii = 1,wlev
    ccv(1:ifull,ii) = ttav(:,ii)*odum ! staggered
  end do
  ! ice
  ! niu and niv hold the free drift solution (staggered).  Wind stress terms are updated in mlo.f90
  niu(1:ifull) = ttau(:,wlev+1)/(1.+(dt*fu(1:ifull))**2) ! staggered
  niv(1:ifull) = ttav(:,wlev+1)/(1.+(dt*fv(1:ifull))**2)

  do ii = 1,wlev
    ! Create arrays to calcuate u and v at t+1, based on pressure gradient at t+1

    ! u^(t+1) = nu = au^(t*) + bu*dpdxu^(t+1)/wrtrho + cu*dpdyu^(t+1)/wrtrho (staggered)
    ! v^(t+1) = nv = av^(t*) + bv*dpdyv^(t+1)/wrtrho + cv*dpdxv^(t+1)/wrtrho (staggered)

    ! u^(t+1)*(1+(0.5*dt*f)^2) = [u + 0.5*dt*f*v - 0.5*dt/wrtrho*dpdx]^(t*) + 0.5*dt*[v - 0.5*dt*f*u - 0.5*dt/wrtrho*dpdy]^(t*)
    !                            - 0.5*dt/wrtrho*[dpdx + 0.5*dt*f*dpdy]^(t+1)     
    ! v^(t+1)*(1+(0.5*dt*f)^2) = [v - 0.5*dt*f*u - 0.5*dt/wrtrho*dpdy]^(t*) - 0.5*dt*[u + 0.5*dt*f*v - 0.5*dt/wrtrho*dpdx]^(t*)
    !                            - 0.5*dt/wrtrho*[dpdy - 0.5*dt*f*dpdx]^(t+1)
            
    au(:) = ccu(1:ifull,ii)
    bu(:) = dumf(:)
    cu(:) =  (1.+ocneps)*0.5*dt*fu(1:ifull)*bu(:)
  
    av(:) = ccv(1:ifull,ii)
    bv(:) = dumg(:)
    cv(:) = -(1.+ocneps)*0.5*dt*fv(1:ifull)*bv(:)

    ! P = grav*rhobar*(z+eta)
    ! dPdx = grav*sig*(dd+eta)*drhobardx + grav*rhobar*detadx
    
    ! Note pressure gradients are along constant z surfaces
    !p = ps + grav*wrtrho*tt + grav*sig*(eta+dd)*wrtrho
    !dpdxu=dpsdxu+grav*wrtrho*dttdxu+grav*sig*ddu*drhobardxu|neta=0+grav*wrtrho*detadxu !+grav*sig*etau*drhobardxu|neta=0
    !dpdyu=dpsdyu+grav*wrtrho*dttdyu+grav*sig*ddu*drhobardyu|neta=0+grav*wrtrho*detadyu !+grav*sig*etau*drhobardyu|neta=0
    !dpdxv=dpsdxv+grav*wrtrho*dttdxv+grav*sig*ddv*drhobardxv|neta=0+grav*wrtrho*detadxv !+grav*sig*etav*drhobardxv|neta=0
    !dpdyv=dpsdyv+grav*wrtrho*dttdyv+grav*sig*ddv*drhobardyv|neta=0+grav*wrtrho*detadyv !+grav*sig*etav*drhobardyv|neta=0

    ! Create arrays for u and v at t+1 in terms of neta gradients
    
    ! nu = kku + llu*ddu + mmu*detadxu + nnu*detadyu !+llu*etau   (staggered)
    ! nv = kkv + llv*ddv + mmv*detadyv + nnv*detadxv !+llv*etav   (staggered)
    
    kku(:,ii) = au + bu*(dpsdxu/wrtrho+grav*dttdxu) + cu*(dpsdyu/wrtrho+grav*dttdyu)
    llu(:,ii) = grav*gosig(ii)*(bu*drhobardxu(:,ii) + cu*drhobardyu(:,ii))/wrtrho
    mmu(:,ii) = bu*grav
    nnu(:,ii) = cu*grav

    kkv(:,ii) = av + bv*(dpsdyv/wrtrho+grav*dttdyv) + cu*(dpsdxv/wrtrho+grav*dttdxv)
    llv(:,ii) = grav*gosig(ii)*(bv*drhobardyv(:,ii) + cv*drhobardxv(:,ii))/wrtrho
    mmv(:,ii) = bv*grav
    nnv(:,ii) = cv*grav
    
    kku(:,ii) = kku(:,ii)*eeu(1:ifull)
    llu(:,ii) = llu(:,ii)*eeu(1:ifull)
    mmu(:,ii) = mmu(:,ii)*eeu(1:ifull)
    nnu(:,ii) = nnu(:,ii)*eeu(1:ifull)

    kkv(:,ii) = kkv(:,ii)*eev(1:ifull)
    llv(:,ii) = llv(:,ii)*eev(1:ifull)
    mmv(:,ii) = mmv(:,ii)*eev(1:ifull)
    nnv(:,ii) = nnv(:,ii)*eev(1:ifull)

    ! Pre-integrate arrays for u and v at t+1 (i.e., for calculating net divergence at t+1)
    ! We can neglect ppu, qqu, nnu, ppv, qqv and nnv terms as they cancel when calculating
    ! the divergence.
  
    !sum nu dz = sou+spu*ddu+squ*detadxu+ssu*detadyu !+spu*etau
    !sum nv dz = sov+spv*ddv+sqv*detadyv+ssv*detadxv !+spv*etav

    sou(1:ifull) = sou(1:ifull) + kku(:,ii)*godsig(ii)
    spu(1:ifull) = spu(1:ifull) + llu(:,ii)*godsig(ii)
    squ(1:ifull) = squ(1:ifull) + mmu(:,ii)*godsig(ii)
    ssu(1:ifull) = ssu(1:ifull) + nnu(:,ii)*godsig(ii)

    sov(1:ifull) = sov(1:ifull) + kkv(:,ii)*godsig(ii)
    spv(1:ifull) = spv(1:ifull) + llv(:,ii)*godsig(ii)
    sqv(1:ifull) = sqv(1:ifull) + mmv(:,ii)*godsig(ii)
    ssv(1:ifull) = ssv(1:ifull) + nnv(:,ii)*godsig(ii)

  end do

  ! calculate terms for ice velocity at t+1

  ! (staggered)
  ! niu(t+1) = niu + ibu*dipdx + icu*dipdy
  ! niv(t+1) = niv + ibv*dipdy + icv*dipdx

  imu=0.5*(imass(1:ifull)+imass_e)
  imv=0.5*(imass(1:ifull)+imass_n)
  ! note missing 0.5 as ip is the average of t and t+1
  ibu(1:ifull)=-dt*eeu(1:ifull)/(imu*(1.+(dt*fu(1:ifull))**2))
  ibv(1:ifull)=-dt*eev(1:ifull)/(imv*(1.+(dt*fv(1:ifull))**2))
  icu(1:ifull)= dt*fu(1:ifull)*ibu(1:ifull)
  icv(1:ifull)=-dt*fv(1:ifull)*ibv(1:ifull)

  ! update boundary
  dumc(1:ifull,1) = sou(1:ifull)
  dumd(1:ifull,1) = sov(1:ifull)
  dumc(1:ifull,2) = spu(1:ifull)
  dumd(1:ifull,2) = spv(1:ifull)
  dumc(1:ifull,3) = squ(1:ifull)
  dumd(1:ifull,3) = sqv(1:ifull)
  dumc(1:ifull,4) = ssu(1:ifull)
  dumd(1:ifull,4) = ssv(1:ifull)
  dumc(1:ifull,5) = ibu(1:ifull)
  dumd(1:ifull,5) = ibv(1:ifull)
  dumc(1:ifull,6) = icu(1:ifull)
  dumd(1:ifull,6) = icv(1:ifull)
  dumc(1:ifull,7) = niu(1:ifull)
  dumd(1:ifull,7) = niv(1:ifull)
  call boundsuv(dumc(:,1:7),dumd(:,1:7),stag=-9) ! stag=-9 updates iwu and isv
  sou(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,1)
  sov(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,1)
  spu(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,2)
  spv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,2)
  squ(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,3)
  sqv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,3)
  ssu(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,4)
  ssv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,4)
  ibu(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,5)
  ibv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,5)
  icu(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,6)
  icv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,6)
  niu(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,7)
  niv(ifull+1:ifull+iextra) = dumd(ifull+1:ifull+iextra,7)

  ! Iteratively solve for free surface height, eta
  ! Iterative loop to estimate ice 'pressure'
  if ( precon<-9999 ) then
    ! Multi-grid
    call mlomg(neta,sou,sov,spu,spv,squ,sqv,ssu,ssv,xps,   &
               ipice,niu,niv,ibu,ibv,icu,icv,              &
               ipmax,totits,maxglobseta,maxglobip,minwater)
  else
    ! Usual SOR
    write(6,*) "ERROR: MLO dynamics requires precon=-10000"
    call ccmpi_abort(-1)
  end if

  dumc(1:ifull,1)=neta(1:ifull)
  dumc(1:ifull,2)=ipice(1:ifull)
  call bounds(dumc(:,1:2),corner=.true.)
  neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
  ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

  call unpack_nsew(neta,neta_n,neta_s,neta_e,neta_w)
  call unpack_nsew(ipice,ipice_n,ipice_s,ipice_e,ipice_w)

!$omp simd
  do iq = 1,ifull
    tnu(iq)=0.5*(neta_n(iq)+neta(ien(iq)))
    tee(iq)=0.5*(neta(iq)+neta_e(iq))
    tsu(iq)=0.5*(neta_s(iq)+neta(ies(iq)))
    tev(iq)=0.5*(neta_e(iq)+neta(ine(iq)))
    tnn(iq)=0.5*(neta(iq)+neta_n(iq))
    twv(iq)=0.5*(neta_w(iq)+neta(inw(iq)))
  end do  
  detadxu=(neta_e-neta(1:ifull))*emu(1:ifull)/ds
  detadyu=0.5*(stwgt(1:ifull,1)*(tnu-tee)+stwgt(1:ifull,2)*(tee-tsu))*emu(1:ifull)/ds
  detadxv=0.5*(stwgt(1:ifull,3)*(tev-tnn)+stwgt(1:ifull,4)*(tnn-twv))*emv(1:ifull)/ds
  detadyv=(neta_n-neta(1:ifull))*emv(1:ifull)/ds

  ! Update currents once neta is calculated
  do ii = 1,wlev
    ! update currents (staggered)
    nu(1:ifull,ii) = kku(:,ii) + llu(:,ii)*ddu(1:ifull) + mmu(:,ii)*detadxu + nnu(:,ii)*detadyu !+llu(:,ii)*oeu
    nv(1:ifull,ii) = kkv(:,ii) + llv(:,ii)*ddv(1:ifull) + mmv(:,ii)*detadyv + nnv(:,ii)*detadxv !+llv(:,ii)*oev 
  end do

  
#ifdef mlodebug
  if ( any(abs(nu(1:ifull,:))>20.) .or. any(abs(nv(1:ifull,:))>20.) ) then
    write(6,*) "ERROR: current out-of-range after solver"
    write(6,*) "u ",minval(nu(1:ifull,:)),maxval(nu(1:ifull,:))
    write(6,*) "v ",minval(nv(1:ifull,:)),maxval(nv(1:ifull,:))
    stop
  end if
#endif

  call START_LOG(wateriadv_begin)

  ! Update ice velocity with internal pressure terms
!$omp simd
  do iq = 1,ifull
    tnu(iq)=0.5*(ipice_n(iq)+ipice(ien(iq)))
    tee(iq)=0.5*(ipice(iq)+ipice_e(iq))
    tsu(iq)=0.5*(ipice_s(iq)+ipice(ies(iq)))
    tev(iq)=0.5*(ipice_e(iq)+ipice(ine(iq)))
    tnn(iq)=0.5*(ipice(iq)+ipice_n(iq))
    twv(iq)=0.5*(ipice_w(iq)+ipice(inw(iq)))
  end do  
  dipdxu=(ipice_e-ipice(1:ifull))*emu(1:ifull)/ds
  dipdyu=0.5*(stwgt(1:ifull,1)*(tnu-tee)+stwgt(1:ifull,2)*(tee-tsu))*emu(1:ifull)/ds
  dipdxv=0.5*(stwgt(1:ifull,3)*(tev-tnn)+stwgt(1:ifull,4)*(tnn-twv))*emv(1:ifull)/ds
  dipdyv=(ipice_n-ipice(1:ifull))*emv(1:ifull)/ds
  niu(1:ifull)=niu(1:ifull)+ibu(1:ifull)*dipdxu+icu(1:ifull)*dipdyu
  niv(1:ifull)=niv(1:ifull)+ibv(1:ifull)*dipdyv+icv(1:ifull)*dipdxv
  niu(1:ifull)=niu(1:ifull)*eeu(1:ifull)
  niv(1:ifull)=niv(1:ifull)*eev(1:ifull)
  call boundsuv(niu,niv,stag=-9)

  ! Normalisation factor for conserving ice flow in and out of gridbox
  call unpack_svwu(niu,niv,ni_isv,ni_iwu)
  spnet(1:ifull)=(-min(ni_iwu*em_iwu,0.)+max(niu(1:ifull)*emu(1:ifull),0.)     &
                  -min(ni_isv*em_isv,0.)+max(niv(1:ifull)*emv(1:ifull),0.))/ds


  ! ADVECT ICE ------------------------------------------------------
  ! use simple upwind scheme

  ! Horizontal advection for ice area
  dumc(1:ifull,1) = fracice/(em(1:ifull)*em(1:ifull))             ! dumc(:,1) is an area
  ! Horizontal advection for ice volume
  dumc(1:ifull,2) = sicedep*fracice/(em(1:ifull)*em(1:ifull))     ! dumc(:,2) is a volume
  ! Horizontal advection for snow volume
  dumc(1:ifull,3) = snowd*0.001*fracice/(em(1:ifull)*em(1:ifull)) ! dumc(:,3) is a volume
  ! Horizontal advection for ice energy store
  dumc(1:ifull,4) = i_sto*fracice/(em(1:ifull)*em(1:ifull))
  ndsn(1:ifull) = snowd(1:ifull)*0.001
  call mloexpgamm(gamm,sicedep,ndsn(1:ifull),0)
  ! Horizontal advection for surface temperature
  dumc(1:ifull,5) = i_it(1:ifull,1)*fracice*gamm(:,1)/(em(1:ifull)*em(1:ifull))
  ! Horizontal advection of snow temperature
  dumc(1:ifull,6) = i_it(1:ifull,2)*fracice*gamm(:,2)/(em(1:ifull)*em(1:ifull))
  ! Horizontal advection of ice temperatures
  dumc(1:ifull,7) = i_it(1:ifull,3)*fracice*gamm(:,3)/(em(1:ifull)*em(1:ifull))
  dumc(1:ifull,8) = i_it(1:ifull,4)*fracice*gamm(:,3)/(em(1:ifull)*em(1:ifull)) 
  ! Conservation
  dumc(1:ifull,9) = spnet(1:ifull)
  call bounds(dumc(:,1:9))
  spnet(ifull+1:ifull+iextra) = dumc(ifull+1:ifull+iextra,9)
  do ii = 1,8
    call upwind_iceadv(dumc(:,ii),niu,niv,spnet)
  end do  
  nfracice(1:ifull) = min( max( dumc(1:ifull,1)*em(1:ifull)*em(1:ifull), 0. ), maxicefrac )
  ndic(1:ifull) = dumc(1:ifull,2)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
  ndsn(1:ifull) = dumc(1:ifull,3)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
  nsto(1:ifull) = dumc(1:ifull,4)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
  call mloexpgamm(gamm,ndic,ndsn,0)
  nit(1:ifull,1) = dumc(1:ifull,5)*em(1:ifull)*em(1:ifull)/max(gamm(:,1)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,2) = dumc(1:ifull,6)*em(1:ifull)*em(1:ifull)/max(gamm(:,2)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,3) = dumc(1:ifull,7)*em(1:ifull)*em(1:ifull)/max(gamm(:,3)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,4) = dumc(1:ifull,8)*em(1:ifull)*em(1:ifull)/max(gamm(:,3)*nfracice(1:ifull),1.E-10)

  ! populate grid points that have no sea ice
  where ( nfracice(1:ifull)<1.E-4 .or. ndic(1:ifull)<1.E-4 )
    nfracice(1:ifull) = 0.
    ndic(1:ifull) = 0.
    ndsn(1:ifull) = 0.
    nsto(1:ifull) = 0.
    nit(1:ifull,1) = 273.16
    nit(1:ifull,2) = 273.16
    nit(1:ifull,3) = 273.16
    nit(1:ifull,4) = 273.16
  elsewhere ( ndsn(1:ifull)<1.E-4 )
    ndsn(1:ifull) = 0.
    nit(1:ifull,2) = 273.16
  end where


  ! ocean
  ! use nu and nv for unstagged currents
  tau(:,1:wlev)=nu(1:ifull,:)
  tav(:,1:wlev)=nv(1:ifull,:)
  ! ice
  ! unstagger ice velocities
  tau(:,wlev+1)=niu(1:ifull)
  tav(:,wlev+1)=niv(1:ifull)
  call mlounstaguv(tau,tav,ttau,ttav)
  ! ocean
  nu(1:ifull,:)=ttau(:,1:wlev)
  nv(1:ifull,:)=ttav(:,1:wlev)
  ! ice
  niu(1:ifull)=ttau(:,wlev+1)
  niv(1:ifull)=ttav(:,wlev+1)

  call END_LOG(wateriadv_end)

  call START_LOG(watermfix_begin)

  ! volume conservation for water
  if ( nud_sfh==0 ) then
    ! this is a mfix=1 method.  mfix=3 may not be viable for neta  
    delpos(1) = 0.
    delneg(1) = 0.
    neta(1:ifull) = max(min(neta(1:ifull), 130.), -130.)
    odum = (neta(1:ifull)-w_e)*ee(1:ifull)
    !odum = neta(1:ifull)*ee(1:ifull)
    call ccglobal_posneg(odum,delpos(1),delneg(1))
    alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-30))
    if ( abs(alph_p)>1.e-20 ) then
      neta(1:ifull) = w_e(1:ifull) + max(0.,odum)*alph_p + min(0.,odum)/alph_p
    end if
  end if

  ! temperature conservation (usually off when nudging SSTs)
  if ( nud_sst==0 ) then
    delpos(1) = 0.
    delneg(1) = 0.
    do ii = 1,wlev
      nt(1:ifull,ii) = min( max( nt(1:ifull,ii), 150.-wrtemp ), 375.-wrtemp )
      where( wtr(1:ifull) )
        mfixdum(:,ii,1)=(nt(1:ifull,ii)-w_t(:,ii))*dd(1:ifull) 
      elsewhere
        mfixdum(:,ii,1)=0.
      end where
    end do
    neg_godsig(1:wlev) = -godsig(1:wlev)
    call ccglobal_posneg(mfixdum(:,:,1),delpos(1),delneg(1),neg_godsig)
    alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
    do ii = 1,wlev
      where(wtr(1:ifull) .and. abs(alph_p)>1.e-20)
        nt(1:ifull,ii)=w_t(1:ifull,ii)                                                  &
                       +(max(0.,mfixdum(:,ii,1))*alph_p+min(0.,mfixdum(:,ii,1))/alph_p) &
                       /dd(1:ifull)
      elsewhere
        nt(1:ifull,ii) = w_t(:,ii)              
      end where
    end do
    
#ifdef mlodebug
    if ( any( nt(1:ifull,:)+wrtemp<100. .or. nt(1:ifull,:)+wrtemp>400. ) ) then
      write(6,*) "ERROR: nt is out of range after conservation fix"
      write(6,*) "minval,maxval ",minval(nt(1:ifull,:)),maxval(nt(1:ifull,:))
      write(6,*) "minloc,maxloc ",minloc(nt(1:ifull,:)),maxloc(nt(1:ifull,:))
      call ccmpi_abort(-1)
    end if
#endif
    
  end if

  ! salinity conservation
  if ( nud_sss==0 ) then
    delpos(1) = 0.
    delneg(1) = 0.
    ndum(:) = 0.
    do ii = 1,wlev
      ndum(:) = ndum(:) + w_s(1:ifull,ii)*godsig(ii)
    end do
    do ii = 1,wlev
      ns(1:ifull,ii) = max( ns(1:ifull,ii), 0. )  
      where( wtr(1:ifull) .and. ndum>0. )
        mfixdum(:,ii,1) = (ns(1:ifull,ii)-w_s(:,ii))*dd(1:ifull)
      elsewhere
        mfixdum(:,ii,1) = 0.
      end where
    end do
    neg_godsig(1:wlev) = -godsig(1:wlev)
    call ccglobal_posneg(mfixdum(:,:,1),delpos(1),delneg(1),neg_godsig)
    alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-30))
    do ii = 1,wlev
      where( wtr(1:ifull) .and. ndum>0. .and. abs(alph_p)>1.e-20 )
        ns(1:ifull,ii) = w_s(1:ifull,ii)                                                &
                       +(max(0.,mfixdum(:,ii,1))*alph_p+min(0.,mfixdum(:,ii,1))/alph_p) &
                       /dd(1:ifull)
      elsewhere
        ns(1:ifull,ii) = w_s(:,ii)              
      end where
    end do
  end if

  if ( myid==0 .and. (ktau<=5.or.maxglobseta>tol.or.maxglobip>itol) ) then
    write(6,*) "MLODYNAMICS ",totits,maxglobseta,maxglobip
  end if

  call END_LOG(watermfix_end)


  ! reset time-step
  dt = dtin_mlo

end do ! mspec_mlo


#ifdef mlodebug
if ( any(abs(nu(1:ifull,:))>20.) .or. any(abs(nv(1:ifull,:))>20.) ) then
  write(6,*) "ERROR: current out-of-range at end of mlodynamics"
  write(6,*) "u ",minval(nu(1:ifull,:)),maxval(nu(1:ifull,:))
  write(6,*) "v ",minval(nv(1:ifull,:)),maxval(nv(1:ifull,:))
  stop
end if
#endif

! IMPORT WATER AND ICE DATA INTO MLO ------------------------------------------
call mloimport(4,neta(1:ifull),0,0)
call mloimport3d(0,nt,0)
call mloimport3d(1,ns,0)
call mloimport3d(2,nu,0)
call mloimport3d(3,nv,0)
do ii = 1,4
  call mloimpice(nit(1:ifull,ii),ii,0)
end do
call mloimpice(nfracice(1:ifull),5,0)
call mloimpice(ndic(1:ifull),6,0)
call mloimpice(ndsn(1:ifull),7,0)
call mloimpice(nsto(1:ifull),8,0)
call mloimpice(niu(1:ifull),9,0)
call mloimpice(niv(1:ifull),10,0)
where ( wtr(1:ifull) )
  fracice = nfracice(1:ifull)
  sicedep = ndic(1:ifull)
  snowd = ndsn(1:ifull)*1000.
end where

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(dt_in,ubar,vbar,nface,xg,yg,x3d,y3d,z3d,wtr)

use cc_mpi
use const_phys
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m
use vecsuv_m
use xyzinfo_m

implicit none

integer iq,i,j,k,n,nn,idel,jdel,intsch,ii
integer, dimension(ifull,wlev), intent(out) :: nface
real, intent(in) :: dt_in
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real, dimension(ifull,wlev), intent(out) :: xg,yg
real(kind=8), dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc
real, dimension(ifull+iextra,wlev,3) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,3) :: sx
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real xxg,yyg
logical, dimension(ifull+iextra), intent(in) :: wtr

call START_LOG(waterdeps_begin)

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do k = 1,wlev
  uc(:,k) = (ax(1:ifull)*ubar(:,k)+bx(1:ifull)*vbar(:,k))*dt_in/rearth ! unit sphere 
  vc(:,k) = (ay(1:ifull)*ubar(:,k)+by(1:ifull)*vbar(:,k))*dt_in/rearth ! unit sphere 
  wc(:,k) = (az(1:ifull)*ubar(:,k)+bz(1:ifull)*vbar(:,k))*dt_in/rearth ! unit sphere 
  x3d(:,k) = x - uc(:,k) ! 1st guess
  y3d(:,k) = y - vc(:,k)
  z3d(:,k) = z - wc(:,k)
end do

! convert to grid point numbering
call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
! Share off processor departure points.
call deptsync(nface,xg,yg)

intsch = mod(ktau,2)
do k = 1,wlev
  where ( .not.wtr(1:ifull) )
    s(1:ifull,k,1) = 0.
    s(1:ifull,k,2) = 0.
    s(1:ifull,k,3) = 0.
  elsewhere
    s(1:ifull,k,1) = uc(1:ifull,k)
    s(1:ifull,k,2) = vc(1:ifull,k)
    s(1:ifull,k,3) = wc(1:ifull,k)
  end where
end do
s(ifull+1:ifull+iextra,:,:) = 0.

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:3) = reshape( s(1:ipan*jpan*npan,1:wlev,1:3), (/ ipan, jpan, npan, wlev, 3 /) )
  do nn = 1,3
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do
      end do
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lwws(n),k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),       k,nn)
      end do               ! n loop
    end do                 ! k loop
  end do                   ! nn loop
 
! Loop over points that need to be calculated for other processes
  do ii = neighnum,1,-1
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop
  
  call intssync_send(3)

  do nn = 1,3
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg  = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg  = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:3) = reshape( s(1:ipan*jpan*npan,1:wlev,1:3), (/ ipan, jpan, npan, wlev, 3 /) )
  do nn = 1,3
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! For other processes
  do ii = neighnum,1,-1
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(3)

  do nn = 1,3
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k = 1,wlev
  where ( wtr(1:ifull) )
    x3d(:,k) = x - 0.5*(uc(:,k)+s(1:ifull,k,1)) ! n+1 guess
    y3d(:,k) = y - 0.5*(vc(:,k)+s(1:ifull,k,2)) ! n+1 guess
    z3d(:,k) = z - 0.5*(wc(:,k)+s(1:ifull,k,3)) ! n+1 guess
  elsewhere
    x3d(:,k) = x
    y3d(:,k) = y
    z3d(:,k) = z
  end where
end do

call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
!     Share off processor departure points.
call deptsync(nface,xg,yg)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then
 
  ! Loop over points that need to be calculated for other processes
  do ii = neighnum,1,-1
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop
  
  call intssync_send(3)

  do nn = 1,3
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg  = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg  = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  ! For other processes
  do ii = neighnum,1,-1
    do nn = 1,3
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(3)

  do nn = 1,3
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k = 1,wlev
  where (wtr(1:ifull))
    x3d(:,k) = x - 0.5*(uc(:,k)+s(1:ifull,k,1)) ! n+1 guess
    y3d(:,k) = y - 0.5*(vc(:,k)+s(1:ifull,k,2)) ! n+1 guess
    z3d(:,k) = z - 0.5*(wc(:,k)+s(1:ifull,k,3)) ! n+1 guess
  elsewhere
    x3d(:,k) = x
    y3d(:,k) = y
    z3d(:,k) = z
  end where
end do

call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
!     Share off processor departure points.
call deptsync(nface,xg,yg)

call END_LOG(waterdeps_end)

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices
! This code is from depts.f90

pure subroutine mlotoij5(x3d,y3d,z3d,nface,xg,yg)

use bigxy4_m
use cc_mpi
use mlo
use newmpar_m
use parm_m
use parmgeom_m
use xyzinfo_m

implicit none

integer loop,iq,i,j,is,js
integer ii
integer, dimension(ifull,wlev), intent(out) :: nface
real, dimension(ifull,wlev), intent(out) :: xg,yg
real, dimension(ifull) :: xstr,ystr,zstr
real, dimension(ifull) :: denxyz,xd,yd,zd
real, dimension(ifull) :: ri,rj
real(kind=8), dimension(ifull,wlev), intent(inout) :: x3d,y3d,z3d
real(kind=8), dimension(ifull) :: den
real(kind=8) alf,alfonsch
real(kind=8) dxx,dxy,dyx,dyy
integer, parameter :: nmaploop = 3

alf = (1._8-schmidt*schmidt)/(1._8+schmidt*schmidt)
alfonsch = 2._8*schmidt/(1._8+schmidt*schmidt)

do ii = 1,wlev

  !     if necessary, transform (x3d, y3d, z3d) to equivalent
  !     coordinates (xstr, ystr, zstr) on regular gnomonic panels
  den=1._8-alf*z3d(:,ii) ! to force real*8
  xstr=real(x3d(:,ii)*(alfonsch/den))
  ystr=real(y3d(:,ii)*(alfonsch/den))
  zstr=real((z3d(:,ii)-alf)/den)

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1
  denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
  xd=xstr/denxyz
  yd=ystr/denxyz
  zd=zstr/denxyz

  where (abs(xstr-denxyz)<1.E-8)
    nface(:,ii)    =0
    xg(:,ii) =      yd
    yg(:,ii) =      zd
  elsewhere (abs(xstr+denxyz)<1.E-8)
    nface(:,ii)    =3
    xg(:,ii) =     -zd
    yg(:,ii) =     -yd
  elsewhere (abs(zstr-denxyz)<1.E-8)
    nface(:,ii)    =1
    xg(:,ii) =      yd
    yg(:,ii) =     -xd
  elsewhere (abs(zstr+denxyz)<1.E-8)
    nface(:,ii)    =4
    xg(:,ii) =      xd
    yg(:,ii) =     -yd
  elsewhere (abs(ystr-denxyz)<1.E-8)
    nface(:,ii)    =2
    xg(:,ii) =     -zd
    yg(:,ii) =     -xd
  elsewhere
    nface(:,ii)    =5
    xg(:,ii) =      xd
    yg(:,ii) =      zd
  end where

  !     use 4* resolution grid il --> 4*il
  xg(:,ii)=min(max(-.999999,xg(:,ii)),.999999)
  yg(:,ii)=min(max(-.999999,yg(:,ii)),.999999)
  !      first guess for ri, rj and nearest i,j
  ri=1.+(1.+xg(:,ii))*real(2*il_g)
  rj=1.+(1.+yg(:,ii))*real(2*il_g)
  do loop=1,nmaploop
    do iq=1,ifull
      i=nint(ri(iq))
      j=nint(rj(iq))
      is=nint(sign(1.,ri(iq)-real(i)))
      js=nint(sign(1.,rj(iq)-real(j)))
      ! predict new value for ri, rj
      dxx=xx4(i+is,j)-xx4(i,j)
      dyx=xx4(i,j+js)-xx4(i,j)
      dxy=yy4(i+is,j)-yy4(i,j)
      dyy=yy4(i,j+js)-yy4(i,j)       
      den(iq)=dxx*dyy-dyx*dxy
      ri(iq) = real(i) + real(is)*real(((xg(iq,ii)-xx4(i,j))*dyy-(yg(iq,ii)-yy4(i,j))*dyx)/den(iq))
      rj(iq) = real(j) + real(js)*real(((yg(iq,ii)-yy4(i,j))*dxx-(xg(iq,ii)-xx4(i,j))*dxy)/den(iq))
    end do
    ri(1:ifull) = min(ri(1:ifull),1.0+1.999999*real(2*il_g))
    ri(1:ifull) = max(ri(1:ifull),1.0+0.000001*real(2*il_g))
    rj(1:ifull) = min(rj(1:ifull),1.0+1.999999*real(2*il_g))
    rj(1:ifull) = max(rj(1:ifull),1.0+0.000001*real(2*il_g))
  end do  ! loop loop
  !      expect xg, yg to range between .5 and il+.5
  xg(:,ii)=0.25*(ri+3.)-0.5  ! -.5 for stag; back to normal ri, rj defn
  yg(:,ii)=0.25*(rj+3.)-0.5  ! -.5 for stag

end do

return
end subroutine mlotoij5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate depature points for semi-Lagrangian advection
! This code is from ints.f

subroutine mlob2ints_uv(s,nface,xg,yg,wtr)

use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m

implicit none

integer idel,iq,jdel
integer i,j,k,n,intsch
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real xxg,yyg
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
logical, dimension(ifull+iextra), intent(in) :: wtr

call START_LOG(waterints_begin)

ntr=size(s,3)
intsch=mod(ktau,2)

do nn=1,ntr
  do k=1,wlev
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn) = 0.
    end where
  end do
end do
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lwws(n),                  k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),                  k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),                  k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),                  k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),                  k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),                  k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),                  k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),                  k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),  k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do nn = 1,ntr
      do iq=1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k=1,wlev      
      do iq=1,ifull
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-real(idel)
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),  k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! For other processes
  do ii=neighnum,1,-1
    do nn = 1,ntr
      do iq=1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k=1,wlev
      do iq=1,ifull
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-real(idel)
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) )
      s(1:ifull,k,nn)=0.
    end where
  end do
end do

call END_LOG(waterints_end)

return
end subroutine mlob2ints_uv


subroutine mlob2ints(s,nface,xg,yg,wtr)

use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m

implicit none

integer idel,iq,jdel
integer i,j,k,n,intsch
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
integer, dimension(ifull) :: s_count
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(ifull+iextra,wlev,size(s,3)) :: s_old
real, dimension(ifull) :: s_tot, s_test_n, s_test_s, s_test_e, s_test_w
real xxg,yyg
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

call START_LOG(waterints_begin)

ntr=size(s,3)
intsch=mod(ktau,2)

do nn=1,ntr
  do k=1,wlev
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn)=cxx-1. ! missing value flag
    end where
  end do
end do

! fill
do ii = 1,3 ! 3 iterations of fill should be enough
  s_old(1:ifull,:,:) = s(1:ifull,:,:)
  call bounds(s_old)
  do nn = 1,ntr
    do k = 1,wlev
      s_tot(:) = 0.
      s_count(:) = 0
      call unpack_nsew(s_old(:,k,nn),s_test_n,s_test_s,s_test_e,s_test_w)
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test_s(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test_s(:)
        s_count(:) = s_count(:) + 1
      end where
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test_w(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test_w(:)
        s_count(:) = s_count(:) + 1
      end where
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test_e(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test_e(:)
        s_count(:) = s_count(:) + 1
      end where
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test_n(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test_n(:)
        s_count(:) = s_count(:) + 1
      end where
      where ( s_count>0 )
        s(1:ifull,k,nn) = s_tot(1:ifull)/real(s_count(1:ifull))
      end where
    end do
  end do
end do

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev      
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lwws(n),                  k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),                  k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),                  k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),                  k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),                  k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),                  k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),                  k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),                  k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),  k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

  ! Loop over points that need to be calculated for other processes
  do ii = neighnum,1,-1
    do nn = 1,ntr
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k = 1,wlev      
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do     ! iq loop
    end do       ! k loop
  end do         ! nn loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! For other processes
  do ii = neighnum,1,-1
    do nn = 1,ntr
      do iq = 1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k = 1,wlev
      do iq = 1,ifull
        idel = int(xg(iq,k))
        xxg = xg(iq,k) - real(idel)
        jdel = int(yg(iq,k))
        yyg = yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4
      end do
    end do
  end do

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn)=0.
    end where
  end do
end do

call END_LOG(waterints_end)

return
end subroutine mlob2ints

subroutine mlob2ints_bs(s,nface,xg,yg,wtr)

use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m
use parmhor_m

implicit none

integer idel,iq,jdel
integer i,j,k,n,intsch
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
integer, dimension(ifull) :: s_count
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(ifull+iextra,wlev,size(s,3)) :: s_old
real, dimension(ifull,wlev,size(s,3)) :: s_store
real, dimension(ifull) :: s_tot, s_test
real cmax, cmin, xxg, yyg
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

call START_LOG(waterints_begin)

ntr = size(s,3)
intsch = mod(ktau,2)

do nn = 1,ntr
  do k = 1,wlev
    s_store(1:ifull,k,nn) = s(1:ifull,k,nn)
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn) = cxx - 1. ! missing value flag
    end where
  end do
end do

! fill
do ii = 1,3 ! 3 iterations of fill should be enough
  s_old(1:ifull,:,:) = s(1:ifull,:,:)
  call bounds(s_old)
  do nn = 1,ntr
    do k = 1,wlev
      s_tot(:) = 0.
      s_count(:) = 0
      s_test(:) = s_old(is,k,nn)
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test(:)
        s_count(:) = s_count(:) + 1
      end where
      s_test(:) = s_old(iw,k,nn)
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test(:)
        s_count(:) = s_count(:) + 1
      end where
      s_test(:) = s_old(ie,k,nn)
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test(:)
        s_count(:) = s_count(:) + 1
      end where
      s_test(:) = s_old(in,k,nn)
      where ( s_old(1:ifull,k,nn)<cxx .and. s_test(1:ifull)>cxx )
        s_tot(:) = s_tot(:) + s_test(:)
        s_count(:) = s_count(:) + 1
      end where
      where ( s_count>0 )
        s(1:ifull,k,nn) = s_tot(1:ifull)/real(s_count(1:ifull))
      end where
    end do
  end do
end do

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if ( intsch==1 ) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i = 1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lwws(n),                  k,nn)
        sx(0,0,n,k,nn)           = s(iws(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),                  k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),                  k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),                  k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),                  k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),                  k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),                  k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),                  k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(1-ipan+n*ipan*jpan),  k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do nn = 1,ntr
      do iq=1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k=1,wlev      
      do iq=1,ifull
        idel=int(xg(iq,k))
        xxg=xg(iq,k) - real(idel)
        jdel=int(yg(iq,k))
        yyg=yg(iq,k) - real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        rmul_1 = sx(idel,  jdel-1,n,k,nn)*dmul_2 + sx(idel+1,jdel-1,n,k,nn)*dmul_3
        rmul_2 = sx(idel-1,jdel,  n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel,  n,k,nn)*cmul_3 + sx(idel+2,jdel,  n,k,nn)*cmul_4
        rmul_3 = sx(idel-1,jdel+1,n,k,nn)*cmul_1 + sx(idel,  jdel+1,n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+2,jdel+1,n,k,nn)*cmul_4
        rmul_4 = sx(idel,  jdel+2,n,k,nn)*dmul_2 + sx(idel+1,jdel+2,n,k,nn)*dmul_3
        s(iq,k,nn) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do       ! iq loop
    end do         ! k loop
  end do           ! nn loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do nn = 1,ntr
    do k = 1,wlev
      do n = 1,npan
!$omp simd
        do j = 1,jpan
          iq = 1+(j-1)*ipan+(n-1)*ipan*jpan
          sx(0,j,n,k,nn)      = s( iw(iq),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(iq),k,nn)
          iq = j*ipan+(n-1)*ipan*jpan
          sx(ipan+1,j,n,k,nn) = s( ie(iq),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(iq),k,nn)
        end do            ! j loop
!$omp simd
        do i=1,ipan
          iq = i+(n-1)*ipan*jpan
          sx(i,0,n,k,nn)      = s( is(iq),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(iq),k,nn)
          iq = i-ipan+n*ipan*jpan
          sx(i,jpan+1,n,k,nn) = s( in(iq),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(iq),k,nn)
        end do            ! i loop
      end do
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!$omp simd
      do n = 1,npan
        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(1+(n-1)*ipan*jpan),   k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(1-ipan+n*ipan*jpan),  k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(n*ipan*jpan),         k,nn)
      end do           ! n loop
    end do             ! k loop
  end do               ! nn loop

! For other processes
  do ii=neighnum,1,-1
    do nn = 1,ntr
      do iq=1,drlen(ii)
        n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(ii)%a(2,iq))
        xxg = dpoints(ii)%a(2,iq) - real(idel)
        jdel = int(dpoints(ii)%a(3,iq))
        yyg = dpoints(ii)%a(3,iq) - real(jdel)
        k = nint(dpoints(ii)%a(4,iq))
        idel = idel - ioff
        jdel = jdel - joff
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        sextra(ii)%a(iq+(nn-1)*drlen(ii)) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do          ! iq loop
    end do            ! nn loop
  end do              ! ii loop

  call intssync_send(ntr)

  do nn = 1,ntr
    do k=1,wlev
      do iq=1,ifull
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-real(idel)
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-real(jdel)
        idel = min( max( idel - ioff, 0), ipan )
        jdel = min( max( jdel - joff, 0), jpan )
        n = min( max( nface(iq,k) + noff, 1), npan )
        ! bi-cubic interpolation
        cmul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        cmul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul_2 = (1.-yyg)
        dmul_3 = yyg
        emul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        emul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        cmin = min(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        cmax = max(sx(idel,  jdel,n,k,nn),sx(idel+1,jdel,  n,k,nn), &
                   sx(idel,jdel+1,n,k,nn),sx(idel+1,jdel+1,n,k,nn))
        rmul_1 = sx(idel-1,jdel,  n,k,nn)*dmul_2 + sx(idel-1,jdel+1,n,k,nn)*dmul_3
        rmul_2 = sx(idel,  jdel-1,n,k,nn)*cmul_1 + sx(idel,  jdel,  n,k,nn)*cmul_2 + &
                 sx(idel,  jdel+1,n,k,nn)*cmul_3 + sx(idel,  jdel+2,n,k,nn)*cmul_4
        rmul_3 = sx(idel+1,jdel-1,n,k,nn)*cmul_1 + sx(idel+1,jdel,  n,k,nn)*cmul_2 + &
                 sx(idel+1,jdel+1,n,k,nn)*cmul_3 + sx(idel+1,jdel+2,n,k,nn)*cmul_4
        rmul_4 = sx(idel+2,jdel,  n,k,nn)*dmul_2 + sx(idel+2,jdel+1,n,k,nn)*dmul_3
        s(iq,k,nn) = min( max( cmin, &
            rmul_1*emul_1 + rmul_2*emul_2 + rmul_3*emul_3 + rmul_4*emul_4 ), cmax ) ! Bermejo & Staniforth
      end do
    end do
  end do

end if                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn) = s_store(1:ifull,k,nn)
    end where
  end do
end do

call END_LOG(waterints_end)

return
end subroutine mlob2ints_bs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

pure subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use mlo
use xyzinfo_m
use newmpar_m

implicit none

integer k
real, dimension(ifull+iextra,wlev), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real(kind=8), dimension(ifull,size(cou,2)), intent(in) :: x3d,y3d,z3d

do k = 1,wlev
!         cross product n1xn2 into vec1
  vec1x = real(y3d(:,k)*z - y*z3d(:,k))
  vec1y = real(z3d(:,k)*x - z*x3d(:,k))
  vec1z = real(x3d(:,k)*y - x*y3d(:,k))
  denb = vec1x*vec1x + vec1y*vec1y + vec1z*vec1z
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
  where ( denb>1.e-10 )
    vecdot = real(x3d(:,k)*x + y3d(:,k)*y + z3d(:,k)*z)
    vec2x = real(x3d(:,k)*vecdot - x)
    vec2y = real(y3d(:,k)*vecdot - y)
    vec2z = real(z3d(:,k)*vecdot - z)
    vec3x = real(x3d(:,k) - vecdot*x)
    vec3y = real(y3d(:,k) - vecdot*y)
    vec3z = real(z3d(:,k) - vecdot*z)
    vdot1 = (vec1x*cou(1:ifull,k) + vec1y*cov(1:ifull,k) + vec1z*cow(1:ifull,k))/denb
    vdot2 = (vec2x*cou(1:ifull,k) + vec2y*cov(1:ifull,k) + vec2z*cow(1:ifull,k))/denb
    cou(1:ifull,k) = vdot1*vec1x + vdot2*vec3x
    cov(1:ifull,k) = vdot1*vec1y + vdot2*vec3y
    cow(1:ifull,k) = vdot1*vec1z + vdot2*vec3z
  end where
end do ! k

return
end subroutine mlorot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stagger u and v
! Modified to include zero velocity next to land points
subroutine mlostaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m

implicit none

integer k,itn,kx,iq
real, dimension(:,:), intent(in) :: u, v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:), allocatable, save :: dtur,dtvr
real, dimension(ifull,0:3) :: nwtu,nwtv
real, dimension(ifull) :: nud,nvd
real, dimension(ifull) :: v_n, v_s, u_e, u_w
real, dimension(:), allocatable, save :: stul,stvl
real, dimension(:), allocatable, save :: stur,stvr
logical, dimension(ifull) :: euetest,euwtest,evntest,evstest
logical, dimension(ifull) :: euewtest,evnstest,eutest,evtest
logical ltest

call START_LOG(ocnstag_begin)

if ( size(u,1)<ifull .or. size(v,1)<ifull ) then
  write(6,*) "ERROR: Argument is too small for mlostaguv"
  call ccmpi_abort(-1)
end if

if ( size(uout,1)<ifull .or. size(vout,1)<ifull ) then
  write(6,*) "ERROR: Argument is too small for mlostaguv"
  call ccmpi_abort(-1)
end if

if (.not.allocated(wtul)) then
  allocate(wtul(ifull+iextra,0:3),wtvl(ifull+iextra,0:3))
  allocate(wtur(ifull+iextra,0:3),wtvr(ifull+iextra,0:3))
  allocate(dtul(ifull,3),dtvl(ifull,3))
  allocate(dtur(ifull,3),dtvr(ifull,3))
  allocate(stul(ifull),stvl(ifull))
  allocate(stur(ifull),stvr(ifull))
  
  ! assign land arrays
  eutest=eeu(1:ifull)>0.5
  evtest=eev(1:ifull)>0.5
  euetest=eutest      .and.eeu(ieu)>0.5
  euwtest=eeu(iwu)>0.5.and.eutest
  evntest=evtest      .and.eev(inv)>0.5
  evstest=eev(isv)>0.5.and.evtest
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest

  ! assign weights (left)

! |   *   | X E   |  EE  |     unstaggered
! W       * X     E            staggered
  where (euewtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.5
    wtul(1:ifull,2)=-0.1
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.1
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.5
    !ud(1:ifull,k)=uin(ieeu,k)*0.1+uin(ieu,k)+uin(1:ifull,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(ieu,k)*0.5-ua(iwu,k)*0.1

! #   *   | X E   |  EE  |     unstaggered
! 0       * X     E            staggered
  elsewhere (euetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.5
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.1
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.5

! |   *   | X E   #  ##  #     unstaggered
!         * X     0  ##  #     staggered
  elsewhere (euwtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=1./3.

! #   *   |   E   #  ##  #     unstaggered
! #       *       #  ##  #     staggered
  elsewhere (eutest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=0.5
    dtul(1:ifull,3)=0.5

  elsewhere
    wtul(1:ifull,0)=0.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=0.
    dtul(1:ifull,3)=0.
            
  end where
  where (evnstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.5
    wtvl(1:ifull,2)=-0.1
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.1
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.5
  elsewhere (evntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.5
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.1
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.5
  elsewhere (evstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0. 
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=1./3.
  elsewhere (evtest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=0.5
    dtvl(1:ifull,3)=0.5
  elsewhere
    wtvl(1:ifull,0)=0.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=0.
    dtvl(1:ifull,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtul,wtvl,stag=-10)
  where (abs(wtul(ieu,0))>1.E-4.and.abs(wtul(1:ifull,1))>1.E-4)
    stul(1:ifull)=-wtul(1:ifull,1)/wtul(ieu,0)
    nwtu(1:ifull,0)=wtul(1:ifull,0)+stul*wtul(ieu,2)
    nwtu(1:ifull,1)=wtul(1:ifull,1)+stul*wtul(ieu,0)
    nwtu(1:ifull,2)=wtul(1:ifull,2)
    nwtu(1:ifull,3)=wtul(1:ifull,3)-stul*wtul(ieu,1)
  elsewhere
    stul(1:ifull)=0.
    nwtu(1:ifull,0)=wtul(1:ifull,0)
    nwtu(1:ifull,1)=wtul(1:ifull,1)
    nwtu(1:ifull,2)=wtul(1:ifull,2)
    nwtu(1:ifull,3)=wtul(1:ifull,3)
  end where
  where (abs(wtvl(inv,0))>1.E-4.and.abs(wtvl(1:ifull,1))>1.E-4)
    stvl(1:ifull)=-wtvl(1:ifull,1)/wtvl(inv,0)
    nwtv(1:ifull,0)=wtvl(1:ifull,0)+stvl*wtvl(inv,2)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)+stvl*wtvl(inv,0)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)-stvl*wtvl(inv,1)
  elsewhere
    stvl(1:ifull)=0.
    nwtv(1:ifull,0)=wtvl(1:ifull,0)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)
  end where
  where (abs(wtul(1:ifull,0))<1.E-4)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
  elsewhere
    wtul(1:ifull,0)=nwtu(1:ifull,0)
    wtul(1:ifull,1)=nwtu(1:ifull,1)
    wtul(1:ifull,2)=nwtu(1:ifull,2)
    wtul(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvl(1:ifull,0))<1.E-4)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
  elsewhere
    wtvl(1:ifull,0)=nwtv(1:ifull,0)
    wtvl(1:ifull,1)=nwtv(1:ifull,1)
    wtvl(1:ifull,2)=nwtv(1:ifull,2)
    wtvl(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtul,wtvl,stag=-10)
  do k=1,3
    dtul(1:ifull,k)=dtul(1:ifull,k)/wtul(1:ifull,0)
    dtvl(1:ifull,k)=dtvl(1:ifull,k)/wtvl(1:ifull,0)
    wtul(1:ifull,k)=wtul(1:ifull,k)/wtul(1:ifull,0)
    wtvl(1:ifull,k)=wtvl(1:ifull,k)/wtvl(1:ifull,0)
  end do
  stul(1:ifull)=stul(1:ifull)*wtul(ieu,0)/wtul(1:ifull,0)
  stvl(1:ifull)=stvl(1:ifull)*wtvl(inv,0)/wtvl(1:ifull,0)
  wtul(1:ifull,0)=1.
  wtvl(1:ifull,0)=1.

  ! assign weights (right)

! |   W   |   * X |  E   |     unstaggered
!         W     X *      E     staggered
  where (euewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.1
    wtur(1:ifull,2)=-0.5
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.1
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.5

! |   W   |   * X |  E   #     unstaggered
!         W     X *      0     staggered
  elsewhere (euwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=-0.5
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.1
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.5

! #  ##   #   * X |   E   |     unstaggered
! #  ##   0     X *             staggered
  elsewhere (euetest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=1./3.

! #  ##   #   *   |  E   #     unstaggered
! #  ##   #       *      #     staggered
  elsewhere (eutest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=0.5
    dtur(1:ifull,3)=0.5
  elsewhere
    wtur(1:ifull,0)=0.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=0.
    dtur(1:ifull,3)=0.
      
  end where
  where (evnstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.1
    wtvr(1:ifull,2)=-0.5
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.1
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.5
  elsewhere (evstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=-0.5
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.1
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.5
  elsewhere (evntest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=1./3.
  elsewhere (evtest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=0.5
    dtvr(1:ifull,3)=0.5
  elsewhere
    wtvr(1:ifull,0)=0.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=0.
    dtvr(1:ifull,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtur,wtvr,stag=-9)
  where (abs(wtur(iwu,0))>1.E-4.and.abs(wtur(1:ifull,2))>1.E-4)
    stur(1:ifull)=-wtur(1:ifull,2)/wtur(iwu,0)
    nwtu(1:ifull,0)=wtur(1:ifull,0)+stur*wtur(iwu,1)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)+stur*wtur(iwu,0)
    nwtu(1:ifull,3)=wtur(1:ifull,3)-stur*wtur(iwu,2)
  elsewhere
    stur(1:ifull)=0.
    nwtu(1:ifull,0)=wtur(1:ifull,0)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)
  end where
  where (abs(wtvr(isv,0))>1.E-4.and.abs(wtvr(1:ifull,2))>1.E-4)
    stvr(1:ifull)=-wtvr(1:ifull,2)/wtvr(isv,0)
    nwtv(1:ifull,0)=wtvr(1:ifull,0)+stvr*wtvr(isv,1)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)+stvr*wtvr(isv,0)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)-stvr*wtvr(isv,2)
  elsewhere
    stvr(1:ifull)=0.
    nwtv(1:ifull,0)=wtvr(1:ifull,0)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)
  end where
  where (abs(wtur(1:ifull,0))<1.E-4)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
  elsewhere
    wtur(1:ifull,0)=nwtu(1:ifull,0)
    wtur(1:ifull,1)=nwtu(1:ifull,1)
    wtur(1:ifull,2)=nwtu(1:ifull,2)
    wtur(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvr(1:ifull,0))<1.E-4)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
  elsewhere
    wtvr(1:ifull,0)=nwtv(1:ifull,0)
    wtvr(1:ifull,1)=nwtv(1:ifull,1)
    wtvr(1:ifull,2)=nwtv(1:ifull,2)
    wtvr(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtur,wtvr,stag=-9)
  do k=1,3
    dtur(1:ifull,k)=dtur(1:ifull,k)/wtur(1:ifull,0)
    dtvr(1:ifull,k)=dtvr(1:ifull,k)/wtvr(1:ifull,0)
    wtur(1:ifull,k)=wtur(1:ifull,k)/wtur(1:ifull,0)
    wtvr(1:ifull,k)=wtvr(1:ifull,k)/wtvr(1:ifull,0)
  end do
  stur(1:ifull)=stur(1:ifull)*wtur(iwu,0)/wtur(1:ifull,0)
  stvr(1:ifull)=stvr(1:ifull)*wtvr(isv,0)/wtvr(1:ifull,0)
  wtur(1:ifull,0)=1.
  wtvr(1:ifull,0)=1.
  
end if

kx=size(u,2)

do k=1,kx
  uin(1:ifull,k)=u(1:ifull,k)*ee(1:ifull)
  vin(1:ifull,k)=v(1:ifull,k)*ee(1:ifull)
end do

if ( mstagf==0 ) then
  ltest = .true.
else if ( mstagf<0 ) then
  ltest = .false.
else
  ! using ktau-1 ensures that the staggering is relative to the C grid
  ltest = mod(ktau-koff-nstagoffmlo,2*mstagf)<mstagf
end if

if ( ltest ) then

  call boundsuv(uin,vin,stag=1)
  do k=1,kx
!$omp simd
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,1)*uin(ieeu(iq),k)+dtul(iq,2)*uin(ieu(iq),k)+dtul(iq,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,1)*vin(innv(iq),k)+dtvl(iq,2)*vin(inv(iq),k)+dtvl(iq,3)*vin(iq,k)
    end do  
  end do
  
  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    call unpack_nveu(ud(:,k),vd(:,k),v_n,u_e)  
    ! Apply JLM's preconditioner
!$omp simd
    do iq = 1,ifull  
      nud(iq)=ud(iq,k)-stul(iq)*u_e(iq)
      nvd(iq)=vd(iq,k)-stvl(iq)*v_n(iq)
    end do
    ud(1:ifull,k) = nud(1:ifull)
    vd(1:ifull,k) = nvd(1:ifull)
    ua(1:ifull,k) = nud(1:ifull)
    va(1:ifull,k) = nvd(1:ifull)
  end do

  ! There are many ways to handle staggering near coastlines.
  ! This version supports the following properties:
  ! - Staggering is exactly reversible for the staggered (C grid) frame
  ! - Zero flux at land boundaries
  ! - The wave amplitude should be preserved

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)  
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,1)*u_e(iq)+wtul(iq,2)*u_w(iq)+wtul(iq,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,1)*v_n(iq)+wtvl(iq,2)*v_s(iq)+wtvl(iq,3)*va(innv(iq),k)
      end do  
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,1)*u_e(iq)+wtul(iq,2)*u_w(iq)+wtul(iq,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,1)*v_n(iq)+wtvl(iq,2)*v_s(iq)+wtvl(iq,3)*vin(innv(iq),k)
      end do  
    end do
  end do                 ! itn=1,itnmax

else

  call boundsuv(uin,vin)
  do k=1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,1)*u_w(iq)+dtur(iq,2)*uin(iq,k)+dtur(iq,3)*u_e(iq)
      vd(iq,k)=dtvr(iq,1)*v_s(iq)+dtvr(iq,2)*vin(iq,k)+dtvr(iq,3)*v_n(iq)
    end do  
  end do

  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    call unpack_svwu(ud(:,k),vd(:,k),v_s,u_w)  
    ! Apply JLM's preconditioner
!$omp simd
    do iq = 1,ifull  
      nud(iq)=ud(iq,k)-stur(iq)*u_w(iq)
      nvd(iq)=vd(iq,k)-stvr(iq)*v_s(iq)
    end do  
    ud(1:ifull,k) = nud(1:ifull)
    vd(1:ifull,k) = nvd(1:ifull)
    ua(1:ifull,k) = nud(1:ifull)
    va(1:ifull,k) = nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)  
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,1)*u_e(iq)+wtur(iq,2)*u_w(iq)+wtur(iq,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,1)*v_n(iq)+wtvr(iq,2)*v_s(iq)+wtvr(iq,3)*va(issv(iq),k)
      end do  
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,1)*u_e(iq)+wtur(iq,2)*u_w(iq)+wtur(iq,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,1)*v_n(iq)+wtvr(iq,2)*v_s(iq)+wtvr(iq,3)*vin(issv(iq),k)
      end do  
    end do
  end do                 ! itn=1,itnmax

end if

uout(1:ifull,1:kx)=ua(1:ifull,1:kx)
vout(1:ifull,1:kx)=va(1:ifull,1:kx)

call END_LOG(ocnstag_end)

return
end subroutine mlostaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unstagger u and v
! Modified to include zero velocity over land points
subroutine mlounstaguv(u,v,uout,vout,toff)

use cc_mpi
use indices_m
use mlo
use newmpar_m
use parm_m

implicit none

integer, intent(in), optional :: toff
integer k,itn,kx,zoff,iq
real, dimension(:,:), intent(in) :: u,v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(ifull,0:3) :: nwtu,nwtv
real, dimension(ifull) :: nud, nvd
real, dimension(ifull) :: v_n, v_s, u_e, u_w
real, dimension(:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:), allocatable, save :: dtur,dtvr
real, dimension(:), allocatable, save :: stul,stvl
real, dimension(:), allocatable, save :: stur,stvr
logical, dimension(ifull) :: eutest,evtest
logical, dimension(ifull) :: eetest,ewtest,entest,estest
logical, dimension(ifull) :: euetest,euwtest,evntest,evstest
logical, dimension(ifull) :: euewtest,evnstest
logical, dimension(ifull) :: eeetest,ewwtest,enntest,esstest
logical ltest

call START_LOG(ocnstag_begin)

if ( size(u,1)<ifull .or. size(v,1)<ifull ) then
  write(6,*) "ERROR: Argument is too small for mlounstaguv"
  call ccmpi_abort(-1)
end if

if ( size(uout,1)<ifull .or. size(vout,1)<ifull ) then
  write(6,*) "ERROR: Argument is too small for mlounstaguv"
  call ccmpi_abort(-1)
end if

if (.not.allocated(wtul)) then
  allocate(wtul(ifull+iextra,0:3),wtvl(ifull+iextra,0:3))
  allocate(wtur(ifull+iextra,0:3),wtvr(ifull+iextra,0:3))
  allocate(dtul(ifull,3),dtvl(ifull,3))
  allocate(dtur(ifull,3),dtvr(ifull,3))
  allocate(stul(ifull),stvl(ifull))
  allocate(stur(ifull),stvr(ifull))

  ! assign land arrays
  eetest=ee(1:ifull)*ee(ie)>0.5
  ewtest=ee(iw)*ee(1:ifull)>0.5
  entest=ee(1:ifull)*ee(in)>0.5
  estest=ee(is)*ee(1:ifull)>0.5

  eutest=eeu(iwu)>0.5
  evtest=eev(isv)>0.5
  euetest=eutest.and.eeu(1:ifull)>0.5
  evntest=evtest.and.eev(1:ifull)>0.5
  euwtest=eutest.and.eeu(iwwu)>0.5
  evstest=evtest.and.eev(issv)>0.5
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest
  eeetest=eetest.and.ee(iee)>0.5
  enntest=entest.and.ee(inn)>0.5

  ! assign weights (left)

!  |   W   | X *   |  E   |     unstaggered
! WW       W X     *            staggered
  where (euewtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.1
    wtul(1:ifull,2)=-0.5
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.1
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.5
    !ud(1:ifull,k)=uin(iwwu,k)*0.1+uin(iwu,k)+uin(1:ifull,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(ieu,k)*0.1-ua(iwu,k)*0.5

! ##   W   | X *   |  E   |     unstaggered
! #0       W X     *            staggered
  elsewhere (euetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.1
    wtul(1:ifull,2)=-0.5
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.5
      
!  |   W   | X *   #  ##  #     unstaggered
!          W X     0  ##  #     staggered
  elsewhere (euwtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=-1./3.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.

! ##   ##  #   * X |   E   |     unstaggered
! ##   ##  0     X *             staggered
  elsewhere (eeetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-1./3.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=0.
    dtul(1:ifull,3)=1.

! ##   ##  #   *   |   E    #     unstaggered
! ##   ##  #       *        #     staggered
  elsewhere (eetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=0.
    dtul(1:ifull,3)=1.

! ##   W   |   *   #  ##  #     unstaggered
! ##       W       #  ##  #     staggered
  elsewhere (eutest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=1.
    dtul(1:ifull,3)=0.

  elsewhere
    wtul(1:ifull,0)=0.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(1:ifull,1)=0.
    dtul(1:ifull,2)=0.
    dtul(1:ifull,3)=0.
      
  end where
  where (evnstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.1
    wtvl(1:ifull,2)=-0.5
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.1
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.5
    !vd(1:ifull,k)=vin(issv,k)*0.1+vin(isv,k)+vin(1:ifull,k)*0.5
    !vin(1:ifull,k)=vd(:,k)-va(inv,k)*0.1-va(isv,k)*0.5
  elsewhere (evntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.1
    wtvl(1:ifull,2)=-0.5
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.5
  elsewhere (evstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=-1./3.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.
  elsewhere (enntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-1./3.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=0.
    dtvl(1:ifull,3)=1.
  elsewhere (entest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=0.
    dtvl(1:ifull,3)=1.
  elsewhere (evtest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=1.
    dtvl(1:ifull,3)=0.
  elsewhere
    wtvl(1:ifull,0)=0.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(1:ifull,1)=0.
    dtvl(1:ifull,2)=0.
    dtvl(1:ifull,3)=0.
  end where
    
  ! Apply JLM's preconditioner
  call boundsuv(wtul,wtvl,stag=-9)
  where (abs(wtul(iwu,0))>1.E-4.and.abs(wtul(1:ifull,2))>1.E-4)
    stul(1:ifull)=-wtul(1:ifull,2)/wtul(iwu,0)
    nwtu(1:ifull,0)=wtul(1:ifull,0)+stul*wtul(iwu,1)
    nwtu(1:ifull,1)=wtul(1:ifull,1)
    nwtu(1:ifull,2)=wtul(1:ifull,2)+stul*wtul(iwu,0)
    nwtu(1:ifull,3)=wtul(1:ifull,3)-stul*wtul(iwu,2)
  elsewhere
    stul(1:ifull)=0.
    nwtu(1:ifull,0)=wtul(1:ifull,0)
    nwtu(1:ifull,1)=wtul(1:ifull,1)
    nwtu(1:ifull,2)=wtul(1:ifull,2)
    nwtu(1:ifull,3)=wtul(1:ifull,3)
  end where
  where (abs(wtvl(isv,0))>1.E-4.and.abs(wtvl(1:ifull,2))>1.E-4)
    stvl(1:ifull)=-wtvl(1:ifull,2)/wtvl(isv,0)
    nwtv(1:ifull,0)=wtvl(1:ifull,0)+stvl*wtvl(isv,1)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)+stvl*wtvl(isv,0)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)-stvl*wtvl(isv,2)
  elsewhere
    stvl(1:ifull)=0.
    nwtv(1:ifull,0)=wtvl(1:ifull,0)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)
  end where
  where (abs(wtul(1:ifull,0))<1.E-4)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
  elsewhere
    wtul(1:ifull,0)=nwtu(1:ifull,0)
    wtul(1:ifull,1)=nwtu(1:ifull,1)
    wtul(1:ifull,2)=nwtu(1:ifull,2)
    wtul(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvl(1:ifull,0))<1.E-4)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
  elsewhere
    wtvl(1:ifull,0)=nwtv(1:ifull,0)
    wtvl(1:ifull,1)=nwtv(1:ifull,1)
    wtvl(1:ifull,2)=nwtv(1:ifull,2)
    wtvl(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtul,wtvl,stag=-9)
  do k=1,3
    wtul(1:ifull,k)=wtul(1:ifull,k)/wtul(1:ifull,0)
    wtvl(1:ifull,k)=wtvl(1:ifull,k)/wtvl(1:ifull,0)
    dtul(1:ifull,k)=dtul(1:ifull,k)/wtul(1:ifull,0)
    dtvl(1:ifull,k)=dtvl(1:ifull,k)/wtvl(1:ifull,0)
  end do
  stul(1:ifull)=stul(1:ifull)*wtul(iwu,0)/wtul(1:ifull,0)
  stvl(1:ifull)=stvl(1:ifull)*wtvl(isv,0)/wtvl(1:ifull,0)
  wtul(1:ifull,0)=1.
  wtvl(1:ifull,0)=1.

  ! assign land arrays
  eutest=eeu(1:ifull)>0.5
  evtest=eev(1:ifull)>0.5
  euetest=eutest.and.eeu(ieu)>0.5
  evntest=evtest.and.eev(inv)>0.5
  euwtest=eutest.and.eeu(iwu)>0.5
  evstest=evtest.and.eev(isv)>0.5
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest
  ewwtest=ewtest.and.ee(iww)>0.5
  esstest=estest.and.ee(iss)>0.5

  ! assign weights (right) 
  
!  |   W   |   * X |  E   |     unstaggered
!          W     X *      E     staggered
  where (euewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.5
    wtur(1:ifull,2)=-0.1
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.1
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.5
    !ud(1:ifull,k)=uin(ieu,k)*0.1+uin(1:ifull,k)+uin(iwu,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(iwu,k)*0.1-ua(ieu,k)*0.5
        
!  |   W   |   * X |  E   #     unstaggered
!          W     X *      0     staggered
  elsewhere (euwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.5
    wtur(1:ifull,2)=-0.1
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.5
      
!  #   ##  #   * X |   E   |     unstaggered
!  #   ##  0     X *             staggered
  elsewhere (euetest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-1./3.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.

!  |   W   | X *   #  ##  #     unstaggered
!          W X     0  ##  #     staggered
  elsewhere (ewwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=-1./3.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=0.
    dtur(1:ifull,3)=1.

!  #   W   |   *   #  ##  #     unstaggered
!  #       W       #  ##  #     staggered
  elsewhere (ewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=0.
    dtur(1:ifull,3)=1.

!  #   ##  #   *   |  E   #     unstaggered
!  #   ##  #       *      #     staggered
  elsewhere (eutest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=1.
    dtur(1:ifull,3)=0.
  elsewhere
    wtur(1:ifull,0)=0.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(1:ifull,1)=0.
    dtur(1:ifull,2)=0.
    dtur(1:ifull,3)=0.

  end where
  where (evnstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.5
    wtvr(1:ifull,2)=-0.1
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.1
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.5
    !vd(1:ifull,k)=vin(inv,k)*0.1+vin(1:ifull,k)+vin(isv,k)*0.5
    !vin(1:ifull,k)=vd(:,k)-va(isv,k)*0.1-va(inv,k)*0.5
  elsewhere (evstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.5
    wtvr(1:ifull,2)=-0.1
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.5
  elsewhere (evntest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-1./3.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.
  elsewhere (esstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=-1./3.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=0.
    dtvr(1:ifull,3)=1.
  elsewhere (estest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=0.
    dtvr(1:ifull,3)=1.
  elsewhere (evtest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=1.
    dtvr(1:ifull,3)=0.
  elsewhere
    wtvr(1:ifull,0)=0.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(1:ifull,1)=0.
    dtvr(1:ifull,2)=0.
    dtvr(1:ifull,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtur,wtvr,stag=-10)
  where (abs(wtur(ieu,0))>1.E-4.and.abs(wtur(1:ifull,1))>1.E-4)
    stur(1:ifull)=-wtur(1:ifull,1)/wtur(ieu,0)
    nwtu(1:ifull,0)=wtur(1:ifull,0)+stur*wtur(ieu,2)
    nwtu(1:ifull,1)=wtur(1:ifull,1)+stur*wtur(ieu,0)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)-stur*wtur(ieu,1)
  elsewhere
    stur(1:ifull)=0.
    nwtu(1:ifull,0)=wtur(1:ifull,0)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)
  end where
  where (abs(wtvr(inv,0))>1.E-4.and.abs(wtvr(1:ifull,1))>1.E-4)
    stvr(1:ifull)=-wtvr(1:ifull,1)/wtvr(inv,0)
    nwtv(1:ifull,0)=wtvr(1:ifull,0)+stvr*wtvr(inv,2)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)+stvr*wtvr(inv,0)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)-stvr*wtvr(inv,1)
  elsewhere
    stvr(1:ifull)=0.
    nwtv(1:ifull,0)=wtvr(1:ifull,0)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)
  end where
  where (abs(wtur(1:ifull,0))<1.E-4)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
  elsewhere
    wtur(1:ifull,0)=nwtu(1:ifull,0)
    wtur(1:ifull,1)=nwtu(1:ifull,1)
    wtur(1:ifull,2)=nwtu(1:ifull,2)
    wtur(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvr(1:ifull,0))<1.E-4)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
  elsewhere
    wtvr(1:ifull,0)=nwtv(1:ifull,0)
    wtvr(1:ifull,1)=nwtv(1:ifull,1)
    wtvr(1:ifull,2)=nwtv(1:ifull,2)
    wtvr(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtur,wtvr,stag=-10)
  do k=1,3
    wtur(1:ifull,k)=wtur(1:ifull,k)/wtur(1:ifull,0)
    wtvr(1:ifull,k)=wtvr(1:ifull,k)/wtvr(1:ifull,0)
    dtur(1:ifull,k)=dtur(1:ifull,k)/wtur(1:ifull,0)
    dtvr(1:ifull,k)=dtvr(1:ifull,k)/wtvr(1:ifull,0)
  end do
  stur(1:ifull)=stur(1:ifull)*wtur(ieu,0)/wtur(1:ifull,0)
  stvr(1:ifull)=stvr(1:ifull)*wtvr(inv,0)/wtvr(1:ifull,0)
  wtur(1:ifull,0)=1.
  wtvr(1:ifull,0)=1.
 
end if

kx=size(u,2)

zoff=0
if (present(toff)) then
  if (toff==1) zoff=koff
end if

do k=1,kx
  uin(1:ifull,k)=u(1:ifull,k)*eeu(1:ifull)
  vin(1:ifull,k)=v(1:ifull,k)*eev(1:ifull)
end do

if (mstagf==0) then
  ltest = .true.
else if (mstagf<0) then
  ltest = .false.
else
  ltest = mod(ktau-zoff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then
  
  call boundsuv(uin,vin,stag=5)
  do k=1,kx
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,1)*uin(iwwu(iq),k)+dtul(iq,2)*u_w(iq)+dtul(iq,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,1)*vin(issv(iq),k)+dtvl(iq,2)*v_s(iq)+dtvl(iq,3)*vin(iq,k)
    end do  
  end do
  
  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    call unpack_svwu(ud(:,k),vd(:,k),v_s,u_w)  
    ! Apply JLM's preconditioner
!$omp simd
    do iq = 1,ifull  
      nud(iq)=ud(iq,k)-stul(iq)*u_w(iq)
      nvd(iq)=vd(iq,k)-stvl(iq)*v_s(iq)
    end do
    ud(1:ifull,k) = nud(1:ifull)
    vd(1:ifull,k) = nvd(1:ifull)
    ua(1:ifull,k) = nud(1:ifull)
    va(1:ifull,k) = nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)  
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,1)*u_e(iq)+wtul(iq,2)*u_w(iq)+wtul(iq,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,1)*v_n(iq)+wtvl(iq,2)*v_s(iq)+wtvl(iq,3)*va(issv(iq),k)
      end do  
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,1)*u_e(iq)+wtul(iq,2)*u_w(iq)+wtul(iq,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,1)*v_n(iq)+wtvl(iq,2)*v_s(iq)+wtvl(iq,3)*vin(issv(iq),k)
      end do  
    end do

  end do                  ! itn=1,itnmax
  
else

  call boundsuv(uin,vin)
  do k=1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,1)*u_e(iq)+dtur(iq,2)*uin(iq,k)+dtur(iq,3)*u_w(iq)
      vd(iq,k)=dtvr(iq,1)*v_n(iq)+dtvr(iq,2)*vin(iq,k)+dtvr(iq,3)*v_s(iq)
    end do  
  end do

  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    call unpack_nveu(ud(:,k),vd(:,k),v_n,u_e)  
    ! Apply JLM's preconditioner
!$omp simd
    do iq = 1,ifull  
      nud(iq)=ud(iq,k)-stur(iq)*u_e(iq)
      nvd(iq)=vd(iq,k)-stvr(iq)*v_n(iq)
    end do
    ud(1:ifull,k) = nud(1:ifull)
    vd(1:ifull,k) = nvd(1:ifull)
    ua(1:ifull,k) = nud(1:ifull)
    va(1:ifull,k) = nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)  
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,1)*u_e(iq)+wtur(iq,2)*u_w(iq)+wtur(iq,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,1)*v_n(iq)+wtvr(iq,2)*v_s(iq)+wtvr(iq,3)*va(innv(iq),k)
      end do  
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,1)*u_e(iq)+wtur(iq,2)*u_w(iq)+wtur(iq,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,1)*v_n(iq)+wtvr(iq,2)*v_s(iq)+wtvr(iq,3)*vin(innv(iq),k)
      end do  
    end do
  end do                  ! itn=1,itnmax

end if

uout(1:ifull,1:kx)=ua(1:ifull,1:kx)
vout(1:ifull,1:kx)=va(1:ifull,1:kx)

call END_LOG(ocnstag_end)

return
end subroutine mlounstaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs vertical advection based on JLMs TVD scheme

subroutine mlovadv(dtin,ww,uu,vv,ss,tt,mm,depdum,idzdum,wtr,cnum)

use cc_mpi
use cc_omp
use mlo
use newmpar_m

implicit none

integer, intent(in) :: cnum
integer ii,iq,its_g,js,je,tile
integer, dimension(imax) :: its
real, intent(in) :: dtin
real, dimension(imax) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depdum,idzdum
real, dimension(:,:), intent(inout) :: uu,vv,ss,tt,mm
real, dimension(imax,0:wlev) :: ww_l
real, dimension(imax,wlev) :: dat_l, depdum_l, dzdum
logical, dimension(ifull), intent(in) :: wtr

call START_LOG(watervadv_begin)

!$omp parallel private(tile,iq,ii,its,js,je,its_g,ww_l)
!$omp do
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  dzdum    = max(idzdum(js:je,:),1.E-10)
  ww_l     = ww(js:je,:)
  depdum_l = depdum(js:je,:)
  
  ! reduce time step to ensure stability
  dtnew(:)=dtin
  do iq = 1,imax
    if (wtr(iq+js-1)) then
      do ii = 1,wlev-1
        ! this trick works if dzdum(iq,ii)<dzdum(iq,ii+1)
        dtnew(iq)=min(dtnew(iq),0.3*dzdum(iq,ii)/max(abs(ww_l(iq,ii)),1.E-12))
      end do
    end if
  end do
  its(1:imax)=int(dtin/(dtnew(:)+0.01))+1
  its_g=maxval(its(:))
  if (its_g>500) then
    write(6,*) "MLOVERT myid,cnum,its_g",myid,cnum,its_g
  end if
  dtnew(:)=dtin/real(its(:))

  dat_l = uu(js:je,:)
  call mlotvd(its,dtnew,ww_l,dat_l,depdum_l,dzdum)
  uu(js:je,:) = dat_l
  dat_l = vv(js:je,:)
  call mlotvd(its,dtnew,ww_l,dat_l,depdum_l,dzdum)
  vv(js:je,:) = dat_l
  dat_l = ss(js:je,:)-34.72
  call mlotvd(its,dtnew,ww_l,dat_l,depdum_l,dzdum)
  ss(js:je,:)=max(dat_l+34.72,0.)
  dat_l = tt(js:je,:)
  call mlotvd(its,dtnew,ww_l,dat_l,depdum_l,dzdum)
  tt(js:je,:)=max(dat_l,-wrtemp)
  dat_l = mm(js:je,:)
  call mlotvd(its,dtnew,ww_l,dat_l,depdum_l,dzdum)
  mm(js:je,:) = dat_l
  
end do
!$omp end do
!$omp end parallel

call END_LOG(watervadv_end)

return
end subroutine mlovadv

pure subroutine mlotvd(its,dtnew,ww,uu,depadj,dzadj)

use mlo
use newmpar_m

implicit none

integer ii,i,iq,kp,kx,ifin
integer, dimension(:), intent(in) :: its
real, dimension(:), intent(in) :: dtnew
real, dimension(:,0:), intent(in) :: ww
real, dimension(:,:), intent(in) :: depadj,dzadj
real, dimension(:,:), intent(inout) :: uu
real, dimension(size(its,1),wlev-1) :: ff
real, dimension(size(its,1),0:wlev) :: delu
real, dimension(size(its,1)) :: fl,fh,cc,rr

! f=(w*u) at half levels
! du/dt = u*dw/dz-df/dz = -w*du/dz

ifin = size(its,1)

delu(1:ifin,0) = 0.
do ii = 1,wlev-1
  delu(1:ifin,ii) = uu(1:ifin,ii+1) - uu(1:ifin,ii)
end do
delu(1:ifin,wlev) = 0.

! TVD part
do ii = 1,wlev-1
  ! +ve ww is downwards to the ocean floor
  where ( ww(1:ifin,ii)>0.)
    rr(1:ifin) = delu(1:ifin,ii-1)/(delu(1:ifin,ii)+sign(1.e-20,delu(1:ifin,ii)))
    fl(1:ifin) = ww(1:ifin,ii)*uu(1:ifin,ii)
  elsewhere
    rr(1:ifin)=delu(1:ifin,ii+1)/(delu(1:ifin,ii)+sign(1.e-20,delu(1:ifin,ii)))  
    fl(1:ifin)=ww(1:ifin,ii)*uu(1:ifin,ii+1)
  end where
  cc(1:ifin)=max(0.,min(1.,2.*rr(1:ifin)),min(2.,rr(1:ifin))) ! superbee
  fh(1:ifin)=ww(1:ifin,ii)*0.5*(uu(1:ifin,ii)+uu(1:ifin,ii+1))         &
   -0.5*(uu(1:ifin,ii+1)-uu(1:ifin,ii))*ww(1:ifin,ii)**2*dtnew(1:ifin) &
   /max(depadj(1:ifin,ii+1)-depadj(1:ifin,ii),1.E-10)
  ff(1:ifin,ii)=fl(1:ifin)+cc(1:ifin)*(fh(1:ifin)-fl(1:ifin))
  !ff(:,ii)=ww(:,ii)*0.5*(uu(1:ifin,ii)+uu(1:ifin,ii+1)) ! explicit
end do
uu(1:ifin,1)=uu(1:ifin,1)+dtnew*(uu(1:ifin,1)*ww(1:ifin,1)-ff(1:ifin,1))/dzadj(1:ifin,1)
do ii=2,wlev-1
  uu(1:ifin,ii)=uu(1:ifin,ii)+dtnew*(uu(1:ifin,ii)*(ww(1:ifin,ii)-ww(1:ifin,ii-1))   &
                                      -ff(1:ifin,ii)+ff(1:ifin,ii-1))/dzadj(1:ifin,ii)
end do
uu(1:ifin,wlev)=uu(1:ifin,wlev)+dtnew*(-uu(1:ifin,wlev)*ww(1:ifin,wlev-1)+ff(1:ifin,wlev-1))/dzadj(1:ifin,wlev)


do iq=1,ifin
  do i=2,its(iq)

    do ii=1,wlev-1
      delu(iq,ii)=uu(iq,ii+1)-uu(iq,ii)
    end do

    ! TVD part
    do ii=1,wlev-1
      ! +ve ww is downwards to the ocean floor
      kp = nint(sign(1.,ww(iq,ii)))
      kx = ii+(1-kp)/2 !  k for ww +ve,  k+1 for ww -ve
      rr(iq)=delu(iq,ii-kp)/(delu(iq,ii)+sign(1.E-20,delu(iq,ii)))
      fl(iq)=ww(iq,ii)*uu(iq,kx)
      cc(iq)=max(0.,min(1.,2.*rr(iq)),min(2.,rr(iq))) ! superbee
      fh(iq)=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1))          &
        -0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq) &
        /max(depadj(iq,ii+1)-depadj(iq,ii),1.E-10)
      ff(iq,ii)=fl(iq)+cc(iq)*(fh(iq)-fl(iq))
    end do
    uu(iq,1)=uu(iq,1)+dtnew(iq)*(uu(iq,1)*ww(iq,1)-ff(iq,1))/dzadj(iq,1)
    do ii=2,wlev-1
      uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1))-ff(iq,ii)+ff(iq,ii-1))/dzadj(iq,ii)
    end do
    uu(iq,wlev)=uu(iq,wlev)+dtnew(iq)*(-uu(iq,wlev)*ww(iq,wlev-1)+ff(iq,wlev-1))/dzadj(iq,wlev)

  end do
end do

return
end subroutine mlotvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use potential temperature and salinity Jacobians to calculate
! density Jacobian

subroutine tsjacobi(nti,nsi,alphabar,betabar,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m, only : f, emu, emv
use mlo, only : wlev,wrtemp
use newmpar_m
use parm_m, only : ds

implicit none

integer ii, jj
real, dimension(ifull+iextra,wlev), intent(in) :: nti, nsi, alphabar, betabar
real, dimension(ifull,wlev), intent(out) :: drhobardxu, drhobardyu, drhobardxv, drhobardyv
real, dimension(ifull+iextra,wlev,2) :: na
real, dimension(ifull) :: absu, bbsu, absv, bbsv
real, dimension(ifull) :: alphabar_n, alphabar_e, betabar_n, betabar_e
real, dimension(ifull) :: m, a, b, c
real, dimension(ifull) :: y_u, y_nu, y_eu, y_su, y_neu, y_seu
real, dimension(ifull) :: y_v, y_nv, y_ev, y_wv, y_env, y_wnv
real, dimension(ifull,wlev,2) :: dnadxu, dnadxv, dnadyu, dnadyv
real, dimension(2) :: na_max, na_min
real, dimension(ifull) :: dd_in
real, dimension(ifull,-1:1) :: na3
real, dimension(ifull,0:1) :: na2

select case( mlojacobi )
  case(0) ! off
    do ii = 1,wlev
      drhobardxu(1:ifull,ii) = 0.
      drhobardxv(1:ifull,ii) = 0.
      drhobardyu(1:ifull,ii) = 0.
      drhobardyv(1:ifull,ii) = 0.
    end do

  case(1,2,3) 
    na(:,:,1) = min(max(271.-wrtemp,nti),373.-wrtemp)
    na(:,:,2) = min(max(0.,  nsi),50. )-34.72
    select case( mlojacobi )
      case(1) ! non-local - spline  
        call seekdelta(na,dnadxu,dnadyu,dnadxv,dnadyv)
      case(2) ! non-local - linear
        call seekdelta_l(na,dnadxu,dnadyu,dnadxv,dnadyv)
      case(3) ! local - Song 1998
        call seekdelta_song(na,dnadxu,dnadyu,dnadxv,dnadyv)
      case default
        write(6,*) "ERROR: unknown mlojacobi option ",mlojacobi
        stop
    end select  
    do ii = 1,wlev
      call unpack_ne(alphabar(:,ii),alphabar_n,alphabar_e)  
      call unpack_ne(betabar(:,ii),betabar_n,betabar_e)
      absu = 0.5*(alphabar(1:ifull,ii)+alphabar_e)*eeu(1:ifull)
      bbsu = 0.5*(betabar(1:ifull,ii) +betabar_e )*eeu(1:ifull)
      absv = 0.5*(alphabar(1:ifull,ii)+alphabar_n)*eev(1:ifull)
      bbsv = 0.5*(betabar(1:ifull,ii) +betabar_n )*eev(1:ifull)

      ! This relationship neglects compression effects due to neta from the EOS.
      drhobardxu(1:ifull,ii) = -absu*dnadxu(1:ifull,ii,1) + bbsu*dnadxu(1:ifull,ii,2)
      drhobardxv(1:ifull,ii) = -absv*dnadxv(1:ifull,ii,1) + bbsv*dnadxv(1:ifull,ii,2)
      drhobardyu(1:ifull,ii) = -absu*dnadyu(1:ifull,ii,1) + bbsu*dnadyu(1:ifull,ii,2)
      drhobardyv(1:ifull,ii) = -absv*dnadyv(1:ifull,ii,1) + bbsv*dnadyv(1:ifull,ii,2)
    end do
    
  case default
    write(6,*) "ERROR: unknown mlojacobi option ",mlojacobi
    stop

end select
    
return
end subroutine tsjacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using an interpolation method - spline

subroutine seekdelta(rhobar,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m
use mlo, only : wlev
use newmpar_m
use parm_m

implicit none

integer ii, jj, iq
real, dimension(ifull,wlev) :: ddux,ddvy
real, dimension(ifull,wlev) :: ddi,dde,ddw,ddn,dds,dden,ddes,ddne,ddnw
real, dimension(ifull,wlev) :: ramp_a,ramp_c,ramp_e
real, dimension(ifull+iextra,wlev) :: dd_i
real, dimension(ifull,wlev,2) :: ri,re,rw,rn,rs,ren,res,rne,rnw
real, dimension(ifull,wlev,2) :: ssi,sse,ssw,ssn,sss,ssen,sses,ssne,ssnw
real, dimension(ifull,wlev,2) :: y2i,y2e,y2w,y2n,y2s,y2en,y2es,y2ne,y2nw
real, dimension(ifull+iextra,wlev,2) :: y2_i
real, dimension(ifull+iextra,wlev,2), intent (in) :: rhobar
real, dimension(ifull,wlev,2), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv

! Here we calculate the slow contribution of the pressure gradient

! dP/dx = g rhobar dneta/dx + g sigma D drhobar/dx + g sigma neta drhobar/dx
!                   (fast)               (slow)          (mixed - neglect?)

! rhobar = int_0^sigma rho dsigma / sigma

! MJT notes - this version fades out extrapolated gradients using ramp_a, etc.
!
! Idealy, we want to separate the neta contribution to drhobar/dx so that it
! can be included in the implicit solution to neta.


do ii = 1,wlev
  dd_i(:,ii) = gosig(ii)*dd(:)
  ddux(:,ii) = gosig(ii)*ddu(1:ifull)
  ddvy(:,ii) = gosig(ii)*ddv(1:ifull)
end do
call mlospline(dd_i,rhobar,y2_i) ! cubic spline

do jj = 1,2
  do ii = 1,wlev
    ssi(:,ii,jj)=rhobar(1:ifull,ii,jj)
    call unpack_ne(rhobar(:,ii,jj),ssn(:,ii,jj),sse(:,ii,jj))
    y2i(:,ii,jj)=y2_i(1:ifull,ii,jj)
    call unpack_ne(y2_i(:,ii,jj),y2n(:,ii,jj),y2e(:,ii,jj))
  end do
end do  
do ii = 1,wlev
  ddi(:,ii)  =dd_i(1:ifull,ii)
  call unpack_ne(dd_i(:,ii),ddn(:,ii),dde(:,ii))
end do  

do jj = 1,2
  do ii = 1,wlev
!$omp simd
    do iq = 1,ifull  
      sss(iq,ii,jj) =rhobar(is(iq),ii,jj)
      ssen(iq,ii,jj)=rhobar(ien(iq),ii,jj)
      sses(iq,ii,jj)=rhobar(ies(iq),ii,jj)
      y2s(iq,ii,jj) =y2_i(is(iq),ii,jj)
      y2en(iq,ii,jj)=y2_i(ien(iq),ii,jj)
      y2es(iq,ii,jj)=y2_i(ies(iq),ii,jj)
    end do  
  end do
end do  
do ii = 1,wlev
!$omp simd
  do iq = 1,ifull
    dds(iq,ii)   =dd_i(is(iq),ii)
    dden(iq,ii)  =dd_i(ien(iq),ii)
    ddes(iq,ii)  =dd_i(ies(iq),ii)
  end do  
end do  
  
! process staggered u locations
ramp_a(:,:)=1.
call seekval(ri,ssi,ddi,ddux,y2i,ramp_a)
call seekval(re,sse,dde,ddux,y2e,ramp_a)
do jj=1,2
  do ii=1,wlev
    drhobardxu(:,ii,jj)=ramp_a(:,ii)*(re(:,ii,jj)-ri(:,ii,jj))*eeu(1:ifull)*emu(1:ifull)/ds
  end do
end do
ramp_c(:,:)=1.
ramp_e(:,:)=1.
call seekval(rn, ssn, ddn, ddux,y2n, ramp_c)
call seekval(ren,ssen,dden,ddux,y2en,ramp_c)
call seekval(rs, sss, dds, ddux,y2s, ramp_e)
call seekval(res,sses,ddes,ddux,y2es,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardyu(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,1)*(rn(:,ii,jj)+ren(:,ii,jj)      &
                                                         -ri(:,ii,jj)-re(:,ii,jj))*emu(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,2)*(ri(:,ii,jj)+re(:,ii,jj)       &
                                                         -rs(:,ii,jj)-res(:,ii,jj))*emu(1:ifull)/ds)
  end do
end do

do jj = 1,2
  do ii = 1,wlev
!$omp simd
    do iq = 1,ifull  
      ssw(iq,ii,jj) =rhobar(iw(iq),ii,jj)
      ssne(iq,ii,jj)=rhobar(ine(iq),ii,jj)
      ssnw(iq,ii,jj)=rhobar(inw(iq),ii,jj)
      y2w(iq,ii,jj) =y2_i(iw(iq),ii,jj)
      y2ne(iq,ii,jj)=y2_i(ine(iq),ii,jj)
      y2nw(iq,ii,jj)=y2_i(inw(iq),ii,jj)
    end do
  end do
end do 
do ii = 1,wlev
!$omp simd
  do iq = 1,ifull  
    ddw(iq,ii) =dd_i(iw(iq),ii)
    ddne(iq,ii)=dd_i(ine(iq),ii)
    ddnw(iq,ii)=dd_i(inw(iq),ii)
  end do
end do  

! now process staggered v locations
ramp_a(:,:)=1.
call seekval(ri,ssi,ddi,ddvy,y2i,ramp_a)
call seekval(rn,ssn,ddn,ddvy,y2n,ramp_a)
do jj=1,2
  do ii=1,wlev
    drhobardyv(:,ii,jj)=ramp_a(:,ii)*(rn(:,ii,jj)-ri(:,ii,jj))*eev(1:ifull)*emv(1:ifull)/ds
  end do
end do
ramp_c(:,:)=1.
ramp_e(:,:)=1.
call seekval(re, sse, dde, ddvy,y2e, ramp_c)
call seekval(rne,ssne,ddne,ddvy,y2ne,ramp_c)
call seekval(rw, ssw, ddw, ddvy,y2w, ramp_e)
call seekval(rnw,ssnw,ddnw,ddvy,y2nw,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardxv(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,3)*(re(:,ii,jj)+rne(:,ii,jj)      &
                                                         -ri(:,ii,jj)-rn(:,ii,jj))*emv(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,4)*(ri(:,ii,jj)+rn(:,ii,jj)       &
                                                         -rw(:,ii,jj)-rnw(:,ii,jj))*emv(1:ifull)/ds)
  end do
end do

return
end subroutine seekdelta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate to common depths - spline

pure subroutine seekval(rout,ssin,ddin,ddseek,y2,ramp)

use cc_mpi
use mlo, only : wlev
use newmpar_m

implicit none

integer iq, ii, jj, ii_min, ii_max
integer, dimension(1) :: pos
integer, dimension(ifull,wlev) :: sindx
real, dimension(ifull,wlev), intent(in) :: ddseek
real, dimension(ifull,wlev), intent(in) :: ddin
real, dimension(ifull,wlev), intent(inout) :: ramp
real, dimension(ifull,wlev,2), intent(in) :: ssin, y2
real, dimension(ifull,wlev,2), intent(out) :: rout
real, dimension(ifull,2) :: ssunpack1, ssunpack0, y2unpack1, y2unpack0
real, dimension(ifull) :: ddunpack1, ddunpack0
real, dimension(ifull) :: h, a, b, tempa, tempb, temph
real, parameter :: dzramp = 0.01 ! extrapolation limit

sindx(:,:) = wlev
do iq = 1,ifull
  ii = 2
  do jj = 1,wlev
    if ( ddseek(iq,jj)<ddin(iq,wlev-1) .and. ii<wlev ) then
      pos = maxloc( ddin(iq,ii:wlev-1), ddseek(iq,jj)<ddin(iq,ii:wlev-1) )
      sindx(iq,jj) = pos(1) + ii - 1
      ii = sindx(iq,jj)
    else
      exit
    end if
  end do
end do
  
do  jj = 1,wlev
  ! MJT notes - This calculation is slow
  ii_min = minval( sindx(:,jj) )
  ii_max = maxval( sindx(:,jj) )
  do ii = ii_min,ii_max
    where ( ii==sindx(:,jj) )
      ddunpack1(:) = ddin(:,ii)
      ddunpack0(:) = ddin(:,ii-1)
      ssunpack1(:,1) = ssin(:,ii,1)
      ssunpack1(:,2) = ssin(:,ii,2)
      ssunpack0(:,1) = ssin(:,ii-1,1)
      ssunpack0(:,2) = ssin(:,ii-1,2)
      y2unpack1(:,1) = y2(:,ii,1)
      y2unpack1(:,2) = y2(:,ii,2)
      y2unpack0(:,1) = y2(:,ii-1,1)    
      y2unpack0(:,2) = y2(:,ii-1,2)    
    end where
  end do

  h(:) = max(ddunpack1(:)-ddunpack0(:), 1.e-8)
  a(:) = (ddunpack1(:)-ddseek(:,jj))/h(:)
  b(:) = 1. - a(:)
  temph(:) = h(:)**2/6.
  tempa(:) = (a(:)**3-a(:))*temph(:)
  tempb(:) = (b(:)**3-b(:))*temph(:)
  
  rout(:,jj,1) = a(:)*ssunpack0(:,1)+b(:)*ssunpack1(:,1)            & ! linear interpolation
                 +tempa(:)*y2unpack0(:,1)+tempb(:)*y2unpack1(:,1)     ! cubic spline terms
  rout(:,jj,2) = a(:)*ssunpack0(:,2)+b(:)*ssunpack1(:,2)            & ! linear interpolation
                 +tempa(:)*y2unpack0(:,2)+tempb(:)*y2unpack1(:,2)     ! cubic spline terms

  
  ! fade out extrapolation
  ramp(:,jj) = ramp(:,jj)*min(max((a(:)+dzramp)/dzramp,0.),1.)*min(max((b(:)+dzramp)/dzramp,0.),1.)
end do

return
end subroutine seekval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate coeff for cubic spline

pure subroutine mlospline(x,y,y2)

use mlo, only : wlev
use newmpar_m

implicit none

integer ii, jj
real, dimension(ifull+iextra,wlev), intent(in) :: x
real, dimension(ifull+iextra,wlev,2), intent(in) :: y
real, dimension(ifull+iextra,wlev,2), intent(out) :: y2
real, dimension(ifull+iextra,wlev) :: u
real, dimension(ifull+iextra,wlev) :: sig
real, dimension(ifull+iextra) :: p

do ii = 2,wlev-1
  sig(:,ii) = (x(:,ii)-x(:,ii-1))/max(x(:,ii+1)-x(:,ii-1), 1.e-8)
end do

do jj = 1,2
  y2(:,1,jj) = 0.
  u(:,1) = 0.
  do ii = 2,wlev-1
    p(:) = sig(:,ii)*y2(:,ii-1,jj) + 2.
    y2(:,ii,jj) = (sig(:,ii)-1.)/p(:)
    u(:,ii) = (6.*((y(:,ii+1,jj)-y(:,ii,jj))/max(x(:,ii+1)-x(:,ii),1.e-8)-(y(:,ii,jj)-y(:,ii-1,jj)) &
             /max(x(:,ii)-x(:,ii-1),1.e-8))/max(x(:,ii+1)-x(:,ii-1),1.e-8)-sig(:,ii)*u(:,ii-1))/p(:)
  end do
  y2(:,wlev,jj) = 0.    
  do ii = wlev-1,1,-1
    y2(:,ii,jj) = y2(:,ii,jj)*y2(:,ii+1,jj) + u(:,ii)
  end do
end do

return
end subroutine mlospline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using an interpolation method - linear

subroutine seekdelta_l(rhobar,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m
use mlo, only : wlev
use newmpar_m
use parm_m

implicit none

integer ii, jj, iq
real, dimension(ifull,wlev) :: ddux,ddvy
real, dimension(ifull,wlev) :: ddi,dde,ddw,ddn,dds,dden,ddes,ddne,ddnw
real, dimension(ifull,wlev) :: ramp_a,ramp_c,ramp_e
real, dimension(ifull+iextra,wlev) :: dd_i
real, dimension(ifull,wlev,2) :: ri,re,rw,rn,rs,ren,res,rne,rnw
real, dimension(ifull,wlev,2) :: ssi,sse,ssw,ssn,sss,ssen,sses,ssne,ssnw
real, dimension(ifull+iextra,wlev,2), intent (in) :: rhobar
real, dimension(ifull,wlev,2), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv

! Here we calculate the slow contribution of the pressure gradient

! dP/dx = g rhobar dneta/dx + g sigma D drhobar/dx + g sigma neta drhobar/dx
!                   (fast)               (slow)          (mixed - neglected?)

! rhobar = int_0^sigma rho dsigma / sigma

! MJT notes - this version fades out extrapolated gradients using ramp_a, etc.
!
! Idealy, we want to separate the neta contribution to drhobar/dx so that it
! can be included in the implicit solution to neta.


do ii = 1,wlev
  dd_i(:,ii) = gosig(ii)*dd(:)
  ddux(:,ii) = gosig(ii)*ddu(1:ifull)
  ddvy(:,ii) = gosig(ii)*ddv(1:ifull)
end do

do jj = 1,2
  do ii = 1,wlev
    ssi(:,ii,jj)=rhobar(1:ifull,ii,jj)
    call unpack_ne(rhobar(:,ii,jj),ssn(:,ii,jj),sse(:,ii,jj))
  end do
end do  
do ii = 1,wlev
  ddi(:,ii)=dd_i(1:ifull,ii)
  call unpack_ne(dd_i(:,ii),ddn(:,ii),dde(:,ii))
end do  

do jj = 1,2
  do ii = 1,wlev
!$omp simd
    do iq = 1,ifull  
      sss(iq,ii,jj) =rhobar(is(iq),ii,jj)
      ssen(iq,ii,jj)=rhobar(ien(iq),ii,jj)
      sses(iq,ii,jj)=rhobar(ies(iq),ii,jj)
    end do  
  end do
end do  
do ii = 1,wlev
!$omp simd
  do iq = 1,ifull
    dds(iq,ii)   =dd_i(is(iq),ii)
    dden(iq,ii)  =dd_i(ien(iq),ii)
    ddes(iq,ii)  =dd_i(ies(iq),ii)
  end do  
end do  

! process staggered u locations
ramp_a(:,:)=1.
call seekval_l(ri,ssi,ddi,ddux,ramp_a)
call seekval_l(re,sse,dde,ddux,ramp_a)
do jj=1,2
  do ii=1,wlev
    drhobardxu(:,ii,jj)=ramp_a(:,ii)*(re(:,ii,jj)-ri(:,ii,jj))*eeu(1:ifull)*emu(1:ifull)/ds
  end do
end do
ramp_c(:,:)=1.
ramp_e(:,:)=1.
call seekval_l(rn, ssn, ddn, ddux,ramp_c)
call seekval_l(ren,ssen,dden,ddux,ramp_c)
call seekval_l(rs, sss, dds, ddux,ramp_e)
call seekval_l(res,sses,ddes,ddux,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardyu(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,1)*(rn(:,ii,jj)+ren(:,ii,jj)      &
                                                         -ri(:,ii,jj)-re(:,ii,jj))*emu(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,2)*(ri(:,ii,jj)+re(:,ii,jj)       &
                                                         -rs(:,ii,jj)-res(:,ii,jj))*emu(1:ifull)/ds)
  end do
end do

do jj = 1,2
  do ii = 1,wlev
!$omp simd
    do iq = 1,ifull  
      ssw(iq,ii,jj) =rhobar(iw(iq),ii,jj)
      ssne(iq,ii,jj)=rhobar(ine(iq),ii,jj)
      ssnw(iq,ii,jj)=rhobar(inw(iq),ii,jj)
    end do
  end do
end do 
do ii = 1,wlev
!$omp simd
  do iq = 1,ifull  
    ddw(iq,ii) =dd_i(iw(iq),ii)
    ddne(iq,ii)=dd_i(ine(iq),ii)
    ddnw(iq,ii)=dd_i(inw(iq),ii)
  end do
end do  

! now process staggered v locations
ramp_a(:,:)=1.
call seekval_l(ri,ssi,ddi,ddvy,ramp_a)
call seekval_l(rn,ssn,ddn,ddvy,ramp_a)
do jj=1,2
  do ii=1,wlev
    drhobardyv(:,ii,jj)=ramp_a(:,ii)*(rn(:,ii,jj)-ri(:,ii,jj))*eev(1:ifull)*emv(1:ifull)/ds
  end do
end do
ramp_c(:,:)=1.
ramp_e(:,:)=1.
call seekval_l(re, sse, dde, ddvy,ramp_c)
call seekval_l(rne,ssne,ddne,ddvy,ramp_c)
call seekval_l(rw, ssw, ddw, ddvy,ramp_e)
call seekval_l(rnw,ssnw,ddnw,ddvy,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardxv(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,3)*(re(:,ii,jj)+rne(:,ii,jj)      &
                                                         -ri(:,ii,jj)-rn(:,ii,jj))*emv(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,4)*(ri(:,ii,jj)+rn(:,ii,jj)       &
                                                         -rw(:,ii,jj)-rnw(:,ii,jj))*emv(1:ifull)/ds)
  end do
end do

return
end subroutine seekdelta_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate to common depths - linear

pure subroutine seekval_l(rout,ssin,ddin,ddseek,ramp)

use cc_mpi
use mlo, only : wlev
use newmpar_m

implicit none

integer iq, ii, jj, ii_min, ii_max
integer, dimension(1) :: pos
integer, dimension(ifull,wlev) :: sindx
real, dimension(ifull,wlev), intent(in) :: ddseek
real, dimension(ifull,wlev), intent(in) :: ddin
real, dimension(ifull,wlev), intent(inout) :: ramp
real, dimension(ifull,wlev,2), intent(in) :: ssin
real, dimension(ifull,wlev,2), intent(out) :: rout
real, dimension(ifull,2) :: ssunpack1, ssunpack0
real, dimension(ifull) :: ddunpack1, ddunpack0
real, dimension(ifull) :: h, a, b
real, parameter :: dzramp = 0.01 ! extrapolation limit

sindx(:,:) = wlev
do iq = 1,ifull
  ii = 2
  do jj = 1,wlev
    if ( ddseek(iq,jj)<ddin(iq,wlev-1) .and. ii<wlev ) then
      pos = maxloc( ddin(iq,ii:wlev-1), ddseek(iq,jj)<ddin(iq,ii:wlev-1) )
      sindx(iq,jj) = pos(1) + ii - 1
      ii = sindx(iq,jj)
    else
      exit
    end if
  end do
end do
  
do  jj = 1,wlev
  ! MJT notes - This calculation is slow
  ii_min = minval( sindx(:,jj) )
  ii_max = maxval( sindx(:,jj) )
  do ii = ii_min,ii_max
    where ( ii==sindx(:,jj) )
      ddunpack1(:) = ddin(:,ii)
      ddunpack0(:) = ddin(:,ii-1)
      ssunpack1(:,1) = ssin(:,ii,1)
      ssunpack1(:,2) = ssin(:,ii,2)
      ssunpack0(:,1) = ssin(:,ii-1,1)
      ssunpack0(:,2) = ssin(:,ii-1,2)
    end where
  end do

  h(:) = max(ddunpack1(:)-ddunpack0(:), 1.e-8)
  a(:) = (ddunpack1(:)-ddseek(:,jj))/h(:)
  b(:) = 1. - a(:)
  
  rout(:,jj,1) = a(:)*ssunpack0(:,1) + b(:)*ssunpack1(:,1)
  rout(:,jj,2) = a(:)*ssunpack0(:,2) + b(:)*ssunpack1(:,2)
  
  ! fade out extrapolation
  ramp(:,jj) = ramp(:,jj)*min(max((a(:)+dzramp)/dzramp,0.),1.)*min(max((b(:)+dzramp)/dzramp,0.),1.)
end do

return
end subroutine seekval_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using an interpolation method - Song

subroutine seekdelta_song(rhobar,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m
use mlo, only : wlev
use newmpar_m
use parm_m

implicit none

integer ii, jj, kk, iq
integer iip1, iim1, kkpm
real, dimension(ifull+iextra,wlev) :: dd_i
real, dimension(ifull,wlev) :: ddi,dde,ddw,ddn,dds,dden,ddes,ddne,ddnw
real, dimension(ifull,wlev) :: ssi,sse,ssw,ssn,sss,ssen,sses,ssne,ssnw
real, dimension(ifull,2) :: tti,tte,ttw,ttn,tts,tten,ttes,ttne,ttnw
real, dimension(ifull+iextra,wlev,2), intent (in) :: rhobar
real, dimension(ifull,wlev,2), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv

! Here we calculate the slow contribution of the pressure gradient

! dP/dx = g rhobar dneta/dx + g sigma D drhobar/dx + g sigma neta drhobar/dx
!                   (fast)               (slow)          (mixed - neglected?)

! rhobar = int_0^sigma rho dsigma / sigma

do ii = 1,wlev
  dd_i(:,ii) = gosig(ii)*dd(:)
  ddi(:,ii) = dd_i(1:ifull,ii)
  call unpack_nsew(dd_i(:,ii),ddn(:,ii),dds(:,ii),dde(:,ii),ddw(:,ii))
!$omp simd
  do iq = 1,ifull
    dden(iq,ii) = dd_i(ien(iq),ii)
    ddes(iq,ii) = dd_i(ies(iq),ii)
    ddne(iq,ii) = dd_i(ine(iq),ii)
    ddnw(iq,ii) = dd_i(inw(iq),ii)
  end do  
end do

do jj = 1,2
  do ii = 1,wlev
    ssi(:,ii) = rhobar(1:ifull,ii,jj)
    call unpack_nsew(rhobar(:,ii,jj),ssn(:,ii),sss(:,ii),sse(:,ii),ssw(:,ii))
!$omp simd
    do iq = 1,ifull
      ssen(iq,ii) = rhobar(ien(iq),ii,jj)
      sses(iq,ii) = rhobar(ies(iq),ii,jj)
      ssne(iq,ii) = rhobar(ine(iq),ii,jj)
      ssnw(iq,ii) = rhobar(inw(iq),ii,jj)
    end do  
  end do
  
  ! process staggered u locations
  do ii = 1,wlev
    iip1 = min(ii + 1,wlev)
    iim1 = max(ii - 1,1)
    
    drhobardxu(:,ii,jj)=ddsong(sse(:,iip1),sse(:,iim1),ssi(:,iip1),ssi(:,iim1), &
                               dde(:,iip1),dde(:,iim1),ddi(:,iip1),ddi(:,iim1)) &
                        *eeu(1:ifull)*emu(1:ifull)/ds

    do kk = 1,2
      kkpm = min( max( ii-3+2*kk, 1 ), wlev )  
      ttn(:,kk) = ssn(:,kkpm)
      tti(:,kk) = ssi(:,kkpm)
      tten(:,kk) = ssen(:,kkpm)
      tte(:,kk) = sse(:,kkpm)
      tts(:,kk) = sss(:,kkpm)
      ttes(:,kk) = sses(:,kkpm)
    end do  
    drhobardyu(:,ii,jj)=0.25*stwgt(1:ifull,1)*emu(1:ifull)/ds                          &
                        *(ddsong(ttn(:,2),ttn(:,1),tti(:,2),tti(:,1),                  &
                                 ddn(:,iip1),ddn(:,iim1),ddi(:,iip1),ddi(:,iim1))      &
                         +ddsong(tten(:,2),tten(:,1),tte(:,2),tte(:,1),                &
                                 dden(:,iip1),dden(:,iim1),dde(:,iip1),dde(:,iim1)))   &
                       +0.25*stwgt(1:ifull,2)*emu(1:ifull)/ds                          &
                        *(ddsong(tti(:,2),tti(:,1),tts(:,2),tts(:,1),                  &
                                 ddi(:,iip1),ddi(:,iim1),dds(:,iip1),dds(:,iim1))      &
                         +ddsong(tte(:,2),tte(:,1),ttes(:,2),ttes(:,1),                &
                                 dde(:,iip1),dde(:,iim1),ddes(:,iip1),ddes(:,iim1)))

  ! now process staggered v locations
    drhobardyv(:,ii,jj)=ddsong(ssn(:,iip1),ssn(:,iim1),ssi(:,iip1),ssi(:,iim1), &
                               ddn(:,iip1),ddn(:,iim1),ddi(:,iip1),ddi(:,iim1)) &
                        *eev(1:ifull)*emv(1:ifull)/ds
    do kk = 1,2  
      kkpm = min( max( ii-3+2*kk, 1 ), wlev )  
      tte(:,kk) = sse(:,kkpm)
      tti(:,kk) = ssi(:,kkpm)
      ttne(:,kk) = ssne(:,kkpm)
      ttn(:,kk) = ssn(:,kkpm)
      ttw(:,kk) = ssw(:,kkpm)
      ttnw(:,kk) = ssnw(:,kkpm)
    end do  
    drhobardxv(:,ii,jj)=0.25*stwgt(1:ifull,3)*emv(1:ifull)/ds                            &
                        *(ddsong(tte(:,2),tte(:,1),tti(:,2),tti(:,1),                    &
                                 dde(:,iip1),dde(:,iim1),ddi(:,iip1),ddi(:,iim1))        &
                         +ddsong(ttne(:,2),ttne(:,1),ttn(:,2),ttn(:,1),                  &
                                 ddne(:,iip1),ddne(:,iim1),ddn(:,iip1),ddn(:,iim1)))     &
                       +0.25*stwgt(1:ifull,4)*emv(1:ifull)/ds                            &
                        *(ddsong(tti(:,2),tti(:,1),ttw(:,2),ttw(:,1),                    &
                                 ddi(:,iip1),ddi(:,iim1),ddw(:,iip1),ddw(:,iim1))        &
                         +ddsong(ttn(:,2),ttn(:,1),ttnw(:,2),ttnw(:,1),                  &
                                 ddn(:,iip1),ddn(:,iim1),ddnw(:,iip1),ddnw(:,iim1)))       
  end do

end do

return
end subroutine seekdelta_song

pure function ddsong(r1,r2,r3,r4,z1,z2,z3,z4) result(ans)

implicit none

real, dimension(:), intent(in) :: r1, r2, r3, r4
real, dimension(:), intent(in) :: z1, z2, z3, z4
real, dimension(size(z1)) :: ans

where ( abs(z2-z1)<1.e-8 .or. abs(z4-z3)<1.e-8 )
  ans = 0.
elsewhere
  ans = 0.5*((r1+r2-r3-r4)*(z2-z1+z4-z3)-(r2-r1+r4-r3)*(z1+z2-z3-z4))/(z2-z1+z4-z3)
end where

return
end function ddsong

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate tidal potential

pure subroutine mlotide(eta,slon,slat,mtimer,jstart)

use newmpar_m

implicit none

integer, intent(in) :: mtimer,jstart
real, dimension(ifull), intent(out) :: eta
real, dimension(ifull), intent(in) :: slon,slat
real ltime,stime,ctime,mn,sn,pn
real, dimension(ifull) :: sinlat,coslat
real, dimension(4) :: aa,ab,ba,bb
real, parameter :: pi = 3.14159265

! amplitude (mom3)
aa(1)=0.141565 ! K1
aa(2)=0.100661 ! O1
aa(3)=0.046848 ! P1
aa(4)=0.019273 ! Q1
ab(1)=0.242334 ! M2
ab(2)=0.112743 ! S2
ab(3)=0.046397 ! N2
ab(4)=0.030684 ! K2
! Love number (mom3)
ba(1)=0.736
ba(2)=0.695
ba(3)=0.706
ba(4)=0.695
bb(1)=0.693
bb(2)=0.693
bb(3)=0.693
bb(4)=0.693
! Frequency (mom3)
!wa(1)=0.7292117
!wa(2)=0.6750774
!wa(3)=0.7252295
!wa(4)=0.6495854
!wb(1)=1.405189
!wb(2)=1.454441 ! exactly twice per day solar day for S2
!wb(3)=1.378797
!wb(4)=1.458423

stime=real(mtimer)/1440.+real(jstart) ! solar time (days)
ctime=stime/36525. ! century (relative to 12Z 31 Dec 1899)

mn=270.43659+481276.89057*ctime+0.00198*ctime*ctime+0.000002*ctime*ctime*ctime ! moon
sn=279.69668+36000.76892*ctime+0.00030*ctime*ctime                             ! sun
pn=334.32956+4069.03404*ctime-0.01032*ctime*ctime-0.00001*ctime*ctime*ctime    ! lunar perigee
mn=mn*pi/180.
sn=sn*pi/180.
pn=pn*pi/180.
stime=(stime-real(int(stime)))*2.*pi ! solar time (rad)
ltime=stime+sn-mn                    ! lunar time (rad)
! 360*stime=360*ltime-9.3+12.2*day

coslat=cos(slat)**2
sinlat=sin(2.*slat)

eta=ba(1)*aa(1)*sinlat*cos(ltime+mn-slon)          &    ! K1 (note sign change)
   -ba(2)*aa(2)*sinlat*cos(ltime-mn-slon)          &    ! O1
   -ba(3)*aa(3)*sinlat*cos(stime-sn-slon)          &    ! P1
   -ba(4)*aa(4)*sinlat*cos(ltime+pn-slon)          &    ! Q1
   -bb(1)*ab(1)*coslat*cos(2.*ltime-2.*slon)       &    ! M2
   -bb(2)*ab(2)*coslat*cos(2.*stime-2.*slon)       &    ! S2
   -bb(3)*ab(3)*coslat*cos(2.*ltime-mn+pn-2.*slon) &    ! N2
   -bb(4)*ab(4)*coslat*cos(2.*stime+2.*sn-2.*slon)      ! K2

!sinlon=sin(2.*slon)
!coslon=cos(2.*slon)
!eta=0.
!do i=1,4
!  eta=eta-ba(i)*aa(i)*coslat*(cos(wa(i)*stime)*coslon-sin(wa(i)*stime)*sinlon)
!  eta=eta-bb(i)*ab(i)*sinlat*(cos(wb(i)*stime)*coslon-sin(wb(i)*stime)*sinlon)
!end do

return
end subroutine mlotide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test for leap year

pure subroutine mloleap(tyear,ttest)

implicit none

integer, intent(in) :: tyear
logical, intent(out) :: ttest

ttest=.false.
if (mod(tyear,4)==0)   ttest=.true.
if (mod(tyear,100)==0) ttest=.false.
if (mod(tyear,400)==0) ttest=.true.

return
end subroutine mloleap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advect sea-ice using simple upwind scheme

subroutine upwind_iceadv(dumc,niu,niv,spnet)

use indices_m
use map_m
use newmpar_m
use parm_m

implicit none

real, dimension(ifull+iextra), intent(inout) :: dumc
real, dimension(ifull+iextra), intent(in) :: niu,niv,spnet
real(kind=8), dimension(ifull) :: odum,dumd,tdum
real, dimension(ifull) :: spnet_n, spnet_s, spnet_e, spnet_w
real, dimension(ifull) :: dumc_n, dumc_s, dumc_e, dumc_w
real, dimension(ifull) :: em_isv, em_iwu, ni_isv, ni_iwu

dumd=dumc(1:ifull)
call unpack_nsew(spnet,spnet_n,spnet_s,spnet_e,spnet_w)
call unpack_nsew(dumc,dumc_n,dumc_s,dumc_e,dumc_w)
call unpack_svwu(niu,niv,ni_isv,ni_iwu)
call unpack_svwu(emu,emv,em_isv,em_iwu)

odum=0.5*dt*(ni_iwu*(dumc(1:ifull)+dumc_w)    -abs(ni_iwu)*(dumc(1:ifull)-dumc_w)    )*emu(iwu)/ds
where ( spnet_w>1.e-10 )
  tdum=dumc_w*max(ni_iwu*em_iwu,0.)/(ds*spnet_w)    
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=dumc(1:ifull)*min(ni_iwu*em_iwu,0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum

odum=-0.5*dt*(niu(1:ifull)*(dumc(1:ifull)+dumc_e)+abs(niu(1:ifull))*(dumc(1:ifull)-dumc_e))*emu(1:ifull)/ds
where ( spnet_e>1.e-10 )
  tdum=-dumc_e*min(niu(1:ifull)*emu(1:ifull),0.)/(ds*spnet_e)
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=-dumc(1:ifull)*max(niu(1:ifull)*emu(1:ifull),0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum  

odum=0.5*dt*(ni_isv*(dumc(1:ifull)+dumc_s)    -abs(ni_isv)*(dumc(1:ifull)-dumc_s)    )*emv(isv)/ds
where ( spnet_s>1.e-10 )
  tdum=dumc_s*max(ni_isv*em_isv,0.)/(ds*spnet_s)
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=dumc(1:ifull)*min(ni_isv*em_isv,0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum

odum=-0.5*dt*(niv(1:ifull)*(dumc(1:ifull)+dumc_n)+abs(niv(1:ifull))*(dumc(1:ifull)-dumc_n))*emv(1:ifull)/ds
where ( spnet_n>1.e-10 )
  tdum=-dumc_n*min(niv(1:ifull)*emv(1:ifull),0.)/(ds*spnet_n)    
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=-dumc(1:ifull)*max(niv(1:ifull)*emv(1:ifull),0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum

dumc(1:ifull)=real(max(dumd,0._8))
  
return
end subroutine upwind_iceadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use multi-grid to solve for free surface and ice pressure

subroutine mlomg(neta,sou,sov,spu,spv,squ,sqv,ssu,ssv,xps,   &
                 ipice,niu,niv,ibu,ibv,icu,icv,              &
                 ipmax,totits,maxglobseta,maxglobip,minwater)

use helmsolve
use indices_m
use map_m
use newmpar_m
use parm_m

implicit none

integer, intent(out) :: totits
real, intent(out) :: maxglobseta, maxglobip
real, intent(in) :: minwater
real, dimension(ifull), intent(in) :: xps
real, dimension(ifull+iextra), intent(inout) :: neta, ipice
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull+iextra), intent(in) :: sou, sov
real, dimension(ifull+iextra), intent(in) :: spu, spv
real, dimension(ifull+iextra), intent(in) :: squ, sqv
real, dimension(ifull+iextra), intent(in) :: ssu, ssv
real, dimension(ifull+iextra), intent(in) :: niu, niv
real, dimension(ifull+iextra), intent(in) :: ibu, ibv
real, dimension(ifull+iextra), intent(in) :: icu, icv
real, dimension(ifull,2) :: zz, zzn, zzs, zze, zzw
real, dimension(ifull,2) :: zzne, zzen, zznw, zzwn, zzse, zzes, zzsw, zzws
real, dimension(ifull,2) :: rhs
real, dimension(ifull) :: sue, svn, suw, svs
real, dimension(ifull) :: que, qvn, quw, qvs
real, dimension(ifull) :: pdiv_d, odiv_d
real, dimension(ifull) :: hh
real, dimension(ifull) :: so_isv, so_iwu
real, dimension(ifull) :: sp_isv, sp_iwu
real, dimension(ifull) :: sq_isv, sq_iwu
real, dimension(ifull) :: ss_isv, ss_iwu
real, dimension(ifull) :: dd_n, dd_s, dd_e, dd_w
real, dimension(ifull) :: ic_isv, ic_iwu, ib_isv, ib_iwu
real, dimension(ifull) :: ni_isv, ni_iwu, em_isv, em_iwu
real, dimension(ifull) :: dd_iwu, dd_isv
real, dimension(ifull) :: stwgt3_isv, stwgt1_iwu
real, dimension(ifull) :: stwgt4_isv, stwgt2_iwu

call mgsor_init

! ocean

call unpack_svwu(stwgt(:,1),stwgt(:,3),stwgt3_isv,stwgt1_iwu)
call unpack_svwu(stwgt(:,2),stwgt(:,4),stwgt4_isv,stwgt2_iwu)
call unpack_svwu(ddu,ddv,dd_isv,dd_iwu)
call unpack_nsew(dd,dd_n,dd_s,dd_e,dd_w)
call unpack_svwu(emu,emv,em_isv,em_iwu)
call unpack_svwu(sou,sov,so_isv,so_iwu)
call unpack_svwu(spu,spv,sp_isv,sp_iwu)
call unpack_svwu(squ,sqv,sq_isv,sq_iwu)
call unpack_svwu(ssu,ssv,ss_isv,ss_iwu)

! prep ocean gradient terms - constant
odiv_d = (sou(1:ifull)*ddu(1:ifull)/emu(1:ifull)-so_iwu*dd_iwu/em_iwu  &
         +sov(1:ifull)*ddv(1:ifull)/emv(1:ifull)-so_isv*dd_isv/em_isv) &
         *em(1:ifull)**2/ds

! prep ocean gradient terms - dd !(+neta)
pdiv_d = (spu(1:ifull)*ddu(1:ifull)**2/emu(1:ifull)-sp_iwu*dd_iwu**2/em_iwu  &
         +spv(1:ifull)*ddv(1:ifull)**2/emv(1:ifull)-sp_isv*dd_isv**2/em_isv) &
         *em(1:ifull)**2/ds

! prep ocean gradient terms - detadxu and detadyv
que = squ(1:ifull)*ddu(1:ifull)/ds
quw = sq_iwu*dd_iwu/ds
qvn = sqv(1:ifull)*ddv(1:ifull)/ds
qvs = sq_isv*dd_isv/ds

! prep ocean gradient terms - detadyu and detadxv
sue = 0.5*ssu(1:ifull)*ddu(1:ifull)/ds
suw = 0.5*ss_iwu*dd_iwu/ds
svn = 0.5*ssv(1:ifull)*ddv(1:ifull)/ds
svs = 0.5*ss_isv*dd_isv/ds


! zz*(d2neta/dx2+d2neta/dy2) + zz*d2neta/dxdy + hh*neta = rhs

zzn(:,1) = (1.+ocneps)*0.5*dt*qvn*em(1:ifull)**2/ds                                           &
          +(1.+ocneps)*0.5*dt*(stwgt(1:ifull,1)*sue+stwgt1_iwu*suw                            &
                              -stwgt(1:ifull,3)*svn-stwgt(1:ifull,4)*svn)*em(1:ifull)**2/ds
zzs(:,1) = (1.+ocneps)*0.5*dt*qvs*em(1:ifull)**2/ds                                           &
          +(1.+ocneps)*0.5*dt*(stwgt(1:ifull,2)*sue+stwgt2_iwu*suw                            &
                              -stwgt3_isv*svs-stwgt4_isv*svs)*em(1:ifull)**2/ds
zze(:,1) = (1.+ocneps)*0.5*dt*que*em(1:ifull)**2/ds                                           &
          +(1.+ocneps)*0.5*dt*(-stwgt(1:ifull,1)*sue-stwgt(1:ifull,2)*sue                     &
                               +stwgt(1:ifull,3)*svn+stwgt3_isv*svs)*em(1:ifull)**2/ds
zzw(:,1) = (1.+ocneps)*0.5*dt*quw*em(1:ifull)**2/ds                                           &
          +(1.+ocneps)*0.5*dt*(-stwgt1_iwu*suw-stwgt2_iwu*suw                                 &
                               +stwgt(1:ifull,4)*svn+stwgt4_isv*svs)*em(1:ifull)**2/ds
zzne(:,1)= (1.+ocneps)*0.5*dt*(stwgt(1:ifull,3)*svn)*em(1:ifull)**2/ds
zzen(:,1)= (1.+ocneps)*0.5*dt*(stwgt(1:ifull,1)*sue)*em(1:ifull)**2/ds
zzse(:,1)= (1.+ocneps)*0.5*dt*(stwgt3_isv*svs)*em(1:ifull)**2/ds 
zzes(:,1)= (1.+ocneps)*0.5*dt*(stwgt(1:ifull,2)*sue)*em(1:ifull)**2/ds
zzsw(:,1)= (1.+ocneps)*0.5*dt*(stwgt4_isv*svs)*em(1:ifull)**2/ds
zzws(:,1)= (1.+ocneps)*0.5*dt*(stwgt2_iwu*suw)*em(1:ifull)**2/ds
zznw(:,1)= (1.+ocneps)*0.5*dt*(stwgt(1:ifull,4)*svn)*em(1:ifull)**2/ds
zzwn(:,1)= (1.+ocneps)*0.5*dt*(stwgt1_iwu*suw)*em(1:ifull)**2/ds
zz(:,1)  = (1.+ocneps)*0.5*dt*(-qvn-qvs-que-quw)*em(1:ifull)**2/ds                            &
          +(1.+ocneps)*0.5*dt*(-stwgt(1:ifull,1)*sue-stwgt(1:ifull,2)*sue                     &
                               -stwgt1_iwu*suw-stwgt2_iwu*suw                                 &
                               -stwgt(1:ifull,3)*svn-stwgt(1:ifull,4)*svn                     &
                               -stwgt3_isv*svs-stwgt4_isv*svs)*em(1:ifull)**2/ds

hh(:) = 1.

rhs(:,1) = xps(1:ifull) - (1.+ocneps)*0.5*dt*(odiv_d+pdiv_d)


! ice
! zz*(d2ipice/dx2 + d2ipice/dy2) + zz*d2ipice/dxdy = rhs

call unpack_svwu(ibu,ibv,ib_isv,ib_iwu)
call unpack_svwu(icu,icv,ic_isv,ic_iwu)
call unpack_svwu(niu,niv,ni_isv,ni_iwu)

zzn(:,2) = -ibv(1:ifull)/ds                                                                        &
           +(-stwgt(1:ifull,1)*icu-stwgt1_iwu*ic_iwu+stwgt(1:ifull,3)*icv+stwgt(1:ifull,4)*icv)/ds
zzs(:,2) = -ib_isv/ds                                                                              &
           +(-stwgt(1:ifull,2)*icu-stwgt2_iwu*ic_iwu+stwgt3_isv*ic_isv+stwgt4_isv*ic_isv)/ds
zze(:,2) = -ibu(1:ifull)/ds                                                                        &
           +(stwgt(1:ifull,1)*icu+stwgt(1:ifull,2)*icu-stwgt(1:ifull,3)*icv-stwgt3_isv*ic_isv)/ds
zzw(:,2) = -ib_iwu/ds                                                                              &
           +(stwgt1_iwu*ic_iwu+stwgt2_iwu*ic_iwu-stwgt(1:ifull,4)*icv-stwgt4_isv*ic_isv)/ds
zzne(:,2)= (-stwgt(1:ifull,3)*icv)/ds
zzen(:,2)= (-stwgt(1:ifull,1)*icu)/ds
zzse(:,2)= (-stwgt3_isv*ic_isv)/ds 
zzes(:,2)= (-stwgt(1:ifull,2)*icu)/ds
zzsw(:,2)= (-stwgt4_isv*ic_isv)/ds
zzws(:,2)= (-stwgt2_iwu*ic_iwu)/ds
zznw(:,2)= (-stwgt(1:ifull,4)*icv)/ds
zzwn(:,2)= (-stwgt1_iwu*ic_iwu)/ds
zz(:,2)  = (ibu(1:ifull)+ib_iwu+ibv(1:ifull)+ib_isv)/ds                                             &
          +(stwgt(1:ifull,1)*icu+stwgt(1:ifull,2)*icu+stwgt1_iwu*ic_iwu+stwgt2_iwu*ic_iwu           &
           +stwgt(1:ifull,3)*icv+stwgt(1:ifull,4)*icv+stwgt3_isv*ic_isv+stwgt4_isv*ic_isv)/ds

rhs(:,2) = min(niu(1:ifull)/emu(1:ifull)-ni_iwu/em_iwu+niv(1:ifull)/emv(1:ifull)-ni_isv/em_isv,0.)

! Ensure that zero is a valid solution for ice free grid points
where ( zz(:,2)>=0. )
  zz(:,2)  = -dt/(ds*100.) ! 100 is the minimum imass
  zzn(:,2) = 0.
  zzs(:,2) = 0.
  zze(:,2) = 0.
  zzw(:,2) = 0.
  zzne(:,2)= 0.
  zzen(:,2)= 0.
  zzse(:,2)= 0.
  zzes(:,2)= 0.
  zzsw(:,2)= 0.
  zzws(:,2)= 0.
  zznw(:,2)= 0.
  zzwn(:,2)= 0.
  rhs(:,2) = 0.
end where

call mgmlo(neta,ipice,zz,zzn,zzs,zze,zzw,                      &
           zzne,zzen,zzse,zzes,zznw,zzwn,zzsw,zzws,            &
           hh,rhs,tol,itol,totits,maxglobseta,maxglobip,ipmax, &
           ee,dd,minwater)

return
end subroutine mlomg

end module mlodynamics
