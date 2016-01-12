! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! - River routing
! - Ocean dynamics
! - Ice dynamics

! This module also links to mlo.f90 which solves for 1D physics of
! ocean and ice processes.  Currently the code assumes the hydrostatic
! approximation which is reasonably valid to 1km resolution.

! Ocean and sea-ice dynamics are based on the R-grid used by CCAM

module mlodynamics

implicit none

private
public mlodiffusion,mlohadv,mlodyninit,ipice
public oldu1,oldv1,oldu2,oldv2,gosig,gosigh,godsig,ocnsmag,ocneps
public mlodiff
public dd
public nstagoffmlo,mstagf

real, dimension(:), allocatable, save :: ipice,ee,eeu,eev,dd,ddu,ddv,dfdyu,dfdxv
real, dimension(:), allocatable, save :: gosig,gosigh,godsig
real, dimension(:,:), allocatable, save :: oldu1,oldu2,oldv1,oldv2
real, dimension(:,:), allocatable, save :: stwgt
integer, save :: nstagoffmlo
integer, parameter :: usetide   = 1       ! tidal forcing (0=off, 1=on)
integer, parameter :: icemode   = 2       ! ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, parameter :: mstagf    = 30      ! alternating staggering (0=off left, -1=off right, >0 alternating)
integer, parameter :: koff      = 1       ! time split stagger relative to A-grid (koff=0) or C-grid (koff=1)
integer, parameter :: nf        = 2       ! power for horizontal diffusion reduction factor
integer, parameter :: itnmax    = 6       ! number of interations for reversible staggering
integer, parameter :: nxtrrho   = 1       ! Estimate rho at t+1 (0=off, 1=on)
integer, parameter :: usepice   = 0       ! include ice in surface pressure (0=without ice, 1=with ice)
integer, save      :: mlodiff   = 0       ! diffusion (0=all, 1=scalars only)
real, parameter :: rhosn      = 330.      ! density snow (kg m^-3)
real, parameter :: rhoic      = 900.      ! density ice  (kg m^-3)
real, parameter :: grav       = 9.80616   ! gravitational constant (m s^-2)
real, parameter :: delphi     = 150.      ! horizontal diffusion reduction factor gradient
real, save      :: ocnsmag    = 1.        ! horizontal diffusion (2. in Griffies (2000), 1.-1.4 in POM (Mellor 2004))
real, save      :: ocneps     = 0.1       ! semi-implicit off-centring term
real, parameter :: maxicefrac = 0.999     ! maximum ice fraction
real, parameter :: tol        = 5.E-4     ! Tolerance for SOR solver (water)
real, parameter :: itol       = 2.E0      ! Tolerance for SOR solver (ice)


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises mlo dynamical arrays
!
subroutine mlodyninit

use cc_mpi
use indices_m
use map_m
use mlo
use soil_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer ii,iq
real, dimension(ifull,0:wlev) :: dephl
real, dimension(ifull,wlev) :: dep,dz
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn
real, dimension(wlev) :: sig,sigh,dsig
real, dimension(3*wlev) :: dumz,gdumz
logical, dimension(ifull+iextra) :: wtr

! staggering
nstagoffmlo=0

! prep land-sea mask
allocate(ee(ifull+iextra))
allocate(eeu(ifull+iextra),eev(ifull+iextra))
ee=0.
eeu=0.
eev=0. 
where(.not.land)
  ee(1:ifull)=1.
end where

call bounds(ee,nrows=2)
wtr=abs(ee-1.)<0.
where (ee>1.5.or.ee<=0.5)
  ee=0.
end where

eeu(1:ifull)=ee(1:ifull)*ee(ie)
eev(1:ifull)=ee(1:ifull)*ee(in)
call boundsuv(eeu,eev,nrows=2)

! The following is for the in-line MLO ocean model ------------------

! Calculate depth arrays (free suface term is included later)
allocate(dd(ifull+iextra))
allocate(ddu(ifull+iextra),ddv(ifull+iextra))
dd=0.
ddu=0.
ddv=0.
dep=0.
dz=0.
do ii=1,wlev
  call mloexpdep(0,dep(:,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
end do
dephl(:,0)=0.
dephl(:,1)=dz(:,1)
do ii=2,wlev
  dephl(:,ii)=dephl(:,ii-1)+dz(:,ii)
end do
dd(1:ifull)=dephl(:,wlev)
dd(1:ifull)=max(dd(1:ifull),1.E-8)

! update bounds values
call bounds(dd,nrows=2)

! update staggered bounds values
ddu(1:ifull)=0.5*(dd(1:ifull)+dd(ie))
ddv(1:ifull)=0.5*(dd(1:ifull)+dd(in))
where (eeu(1:ifull)==0.)
  ddu(1:ifull)=1.E-8
end where
where (eev(1:ifull)==0.)
  ddv(1:ifull)=1.E-8
end where
call boundsuv(ddu,ddv,nrows=2)

if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then

  ! Precompute weights for calculating staggered gradients
  allocate(stwgt(ifull,4))
  stwgt=0.
  where (wtr(in).and.wtr(ine).and.wtr(ien).and.wtr(ie).and.wtr(1:ifull))
    stwgt(:,1)=1.
  end where
  where (wtr(is).and.wtr(ise).and.wtr(ies).and.wtr(ie).and.wtr(1:ifull))
    stwgt(:,2)=1.
  end where
  stwgt(:,3)=stwgt(:,1)
  where (wtr(iw).and.wtr(iwn).and.wtr(inw).and.wtr(in).and.wtr(1:ifull))
    stwgt(:,4)=1.
  end where

  onedice = 0 ! Turn off 1D ice model

  ! dynamics save arrays
  allocate(oldu1(ifull,wlev),oldv1(ifull,wlev))
  allocate(oldu2(ifull,wlev),oldv2(ifull,wlev))
  allocate(ipice(ifull+iextra))
  oldu1 = 0.
  oldv1 = 0.
  oldu2 = oldu1
  oldv2 = oldv1
  ipice = 0.

  ! special gradient arrays
  allocate(dfdyu(ifull),dfdxv(ifull))
  tnu=0.5*(f(in)+f(ine))
  tee=0.5*(f(1:ifull)+f(ie))
  tsu=0.5*(f(is)+f(ise))
  tev=0.5*(f(ie)+f(ien))
  tnn=0.5*(f(1:ifull)+f(in))
  twv=0.5*(f(iw)+f(iwn))
  dfdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
  dfdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
  
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
      if (abs(sig(ii)*dd(iq)-dep(iq,ii))>0.01.or.    &
          abs(sigh(ii)*dd(iq)-dephl(iq,ii))>0.01.or. &
          abs(dsig(ii)*dd(iq)-dz(iq,ii))>0.01) then
        write(6,*) "ERROR: MLO not configured for sigma levels"
        call ccmpi_abort(-1)
      end if
    end do
  end if
end do
dumz(1:wlev)=sig
dumz(wlev+1:2*wlev)=sigh
dumz(2*wlev+1:3*wlev)=dsig
call ccmpi_allreduce(dumz(1:3*wlev),gdumz(1:3*wlev),"max",comm_world)
gosig=gdumz(1:wlev)
gosigh=gdumz(wlev+1:2*wlev)
godsig=gdumz(2*wlev+1:3*wlev)

return
end subroutine mlodyninit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, if mlohadv is
! not called.  Calling diffusion from mlodadv avoids additional
! unpacking and packing
subroutine mlodiffusion

use mlo

implicit none

include 'newmpar.h'

real, dimension(ifull,wlev) :: u,v,tt,ss
real, dimension(ifull) :: eta

! extract data from MLO
u=0.
v=0.
eta=0.
tt=0.
ss=0.
call mloexport3d(0,tt,0)
call mloexport3d(1,ss,0)
call mloexport3d(2,u,0)
call mloexport3d(3,v,0)
call mloexport(4,eta(1:ifull),0,0)

call mlodiffusion_main(u,v,u,v,eta,tt,ss)

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion_main(uauin,uavin,u,v,etain,tt,ss)

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

integer k
real hdif
real, dimension(ifull,wlev), intent(in) :: uauin,uavin
real, dimension(:,:), intent(in) :: u,v,tt,ss
real, dimension(:), intent(in) :: etain
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull,wlev) :: ft,fs,base,outu,outv
real, dimension(ifull+iextra) :: depadj,eta
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: nu,nv,nw
real, dimension(ifull) :: tx_fact,ty_fact
real, dimension(ifull) :: emi

call START_LOG(waterdiff_begin)

! Define diffusion scale and grid spacing
hdif=dt*(ocnsmag/pi)**2
emi=max(dd(1:ifull)+etain(1:ifull),minwater)/em(1:ifull)

! extract data from MLO
do k=1,wlev
  uau(1:ifull,k)=uauin(:,k)*ee(1:ifull)
  uav(1:ifull,k)=uavin(:,k)*ee(1:ifull)
end do
call boundsuv_allvec(uau,uav)

! calculate diffusion following Smagorinsky
do k=1,wlev
  dudx=0.5*((uau(ieu,k)-uau(1:ifull,k))*emu(1:ifull)*eeu(1:ifull) &
           +(uau(1:ifull,k)-uau(iwu,k))*emu(iwu)*eeu(iwu))/ds
  dudy=0.5*((uau(inu,k)-uau(1:ifull,k))*emv(1:ifull)*eev(1:ifull) &
           +(uau(1:ifull,k)-uau(isu,k))*emv(isv)*eev(isv))/ds
  dvdx=0.5*((uav(iev,k)-uav(1:ifull,k))*emu(1:ifull)*eeu(1:ifull) &
           +(uav(1:ifull,k)-uav(iwv,k))*emu(iwu)*eeu(iwu))/ds
  dvdy=0.5*((uav(inv,k)-uav(1:ifull,k))*emv(1:ifull)*eev(1:ifull) &
           +(uav(1:ifull,k)-uav(isv,k))*emv(isv)*eev(isv))/ds

  t_kh(1:ifull,k)=sqrt((dudx-dvdy)**2+(dudy+dvdx)**2)*hdif*emi
end do
t_kh(1:ifull,wlev+1)=etain(1:ifull)
call bounds(t_kh,nehalf=.true.)
eta(:)=t_kh(:,wlev+1)

! reduce diffusion errors where bathymetry gradients are strong
do k=1,wlev
  depadj=gosig(k)*max(dd+eta,minwater)
  tx_fact=1./(1.+(abs(depadj(ie)-depadj(1:ifull))/delphi)**nf)
  ty_fact=1./(1.+(abs(depadj(in)-depadj(1:ifull))/delphi)**nf)

  xfact(1:ifull,k)=0.5*(t_kh(1:ifull,k)+t_kh(ie,k))*tx_fact*eeu(1:ifull) ! reduction factor
  yfact(1:ifull,k)=0.5*(t_kh(1:ifull,k)+t_kh(in,k))*ty_fact*eev(1:ifull) ! reduction factor
end do
call boundsuv(xfact,yfact,stag=-9)

do k=1,wlev
  base(:,k)=emi+xfact(1:ifull,k)+xfact(iwu,k)+yfact(1:ifull,k)+yfact(isv,k)
end do

if (mlodiff==0) then
  ! Laplacian diffusion terms (closure #1)
  do k=1,wlev
    duma(1:ifull,k,1)=ax(1:ifull)*u(1:ifull,k)+bx(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,2)=ay(1:ifull)*u(1:ifull,k)+by(1:ifull)*v(1:ifull,k)
    duma(1:ifull,k,3)=az(1:ifull)*u(1:ifull,k)+bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(duma(:,:,1:3))
  do k=1,wlev
    nu = ( duma(1:ifull,k,1)*emi +                      &
           xfact(1:ifull,k)*duma(ie,k,1) +              &
           xfact(iwu,k)*duma(iw,k,1) +                  &
           yfact(1:ifull,k)*duma(in,k,1) +              &
           yfact(isv,k)*duma(is,k,1) ) / base(:,k)
    nv = ( duma(1:ifull,k,2)*emi +                      &
           xfact(1:ifull,k)*duma(ie,k,2) +              &
           xfact(iwu,k)*duma(iw,k,2) +                  &
           yfact(1:ifull,k)*duma(in,k,2) +              &
           yfact(isv,k)*duma(is,k,2) ) / base(:,k)
    nw = ( duma(1:ifull,k,3)*emi +                      &
           xfact(1:ifull,k)*duma(ie,k,3) +              &
           xfact(iwu,k)*duma(iw,k,3) +                  &
           yfact(1:ifull,k)*duma(in,k,3) +              &
           yfact(isv,k)*duma(is,k,3) ) / base(:,k)
    outu(1:ifull,k)=ax(1:ifull)*nu+ay(1:ifull)*nv+az(1:ifull)*nw
    outv(1:ifull,k)=bx(1:ifull)*nu+by(1:ifull)*nv+bz(1:ifull)*nw
  end do

!  Laplacian diffusion and viscosity terms (closure #2)
!  duma(1:ifull,:,1)=u(1:ifull,:)
!  duma(1:ifull,:,2)=v(1:ifull,:)
!  call boundsuv_allvec(duma(:,:,1),duma(:,:,2))
!  do k=1,wlev
!    outu(:,k)=(duma(1:ifull,k,1)*emi+2.*xfact(1:ifull,k)*duma(ieu,k,1)+2.*xfact(iwu,k)*duma(iwu,k,1) &
!      +yfact(1:ifull,k)*duma(inu,k,1)+yfact(isv,k)*duma(isu,k,1)                                     &
!      +(yfact(1:ifull,k)-yfact(isv,k))*0.5*(duma(iev,k,2)-duma(iwv,k,2))                             &
!      +t_kh(1:ifull,k)*0.5*(duma(inv,k,2)+duma(iev,k,2)-duma(isv,k,2)-duma(iwv,k,2)))                &
!      /(emi+2.*xfact(1:ifull,k)+2.*xfact(iwu,k)+yfact(1:ifull,k)+yfact(isv,k))
!    outv(:,k)=(duma(1:ifull,k,2)*emi+2.*yfact(1:ifull,k)*duma(inv,k,2)+2.*yfact(isv,k)*duma(isv,k,2) &
!      +xfact(1:ifull,k)*duma(iev,k,2)+xfact(iwu,k)*duma(iwv,k,2)                                     &
!      +(xfact(1:ifull,k)-xfact(iwu,k))*0.5*(duma(inu,k,1)-duma(isu,k,1))                             &
!      +t_kh(1:ifull,k)*0.5*(duma(inu,k,1)+duma(ieu,k,1)-duma(isu,k,1)-duma(iwu,k,1)))                &
!      /(emi+2.*yfact(1:ifull,k)+2.*yfact(isv,k)+xfact(1:ifull,k)+xfact(iwu,k))
!  end do

else
  outu(1:ifull,:)=u(1:ifull,:)
  outv(1:ifull,:)=v(1:ifull,:)
end if
  
! Potential temperature and salinity
duma(1:ifull,:,1)=tt(1:ifull,:)
duma(1:ifull,:,2)=ss(1:ifull,:)-34.72
call bounds(duma(:,:,1:2))
do k=1,wlev
  ft(:,k) = ( duma(1:ifull,k,1)*emi +                      &
              xfact(1:ifull,k)*duma(ie,k,1) +              &
              xfact(iwu,k)*duma(iw,k,1) +                  &
              yfact(1:ifull,k)*duma(in,k,1) +              &
              yfact(isv,k)*duma(is,k,1) ) / base(:,k)
  fs(:,k) = ( duma(1:ifull,k,2)*emi +                      &
              xfact(1:ifull,k)*duma(ie,k,2) +              &
              xfact(iwu,k)*duma(iw,k,2) +                  &
              yfact(1:ifull,k)*duma(in,k,2) +              &
              yfact(isv,k)*duma(is,k,2) ) / base(:,k)
end do
fs=max(fs+34.72,0.)

call mloimport3d(0,ft,0)
call mloimport3d(1,fs,0)
call mloimport3d(2,outu,0)
call mloimport3d(3,outv,0)

call END_LOG(waterdiff_end)

return
end subroutine mlodiffusion_main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements some basic hydrostatic dynamics for the
! ocean and ice.  The ocean component employs the R-grid design used
! in CCAM semi-Lagragain dynamics, but uses sigma z vertical
! coordinates.  The staggered/unstaggered pivoting has been modified
! for the land-sea boundary.  Sea-ice is advected using an upwind
! scheme.  Internal sea-ice pressure follows a cavitating fluid
! apporximation.
subroutine mlohadv

use arrays_m
use cc_mpi
use infile
use indices_m
use latlong_m
use map_m
use mlo
use nharrs_m, only : lrestart
use soil_m
use soilsnow_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'parmdyn.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer ii,totits
integer jyear,jmonth,jday,jhour,jmin,mins
integer tyear,jstart
integer itc
integer, dimension(ifull,wlev) :: nface
real maxglobseta,maxglobip
real delpos,delneg,alph_p,fjd
real, dimension(ifull+iextra) :: neta,pice,imass,xodum
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv,nis
real, dimension(ifull+iextra) :: snu,sou,spu,squ,ssu,snv,sov,spv,sqv,ssv
real, dimension(ifull+iextra) :: ibu,ibv,icu,icv,idu,idv,spnet,oeu,oev,tide
real, dimension(ifull+iextra) :: ipmax
real, dimension(ifull) :: i_u,i_v,i_sto,i_sal,rhobaru,rhobarv,ndum
real, dimension(ifull) :: pdiv,qdiv,sdiv,odiv,w_e
real, dimension(ifull) :: pdivb,qdivb,sdivb,odivb,xps
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn,rhou,rhov
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: dipdxu,dipdyu,dipdxv,dipdyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv,odum
real, dimension(ifull) :: imu,imv
real, dimension(ifull) :: sue,suw,svn,svs
real, dimension(ifull) :: pue,puw,pvn,pvs
real, dimension(ifull) :: que,quw,qvn,qvs
real, dimension(ifull) :: piceu,picev,tideu,tidev,ipiceu,ipicev
real, dimension(ifull) :: dumf,dumg
real, dimension(ifull+iextra,wlev,3) :: cou
real, dimension(ifull+iextra,wlev+1) :: eou,eov
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns,mps,dzdum_rho
real, dimension(ifull+iextra,wlev) :: rhobar,rho,dalpha,dbeta
real, dimension(ifull+iextra,wlev) :: ccu,ccv
real, dimension(ifull+iextra,10) :: dumc,dumd
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull,wlev+1) :: tau,tav,ttau,ttav
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,dum
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu,oou,ppu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv,oov,ppv
real, dimension(ifull,wlev) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull,wlev) :: depdum,dzdum
real, dimension(ifull,0:wlev) :: nw
real, dimension(ifull,4) :: i_it
real, dimension(ifull,3) :: gamm
real(kind=8), dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical lleap

! new z levels for including free surface eta (effectively sigma-depth levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! newdz=olddz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth
! depth below free suface is newz+eta=oldz*(1+eta/maxdepth)

! We use a modified form of JLM's trick where:
!   dP/dx+f dP/dy = dP/dx + d(fP)/dy - P df/dy
!   dP/dy-f dP/dx = dP/dy - d(fP)/dy + P df/dx
! which produces a 5-point stencil when solving for
! neta and ip.

! Instead of packing ocean points, we keep the original
! grid and define a land/sea mask.  This allows us to
! use the same index and map factor arrays from the
! atmospheric dynamical core.

call START_LOG(watermisc_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Start"
end if
#endif

! Define land/sea mask
wtr=ee>0.5

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
i_sal = 0.            ! ice salinity
rho   = 0.            ! water density
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

! IMPORT WATER AND ICE DATA -----------------------------------------
#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Import"
end if
#endif
call mloexport3d(0,w_t,0)
call mloexport3d(1,w_s,0)
call mloexport3d(2,w_u,0)
call mloexport3d(3,w_v,0)
call mloexport(4,w_e,0,0)
do ii=1,4
  call mloexpice(i_it(:,ii),ii,0)
end do
call mloexpice(fracice,5,0)
call mloexpice(sicedep,6,0)
call mloexpice(snowd,7,0)
call mloexpice(i_sto,8,0)
call mloexpice(i_u,9,0)
call mloexpice(i_v,10,0)
call mloexpice(i_sal,11,0)

! estimate tidal forcing (assumes leap days)
#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Tides"
end if
#endif
if ( usetide==1 ) then
  call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins,allleap=.true.)
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
nfracice(ifull+1:ifull+iextra)=0.
ndic(ifull+1:ifull+iextra)=0.
ndsn(ifull+1:ifull+iextra)=0.
where (wtr(1:ifull))
  nfracice(1:ifull)=fracice
  ndic(1:ifull)=sicedep
  ndsn(1:ifull)=snowd*0.001
elsewhere
  nfracice(1:ifull)=0.
  ndic(1:ifull)=0.
  ndsn(1:ifull)=0.
end where
nit(1:ifull,:)=i_it
nsto(1:ifull)=i_sto
niu(1:ifull)=i_u
niv(1:ifull)=i_v
nis(1:ifull)=i_sal

! surface pressure and ice mass
! (assume ice velocity is 'slow' compared to 'fast' change in neta)
imass(1:ifull)=ndic(1:ifull)*rhoic+ndsn(1:ifull)*rhosn  ! ice mass per unit area (kg/m^2), unstaggered at time t
pice(1:ifull)=ps(1:ifull)                               ! pressure due to atmosphere at top of water column (unstaggered at t)
if ( usepice==1 ) then
  pice(1:ifull) = pice(1:ifull) + grav*nfracice(1:ifull)*imass(1:ifull) ! include ice in surface pressure
end if

! Limit minimum ice mass for ice velocity calculation.  Hence, we can solve for the ice velocity at
! grid points where the ice is not yet present.  
imass(1:ifull)=max(imass(1:ifull),10.)

! maximum pressure for cavitating fluid
select case(icemode)
  case(2)
    ! cavitating fluid
    ipmax(1:ifull)=27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))*ee(1:ifull)
  case(1)
    ! incompressible fluid
    ipmax(1:ifull)=9.E9*ee(1:ifull)
  case DEFAULT
    ! free drift
    ipmax(1:ifull)=0.
end select

! update scalar bounds and various gradients (aggregate fields for MPI)
dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=pice(1:ifull)
dumc(1:ifull,3)=tide(1:ifull)
dumc(1:ifull,4)=imass(1:ifull)
dumc(1:ifull,5)=ipmax(1:ifull)
call bounds(dumc(:,1:5),corner=.true.)
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
pice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
tide(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,3)
imass(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,4)
ipmax(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,5)
! surface pressure
tnu=0.5*(pice(in)*f(in)+pice(ine)*f(ine))
tee=0.5*(pice(1:ifull)*f(1:ifull)+pice(ie)*f(ie))
tsu=0.5*(pice(is)*f(is)+pice(ise)*f(ise))
tev=0.5*(pice(ie)*f(ie)+pice(ien)*f(ien))
tnn=0.5*(pice(1:ifull)*f(1:ifull)+pice(in)*f(in))
twv=0.5*(pice(iw)*f(iw)+pice(iwn)*f(iwn))
dpsdxu=(pice(ie)-pice(1:ifull))*emu(1:ifull)/ds ! staggered at time t
dpsdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dpsdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dpsdyv=(pice(in)-pice(1:ifull))*emv(1:ifull)/ds
piceu=0.5*(pice(1:ifull)+pice(ie))
picev=0.5*(pice(1:ifull)+pice(in))
! tides
tnu=0.5*(tide(in)*f(in)+tide(ine)*f(ine))
tee=0.5*(tide(1:ifull)*f(1:ifull)+tide(ie)*f(ie))
tsu=0.5*(tide(is)*f(is)+tide(ise)*f(ise))
tev=0.5*(tide(ie)*f(ie)+tide(ien)*f(ien))
tnn=0.5*(tide(1:ifull)*f(1:ifull)+tide(in)*f(in))
twv=0.5*(tide(iw)*f(iw)+tide(iwn)*f(iwn))
dttdxu=(tide(ie)-tide(1:ifull))*emu(1:ifull)/ds ! staggered
dttdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dttdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dttdyv=(tide(in)-tide(1:ifull))*emv(1:ifull)/ds
tideu=0.5*(tide(1:ifull)+tide(ie))
tidev=0.5*(tide(1:ifull)+tide(in))

call END_LOG(watermisc_end)

call START_LOG(watereos_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: density EOS 1"
end if
#endif

! Calculate adjusted depths and thicknesses
xodum=max(dd(:)+neta(:),minwater)
do ii=1,wlev
  depdum(:,ii)=gosig(ii)*xodum(1:ifull)
  dzdum(:,ii) =godsig(ii)*xodum(1:ifull)
  dzdum_rho(:,ii)=godsig(ii)*dd(:) ! MJT suggestion to avoid density oscillations  
end do

! Calculate normalised density rhobar (unstaggered at time t)
! (Assume free surface correction is small so that changes in the compression 
! effect due to neta can be neglected.  Consequently, the neta dependence is 
! separable in the iterative loop)
cou(1:ifull,:,1)=nt(1:ifull,:)
cou(1:ifull,:,2)=ns(1:ifull,:)
call bounds(cou(:,:,1:2),corner=.true.)
nt(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,1)
ns(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,2)
call mloexpdensity(rho,dalpha,dbeta,nt,ns,dzdum_rho,pice,0)
rhobar(:,1)=(rho(:,1)-1030.)*godsig(1)
do ii=2,wlev
  rhobar(:,ii)=rhobar(:,ii-1)+(rho(:,ii)-1030.)*godsig(ii)
end do
do ii=1,wlev
  rhobar(:,ii)=rhobar(:,ii)/gosigh(ii)+1030.
end do

! Calculate normalised density gradients
! method 2: Use potential temperature and salinity Jacobians (see Shchepetkin and McWilliams 2003)
call tsjacobi(nt,ns,dalpha,dbeta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
drhobardxu(:,1)=drhobardxu(:,1)*godsig(1)
drhobardxv(:,1)=drhobardxv(:,1)*godsig(1)
drhobardyu(:,1)=drhobardyu(:,1)*godsig(1)
drhobardyv(:,1)=drhobardyv(:,1)*godsig(1)
do ii=2,wlev
  drhobardxu(:,ii)=drhobardxu(:,ii-1)+drhobardxu(:,ii)*godsig(ii)
  drhobardxv(:,ii)=drhobardxv(:,ii-1)+drhobardxv(:,ii)*godsig(ii)
  drhobardyu(:,ii)=drhobardyu(:,ii-1)+drhobardyu(:,ii)*godsig(ii)
  drhobardyv(:,ii)=drhobardyv(:,ii-1)+drhobardyv(:,ii)*godsig(ii)
end do
do ii=1,wlev
  drhobardxu(:,ii)=drhobardxu(:,ii)/gosigh(ii)
  drhobardxv(:,ii)=drhobardxv(:,ii)/gosigh(ii)
  drhobardyu(:,ii)=drhobardyu(:,ii)/gosigh(ii)
  drhobardyv(:,ii)=drhobardyv(:,ii)/gosigh(ii)
end do

call END_LOG(watereos_end)

call START_LOG(waterdeps_begin)

! ADVECT WATER AND ICE ----------------------------------------------
! Water currents are advected using semi-Lagrangian advection
! based on McGregor's CCAM advection routines.
! Velocity is set to zero at ocean boundaries.

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Departure points"
end if
#endif

! save arrays
if (ktau==1.and..not.lrestart) then
  oldu1=nu(1:ifull,:)
  oldv1=nv(1:ifull,:)
  oldu2=nu(1:ifull,:)
  oldv2=nv(1:ifull,:)
  ipice=0. ! ice free drift solution
end if

! Estimate currents at t+1/2 for semi-Lagrangian advection
do ii=1,wlev
  nuh(:,ii)=(15.*nu(1:ifull,ii)-10.*oldu1(:,ii)+3.*oldu2(:,ii))*ee(1:ifull)/8. ! U at t+1/2
  nvh(:,ii)=(15.*nv(1:ifull,ii)-10.*oldv1(:,ii)+3.*oldv2(:,ii))*ee(1:ifull)/8. ! V at t+1/2
end do

! Calculate depature points
call mlodeps(dt,nuh,nvh,nface,xg,yg,x3d,y3d,z3d,wtr)

oldu2=oldu1
oldv2=oldv1
oldu1=nu(1:ifull,:)
oldv1=nv(1:ifull,:)

call END_LOG(waterdeps_end)

call START_LOG(waterhadv_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: continuity equation"
end if
#endif

! Advect continuity equation to tstar
! Calculate velocity on C-grid for consistancy with iterative free surface calculation
! nw is at half levels with nw(:,1) at the bottom of level 1
! positive nw is moving downwards to the ocean floor
! For now, assume Boussinesq fluid and treat density in continuity equation as constant
! true vertical velocity = nw-u*((1-sig)*deta/dx-sig*d(dd)/dx)-v*((1-sig)*deta/dy-sig*d(dd)/dy)-(1-sig)*deta/dt
call mlostaguv(nu,nv,eou,eov)
! surface height at staggered coordinate
oeu(1:ifull)=0.5*(neta(1:ifull)+neta(ie))*eeu(1:ifull) ! height at staggered coordinate
oev(1:ifull)=0.5*(neta(1:ifull)+neta(in))*eev(1:ifull) ! height at staggered coordinate
eou(1:ifull,wlev+1)=oeu(1:ifull)
eov(1:ifull,wlev+1)=oev(1:ifull)
call boundsuv(eou,eov,stag=-9)
oeu(ifull+1:ifull+iextra)=eou(ifull+1:ifull+iextra,wlev+1)
oev(ifull+1:ifull+iextra)=eov(ifull+1:ifull+iextra,wlev+1)
ccu(:,1)=eou(:,1)*godsig(1)
ccv(:,1)=eov(:,1)*godsig(1)
do ii=2,wlev
  ccu(:,ii)=(ccu(:,ii-1)+eou(:,ii)*godsig(ii))
  ccv(:,ii)=(ccv(:,ii-1)+eov(:,ii)*godsig(ii))
end do
sdiv=(ccu(1:ifull,wlev)*max(ddu(1:ifull)+oeu(1:ifull),0.)/emu(1:ifull)-ccu(iwu,wlev)*max(ddu(iwu)+oeu(iwu),0.)/emu(iwu)  &
     +ccv(1:ifull,wlev)*max(ddv(1:ifull)+oev(1:ifull),0.)/emv(1:ifull)-ccv(isv,wlev)*max(ddv(isv)+oev(isv),0.)/emv(isv)) &
     *em(1:ifull)*em(1:ifull)/ds
do ii=1,wlev-1
  nw(:,ii)=ee(1:ifull)*(sdiv*gosigh(ii)-                                                                                 &
      (ccu(1:ifull,ii)*max(ddu(1:ifull)+oeu(1:ifull),0.)/emu(1:ifull)-ccu(iwu,ii)*max(ddu(iwu)+oeu(iwu),0.)/emu(iwu)     &
      +ccv(1:ifull,ii)*max(ddv(1:ifull)+oev(1:ifull),0.)/emv(1:ifull)-ccv(isv,ii)*max(ddv(isv)+oev(isv),0.)/emv(isv))    &
      *em(1:ifull)*em(1:ifull)/ds)
end do
! compute contunity equation horizontal transport terms
do ii=1,wlev
  mps(1:ifull,ii)=neta(1:ifull)-(1.-ocneps)*0.5*dt*ee(1:ifull)*                                                          &
     ((eou(1:ifull,ii)*(ddu(1:ifull)+neta(1:ifull))/emu(1:ifull)-eou(iwu,ii)*(ddu(iwu)+neta(1:ifull))/emu(iwu)           &
      +eov(1:ifull,ii)*(ddv(1:ifull)+neta(1:ifull))/emv(1:ifull)-eov(isv,ii)*(ddv(isv)+neta(1:ifull))/emv(isv))          &
      *em(1:ifull)*em(1:ifull)/ds+(nw(:,ii)-nw(:,ii-1))/godsig(ii)) ! nw is at half levels
end do

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: pressure gradient terms"
end if
#endif

! ocean
! Prepare pressure gradient terms at t=t and incorporate into velocity field
do ii=1,wlev
  rhou=0.5*(rho(1:ifull,ii)+rho(ie,ii))
  rhov=0.5*(rho(1:ifull,ii)+rho(in,ii))
  rhobaru=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
  rhobarv=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))
  tau(:,ii)=grav*(gosig(ii)*max(oeu(1:ifull)+ddu(1:ifull),0.)*drhobardxu(:,ii)                                           &
           +(rhobaru+(1.-gosig(ii))*rhou)*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds)/rhou                                 &
           +(dpsdxu+grav*rhou*dttdxu)/rhou ! staggered
  tav(:,ii)=grav*(gosig(ii)*max(oev(1:ifull)+ddv(1:ifull),0.)*drhobardyv(:,ii)                                           &
           +(rhobarv+(1.-gosig(ii))*rhov)*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds)/rhov                                 &
           +(dpsdyv+grav*rhov*dttdyv)/rhov
end do
! ice
!tau(:,wlev+1)=grav*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds ! staggered
!tav(:,wlev+1)=grav*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
!call mlounstaguv(tau,tav,ttau,ttav,toff=1)
call mlounstaguv(tau(:,1:wlev),tav(:,1:wlev),ttau(:,1:wlev),ttav(:,1:wlev),toff=1)
! ocean
do ii=1,wlev
  uau(:,ii)=nu(1:ifull,ii)+(1.-ocneps)*0.5*dt*( f(1:ifull)*nv(1:ifull,ii)-ttau(:,ii)) ! unstaggered
  uav(:,ii)=nv(1:ifull,ii)+(1.-ocneps)*0.5*dt*(-f(1:ifull)*nu(1:ifull,ii)-ttav(:,ii))
  uau(:,ii)=uau(:,ii)*ee(1:ifull)
  uav(:,ii)=uav(:,ii)*ee(1:ifull)
end do
! ice
snu(1:ifull)=i_u !-dt*ttau(:,wlev+1)
snv(1:ifull)=i_v !-dt*ttav(:,wlev+1)

call END_LOG(waterhadv_end)

call START_LOG(watervadv_begin)

#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "mlohadv: water vertical advection 1"
end if
#endif

! Vertical advection (first call for 0.5*dt)
call mlovadv(0.5*dt,nw,uau,uav,ns,nt,mps,depdum,dzdum,wtr(1:ifull),1)

call END_LOG(watervadv_end)

call START_LOG(waterhadv_begin)

#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "mlohadv: horizontal advection"
end if
#endif

! Convert (u,v) to cartesian coordinates (U,V,W)
do ii=1,wlev
  cou(1:ifull,ii,1)=ax(1:ifull)*uau(:,ii)+bx(1:ifull)*uav(:,ii)
  cou(1:ifull,ii,2)=ay(1:ifull)*uau(:,ii)+by(1:ifull)*uav(:,ii)
  cou(1:ifull,ii,3)=az(1:ifull)*uau(:,ii)+bz(1:ifull)*uav(:,ii)
end do

! Horizontal advection for U, V, W
call mlob2ints_uv(cou(:,:,1:3),nface,xg,yg,wtr)

! Rotate vector to arrival point
call mlorot(cou(:,:,1),cou(:,:,2),cou(:,:,3),x3d,y3d,z3d)

! Convert (U,V,W) back to conformal cubic coordinates
do ii=1,wlev
  uau(:,ii)=ax(1:ifull)*cou(1:ifull,ii,1)+ay(1:ifull)*cou(1:ifull,ii,2)+az(1:ifull)*cou(1:ifull,ii,3)
  uav(:,ii)=bx(1:ifull)*cou(1:ifull,ii,1)+by(1:ifull)*cou(1:ifull,ii,2)+bz(1:ifull)*cou(1:ifull,ii,3)
  uau(:,ii)=uau(:,ii)*ee(1:ifull)
  uav(:,ii)=uav(:,ii)*ee(1:ifull)
end do

! Horizontal advection for continuity
cou(1:ifull,1:wlev,1)=mps(1:ifull,1:wlev)
call mlob2ints(cou(:,:,1:1),nface,xg,yg,wtr)
mps(1:ifull,1:wlev)=cou(1:ifull,1:wlev,1)

do ii=1,wlev
  cou(1:ifull,ii,1)=nt(1:ifull,ii)
  cou(1:ifull,ii,2)=ns(1:ifull,ii)-34.72
end do

! Horizontal advection for T and S
call mlob2ints_bs(cou(:,:,1:2),nface,xg,yg,wtr)

! Convert (U,V,W) back to conformal cubic coordinates
do ii=1,wlev
  nt(1:ifull,ii) =cou(1:ifull,ii,1)
  ns(1:ifull,ii) =cou(1:ifull,ii,2)+34.72
end do

call END_LOG(waterhadv_end)

call START_LOG(watervadv_begin)

#ifdef debug
if ( myid==0 .and. nmaxpr==1 ) then
  write(6,*) "mlohadv: water vertical advection 2"
end if
#endif

! Vertical advection (second call for 0.5*dt)
! use explicit nw and depdum,dzdum from t=tau step (i.e., following JLM in CCAM atmospheric dynamics)
! Could use nuh and nvh to estimate nw at t+1/2, but also require an estimate of neta at t+1/2
call mlovadv(0.5*dt,nw,uau,uav,ns,nt,mps,depdum,dzdum,wtr(1:ifull),2)

! Integrate advected mps
xps(:)=mps(1:ifull,1)*godsig(1)
do ii=2,wlev
  xps(:)=xps(:)+mps(1:ifull,ii)*godsig(ii)
end do
xps(:)=xps(:)*ee(1:ifull)

call END_LOG(watervadv_end)

call START_LOG(watereos_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: density EOS 2"
end if
#endif

! Approximate normalised density rhobar at t+1 (unstaggered, using T and S at t+1)
if (nxtrrho==1) then
  cou(1:ifull,:,1)=nt(1:ifull,:)
  cou(1:ifull,:,2)=ns(1:ifull,:)
  call bounds(cou(:,:,1:2),corner=.true.)
  nt(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,1)
  ns(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,2)
  call mloexpdensity(rho,dalpha,dbeta,nt,ns,dzdum_rho,pice,0)
  rhobar(:,1)=(rho(:,1)-1030.)*godsig(1)
  do ii=2,wlev
    rhobar(:,ii)=rhobar(:,ii-1)+(rho(:,ii)-1030.)*godsig(ii)
  end do
  do ii=1,wlev
    rhobar(:,ii)=rhobar(:,ii)/gosigh(ii)+1030.
  end do

  ! update normalised density gradients
  ! method 2
  call tsjacobi(nt,ns,dalpha,dbeta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
  drhobardxu(:,1)=drhobardxu(:,1)*godsig(1)
  drhobardxv(:,1)=drhobardxv(:,1)*godsig(1)
  drhobardyu(:,1)=drhobardyu(:,1)*godsig(1)
  drhobardyv(:,1)=drhobardyv(:,1)*godsig(1)
  do ii=2,wlev
    drhobardxu(:,ii)=drhobardxu(:,ii-1)+drhobardxu(:,ii)*godsig(ii)
    drhobardxv(:,ii)=drhobardxv(:,ii-1)+drhobardxv(:,ii)*godsig(ii)
    drhobardyu(:,ii)=drhobardyu(:,ii-1)+drhobardyu(:,ii)*godsig(ii)
    drhobardyv(:,ii)=drhobardyv(:,ii-1)+drhobardyv(:,ii)*godsig(ii)
  end do
  do ii=1,wlev
    drhobardxu(:,ii)=drhobardxu(:,ii)/gosigh(ii)
    drhobardxv(:,ii)=drhobardxv(:,ii)/gosigh(ii)
    drhobardyu(:,ii)=drhobardyu(:,ii)/gosigh(ii)
    drhobardyv(:,ii)=drhobardyv(:,ii)/gosigh(ii)
  end do
end if

call END_LOG(watereos_end)

! FREE SURFACE CALCULATION ----------------------------------------

call START_LOG(waterhelm_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: free surface"
end if
#endif

! Prepare integral terms
sou=0.
spu=0.
squ=0.
ssu=0.
sov=0.
spv=0.
sqv=0.
ssv=0.

! ocean
! Precompute U,V current and integral terms at t+1
do ii=1,wlev
  au=uau(:,ii)
  av=uav(:,ii)
  tau(:,ii)=uau(:,ii)+(1.+ocneps)*0.5*dt*f(1:ifull)*av
  tav(:,ii)=uav(:,ii)-(1.+ocneps)*0.5*dt*f(1:ifull)*au
end do
! ice
tau(:,wlev+1)=snu(1:ifull)+dt*f(1:ifull)*snv(1:ifull) ! unstaggered
tav(:,wlev+1)=snv(1:ifull)-dt*f(1:ifull)*snu(1:ifull)
call mlostaguv(tau,tav,ttau,ttav)
! ocean
odum=1./(1.+(1.+ocneps)*(1.+ocneps)*0.25*dt*dt*fu(1:ifull)*fu(1:ifull))
odum=odum*eeu(1:ifull)
dumf=-(1.+ocneps)*0.5*dt*odum
do ii=1,wlev
  ccu(1:ifull,ii)=ttau(:,ii)*odum ! staggered
end do
odum=1./(1.+(1.+ocneps)*(1.+ocneps)*0.25*dt*dt*fv(1:ifull)*fv(1:ifull))
odum=odum*eev(1:ifull)
dumg=-(1.+ocneps)*0.5*dt*odum
do ii=1,wlev
  ccv(1:ifull,ii)=ttav(:,ii)*odum ! staggered
end do
! ice
! niu and niv hold the free drift solution (staggered).  Wind stress terms are updated in mlo.f90
niu(1:ifull)=ttau(:,wlev+1)/(1.+dt*dt*fu(1:ifull)*fu(1:ifull)) ! staggered
niv(1:ifull)=ttav(:,wlev+1)/(1.+dt*dt*fv(1:ifull)*fv(1:ifull))

do ii=1,wlev
  rhou   =0.5*(rho(1:ifull,ii)+rho(ie,ii))
  rhov   =0.5*(rho(1:ifull,ii)+rho(in,ii))
  rhobaru=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
  rhobarv=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))

  ! Create arrays to calcuate u and v at t+1, based on pressure gradient at t+1

  ! u^(t+1) = nu = au^(t*) + bu*dpdx^(t+1) + cu*dpdy^(t+1) (staggered)
  ! v^(t+1) = nv = av^(t*) + bv*dpdy^(t+1) + cv*dpdx^(t+1) (staggered)

  au=ccu(1:ifull,ii)
  bu=dumf/rhou
  cu= (1.+ocneps)*0.5*dt*bu ! fu now included in dpdy

  av=ccv(1:ifull,ii)
  bv=dumg/rhov
  cv=-(1.+ocneps)*0.5*dt*bv ! fv now included in dpdx

  ! Note pressure gradients are along constant z surfaces
  !dppdxu=dpsdxu+grav*sig*(etau+ddu)*drhobardxu|neta=0+grav*(rhobaru+(1-sig)*rhou)*detadxu
  !dppdyu=dpsdyu+grav*sig*(etau+ddu)*drhobardyu|neta=0+grav*(rhobaru+(1-sig)*rhou)*detadyu
  !dppdxv=dpsdxv+grav*sig*(etav+ddv)*drhobardxv|neta=0+grav*(rhobarv+(1-sig)*rhov)*detadxv
  !dppdyv=dpsdyv+grav*sig*(etav+ddv)*drhobardyv|neta=0+grav*(rhobarv+(1-sig)*rhov)*detadyv

  ! Create arrays for u and v at t+1 in terms of neta gradients
    
  !nu=kku+ppu+oou*etau+llu*(etau+ddu)+mmu*detadxu+nnu*detadyu (staggered)
  !nv=kkv+ppv+oov*etav+llv*(etav+ddv)+mmv*detadyv+nnv*detadxv (staggered)

  llu(:,ii)=grav*gosig(ii)*(bu*drhobardxu(:,ii)+cu*drhobardyu(:,ii)) ! drhobardyu contains dfdyu term
  mmu(:,ii)=bu*grav*(rhobaru+(1.-gosig(ii))*rhou)
  nnu(:,ii)=cu*grav*(rhobaru+(1.-gosig(ii))*rhou)
  oou(:,ii)=-nnu(:,ii)*dfdyu
  kku(:,ii)=au+bu*(dpsdxu+grav*rhou*dttdxu)-cu*dfdyu*(piceu+grav*rhou*tideu)
  ppu(:,ii)=cu*(dpsdyu+grav*rhou*dttdyu)

  llv(:,ii)=grav*gosig(ii)*(bv*drhobardyv(:,ii)+cv*drhobardxv(:,ii)) ! drhobardxv contains dfdxv term
  mmv(:,ii)=bv*grav*(rhobarv+(1.-gosig(ii))*rhov)
  nnv(:,ii)=cv*grav*(rhobarv+(1.-gosig(ii))*rhov)
  oov(:,ii)=-nnv(:,ii)*dfdxv
  kkv(:,ii)=av+bv*(dpsdyv+grav*rhov*dttdyv)-cv*dfdxv*(picev+grav*rhov*tidev)
  ppv(:,ii)=cv*(dpsdxv+grav*rhov*dttdxv)

  llu(:,ii)=llu(:,ii)*eeu(1:ifull)
  mmu(:,ii)=mmu(:,ii)*eeu(1:ifull)
  nnu(:,ii)=nnu(:,ii)*eeu(1:ifull)
  oou(:,ii)=oou(:,ii)*eeu(1:ifull)
  kku(:,ii)=kku(:,ii)*eeu(1:ifull)
  ppu(:,ii)=ppu(:,ii)*eeu(1:ifull)

  llv(:,ii)=llv(:,ii)*eev(1:ifull)
  mmv(:,ii)=mmv(:,ii)*eev(1:ifull)
  nnv(:,ii)=nnv(:,ii)*eev(1:ifull)
  oov(:,ii)=oov(:,ii)*eev(1:ifull)
  kkv(:,ii)=kkv(:,ii)*eev(1:ifull)
  ppv(:,ii)=ppv(:,ii)*eev(1:ifull)

  ! Pre-integrate arrays for u and v at t+1 (i.e., for calculating net divergence at t+1)

  !int nu dz = sou+ssu*etau+spu*(etau+ddu)+squ*detadxu
  !int nv dz = sov+ssv*etav+spv*(etav+ddv)+sqv*detadyv

  sou(1:ifull)=sou(1:ifull)+kku(:,ii)*godsig(ii)
  spu(1:ifull)=spu(1:ifull)+llu(:,ii)*godsig(ii)
  squ(1:ifull)=squ(1:ifull)+mmu(:,ii)*godsig(ii)
  ssu(1:ifull)=ssu(1:ifull)+oou(:,ii)*godsig(ii)

  sov(1:ifull)=sov(1:ifull)+kkv(:,ii)*godsig(ii)
  spv(1:ifull)=spv(1:ifull)+llv(:,ii)*godsig(ii)
  sqv(1:ifull)=sqv(1:ifull)+mmv(:,ii)*godsig(ii)
  ssv(1:ifull)=ssv(1:ifull)+oov(:,ii)*godsig(ii)

end do

! calculate terms for ice velocity at t+1

! (staggered)
! niu(t+1) = niu + idu*ipu + ibu*dipdx + icu*dipdy
! niv(t+1) = niv + idv*ipv + ibv*dipdy + icv*dipdx

imu=0.5*(imass(1:ifull)+imass(ie))
imv=0.5*(imass(1:ifull)+imass(in))
! note missing 0.5 as ip is the average of t and t+1
ibu(1:ifull)=-dt/(imu*(1.+dt*dt*fu(1:ifull)*fu(1:ifull)))
ibv(1:ifull)=-dt/(imv*(1.+dt*dt*fv(1:ifull)*fv(1:ifull)))
ibu(1:ifull)=ibu(1:ifull)*eeu(1:ifull)
ibv(1:ifull)=ibv(1:ifull)*eev(1:ifull)
icu(1:ifull)= dt*ibu(1:ifull) ! fu now included in dipdy
icv(1:ifull)=-dt*ibv(1:ifull) ! fv now includex in dipdx
idu(1:ifull)=-icu(1:ifull)*dfdyu
idv(1:ifull)=-icv(1:ifull)*dfdxv

dumc(1:ifull,1)=sou(1:ifull)
dumd(1:ifull,1)=sov(1:ifull)
dumc(1:ifull,2)=spu(1:ifull)
dumd(1:ifull,2)=spv(1:ifull)
dumc(1:ifull,3)=squ(1:ifull)
dumd(1:ifull,3)=sqv(1:ifull)
dumc(1:ifull,4)=ssu(1:ifull)
dumd(1:ifull,4)=ssv(1:ifull)
dumc(1:ifull,5)=ibu(1:ifull)
dumd(1:ifull,5)=ibv(1:ifull)
dumc(1:ifull,6)=idu(1:ifull)
dumd(1:ifull,6)=idv(1:ifull)
dumc(1:ifull,7)=niu(1:ifull)
dumd(1:ifull,7)=niv(1:ifull)
call boundsuv(dumc(:,1:7),dumd(:,1:7),stag=-9) ! stag=-9 updates iwu and isv
sou(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
sov(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,1)
spu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
spv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,2)
squ(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,3)
sqv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,3)
ssu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,4)
ssv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,4)
ibu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,5)
ibv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,5)
idu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,6)
idv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,6)
niu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,7)
niv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,7)

! prep ocean gradient terms
odiv=(sou(1:ifull)/emu(1:ifull)-sou(iwu)/emu(iwu)  &
     +sov(1:ifull)/emv(1:ifull)-sov(isv)/emv(isv)) &
     *em(1:ifull)*em(1:ifull)/ds
odivb=(sou(1:ifull)*ddu(1:ifull)/emu(1:ifull)-sou(iwu)*ddu(iwu)/emu(iwu)  &
      +sov(1:ifull)*ddv(1:ifull)/emv(1:ifull)-sov(isv)*ddv(isv)/emv(isv)) &
      *em(1:ifull)*em(1:ifull)/ds

pue=spu(1:ifull)*0.5/emu(1:ifull)
puw=spu(iwu)*0.5/emu(iwu)
pvn=spv(1:ifull)*0.5/emv(1:ifull)
pvs=spv(isv)*0.5/emv(isv)
pdiv=(pue-puw+pvn-pvs)*em(1:ifull)*em(1:ifull)/ds
pdivb=(pue*ddu(1:ifull)-puw*ddu(iwu)+pvn*ddv(1:ifull)-pvs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

sue=squ(1:ifull)/ds
suw=-squ(iwu)/ds
svn=sqv(1:ifull)/ds
svs=-sqv(isv)/ds
sdiv=(-sue+suw-svn+svs)*em(1:ifull)*em(1:ifull)/ds
sdivb=(-sue*ddu(1:ifull)+suw*ddu(iwu)-svn*ddv(1:ifull)+svs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

que=ssu(1:ifull)*0.5/emu(1:ifull)
quw=ssu(iwu)*0.5/emu(iwu)
qvn=ssv(1:ifull)*0.5/emv(1:ifull)
qvs=ssv(isv)*0.5/emv(isv)
qdiv=(que-quw+qvn-qvs)*em(1:ifull)*em(1:ifull)/ds
qdivb=(que*ddu(1:ifull)-quw*ddu(iwu)+qvn*ddv(1:ifull)-qvs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

! Iteratively solve for free surface height, eta
! Iterative loop to estimate ice 'pressure'
if (precon<-9999) then
  ! Multi-grid
  call mlomg(neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                               &
             pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                    &
             ipice,ibu,ibv,idu,idv,niu,niv,ipmax,totits,itc,maxglobseta,maxglobip)
else
  ! Usual SOR
  call mlosor(neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                              &
              pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                   &
              ipice,ibu,ibv,idu,idv,niu,niv,ipmax,totits,itc,maxglobseta,maxglobip)
end if

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: update currents"
end if
#endif

dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc(:,1:2),corner=.true.)
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

tnu=0.5*(neta(in)*f(in)+neta(ine)*f(ine))
tee=0.5*(neta(1:ifull)*f(1:ifull)+neta(ie)*f(ie))
tsu=0.5*(neta(is)*f(is)+neta(ise)*f(ise))
tev=0.5*(neta(ie)*f(ie)+neta(ien)*f(ien))
tnn=0.5*(neta(1:ifull)*f(1:ifull)+neta(in)*f(in))
twv=0.5*(neta(iw)*f(iw)+neta(iwn)*f(iwn))
detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
detadyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
detadxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
oeu(1:ifull)=0.5*(neta(1:ifull)+neta(ie))
oev(1:ifull)=0.5*(neta(1:ifull)+neta(in))

! Update currents once neta is calculated
do ii=1,wlev
  ! update currents (staggered)
  nu(1:ifull,ii)=kku(:,ii)+ppu(:,ii)+oou(:,ii)*oeu(1:ifull)+llu(:,ii)*max(oeu(1:ifull)+ddu(1:ifull),0.) &
                 +mmu(:,ii)*detadxu+nnu(:,ii)*detadyu
  nv(1:ifull,ii)=kkv(:,ii)+ppv(:,ii)+oov(:,ii)*oev(1:ifull)+llv(:,ii)*max(oev(1:ifull)+ddv(1:ifull),0.) &
                 +mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
end do

call END_LOG(waterhelm_end)

call START_LOG(wateriadv_begin)

! Update ice velocity with internal pressure terms
tnu=0.5*(ipice(in)*f(in)+ipice(ine)*f(ine))
tee=0.5*(ipice(1:ifull)*f(1:ifull)+ipice(ie)*f(ie))
tsu=0.5*(ipice(is)*f(is)+ipice(ise)*f(ise))
tev=0.5*(ipice(ie)*f(ie)+ipice(ien)*f(ien))
tnn=0.5*(ipice(1:ifull)*f(1:ifull)+ipice(in)*f(in))
twv=0.5*(ipice(iw)*f(iw)+ipice(iwn)*f(iwn))
dipdxu=(ipice(ie)-ipice(1:ifull))*emu(1:ifull)/ds
dipdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dipdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dipdyv=(ipice(in)-ipice(1:ifull))*emv(1:ifull)/ds
ipiceu=0.5*(ipice(1:ifull)+ipice(ie))
ipicev=0.5*(ipice(1:ifull)+ipice(in))
niu(1:ifull)=niu(1:ifull)+idu(1:ifull)*ipiceu+ibu(1:ifull)*dipdxu+icu(1:ifull)*dipdyu
niv(1:ifull)=niv(1:ifull)+idv(1:ifull)*ipicev+ibv(1:ifull)*dipdyv+icv(1:ifull)*dipdxv
niu(1:ifull)=niu(1:ifull)*eeu(1:ifull)
niv(1:ifull)=niv(1:ifull)*eev(1:ifull)
call boundsuv(niu,niv,stag=-9)

! Normalisation factor for conserving ice flow in and out of gridbox
spnet(1:ifull)=(-min(niu(iwu)*emu(iwu),0.)+max(niu(1:ifull)*emu(1:ifull),0.)     &
                -min(niv(isv)*emv(isv),0.)+max(niv(1:ifull)*emv(1:ifull),0.))/ds

! ADVECT ICE ------------------------------------------------------
! use simple upwind scheme

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Ice advection"
end if
#endif


! Horizontal advection for ice area
dumc(1:ifull,1)=fracice/(em(1:ifull)*em(1:ifull))             ! dumc(:,1) is an area
! Horizontal advection for ice volume
dumc(1:ifull,2)=sicedep*fracice/(em(1:ifull)*em(1:ifull))     ! dumc(:,2) is a volume
! Horizontal advection for snow volume
dumc(1:ifull,3)=snowd*0.001*fracice/(em(1:ifull)*em(1:ifull)) ! dumc(:,3) is a volume
! Horizontal advection for ice energy store
dumc(1:ifull,4)=i_sto*fracice/(em(1:ifull)*em(1:ifull))
! Horizontal advection for ice salinity
dumc(1:ifull,5)=i_sal*fracice*sicedep/(em(1:ifull)*em(1:ifull))
call mloexpgamm(gamm,sicedep,dumd(1:ifull,1),0)
! Horizontal advection for surface temperature
dumc(1:ifull,6)=i_it(1:ifull,1)*fracice*gamm(:,1)/(em(1:ifull)*em(1:ifull))
! Horizontal advection of snow temperature
dumc(1:ifull,7)=i_it(1:ifull,2)*fracice*gamm(:,2)/(em(1:ifull)*em(1:ifull))
! Horizontal advection of ice temperature
dumc(1:ifull,8)=i_it(1:ifull,3)*fracice*gamm(:,3)/(em(1:ifull)*em(1:ifull))
dumc(1:ifull,9)=i_it(1:ifull,4)*fracice*gamm(:,3)/(em(1:ifull)*em(1:ifull)) 
! Conservation
dumc(1:ifull,10)=spnet(1:ifull)
call bounds(dumc(:,1:10))
spnet(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,10)
do ii=1,9
  call upwind_iceadv(dumc(:,ii),niu,niv,spnet)
end do  
nfracice(1:ifull)=min(max(dumc(1:ifull,1)*em(1:ifull)*em(1:ifull),0.),maxicefrac)
ndic(1:ifull)=dumc(1:ifull,2)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
ndsn(1:ifull)=dumc(1:ifull,3)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
nsto(1:ifull)=dumc(1:ifull,4)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
nis(1:ifull) =dumc(1:ifull,5)*em(1:ifull)*em(1:ifull)/max(ndic(1:ifull)*nfracice(1:ifull),1.E-10)
call mloexpgamm(gamm,ndic,ndsn,0)
nit(1:ifull,1)=dumc(1:ifull,6)*em(1:ifull)*em(1:ifull)/max(gamm(:,1)*nfracice(1:ifull),1.E-10)
nit(1:ifull,2)=dumc(1:ifull,7)*em(1:ifull)*em(1:ifull)/max(gamm(:,2)*nfracice(1:ifull),1.E-10)
nit(1:ifull,3)=dumc(1:ifull,8)*em(1:ifull)*em(1:ifull)/max(gamm(:,3)*nfracice(1:ifull),1.E-10)
nit(1:ifull,4)=dumc(1:ifull,9)*em(1:ifull)*em(1:ifull)/max(gamm(:,3)*nfracice(1:ifull),1.E-10)

! populate grid points that have no sea ice
where ( nfracice(1:ifull)<1.E-4 .or. ndic(1:ifull)<1.E-4 )
  nfracice(1:ifull)=0.
  ndic(1:ifull)=0.
  ndsn(1:ifull)=0.
  nsto(1:ifull)=0.
  nis(1:ifull)=0.
  nit(1:ifull,1)=273.16
  nit(1:ifull,2)=273.16
  nit(1:ifull,3)=273.16
  nit(1:ifull,4)=273.16
elsewhere ( ndsn(1:ifull)<1.E-4 )
  ndsn(1:ifull)=0.
  nit(1:ifull,2)=273.16
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

call START_LOG(watermisc_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: conserve free surface and salinity"
end if
#endif

! volume conservation for water
if (nud_sfh==0) then
  neta(1:ifull)=max(min(neta(1:ifull),120.),-120.)
  odum=(neta(1:ifull)-w_e)*ee(1:ifull)
  call ccglobal_posneg(odum,delpos,delneg)
  alph_p = -delneg/max(delpos,1.E-20)
  alph_p = min(max(sqrt(alph_p),1.E-20),1.E20)
  neta(1:ifull)=w_e+max(0.,odum)*alph_p+min(0.,odum)/alph_p
end if

! temperature conservation (usually off when nudging SSTs)
if (nud_sst==0) then
  delpos=0.
  delneg=0.
  do ii=1,wlev
    where(wtr(1:ifull))
      dum(:,ii)=nt(1:ifull,ii)*max(dd(1:ifull)+neta(1:ifull),minwater)-w_t(:,ii)*max(dd(1:ifull)+w_e(1:ifull),minwater)
    elsewhere
      dum(:,ii)=0.
    end where
  end do
  call ccglobal_posneg(dum,delpos,delneg,dsigin=godsig)
  alph_p = -delneg/max(delpos,1.E-20)
  alph_p = min(sqrt(alph_p),alph_p)
  do ii=1,wlev
    where(wtr(1:ifull))
      nt(1:ifull,ii)=(w_t(:,ii)*max(dd(1:ifull)+w_e(1:ifull),minwater)           &
                     +max(0.,dum(:,ii))*alph_p+min(0.,dum(:,ii))/max(1.,alph_p)) &
                     /max(dd(1:ifull)+neta(1:ifull),minwater)
    end where
  end do
end if

! salinity conservation
if (nud_sss==0) then
  delpos=0.
  delneg=0.
  ndum=0.
  do ii=1,wlev
    ndum=ndum+w_s(1:ifull,ii)*godsig(ii)
  end do
  do ii=1,wlev
    where(wtr(1:ifull).and.ndum>0.)
      dum(:,ii)=ns(1:ifull,ii)*max(dd(1:ifull)+neta(1:ifull),minwater)-w_s(:,ii)*max(dd(1:ifull)+w_e(1:ifull),minwater)
    elsewhere
      dum(:,ii)=0.
    end where
  end do
  call ccglobal_posneg(dum,delpos,delneg,dsigin=godsig)
  alph_p = -delneg/max(delpos,1.E-20)
  alph_p = min(sqrt(alph_p),alph_p)
  do ii=1,wlev
    where(wtr(1:ifull).and.ndum>0.)
      ns(1:ifull,ii)=(w_s(:,ii)*max(dd(1:ifull)+w_e(1:ifull),minwater)           &
                     +max(0.,dum(:,ii))*alph_p+min(0.,dum(:,ii))/max(1.,alph_p)) &
                     /max(dd(1:ifull)+neta(1:ifull),minwater)
    end where
  end do
end if

if ( myid==0 .and. (ktau<=5.or.maxglobseta>tol.or.maxglobip>itol) ) then
  write(6,*) "MLODYNAMICS ",totits,itc,maxglobseta,maxglobip
end if

call END_LOG(watermisc_end)

! DIFFUSION -------------------------------------------------------------------
#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: diffusion"
end if
#endif

uau=av_vmod*nu(1:ifull,:)+(1.-av_vmod)*oldu1
uav=av_vmod*nv(1:ifull,:)+(1.-av_vmod)*oldv1
call mlodiffusion_main(uau,uav,nu,nv,neta,nt,ns)
call mloimport(4,neta(1:ifull),0,0) ! neta

! EXPORT ----------------------------------------------------------------------
! Water data is exported in mlodiffusion

call START_LOG(watermisc_begin)

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Export"
end if
#endif

do ii=1,4
  call mloimpice(nit(1:ifull,ii),ii,0)
end do
call mloimpice(nfracice(1:ifull),5,0)
call mloimpice(ndic(1:ifull),6,0)
call mloimpice(ndsn(1:ifull),7,0)
call mloimpice(nsto(1:ifull),8,0)
call mloimpice(niu(1:ifull),9,0)
call mloimpice(niv(1:ifull),10,0)
call mloimpice(nis(1:ifull),11,0)
where (wtr(1:ifull))
  fracice=nfracice(1:ifull)
  sicedep=ndic(1:ifull)
  snowd=ndsn(1:ifull)*1000.
end where

#ifdef debug
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Finish"
end if
#endif

call END_LOG(watermisc_end)

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(dt_in,ubar,vbar,nface,xg,yg,x3d,y3d,z3d,wtr)

use cc_mpi
use indices_m
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'parmhor.h'

integer iq,i,j,k,n,nn,idel,jdel,intsch,ncount,ii
integer, dimension(ifull,wlev), intent(out) :: nface
real, intent(in) :: dt_in
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real, dimension(ifull,wlev), intent(out) :: xg,yg
real(kind=8), dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc
real, dimension(ifull+iextra,wlev,3) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,3) :: sx
real, dimension(-1:2,-1:2,3) :: sc
real, dimension(4) :: cmul, emul, rmul
real, dimension(2:3) :: dmul
real xxg,yyg
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do k=1,wlev
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

s(1:ifull,:,1) = uc
s(1:ifull,:,2) = vc
s(1:ifull,:,3) = wc

intsch = mod(ktau,2)
sc = cxx-1.

do k=1,wlev
  where (.not.wtr(1:ifull))
    s(1:ifull,k,1) = cxx-1.
    s(1:ifull,k,2) = cxx-1.
    s(1:ifull,k,3) = cxx-1.
  end where
end do
s(ifull+1:,:,:) = cxx-1.
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:3) = reshape( s(1:ipan*jpan*npan,1:wlev,1:3), (/ ipan, jpan, npan, wlev, 3 /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do
    sx(-1,0,n,:,:)          = s(lwws(n),:,:)
    sx(0,0,n,:,:)           = s(iws(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lwss(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),:,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(1-ipan+n*ipan*jpan),  :,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(n*ipan*jpan),         :,:)
  end do               ! n loop
 
! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff
      
      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2)=(1.-xxg)
        dmul(3)=xxg
        emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn=1,3
          rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
          rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
          rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
          rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*3) = sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,3
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*3)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop
  
  call intssync_send(3)

  do k=1,wlev
    do iq=1,ifull
      idel = int(xg(iq,k))
      xxg  = xg(iq,k) - idel
      jdel = int(yg(iq,k))
      yyg  = yg(iq,k) - jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
           
        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2)=(1.-xxg)
          dmul(3)=xxg
          emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn = 1,3
            rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
            rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
            rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
            rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,3
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
    
      end if
    end do       ! iq loop
  end do         ! k loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:3) = reshape( s(1:ipan*jpan*npan,1:wlev,1:3), (/ ipan, jpan, npan, wlev, 3 /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(1-ipan+n*ipan*jpan),  :,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(n*ipan*jpan),         :,:)
 end do               ! n loop

! For other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn=1,3
          rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
          rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
          rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
          rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*3) = sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,3
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*3)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(3)

  do k=1,wlev
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn=1,3
            rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
            rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
            rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
            rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
            s(iq,k,nn) = sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,3
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k=1,wlev
  where (wtr(1:ifull))
    x3d(:,k) = x - 0.5*(uc(:,k)+s(1:ifull,k,1)) ! n+1 guess
    y3d(:,k) = y - 0.5*(vc(:,k)+s(1:ifull,k,2)) ! n+1 guess
    z3d(:,k) = z - 0.5*(wc(:,k)+s(1:ifull,k,3)) ! n+1 guess
  elsewhere
    x3d(:,k) = x
    y3d(:,k) = y
    z3d(:,k) = z
    s(1:ifull,k,1) = cxx-1.
    s(1:ifull,k,2) = cxx-1.
    s(1:ifull,k,3) = cxx-1.
  end where
end do
call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
!     Share off processor departure points.
call deptsync(nface,xg,yg)

!======================== start of intsch=1 section ====================
if(intsch==1)then
 
  ! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff
      
      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2)=(1.-xxg)
        dmul(3)=xxg
        emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn=1,3
          rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
          rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
          rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
          rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*3) = sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,3
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*3)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop
  
  call intssync_send(3)

  do k=1,wlev
    do iq=1,ifull
      idel = int(xg(iq,k))
      xxg  = xg(iq,k) - idel
      jdel = int(yg(iq,k))
      yyg  = yg(iq,k) - jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then
           
        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2)=(1.-xxg)
          dmul(3)=xxg
          emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn = 1,3
            rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
            rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
            rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
            rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,3
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do       ! iq loop
  end do         ! k loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================

  ! For other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn=1,3
          rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
          rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
          rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
          rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*3) = sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,3
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*3)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(3)

  do k=1,wlev
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn=1,3
            rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
            rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
            rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
            rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
            s(iq,k,nn) = sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,3
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do k=1,wlev
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

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices
! This code is from depts.f

subroutine mlotoij5(x3d,y3d,z3d,nface,xg,yg)

use bigxy4_m
use cc_mpi
use mlo
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmgeom.h'

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

alf=(1._8-schmidt*schmidt)/(1._8+schmidt*schmidt)
alfonsch=2._8*schmidt/(1._8+schmidt*schmidt)

do ii=1,wlev

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
      ri(iq)=real(real(i)+real(is)*((xg(iq,ii)-xx4(i,j))*dyy-(yg(iq,ii)-yy4(i,j))*dyx)/den(iq))
      rj(iq)=real(real(j)+real(js)*((yg(iq,ii)-yy4(i,j))*dxx-(xg(iq,ii)-xx4(i,j))*dxy)/den(iq))
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

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'

integer idel,iq,jdel
integer i,j,k,n,intsch,ncount
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(-1:2,-1:2,size(s,3)) :: sc
real, dimension(4) :: cmul, emul, rmul
real, dimension(2:3) :: dmul
real xxg,yyg
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

ntr=size(s,3)
intsch=mod(ktau,2)
sc=cxx-1.

do nn=1,ntr
  do k=1,wlev
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn)=cxx-1. ! missing value flag
    end where
  end do
end do
s(ifull+1:,:,:)=cxx-1.
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if (intsch==1) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lwws(n),                  :,:)
    sx(0,0,n,:,:)           = s(iws(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lwss(n),                  :,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),                  :,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),                  :,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),                  :,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),                  :,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),                  :,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),                  :,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(1-ipan+n*ipan*jpan),  :,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(n*ipan*jpan),         :,:)
  end do               ! n loop

! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2)=(1.-xxg)
        dmul(3)=xxg
        emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
          rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
          rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
          rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev      
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2)=(1.-xxg)
          dmul(3)=xxg
          emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
            rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
            rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
            rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do       ! iq loop
  end do         ! k loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(1-ipan+n*ipan*jpan),  :,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(n*ipan*jpan),         :,:)
  end do               ! n loop

! For other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
          rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
          rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
          rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        where (sc(0:1,0:1,:)<=cxx)
          sc(0:1,0:1,:)=0.
        end where
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
            rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
            rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
            rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          where (sc(0:1,0:1,:)<=cxx)
            sc(0:1,0:1,:)=0.
          end where
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn)=0.
    end where
  end do
end do

return
end subroutine mlob2ints_uv


subroutine mlob2ints(s,nface,xg,yg,wtr)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'

integer idel,iq,jdel
integer i,j,k,n,intsch,ncount
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(-1:2,-1:2,size(s,3)) :: sc
real, dimension(4) :: cmul, emul, rmul
real, dimension(2:3) :: dmul
real xxg,yyg
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

ntr=size(s,3)
intsch=mod(ktau,2)
sc=cxx-1.

do nn=1,ntr
  do k=1,wlev
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn)=cxx-1. ! missing value flag
    end where
  end do
end do
s(ifull+1:,:,:)=cxx-1.
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if (intsch==1) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lwws(n),                  :,:)
    sx(0,0,n,:,:)           = s(iws(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lwss(n),                  :,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),                  :,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),                  :,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),                  :,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),                  :,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),                  :,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),                  :,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(1-ipan+n*ipan*jpan),  :,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(n*ipan*jpan),         :,:)
  end do               ! n loop

! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2)=(1.-xxg)
        dmul(3)=xxg
        emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
          rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
          rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
          rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        call lfill(sc,cxx)
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev      
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2)=(1.-xxg)
          dmul(3)=xxg
          emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
            rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
            rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
            rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          call lfill(sc,cxx)
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do       ! iq loop
  end do         ! k loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(1-ipan+n*ipan*jpan),  :,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(n*ipan*jpan),         :,:)
  end do               ! n loop

! For other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
          rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
          rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
          rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        call lfill(sc,cxx)
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sextra(ii)%a(nn+(iq-1)*ntr)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
            rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
            rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
            rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
            s(iq,k,nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          call lfill(sc,cxx)
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            s(iq,k,nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        
      end if
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn)=0.
    end where
  end do
end do

return
end subroutine mlob2ints

subroutine mlob2ints_bs(s,nface,xg,yg,wtr)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'

integer idel,iq,jdel
integer i,j,k,n,intsch,ncount
integer ii,ntr,nn
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(ifull,wlev,size(s,3)) :: ssav
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev,size(s,3)) :: sx
real, dimension(-1:2,-1:2,size(s,3)) :: sc
real, dimension(4) :: cmul, emul, rmul
real, dimension(2:3) :: dmul
real, dimension(size(s,3)) :: sans, cmax, cmin
real xxg,yyg
real, parameter :: cxx = -9999. ! missing value flag
logical, dimension(ifull+iextra), intent(in) :: wtr

ntr=size(s,3)
intsch=mod(ktau,2)
sc=cxx-1.
ssav=s(1:ifull,:,:)

do nn=1,ntr
  do k=1,wlev
    where (.not.wtr(1:ifull))
      s(1:ifull,k,nn)=cxx-1. ! missing value flag
    end where
  end do
end do
s(ifull+1:,:,:)=cxx-1.
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if (intsch==1) then

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lwws(n),                  :,:)
    sx(0,0,n,:,:)           = s(iws(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lwss(n),                  :,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),                  :,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),                  :,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),                  :,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),                  :,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),                  :,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),                  :,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(1-ipan+n*ipan*jpan),  :,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(n*ipan*jpan),         :,:)
  end do               ! n loop

! Loop over points that need to be calculated for other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul(2)=(1.-xxg)
        dmul(3)=xxg
        emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
          rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
          rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
          rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
          sans(nn)=sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        call lfill(sc,cxx)
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sans(nn)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
      cmax(:)=max(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
      cmin(:)=min(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
      sextra(ii)%a(1+(iq-1)*ntr:iq*ntr) = min(max(sans(:),cmin(:)),cmax(:))
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev      
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)
     
        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          cmul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          cmul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          cmul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          dmul(2)=(1.-xxg)
          dmul(3)=xxg
          emul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          emul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          emul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          emul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(0:1,-1,nn)*dmul(:))
            rmul(2)=sum(sc(-1:2,0,nn)*cmul(:))
            rmul(3)=sum(sc(-1:2,1,nn)*cmul(:))
            rmul(4)=sum(sc(0:1,2,nn)*dmul(:))
            sans(nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          call lfill(sc,cxx)
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            sans(nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        cmax(:)=max(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
        cmin(:)=min(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
        s(iq,k,:) = min(max(sans(:),cmin(:)),cmax(:))
        
      end if    
    end do       ! iq loop
  end do         ! k loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  sx(1:ipan,1:jpan,1:npan,1:wlev,1:ntr) = reshape( s(1:ipan*jpan*npan,1:wlev,1:ntr), (/ ipan, jpan, npan, wlev, ntr /) )
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,:,:)      = s( iw(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(-1,j,n,:,:)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),:,:)
      sx(ipan+1,j,n,:,:) = s( ie(j*ipan+(n-1)*ipan*jpan),      :,:)
      sx(ipan+2,j,n,:,:) = s(iee(j*ipan+(n-1)*ipan*jpan),      :,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s( is(i+(n-1)*ipan*jpan), :,:)
      sx(i,-1,n,:,:)     = s(iss(i+(n-1)*ipan*jpan), :,:)
      sx(i,jpan+1,n,:,:) = s( in(i-ipan+n*ipan*jpan),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(i-ipan+n*ipan*jpan),:,:)
    end do            ! i loop
!   for ns interpolation, sometimes need (different from ew):
!        (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!      (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(1+(n-1)*ipan*jpan),   :,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(1-ipan+n*ipan*jpan),  :,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ipan+(n-1)*ipan*jpan),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(n*ipan*jpan),         :,:)
  end do               ! n loop

! For other processes
  do ii=neighnum,1,-1
    do iq=1,drlen(ii)
      n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
      !  Need global face index in fproc call
      idel = int(dpoints(ii)%a(2,iq))
      xxg = dpoints(ii)%a(2,iq) - idel
      jdel = int(dpoints(ii)%a(3,iq))
      yyg = dpoints(ii)%a(3,iq) - jdel
      k = nint(dpoints(ii)%a(4,iq))
      idel = idel - ioff
      jdel = jdel - joff

      sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
      sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
      sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12) then
        ! bi-cubic interpolation
        cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
        cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
        cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
        emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
        emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
        do nn=1,ntr
          rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
          rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
          rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
          rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
          sans(nn) = sum(rmul(:)*emul(:))
        end do
      else
        ! bi-linear interpolation
        call lfill(sc,cxx)
        cmul(2)=(1.-xxg)
        cmul(3)=xxg
        dmul(2)=(1.-yyg)
        dmul(3)=yyg
        do nn=1,ntr
          rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
          rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
          sans(nn)=sum(rmul(2:3)*dmul(2:3))
        end do
      end if
      cmax(:)=max(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
      cmin(:)=min(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
      sextra(ii)%a(1+(iq-1)*ntr:iq*ntr)=min(max(sans(:),cmin(:)),cmax(:))
    end do            ! iq loop
  end do              ! ii loop

  call intssync_send(ntr)

  do k=1,wlev
    do iq=1,ifull
      idel=int(xg(iq,k))
      xxg=xg(iq,k)-idel
      jdel=int(yg(iq,k))
      yyg=yg(iq,k)-jdel
      idel = idel - ioff
      jdel = jdel - joff
      n = nface(iq,k) + noff ! Make this a local index
      if ( idel>=0 .and. idel<=ipan .and. jdel>=0 .and. jdel<=jpan .and. n>=1 .and. n<=npan ) then

        sc(0:1,-1,:)   = sx(idel:idel+1,jdel-1,n,k,:)
        sc(-1:2,0:1,:) = sx(idel-1:idel+2,jdel:jdel+1,n,k,:)
        sc(0:1,2,:)    = sx(idel:idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12) then
          ! bi-cubic interpolation
          cmul(1)=(1.-yyg)*(2.-yyg)*(-yyg)/6.
          cmul(2)=(1.-yyg)*(2.-yyg)*(1.+yyg)/2.
          cmul(3)=yyg*(1.+yyg)*(2.-yyg)/2.
          cmul(4)=(1.-yyg)*(-yyg)*(1.+yyg)/6.
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          emul(1)=(1.-xxg)*(2.-xxg)*(-xxg)/6.
          emul(2)=(1.-xxg)*(2.-xxg)*(1.+xxg)/2.
          emul(3)=xxg*(1.+xxg)*(2.-xxg)/2.
          emul(4)=(1.-xxg)*(-xxg)*(1.+xxg)/6.
          do nn=1,ntr
            rmul(1)=sum(sc(-1,0:1,nn)*dmul(:))
            rmul(2)=sum(sc(0,-1:2,nn)*cmul(:))
            rmul(3)=sum(sc(1,-1:2,nn)*cmul(:))
            rmul(4)=sum(sc(2,0:1,nn)*dmul(:))
            sans(nn)=sum(rmul(:)*emul(:))
          end do
        else
          ! bi-linear interpolation
          call lfill(sc,cxx)
          cmul(2)=(1.-xxg)
          cmul(3)=xxg
          dmul(2)=(1.-yyg)
          dmul(3)=yyg
          do nn=1,ntr
            rmul(2)=sum(sc(0:1,0,nn)*cmul(2:3))
            rmul(3)=sum(sc(0:1,1,nn)*cmul(2:3))
            sans(nn)=sum(rmul(2:3)*dmul(2:3))
          end do
        end if
        cmax(:)=max(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
        cmin(:)=min(sc(0,0,:),sc(1,0,:),sc(0,1,:),sc(1,1,:))
        s(iq,k,:)=min(max(sans(:),cmin(:)),cmax(:))
        
      end if
    end do
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do nn=1,ntr
  do k=1,wlev
    where ( .not.wtr(1:ifull) .or. s(1:ifull,k,nn)<cxx+10. )
      s(1:ifull,k,nn)=ssav(1:ifull,k,nn)
    end where
  end do
end do

return
end subroutine mlob2ints_bs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use mlo
use xyzinfo_m

implicit none

include 'newmpar.h'

integer k
real, dimension(ifull+iextra,wlev), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real(kind=8), dimension(ifull,size(cou,2)), intent(in) :: x3d,y3d,z3d

do k=1,wlev
!         cross product n1xn2 into vec1
  vec1x = real(y3d(:,k)*z - y*z3d(:,k))
  vec1y = real(z3d(:,k)*x - z*x3d(:,k))
  vec1z = real(x3d(:,k)*y - x*y3d(:,k))
  denb = vec1x**2 + vec1y**2 + vec1z**2
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
  where (denb>1.e-4)
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

implicit none

include 'newmpar.h'
include 'parm.h'

integer k,itn,kx
real, dimension(:,:), intent(in) :: u, v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:), allocatable, save :: dtur,dtvr
real, dimension(ifull) :: nud,nvd
real, dimension(ifull,0:3) :: nwtu,nwtv
real, dimension(:), allocatable, save :: stul,stvl
real, dimension(:), allocatable, save :: stur,stvr
logical, dimension(ifull) :: euetest,euwtest,evntest,evstest
logical, dimension(ifull) :: euewtest,evnstest,eutest,evtest
logical ltest

call START_LOG(ocnstag_begin)

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

if (mstagf==0) then
  ltest=.true.
else if (mstagf<0) then
  ltest=.false.
else
  ! using ktau-1 ensures that the staggering is relative to the C grid
  ltest=mod(ktau-koff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then

  call boundsuv(uin,vin,stag=1)
  do k=1,kx
    ud(1:ifull,k)=dtul(:,1)*uin(ieeu,k)+dtul(:,2)*uin(ieu,k)+dtul(:,3)*uin(1:ifull,k)
    vd(1:ifull,k)=dtvl(:,1)*vin(innv,k)+dtvl(:,2)*vin(inv,k)+dtvl(:,3)*vin(1:ifull,k)
  end do
  
  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud=ud(1:ifull,k)-stul*ud(ieu,k)
    nvd=vd(1:ifull,k)-stvl*vd(inv,k)
    ud(1:ifull,k)=nud
    vd(1:ifull,k)=nvd

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(ieu,k))*0.5*eeu(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(inv,k))*0.5*eev(1:ifull)
    ua(1:ifull,k)=nud
    va(1:ifull,k)=nvd
  end do

  ! There are many ways to handle staggering near coastlines.
  ! This version supports the following properties:
  ! - Staggering is exactly reversible for the staggered (C grid) frame
  ! - Zero flux at land boundaries
  ! - The wave amplitude should be preserved

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtul(1:ifull,1)*ua(ieu,k)+wtul(1:ifull,2)*ua(iwu,k)+wtul(1:ifull,3)*ua(ieeu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvl(1:ifull,1)*va(inv,k)+wtvl(1:ifull,2)*va(isv,k)+wtvl(1:ifull,3)*va(innv,k)
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtul(1:ifull,1)*uin(ieu,k)+wtul(1:ifull,2)*uin(iwu,k)+wtul(1:ifull,3)*uin(ieeu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvl(1:ifull,1)*vin(inv,k)+wtvl(1:ifull,2)*vin(isv,k)+wtvl(1:ifull,3)*vin(innv,k)
    end do
  end do                 ! itn=1,itnmax

else

  call boundsuv(uin,vin)
  do k=1,kx
    ud(1:ifull,k)=dtur(1:ifull,1)*uin(iwu,k)+dtur(1:ifull,2)*uin(1:ifull,k)+dtur(1:ifull,3)*uin(ieu,k)
    vd(1:ifull,k)=dtvr(1:ifull,1)*vin(isv,k)+dtvr(1:ifull,2)*vin(1:ifull,k)+dtvr(1:ifull,3)*vin(inv,k)
  end do

  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud(1:ifull)=ud(1:ifull,k)-stur*ud(iwu,k)
    nvd(1:ifull)=vd(1:ifull,k)-stvr*vd(isv,k)
    ud(1:ifull,k)=nud(1:ifull)
    vd(1:ifull,k)=nvd(1:ifull)

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(ieu,k))*0.5*eeu(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(inv,k))*0.5*eev(1:ifull)
    ua(1:ifull,k)=nud(1:ifull)
    va(1:ifull,k)=nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtur(1:ifull,1)*ua(ieu,k)+wtur(1:ifull,2)*ua(iwu,k)+wtur(1:ifull,3)*ua(iwwu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvr(1:ifull,1)*va(inv,k)+wtvr(1:ifull,2)*va(isv,k)+wtvr(1:ifull,3)*va(issv,k)
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtur(1:ifull,1)*uin(ieu,k)+wtur(1:ifull,2)*uin(iwu,k)+wtur(1:ifull,3)*uin(iwwu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvr(1:ifull,1)*vin(inv,k)+wtvr(1:ifull,2)*vin(isv,k)+wtvr(1:ifull,3)*vin(issv,k)
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

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(in), optional :: toff
integer k,itn,kx,zoff
real, dimension(:,:), intent(in) :: u,v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(ifull) :: nud,nvd
real, dimension(ifull,0:3) :: nwtu,nwtv
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
  ltest=.true.
else if (mstagf<0) then
  ltest=.false.
else
  ltest=mod(ktau-zoff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then
  
  call boundsuv(uin,vin,stag=5)
  do k=1,kx
    ud(1:ifull,k)=dtul(1:ifull,1)*uin(iwwu,k)+dtul(1:ifull,2)*uin(iwu,k)+dtul(1:ifull,3)*uin(1:ifull,k)
    vd(1:ifull,k)=dtvl(1:ifull,1)*vin(issv,k)+dtvl(1:ifull,2)*vin(isv,k)+dtvl(1:ifull,3)*vin(1:ifull,k)
  end do
  
  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud(1:ifull)=ud(1:ifull,k)-stul*ud(iwu,k)
    nvd(1:ifull)=vd(1:ifull,k)-stvl*vd(isv,k)
    ud(1:ifull,k)=nud(1:ifull)
    vd(1:ifull,k)=nvd(1:ifull)

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(iwu,k))*0.5*ee(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(isv,k))*0.5*ee(1:ifull)
    ua(1:ifull,k)=nud(1:ifull)
    va(1:ifull,k)=nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtul(1:ifull,1)*ua(ieu,k)+wtul(1:ifull,2)*ua(iwu,k)+wtul(1:ifull,3)*ua(iwwu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvl(1:ifull,1)*va(inv,k)+wtvl(1:ifull,2)*va(isv,k)+wtvl(1:ifull,3)*va(issv,k)
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtul(1:ifull,1)*uin(ieu,k)+wtul(1:ifull,2)*uin(iwu,k)+wtul(1:ifull,3)*uin(iwwu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvl(1:ifull,1)*vin(inv,k)+wtvl(1:ifull,2)*vin(isv,k)+wtvl(1:ifull,3)*vin(issv,k)
    end do

  end do                  ! itn=1,itnmax
  
else

  call boundsuv(uin,vin)
  do k=1,kx
    ud(1:ifull,k)=dtur(1:ifull,1)*uin(ieu,k)+dtur(1:ifull,2)*uin(1:ifull,k)+dtur(1:ifull,3)*uin(iwu,k)
    vd(1:ifull,k)=dtvr(1:ifull,1)*vin(inv,k)+dtvr(1:ifull,2)*vin(1:ifull,k)+dtvr(1:ifull,3)*vin(isv,k)
  end do

  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud(1:ifull)=ud(1:ifull,k)-stur*ud(ieu,k)
    nvd(1:ifull)=vd(1:ifull,k)-stvr*vd(inv,k)
    ud(1:ifull,k)=nud(1:ifull)
    vd(1:ifull,k)=nvd(1:ifull)

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(iwu,k))*0.5*ee(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(isv,k))*0.5*ee(1:ifull)
    ua(1:ifull,k)=nud(1:ifull)
    va(1:ifull,k)=nvd(1:ifull)
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtur(1:ifull,1)*ua(ieu,k)+wtur(1:ifull,2)*ua(iwu,k)+wtur(1:ifull,3)*ua(ieeu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvr(1:ifull,1)*va(inv,k)+wtvr(1:ifull,2)*va(isv,k)+wtvr(1:ifull,3)*va(innv,k)
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtur(1:ifull,1)*uin(ieu,k)+wtur(1:ifull,2)*uin(iwu,k)+wtur(1:ifull,3)*uin(ieeu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvr(1:ifull,1)*vin(inv,k)+wtvr(1:ifull,2)*vin(isv,k)+wtvr(1:ifull,3)*vin(innv,k)
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
use mlo

implicit none

include 'newmpar.h'

integer, intent(in) :: cnum
integer ii,iq,its_g
integer, dimension(ifull) :: its
real, intent(in) :: dtin
real, dimension(ifull) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depdum,idzdum
real, dimension(:,:), intent(inout) :: uu,vv,ss,tt,mm
real, dimension(ifull,wlev) :: dzdum
logical, dimension(ifull), intent(in) :: wtr

dzdum=max(idzdum,1.E-10)

! reduce time step to ensure stability
dtnew(:)=dtin
do iq=1,ifull
  if (wtr(iq)) then
    do ii=1,wlev-1
      ! this trick works if dzdum(iq,ii)<dzdum(iq,ii+1)
      dtnew(iq)=min(dtnew(iq),0.3*dzdum(iq,ii)/max(abs(ww(iq,ii)),1.E-12))
    end do
  end if
end do
its(:)=int(dtin/(dtnew(:)+0.01))+1
its_g=maxval(its(:))
if (its_g>500) then
  write(6,*) "MLOVERT myid,cnum,its_g",myid,cnum,its_g
end if
dtnew(:)=dtin/real(its(:))

ss=ss-34.72
call mlotvd(its,dtnew,ww,uu,depdum,dzdum)
call mlotvd(its,dtnew,ww,vv,depdum,dzdum)
call mlotvd(its,dtnew,ww,ss,depdum,dzdum)
call mlotvd(its,dtnew,ww,tt,depdum,dzdum)
call mlotvd(its,dtnew,ww,mm,depdum,dzdum)
ss=ss+34.72
ss=max(ss,0.)

return
end subroutine mlovadv

subroutine mlotvd(its,dtnew,ww,uu,depadj,dzadj)

use mlo

implicit none

include 'newmpar.h'

integer ii,i,iq,kp,kx
integer, dimension(ifull), intent(in) :: its
real, dimension(ifull), intent(in) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depadj,dzadj
real, dimension(:,:), intent(inout) :: uu
real, dimension(ifull,wlev-1) :: ff
real, dimension(ifull,0:wlev) :: delu
real, dimension(ifull) :: fl,fh,cc,rr

! f=(w*u) at half levels
! du/dt = u*dw/dz-df/dz = -w*du/dz

delu(1:ifull,0)=0.
do ii=1,wlev-1
  delu(1:ifull,ii)=uu(1:ifull,ii+1)-uu(1:ifull,ii)
end do
delu(1:ifull,wlev)=0.

! TVD part
do ii=1,wlev-1
  ! +ve ww is downwards to the ocean floor
  where ( ww(1:ifull,ii)>0.)
    rr(1:ifull)=delu(1:ifull,ii-1)/(delu(1:ifull,ii)+sign(1.e-20,delu(1:ifull,ii)))
    fl(1:ifull)=ww(1:ifull,ii)*uu(1:ifull,ii)
  elsewhere
    rr(1:ifull)=delu(1:ifull,ii+1)/(delu(1:ifull,ii)+sign(1.e-20,delu(1:ifull,ii)))  
    fl(1:ifull)=ww(1:ifull,ii)*uu(1:ifull,ii+1)
  end where
  cc(1:ifull)=max(0.,min(1.,2.*rr(1:ifull)),min(2.,rr(1:ifull))) ! superbee
  fh(1:ifull)=ww(1:ifull,ii)*0.5*(uu(1:ifull,ii)+uu(1:ifull,ii+1))         &
   -0.5*(uu(1:ifull,ii+1)-uu(1:ifull,ii))*ww(1:ifull,ii)**2*dtnew(1:ifull) &
   /max(depadj(1:ifull,ii+1)-depadj(1:ifull,ii),1.E-10)
  ff(1:ifull,ii)=fl(1:ifull)+cc(1:ifull)*(fh(1:ifull)-fl(1:ifull))
  !ff(:,ii)=ww(:,ii)*0.5*(uu(1:ifull,ii)+uu(1:ifull,ii+1)) ! explicit
end do
uu(1:ifull,1)=uu(1:ifull,1)+dtnew*(uu(1:ifull,1)*ww(1:ifull,1)-ff(1:ifull,1))/dzadj(1:ifull,1)
do ii=2,wlev-1
  uu(1:ifull,ii)=uu(1:ifull,ii)+dtnew*(uu(1:ifull,ii)*(ww(1:ifull,ii)-ww(1:ifull,ii-1))   &
                                      -ff(1:ifull,ii)+ff(1:ifull,ii-1))/dzadj(1:ifull,ii)
end do
uu(1:ifull,wlev)=uu(1:ifull,wlev)+dtnew*(-uu(1:ifull,wlev)*ww(1:ifull,wlev-1)+ff(1:ifull,wlev-1))/dzadj(1:ifull,wlev)


do iq=1,ifull
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
use mlo, only : wlev,wrtemp

implicit none

include 'newmpar.h'

integer ii
real, dimension(ifull+iextra,wlev), intent(in) :: nti, nsi, alphabar, betabar
real, dimension(ifull,wlev), intent(out) :: drhobardxu, drhobardyu, drhobardxv, drhobardyv
real, dimension(ifull+iextra,wlev,2) :: na
real, dimension(ifull) :: absu, bbsu, absv, bbsv
real, dimension(ifull,wlev,2) :: dnadxu, dnadxv, dnadyu, dnadyv

na(:,:,1)=min(max(271.-wrtemp,nti),373.-wrtemp)
na(:,:,2)=min(max(0.,  nsi),50. )-34.72

call seekdelta(na,dnadxu,dnadyu,dnadxv,dnadyv)

do ii=1,wlev
  absu=0.5*(alphabar(1:ifull,ii)+alphabar(ie,ii))*eeu(1:ifull)
  bbsu=0.5*(betabar(1:ifull,ii) +betabar(ie,ii) )*eeu(1:ifull)
  absv=0.5*(alphabar(1:ifull,ii)+alphabar(in,ii))*eev(1:ifull)
  bbsv=0.5*(betabar(1:ifull,ii) +betabar(in,ii) )*eev(1:ifull)

  ! This relationship neglects compression effects due to neta from the EOS.
  drhobardxu(1:ifull,ii)=-absu*dnadxu(1:ifull,ii,1)+bbsu*dnadxu(1:ifull,ii,2)
  drhobardxv(1:ifull,ii)=-absv*dnadxv(1:ifull,ii,1)+bbsv*dnadxv(1:ifull,ii,2)
  drhobardyu(1:ifull,ii)=-absu*dnadyu(1:ifull,ii,1)+bbsu*dnadyu(1:ifull,ii,2)
  drhobardyv(1:ifull,ii)=-absv*dnadyv(1:ifull,ii,1)+bbsv*dnadyv(1:ifull,ii,2)
end do

return
end subroutine tsjacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using an interpolation method

subroutine seekdelta(rhobar,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use cc_mpi
use indices_m
use map_m
use mlo, only : wlev

implicit none

include 'newmpar.h'
include 'parm.h'

integer ii, jj
real, dimension(ifull,wlev) :: ddux,ddvy
real, dimension(ifull,wlev) :: ddi,dde,ddw,ddn,dds,dden,ddse,ddne,ddwn
real, dimension(ifull,wlev) :: ramp_a,ramp_c,ramp_e
real, dimension(ifull+iextra,wlev) :: dd_i
real, dimension(ifull,wlev,2) :: ri,re,rw,rn,rs,ren,rse,rne,rwn
real, dimension(ifull,wlev,2) :: ssi,sse,ssw,ssn,sss,ssen,ssse,ssne,sswn
real, dimension(ifull,wlev,2) :: y2i,y2e,y2w,y2n,y2s,y2en,y2se,y2ne,y2wn
real, dimension(ifull+iextra,wlev,2) :: ss_i, y2_i
real, dimension(ifull+iextra,wlev,2), intent (in) :: rhobar
real, dimension(ifull,wlev,2), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull) :: f_in,f_ine,f_ie,f_is,f_ise,f_ien,f_iw,f_iwn

! Here we calculate the slow contribution of the pressure gradient

! dP/dx = g rhobar dneta/dx + g sigma (D+neta) drhobar/dx
!                   (fast)                       (slow)

! rhobar = int_0^(z+neta) rho dz / (z+neta) = int_0^sigma rho dsigma / sigma

! MJT notes - this version fades out extrapolated gradients using ramp_a, etc.
!
! Idealy, we want to separate the neta contribution to drhobar/dx so that it
! can be included in the implicit solution to neta.

! drhobar/dx = drhobar/dx|neta=0 + (1-sigma) / (sigma*(D+neta)) rho dneta/dx
!
! dP/dx = g ( rhobar + (1-sigma) rho ) dneta/dx + g sigma (D+neta) drhobar/dx|neta=0

ss_i=rhobar(:,:,:)
do ii=1,wlev
  dd_i(:,ii)=gosig(ii)*dd(:)
  ddux(:,ii)=gosig(ii)*ddu(1:ifull)
  ddvy(:,ii)=gosig(ii)*ddv(1:ifull)
end do
call mlospline(dd_i,ss_i,y2_i) ! cubic spline

ssi(:,:,:)=rhobar(1:ifull,:,:)
sse(:,:,:)=rhobar(ie,:,:)
ssn(:,:,:)=rhobar(in,:,:)
ddi(:,:)  =dd_i(1:ifull,:)
dde(:,:)  =dd_i(ie,:)
ddn(:,:)  =dd_i(in,:)
y2i(:,:,:)=y2_i(1:ifull,:,:)
y2e(:,:,:)=y2_i(ie,:,:)
y2n(:,:,:)=y2_i(in,:,:)

sss(:,:,:) =rhobar(is,:,:)
ssne(:,:,:)=rhobar(ine,:,:)
ssse(:,:,:)=rhobar(ise,:,:)
dds(:,:)   =dd_i(is,:)
ddne(:,:)  =dd_i(ine,:)
ddse(:,:)  =dd_i(ise,:)
y2s(:,:,:) =y2_i(is,:,:)
y2ne(:,:,:)=y2_i(ine,:,:)
y2se(:,:,:)=y2_i(ise,:,:)

f_in(:) =f(in)
f_ine(:)=f(ine)
f_ie(:) =f(ie)
f_is(:) =f(is)
f_ise(:)=f(ise)
f_ien(:)=f(ien)
f_iw(:) =f(iw)
f_iwn(:)=f(iwn)
  
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
call seekval(rne,ssne,ddne,ddux,y2ne,ramp_c)
call seekval(rs, sss, dds, ddux,y2s, ramp_e)
call seekval(rse,ssse,ddse,ddux,y2se,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardyu(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,1)*(rn(:,ii,jj)*f_in+rne(:,ii,jj)*f_ine               &
                                                         -ri(:,ii,jj)*f(1:ifull)-re(:,ii,jj)*f_ie)*emu(1:ifull)/ds          &
                        -0.125*stwgt(1:ifull,1)*(ri(:,ii,jj)+re(:,ii,jj))*(f_in+f_ine-f(1:ifull)-f_ie)*emu(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,2)*(ri(:,ii,jj)*f(1:ifull)+re(:,ii,jj)*f_ie           &
                                                         -rs(:,ii,jj)*f_is-rse(:,ii,jj)*f_ise)*emu(1:ifull)/ds              &
                        -0.125*stwgt(1:ifull,2)*(ri(:,ii,jj)+re(:,ii,jj))*(f(1:ifull)+f_ie-f_is-f_ise)*emu(1:ifull)/ds)
  end do
end do

ssw(:,:,:) =rhobar(iw,:,:)
ssen(:,:,:)=rhobar(ien,:,:)
sswn(:,:,:)=rhobar(iwn,:,:)
ddw(:,:)   =dd_i(iw,:)
dden(:,:)  =dd_i(ien,:)
ddwn(:,:)  =dd_i(iwn,:)
y2w(:,:,:) =y2_i(iw,:,:)
y2en(:,:,:)=y2_i(ien,:,:)
y2wn(:,:,:)=y2_i(iwn,:,:)

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
call seekval(ren,ssen,dden,ddvy,y2en,ramp_c)
call seekval(rw, ssw, ddw, ddvy,y2W, ramp_e)
call seekval(rwn,sswn,ddwn,ddvy,y2wn,ramp_e)
do jj=1,2
  do ii=1,wlev
    drhobardxv(:,ii,jj)=ramp_a(:,ii)*ramp_c(:,ii)*(0.25*stwgt(1:ifull,3)*(re(:,ii,jj)*f_ie+ren(:,ii,jj)*f_ien               &
                                                         -ri(:,ii,jj)*f(1:ifull)-rn(:,ii,jj)*f_in)*emv(1:ifull)/ds          &
                        -0.125*stwgt(1:ifull,3)*(ri(:,ii,jj)+rn(:,ii,jj))*(f_ie+f_ien-f(1:ifull)-f_in)*emv(1:ifull)/ds)     &
                       +ramp_a(:,ii)*ramp_e(:,ii)*(0.25*stwgt(1:ifull,4)*(ri(:,ii,jj)*f(1:ifull)+rn(:,ii,jj)*f_in           &
                                                         -rw(:,ii,jj)*f_iw-rwn(:,ii,jj)*f_iwn)*emv(1:ifull)/ds              &
                        -0.125*stwgt(1:ifull,4)*(ri(:,ii,jj)+rn(:,ii,jj))*(f(1:ifull)+f_in-f_iw-f_iwn)*emv(1:ifull)/ds)
  end do
end do

return
end subroutine seekdelta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate to common depths

subroutine seekval(rout,ssin,ddin,ddseek,y2,ramp)

use cc_mpi
use mlo, only : wlev

implicit none

include 'newmpar.h'

integer iq, ii, jj, kk
integer, dimension(ifull,wlev) :: sindx
real, dimension(ifull,wlev), intent(in) :: ddseek
real, dimension(ifull,wlev), intent(in) :: ddin
real, dimension(ifull,wlev), intent(inout) :: ramp
real, dimension(ifull,wlev,2), intent(in) :: ssin, y2
real, dimension(ifull,wlev,2), intent(out) :: rout
real, dimension(ifull,2) :: ssunpack1, ssunpack0, y2unpack1, y2unpack0
real, dimension(ifull) :: h, a, b, tempa, tempb, temph
real, parameter :: dzramp = 0.1 ! extrapolation limit

sindx(:,:)=wlev
do iq = 1,ifull
  ii = 2
  do jj = 1,wlev
    do kk = ii,wlev-1
      if ( ddseek(iq,jj)<=ddin(iq,kk) ) then
        sindx(iq,jj) = kk
        exit
      end if
    end do
    ii = sindx(iq,jj)
  end do
end do
  
do  jj = 1,wlev
  ! MJT notes - This calculation is slow
  do iq = 1,ifull
    ii = sindx(iq,jj)
    h(iq) = max(ddin(iq,ii)-ddin(iq,ii-1),1.e-8)
    a(iq) = (ddin(iq,ii)-ddseek(iq,jj))/h(iq)
    !b(iq) = (ddseek(iq,jj)-ddin(iq,ii-1))/h(iq,jj)
    ssunpack1(iq,:) = ssin(iq,ii,:)
    ssunpack0(iq,:) = ssin(iq,ii-1,:)
    y2unpack1(iq,:) = y2(iq,ii,:)
    y2unpack0(iq,:) = y2(iq,ii-1,:)    
  end do

  b(:) = 1.-a(:)
  temph(:) = h(:)*h(:)/6.
  tempa(:) = (a(:)**3-a(:))*temph(:)
  tempb(:) = (b(:)**3-b(:))*temph(:)
  
  do ii = 1,2
    rout(:,jj,ii) = a(:)*ssunpack0(:,ii)+b(:)*ssunpack1(:,ii)           & ! linear interpolation
                   +tempa(:)*y2unpack0(:,ii)+tempb(:)*y2unpack1(:,ii)     ! cubic spline terms
  end do

  ! fade out extrapolation
  ramp(:,jj) = ramp(:,jj)*min(max((a(:)+dzramp)/dzramp,0.),1.)*min(max((b(:)+dzramp)/dzramp,0.),1.)
end do

return
end subroutine seekval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate coeff for cubic spline

subroutine mlospline(x,y,y2)

use mlo, only : wlev

implicit none

include 'newmpar.h'

integer ii, jj
real, dimension(ifull+iextra,wlev), intent(in) :: x
real, dimension(ifull+iextra,wlev,2), intent(in) :: y
real, dimension(ifull+iextra,wlev,2), intent(out) :: y2
real, dimension(ifull+iextra,wlev) :: u
real, dimension(ifull+iextra,wlev) :: sig
real, dimension(ifull+iextra) :: p

do ii=2,wlev-1
  sig(:,ii)=(x(:,ii)-x(:,ii-1))/max(x(:,ii+1)-x(:,ii-1),1.e-8)
end do

do jj=1,2
  y2(:,1,jj)=0.
  u(:,1)=0.
  do ii=2,wlev-1
    p(:)=sig(:,ii)*y2(:,ii-1,jj)+2.
    y2(:,ii,jj)=(sig(:,ii)-1.)/p(:)
    u(:,ii)=(6.*((y(:,ii+1,jj)-y(:,ii,jj))/max(x(:,ii+1)-x(:,ii),1.e-8)-(y(:,ii,jj)-y(:,ii-1,jj)) &
           /max(x(:,ii)-x(:,ii-1),1.e-8))/max(x(:,ii+1)-x(:,ii-1),1.e-8)-sig(:,ii)*u(:,ii-1))/p(:)
  end do
  y2(:,wlev,jj)=0.    
  do ii=wlev-1,1,-1
    y2(:,ii,jj)=y2(:,ii,jj)*y2(:,ii+1,jj)+u(:,ii)
  end do
end do

return
end subroutine mlospline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate tidal potential

subroutine mlotide(eta,slon,slat,mtimer,jstart)

implicit none

include 'newmpar.h'

integer, intent(in) :: mtimer,jstart
real, dimension(ifull), intent(out) :: eta
real, dimension(ifull), intent(in) :: slon,slat
real ltime,stime,ctime,mn,sn,pn
real, dimension(ifull) :: sinlat,coslat
real, dimension(4) :: aa,ab,ba,bb
real, parameter :: pi = 3.1415927

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
! this fill version is for a simply connected 2x2 grid
subroutine lfill(sc,sx)

implicit none

integer ncount
real, intent(in) :: sx
real, dimension(:,:,:), intent(inout) :: sc
real, dimension(1:4,1:4,size(sc,3)) :: nc
logical globc

ncount=count(sc(2:3,2:3,1)>sx)

! no missing values or all missing values
if (ncount>=4.or.ncount<=0) return

nc=sc
globc=.false.
call trial_sw(nc,sc,globc)
call trial_se(nc,sc,globc)
call trial_nw(nc,sc,globc)
call trial_ne(nc,sc,globc)
sc=nc
! if some points failed to fill, try again after first pass
if (globc) then
  call trial_sw(nc,sc,globc)
  call trial_se(nc,sc,globc)
  call trial_nw(nc,sc,globc)
  call trial_ne(nc,sc,globc)
  sc=nc
end if

return
end subroutine lfill

! This routine averages adjacent points for an average fill
subroutine trial_sw(nc,sc,globc)

implicit none

integer nec
real, dimension(:,:,:), intent(inout) :: nc
real, dimension(:,:,:), intent(in) :: sc
real, dimension(size(sc,3)) :: new
logical, intent(inout) :: globc

! this point has valid data
if (nc(2,2,1)>-990.) return

! attempt to fill data from adjacent point
new(:)=0.
nec=0
if (sc(3,2,1)>-990.) then
  new(:)=new(:)+sc(3,2,:)
  nec=nec+1
end if
if (sc(2,3,1)>-990.) then
  new(:)=new(:)+sc(2,3,:)
  nec=nec+1
end if
if (nec>0) then
  nc(2,2,:)=new(:)/real(nec)
else
  globc=.true. ! cannot find data flag
end if
  
return
end subroutine trial_sw

subroutine trial_se(nc,sc,globc)

implicit none

integer nec
real, dimension(:,:,:), intent(inout) :: nc
real, dimension(:,:,:), intent(in) :: sc
real, dimension(size(sc,3)) :: new
logical, intent(inout) :: globc

! this point has valid data
if (nc(3,2,1)>-990.) return

! attempt to fill data from adjacent point
new(:)=0.
nec=0
if (sc(3,3,1)>-990.) then
  new(:)=new(:)+sc(3,3,:)
  nec=nec+1
end if
if (sc(2,2,1)>-990.) then
  new(:)=new(:)+sc(2,2,:)
  nec=nec+1
end if
if (nec>0) then
  nc(3,2,:)=new(:)/real(nec)
else
  globc=.true. ! cannot find data flag
end if
  
return
end subroutine trial_se

subroutine trial_nw(nc,sc,globc)

implicit none

integer nec
real, dimension(:,:,:), intent(inout) :: nc
real, dimension(:,:,:), intent(in) :: sc
real, dimension(size(sc,3)) :: new
logical, intent(inout) :: globc

! this point has valid data
if (nc(2,3,1)>-990.) return

! attempt to fill data from adjacent point
new(:)=0.
nec=0
if (sc(2,2,1)>-990.) then
  new(:)=new(:)+sc(2,2,:)
  nec=nec+1
end if
if (sc(3,3,1)>-990.) then
  new(:)=new(:)+sc(3,3,:)
  nec=nec+1
end if
if (nec>0) then
  nc(2,3,:)=new(:)/real(nec)
else
  globc=.true. ! cannot find data flag
end if
  
return
end subroutine trial_nw

subroutine trial_ne(nc,sc,globc)

implicit none

integer nec
real, dimension(:,:,:), intent(inout) :: nc
real, dimension(:,:,:), intent(in) :: sc
real, dimension(size(sc,3)) :: new
logical, intent(inout) :: globc

! this point has valid data
if (nc(3,3,1)>-990.) return

! attempt to fill data from adjacent point
new(:)=0.
nec=0
if (sc(3,2,1)>-990.) then
  new(:)=new(:)+sc(3,2,:)
  nec=nec+1
end if
if (sc(2,3,1)>-990.) then
  new(:)=new(:)+sc(2,3,:)
  nec=nec+1
end if
if (nec>0) then
  nc(3,3,:)=new(:)/real(nec)
else
  globc=.true. ! cannot find data flag
end if
  
return
end subroutine trial_ne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test for leap year

subroutine mloleap(tyear,ttest)

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
! Advect ice using upwind scheme

subroutine upwind_iceadv(dumc,niu,niv,spnet)

use indices_m
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

real, dimension(ifull+iextra), intent(inout) :: dumc
real, dimension(ifull+iextra), intent(in) :: niu,niv,spnet
real(kind=8), dimension(ifull) :: odum,dumd,tdum

dumd=dumc(1:ifull)

odum=0.5*dt*(niu(iwu)*(dumc(1:ifull)+dumc(iw))    -abs(niu(iwu))*(dumc(1:ifull)-dumc(iw))    )*emu(iwu)/ds
where ( spnet(iw)>1.e-10 )
  tdum=dumc(iw)*max(niu(iwu)*emu(iwu),0.)/(ds*spnet(iw))    
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=dumc(1:ifull)*min(niu(iwu)*emu(iwu),0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum

odum=-0.5*dt*(niu(1:ifull)*(dumc(1:ifull)+dumc(ie))+abs(niu(1:ifull))*(dumc(1:ifull)-dumc(ie)))*emu(1:ifull)/ds
where ( spnet(ie)>1.e-10 )
  tdum=-dumc(ie)*min(niu(1:ifull)*emu(1:ifull),0.)/(ds*spnet(ie))
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

odum=0.5*dt*(niv(isv)*(dumc(1:ifull)+dumc(is))    -abs(niv(isv))*(dumc(1:ifull)-dumc(is))    )*emv(isv)/ds
where ( spnet(is)>1.e-10 )
  tdum=dumc(is)*max(niv(isv)*emv(isv),0.)/(ds*spnet(is))
  odum=min(odum,tdum)
elsewhere
  odum=min(odum,0._8)
end where
where ( spnet(1:ifull)>1.e-10 )
  tdum=dumc(1:ifull)*min(niv(isv)*emv(isv),0.)/(ds*spnet(1:ifull))
  odum=max(odum,tdum)
elsewhere
  odum=max(odum,0._8)
end where
dumd=dumd+odum

odum=-0.5*dt*(niv(1:ifull)*(dumc(1:ifull)+dumc(in))+abs(niv(1:ifull))*(dumc(1:ifull)-dumc(in)))*emv(1:ifull)/ds
where ( spnet(in)>1.e-10 )
  tdum=-dumc(in)*min(niv(1:ifull)*emv(1:ifull),0.)/(ds*spnet(in))    
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
! Use SOR to solve for free surface and ice pressure

subroutine mlosor(neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                             &
                  pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                  &
                  ipice,ibu,ibv,idu,idv,niu,niv,ipmax,totits,itc,maxglobseta,maxglobip)

use cc_mpi
use indices_m
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(out) :: totits,itc
!integer itstest,itsave1,itsave2
integer nx,ll
integer isc,iec
real, intent(out) :: maxglobseta,maxglobip
real maxloclseta,maxloclip
!real gd,ci,itserr1,itserr2
real, dimension(ifull+iextra), intent(inout) :: neta
real, dimension(ifull), intent(in) :: sue,svn,suw,svs
real, dimension(ifull), intent(in) :: pue,pvn,puw,pvs
real, dimension(ifull), intent(in) :: que,qvn,quw,qvs
real, dimension(ifull), intent(in) :: pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps
real, dimension(ifull+iextra), intent(inout) :: ipice
real, dimension(ifull+iextra), intent(in) :: ibu,ibv,idu,idv
real, dimension(ifull+iextra), intent(in) :: niu,niv
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull,2) :: zz,zzn,zzs,zze,zzw,rhs
real, dimension(ifull) :: yy,yyn,yys,yye,yyw
real, dimension(ifull) :: hh
real, dimension(ifull) :: seta,setab
real, dimension(ifullmaxcol) :: au,bu,cu,setac,nip
real, dimension(ifullmaxcol,maxcolour,2) :: zzc,zznc,zzsc,zzec,zzwc,rhsc
real, dimension(ifullmaxcol,maxcolour) :: yyc,yync,yysc,yyec,yywc
real, dimension(ifullmaxcol,maxcolour) :: hhc
real, dimension(ifull+iextra,2) :: dumc
real, dimension(2) :: dume,dumf

integer, parameter :: llmax      = 400 ! Iterations for calculating surface height

! Use SOR as it converge faster that Conjugate Gradient (problems with coastlines?)
! and it is easy to include the non-linear terms and constraints.  Possibly combine
! with a multi-grid method to improve convergence at large scales, although it is
! more likely that the dynamical core will be replaced with JLM's mass flux scheme
!itsave2=0
!itserr2=9.E9
!itstest=1
!itc=0

dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc(:,1:2))
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

! zz*(DIV^2 neta) + yy*neta*(DIV^2 neta) + hh*neta = rhs

! ocean
zzn(:,1)= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*ddv(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzs(:,1)=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*ddv(isv)    *em(1:ifull)*em(1:ifull)/ds
zze(:,1)= (1.+ocneps)*0.5*dt*(que+sue+pue)*ddu(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzw(:,1)=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*ddu(iwu)    *em(1:ifull)*em(1:ifull)/ds
zz(:,1) = (1.+ocneps)*0.5*dt*(qdivb+sdivb+pdivb)

yyn= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*em(1:ifull)*em(1:ifull)/ds
yys=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*em(1:ifull)*em(1:ifull)/ds
yye= (1.+ocneps)*0.5*dt*(que+sue+pue)*em(1:ifull)*em(1:ifull)/ds
yyw=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*em(1:ifull)*em(1:ifull)/ds
yy = (1.+ocneps)*0.5*dt*(qdiv+sdiv+pdiv)

hh     =1.+(1.+ocneps)*0.5*dt*(odiv+(pvn*dd(in)-pvs*dd(is)+pue*dd(ie)-puw*dd(iw))*em(1:ifull)**2/ds    &
                              +pdiv*dd(1:ifull))
rhs(:,1)=xps-(1.+ocneps)*0.5*dt*(odivb+(pvn*ddv(1:ifull)*dd(in)-pvs*ddv(isv)*dd(is)                    &
                                       +pue*ddu(1:ifull)*dd(ie)-puw*ddu(iwu)*dd(iw))*em(1:ifull)**2/ds &
                                 +pdivb*dd(1:ifull))

! ice
zzn(:,2)=(-idv(1:ifull)*0.5-ibv(1:ifull))/ds
zzs(:,2)=( idv(isv)*0.5    -ibv(isv)    )/ds
zze(:,2)=(-idu(1:ifull)*0.5-ibu(1:ifull))/ds
zzw(:,2)=( idu(iwu)*0.5    -ibu(iwu)    )/ds
zz(:,2) =(ibu(1:ifull)+ibu(iwu)+ibv(1:ifull)+ibv(isv)-0.5*(idu(1:ifull)-idu(iwu)+idv(1:ifull)-idv(isv)))/ds

rhs(:,2)=min(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu)+niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv),0.)

do nx=1,maxcolour
  zznc(1:ifullcol(nx),nx,:)=zzn(iqx(1:ifullcol(nx),nx),:)
  zzsc(1:ifullcol(nx),nx,:)=zzs(iqx(1:ifullcol(nx),nx),:)
  zzec(1:ifullcol(nx),nx,:)=zze(iqx(1:ifullcol(nx),nx),:)
  zzwc(1:ifullcol(nx),nx,:)=zzw(iqx(1:ifullcol(nx),nx),:)
  zzc(1:ifullcol(nx),nx,:) = zz(iqx(1:ifullcol(nx),nx),:)
  yync(1:ifullcol(nx),nx)  =yyn(iqx(1:ifullcol(nx),nx))
  yysc(1:ifullcol(nx),nx)  =yys(iqx(1:ifullcol(nx),nx))
  yyec(1:ifullcol(nx),nx)  =yye(iqx(1:ifullcol(nx),nx))
  yywc(1:ifullcol(nx),nx)  =yyw(iqx(1:ifullcol(nx),nx))
  yyc(1:ifullcol(nx),nx)   = yy(iqx(1:ifullcol(nx),nx))
  hhc(1:ifullcol(nx),nx)   = hh(iqx(1:ifullcol(nx),nx))
  rhsc(1:ifullcol(nx),nx,:)=rhs(iqx(1:ifullcol(nx),nx),:)
end do

do ll=1,llmax

  do nx=1,maxcolour

    isc = 1
    iec = ifullcol_border(nx)
      
    ! ocean
    ! 5-point version -----------------------------------------------

    ! For now, assume Boussinesq fluid and treat density in continuity equation as constant
    au(isc:iec)=yyc(isc:iec,nx)
    bu(isc:iec)=zzc(isc:iec,nx,1)+hhc(isc:iec,nx)+yync(isc:iec,nx)*neta(iqn(isc:iec,nx)) &
                                                 +yysc(isc:iec,nx)*neta(iqs(isc:iec,nx)) &
                                                 +yyec(isc:iec,nx)*neta(iqe(isc:iec,nx)) &
                                                 +yywc(isc:iec,nx)*neta(iqw(isc:iec,nx))
    cu(isc:iec)=-rhsc(isc:iec,nx,1)+zznc(isc:iec,nx,1)*neta(iqn(isc:iec,nx)) &
                                   +zzsc(isc:iec,nx,1)*neta(iqs(isc:iec,nx)) &
                                   +zzec(isc:iec,nx,1)*neta(iqe(isc:iec,nx)) &
                                   +zzwc(isc:iec,nx,1)*neta(iqw(isc:iec,nx))
    
    ! alternative form
    setac(isc:iec)=-neta(iqx(isc:iec,nx))-2.*cu(isc:iec)/(bu(isc:iec) &
           +sqrt(bu(isc:iec)*bu(isc:iec)-4.*au(isc:iec)*cu(isc:iec)))
    
    ! The following expression limits the minimum depth
    ! (should not occur for typical eta values)
    seta(iqx(isc:iec,nx))=max(setac(isc:iec),-(neta(iqx(isc:iec,nx))+dd(iqx(isc:iec,nx)))) 
    seta(iqx(isc:iec,nx))=seta(iqx(isc:iec,nx))*ee(iqx(isc:iec,nx))
    dumc(iqx(isc:iec,nx),1)=neta(iqx(isc:iec,nx))+seta(iqx(isc:iec,nx))

    ! ice
    ! 5-point version -------------------------------------------------
 
    bu(isc:iec)=zzc(isc:iec,nx,2)
    cu(isc:iec)=-rhsc(isc:iec,nx,2)+zznc(isc:iec,nx,2)*ipice(iqn(isc:iec,nx)) &
                                   +zzsc(isc:iec,nx,2)*ipice(iqs(isc:iec,nx)) &
                                   +zzec(isc:iec,nx,2)*ipice(iqe(isc:iec,nx)) &
                                   +zzwc(isc:iec,nx,2)*ipice(iqw(isc:iec,nx))
    nip(isc:iec)=0.
    where (bu(isc:iec)/=0.)
      nip(isc:iec)=-cu(isc:iec)/bu(isc:iec)
    end where
 
    ! cavitating fluid
    nip(isc:iec)=max(min(nip(isc:iec),ipmax(iqx(isc:iec,nx))),0.)

    setab(iqx(isc:iec,nx))=nip(isc:iec)-ipice(iqx(isc:iec,nx))
    dumc(iqx(isc:iec,nx),2)=ipice(iqx(isc:iec,nx))+setab(iqx(isc:iec,nx))

    call bounds_colour_send(dumc(:,1:2),nx)
    
    isc = ifullcol_border(nx) + 1
    iec = ifullcol(nx)
      
    ! ocean
    ! 5-point version -----------------------------------------------

    ! For now, assume Boussinesq fluid and treat density in continuity equation as constant
    au(isc:iec)=yyc(isc:iec,nx)
    bu(isc:iec)=zzc(isc:iec,nx,1)+hhc(isc:iec,nx)+yync(isc:iec,nx)*neta(iqn(isc:iec,nx)) &
                                                 +yysc(isc:iec,nx)*neta(iqs(isc:iec,nx)) &
                                                 +yyec(isc:iec,nx)*neta(iqe(isc:iec,nx)) &
                                                 +yywc(isc:iec,nx)*neta(iqw(isc:iec,nx))
    cu(isc:iec)=-rhsc(isc:iec,nx,1)+zznc(isc:iec,nx,1)*neta(iqn(isc:iec,nx)) &
                                   +zzsc(isc:iec,nx,1)*neta(iqs(isc:iec,nx)) &
                                   +zzec(isc:iec,nx,1)*neta(iqe(isc:iec,nx)) &
                                   +zzwc(isc:iec,nx,1)*neta(iqw(isc:iec,nx))
    
    ! alternative form
    setac(isc:iec)=-neta(iqx(isc:iec,nx))-2.*cu(isc:iec)/(bu(isc:iec) &
           +sqrt(bu(isc:iec)*bu(isc:iec)-4.*au(isc:iec)*cu(isc:iec)))
    
    ! The following expression limits the minimum depth
    ! (should not occur for typical eta values)
    seta(iqx(isc:iec,nx))=max(setac(isc:iec),-(neta(iqx(isc:iec,nx))+dd(iqx(isc:iec,nx)))) 
    seta(iqx(isc:iec,nx))=seta(iqx(isc:iec,nx))*ee(iqx(isc:iec,nx))
    neta(iqx(isc:iec,nx))=neta(iqx(isc:iec,nx))+seta(iqx(isc:iec,nx))

    ! ice
    ! 5-point version -------------------------------------------------
 
    bu(isc:iec)=zzc(isc:iec,nx,2)
    cu(isc:iec)=-rhsc(isc:iec,nx,2)+zznc(isc:iec,nx,2)*ipice(iqn(isc:iec,nx)) &
                                   +zzsc(isc:iec,nx,2)*ipice(iqs(isc:iec,nx)) &
                                   +zzec(isc:iec,nx,2)*ipice(iqe(isc:iec,nx)) &
                                   +zzwc(isc:iec,nx,2)*ipice(iqw(isc:iec,nx))
    nip(isc:iec)=0.
    where (bu(isc:iec)/=0.)
      nip(isc:iec)=-cu(isc:iec)/bu(isc:iec)
    end where
 
    ! cavitating fluid
    nip(isc:iec)=max(min(nip(isc:iec),ipmax(iqx(isc:iec,nx))),0.)

    setab(iqx(isc:iec,nx))=nip(isc:iec)-ipice(iqx(isc:iec,nx))
    ipice(iqx(isc:iec,nx))=ipice(iqx(isc:iec,nx))+setab(iqx(isc:iec,nx))

    dumc(iqx(isc:iec,nx),1)=neta(iqx(isc:iec,nx))
    dumc(iqx(isc:iec,nx),2)=ipice(iqx(isc:iec,nx))
    iec = ifullcol_border(nx)
    neta(iqx(1:iec,nx)) =dumc(iqx(1:iec,nx),1)
    ipice(iqx(1:iec,nx))=dumc(iqx(1:iec,nx),2)
    call bounds_colour_recv(dumc(:,1:2),nx)    
    neta(ifull+1:ifull+iextra) =dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  
  end do
  
  maxloclseta=maxval(abs(seta))
  maxloclip  =maxval(abs(setab))

  !if (ll>=itstest) then
    dume(1)=maxloclseta
    dume(2)=maxloclip
    call ccmpi_allreduce(dume(1:2),dumf(1:2),"max",comm_world)
    maxglobseta=dumf(1)
    maxglobip  =dumf(2)

    if (maxglobseta<tol.and.maxglobip<itol) exit
    
  !  itsave1=itsave2
  !  itsave2=ll
  !  itserr1=itserr2
  !  itserr2=log10(maxglobseta)
  !  
  !  gd=(itserr2-itserr1)/real(itsave2-itsave1)
  !  ci=itserr2-gd*real(itsave2)
  !  if (gd/=0.) then
  !    itstest=nint((log10(tol)-ci)/gd)
  !    itstest=max(itstest,ll+1)
  !  else
  !    itstest=ll+1
  !  end if
  !  itc=itc+1    
  !end if

end do
totits=ll
itc=totits

return
end subroutine mlosor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use multi-grid to solve for free surface and ice pressure

subroutine mlomg(neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                              &
                  pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                  &
                  ipice,ibu,ibv,idu,idv,niu,niv,ipmax,totits,itc,maxglobseta,maxglobip)

use helmsolve
use indices_m
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(out) :: totits,itc
real, intent(out) :: maxglobseta,maxglobip
real, dimension(ifull+iextra), intent(inout) :: neta,ipice
real, dimension(ifull), intent(in) :: sue,svn,suw,svs
real, dimension(ifull), intent(in) :: pue,pvn,puw,pvs
real, dimension(ifull), intent(in) :: que,qvn,quw,qvs
real, dimension(ifull), intent(in) :: pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps
real, dimension(ifull+iextra), intent(in) :: ibu,ibv,idu,idv
real, dimension(ifull+iextra), intent(in) :: niu,niv
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull,2) :: zz,zzn,zzs,zze,zzw,rhs
real, dimension(ifull) :: yy,yyn,yys,yye,yyw
real, dimension(ifull) :: hh

call mgsor_init

! zz*(DIV^2 neta) + yy*neta*(DIV^2 neta) + hh*neta = rhs

! ocean
zzn(:,1) =  (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*ddv(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzs(:,1) = -(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*ddv(isv)    *em(1:ifull)*em(1:ifull)/ds
zze(:,1) =  (1.+ocneps)*0.5*dt*(que+sue+pue)*ddu(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzw(:,1) = -(1.+ocneps)*0.5*dt*(quw+suw+puw)*ddu(iwu)    *em(1:ifull)*em(1:ifull)/ds
zz(:,1)  =  (1.+ocneps)*0.5*dt*(qdivb+sdivb+pdivb)

yyn =  (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*em(1:ifull)*em(1:ifull)/ds
yys = -(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*em(1:ifull)*em(1:ifull)/ds
yye =  (1.+ocneps)*0.5*dt*(que+sue+pue)*em(1:ifull)*em(1:ifull)/ds
yyw = -(1.+ocneps)*0.5*dt*(quw+suw+puw)*em(1:ifull)*em(1:ifull)/ds
yy  =  (1.+ocneps)*0.5*dt*(qdiv+sdiv+pdiv)

hh      = 1.+(1.+ocneps)*0.5*dt*(odiv                                                                                           &
          +(pvn*dd(in)-pvs*dd(is)+pue*dd(ie)-puw*dd(iw))*em(1:ifull)*em(1:ifull)/ds+pdiv*dd(1:ifull))
rhs(:,1) = xps(1:ifull)-(1.+ocneps)*0.5*dt*(odivb                                                                               &
          +(pvn*ddv(1:ifull)*dd(in)-pvs*ddv(isv)*dd(is)+pue*ddu(1:ifull)*dd(ie)-puw*ddu(iwu)*dd(iw))*em(1:ifull)*em(1:ifull)/ds &
          +pdivb*dd(1:ifull))

! zz*(DIV^2 ipice) = rhs

! ice
zzn(:,2) = (-idv(1:ifull)*0.5-ibv(1:ifull))/ds
zzs(:,2) = ( idv(isv)*0.5    -ibv(isv)    )/ds
zze(:,2) = (-idu(1:ifull)*0.5-ibu(1:ifull))/ds
zzw(:,2) = ( idu(iwu)*0.5    -ibu(iwu)    )/ds
zz(:,2)  = (-0.5*(idu(1:ifull)-idu(iwu)+idv(1:ifull)-idv(isv))+ibu(1:ifull)+ibu(iwu)+ibv(1:ifull)+ibv(isv))/ds

rhs(:,2) = min(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu)+niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv),0.)

! Ensure that zero is a valid solution for ice free grid points
where ( zz(:,2) >= 0. )
  zz(:,2)  = -dt/(ds*10.) ! 10 is the minimum imass
  zzn(:,2) = 0.
  zzs(:,2) = 0.
  zze(:,2) = 0.
  zzw(:,2) = 0.
  rhs(:,2) = 0.
end where

call mgmlo(neta,ipice,yy,yyn,yys,yye,yyw,zz,zzn,zzs,zze,zzw,hh,rhs,tol,itol,totits,maxglobseta,maxglobip, &
           ipmax,ee,dd)
itc = totits

return
end subroutine mlomg

end module mlodynamics
