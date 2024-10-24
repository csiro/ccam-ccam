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
    
! These subroutines handle dynamics for the Mixed-Layer-Ocean model

! - Ocean dynamics
! - Ice dynamics

! This module links to mlo.f90 which solves for 1D physics of ocean
! and ice processes, as well as mlodiffg.f90 for horizontal
! diffusion.  Currently the code assumes the hydrostatic
! approximation which is reasonably valid to 1km resolution.

! Ocean and sea-ice dynamics are based on the R-grid used by CCAM.
! z* coordinates are used by the ocean to improve density gradients.
! Flexible nudging options are used for error correction (see
! nesting.f90).

module mlodynamics

use mlo_ctrl                                ! Ocean physics control layer
use mlodepts                                ! Ocean depature points
use mlodiffg                                ! Ocean dynamics horizontal diffusion
use mlodynamicsarrays_m                     ! Ocean dynamics data
use mloints                                 ! Ocean horizontal advection
use mlostag                                 ! Ocean reversible staggering
use mlovadvtvd                              ! Ocean vertical advection

implicit none

private
public mlohadv,mlodyninit
public ocneps,ocnepr,nxtrrho
public usetide,mlojacobi,mlomfix,nodrift,mlodps,mlo_bs,mlointschf
public mloiceadv

public mlodiffusion,mlodiffusion_work
public mlodiff,ocnsmag,ocnlap,mlodiff_numits

public ee, eeu, eev
public dd, ddu, ddv
public stwgtu, stwgtv
public gosig, gosigh, godsig
public godsigu, godsigv, gosighu, gosighv
public oldu1, oldu2, oldv1, oldv2
public ipice
public olddrhobardxu, olddrhobardyu, olddrhobardxv, olddrhobardyv
public oldrhobar_dash, oldrhobaru_dash, oldrhobarv_dash
public mlodynamicsarrays_init, mlodynamicsarrays_end

public mlostaguv, mlounstaguv
public mstagf, koff, nstagoffmlo

public mlovadv, mlontvd

integer, save      :: usetide     = 1       ! tidal forcing (0=off, 1=on)
integer, parameter :: icemode     = 2       ! ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, save      :: nxtrrho     = 1       ! estimate rho at t+1 (0=off, 1=on)
integer, save      :: mlojacobi   = 7       ! density gradient method (0=off, 6=dDdx=0, 7=AC2003)
integer, save      :: nodrift     = 0       ! remove drift from eta (0=off, 1=on)
integer, save      :: mlomfix     = 2       ! conserve T & S (0=off, 1=no free surface, 2=free surface)
integer, save      :: mlodps      = 1       ! include atmosphere pressure gradient terms (0=off, 1=on)
integer, save      :: mlo_bs      = 3       ! category to apply B-S advection (1=all, 2=mps, 3=uv, 4=t, 5=sal)
integer, save      :: mlointschf  = 2       ! direction for interpolation at depature points (0=off, 2=usual)
integer, save      :: mloiceadv   = 0       ! ice advection method (0=Upwind-orig, 1=Upwind-revised)
real, parameter :: rhosn          = 330.    ! density snow (kg m^-3)
real, parameter :: rhoic          = 900.    ! density ice  (kg m^-3)
real, parameter :: grav           = 9.80616 ! gravitational constant (m s^-2)
real, save      :: ocneps         = 0.2     ! lagrangian off-centring term for momentum
real, save      :: ocnepr         = 1.      ! lagrangian off-centring term for density gradient
real, parameter :: maxicefrac     = 0.999   ! maximum ice fraction
real, parameter :: tol            = 5.E-5   ! tolerance for SOR solver (water)
real, parameter :: itol           = 1.E1    ! tolerance for SOR solver (ice)
complex, save   :: emsum                    ! conservation reference value

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises mlo dynamical arrays
!
subroutine mlodyninit

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m
use soil_m
use sumdd_m

implicit none

integer ii, iq
real, dimension(ifull) :: odum
real, dimension(ifull,0:wlev) :: dephl
real, dimension(ifull,wlev) :: dep,dz
real, dimension(3*wlev) :: dumz,gdumz
real, dimension(ifull+iextra,wlev,2) :: dumf
logical, dimension(ifull+iextra,wlev) :: wtr
complex lsum

! calculate depths
allocate( dd(ifull+iextra) )
dd(:) = 0.
dep(:,:) = 0.
dz(:,:) = 0.
do ii = 1,wlev
  call mloexpdep("depth_fl",dep(:,ii),ii,0)
  call mloexpdep("depth_dz_fl",dz(:,ii),ii,0)
end do
dephl(:,0) = 0.
dephl(:,1) = dz(:,1)
do ii = 2,wlev
  dephl(:,ii) = dephl(:,ii-1) + dz(:,ii)
end do
dd(1:ifull) = max(dephl(:,wlev), 1.E-8)
call bounds(dd,nrows=2)

! prep land-sea mask
allocate( ee(ifull+iextra,wlev) )
allocate( eeu(ifull+iextra,wlev), eev(ifull+iextra,wlev) )
ee  = 0.
eeu = 0.
eev = 0. 
do ii = 1,wlev
  where ( dz(:,ii)>1.e-4 .and. .not.land(1:ifull) )
    ee(1:ifull,ii) = 1.
  end where
end do  
call bounds(ee,nrows=2)
do ii = 1,wlev
  where ( ee(:,ii)>1.5 .or. ee(:,ii)<=0.5 )
    ee(:,ii) = 0.
  end where 
  eeu(1:ifull,ii) = ee(1:ifull,ii)*ee(ie,ii)
  eev(1:ifull,ii) = ee(1:ifull,ii)*ee(in,ii)
end do  
call boundsuv(eeu,eev)

wtr = abs(ee-1.)<1.e-20

! Calculate staggered depth arrays
allocate( ddu(ifull+iextra), ddv(ifull+iextra) )
ddu = 0.
ddv = 0.
if ( mlosigma>=0 .and. mlosigma<=3 ) then
  write(6,*) "ERROR: mlosigma levels are not supported with mlosigma = ",mlosigma
  call ccmpi_abort(-1)
else
  ! z* levels
  ddu(1:ifull) = min(dd(1:ifull),dd(ie))
  ddv(1:ifull) = min(dd(1:ifull),dd(in))
end if
where ( abs(eeu(1:ifull,1))<1.e-20 )
  ddu(1:ifull) = 1.E-8
end where
where ( abs(eev(1:ifull,1))<1.e-20 )
  ddv(1:ifull) = 1.E-8
end where
call boundsuv(ddu,ddv)


! allocate memory for mlo dynamics arrays
call mlodynamicsarrays_init(ifull,iextra,wlev)

if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
  ! Precompute weights for calculating staggered gradients
  allocate(stwgtu(ifull,wlev),stwgtv(ifull,wlev))
  stwgtu = 0.
  stwgtv = 0.
  do ii = 1,wlev
    where ( wtr(in,ii) .and. wtr(ine,ii) .and. wtr(is,ii) .and. wtr(ise,ii) .and. &
            wtr(1:ifull,ii) .and. wtr(ie,ii) )
      stwgtu(1:ifull,ii) = 1.
    end where
    where ( wtr(ie,ii) .and. wtr(ien,ii) .and. wtr(iw,ii) .and. wtr(iwn,ii) .and. &
            wtr(1:ifull,ii) .and. wtr(in,ii) )
      stwgtv(1:ifull,ii) = 1.
    end where
            
  end do      
end if ! abs(nmlo)>=3.and.abs(nmlo)<=9

! vertical levels (using sigma notation)
allocate(gosig(1:ifull+iextra,wlev),gosigh(1:ifull,wlev),godsig(1:ifull+iextra,wlev))
allocate(godsigu(1:ifull+iextra,wlev),godsigv(1:ifull+iextra,wlev))
allocate(gosighu(1:ifull,wlev),gosighv(1:ifull,wlev))
gosig = 0.
gosigh = 0.
godsig = 0.
godsigu = 0.
godsigv = 0.
do ii = 1,wlev
  gosig(1:ifull,ii)  = max(dep(1:ifull,ii),1.e-8)/dd(1:ifull)
  gosigh(1:ifull,ii) = max(dephl(1:ifull,ii),1.e-8)/dd(1:ifull)
  godsig(1:ifull,ii) = dz(1:ifull,ii)/dd(1:ifull)*ee(1:ifull,ii)
end do
if ( mlosigma>=0 .and. mlosigma<=3 ) then
  write(6,*) "ERROR: mlo sigma levels not supported with mlosigma = ",mlosigma
  call ccmpi_abort(-1)
else
  ! z* levels
  call bounds(gosig,nrows=2)
  call bounds(godsig,nrows=2)
  do ii = 1,wlev
    godsigu(1:ifull,ii) = min(dd(1:ifull)*godsig(1:ifull,ii),dd(ie)*godsig(ie,ii))/ddu(1:ifull)
    godsigv(1:ifull,ii) = min(dd(1:ifull)*godsig(1:ifull,ii),dd(in)*godsig(in,ii))/ddv(1:ifull)
  end do  
end if
call boundsuv(godsigu,godsigv)
gosighu(:,1) = godsigu(1:ifull,1)
gosighv(:,1) = godsigv(1:ifull,1)
do ii = 2,wlev
  gosighu(1:ifull,ii) = gosighu(1:ifull,ii-1) + godsigu(1:ifull,ii)
  gosighv(1:ifull,ii) = gosighv(1:ifull,ii-1) + godsigv(1:ifull,ii)
end do
do ii = 1,wlev
  gosighu(1:ifull,ii) = max(gosighu(1:ifull,ii),1.e-8/ddu(1:ifull))
  gosighv(1:ifull,ii) = max(gosighv(1:ifull,ii),1.e-8/ddv(1:ifull))
end do

! store correction for nodrift==1
odum = ee(1:ifull,1)*em(1:ifull)**2
lsum = cmplx( 0., 0. )
call drpdr_local(odum,lsum)
call ccmpi_allreduce(lsum,emsum,"sumdr",comm_world)

return
end subroutine mlodyninit

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
use newmpar_m
use nharrs_m, only : lrestart
use parm_m
use parmdyn_m
use soil_m
use soilsnow_m
use sumdd_m
use vecsuv_m

implicit none

integer ii,totits
integer jyear,jmonth,jday,jhour,jmin,mins
integer tyear,jstart, iq, mspec_mlo, mspeca_mlo
integer, dimension(ifull,wlev) :: nface
real maxglobseta,maxglobip,hdt,dtin_mlo
real alph_p, delta, rone
real, dimension(2) :: delpos, delneg
real, dimension(ifull+iextra) :: neta,pice,imass
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv
real, dimension(ifull+iextra) :: sou,sov,snu,snv,squ,sqv
real, dimension(ifull+iextra) :: tide, depdum_rho
real, dimension(ifull+iextra) :: ipmax, spnet
real, dimension(ifull+iextra) :: ibu, ibv
real, dimension(ifull+iextra) :: bb, ff, bb3u, bb4v, ibb, ibb3u, ibb4v
real, dimension(ifull) :: ddddxu, ddddxv, ddddyu, ddddyv
real, dimension(ifull) :: i_u,i_v,i_sto
real, dimension(ifull) :: w_e,xps,oeu,oev,ipu,ipv
real, dimension(ifull) :: tnu,tsu,tev,twv
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyv
real, dimension(ifull) :: dipdxu,dipdyv
real, dimension(ifull) :: odum
real, dimension(ifull) :: dnetadx, dnetady, ddddx, ddddy
real, dimension(ifull) :: sdiv, imu, imv
real, dimension(ifull) :: oeu_iwu, oev_isv
real, dimension(ifull) :: bu, bv, cu, cv
real, dimension(ifull+iextra,wlev,3) :: cou
real, dimension(ifull+iextra,wlev) :: cc
real, dimension(ifull+iextra,wlev) :: eou, eov, ccu, ccv
real, dimension(ifull+iextra,wlev) :: nu, nv, nt, ns, mps
real, dimension(ifull+iextra,9) :: data_c, data_d
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull,wlev+1) :: tau, tav, ttau, ttav
real, dimension(ifull,wlev) :: u_tstar, v_tstar
real, dimension(ifull,wlev) :: cc3u, cc4v
real, dimension(ifull,wlev) :: w_u, w_v, w_t, w_s
real, dimension(ifull,wlev) :: nuh, nvh, xg, yg, uau, uav
real, dimension(ifull,wlev) :: kku, kkv, oou, oov, mmu, mmv
real, dimension(ifull,wlev) :: drhobardxu, drhobardyu, drhobardxv, drhobardyv
real, dimension(ifull,wlev) :: rhou_dash, rhov_dash, rho_dash
real, dimension(ifull,wlev) :: depdum,dzdum,mfixdum
real, dimension(ifull,wlev) :: workdata, workdata2, worku, workv
real, dimension(ifull,wlev) :: w_ocn
real, dimension(ifull,wlev) :: gosigu, gosigv
real, dimension(ifull,0:wlev) :: nw
real, dimension(ifull,4) :: i_it
real, dimension(ifull,3) :: gamm
real(kind=8), dimension(ifull,wlev) :: x3d,y3d,z3d
logical lleap, bs_test
logical, dimension(ifull+iextra,wlev) :: wtr
complex lsum, gsum

! Vertical coordinates are defined as:
!  mlosigma=4,5,6,7
!   z* = D*(z+eta)/(D+eta) = sig*D
!   from Adcroft and Campin 2003 (see also MPAS)

! Instead of packing ocean points, we keep the original
! grid and define a land/sea mask.  This allows us to
! use the same index and map factor arrays from the
! atmospheric dynamical core.

! save time-step
dtin_mlo = dt
mspeca_mlo = 1

! Define land/sea mask
wtr(:,:) = ee(:,:)>0.5

! Default values
w_t      = 293.16-wrtemp ! potential water temperature delta at tau=t
w_s      = 34.72         ! water salinity at tau=t
w_u      = 0.            ! u component of water current at tau=t
w_v      = 0.            ! v component of water current at tau=t
w_e      = 0.            ! free surface height at tau=t
i_it     = 273.16        ! ice temperature
i_sto    = 0.            ! ice brine storage
i_u      = 0.            ! u component of ice velocity
i_v      = 0.            ! v component of ice velocity
nw       = 0.            ! water vertical velocity
pice     = 0.            ! ice pressure for cavitating fluid
imass    = 0.            ! ice mass
tide     = 0.            ! tidal forcing
nt       = 293.16-wrtemp ! new water temperature
ns       = 0.            ! new water salinity
nu       = 0.            ! new u component of water current
nv       = 0.            ! new v component of water current
neta     = 0.            ! new free surface height
niu      = 0.            ! new u component of ice velocity
niv      = 0.            ! new v component of ice velocity
nfracice = 0.            ! new ice fraction
ndic     = 0.            ! new ice thickness 
ndsn     = 0.            ! new ice snow depth 
cou      = 0.            ! working array
data_c   = 0.            ! working array
data_d   = 0.            ! working array
eou      = 0.            ! working array
eov      = 0.            ! working array

! EXPORT WATER AND ICE DATA FROM MLO ------------------------------------------
call mloexport3d(0,w_t,0)
call mloexport3d(1,w_s,0)
call mloexport3d(2,w_u,0)
call mloexport3d(3,w_v,0)
call mloexport("eta",w_e,0,0)
call mloexpice("tsurf",i_it(:,1),0)
call mloexpice("temp0",i_it(:,2),0)
call mloexpice("temp1",i_it(:,3),0)
call mloexpice("temp2",i_it(:,4),0)
call mloexpice("fracice",nfracice,0)
call mloexpice("thick",ndic,0)
call mloexpice("snowd",ndsn,0)
call mloexpice("store",i_sto,0)
call mloexpice("u",i_u,0)
call mloexpice("v",i_v,0)
where ( wtr(1:ifull,1) )
  fracice(1:ifull) = nfracice(1:ifull)  
  sicedep(1:ifull) = ndic(1:ifull)
  snowd(1:ifull)   = ndsn(1:ifull)*1000.
end where  

! estimate tidal forcing
if ( usetide==1 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  if ( leap==0 ) then ! 365 day calendar
    jstart = 365*(jyear-1900)  
  else if ( leap==1 ) then ! 365/366 day calendar
    jstart = 0
    if ( jyear>1900 ) then
      do tyear = 1900,jyear-1
        call mloleap(tyear,lleap)
        if ( lleap ) jstart = jstart + 1
        jstart = jstart + 365
      end do
    else if ( jyear<1900 ) then
      do tyear = 1899,jyear,-1
        call mloleap(tyear,lleap)
        if ( lleap ) jstart = jstart - 1
        jstart = jstart - 365
      end do
    end if
  else if ( leap==2 ) then ! 360 day calendar  
    jstart = 360*(jyear-1900)
  else
    jstart = 0
    write(6,*) "ERROR: Unknown option leap = ",leap
    call ccmpi_abort(-1)
  end if    
  mins = mins + 720 ! base time is 12Z 31 Dec 1899
  call mlotide(tide,rlongg,rlatt,mins,jstart)
else if ( usetide==0 ) then  
  tide = 0.
else
  write(6,*) "ERROR: Unknown option usetide = ",usetide
  call ccmpi_abort(-1)
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


workdata = nt(1:ifull,:)
workdata2 = ns(1:ifull,:)
worku = nu(1:ifull,:)
workv = nv(1:ifull,:)
call mlocheck("start of mlodynamics",water_temp=workdata,water_sal=workdata2,water_u=worku, &
              water_v=workv,ice_tsurf=nit(1:ifull,1),ice_u=niu(1:ifull),ice_v=niv(1:ifull), &
              ice_thick=ndic(1:ifull))


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
  ! Limit minimum ice mass for ice velocity calculation.  Hence, we can solve for the ice velocity at
  ! grid points where the ice is not yet present.  
  imass(1:ifull) = max( imass(1:ifull), minicemass )
  pice(1:ifull) = ps(1:ifull)                                 ! pressure due to atmosphere at top of water column (unstaggered at t)
  if ( usepice==1 ) then
    pice(1:ifull) = pice(1:ifull) + grav*nfracice(1:ifull)*imass(1:ifull) ! include ice in surface pressure
  end if
  ! maximum pressure for ice
  select case(icemode)
    case(2)
      ! cavitating fluid
      ipmax(1:ifull) = 27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))*ee(1:ifull,1)
    case(1)
      ! incompressible fluid
      ipmax(1:ifull) = 9.E9*ee(1:ifull,1)
    case DEFAULT
      ! free drift
      ipmax(1:ifull) = 0.
  end select

  ! update scalar bounds and various gradients (aggregate fields for MPI)
  data_c(1:ifull,1) = neta(1:ifull)
  data_c(1:ifull,2) = pice(1:ifull)
  data_c(1:ifull,3) = tide(1:ifull)
  data_c(1:ifull,4) = imass(1:ifull)
  data_c(1:ifull,5) = ipmax(1:ifull)
  call bounds(data_c(:,1:5),corner=.true.)
  neta(ifull+1:ifull+iextra)  = data_c(ifull+1:ifull+iextra,1)
  pice(ifull+1:ifull+iextra)  = data_c(ifull+1:ifull+iextra,2)
  tide(ifull+1:ifull+iextra)  = data_c(ifull+1:ifull+iextra,3)
  imass(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,4)
  ipmax(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,5)

  ! surface pressure
  if ( mlodps==1 ) then
    do iq = 1,ifull
      tnu(iq) = 0.5*( pice(in(iq)) + pice(ine(iq)) )
      tsu(iq) = 0.5*( pice(is(iq)) + pice(ise(iq)) )
      tev(iq) = 0.5*( pice(ie(iq)) + pice(ien(iq)) )
      twv(iq) = 0.5*( pice(iw(iq)) + pice(iwn(iq)) )
    end do  
    dpsdxu = (pice(ie)-pice(1:ifull))*emu(1:ifull)/ds ! staggered at time t
    dpsdyu = 0.5*(tnu-tsu)*emu(1:ifull)/ds
    dpsdxv = 0.5*(tev-twv)*emv(1:ifull)/ds
    dpsdyv = (pice(in)-pice(1:ifull))*emv(1:ifull)/ds
  else if ( mlodps==0 ) then
    dpsdxu = 0.
    dpsdyu = 0.
    dpsdxv = 0.
    dpsdyv = 0.
  else
    write(6,*) "ERROR: Unknown option mlodps =",mlodps
    call ccmpi_abort(-1)
  end if
 
  ! tides
  if ( usetide==1 ) then
    do iq = 1,ifull
      tnu(iq) = 0.5*( tide(in(iq)) + tide(ine(iq)) )
      tsu(iq) = 0.5*( tide(is(iq)) + tide(ise(iq)) )
      tev(iq) = 0.5*( tide(ie(iq)) + tide(ien(iq)) )
      twv(iq) = 0.5*( tide(iw(iq)) + tide(iwn(iq)) )
    end do  
    dttdxu = (tide(ie)-tide(1:ifull))*emu(1:ifull)/ds ! staggered
    dttdyu = 0.5*(tnu-tsu)*emu(1:ifull)/ds
    dttdxv = 0.5*(tev-twv)*emv(1:ifull)/ds
    dttdyv = (tide(in)-tide(1:ifull))*emv(1:ifull)/ds
  else if ( usetide==0 ) then
    dttdxu = 0.
    dttdyu = 0.
    dttdxv = 0.
    dttdyv = 0.
  else    
    write(6,*) "ERROR: Unknown option usetide = ",usetide
    call ccmpi_abort(-1)
  end if    

  ! for 5-point stencil
  bu(1:ifull) = -eeu(1:ifull,1)*(1.+ocneps)*0.5*dt/(1.+((1.+ocneps)*0.5*dt*fu(1:ifull))**2)
  bv(1:ifull) = -eev(1:ifull,1)*(1.+ocneps)*0.5*dt/(1.+((1.+ocneps)*0.5*dt*fv(1:ifull))**2)
  cu(1:ifull) =  bu(1:ifull)*(1.+ocneps)*0.5*dt*fu(1:ifull)
  cv(1:ifull) = -bv(1:ifull)*(1.+ocneps)*0.5*dt*fv(1:ifull)

  ! Calculate adjusted depths and thicknesses
  do ii = 1,wlev
    do iq = 1,ifull
      depdum(iq,ii) = gosig(iq,ii)*max( dd(iq)+neta(iq), minwater ) !-neta(iq)
      dzdum(iq,ii)  = godsig(iq,ii)*max( dd(iq)+neta(iq), minwater )
      gosigu(iq,ii) = min(gosig(iq,ii)*dd(iq),gosig(ie(iq),ii)*dd(ie(iq)))/max(ddu(iq),1.e-8)
      gosigv(iq,ii) = min(gosig(iq,ii)*dd(iq),gosig(in(iq),ii)*dd(in(iq)))/max(ddv(iq),1.e-8)
    end do
  end do

  
  call START_LOG(watereos_begin)

  ! Calculate normalised density coeffs dalpha and dbeta (unstaggered at time t)
  ! (Assume free surface correction is small so that changes in the compression 
  ! effect due to neta can be neglected.  Consequently, the neta dependence is 
  ! separable in the iterative loop)
  call bounds(nt,corner=.true.)
  call bounds(ns,corner=.true.)

  ! rho(x) = wrtrho + rho_dash(x), where wrtrho is a constant
  
  ! Calculate normalised density gradients
  ! method 2: Use potential temperature and salinity Jacobians (see Shchepetkin and McWilliams 2003)
  call tsjacobi(nt,ns,pice,drhobardxu,drhobardyu,drhobardxv,drhobardyv, &
                rho_dash,rhou_dash,rhov_dash)
  
  call END_LOG(watereos_end)


  ! ADVECT WATER AND ICE ----------------------------------------------
  ! Water currents are advected using semi-Lagrangian advection
  ! based on McGregor's CCAM advection routines.
  

  ! Define horizontal transport
  ! dH(phi)/dt = d(phi)/dt + u*d(phi)/dx + v*d(phi)/dy + w*d(phi)/dz
  
  ! Continuity equation
  ! 1/(neta+D) dH(neta+D)/dt + du/dx + dv/dy + dw/dz = 0
  ! dH(neta+D)/dt + (neta+D)*(du/dx+dv/dy) + (neta+D)*dw/dz = 0
  ! dH(neta)/dt + neta*(du/dx+dv/dy) + d(D*u)/dx + d(D*v)/dy + dw/dsig = 0  ! where D is now included in flux form
  
  ! Lagrangian version of the contunity equation with D included in flux form
  !   [neta + (1+eps)*0.5*dt*(neta*(du/dx+dv/dy)+d(D*u)/dx+d(D*v)/dy+dw/dsig)]^(t+1)
  ! = [neta - (1-eps)*0.5*dt*(neta*(du/dx+dv/dy)+d(D*u)/dx+d(D*v)/dy+dw/dsig)]^(t*)
  ! = mps/dsig
  
  ! Vertical velocity is then (sum_0_sig is the integral from 0 to sig)

  ! dw/dsig = -(d(u*(neta+D))/dx+d(v*(neta+D))/dy) - d(neta)/dt
  ! d(neta)/dt = - sum_0_1(d(u*(neta+D))/dx+d(v*(neta+D))/dy,dsig)
  ! w = -sum_0_sig(d(u*(neta+D))/dx+d(v*(neta+D))/dy,dsig)
  !    + sig*sum_0_1(d(u*(neta+D))/dx+d(v*(neta+D))/dy,dsig)

  ! Advect continuity equation to tstar
  ! Calculate velocity on C-grid for consistancy with iterative free surface calculation
  ! nw is at half levels with nw(:,1) at the bottom of level 1
  ! positive nw is moving downwards to the ocean floor
  ! Assume Boussinesq fluid and treat density in continuity equation as constant
  ! true vertical velocity = nw-u*((1-sig)*deta/dx-sig*d(dd)/dx)-v*((1-sig)*deta/dy-sig*d(dd)/dy)-(1-sig)*deta/dt
  call mlostaguv(nu(:,1:wlev),nv(:,1:wlev),eou(:,1:wlev),eov(:,1:wlev))
  call boundsuv(eou(:,1:wlev),eov(:,1:wlev),stag=-9)
  ! vertically integrate currents
  ccu(:,1) = eou(:,1)*godsigu(:,1)
  ccv(:,1) = eov(:,1)*godsigv(:,1)
  do ii = 2,wlev
    ccu(:,ii) = ccu(:,ii-1) + eou(:,ii)*godsigu(:,ii)
    ccv(:,ii) = ccv(:,ii-1) + eov(:,ii)*godsigv(:,ii)
  end do
  ! surface height at staggered coordinate
  oeu(1:ifull) = 0.5*(neta(ie)+neta(1:ifull))*eeu(1:ifull,1)
  oev(1:ifull) = 0.5*(neta(in)+neta(1:ifull))*eev(1:ifull,1)
  ! trick to avoid additional call to bounds
  oeu_iwu(1:ifull) = 0.5*(neta(iw)+neta(1:ifull))*eeu(iwu,1)
  oev_isv(1:ifull) = 0.5*(neta(is)+neta(1:ifull))*eev(isv,1)
  ! calculate divergence for use with vertical velocity
  sdiv(:) = (ccu(1:ifull,wlev)*(ddu(1:ifull)+oeu(1:ifull))/emu(1:ifull)    &
            -ccu(iwu,wlev)*(ddu(iwu)+oeu_iwu)/emu(iwu)                     &
            +ccv(1:ifull,wlev)*(ddv(1:ifull)+oev(1:ifull))/emv(1:ifull)    &
            -ccv(isv,wlev)*(ddv(isv)+oev_isv)/emv(isv))*em(1:ifull)**2/ds

  ! calculate gradients for neta and dd
  detadxu = (neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
  detadyv = (neta(in)-neta(1:ifull))*emv(1:ifull)/ds
  dnetadx = (oeu(1:ifull)-oeu_iwu)*em(1:ifull)/ds
  dnetady = (oev(1:ifull)-oev_isv)*em(1:ifull)/ds
  do iq = 1,ifull
    tnu(iq) = 0.5*( dd(in(iq)) + dd(ine(iq)) )
    tsu(iq) = 0.5*( dd(is(iq)) + dd(ise(iq)) )
    tev(iq) = 0.5*( dd(ie(iq)) + dd(ien(iq)) )
    twv(iq) = 0.5*( dd(iw(iq)) + dd(iwn(iq)) )
  end do
  ddddxu = (dd(ie)-dd(1:ifull))*emu(1:ifull)/ds
  ddddyu = 0.5*(tnu-tsu)*emu(1:ifull)/ds
  ddddxv = 0.5*(tev-twv)*emv(1:ifull)/ds
  ddddyv = (dd(in)-dd(1:ifull))*emv(1:ifull)/ds
  ddddx = (ddu(1:ifull)-ddu(iwu))*em(1:ifull)/ds
  ddddy = (ddv(1:ifull)-ddv(isv))*em(1:ifull)/ds

  ! calculate vertical velocity (use flux form)
  do ii = 1,wlev-1
    nw(:,ii) = ee(1:ifull,ii)*ee(1:ifull,ii+1)*(sdiv(:)*gosigh(1:ifull,ii)  &
               - (ccu(1:ifull,ii)*(ddu(1:ifull)+oeu(1:ifull))/emu(1:ifull)  &
                 -ccu(iwu,ii)*(ddu(iwu)+oeu_iwu(:))/emu(iwu)                &
                 +ccv(1:ifull,ii)*(ddv(1:ifull)+oev(1:ifull))/emv(1:ifull)  &
                 -ccv(isv,ii)*(ddv(isv)+oev_isv(:))/emv(isv)                &
               )*em(1:ifull)**2/ds)
  end do
  !nw(:,0) = 0.
  !nw(:,wlev) = 0.

  ! compute contunity equation horizontal transport terms
  ! Include dsig terms before advection
  do ii = 1,wlev
    mps(1:ifull,ii) = neta(1:ifull)*godsig(1:ifull,ii) - (1.-ocneps)*0.5*dt                              &
                        *((eou(1:ifull,ii)*(ddu(1:ifull)+neta(1:ifull))/emu(1:ifull)*godsigu(1:ifull,ii) &
                          -eou(iwu,ii)*(ddu(iwu)+neta(1:ifull))/emu(iwu)*godsigu(iwu,ii)                 &
                          +eov(1:ifull,ii)*(ddv(1:ifull)+neta(1:ifull))/emv(1:ifull)*godsigv(1:ifull,ii) &
                          -eov(isv,ii)*(ddv(isv)+neta(1:ifull))/emv(isv)*godsigv(isv,ii)                 &
                         )*em(1:ifull)**2/ds + (nw(:,ii)-nw(:,ii-1)))
  end do

  ! calculate vertical velocity at full levels for diagnostic output
  do ii = 1,wlev
    w_ocn(:,ii) = ee(1:ifull,ii)*(0.5*(nw(:,ii-1)+nw(:,ii))                                         &
                - nu(1:ifull,ii)*((1.-gosig(1:ifull,ii))*dnetadx(:) - gosig(1:ifull,ii)*ddddx(:))   &
                - nv(1:ifull,ii)*((1.-gosig(1:ifull,ii))*dnetady(:) - gosig(1:ifull,ii)*ddddy(:)))  &
                - (1.-gosig(1:ifull,ii))*sdiv(:) ! sdiv=d(eta)/dt 
  end do


  ! ocean
  ! Prepare pressure gradient terms at time t and incorporate into velocity field
  do ii = 1,wlev
    tau(:,ii) = (1.-ocneps)*0.5*dt*grav*detadxu                                                           & ! staggered
              + (1.-ocnepr)*0.5*dt*grav*gosigu(:,ii)*(ddu(1:ifull)+oeu(1:ifull))*drhobardxu(:,ii)/wrtrho  &
              + (1.-ocnepr)*0.5*dt*dpsdxu/wrtrho + (1.-ocnepr)*0.5*dt*grav*dttdxu
    tav(:,ii) = (1.-ocneps)*0.5*dt*grav*detadyv                                                           &
              + (1.-ocnepr)*0.5*dt*grav*gosigv(:,ii)*(ddv(1:ifull)+oev(1:ifull))*drhobardyv(:,ii)/wrtrho  &
              + (1.-ocnepr)*0.5*dt*dpsdyv/wrtrho + (1.-ocnepr)*0.5*dt*grav*dttdyv
    tau(:,ii) = tau(:,ii)*eeu(1:ifull,ii)
    tav(:,ii) = tav(:,ii)*eev(1:ifull,ii)
  end do
  
  ! include rho_dash/wrtrho terms (i.e., sig*dD/dx and (1-sig)*d(eta)/dx )
  select case(mlojacobi)
    case(7)
      do ii = 1,wlev
        tau(:,ii) = tau(:,ii) + (1.-ocneps)*0.5*dt*grav*rhou_dash(:,ii)/wrtrho*detadxu*(1.-gosigu(:,ii)) &
                  + (1.-ocneps)*0.5*dt*grav*rhou_dash(:,ii)/wrtrho*ddddxu*gosigu(:,ii)*oeu(1:ifull)/ddu(1:ifull) ! staggered
        tav(:,ii) = tav(:,ii) + (1.-ocneps)*0.5*dt*grav*rhov_dash(:,ii)/wrtrho*detadyv*(1.-gosigv(:,ii)) &
                  + (1.-ocneps)*0.5*dt*grav*rhov_dash(:,ii)/wrtrho*ddddyv*gosigv(:,ii)*oev(1:ifull)/ddv(1:ifull)
        tau(:,ii) = tau(:,ii)*eeu(1:ifull,ii)
        tav(:,ii) = tav(:,ii)*eev(1:ifull,ii)
      end do     
  end select    
    
  ! ice
  !tau(:,wlev+1)=grav*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds ! staggered
  !tav(:,wlev+1)=grav*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
  call mlounstaguv(tau(:,1:wlev),tav(:,1:wlev),ttau(:,1:wlev),ttav(:,1:wlev),toff=1)
  ! ocean
  do ii = 1,wlev
    u_tstar(:,ii) = nu(1:ifull,ii) + (1.-ocneps)*0.5*dt*f(1:ifull)*nv(1:ifull,ii) - ttau(:,ii) ! unstaggered
    u_tstar(:,ii) = u_tstar(:,ii)*ee(1:ifull,ii)
    v_tstar(:,ii) = nv(1:ifull,ii) - (1.-ocneps)*0.5*dt*f(1:ifull)*nu(1:ifull,ii) - ttav(:,ii)
    v_tstar(:,ii) = v_tstar(:,ii)*ee(1:ifull,ii)
  end do
  ! ice
  snu(1:ifull) = i_u !-dt*ttau(:,wlev+1)
  snv(1:ifull) = i_v !-dt*ttav(:,wlev+1)
  

  ! Vertical advection (first call for 0.5*dt)
  call mlovadv(hdt,nw,u_tstar,v_tstar,ns,nt,mps,depdum,dzdum,ee,1)

  
  workdata = nt(1:ifull,:)
  workdata2 = ns(1:ifull,:)
  call mlocheck("first vertical advection",water_temp=workdata,water_sal=workdata2)

  ! Calculate depature points for horizontal advection
  do ii = 1,wlev
    ! Estimate currents at t+1/2 for semi-Lagrangian advection
    nuh(:,ii) = (15.*nu(1:ifull,ii)-10.*oldu1(:,ii)+3.*oldu2(:,ii))*ee(1:ifull,ii)/8. ! U at t+1/2
    nvh(:,ii) = (15.*nv(1:ifull,ii)-10.*oldv1(:,ii)+3.*oldv2(:,ii))*ee(1:ifull,ii)/8. ! V at t+1/2
    nuh(:,ii) = min( max( nuh(:,ii), -50. ), 50. )
    nvh(:,ii) = min( max( nvh(:,ii), -50. ), 50. )

    ! save arrays for extrapolating currents at next time-step
    oldu2(:,ii) = oldu1(:,ii)
    oldv2(:,ii) = oldv1(:,ii)
    oldu1(:,ii) = nu(1:ifull,ii)
    oldv1(:,ii) = nv(1:ifull,ii)
  end do

  call mlodeps(nuh,nvh,nface,xg,yg,x3d,y3d,z3d,wtr,mlointschf)
  
  do ii = 1,wlev
    ! Convert (u,v) to cartesian coordinates (U,V,W)  
    cou(1:ifull,ii,1) = ax(1:ifull)*u_tstar(:,ii) + bx(1:ifull)*v_tstar(:,ii)
    cou(1:ifull,ii,2) = ay(1:ifull)*u_tstar(:,ii) + by(1:ifull)*v_tstar(:,ii)
    cou(1:ifull,ii,3) = az(1:ifull)*u_tstar(:,ii) + bz(1:ifull)*v_tstar(:,ii)
  end do
  
  
#ifdef GPU
  !$acc data create(xg,yg,nface)
  !$acc update device(xg,yg,nface)
#endif
  
  
  ! Horizontal advection for U, V, W, T, continuity and S
  bs_test = mlo_bs<=2
  call mlob2ints_bs(cou(:,:,1:3),nface,xg,yg,wtr,0,bs_test,mlointschf)
  bs_test = mlo_bs<=3
  call mlob2ints_bs(nt,nface,xg,yg,wtr,0,bs_test,mlointschf)
  bs_test = mlo_bs<=1
  call mlob2ints_bs(mps,nface,xg,yg,wtr,2,bs_test,mlointschf)
  bs_test = mlo_bs<=4
  call mlob2ints_bs(ns,nface,xg,yg,wtr,1,bs_test,mlointschf)
  
  
#ifdef GPU
  !$acc end data
#endif

  
  ! Rotate vector to arrival point
  call mlorot(cou(:,:,1),cou(:,:,2),cou(:,:,3),x3d,y3d,z3d)
  
  do ii = 1,wlev
    ! Convert (U,V,W) back to conformal cubic coordinates  
    u_tstar(:,ii) = ax(1:ifull)*cou(1:ifull,ii,1) + ay(1:ifull)*cou(1:ifull,ii,2) + az(1:ifull)*cou(1:ifull,ii,3)
    u_tstar(:,ii) = u_tstar(:,ii)*ee(1:ifull,ii)    
    v_tstar(:,ii) = bx(1:ifull)*cou(1:ifull,ii,1) + by(1:ifull)*cou(1:ifull,ii,2) + bz(1:ifull)*cou(1:ifull,ii,3)
    v_tstar(:,ii) = v_tstar(:,ii)*ee(1:ifull,ii)
    ! Update T, S and continuity
    nt(1:ifull,ii) = max( nt(1:ifull,ii), -wrtemp )
    ! MJT notes - only advect salinity for salt water
    ns(1:ifull,ii) = max( ns(1:ifull,ii), 0. )
  end do

  workdata = nt(1:ifull,:)
  workdata2 = ns(1:ifull,:)
  call mlocheck("horizontal advection",water_temp=workdata,water_sal=workdata2)


  ! Vertical advection (second call for 0.5*dt)
  ! use explicit nw and depdum,dzdum from t=tau step (i.e., following JLM in CCAM atmospheric dynamics)
  ! Could use nuh and nvh to estimate nw at t+1/2, but also require an estimate of neta at t+1/2
  call mlovadv(hdt,nw,u_tstar,v_tstar,ns,nt,mps,depdum,dzdum,ee,2)


  workdata = nt(1:ifull,:)
  workdata2 = ns(1:ifull,:)
  call mlocheck("second vertical advection",water_temp=workdata,water_sal=workdata2)


  ! dsig terms are included before advection 
  xps(:) = mps(1:ifull,1)
  do ii = 2,wlev
    xps(:) = xps(:) + mps(1:ifull,ii)
  end do
  xps(:) = xps(:)*ee(1:ifull,1) ! only needed for sigma levels
  
  
  call START_LOG(watereos_begin)

  ! Approximate normalised density rhobar at t+1 (unstaggered, using T and S at t+1)
  if ( nxtrrho==1 ) then
    call bounds(nt,corner=.true.)
    call bounds(ns,corner=.true.)
    ! update normalised density gradients
    call tsjacobi(nt,ns,pice,drhobardxu,drhobardyu,drhobardxv,drhobardyv, &
                  rho_dash,rhou_dash,rhov_dash)
  end if

  call END_LOG(watereos_end)


  ! FREE SURFACE CALCULATION ----------------------------------------

  ! Precompute U,V current and integral terms at t+1
  ! ocean
  do ii = 1,wlev
    tau(:,ii) = u_tstar(:,ii) + (1.+ocneps)*0.5*dt*f(1:ifull)*v_tstar(:,ii) ! unstaggered
    tav(:,ii) = v_tstar(:,ii) - (1.+ocneps)*0.5*dt*f(1:ifull)*u_tstar(:,ii)
  end do
  ! ice
  tau(:,wlev+1) = snu(1:ifull) + dt*f(1:ifull)*snv(1:ifull) ! unstaggered
  tav(:,wlev+1) = snv(1:ifull) - dt*f(1:ifull)*snu(1:ifull)
  call mlostaguv(tau(:,1:wlev+1),tav(:,1:wlev+1),ttau(:,1:wlev+1),ttav(:,1:wlev+1))
  
  ! Set-up calculation of ocean and ice at t+1
  ! ocean
  ! ccu and ccv have the correct magnitude for the advected current
  odum = 1./(1.+((1.+ocneps)*0.5*dt*fu(1:ifull))**2)
  do ii = 1,wlev
    ccu(1:ifull,ii) = ttau(:,ii)*odum*eeu(1:ifull,ii) ! staggered
    !ccu(1:ifull,ii) = max(min( ccu(1:ifull,ii), mlomaxuv),-mlomaxuv)
  end do
  odum = 1./(1.+((1.+ocneps)*0.5*dt*fv(1:ifull))**2)
  do ii = 1,wlev
    ccv(1:ifull,ii) = ttav(:,ii)*odum*eev(1:ifull,ii) ! staggered
    !ccv(1:ifull,ii) = max(min( ccv(1:ifull,ii), mlomaxuv),-mlomaxuv)
  end do
  ! ice
  ! niu and niv hold the free drift solution (staggered).  Wind stress terms are updated in mlo.f90
  niu(1:ifull) = ttau(:,wlev+1)/(1.+(dt*fu(1:ifull))**2) ! staggered
  niv(1:ifull) = ttav(:,wlev+1)/(1.+(dt*fv(1:ifull))**2)
  

  ! calculate 5-point stencil for free surface
  ! (i.e., terms for free surface gradient contribution to currents)

  ! temporary terms to be updated later
  bb(:) = (1.+ocneps)*0.5*dt/(1.+((1.+ocneps)*0.5*dt*f(:))**2) ! unstaggered
  ff(:) = (1.+ocneps)*0.5*dt*f(:)

  do ii = 1,wlev
    cc(:,ii) = grav*bb(:)
  end do

  ! include rho_dash/wrtrho terms (i.e., sig*dD/dx and (1-sig)*d(eta)/dx )
  select case(mlojacobi)
    case(7)
      do ii = 1,wlev  
        cc(1:ifull,ii) = cc(1:ifull,ii) + grav*ee(1:ifull,ii)*bb(1:ifull) &
            *(1.-gosig(1:ifull,ii))*rho_dash(:,ii)/wrtrho
      end do
      call bounds(cc(:,1:wlev),corner=.true.)
  end select

  ! update cc, cc3u, cc4v terms with f terms
  do ii = 1,wlev
    cc3u(:,ii) = 0.
    cc4v(:,ii) = 0.
    do iq = 1,ifull
      tnu(iq) = 0.5*( cc(in(iq),ii)*ff(in(iq)) + cc(ine(iq),ii)*ff(ine(iq)) )
      tsu(iq) = 0.5*( cc(is(iq),ii)*ff(is(iq)) + cc(ise(iq),ii)*ff(ise(iq)) )
      tev(iq) = 0.5*( cc(ie(iq),ii)*ff(ie(iq)) + cc(ien(iq),ii)*ff(ien(iq)) )
      twv(iq) = 0.5*( cc(iw(iq),ii)*ff(iw(iq)) + cc(iwn(iq),ii)*ff(iwn(iq)) )
    end do
    cc3u(1:ifull,ii) = 0.5*stwgtu(1:ifull,ii)*( tnu-tsu )*emu(1:ifull)/ds
    cc4v(1:ifull,ii) = 0.5*stwgtv(1:ifull,ii)*( tev-twv )*emv(1:ifull)/ds
  end do

  ! ocean - calculate bb, bb3u, bb4v now possibly including rho_dash contribution to cc
  bb(:) = 0.
  bb3u(:) = 0.
  bb4v(:) = 0.
  do ii = 1,wlev
    bb(:) = bb(:) + cc(:,ii)*godsig(:,ii)
    bb3u(1:ifull) = bb3u(1:ifull) + cc3u(1:ifull,ii)*godsigu(1:ifull,ii)
    bb4v(1:ifull) = bb4v(1:ifull) + cc4v(1:ifull,ii)*godsigv(1:ifull,ii)
  end do

  ! ice - calculate ibb, ibb3u, ibb4v
  ibb(:) = ee(:,1)*dt/(max(imass(:),minicemass)*(1.+(dt*f(:))**2)) ! unstaggered
  ff(:) = dt*f(:)  
  do iq = 1,ifull
    tnu(iq) = 0.5*( ibb(in(iq))*ff(in(iq)) + ibb(ine(iq))*ff(ine(iq)) )
    tsu(iq) = 0.5*( ibb(is(iq))*ff(is(iq)) + ibb(ise(iq))*ff(ise(iq)) )
    tev(iq) = 0.5*( ibb(ie(iq))*ff(ie(iq)) + ibb(ien(iq))*ff(ien(iq)) )
    twv(iq) = 0.5*( ibb(iw(iq))*ff(iw(iq)) + ibb(iwn(iq))*ff(iwn(iq)) )
  end do
  ibb3u(1:ifull) = 0.5*stwgtu(1:ifull,1)*( tnu-tsu )*emu(1:ifull)/ds
  ibb4v(1:ifull) = 0.5*stwgtv(1:ifull,1)*( tev-twv )*emv(1:ifull)/ds
 
  
  rone = (1.+ocnepr)/(1.+ocneps)
  do ii = 1,wlev
    ! Create arrays to calcuate u and v at t+1, based on pressure gradient at t+1
    
    ! note explicit drhobardx, dttdx, dpsdx terms are treated as t+1 with ocneps=1.
    ! (not avergaed over advection)
      
    !   u^(t+1) - 0.5*(1+eps)*dt*f*v^(t+1) + 0.5*(1+eps)*dt*(grav/rho)*dpdx^(t+1)
    ! = u^t + 0.5*(1-eps)*dt*f*v^t - 0.5*(1-eps)*dt*(grav/rho)*dpdx^t
    ! = u^(t*)

    !   v^(t+1) + 0.5*(1+eps)*dt*f*u^(t+1) + 0.5*(1+eps)*dt*(grav/rho)*dpdy^(t+1)
    ! = v^t - 0.5*(1-eps)*dt*f*u^t - 0.5*(1-eps)*dt*(grav/rho)*dpdy^t
    ! = v^(t*)

    !   u^(t+1) + 0.5*(1+eps)*dt*f*( 0.5*(1+eps)*dt*f*u^(t+1) + 0.5*(1+eps)*dt*(grav/rho)*dpdy^(t+1) )
    !     + 0.5*(1+eps)*dt*(grav/rho)*dpdx^(t+1) = u^(t*) + 0.5*(1+eps)*dt*f*v^(t*)

    !   v^(t+1) + 0.5*(1+eps)*dt*f*( 0.5*(1+eps)*dt*f*v^(t+1) - 0.5*(1+eps)*dt*(grav/rho)*dpdx^(t+1) )
    !     + 0.5*(1+eps)*dt*(grav/rho)*dpdy^(t+1) = v^(t*) - 0.5*(1+eps)*dt*f*u^(t*)

    ! u^(t+1)*(1+(0.5*(1+eps)*dt*f)**2)
    !   + 0.5*(1+eps)*dt*(grav/rho)*dpdx^(t+1)
    !   + (0.5*(1+eps)*dt)**2*f*(grav/rho)*dpdy^(t+1)
    ! = u^(t*) + 0.5*(1+eps)*dt*f*v^(t*)

    ! v^(t+1)*(1+(0.5*(1+eps)*dt*f)**2)
    !   + 0.5*(1+eps)*dt*(grav/rho)*dpdy^(t+1)
    !   - (0.5*(1+eps)*dt)**2*f*(grav/rho)*dpdx^(t+1) 
    ! = v^(t*) - 0.5*(1+eps)*dt*f*u^(t*)
  
    ! ffu = (1.+ocneps)*0.5*dt*fu
    ! ffv = (1.+ocneps)*0.5*dt*fv     
    ! u^(t+1) = nu = ccu^(t*) + bu*dpdxu^(t+1)/wrtrho - bu*ffu*dpdyu^(t+1)/wrtrho (staggered)
    ! v^(t+1) = nv = ccv^(t*) + bv*dpdyv^(t+1)/wrtrho + bv*ffv*dpdxv^(t+1)/wrtrho (staggered)
    
    ! rho = wrtrho + rho_dash
    ! rhobar = (int_0_z rho dz)/z
      
    ! dPdz = grav*rho
    ! dPdz* = grav*rho*(1+eta/D)
    ! P = grav*rhobar*(eta+D)*(z*)/D = grav*rhobar*(z+eta) = grav*rhobar*sig*(eta+D)
    ! dPdx = grav*sig*(eta+D)*d(rhobar)dx + grav*rhobar*d(eta)dx
      
    ! dGdz = grav
    ! dGdz* = grav*(D+eta)/D
    ! G = grav*z
    ! dGdz* * dz*dx = grav*( dDdx*sig*eta/D + d(eta)dx*(1-sig) )
    
    ! dz*dx = dDdx*sig*eta/(D+eta) + detadx*D/(D+eta)
      
    ! -(1/wrtrho)*dPdz + (rho/wrtrho)*dGdz = 0   ! hydrostatic approximation

    ! -(1/wrtrho)*dPdx + dGdx
    ! = -grav*sig*(eta+D)*d(rhobar)dx/wrtrho - grav*d(eta)dx
    !   -grav*(rho_dash)*( dDdx*sig*eta/D + d(eta)dx*(1-sig) )
      
    ! (mlojacobi=6 neglects rho_dash/wrtrho and mlojacobi=7 includes rho_dash/wrtrho terms)
    
    ! We use a modified form of JLM's trick where:
    !   -B*F*dP/dy = - d(B*F*P)/dy + (d(B*F)/dy)*P
    !    B*F*dP/dx =   d(B*F*P)/dx - (d(B*F)/dx)*P
    !   B = dt/(1.+((1+eps)*0.5*dt*f)**2)*(1+(1-sig)*rho_dash/wrtrho)
    !   F = (1.+ocneps)*0.5*dt*f
    ! which produces a 5-point stencil when solving for neta and ip.      
      
    !dpdxu=dpsdxu+grav*wrtrho*dttdxu+grav*sig*(ddu+etau)*drhobardxu+grav*rhobar*detadxu
    !      +grav*rho_dash*sig*etau/ddu*ddddxu + grav*rho_dash*(1-sig)*detadxu
    !dpdyu=dpsdyu+grav*wrtrho*dttdyu+grav*sig*(ddu+etau)*drhobardyu+grav*rhobar*detadyu
    !      +grav*rho_dash*sig*etau/ddu*ddddyu + grav*rho_dash*(1-sig)*detadyu
    !dpdxv=dpsdxv+grav*wrtrho*dttdxv+grav*sig*(ddv+etav)*drhobardxv+grav*rhobar*detadxv
    !      +grav*rho_dash*sig*etav/ddv*ddddxv + grav*rho_dash*(1-sig)*detadxv
    !dpdyv=dpsdyv+grav*wrtrho*dttdyv+grav*sig*(ddv+etav)*drhobardyv+grav*rhobar*detadyv
    !      +grav*rho_dash*sig*etav/ddv*ddddyv + grav*rho_dash*(1-sig)*detadyv

    ! Create arrays for u and v at t+1 in terms of neta gradients
    
    ! nu = kku + oou*etau + mmu*detadxu + nnu*detadyu   (staggered)
    ! nv = kkv + oov*etav + mmv*detadyv + nnv*detadxv   (staggered)

    ! nu = kku + oou*etau + mmu*detadxu - d(cc*eta)dyu + cc3u*eta   (staggered)
    ! nv = kkv + oov*etav + mmv*detadyv + d(cc*eta)dxv - cc4v*eta   (staggered)      
     
    kku(:,ii) = ccu(1:ifull,ii) + bu*(dpsdxu/wrtrho+grav*dttdxu)*rone      &
              + cu*(dpsdyu/wrtrho+grav*dttdyu)*rone                        &
              + grav*gosigu(:,ii)*ddu(1:ifull)*(bu*drhobardxu(:,ii)/wrtrho &
              + cu*drhobardyu(:,ii)/wrtrho)*rone
    oou(:,ii) = grav*gosigu(:,ii)*(bu*drhobardxu(:,ii)/wrtrho              &
              + cu*drhobardyu(:,ii)/wrtrho)*rone
    mmu(:,ii) = grav*bu
    
    kkv(:,ii) = ccv(1:ifull,ii) + bv*(dpsdyv/wrtrho+grav*dttdyv)*rone      &
              + cv*(dpsdxv/wrtrho+grav*dttdxv)*rone                        &
              + grav*gosigv(:,ii)*ddv(1:ifull)*(bv*drhobardyv(:,ii)/wrtrho &
              + cv*drhobardxv(:,ii)/wrtrho)*rone   
    oov(:,ii) = grav*gosigv(:,ii)*(bv*drhobardyv(:,ii)/wrtrho              &
              + cv*drhobardxv(:,ii)/wrtrho)*rone
    mmv(:,ii) = grav*bv
    
    kku(:,ii) = kku(:,ii)*eeu(1:ifull,ii)
    oou(:,ii) = oou(:,ii)*eeu(1:ifull,ii)
    mmu(:,ii) = mmu(:,ii)*eeu(1:ifull,ii)
    kkv(:,ii) = kkv(:,ii)*eev(1:ifull,ii)
    oov(:,ii) = oov(:,ii)*eev(1:ifull,ii)
    mmv(:,ii) = mmv(:,ii)*eev(1:ifull,ii)
  end do
  
  ! include rho_dash/wrtrho terms (i.e., sig*(eta/D)*dD/dx and (1-sig)*d(eta)/dx )
  select case(mlojacobi)
    case(7)
      do ii = 1,wlev
        oou(:,ii) = oou(:,ii) + bu*grav*rhou_dash(:,ii)/wrtrho*ddddxu*gosigu(:,ii)/ddu(1:ifull) &
                  + cu*grav*rhou_dash(:,ii)/wrtrho*ddddyu*gosigu(:,ii)/ddu(1:ifull)
        mmu(:,ii) = mmu(:,ii) + grav*bu*(1.-gosigu(:,ii))*rhou_dash(:,ii)/wrtrho
        oov(:,ii) = oov(:,ii) + bv*grav*rhov_dash(:,ii)/wrtrho*ddddyv*gosigv(:,ii)/ddv(1:ifull) &
                  + cv*grav*rhov_dash(:,ii)/wrtrho*ddddxv*gosigv(:,ii)/ddv(1:ifull)
        mmv(:,ii) = mmv(:,ii) + grav*bv*(1.-gosigv(:,ii))*rhov_dash(:,ii)/wrtrho

        oou(:,ii) = oou(:,ii)*eeu(1:ifull,ii)
        mmu(:,ii) = mmu(:,ii)*eeu(1:ifull,ii) 
        oov(:,ii) = oov(:,ii)*eev(1:ifull,ii) 
        mmv(:,ii) = mmv(:,ii)*eev(1:ifull,ii) 
      end do    
  end select

  ! partial step correction
  ! -(1/wrtho)dPdz* + dGdz* = grav*rho_dash*(1+eta/D)
  if ( mlo_step==2 ) then
    do ii = 1,wlev
      do iq = 1,ifull
        tnu(iq) = 0.5*( gosig(in(iq),ii)*dd(in(iq)) + gosig(ine(iq),ii)*dd(ine(iq)) )
        tsu(iq) = 0.5*( gosig(is(iq),ii)*dd(is(iq)) + gosig(ise(iq),ii)*dd(ise(iq)) )
        tev(iq) = 0.5*( gosig(ie(iq),ii)*dd(ie(iq)) + gosig(ien(iq),ii)*dd(ien(iq)) )
        twv(iq) = 0.5*( gosig(iw(iq),ii)*dd(iw(iq)) + gosig(iwn(iq),ii)*dd(iwn(iq)) )
      end do
      kku(:,ii) = kku(:,ii) + bu(:)*grav*(dd(ie)*gosig(ie,ii)-dd(1:ifull)*gosig(1:ifull,ii))*emu(1:ifull)/ds              &
                              *rhou_dash(:,ii)/wrtrho*rone                                                                &
                            + cu(:)*grav*0.5*stwgtu(1:ifull,ii)*( tnu-tsu )*emu(1:ifull)/ds                               &
                              *rhou_dash(:,ii)/wrtrho*rone
      oou(:,ii) = oou(:,ii) + bu(:)*grav/ddu(1:ifull)*(dd(ie)*gosig(ie,ii)-dd(1:ifull)*gosig(1:ifull,ii))*emu(1:ifull)/ds &
                              *rhou_dash(:,ii)/wrtrho*rone                                                                &
                            + cu(:)*grav/ddu(1:ifull)*0.5*stwgtu(1:ifull,ii)*( tnu-tsu )*emu(1:ifull)/ds                  &
                              *rhou_dash(:,ii)/wrtrho*rone
      kkv(:,ii) = kkv(:,ii) + bv(:)*grav*(dd(in)*gosig(in,ii)-dd(1:ifull)*gosig(1:ifull,ii))*emv(1:ifull)/ds              &
                              *rhov_dash(:,ii)/wrtrho*rone                                                                &
                            + cv(:)*grav*0.5*stwgtv(1:ifull,ii)*( tev-twv )*emv(1:ifull)/ds                               &
                              *rhov_dash(:,ii)/wrtrho*rone
      oov(:,ii) = oov(:,ii) + bv(:)*grav/ddv(1:ifull)*(dd(in)*gosig(in,ii)-dd(1:ifull)*gosig(1:ifull,ii))*emv(1:ifull)/ds &
                              *rhov_dash(:,ii)/wrtrho*rone                                                                &
                            + cv(:)*grav/ddv(1:ifull)*0.5*stwgtv(1:ifull,ii)*( tev-twv )*emv(1:ifull)/ds                  &
                              *rhov_dash(:,ii)/wrtrho*rone

      kku(:,ii) = kku(:,ii)*eeu(1:ifull,ii)
      oou(:,ii) = oou(:,ii)*eeu(1:ifull,ii)
      kkv(:,ii) = kkv(:,ii)*eev(1:ifull,ii)
      oov(:,ii) = oov(:,ii)*eev(1:ifull,ii)
    end do
  end if
    
  
  ! Pre-integrate arrays for u and v at t+1 (i.e., for calculating net divergence at t+1)
  
  !sum nu dz = sou+snu*etau+squ*detadxu-d(bb*eta)dyu+bb3u*eta
  !sum nv dz = sov+snv*etav+sqv*detadyv+d(bb*eta)dxv-bb4v*eta

  sou = 0.
  sov = 0.
  snu = 0.
  snv = 0.
  squ = 0.
  sqv = 0.
  do ii = 1,wlev
    sou(1:ifull) = sou(1:ifull) + kku(:,ii)*godsigu(1:ifull,ii)
    snu(1:ifull) = snu(1:ifull) + oou(:,ii)*godsigu(1:ifull,ii)
    squ(1:ifull) = squ(1:ifull) + mmu(:,ii)*godsigu(1:ifull,ii)
    sov(1:ifull) = sov(1:ifull) + kkv(:,ii)*godsigv(1:ifull,ii)
    snv(1:ifull) = snv(1:ifull) + oov(:,ii)*godsigv(1:ifull,ii)
    sqv(1:ifull) = sqv(1:ifull) + mmv(:,ii)*godsigv(1:ifull,ii)
  end do
  
  
  ! calculate terms for ice velocity at t+1

  ! (staggered)
  ! niu(t+1) = niu + ibu*dipdx + ibu*dt*fu*dipdy
  ! niv(t+1) = niv + ibv*dipdy - ibv*dt*fv*dipdx

  ! niu(t+1) = niu + ibu*dipdx - d(BI*F*ip)dy + d(BI*F)dy*ip
  ! niv(t+1) = niv + ibv*dipdy + d(BI*F*ip)dx - d(BI*F)dx*ip
  ! BI = -dt/(imass*(1.+f**2))
  ! F = dt*f
  
  imu = 0.5*(imass(1:ifull)+imass(ie))
  imv = 0.5*(imass(1:ifull)+imass(in))
  ! note missing 0.5 as ip is the average of t and t+1
  ibu(1:ifull) = -dt*eeu(1:ifull,1)/(imu*(1.+(dt*fu(1:ifull))**2))
  ibv(1:ifull) = -dt*eev(1:ifull,1)/(imv*(1.+(dt*fv(1:ifull))**2))

  
  ! update boundary
  data_c(1:ifull,1) = sou(1:ifull)
  data_d(1:ifull,1) = sov(1:ifull)
  data_c(1:ifull,2) = snu(1:ifull)
  data_d(1:ifull,2) = snv(1:ifull)
  data_c(1:ifull,3) = squ(1:ifull)
  data_d(1:ifull,3) = sqv(1:ifull)
  data_c(1:ifull,4) = bb3u(1:ifull)
  data_d(1:ifull,4) = bb4v(1:ifull)
  data_c(1:ifull,5) = ibu(1:ifull)
  data_d(1:ifull,5) = ibv(1:ifull)
  data_c(1:ifull,6) = ibb3u(1:ifull)
  data_d(1:ifull,6) = ibb4v(1:ifull)
  data_c(1:ifull,7) = niu(1:ifull)
  data_d(1:ifull,7) = niv(1:ifull)
  call boundsuv(data_c(:,1:7),data_d(:,1:7),stag=-9) ! stag=-9 updates iwu and isv
  sou(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,1)
  sov(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,1)
  snu(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,2)
  snv(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,2)
  squ(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,3)
  sqv(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,3)
  bb3u(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,4)
  bb4v(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,4)
  ibu(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,5)
  ibv(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,5)
  ibb3u(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,6)
  ibb4v(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,6)
  niu(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,7)
  niv(ifull+1:ifull+iextra) = data_d(ifull+1:ifull+iextra,7)

  ! Iteratively solve for free surface height, eta
  ! Iterative loop to estimate ice 'pressure'
  if ( precon<-9999 ) then
    ! Multi-grid
    call mlomg(neta,sou,sov,snu,snv,squ,sqv,bb,bb3u,bb4v,xps, &
               niu,niv,ibu,ibv,ibb,ibb3u,ibb4v,               &
               ipmax,totits,maxglobseta,maxglobip)
  else
    ! Usual SOR
    write(6,*) "ERROR: MLO dynamics requires precon=-10000"
    call ccmpi_abort(-1)
  end if
  
  
  ! update bounds for solution to neta and ipice
  data_c(1:ifull,1) = neta(1:ifull)
  data_c(1:ifull,2) = ipice(1:ifull)
  call bounds(data_c(:,1:2),corner=.true.)
  neta(ifull+1:ifull+iextra)  = data_c(ifull+1:ifull+iextra,1)
  ipice(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,2)
  
  
  ! Update currents once neta is calculated
  detadxu = (neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
  detadyv = (neta(in)-neta(1:ifull))*emv(1:ifull)/ds
  oeu(1:ifull) = 0.5*(neta(ie)+neta(1:ifull))
  oev(1:ifull) = 0.5*(neta(in)+neta(1:ifull))
  ff(:) = (1.+ocneps)*0.5*dt*f(:)
  do ii = 1,wlev
    ! update currents (staggered)
    do iq = 1,ifull
      tnu(iq) = 0.5*( cc(in(iq),ii)*ff(in(iq))*neta(in(iq)) + cc(ine(iq),ii)*ff(ine(iq))*neta(ine(iq)) )
      tsu(iq) = 0.5*( cc(is(iq),ii)*ff(is(iq))*neta(is(iq)) + cc(ise(iq),ii)*ff(ise(iq))*neta(ise(iq)) )
      tev(iq) = 0.5*( cc(ie(iq),ii)*ff(ie(iq))*neta(ie(iq)) + cc(ien(iq),ii)*ff(ien(iq))*neta(ien(iq)) )
      twv(iq) = 0.5*( cc(iw(iq),ii)*ff(iw(iq))*neta(iw(iq)) + cc(iwn(iq),ii)*ff(iwn(iq))*neta(iwn(iq)) )
    end do  
    nu(1:ifull,ii) = kku(:,ii) + oou(:,ii)*oeu(1:ifull) + mmu(:,ii)*detadxu &
                   - 0.5*stwgtu(1:ifull,ii)*(tnu-tsu)*emu(1:ifull)/ds       &
                   + cc3u(1:ifull,ii)*oeu(1:ifull)
    nv(1:ifull,ii) = kkv(:,ii) + oov(:,ii)*oev(1:ifull) + mmv(:,ii)*detadyv &
                   + 0.5*stwgtv(1:ifull,ii)*(tev-twv)*emv(1:ifull)/ds       &
                   - cc4v(1:ifull,ii)*oev(1:ifull)
  end do
  

  ! error checking
  do ii = 1,wlev
    do iq = 1,ifull
      if ( abs(nu(iq,ii))>40. ) then
        tnu(iq) = 0.5*( cc(in(iq),ii)*ff(in(iq))*neta(in(iq)) + cc(ine(iq),ii)*ff(ine(iq))*neta(ine(iq)) )
        tsu(iq) = 0.5*( cc(is(iq),ii)*ff(is(iq))*neta(is(iq)) + cc(ise(iq),ii)*ff(ise(iq))*neta(ise(iq)) )
        write(6,*) "WARN: nu,iq,ii ",nu(iq,ii),iq,ii
        write(6,*) " kku,oou*oeu,mmu*detadxu,nnu*detadyu ",        &
          kku(iq,ii),oou(iq,ii)*oeu(iq),mmu(iq,ii)*detadxu(iq), &
          - 0.5*stwgtu(iq,ii)*(tnu(iq)-tsu(iq))*emu(iq)/ds      &
          + cc3u(iq,ii)*oeu(iq)
      end if
      if ( abs(nv(iq,ii))>40. ) then
        tev(iq) = 0.5*( cc(ie(iq),ii)*ff(ie(iq))*neta(ie(iq)) + cc(ien(iq),ii)*ff(ien(iq))*neta(ien(iq)) )
        twv(iq) = 0.5*( cc(iw(iq),ii)*ff(iw(iq))*neta(iw(iq)) + cc(iwn(iq),ii)*ff(iwn(iq))*neta(iwn(iq)) )
        write(6,*) "WARN: nv,iq,ii ",nv(iq,ii),iq,ii
        write(6,*) " kkv,oov*oev,mmv*detadyv,nnv*detadxv ",        &
          kkv(iq,ii),oov(iq,ii)*oev(iq),mmv(iq,ii)*detadyv(iq), &
          + 0.5*stwgtv(iq,ii)*(tev(iq)-twv(iq))*emv(iq)/ds      &
          + cc4v(iq,ii)*oev(iq)
      end if
    end do
  end do


  worku = nu(1:ifull,:)
  workv = nv(1:ifull,:)
  call mlocheck("solver",water_u=worku,water_v=workv)


  call START_LOG(wateriadv_begin)

  ! Update ice velocity with internal pressure terms
  dipdxu = (ipice(ie)-ipice(1:ifull))*emu(1:ifull)/ds
  dipdyv = (ipice(in)-ipice(1:ifull))*emv(1:ifull)/ds
  ipu = 0.5*(ipice(ie)+ipice(1:ifull))
  ipv = 0.5*(ipice(in)+ipice(1:ifull))
  ff(:) = dt*f(:)
  do iq = 1,ifull
    tnu(iq) = 0.5*( ibb(in(iq))*ff(in(iq))*ipice(in(iq))+ibb(ine(iq))*ff(ine(iq))*ipice(ine(iq)) )
    tsu(iq) = 0.5*( ibb(is(iq))*ff(is(iq))*ipice(is(iq))+ibb(ise(iq))*ff(ise(iq))*ipice(ise(iq)) )
    tev(iq) = 0.5*( ibb(ie(iq))*ff(ie(iq))*ipice(ie(iq))+ibb(ien(iq))*ff(ien(iq))*ipice(ien(iq)) )
    twv(iq) = 0.5*( ibb(iw(iq))*ff(iw(iq))*ipice(iw(iq))+ibb(iwn(iq))*ff(iwn(iq))*ipice(iwn(iq)) )
  end do  
  niu(1:ifull) = niu(1:ifull) + ibu(1:ifull)*dipdxu              &
               - 0.5*stwgtu(1:ifull,1)*(tnu-tsu)*emu(1:ifull)/ds &
               + ibb3u(1:ifull)*ipu(1:ifull)
  niv(1:ifull) = niv(1:ifull) + ibv(1:ifull)*dipdyv              &
               + 0.5*stwgtv(1:ifull,1)*(tev-twv)*emv(1:ifull)/ds &
               - ibb4v(1:ifull)*ipv(1:ifull)
  
  
  ! Normalisation factor for conserving ice flow in and out of gridbox
  call boundsuv(niu,niv,stag=-9)  
  spnet(1:ifull) = (-min(niu(iwu)*emu(iwu),0.)+max(niu(1:ifull)*emu(1:ifull),0.)     &
                    -min(niv(isv)*emv(isv),0.)+max(niv(1:ifull)*emv(1:ifull),0.))/ds

  ! ADVECT ICE ------------------------------------------------------
  ! use simple upwind scheme

  ! Horizontal advection for ice area
  data_c(1:ifull,1) = nfracice(1:ifull)/(em(1:ifull)**2)               ! data_c(:,1) is an area
  ! Horizontal advection for ice volume
  data_c(1:ifull,2) = ndic(1:ifull)*nfracice(1:ifull)/(em(1:ifull)**2) ! data_c(:,2) is a volume
  ! Horizontal advection for snow volume
  data_c(1:ifull,3) = ndsn(1:ifull)*nfracice(1:ifull)/(em(1:ifull)**2) ! data_c(:,3) is a volume
  ! Horizontal advection for ice energy store
  data_c(1:ifull,4) = nsto(1:ifull)*nfracice(1:ifull)/(em(1:ifull)**2)
  call mloexpgamm(gamm,ndic(1:ifull),ndsn(1:ifull),0)
  ! Horizontal advection for surface temperature
  data_c(1:ifull,5) = i_it(1:ifull,1)*nfracice(1:ifull)*gamm(:,1)/(em(1:ifull)**2)
  ! Horizontal advection of snow temperature
  data_c(1:ifull,6) = i_it(1:ifull,2)*nfracice(1:ifull)*gamm(:,2)/(em(1:ifull)**2)
  ! Horizontal advection of ice temperatures
  data_c(1:ifull,7) = i_it(1:ifull,3)*nfracice(1:ifull)*gamm(:,3)/(em(1:ifull)**2)
  data_c(1:ifull,8) = i_it(1:ifull,4)*nfracice(1:ifull)*gamm(:,3)/(em(1:ifull)**2) 
  if ( mloiceadv==0 ) then
    ! Conservation (original method)
    data_c(1:ifull,9) = spnet(1:ifull)
    call bounds(data_c(:,1:9))
    spnet(ifull+1:ifull+iextra) = data_c(ifull+1:ifull+iextra,9)
  else
    call bounds(data_c(:,1:8))
    spnet=0.  
  end if  
  call upwind_iceadv(data_c(:,1:8),niu,niv,spnet)
  nfracice(1:ifull) = min( max( data_c(1:ifull,1)*em(1:ifull)**2, 0. ), maxicefrac )
  ndic(1:ifull) = data_c(1:ifull,2)*em(1:ifull)**2/max(nfracice(1:ifull),1.E-10)
  ndsn(1:ifull) = data_c(1:ifull,3)*em(1:ifull)**2/max(nfracice(1:ifull),1.E-10)
  nsto(1:ifull) = data_c(1:ifull,4)*em(1:ifull)**2/max(nfracice(1:ifull),1.E-10)
  call mloexpgamm(gamm,ndic,ndsn,0)
  nit(1:ifull,1) = data_c(1:ifull,5)*em(1:ifull)**2/max(gamm(:,1)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,2) = data_c(1:ifull,6)*em(1:ifull)**2/max(gamm(:,2)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,3) = data_c(1:ifull,7)*em(1:ifull)**2/max(gamm(:,3)*nfracice(1:ifull),1.E-10)
  nit(1:ifull,4) = data_c(1:ifull,8)*em(1:ifull)**2/max(gamm(:,3)*nfracice(1:ifull),1.E-10)

  ! populate grid points that have no sea ice
  nit(1:ifull,:) = min( nit(1:ifull,:), 399. )
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
  
  call mlocheck("seaice advection",ice_tsurf=nit(1:ifull,1), &
                ice_u=niu(1:ifull),ice_v=niv(1:ifull))
  
  
  ! unstagged currents and ice velocity
  ttau(:,1:wlev) = nu(1:ifull,:)
  ttav(:,1:wlev) = nv(1:ifull,:)
  ttau(:,wlev+1) = niu(1:ifull)
  ttav(:,wlev+1) = niv(1:ifull)
  call mlounstaguv(ttau(:,1:wlev+1),ttav(:,1:wlev+1),tau(:,1:wlev+1),tav(:,1:wlev+1))
  nu(1:ifull,:) = tau(:,1:wlev)
  nv(1:ifull,:) = tav(:,1:wlev)
  niu(1:ifull) = tau(:,wlev+1)
  niv(1:ifull) = tav(:,wlev+1)
  
  worku = nu(1:ifull,:)
  workv = nv(1:ifull,:)
  call mlocheck("unstagger",water_u=worku,water_v=workv,ice_u=niu(1:ifull),ice_v=niv(1:ifull))
  
  call END_LOG(wateriadv_end)

  
  call START_LOG(watermfix_begin)

  ! volume conservation for water
  if ( nud_sfh==0 ) then
    ! this is a mfix=1 method.  mfix=3 may not be viable for neta  
    delpos(1) = 0.
    delneg(1) = 0.
    select case(mlomfix)
      case(1,2)
        neta(1:ifull) = max(min(neta(1:ifull), 130.), -130.)
        odum = (neta(1:ifull)-w_e)*ee(1:ifull,1)
        call ccglobal_posneg(odum,delpos(1),delneg(1))
        alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
        if ( abs(alph_p)>1.e-20 ) then
          neta(1:ifull) = w_e(1:ifull) + max(0.,odum)*alph_p + min(0.,odum)/alph_p
        end if
      case default
        neta(1:ifull) = max(min(neta(1:ifull), 130.), -130.)
    end select
    if ( nodrift==1 ) then
      odum = neta(1:ifull)*ee(1:ifull,1)*em(1:ifull)**2
      lsum = cmplx( 0., 0. )
      call drpdr_local(odum,lsum)
      call ccmpi_allreduce(lsum,gsum,"sumdr",comm_world)
      delta = real(gsum)/real(emsum)
      neta(1:ifull) = (neta(1:ifull) - delta)*ee(1:ifull,1)
    end if
  end if

  ! temperature conservation (usually off when nudging SSTs)
  if ( nud_sst==0 ) then
    delpos(1) = 0.
    delneg(1) = 0.
    select case( mlomfix )
      case(1)  ! no free surface
        do ii = 1,wlev
          nt(1:ifull,ii) = min( max( nt(1:ifull,ii), 150.-wrtemp ), 375.-wrtemp )
          mfixdum(:,ii) = (nt(1:ifull,ii)-w_t(:,ii))*dd(1:ifull)*ee(1:ifull,ii)
        end do
        call ccglobal_posneg(mfixdum,delpos(1),delneg(1),godsig)
        alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
        do ii = 1,wlev
          where(wtr(1:ifull,ii) .and. abs(alph_p)>1.e-20)
            nt(1:ifull,ii) = w_t(1:ifull,ii)                                              &
                             +(max(0.,mfixdum(:,ii))*alph_p+min(0.,mfixdum(:,ii))/alph_p) &
                             /dd(1:ifull)
          elsewhere
            nt(1:ifull,ii) = w_t(:,ii)              
          end where
        end do  
      case(2)  ! include free surface
        do ii = 1,wlev
          nt(1:ifull,ii) = min( max( nt(1:ifull,ii), 150.-wrtemp ), 375.-wrtemp )
          mfixdum(:,ii) = (nt(1:ifull,ii)-w_t(:,ii))*max(dd(1:ifull)+w_e(1:ifull),minwater)*ee(1:ifull,ii)
        end do
        call ccglobal_posneg(mfixdum,delpos(1),delneg(1),godsig)
        alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
        do ii = 1,wlev
          where(wtr(1:ifull,ii) .and. abs(alph_p)>1.e-20)
            nt(1:ifull,ii) = w_t(1:ifull,ii)                                              &
                             +(max(0.,mfixdum(:,ii))*alph_p+min(0.,mfixdum(:,ii))/alph_p) &
                             /max(dd(1:ifull)+w_e(1:ifull),minwater)
          elsewhere
            nt(1:ifull,ii) = w_t(:,ii)              
          end where
        end do  
      case default  
        do ii = 1,wlev
          nt(1:ifull,ii) = min( max( nt(1:ifull,ii), 150.-wrtemp ), 375.-wrtemp )
        end do  
    end select
    workdata = nt(1:ifull,:)
    call mlocheck("conservation fix",water_temp=workdata)
  end if ! nud_sst==0

  ! salinity conservation
  if ( nud_sss==0 ) then
    delpos(1) = 0.
    delneg(1) = 0.
    select case( mlomfix )
      case(1) ! no free surface
        do ii = 1,wlev
          ns(1:ifull,ii) = max( ns(1:ifull,ii), 0. )
          mfixdum(:,ii) = (ns(1:ifull,ii)-w_s(:,ii))*dd(1:ifull)*ee(1:ifull,ii)
        end do
        call ccglobal_posneg(mfixdum,delpos(1),delneg(1),godsig)
        alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
        do ii = 1,wlev
          where ( wtr(1:ifull,ii) .and. abs(alph_p)>1.e-20 )   
            ns(1:ifull,ii) = w_s(1:ifull,ii)                                            &
                           +(max(0.,mfixdum(:,ii))*alph_p+min(0.,mfixdum(:,ii))/alph_p) &
                           /dd(1:ifull)
          elsewhere
            ns(1:ifull,ii) = w_s(1:ifull,ii)  
          end where
          ns(1:ifull,ii) = max( ns(1:ifull,ii), 0. )
        end do  
      case(2) ! include free surface
        do ii = 1,wlev
          mfixdum(:,ii) = (ns(1:ifull,ii)-w_s(:,ii))*max(dd(1:ifull)+w_e(1:ifull),minwater)*ee(1:ifull,ii)
        end do
        call ccglobal_posneg(mfixdum,delpos(1),delneg(1),godsig)
        alph_p = sqrt(-delneg(1)/max(delpos(1),1.e-20))
        do ii = 1,wlev
          where ( wtr(1:ifull,ii) .and. abs(alph_p)>1.e-20 )
            ns(1:ifull,ii) = w_s(1:ifull,ii)                                            &
                           +(max(0.,mfixdum(:,ii))*alph_p+min(0.,mfixdum(:,ii))/alph_p) &
                           /max(dd(1:ifull)+w_e(1:ifull),minwater)
          elsewhere
            ns(1:ifull,ii) = w_s(1:ifull,ii)    
          end where
        end do
      case default
        do ii = 1,wlev
          ns(1:ifull,ii) = max( ns(1:ifull,ii), 0. )  
        end do  
    end select
    workdata2 = ns(1:ifull,:)
    call mlocheck("conservation fix",water_sal=workdata2)
  end if ! nud_sss==0

  if ( myid==0 .and. (ktau<=5.or.maxglobseta>tol.or.maxglobip>itol) ) then
    write(6,*) "MLODYNAMICS ",totits,maxglobseta,maxglobip
  end if

  call END_LOG(watermfix_end)


  ! reset time-step
  dt = dtin_mlo

end do ! mspec_mlo


w_u = nu(1:ifull,:)
w_v = nv(1:ifull,:)
w_t = nt(1:ifull,:)
w_s = ns(1:ifull,:)
call mlocheck("end of mlodynamics",water_temp=w_t,water_sal=w_s,water_u=w_u,water_v=w_v, &
              ice_tsurf=nit(1:ifull,1),ice_u=niu(1:ifull),ice_v=niv(1:ifull))

call mlodiffusion_work(w_u,w_v,w_t,w_s)

call mlocheck("after diffusion",water_temp=w_t,water_sal=w_s,water_u=w_u,water_v=w_v)


! STORE WATER AND ICE DATA IN MLO ------------------------------------------
do ii = 1,wlev
  call mloimport("temp",w_t(:,ii),ii,0)
  call mloimport("sal",w_s(:,ii),ii,0)
  call mloimport("u",w_u(:,ii),ii,0)
  call mloimport("v",w_v(:,ii),ii,0)
  call mloimport("w",w_ocn(:,ii),ii,0)
end do    
call mloimport("eta",neta(1:ifull),0,0)
call mloimpice("tsurf",nit(1:ifull,1),0)
call mloimpice("temp0",nit(1:ifull,2),0)
call mloimpice("temp1",nit(1:ifull,3),0)
call mloimpice("temp2",nit(1:ifull,4),0)
call mloimpice("fracice",nfracice(1:ifull),0)
call mloimpice("thick",ndic(1:ifull),0)
call mloimpice("snowd",ndsn(1:ifull),0)
call mloimpice("store",nsto(1:ifull),0)
call mloimpice("u",niu(1:ifull),0)
call mloimpice("v",niv(1:ifull),0)
where ( wtr(1:ifull,1) )
  fracice = nfracice(1:ifull)
  sicedep = ndic(1:ifull)
  snowd = ndsn(1:ifull)*1000.
end where

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use xyzinfo_m
use newmpar_m

implicit none

integer iq, k
real, dimension(ifull+iextra,wlev), intent(inout) :: cou,cov,cow
real vec1x,vec1y,vec1z,denb
real vec2x,vec2y,vec2z,vecdot
real vec3x,vec3y,vec3z,vdot1,vdot2
real(kind=8), dimension(ifull,size(cou,2)), intent(in) :: x3d,y3d,z3d

do k = 1,wlev
  do iq = 1,ifull
    !         cross product n1xn2 into vec1
    vec1x = real(y3d(iq,k)*z(iq) - y(iq)*z3d(iq,k))
    vec1y = real(z3d(iq,k)*x(iq) - z(iq)*x3d(iq,k))
    vec1z = real(x3d(iq,k)*y(iq) - x(iq)*y3d(iq,k))
    !         N.B. rotation formula is singular for small denb,
    !         but the rotation is unnecessary in this case
    denb = vec1x**2 + vec1y**2 + vec1z**2
    if ( denb>1.e-10 ) then
      vecdot = real(x3d(iq,k)*x(iq) + y3d(iq,k)*y(iq) + z3d(iq,k)*z(iq))
      vec2x = real(x3d(iq,k)*vecdot - x(iq))
      vec2y = real(y3d(iq,k)*vecdot - y(iq))
      vec2z = real(z3d(iq,k)*vecdot - z(iq))
      vec3x = real(x3d(iq,k) - vecdot*x(iq))
      vec3y = real(y3d(iq,k) - vecdot*y(iq))
      vec3z = real(z3d(iq,k) - vecdot*z(iq))
      vdot1 = (vec1x*cou(iq,k) + vec1y*cov(iq,k) + vec1z*cow(iq,k))/denb
      vdot2 = (vec2x*cou(iq,k) + vec2y*cov(iq,k) + vec2z*cow(iq,k))/denb
      cou(iq,k) = vdot1*vec1x + vdot2*vec3x
      cov(iq,k) = vdot1*vec1y + vdot2*vec3y
      cow(iq,k) = vdot1*vec1z + vdot2*vec3z
    end if
  end do ! iq
end do   ! k

return
end subroutine mlorot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use potential temperature and salinity Jacobians to calculate
! density Jacobian

subroutine tsjacobi(nti,nsi,pice,drhobardxu,drhobardyu,drhobardxv,drhobardyv, &
                    rho_dash,rhou_dash,rhov_dash)

use indices_m
use map_m, only : emu, emv
use newmpar_m
use parm_m, only : ds

implicit none

integer iq, ii
real, dimension(ifull+iextra,wlev), intent(in) :: nti, nsi
real, dimension(ifull,wlev), intent(out) :: drhobardxu, drhobardyu, drhobardxv, drhobardyv
real, dimension(ifull,wlev), intent(out) :: rhou_dash, rhov_dash, rho_dash
real, dimension(ifull,wlev) :: drhodxu, drhodyu, drhodxv, drhodyv
real, dimension(ifull+iextra), intent(in) :: pice
real, dimension(ifull+iextra,wlev,2) :: na
real, dimension(ifull+iextra,wlev) :: alpha, beta, lrho_dash, dzdum_rho
real absu, bbsu, absv, bbsv
real dnadxu1, dnadxv1, dnadyu1, dnadyv1
real dnadxu2, dnadxv2, dnadyu2, dnadyv2

! rho(x) = rho0 + rho_dash(x), where rho0 is a constant equal to wrtrho

do ii = 1,wlev
  ! neglect neta for calculating density  
  dzdum_rho(1:ifull+iextra,ii) = godsig(1:ifull+iextra,ii)*dd(1:ifull+iextra)
end do
call mloexpdensity(lrho_dash,alpha,beta,nti,nsi,dzdum_rho,pice,0,rawrho=.true.)

na(:,:,1) = min(max(271.-wrtemp,nti),373.-wrtemp)
where ( nsi<2. )
  na(:,:,2) = 0. ! 34.72 PSU with offset
elsewhere  
  na(:,:,2) = min(max(minsal, nsi),maxsal)-34.72
end where  

if ( mlojacobi==0 ) then !off
  do ii = 1,wlev
    drhobardxu(1:ifull,ii) = 0.
    drhobardxv(1:ifull,ii) = 0.
    drhobardyu(1:ifull,ii) = 0.
    drhobardyv(1:ifull,ii) = 0.
    rho_dash(:,ii) = 0.
    rhou_dash(:,ii) = 0.
    rhov_dash(:,ii) = 0.
  end do
else
  ! Here we calculate the slow contribution of the pressure gradient

  ! dP/dx = g rhobar dneta/dx + g sig D drhobar/dx + g sig neta drhobar/dx
  !                   (fast)               (slow)          (mixed)

  ! rhobar = int_0^sigma rho dsigma / sigma

  !$omp parallel do schedule(static) private(ii,iq,absu,bbsu,absv,bbsv)            &
  !$omp   private(dnadxu1,dnadyu1,dnadyv1,dnadxv1,dnadxu2,dnadyu2,dnadyv2,dnadxv2)
  do ii = 1,wlev
    do iq = 1,ifull
      absu = 0.5*(alpha(iq,ii)+alpha(ie(iq),ii))
      bbsu = 0.5*(beta(iq,ii) +beta(ie(iq),ii) )
      absv = 0.5*(alpha(iq,ii)+alpha(in(iq),ii))
      bbsv = 0.5*(beta(iq,ii) +beta(in(iq),ii) )
      ! process staggered u locations  
      dnadxu1=(na(ie(iq),ii,1)-na(iq,ii,1))*emu(iq)/ds
      dnadyu1=0.25*emu(iq)/ds*(na(in(iq),ii,1)-na(is(iq),ii,1) &
                              +na(ine(iq),ii,1)-na(ise(iq),ii,1))
      ! process staggered v locations
      dnadyv1=(na(in(iq),ii,1)-na(iq,ii,1))*emv(iq)/ds
      dnadxv1=0.25*emv(iq)/ds*(na(ie(iq),ii,1)-na(iw(iq),ii,1) &
                              +na(ien(iq),ii,1)-na(iwn(iq),ii,1))
      ! process staggered u locations  
      dnadxu2=(na(ie(iq),ii,2)-na(iq,ii,2))*emu(iq)/ds
      dnadyu2=0.25*emu(iq)/ds*(na(in(iq),ii,2)-na(is(iq),ii,2) &
                              +na(ine(iq),ii,2)-na(ise(iq),ii,2))
      ! process staggered v locations
      dnadyv2=(na(in(iq),ii,2)-na(iq,ii,2))*emv(iq)/ds
      dnadxv2=0.25*emv(iq)/ds*(na(ie(iq),ii,2)-na(iw(iq),ii,2) &
                              +na(ien(iq),ii,2)-na(iwn(iq),ii,2))
      ! This relationship neglects compression effects due to neta from the EOS.
      drhodxu(iq,ii) = -absu*dnadxu1 + bbsu*dnadxu2
      drhodxv(iq,ii) = -absv*dnadxv1 + bbsv*dnadxv2
      drhodyu(iq,ii) = -absu*dnadyu1 + bbsu*dnadyu2
      drhodyv(iq,ii) = -absv*dnadyv1 + bbsv*dnadyv2
    end do
  end do
  !$omp end parallel do

  ! integrate density gradient  
  drhobardxu(:,1) = drhodxu(:,1)*godsigu(1:ifull,1)
  drhobardxv(:,1) = drhodxv(:,1)*godsigv(1:ifull,1)
  drhobardyu(:,1) = drhodyu(:,1)*godsigu(1:ifull,1)
  drhobardyv(:,1) = drhodyv(:,1)*godsigv(1:ifull,1)
  rho_dash(:,1) = lrho_dash(1:ifull,1)
  rhou_dash(:,1) = 0.5*(lrho_dash(1:ifull,1)+lrho_dash(ie,1))
  rhov_dash(:,1) = 0.5*(lrho_dash(1:ifull,1)+lrho_dash(in,1))
  do ii = 2,wlev
    drhobardxu(:,ii) = drhobardxu(:,ii-1) + drhodxu(:,ii)*godsigu(1:ifull,ii)
    drhobardxv(:,ii) = drhobardxv(:,ii-1) + drhodxv(:,ii)*godsigv(1:ifull,ii)
    drhobardyu(:,ii) = drhobardyu(:,ii-1) + drhodyu(:,ii)*godsigu(1:ifull,ii)
    drhobardyv(:,ii) = drhobardyv(:,ii-1) + drhodyv(:,ii)*godsigv(1:ifull,ii)
    rho_dash(:,ii) = lrho_dash(1:ifull,ii)
    rhou_dash(:,ii) = 0.5*(lrho_dash(1:ifull,ii)+lrho_dash(ie,ii))
    rhov_dash(:,ii) = 0.5*(lrho_dash(1:ifull,ii)+lrho_dash(in,ii))
  end do
  do ii = 1,wlev
    drhobardxu(:,ii) = drhobardxu(:,ii)/gosighu(:,ii)
    drhobardxv(:,ii) = drhobardxv(:,ii)/gosighv(:,ii)
    drhobardyu(:,ii) = drhobardyu(:,ii)/gosighu(:,ii)
    drhobardyv(:,ii) = drhobardyv(:,ii)/gosighv(:,ii)
  end do
end if
  
return
end subroutine tsjacobi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate tidal potential

pure subroutine mlotide(eta,slon,slat,mtimer,jstart)

use parm_m, only : leap
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

stime = real(mtimer)/1440. + real(jstart) ! solar time (days)
if ( leap==0 ) then
  ctime = stime/36500. ! century (relative to 12Z 31 Dec 1899)
else if ( leap==1 ) then
  ctime = stime/36525. ! century (relative to 12Z 31 Dec 1899)
else if ( leap==2 ) then
  ctime = stime/36000.  
else
  ctime = 0. ! error
end if

mn = 270.43659+481276.89057*ctime+0.00198*ctime*ctime+0.000002*ctime*ctime*ctime ! moon
sn = 279.69668+36000.76892*ctime+0.00030*ctime*ctime                             ! sun
pn = 334.32956+4069.03404*ctime-0.01032*ctime*ctime-0.00001*ctime*ctime*ctime    ! lunar perigee
mn = mn*pi/180.
sn = sn*pi/180.
pn = pn*pi/180.
stime = (stime-real(int(stime)))*2.*pi ! solar time (rad)
ltime = stime+sn-mn                    ! lunar time (rad)
! 360*stime=360*ltime-9.3+12.2*day

coslat = cos(slat)**2
sinlat = sin(2.*slat)

eta = ba(1)*aa(1)*sinlat*cos(ltime+mn-slon)          &    ! K1 (note sign change)
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

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m

implicit none


integer iq, n, ntr
real, dimension(:,:), intent(inout) :: dumc
real, dimension(:), intent(in) :: niu,niv,spnet
real, dimension(ifull,size(dumc,2)) :: dum_out
real odum, tdum
real(kind=8) dumd

ntr = size(dumc,2)

select case( mloiceadv )
  case(0) ! Original 

    !$omp parallel do schedule(static) private(n,iq,dumd,odum,tdum)
    do n = 1,ntr
      do iq = 1,ifull
        dumd=real(dumc(iq,n),8)
        odum=0.5*dt*(niu(iwu(iq))*(dumc(iq,n)+dumc(iw(iq),n))    -abs(niu(iwu(iq)))*(dumc(iq,n)-dumc(iw(iq),n))    )*emu(iwu(iq))/ds
        if ( spnet(iw(iq))>1.e-10 ) then
          tdum=dumc(iw(iq),n)*max(niu(iwu(iq))*emu(iwu(iq)),0.)/(ds*spnet(iw(iq)))    
          odum=min(odum,tdum)
        else
          odum=min(odum,0.)
        end if
        if ( spnet(iq)>1.e-10 ) then
          tdum=dumc(iq,n)*min(niu(iwu(iq))*emu(iwu(iq)),0.)/(ds*spnet(iq))
          odum=max(odum,tdum)
        else
          odum=max(odum,0.)
        end if
        dumd=dumd+real(odum,8)

        odum=-0.5*dt*(niu(iq)*(dumc(iq,n)+dumc(ie(iq),n))+abs(niu(iq))*(dumc(iq,n)-dumc(ie(iq),n)))*emu(iq)/ds
        if ( spnet(ie(iq))>1.e-10 ) then
          tdum=-dumc(ie(iq),n)*min(niu(iq)*emu(iq),0.)/(ds*spnet(ie(iq)))
          odum=min(odum,tdum)
        else
          odum=min(odum,0.)
        end if
        if ( spnet(iq)>1.e-10 ) then
          tdum=-dumc(iq,n)*max(niu(iq)*emu(iq),0.)/(ds*spnet(iq))
          odum=max(odum,tdum)
        else
          odum=max(odum,0.)
        end if
        dumd=dumd+real(odum,8)

        odum=0.5*dt*(niv(isv(iq))*(dumc(iq,n)+dumc(is(iq),n))    -abs(niv(isv(iq)))*(dumc(iq,n)-dumc(is(iq),n))    )*emv(isv(iq))/ds
        if ( spnet(is(iq))>1.e-10 ) then
          tdum=dumc(is(iq),n)*max(niv(isv(iq))*emv(isv(iq)),0.)/(ds*spnet(is(iq)))
          odum=min(odum,tdum)
        else
          odum=min(odum,0.)
        end if
        if ( spnet(iq)>1.e-10 ) then
          tdum=dumc(iq,n)*min(niv(isv(iq))*emv(isv(iq)),0.)/(ds*spnet(iq))
          odum=max(odum,tdum)
        else
          odum=max(odum,0.)
        end if
        dumd=dumd+real(odum,8)

        odum=-0.5*dt*(niv(iq)*(dumc(iq,n)+dumc(in(iq),n))+abs(niv(iq))*(dumc(iq,n)-dumc(in(iq),n)))*emv(iq)/ds
        if ( spnet(in(iq))>1.e-10 ) then
          tdum=-dumc(in(iq),n)*min(niv(iq)*emv(iq),0.)/(ds*spnet(in(iq)))    
          odum=min(odum,tdum)
        else
          odum=min(odum,0.)
        end if
        if ( spnet(iq)>1.e-10 ) then
          tdum=-dumc(iq,n)*max(niv(iq)*emv(iq),0.)/(ds*spnet(iq))
          odum=max(odum,tdum)
        else
          odum=max(odum,0.)
        end if
        dumd=dumd+real(odum,8)

        dum_out(iq,n)=real(max(dumd,0._8))
      end do
    end do
    !$omp end parallel do
    dumc(1:ifull,1:ntr) = dum_out(1:ifull,1:ntr)
  
  case(1) ! revised
      
    !$omp parallel do schedule(static) private(n,iq,dumd,odum)
    do n = 1,ntr
      do iq = 1,ifull
        dumd = real(dumc(iq,n),8)
        odum = 0.5*dt*(niu(iwu(iq))*(dumc(iq,n)+dumc(iw(iq),n))-abs(niu(iwu(iq)))*(dumc(iq,n)-dumc(iw(iq),n)))*emu(iwu(iq))/ds
        dumd = dumd + real(odum,8)
        odum = -0.5*dt*(niu(iq)*(dumc(iq,n)+dumc(ie(iq),n))+abs(niu(iq))*(dumc(iq,n)-dumc(ie(iq),n)))*emu(iq)/ds
        dumd = dumd + real(odum,8)
        odum = 0.5*dt*(niv(isv(iq))*(dumc(iq,n)+dumc(is(iq),n))-abs(niv(isv(iq)))*(dumc(iq,n)-dumc(is(iq),n)))*emv(isv(iq))/ds
        dumd = dumd + real(odum,8)
        odum = -0.5*dt*(niv(iq)*(dumc(iq,n)+dumc(in(iq),n))+abs(niv(iq))*(dumc(iq,n)-dumc(in(iq),n)))*emv(iq)/ds
        dumd = dumd + real(odum,8)
        dum_out(iq,n) = real(max(dumd,0._8))
      end do
    end do
    !$omp end parallel do
    dumc(1:ifull,1:ntr) = dum_out(1:ifull,1:ntr)
      
  case default
    write(6,*) "ERROR: Unknown option mloiceadv = ",mloiceadv
    call ccmpi_abort(-1)

end select
      
return
end subroutine upwind_iceadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use multi-grid to solve for free surface and ice pressure

subroutine mlomg(neta,sou,sov,snu,snv,squ,sqv,bb,bb3u,bb4v,xps, &
                 niu,niv,ibu,ibv,ibb,ibb3u,ibb4v,               &
                 ipmax,totits,maxglobseta,maxglobip)

use cc_mpi
use helmsolve
use indices_m
use map_m
use newmpar_m
use parm_m

implicit none

integer iq
integer, intent(out) :: totits
real, intent(out) :: maxglobseta, maxglobip
real, dimension(ifull), intent(in) :: xps
real, dimension(ifull+iextra), intent(inout) :: neta
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull+iextra), intent(in) :: sou, sov, snu, snv, squ, sqv, bb, bb3u, bb4v
real, dimension(ifull+iextra), intent(in) :: niu, niv, ibu, ibv, ibb, ibb3u, ibb4v
real, dimension(ifull,2) :: zz, zzn, zzs, zze, zzw, rhs
real :: odiv_d, odiv_n
real, dimension(ifull) :: hh
real, dimension(ifull) :: yy, yyn, yys, yye, yyw

call mgsor_init

! ocean

!sum nu dz = sou + snu*etau + squ*detadxu + d(bb*ff*eta)dyu - bb3u*eta
!sum nv dz = sov + snv*etav + sqv*detadyv - d(bb*ff*eta)dxv + bb4v*eta

! int_0_1{dw/dsig,dsig} = w(sig=0) - w(sig=1) = 0.

! Lagrangian version of the contunity equation with D included in flux form
!   [neta + (1+eps)*0.5*dt*(neta*(du/dx+dv/dy)+d(D*u)/dx+d(D*v)/dy+dw/dsig)]^(t+1) = xps

! neta + (1+eps)*0.5*dt( neta*(d(sou)/dx+d(sov)/dy) + neta*(d(snu*eta)/dx+d(snv*eta)/dy)
!                      + neta*(d(squ*d(eta)/dx)/dx+d(sqv*d(eta)/dy)/dy)
!                      + d(D*sou)/dx+d(D*sov)/dy + d(D*snu*eta)/dx+d(D*snv*eta)/dy
!                      + d(D*squ*d(eta)/dx)/dx+d(D*sqv*d(eta)/dy)/dy)
!                      + neta*(-d(bb3u*eta)/dx+d(bb4v*eta)/dy)
!                      - d(D*bb3u*eta)/dx + d(D*bb4v*eta)/dy )
!                      = xps
!                        
! neta + (1+eps)*0.5*dt*( neta*odiv_n + odiv_d
!                         + yy*neta*(d2neta/dx2+d2neta/dy2+dneta/dx+dneta/dy)
!                         + zz*(d2neta/dx2+d2neta/dy2+dneta/dx+dneta/dy) ) = xps

! yy*neta*(d2neta/dx2+d2neta/dy2+dneta/dx+dneta/dy) + zz*(d2neta/dx2+d2neta/dy2+dneta/dx+dneta/dy) + hh*neta = rhs

yyn(:) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(sqv(1:ifull)*emv(1:ifull)/ds + 0.5*snv(1:ifull) - 0.5*bb4v(1:ifull))/emv(1:ifull)
yys(:) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(sqv(isv)*emv(isv)/ds - 0.5*snv(isv) + 0.5*bb4v(isv))/emv(isv)
yye(:) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(squ(1:ifull)*emu(1:ifull)/ds + 0.5*snu(1:ifull) + 0.5*bb3u(1:ifull))/emu(1:ifull)
yyw(:) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(squ(iwu)*emu(iwu)/ds - 0.5*snu(iwu) - 0.5*bb3u(iwu))/emu(iwu)
yy(:) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*((-sqv(1:ifull)-sqv(isv)-squ(1:ifull)-squ(iwu))/ds                           &
          +0.5*snv(1:ifull)/emv(1:ifull)-0.5*snv(isv)/emv(isv)+0.5*snu(1:ifull)/emu(1:ifull)-0.5*snu(iwu)/emu(iwu)        &
          -0.5*bb4v(1:ifull)/emv(1:ifull)+0.5*bb4v(isv)/emv(isv)+0.5*bb3u(1:ifull)/emu(1:ifull)-0.5*bb3u(iwu)/emu(iwu))

zzn(:,1) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(ddv(1:ifull)*sqv(1:ifull)*emv(1:ifull)/ds + 0.5*ddv(1:ifull)*snv(1:ifull) &
                                                     - 0.5*ddv(1:ifull)*bb4v(1:ifull))/emv(1:ifull)
zzs(:,1) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(ddv(isv)*sqv(isv)*emv(isv)/ds - 0.5*ddv(isv)*snv(isv)                     &
                                                     + 0.5*ddv(isv)*bb4v(isv))/emv(isv)
zze(:,1) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(ddu(1:ifull)*squ(1:ifull)*emu(1:ifull)/ds + 0.5*ddu(1:ifull)*snu(1:ifull) &
                                                     + 0.5*ddu(1:ifull)*bb3u(1:ifull))/emu(1:ifull)
zzw(:,1) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(ddu(iwu)*squ(iwu)*emu(iwu)/ds - 0.5*ddu(iwu)*snu(iwu)                     &
                                                     - 0.5*ddu(iwu)*bb3u(iwu))/emu(iwu)
zz(:,1) = (1.+ocneps)*0.5*dt*em(1:ifull)**2/ds*(-ddv(1:ifull)*sqv(1:ifull)/ds-ddv(isv)*sqv(isv)/ds &
             -ddu(1:ifull)*squ(1:ifull)/ds-ddu(iwu)*squ(iwu)/ds                                    &
             +0.5*ddv(1:ifull)*snv(1:ifull)/emv(1:ifull)-0.5*ddv(isv)*snv(isv)/emv(isv)            &
             +0.5*ddu(1:ifull)*snu(1:ifull)/emu(1:ifull)-0.5*ddu(iwu)*snu(iwu)/emu(iwu)            &
             -0.5*ddv(1:ifull)*bb4v(1:ifull)/emv(1:ifull)+0.5*ddv(isv)*bb4v(isv)/emv(isv)          &
             +0.5*ddu(1:ifull)*bb3u(1:ifull)/emu(1:ifull)-0.5*ddu(iwu)*bb3u(iwu)/emu(iwu))


yyn(:) = yyn(:)*ee(1:ifull,1)
yys(:) = yys(:)*ee(1:ifull,1)
yye(:) = yye(:)*ee(1:ifull,1)
yyw(:) = yyw(:)*ee(1:ifull,1)
yy(:) = yy(:)*ee(1:ifull,1)

zzn(:,1) = zzn(:,1)*ee(1:ifull,1)
zzs(:,1) = zzs(:,1)*ee(1:ifull,1)
zze(:,1) = zze(:,1)*ee(1:ifull,1)
zzw(:,1) = zzw(:,1)*ee(1:ifull,1)
zz(:,1) = zz(:,1)*ee(1:ifull,1)

do iq = 1,ifull

  ! prep ocean gradient terms - constant
  odiv_d = (sou(iq)*ddu(iq)/emu(iq)-sou(iwu(iq))*ddu(iwu(iq))/emu(iwu(iq))  &
           +sov(iq)*ddv(iq)/emv(iq)-sov(isv(iq))*ddv(isv(iq))/emv(isv(iq))) &
           *em(iq)**2/ds

  odiv_n = (sou(iq)/emu(iq)-sou(iwu(iq))/emu(iwu(iq))+sov(iq)/emv(iq)-sov(isv(iq))/emv(isv(iq))) &
           *em(iq)**2/ds

  odiv_d = odiv_d*ee(iq,1)
  odiv_n = odiv_n*ee(iq,1)

  hh(iq) = 1. + (1.+ocneps)*0.5*dt*odiv_n

  rhs(iq,1) = xps(iq) - (1.+ocneps)*0.5*dt*odiv_d
  rhs(iq,1) = rhs(iq,1)*ee(iq,1)

end do

! ice
! zz*(d2ipice/dx2 + d2ipice/dy2) + zz*d2ipice/dxdy = rhs

zzn(:,2) = -ibv(1:ifull)/ds + 0.5*ibb4v(1:ifull)/emv(1:ifull)
zzs(:,2) = -ibv(isv)/ds - 0.5*ibb4v(isv)/emv(isv)
zze(:,2) = -ibu(1:ifull)/ds - 0.5*ibb3u(1:ifull)/emu(1:ifull)
zzw(:,2) = -ibu(iwu)/ds + 0.5*ibb3u(iwu)/emu(iwu)
zz(:,2)  = (ibv(1:ifull)+ibv(isv)+ibu(1:ifull)+ibu(iwu))/ds          &
           +0.5*ibb4v(1:ifull)/emv(1:ifull)-0.5*ibb4v(isv)/emv(isv)  &
           -0.5*ibb3u(1:ifull)/emu(1:ifull)+0.5*ibb3u(iwu)/emu(iwu)


zzn(:,2) = zzn(:,2)*ee(1:ifull,1)
zzs(:,2) = zzs(:,2)*ee(1:ifull,1)
zze(:,2) = zze(:,2)*ee(1:ifull,1)
zzw(:,2) = zzw(:,2)*ee(1:ifull,1)
zz(:,2) = zz(:,2)*ee(1:ifull,1)

rhs(:,2) = min(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu) &
              +niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv),0.)
rhs(:,2) = rhs(:,2)*ee(1:ifull,1)

! Ensure that zero is a valid solution for ice free grid points
where ( zz(:,2)>=0. )
  zz(:,2)  = -dt/(ds*100.) ! 100 is the minimum imass
  zzn(:,2) = 0.
  zzs(:,2) = 0.
  zze(:,2) = 0.
  zzw(:,2) = 0.
  rhs(:,2) = 0.
end where

call mgmlo(neta,ipice,yy,yyn,yys,yye,yyw,zz,zzn,zzs,zze,zzw,   &
           hh,rhs,tol,itol,totits,maxglobseta,maxglobip,ipmax, &
           ee(:,1),dd)

return
end subroutine mlomg

end module mlodynamics
