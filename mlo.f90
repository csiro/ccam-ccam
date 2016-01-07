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
    
! This is a 1D, mixed layer ocean model based on Large, et al (1994), for ensemble regional climate simulations.
! This code is also used for modelling lakes in ACCESS and CCAM.  In CCAM this module interfaces with
! mlodynamics.f90 for river routing, diffusion and advection routines.

! This version has a relatively thin 1st layer (e.g, 0.5m) so as to better reproduce a diurnal cycle in SST.  It also
! supports a thermodynamic model of sea ice based on O'Farrell's from Mk3.5.  We have included a free surface so that
! lakes can change depth, etc.

! This version can assimilate SSTs from GCMs, using a convolution based digital filter (see nestin.f),
! which avoids problems with complex land-sea boundary conditions.

! Order of calls:
!    call mloinit
!    call mloload
!    ...
!    main loop
!       ...
!       call mloalb4
!       ...
!       call mloeval
!       ...
!       call mloscrnout
!       ...
!    end loop
!    ...
!    call mlosave
!    call mloend

module mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,mlodiag,mloalb2,mloalb4, &
       mloscrnout,mloextra,mloimpice,mloexpice,mloexpdep,mloexpdensity,mloexpmelt,mloexpgamm,        &
       mloimport3d,mloexport3d
public micdwn
public wlev,zomode,wrtemp,onedice,mxd,mindep,minwater,zoseaice,factchseaice

! parameters
integer, save      :: wlev = 20                                        ! Number of water layers
integer, parameter :: iqx = 4148                                       ! Diagnostic grid point (host index)
! model arrays
integer, save :: wfull, ifull, iqwt                                    ! Grid size and dignostic point (local index)
logical, dimension(:), allocatable, save :: wpack                      ! Map for packing/unpacking water points
real, dimension(:,:), allocatable, save :: depth, dz, depth_hl, dz_hl  ! Column depth and thickness (m)
real, dimension(:,:), allocatable, save :: micdwn                      ! This variable is for CCAM onthefly.f

type waterdata
  real, dimension(:,:), allocatable :: temp         ! water layer temperature delta
  real, dimension(:,:), allocatable :: sal          ! water layer salinity (PSU)
  real, dimension(:,:), allocatable :: u            ! water u-current (m/s)
  real, dimension(:,:), allocatable :: v            ! water v-current (m/s)
  real, dimension(:), allocatable :: eta            ! water free surface height (m)
end type waterdata

type icedata
  real, dimension(:), allocatable :: thick          ! ice thickness (m)
  real, dimension(:), allocatable :: snowd          ! snow thickness (m)
  real, dimension(:), allocatable :: fracice        ! ice area cover fraction
  real, dimension(:), allocatable :: tsurf          ! surface temperature (K)
  real, dimension(:,:), allocatable :: temp         ! layer temperature (K)
  real, dimension(:), allocatable :: store          ! brine energy storage (J/m^2)
  real, dimension(:), allocatable :: u              ! ice u-velocity (m/s)
  real, dimension(:), allocatable :: v              ! ice v-velocity (m/s)
  real, dimension(:), allocatable :: sal            ! ice salinity (PSU)
end type icedata

type dgwaterdata
  real, dimension(:), allocatable :: mixdepth       ! mixed layer depth (m)
  integer, dimension(:), allocatable :: mixind      ! index for mixed layer depth
  real, dimension(:), allocatable :: bf             ! bouyancy forcing
  real, dimension(:), allocatable :: visdiralb      ! water VIS direct albedo
  real, dimension(:), allocatable :: visdifalb      ! water VIS difuse albedo
  real, dimension(:), allocatable :: nirdiralb      ! water NIR direct albedo
  real, dimension(:), allocatable :: nirdifalb      ! water NIR direct albedo
  real, dimension(:), allocatable :: zo             ! water roughness length for momentum (m)
  real, dimension(:), allocatable :: zoh            ! water roughness length for heat (m)
  real, dimension(:), allocatable :: zoq            ! water roughness length for moisture (m)
  real, dimension(:), allocatable :: cd             ! water drag coeff for momentum
  real, dimension(:), allocatable :: cdh            ! water drag coeff for heat
  real, dimension(:), allocatable :: cdq            ! water drag coeff for moisture
  real, dimension(:), allocatable :: fg             ! water sensible heat flux (W/m2)
  real, dimension(:), allocatable :: eg             ! water latent heat flux (W/m2)
  real, dimension(:), allocatable :: taux           ! wind/water u-component stress (N/m2)
  real, dimension(:), allocatable :: tauy           ! wind/water v-component stress (N/m2)
  !real(kind=8), dimension(:), allocatable :: deleng         ! Change in energy stored
end type dgwaterdata

type dgicedata
  real, dimension(:), allocatable :: wetfrac        ! ice area cover of liquid water
  real, dimension(:), allocatable :: visdiralb      ! ice VIS direct albedo
  real, dimension(:), allocatable :: visdifalb      ! ice VIS difuse albedo
  real, dimension(:), allocatable :: nirdiralb      ! ice NIR direct albedo
  real, dimension(:), allocatable :: nirdifalb      ! ice NIR direct albedo
  real, dimension(:), allocatable :: zo             ! ice roughness length for momentum (m)
  real, dimension(:), allocatable :: zoh            ! ice roughness length for heat (m)
  real, dimension(:), allocatable :: zoq            ! ice roughness length for moisture (m)
  real, dimension(:), allocatable :: cd             ! ice drag coeff for momentum
  real, dimension(:), allocatable :: cdh            ! ice drag coeff for heat
  real, dimension(:), allocatable :: cdq            ! ice drag coeff for moisture
  real, dimension(:), allocatable :: fg             ! ice sensible heat flux (W/m2)
  real, dimension(:), allocatable :: eg             ! ice latent heat flux (W/m2)
  real, dimension(:), allocatable :: tauxica        ! water/ice u-component stress (N/m2)
  real, dimension(:), allocatable :: tauyica        ! water/ice v-component stress (N/m2)
  real, dimension(:), allocatable :: tauxicw        ! water/ice u-component stress (N/m2)
  real, dimension(:), allocatable :: tauyicw        ! water/ice v-component stress (N/m2)
  !real(kind=8), dimension(:), allocatable :: deleng         ! Change in energy stored
end type dgicedata

type dgscrndata
  real, dimension(:), allocatable :: temp           ! 2m air temperature (K)    
  real, dimension(:), allocatable :: qg             ! 2m water vapor mixing ratio (kg/kg)
  real, dimension(:), allocatable :: u2             ! 2m wind speed (m/s)
  real, dimension(:), allocatable :: u10            ! 10m wind speed (m/s)
end type dgscrndata

type(waterdata), save :: water
type(icedata), save :: ice
type(dgwaterdata), save :: dgwater
type(dgicedata), save :: dgice
type(dgscrndata), save :: dgscrn
  
! mode
integer, parameter :: incradbf  = 1       ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 1       ! include shortwave in non-local term
integer, save      :: zomode    = 2       ! roughness calculation (0=Charnock (CSIRO9), 1=Charnock (zot=zom), 2=Beljaars)
integer, parameter :: mixmeth   = 1       ! Refine mixed layer depth calculation (0=None, 1=Iterative)
integer, parameter :: deprelax  = 0       ! surface height (0=vary, 1=relax, 2=set to zero)
integer, save      :: onedice   = 1       ! use 1D ice model (0=Off, 1=On)
! model parameters
real, save :: mxd      = 5002.18          ! Max depth (m)
real, save :: mindep   = 1.               ! Thickness of first layer (m)
real, save :: minwater = 5.               ! Minimum water depth (m)
real, parameter :: ric     = 0.3          ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1          ! Ratio of surface layer and mixed layer thickness
real, parameter :: minsfc  = 1.           ! Minimum thickness to average surface layer properties (m)
real, parameter :: maxsal  = 50.          ! Maximum salinity used in density and melting point calculations (PSU)
real, parameter :: mu_1    = 23.          ! VIS depth (m) - Type I
real, parameter :: mu_2    = 0.35         ! NIR depth (m) - Type I
real, parameter :: fluxwgt = 0.7          ! Time filter for flux calculations
real, parameter :: wrtemp  = 290.         ! Water reference temperature (K)
! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5             ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf               ! Latent heat of sublimation (J kg^-1)
real, parameter :: grav=9.80              ! graviational constant (m s^-2)
real, parameter :: sbconst=5.67e-8        ! Stefan-Boltzmann constant
real, parameter :: cdbot=2.4e-3           ! bottom drag coefficent
real, parameter :: cpair=1004.64          ! Specific heat of dry air at const P
real, parameter :: rdry=287.04            ! Specific gas const for dry air
real, parameter :: rvap=461.5             ! Gas constant for water vapor
! water parameters
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: rhowt=1025.            ! reference water density (Boussinesq fluid) (kg/m3)
real, parameter :: salwt=34.72            ! reference water salinity (PSU)
! ice parameters
real, save      :: zoseaice=0.0005        ! roughnes length for sea-ice (m)
real, save      :: factchseaice=1.        ! =sqrt(zo/zoh) for sea-ice
real, parameter :: himin=0.1              ! minimum ice thickness for multiple layers (m)
real, parameter :: icemin=0.01            ! minimum ice thickness (m)
real, parameter :: icemax=6.              ! maximum ice thickness (m)
real, parameter :: rhoic=900.             ! ice density (kg/m3)
real, parameter :: rhosn=330.             ! snow density (kg/m3)
real, parameter :: qice=lf*rhoic          ! latent heat of fusion for ice (J m^-3)
real, parameter :: qsnow=lf*rhosn         ! latent heat of fusion for snow (J m^-3)
real, parameter :: cpi=1.8837e6           ! specific heat ice  (J/m**3/K)
real, parameter :: cps=6.9069e5           ! specific heat snow (J/m**3/K)
real, parameter :: condice=2.03439        ! conductivity ice
real, parameter :: condsnw=0.30976        ! conductivity snow
real, parameter :: gammi=0.5*cpi*himin    ! specific heat*depth (for ice/snow) (J m^-2 K^-1)
real, parameter :: emisice=1.             ! emissivity of ice (0.95?)
real, parameter :: maxicesal=4.           ! Maximum salinity for sea-ice (PSU)
! 1-D ice model
!real, parameter :: fracbreak=0.05        ! Minimum ice fraction
! stability function parameters
real, parameter :: bprm=5.                ! 4.7 in rams
real, parameter :: chs=2.6                ! 5.3 in rams
real, parameter :: cms=5.                 ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialise MLO

subroutine mloinit(ifin,depin,diag)

implicit none

integer, intent(in) :: ifin,diag
integer iq,iqw,ii
real, dimension(ifin), intent(in) :: depin
real, dimension(ifin) :: deptmp
real, dimension(wlev) :: dumdf
real, dimension(wlev+1) :: dumdh

if (diag>=1) write(6,*) "Initialising MLO"

ifull=ifin
allocate(wpack(ifull))
wpack=depin>minwater
wfull=count(wpack)
if (wfull==0) then
  deallocate(wpack)
  return
end if

allocate(water%temp(wfull,wlev),water%sal(wfull,wlev))
allocate(water%u(wfull,wlev),water%v(wfull,wlev))
allocate(water%eta(wfull))
allocate(ice%temp(wfull,0:2),ice%thick(wfull),ice%snowd(wfull))
allocate(ice%fracice(wfull),ice%tsurf(wfull),ice%store(wfull))
allocate(ice%u(wfull),ice%v(wfull),ice%sal(wfull))
allocate(dgwater%visdiralb(wfull),dgwater%visdifalb(wfull))
allocate(dgwater%nirdiralb(wfull),dgwater%nirdifalb(wfull))
allocate(dgice%visdiralb(wfull),dgice%visdifalb(wfull))
allocate(dgice%nirdiralb(wfull),dgice%nirdifalb(wfull))
allocate(dgwater%zo(wfull),dgwater%zoh(wfull),dgwater%zoq(wfull))
allocate(dgwater%cd(wfull),dgwater%cdh(wfull),dgwater%cdq(wfull))
allocate(dgwater%fg(wfull),dgwater%eg(wfull))
allocate(dgice%zo(wfull),dgice%zoh(wfull),dgice%zoq(wfull))
allocate(dgice%cd(wfull),dgice%cdh(wfull),dgice%cdq(wfull))
allocate(dgice%fg(wfull),dgice%eg(wfull))
allocate(dgwater%mixdepth(wfull),dgwater%bf(wfull))
allocate(dgice%wetfrac(wfull),dgwater%mixind(wfull))
allocate(dgscrn%temp(wfull),dgscrn%u2(wfull),dgscrn%qg(wfull),dgscrn%u10(wfull))
allocate(dgwater%taux(wfull),dgwater%tauy(wfull))
allocate(dgice%tauxica(wfull),dgice%tauyica(wfull))
allocate(dgice%tauxicw(wfull),dgice%tauyicw(wfull))
allocate(depth(wfull,wlev),dz(wfull,wlev))
allocate(depth_hl(wfull,wlev+1))
allocate(dz_hl(wfull,2:wlev))
!allocate(dgwater%deleng(wfull),dgice%deleng(wfull))

iqw=0
iqwt=0
do iq=1,ifull
  if (wpack(iq)) then
    iqw=iqw+1
    if (iq>=iqx) then
      iqwt=iqw
      exit
    end if
  end if
end do

water%temp=288.-wrtemp    ! K
water%sal=35.             ! PSU
water%u=0.                ! m/s
water%v=0.                ! m/s
water%eta=0.              ! m

ice%thick=0.              ! m
ice%snowd=0.              ! m
ice%fracice=0.            ! %
ice%tsurf=271.2           ! K
ice%temp=271.2            ! K
ice%store=0.
ice%u=0.                  ! m/s
ice%v=0.                  ! m/s
ice%sal=0.                ! PSU

dgwater%mixdepth=-1.      ! m
dgwater%bf=0.
dgwater%visdiralb=0.06
dgwater%visdifalb=0.06
dgwater%nirdiralb=0.06
dgwater%nirdifalb=0.06
dgice%visdiralb=0.65
dgice%visdifalb=0.65
dgice%nirdiralb=0.65
dgice%nirdifalb=0.65
dgwater%zo=0.001
dgwater%zoh=0.001
dgwater%zoq=0.001
dgwater%cd=0.
dgwater%cdh=0.
dgwater%cdq=0.
dgwater%fg=0.
dgwater%eg=0.
dgice%zo=0.001
dgice%zoh=0.001
dgice%zoq=0.001
dgice%cd=0.
dgice%cdh=0.
dgice%cdq=0.
dgice%fg=0.
dgice%eg=0.
dgice%wetfrac=1.
dgscrn%temp=273.2
dgscrn%qg=0.
dgscrn%u2=0.
dgscrn%u10=0.
dgwater%mixind=wlev-1
dgwater%taux=0.
dgwater%tauy=0.
dgice%tauxica=0.
dgice%tauyica=0.
dgice%tauxicw=0.
dgice%tauyicw=0.
!dgwater%deleng=0._8
!dgice%deleng=0._8

! MLO - 20 level
!depth = (/   0.5,   3.4,  11.9,  29.7,  60.7, 108.5, 176.9, 269.7, 390.6, 543.3, &
!           731.6, 959.2,1230.0,1547.5,1915.6,2338.0,2818.5,3360.8,3968.7,4645.8 /)
! MLO - 30 level
!depth = (/   0.5,   2.1,   5.3,  11.3,  21.1,  36.0,  56.9,  85.0, 121.5, 167.3, &
!           223.7, 291.7, 372.4, 467.0, 576.5, 702.1, 844.8,1005.8,1186.2,1387.1, &
!          1609.6,1854.7,2123.7,2417.6,2737.5,3084.6,3456.9,3864.5,4299.6,4766.2 /)
! MLO - 40 level
!depth = (/   0.5,   1.7,   3.7,   6.8,  11.5,  18.3,  27.7,  40.1,  56.0,  75.9, &
!           100.2, 129.4, 164.0, 204.3, 251.0, 304.4, 365.1, 433.4, 509.9, 595.0, &
!           689.2, 793.0, 906.7,1031.0,1166.2,1312.8,1471.3,1642.2,1825.9,2022.8, &
!          2233.5,2458.4,2698.0,2952.8,3233.1,3509.5,3812.5,4132.4,4469.9,4825.3 /)
! MLO - 50 level
!depth = (/   0.5,   1.7,   3.1,   5.2,   8.1,  12.1,  17.3,  24.2,  32.8,  43.4, &
!            56.3,  71.7,  89.9, 111.0, 135.3, 163.1, 194.5, 229.9, 269.5, 313.4, &
!           362.0, 415.5, 474.1, 538.1, 607.6, 683.0, 764.4, 852.2, 946.5,1047.6, &
!          1155.7,1271.0,1393.9,1524.5,1663.0,1809.8,1965.0,2128.9,2301.8,2483.8, &
!          2675.2,2876.2,3087.1,3308.1,3539.5,3781.5,4034.3,4298.2,4573.4,4860.1 /)
! Mk3.5 - 31 level
!depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
!           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1,1095.2, &
!          1284.7,1502.9,1758.8,2054.3,2400.0,2800.0,3200.0,3600.0,4000.0,4400.0, &
!          4800.0 /)
!ACCESS 1.3 - 50 level
!depth = (/   5.0,  15.0,  25.0,  35.0,  45.0,  55.0,  65.0,  75.0,  85.0,  95.0, &
!           105.0, 115.0, 125.0, 135.0, 145.0, 155.0, 165.0, 175.0, 185.0, 195.0, &
!           205.0, 216.8, 241.4, 280.8, 343.3, 427.3, 536.7, 665.4, 812.8, 969.1, &
!          1130.9,1289.6,1455.8,1622.9,1801.6,1984.9,2182.9,2388.4,2610.9,2842.6, &
!          3092.2,3351.3,3628.1,3913.3,4214.5,4521.9,4842.6,5166.1,5499.2,5831.3 /)

deptmp(1:wfull)=pack(depin,wpack)

!smxd=maxval(deptmp(1:wfull))
!smnd=minval(deptmp(1:wfull))

do iqw=1,wfull
  call vgrid(wlev,deptmp(iqw),dumdf,dumdh)
  depth(iqw,:)=dumdf
  depth_hl(iqw,:)=dumdh
  !if (smxd==deptmp(iqw)) then
  !  write(6,*) "MLO max depth ",depth(iqw,:)
  !  smxd=smxd+10.
  !end if
  !if (smnd==deptmp(iqw)) then
  !  write (6,*) "MLO min depth ",depth(iqw,:)
  !  smnd=smnd-10.
  !end if
end do
do ii=1,wlev
  dz(:,ii)=depth_hl(:,ii+1)-depth_hl(:,ii)
end do
do ii=2,wlev
  dz_hl(:,ii)=depth(:,ii)-depth(:,ii-1)
end do

if ( minsfc>minwater ) then
  write(6,*) "ERROR: MLO parameters are invalid.  minsfc>minwater"
  stop
end if

return
end subroutine mloinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define vertical grid

subroutine vgrid(wlin,depin,depthout,depth_hlout)

implicit none

integer, intent(in) :: wlin
integer ii
real, intent(in) :: depin
real, dimension(wlin), intent(out) :: depthout
real, dimension(wlin+1), intent(out) :: depth_hlout
real dd,x,al,bt

dd=min(mxd,max(mindep,depin))
x=real(wlin)
!if (dd>x) then
!  al=(mindep*x-dd)/(x-x*x*x)           ! hybrid levels
!  bt=(mindep*x*x*x-dd)/(x*x*x-x)       ! hybrid levels
  al=dd*(mindep*x/mxd-1.)/(x-x*x*x)     ! sigma levels
  bt=dd*(mindep*x*x*x/mxd-1.)/(x*x*x-x) ! sigma levels 
  do ii=1,wlin+1
    x=real(ii-1)
    depth_hlout(ii)=al*x*x*x+bt*x ! ii is for half level ii-0.5
  end do
!else
!  depth_hlout(1)=0.
!  depth_hlout(2)=mindep       ! hybrid levels
!  do ii=3,wlev+1
!    x=(dd-mindep)*real(ii-2)/real(wlev-1)+mindep             ! hybrid levels
!    depth_hlout(ii)=x
!  end do
!end if
do ii=1,wlin
  depthout(ii)=0.5*(depth_hlout(ii)+depth_hlout(ii+1))
end do

return
end subroutine vgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate MLO arrays

subroutine mloend

implicit none

if (wfull==0) return

deallocate(wpack)
deallocate(water%temp,water%sal,water%u,water%v)
deallocate(water%eta)
deallocate(ice%temp,ice%thick,ice%snowd)
deallocate(ice%fracice,ice%tsurf,ice%store)
deallocate(ice%u,ice%v,ice%sal)
deallocate(dgwater%mixdepth,dgwater%bf)
deallocate(dgwater%visdiralb,dgwater%visdifalb)
deallocate(dgwater%nirdiralb,dgwater%nirdifalb)
deallocate(dgice%visdiralb,dgice%visdifalb)
deallocate(dgice%nirdiralb,dgice%nirdifalb)
deallocate(dgwater%zo,dgwater%zoh,dgwater%zoq,dgwater%cd,dgwater%cdh,dgwater%fg,dgwater%eg)
deallocate(dgice%zo,dgice%zoh,dgice%zoq,dgice%cd,dgice%cdh,dgice%fg,dgice%eg)
deallocate(dgwater%cdq,dgice%cdq)
deallocate(dgice%wetfrac,dgwater%mixind)
deallocate(dgscrn%temp,dgscrn%u2,dgscrn%qg,dgscrn%u10)
deallocate(dgwater%taux,dgwater%tauy,dgice%tauxica,dgice%tauyica)
deallocate(dgice%tauxicw,dgice%tauyicw)
deallocate(depth,dz,depth_hl,dz_hl)
!deallocate(dgwater%deleng,dgice%deleng)

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(datain,shin,icein,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(in) :: datain
real, dimension(ifull,11), intent(in) :: icein
real, dimension(ifull), intent(in) :: shin

if (wfull==0) return

do ii=1,wlev
  water%temp(:,ii)=pack(datain(:,ii,1),wpack)
  water%sal(:,ii) =pack(datain(:,ii,2),wpack)
  water%u(:,ii)   =pack(datain(:,ii,3),wpack)
  water%v(:,ii)   =pack(datain(:,ii,4),wpack)
end do
ice%tsurf(:)   =pack(icein(:,1),wpack)
ice%temp(:,0)  =pack(icein(:,2),wpack)
ice%temp(:,1)  =pack(icein(:,3),wpack)
ice%temp(:,2)  =pack(icein(:,4),wpack)
ice%fracice(:) =pack(icein(:,5),wpack)
ice%thick(:)   =pack(icein(:,6),wpack)
ice%snowd(:)   =pack(icein(:,7),wpack)
ice%store(:)   =pack(icein(:,8),wpack)
ice%u(:)       =pack(icein(:,9),wpack)
ice%v(:)       =pack(icein(:,10),wpack)
ice%sal(:)     =pack(icein(:,11),wpack)
water%eta(:)   =pack(shin(:),wpack)

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(dataout,depout,shout,iceout,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(inout) :: dataout
real, dimension(ifull,11), intent(inout) :: iceout
real, dimension(ifull), intent(inout) :: depout,shout

iceout(:,8)=0.
depout=0.
shout=0.

if (wfull==0) return

do ii=1,wlev
  dataout(:,ii,1)=unpack(water%temp(:,ii),wpack,dataout(:,ii,1))
  dataout(:,ii,2)=unpack(water%sal(:,ii),wpack,dataout(:,ii,2))
  dataout(:,ii,3)=unpack(water%u(:,ii),wpack,dataout(:,ii,3))
  dataout(:,ii,4)=unpack(water%v(:,ii),wpack,dataout(:,ii,4))
end do
iceout(:,1)   =unpack(ice%tsurf(:),wpack,iceout(:,1))
iceout(:,2)   =unpack(ice%temp(:,0),wpack,iceout(:,2))
iceout(:,3)   =unpack(ice%temp(:,1),wpack,iceout(:,3))
iceout(:,4)   =unpack(ice%temp(:,2),wpack,iceout(:,4))
iceout(:,5)   =unpack(ice%fracice(:),wpack,iceout(:,5))
iceout(:,6)   =unpack(ice%thick(:),wpack,iceout(:,6))
iceout(:,7)   =unpack(ice%snowd(:),wpack,iceout(:,7))
iceout(:,8)   =unpack(ice%store(:),wpack,iceout(:,8))
iceout(:,9)   =unpack(ice%u(:),wpack,iceout(:,9))
iceout(:,10)  =unpack(ice%v(:),wpack,iceout(:,10))
iceout(:,11)  =unpack(ice%sal(:),wpack,iceout(:,11))
depout(:)     =unpack(depth_hl(:,wlev+1),wpack,depout)
shout(:)      =unpack(water%eta(:),wpack,shout)

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mloimport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(in) :: sst

if (wfull==0) return

select case(mode)
  case(0)
    water%temp(:,ilev)=pack(sst,wpack)
  case(1)
    water%sal(:,ilev)=pack(sst,wpack)
  case(2)
    water%u(:,ilev)=pack(sst,wpack)
  case(3)
    water%v(:,ilev)=pack(sst,wpack)
  case(4)
    water%eta=pack(sst,wpack)
end select

return
end subroutine mloimport

subroutine mloimport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
integer ii
real, dimension(ifull,wlev), intent(in) :: sst

if (wfull==0) return

select case(mode)
  case(0)
    do ii=1,wlev
      water%temp(:,ii)=pack(sst(:,ii),wpack)
    end do
  case(1)
    do ii=1,wlev
      water%sal(:,ii) =pack(sst(:,ii),wpack)
    end do
  case(2)
    do ii=1,wlev
      water%u(:,ii)   =pack(sst(:,ii),wpack)
    end do
  case(3)
    do ii=1,wlev
      water%v(:,ii)   =pack(sst(:,ii),wpack)
    end do
end select

return
end subroutine mloimport3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloimpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull==0) return

select case(ilev)
  case(1)
    ice%tsurf=pack(tsn,wpack)
  case(2)
    ice%temp(:,0)=pack(tsn,wpack)
  case(3)
    ice%temp(:,1)=pack(tsn,wpack)
  case(4)
    ice%temp(:,2)=pack(tsn,wpack)
  case(5)
    ice%fracice=pack(tsn,wpack)
  case(6)
    ice%thick=pack(tsn,wpack)
  case(7)
    ice%snowd=pack(tsn,wpack)
  case(8)
    ice%store=pack(tsn,wpack)
  case(9)
    ice%u=pack(tsn,wpack)
  case(10)
    ice%v=pack(tsn,wpack)
  case(11)
    ice%sal=pack(tsn,wpack)
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",ilev
    stop
end select

return
end subroutine mloimpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mloexport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(inout) :: sst

if (wfull==0) return

select case(mode)
  case(0)
    sst=unpack(water%temp(:,ilev),wpack,sst)
  case(1)
    sst=unpack(water%sal(:,ilev),wpack,sst)
  case(2)
    sst=unpack(water%u(:,ilev),wpack,sst)
  case(3)
    sst=unpack(water%v(:,ilev),wpack,sst)
  case(4)
    sst=unpack(water%eta,wpack,sst)
end select

return
end subroutine mloexport

subroutine mloexport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
integer ii
real, dimension(ifull,wlev), intent(inout) :: sst

if (wfull==0) return

select case(mode)
  case(0)
    do ii=1,wlev
      sst(:,ii)=unpack(water%temp(:,ii),wpack,sst(:,ii))
    end do
  case(1)
    do ii=1,wlev
      sst(:,ii)=unpack(water%sal(:,ii),wpack,sst(:,ii))
    end do
  case(2)
    do ii=1,wlev
      sst(:,ii)=unpack(water%u(:,ii),wpack,sst(:,ii))
    end do
  case(3)
    do ii=1,wlev
      sst(:,ii)=unpack(water%v(:,ii),wpack,sst(:,ii))
    end do
end select

return
end subroutine mloexport3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloexpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull==0) return

select case(ilev)
  case(1)
    tsn=unpack(ice%tsurf,wpack,tsn)
  case(2)
    tsn=unpack(ice%temp(:,0),wpack,tsn)
  case(3)
    tsn=unpack(ice%temp(:,1),wpack,tsn)
  case(4)
    tsn=unpack(ice%temp(:,2),wpack,tsn)
  case(5)
    tsn=unpack(ice%fracice,wpack,tsn)
  case(6)
    tsn=unpack(ice%thick,wpack,tsn)
  case(7)
    tsn=unpack(ice%snowd,wpack,tsn)
  case(8)
    tsn=unpack(ice%store,wpack,tsn)
  case(9)
    tsn=unpack(ice%u,wpack,tsn)
  case(10)
    tsn=unpack(ice%v,wpack,tsn)
  case(11)
    tsn=unpack(ice%sal,wpack,tsn)
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",ilev
    stop
end select

return
end subroutine mloexpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return mixed layer depth

subroutine mlodiag(mld,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(out) :: mld

mld=0.
if (wfull==0) return
mld=unpack(dgwater%mixdepth,wpack,mld)

return
end subroutine mlodiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return roughness length for heat

subroutine mloextra(mode,zoh,zmin,diag)

implicit none

integer, intent(in) :: mode,diag
real, dimension(ifull), intent(out) :: zoh
real, dimension(ifull), intent(in) :: zmin
real, dimension(wfull) :: atm_zmin
real, dimension(wfull) :: workb,workc

zoh=0.
if (wfull==0) return

select case(mode)
  case(0) ! zoh
    atm_zmin=pack(zmin,wpack)
    workb=(1.-ice%fracice)/(log(atm_zmin/dgwater%zo)*log(atm_zmin/dgwater%zoh)) &
         +ice%fracice/(log(atm_zmin/dgice%zo)*log(atm_zmin/dgice%zoh))
    workc=(1.-ice%fracice)/log(atm_zmin/dgwater%zo)**2+ice%fracice/log(atm_zmin/dgice%zo)**2
    workc=sqrt(workc)
    zoh=unpack(atm_zmin*exp(-workc/workb),wpack,zoh)
  case(1) ! taux
    workb=(1.-ice%fracice)*dgwater%taux+ice%fracice*dgice%tauxica
    zoh=unpack(workb,wpack,zoh)
  case(2) ! tauy
    workb=(1.-ice%fracice)*dgwater%tauy+ice%fracice*dgice%tauyica
    zoh=unpack(workb,wpack,zoh)
  case(3) ! zoq
    atm_zmin=pack(zmin,wpack)
    workb=(1.-ice%fracice)/(log(atm_zmin/dgwater%zo)*log(atm_zmin/dgwater%zoq)) &
         +ice%fracice/(log(atm_zmin/dgice%zo)*log(atm_zmin/dgice%zoq))
    workc=(1.-ice%fracice)/log(atm_zmin/dgwater%zo)**2+ice%fracice/log(atm_zmin/dgice%zo)**2
    workc=sqrt(workc)
    zoh=unpack(atm_zmin*exp(-workc/workb),wpack,zoh)
  case default
    write(6,*) "ERROR: Invalid mode ",mode
    stop
end select

return
end subroutine mloextra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine mloscrnout(tscrn,qgscrn,uscrn,u10,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: tscrn,qgscrn,uscrn,u10

if (wfull==0) return

tscrn =unpack(dgscrn%temp,wpack,tscrn)
qgscrn=unpack(dgscrn%qg,wpack,qgscrn)
uscrn =unpack(dgscrn%u2,wpack,uscrn)
u10   =unpack(dgscrn%u10,wpack,u10)

return
end subroutine mloscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS + NIR)

subroutine mloalb2(istart,ifin,coszro,ovisalb,oniralb,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisalb,oniralb
real, dimension(wfull) :: watervis,waternir,icevis,icenir
real, dimension(wfull) :: costmp,pond,snow

if (wfull==0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie<ib) return

costmp(ib:ie)=pack(coszro,wpack(istart:ifinish))

!pond(ib:ie)=max(1.+.008*min(ice%tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
!snow(ib:ie)=min(max(ice%snowd(ib:ie)/0.05,0.),1.)
snow(ib:ie)=0.

watervis(ib:ie)=.05/(costmp(ib:ie)+0.15)
waternir(ib:ie)=.05/(costmp(ib:ie)+0.15)
! need to factor snow age into albedo
icevis(ib:ie)=(0.85*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie))+(0.95*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
icenir(ib:ie)=(0.45*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie))+(0.65*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
ovisalb(:)=unpack(icevis(ib:ie)*ice%fracice(ib:ie)+(1.-ice%fracice(ib:ie))*watervis(ib:ie),wpack(istart:ifinish),ovisalb)
oniralb(:)=unpack(icenir(ib:ie)*ice%fracice(ib:ie)+(1.-ice%fracice(ib:ie))*waternir(ib:ie),wpack(istart:ifinish),oniralb)

dgwater%visdiralb(ib:ie)=watervis(ib:ie)
dgwater%visdifalb(ib:ie)=watervis(ib:ie)
dgwater%nirdiralb(ib:ie)=waternir(ib:ie)
dgwater%nirdifalb(ib:ie)=waternir(ib:ie)
dgice%visdiralb(ib:ie)=icevis(ib:ie)
dgice%visdifalb(ib:ie)=icevis(ib:ie)
dgice%nirdiralb(ib:ie)=icenir(ib:ie)
dgice%nirdifalb(ib:ie)=icenir(ib:ie)

return
end subroutine mloalb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS-DIR, VIS-DIF, NIR-DIR & NIR-DIF)

subroutine mloalb4(istart,ifin,coszro,ovisdir,ovisdif,onirdir,onirdif,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisdir,ovisdif,onirdir,onirdif
real, dimension(wfull) :: costmp,pond,snow

if (wfull==0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie<ib) return

costmp(ib:ie)=pack(coszro,wpack(istart:ifinish))

!pond(ib:ie)=max(1.+.008*min(ice%tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
!snow(ib:ie)=min(max(ice%snowd(ib:ie)/0.05,0.),1.)
snow(ib:ie)=0.

where (costmp(ib:ie)>0.)
  dgwater%visdiralb(ib:ie)=0.026/(costmp(ib:ie)**1.7+0.065)+0.15*(costmp(ib:ie)-0.1)* &
                          (costmp(ib:ie)-0.5)*(costmp(ib:ie)-1.)
elsewhere
  dgwater%visdiralb(ib:ie)=0.3925
end where
dgwater%visdifalb(ib:ie)=0.06
dgwater%nirdiralb(ib:ie)=dgwater%visdiralb(ib:ie)
dgwater%nirdifalb(ib:ie)=0.06
! need to factor snow age into albedo
dgice%visdiralb(ib:ie)=(0.85*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie)) &
                        +(0.95*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
dgice%visdifalb(ib:ie)=dgice%visdiralb(ib:ie)
dgice%nirdiralb(ib:ie)=(0.45*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie)) &
                        +(0.65*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
dgice%nirdifalb(ib:ie)=dgice%nirdiralb(ib:ie)
ovisdir=unpack(ice%fracice(ib:ie)*dgice%visdiralb(ib:ie)+(1.-ice%fracice(ib:ie))*dgwater%visdiralb(ib:ie), &
    wpack(istart:ifinish),ovisdir)
ovisdif=unpack(ice%fracice(ib:ie)*dgice%visdifalb(ib:ie)+(1.-ice%fracice(ib:ie))*dgwater%visdifalb(ib:ie), &
    wpack(istart:ifinish),ovisdif)
onirdir=unpack(ice%fracice(ib:ie)*dgice%nirdiralb(ib:ie)+(1.-ice%fracice(ib:ie))*dgwater%nirdiralb(ib:ie), &
    wpack(istart:ifinish),onirdir)
onirdif=unpack(ice%fracice(ib:ie)*dgice%nirdifalb(ib:ie)+(1.-ice%fracice(ib:ie))*dgwater%nirdifalb(ib:ie), &
    wpack(istart:ifinish),onirdif)

return
end subroutine mloalb4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data

subroutine mloregrid(wlin,depin,mloin,mlodat,mode)

implicit none

integer, intent(in) :: wlin,mode
integer iqw,ii,pos(1)
real, dimension(ifull), intent(in) :: depin
real, dimension(ifull,wlin), intent(in) :: mloin
real, dimension(ifull,wlev), intent(inout) :: mlodat
real, dimension(wfull,wlin) :: newdata
real, dimension(wfull,wlev) :: newdatb
real, dimension(wfull) :: deptmp
real, dimension(wlin) :: dpin,sgin
real, dimension(wlin+1) :: dp_hlin
real, dimension(wlev) :: sig
real x

if ( wfull==0 ) return

deptmp=pack(depin,wpack)
do ii=1,wlin
  newdata(:,ii)=pack(mloin(:,ii),wpack)
end do

select case(mode)
  case(0,1)
    do iqw=1,wfull
      call vgrid(wlin,deptmp(iqw),dpin,dp_hlin)
      do ii=1,wlev
        if (depth(iqw,ii)>=dpin(wlin)) then
          newdatb(iqw,ii)=newdata(iqw,wlin)
        else if (depth(iqw,ii)<=dpin(1)) then
          newdatb(iqw,ii)=newdata(iqw,1)
        else
          pos=maxloc(dpin,dpin<depth(iqw,ii))
          pos(1)=max(1,min(wlin-1,pos(1)))
          x=(depth(iqw,ii)-dpin(pos(1)))/(dpin(pos(1)+1)-dpin(pos(1)))
          x=max(0.,min(1.,x))
          newdatb(iqw,ii)=newdata(iqw,pos(1)+1)*x+newdata(iqw,pos(1))*(1.-x)
        end if
      end do
    end do
  case(2,3)
    do iqw=1,wfull
      sig=depth(iqw,:)/depth_hl(iqw,wlev)
      call vgrid(wlin,deptmp(iqw),dpin,dp_hlin)
      sgin=dpin/dp_hlin(wlin)
      do ii=1,wlev
        if (sig(ii)>=sgin(wlin)) then
          newdatb(iqw,ii)=newdata(iqw,wlin)
        else if (sig(ii)<=sgin(1)) then
          newdatb(iqw,ii)=newdata(iqw,1)
        else
          pos=maxloc(sgin,sgin<sig(ii))
          pos(1)=max(1,min(wlin-1,pos(1)))
          x=(sig(ii)-sgin(pos(1)))/(sgin(pos(1)+1)-sgin(pos(1)))
          x=max(0.,min(1.,x))
          newdatb(iqw,ii)=newdata(iqw,pos(1)+1)*x+newdata(iqw,pos(1))*(1.-x)
        end if
      end do
    end do
  case default
    write(6,*) "ERROR: Unknown mode for mloregrid ",mode
    stop
end select

do ii=1,wlev
  mlodat(:,ii)=unpack(newdatb(:,ii),wpack,mlodat(:,ii))
end do

return
end subroutine mloregrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ocean depth data

subroutine mloexpdep(mode,odep,ii,diag)

implicit none

integer, intent(in) :: mode,ii,diag
real, dimension(ifull), intent(out) :: odep

odep=0.
if (wfull==0) return

select case(mode)
  case(0)
    odep=unpack(depth(:,ii),wpack,odep)
  case(1)
    odep=unpack(dz(:,ii),wpack,odep)
  case default
    write(6,*) "ERROR: Unknown mloexpdep mode ",mode
    stop
end select
  
return
end subroutine mloexpdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract density

subroutine mloexpdensity(odensity,alpha,beta,tt,ss,ddz,pxtr,diag)

implicit none

integer, intent(in) :: diag
real, dimension(:,:), intent(in) :: tt
real, dimension(size(tt,1)), intent(in) :: pxtr
real, dimension(size(tt,1),size(tt,2)), intent(in) :: ss,ddz
real, dimension(size(tt,1),size(tt,2)), intent(out) :: odensity,alpha,beta
real, dimension(size(tt,1)) :: rho0

call calcdensity(odensity,alpha,beta,rho0,tt,ss,ddz,pxtr)

return
end subroutine mloexpdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract melting temperature

subroutine mloexpmelt(omelt)

implicit none

real, dimension(ifull), intent(out) :: omelt
real, dimension(wfull) :: tmelt,d_zcr

omelt=273.16
if (wfull==0) return
d_zcr=max(1.+water%eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))
call calcmelt(tmelt,d_zcr)
omelt=unpack(tmelt,wpack,omelt)

return
end subroutine mloexpmelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract factors for energy conservation

subroutine mloexpgamm(gamm,ip_dic,ip_dsn,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(in) :: ip_dic, ip_dsn
real, dimension(ifull,3), intent(out) :: gamm

gamm(:,1)=gammi
gamm(:,2)=max(ip_dsn,0.)*cps
gamm(:,3)=max(ip_dic,0.)*0.5*cpi

return
end subroutine mloexpgamm

!subroutine mloexpenergy(engout,diag)
!
!implicit none
!
!integer, intent(in) :: diag
!real, dimension(ifull), intent(out) :: engout
!
!engout=unpack(real(dgwater%deleng+dgice%deleng),wpack,0.)
!
!return
!end subroutine mloexpenergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

subroutine mloeval(sst,zo,cd,cds,fg,eg,wetfac,epot,epan,fracice,siced,snowd, &
                   dt,zmin,zmins,sg,rg,precp,precs,uatm,vatm,temp,qg,ps,f,   &
                   visnirratio,fbvis,fbnir,inflow,diag,calcprog,oldu,oldv)

implicit none

integer, intent(in) :: diag
real, intent(in) :: dt
real, dimension(ifull), intent(in) :: sg,rg,precp,precs,f,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,inflow,zmin,zmins
real, dimension(ifull), intent(inout) :: sst,zo,cd,cds,fg,eg,wetfac,fracice,siced,epot,epan,snowd
real, dimension(ifull), intent(in), optional :: oldu,oldv
real, dimension(wfull) :: workb,workc
real, dimension(wfull) :: atm_sg,atm_rg,atm_rnd,atm_snd,atm_f,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v,atm_temp,atm_qg
real, dimension(wfull) :: atm_ps,atm_zmin,atm_zmins,atm_inflow,atm_oldu,atm_oldv
real, dimension(wfull,wlev) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_zcr
real, dimension(wfull) :: d_fb,d_timelt,d_neta,d_ndsn
real, dimension(wfull) :: d_ndic,d_nsto,d_delstore,d_delinflow,d_imass
!real, dimension(wfull,wlev) :: oldwatertemp
!real, dimension(wfull,0:3) :: oldicetemp
!real, dimension(wfull) :: oldicestore,oldicesnowd,oldicethick,oldicefrac,oldzcr
integer, dimension(wfull) :: d_nk
logical, intent(in) :: calcprog ! flag to update prognostic variables (or just calculate fluxes)

if (wfull==0) return

atm_sg     =pack(sg,wpack)
atm_rg     =pack(rg,wpack)
atm_f      =pack(f,wpack)
atm_vnratio=pack(visnirratio,wpack)
atm_fbvis  =pack(fbvis,wpack)
atm_fbnir  =pack(fbnir,wpack)
atm_u      =pack(uatm,wpack)
atm_v      =pack(vatm,wpack)
atm_temp   =pack(temp,wpack)
atm_qg     =pack(qg,wpack)
atm_ps     =pack(ps,wpack)
atm_zmin   =pack(zmin,wpack)
atm_zmins  =pack(zmins,wpack)
atm_inflow =pack(inflow,wpack)
atm_rnd    =pack(precp-precs,wpack)
atm_snd    =pack(precs,wpack)
if (present(oldu).and.present(oldv)) then
  atm_oldu=pack(oldu,wpack)
  atm_oldv=pack(oldv,wpack)
else
  atm_oldu=water%u(:,1)
  atm_oldv=water%v(:,1)
end if

! adjust levels for free surface
d_zcr=max(1.+water%eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))

! store state variables for energy conservation check
!oldwatertemp=water%temp
!oldicetemp(:,0)=ice%tsurf
!oldicetemp(:,1:3)=ice%temp(:,0:2)
!oldicestore=ice%store
!oldicesnowd=ice%snowd
!oldicethick=ice%thick
!oldicefrac=ice%fracice
!oldzcr=d_zcr

! calculate melting temperature
call calcmelt(d_timelt,d_zcr)

! equation of state
call getrho(atm_ps,d_rho,d_alpha,d_beta,d_zcr)

! ice mass per unit area
! MJT notes - a limit of 10 can cause reproducibility issues with
! single precision and multple processes
d_imass=max(rhoic*ice%thick+rhosn*ice%snowd,100.) 

! split adjustment of free surface and ice thickness to ensure conservation
d_ndsn=ice%snowd
d_ndic=ice%thick
d_nsto=ice%store
where ( d_ndic>icemin )
  d_delinflow=qsnow*0.001*(atm_rnd+atm_snd)
  d_ndsn=d_ndsn+0.001*dt*(atm_rnd+atm_snd)
  d_neta=water%eta+0.001*dt*(atm_inflow+(1.-ice%fracice)*(atm_rnd+atm_snd))
elsewhere
  d_neta=water%eta+0.001*dt*(atm_inflow+atm_rnd+atm_snd)
  d_delinflow=0.
end where

! water fluxes
call fluxcalc(dt,atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,d_neta,atm_oldu,atm_oldv,     &
              diag)

! boundary conditions
call getwflux(atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow,d_rho,d_nsq,  &
              d_rad,d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr)

! ice fluxes
call iceflux(dt,atm_sg,atm_rg,atm_rnd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v,atm_temp,atm_qg,   &
             atm_ps,atm_zmin,atm_zmins,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_ndsn,d_ndic,d_nsto,          &
             d_delstore,d_imass,atm_oldu,atm_oldv,diag)

if (calcprog) then

  ! update ice
  ice%tsurf=ice%tsurf+ice%temp(:,0)*cps*(ice%snowd-d_ndsn)/gammi                     ! close energy budget
  ice%tsurf=ice%tsurf+0.5*(ice%temp(:,1)+ice%temp(:,2))*cpi*(ice%thick-d_ndic)/gammi ! close energy budget
  ice%thick=d_ndic
  ice%snowd=d_ndsn
  ice%store=d_nsto
  call mloice(dt,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt,       &
              d_ustar,d_nk,d_neta,d_imass,diag)

  ! create or destroy ice
  ! MJT notes - this is done after the flux calculations to agree with the albedo passed to the radiation
  call mlonewice(d_timelt,d_zcr,diag)
  
  ! update water
  call mlocalc(dt,atm_f,d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr, &
               d_neta,diag)

end if
! screen diagnostics
call scrncalc(atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,diag)

! energy conservation check
!d_zcr=max(1.+water%eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))
!dgwater%deleng=0._8
!do ii=1,wlev
!  dgwater%deleng=dgwater%deleng+real((water%temp(:,ii)*d_zcr-oldwatertemp(:,ii)*oldzcr)*dz(:,ii),8)
!end do
!dgwater%deleng=dgwater%deleng*real(rhowt*cp0/dt,8)
!dgice%deleng=real(ice%store-oldicestore,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%tsurf-oldicefrac*oldicetemp(:,0))*gammi,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,0)*ice%snowd-oldicefrac*oldicetemp(:,1)*oldicesnowd)*cps,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,1)*ice%thick-oldicefrac*oldicetemp(:,2)*oldicethick)*0.5*cpi,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,2)*ice%thick-oldicefrac*oldicetemp(:,3)*oldicethick)*0.5*cpi,8)
!dgice%deleng=dgice%deleng/real(dt,8)

workb=emisice**0.25*ice%tsurf
sst    =unpack((1.-ice%fracice)*(water%temp(:,1)+wrtemp)+ice%fracice*workb,wpack,sst)
workc=(1.-ice%fracice)/log(atm_zmin/dgwater%zo)**2+ice%fracice/log(atm_zmin/dgice%zo)**2
zo     =unpack(atm_zmin*exp(-1./sqrt(workc)),wpack,zo)
cd     =unpack((1.-ice%fracice)*dgwater%cd +ice%fracice*dgice%cd,wpack,cd)
cds    =unpack((1.-ice%fracice)*dgwater%cdh+ice%fracice*dgice%cdh,wpack,cds)
fg     =unpack((1.-ice%fracice)*dgwater%fg +ice%fracice*dgice%fg,wpack,fg)
eg     =unpack((1.-ice%fracice)*dgwater%eg +ice%fracice*dgice%eg,wpack,eg)
wetfac =unpack((1.-ice%fracice)            +ice%fracice*dgice%wetfrac,wpack,wetfac)
epan   =unpack(dgwater%eg,wpack,epan)
epot   =unpack((1.-ice%fracice)*dgwater%eg +ice%fracice*dgice%eg/max(dgice%wetfrac,1.e-20),wpack,epot)
fracice=0.
fracice=unpack(ice%fracice,wpack,fracice)
siced=0.
siced  =unpack(ice%thick,wpack,siced)
snowd  =unpack(ice%snowd,wpack,snowd)

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,atm_f,d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr, &
                   d_neta,diag)

implicit none

integer, intent(in) :: diag
integer ii, iqw
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km, ks, gammas, rhs
real(kind=8), dimension(wfull,2:wlev) :: aa
real(kind=8), dimension(wfull,wlev) :: bb, dd
real(kind=8), dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull,wlev), intent(in) :: d_rho, d_nsq, d_rad, d_alpha
real, dimension(wfull) :: dumt0, umag, avearray
real, dimension(wfull), intent(in) :: atm_f
real, dimension(wfull), intent(inout) :: d_b0, d_ustar, d_wu0, d_wv0, d_wt0, d_ws0, d_zcr, d_neta

! solve for mixed layer depth (calculated at full levels)
call getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,atm_f,d_zcr) 
! solve for stability functions and non-local term (calculated at half levels)
call getstab(km,ks,gammas,d_nsq,d_ustar,d_zcr) 


! Counter-gradient term for scalars (rhs)
! +ve sign for rhs terms since z +ve is down
rhs(:,1)=ks(:,2)*gammas(:,2)/(dz(:,1)*d_zcr)
do ii=2,wlev-1
  rhs(:,ii)=(ks(:,ii+1)*gammas(:,ii+1)-ks(:,ii)*gammas(:,ii))/(dz(:,ii)*d_zcr)
end do
rhs(:,wlev)=-ks(:,wlev)*gammas(:,wlev)/(dz(:,wlev)*d_zcr)


! Diffusion term for scalars (aa,bb,cc)
cc(:,1)=-dt*ks(:,2)/(dz_hl(:,2)*dz(:,1)*d_zcr*d_zcr)
bb(:,1)=1._8-cc(:,1)
do ii=2,wlev-1
  aa(:,ii)=-dt*ks(:,ii)/(dz_hl(:,ii)*dz(:,ii)*d_zcr*d_zcr)
  cc(:,ii)=-dt*ks(:,ii+1)/(dz_hl(:,ii+1)*dz(:,ii)*d_zcr*d_zcr)
  bb(:,ii)=1._8-aa(:,ii)-cc(:,ii)
end do
aa(:,wlev)=-dt*ks(:,wlev)/(dz_hl(:,wlev)*dz(:,wlev)*d_zcr*d_zcr)
bb(:,wlev)=1._8-aa(:,wlev)


! POTENTIAL TEMPERATURE
if ( incradgam>0 ) then
  ! include radiation in counter-gradient term
  do iqw=1,wfull
    dumt0(iqw)=d_wt0(iqw)+sum(d_rad(iqw,1:dgwater%mixind(iqw)))
  end do
else
  dumt0=d_wt0
end if
avearray=0.5*(minval(water%temp,dim=2)+maxval(water%temp,dim=2))
do ii=1,wlev
  water%temp(:,ii)=water%temp(:,ii)-avearray
  dd(:,ii)=water%temp(:,ii)+dt*rhs(:,ii)*dumt0-dt*d_rad(:,ii)/(dz(:,ii)*d_zcr)
end do
dd(:,1)=dd(:,1)-dt*d_wt0/(dz(:,1)*d_zcr)
call thomas(water%temp,aa,bb,cc,dd)
do ii=1,wlev
  water%temp(:,ii)=water%temp(:,ii)+avearray
end do


! SALINITY
avearray=0.5*(maxval(water%sal,dim=2)+minval(water%sal,dim=2))
do ii=1,wlev
  water%sal(:,ii)=water%sal(:,ii)-avearray
end do
do ii=1,wlev
  dd(:,ii)=water%sal(:,ii)+dt*rhs(:,ii)*d_ws0
end do
dd(:,1)=dd(:,1)-dt*d_ws0/(dz(:,1)*d_zcr)
call thomas(water%sal,aa,bb,cc,dd)
do ii=1,wlev
  water%sal(:,ii)=water%sal(:,ii)+avearray
end do
water%sal=max(0.,water%sal)


! Diffusion term for momentum (aa,bb,cc)
cc(:,1)=-dt*km(:,2)/(dz_hl(:,2)*dz(:,1)*d_zcr*d_zcr)
bb(:,1)=1._8-cc(:,1)
do ii=2,wlev-1
  aa(:,ii)=-dt*km(:,ii)/(dz_hl(:,ii)*dz(:,ii)*d_zcr*d_zcr)
  cc(:,ii)=-dt*km(:,ii+1)/(dz_hl(:,ii+1)*dz(:,ii)*d_zcr*d_zcr)
  bb(:,ii)=1._8-aa(:,ii)-cc(:,ii)
end do
aa(:,wlev)=-dt*km(:,wlev)/(dz_hl(:,wlev)*dz(:,wlev)*d_zcr*d_zcr)
bb(:,wlev)=1._8-aa(:,wlev)
umag=sqrt(water%u(:,wlev)*water%u(:,wlev)+water%v(:,wlev)*water%v(:,wlev))
! bottom drag
where ( depth_hl(:,wlev+1)<mxd )
  bb(:,wlev)=bb(:,wlev)+dt*cdbot*umag/(dz(:,wlev)*d_zcr)
end where


! U diffusion term
do ii=1,wlev
  dd(:,ii)=water%u(:,ii)
end do
dd(:,1)=dd(:,1)-dt*d_wu0/(dz(:,1)*d_zcr) ! explicit
call thomas(water%u,aa,bb,cc,dd)


! V diffusion term
do ii=1,wlev
  dd(:,ii)=water%v(:,ii)
end do
dd(:,1)=dd(:,1)-dt*d_wv0/(dz(:,1)*d_zcr) ! explicit
call thomas(water%v,aa,bb,cc,dd)


! --- Turn off coriolis terms as this is processed in mlodynamics.f90 ---
!! Split U and V coriolis terms
!xp=1.+(0.5*dt*atm_f)**2
!xm=1.-(0.5*dt*atm_f)**2
!do ii=1,wlev
!  newa=(water%u(:,ii)*xm+water%v(:,ii)*dt*atm_f)/xp
!  newb=(water%v(:,ii)*xm-water%u(:,ii)*dt*atm_f)/xp
!  water%u(:,ii)=newa
!  water%v(:,ii)=newb
!end do


! adjust surface height
select case(deprelax)
  case(0) ! free surface height
    water%eta=d_neta
  case(1) ! relax surface height
    water%eta=d_neta-dt*d_neta/(3600.*24.)
  case(2) ! fix surface height
    water%eta=0.
  case DEFAULT
    write(6,*) "ERROR: Invalid deprelax ",deprelax
    stop
end select

return
end subroutine mlocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for stability functions
! GFDL use a look-up table to speed-up the code...

subroutine getstab(km,ks,gammas,d_nsq,d_ustar,d_zcr)

implicit none

integer ii,jj,iqw
integer, dimension(wfull) :: mixind_hl
real, dimension(wfull,wlev), intent(out) :: km,ks,gammas
real, dimension(wfull,wlev) :: num,nus,wm,ws,ri
real, dimension(wfull) :: sigma,d_depth_hl,d_depth_hlp1
real, dimension(wfull) :: a2m,a3m,a2s,a3s
real, dimension(wfull) :: numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real, dimension(wfull) :: nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real xp,cg
real, dimension(wfull,wlev), intent(in) :: d_nsq
real, dimension(wfull), intent(in) :: d_ustar
real, dimension(wfull), intent(in) :: d_zcr
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.0E-4
real, parameter :: nusw = 0.1E-4

ri=0. ! ri(:,1) is not used
do ii=2,wlev
  ri(:,ii)=d_nsq(:,ii)*(dz_hl(:,ii)*d_zcr)**2    & 
           /max((water%u(:,ii-1)-water%u(:,ii))**2       &
               +(water%v(:,ii-1)-water%v(:,ii))**2,1.E-20)
end do

! diffusion ---------------------------------------------------------
do ii=2,wlev
  ! s - shear
  where (ri(:,ii)<0.)
    num(:,ii)=nu0
  elsewhere (ri(:,ii)<ri0)
    num(:,ii)=nu0*(1.-(ri(:,ii)/ri0)**2)**3
  elsewhere
    num(:,ii)=0.
  endwhere
  ! d - double-diffusive
  ! double-diffusive mixing is neglected for now
  ! w - internal ave
  num(:,ii)=num(:,ii)+numw
  nus(:,ii)=num(:,ii)+nusw ! note that nus is a copy of num
end do
num(:,1)=num(:,2) ! to avoid problems calculating shallow mixed layer
nus(:,1)=nus(:,2)
!--------------------------------------------------------------------

! stability ---------------------------------------------------------
wm=0.
ws=0.
do ii=2,wlev
  d_depth_hl=depth_hl(:,ii)*d_zcr
  call getwx(wm(:,ii),ws(:,ii),d_depth_hl,dgwater%bf,d_ustar,dgwater%mixdepth)
end do
wm(:,1)=wm(:,2) ! to avoid problems calculating shallow mixed layer
ws(:,1)=ws(:,2)
!--------------------------------------------------------------------

! calculate G profile -----------------------------------------------
! -ve as z is down
do iqw=1,wfull
  mixind_hl(iqw)=dgwater%mixind(iqw)
  jj=min(mixind_hl(iqw)+1,wlev-1)
  d_depth_hl(iqw)=depth_hl(iqw,jj)*d_zcr(iqw)
  if (dgwater%mixdepth(iqw)>d_depth_hl(iqw)) mixind_hl(iqw)=jj
  d_depth_hl(iqw)=depth_hl(iqw,mixind_hl(iqw))*d_zcr(iqw)
  d_depth_hlp1(iqw)=depth_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw)
  xp=(dgwater%mixdepth(iqw)-d_depth_hl(iqw))/(d_depth_hlp1(iqw)-d_depth_hl(iqw))
  xp=max(0.,min(1.,xp))
  numh(iqw)=(1.-xp)*num(iqw,mixind_hl(iqw))+xp*num(iqw,mixind_hl(iqw)+1)
  nush(iqw)=(1.-xp)*nus(iqw,mixind_hl(iqw))+xp*nus(iqw,mixind_hl(iqw)+1)
  wm1(iqw)=(1.-xp)*wm(iqw,mixind_hl(iqw))+xp*wm(iqw,mixind_hl(iqw)+1)
  ws1(iqw)=(1.-xp)*ws(iqw,mixind_hl(iqw))+xp*ws(iqw,mixind_hl(iqw)+1)
  dnumhdz(iqw)=min(num(iqw,mixind_hl(iqw)+1)-num(iqw,mixind_hl(iqw)),0.)/(dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
  dnushdz(iqw)=min(nus(iqw,mixind_hl(iqw)+1)-nus(iqw,mixind_hl(iqw)),0.)/(dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
  !dwm1ds and dws1ds are now multipled by 1/(mixdepth*wx1*wx1)
  dwm1ds(iqw)=-5.*max(dgwater%bf(iqw),0.)/d_ustar(iqw)**4
  dws1ds(iqw)=dwm1ds(iqw)
end do

g1m=numh/(dgwater%mixdepth*wm1)
g1s=nush/(dgwater%mixdepth*ws1)
dg1mds=dnumhdz/wm1-numh*dwm1ds
dg1sds=dnushdz/ws1-nush*dws1ds
  
a2m=-2.+3.*g1m-dg1mds
a2s=-2.+3.*g1s-dg1sds
a3m=1.-2.*g1m+dg1mds
a3s=1.-2.*g1s+dg1sds

!--------------------------------------------------------------------
! calculate diffusion terms
km(:,1)=0.
ks(:,1)=0.
do ii=2,wlev
  where (ii<=mixind_hl)
    sigma=depth_hl(:,ii)*d_zcr/dgwater%mixdepth
    km(:,ii)=max(dgwater%mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma)),num(:,ii))
    ks(:,ii)=max(dgwater%mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2s+a3s*sigma)),nus(:,ii))
  elsewhere
    km(:,ii)=num(:,ii)
    ks(:,ii)=nus(:,ii)
  end where
end do

!--------------------------------------------------------------------
! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Large (1994)
!cg=5.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Bernie (2004)
gammas(:,1)=0.
do ii=2,wlev
  where (dgwater%bf<0..and.ii<=mixind_hl) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*dgwater%mixdepth,1.E-20)
  elsewhere
    gammas(:,ii)=0.
  end where
end do

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixed layer depth

subroutine getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,atm_f,d_zcr)

implicit none

integer ii,jj,kk,iqw
integer, dimension(wfull) :: isf
real vtc,dvsq,vtsq,xp
real tnsq,tws,twu,twv,tdepth,tbuoy,trho
real oldxp,oldtrib,trib,newxp,deldz,aa,bb,dsf
real, dimension(wfull,wlev) :: ws,wm,dumbuoy,rib
real, dimension(wfull) :: dumbf,l,d_depth,usf,vsf,rsf
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha
real, dimension(wfull), intent(in) :: d_b0,d_ustar,d_zcr
real, dimension(wfull), intent(in) :: atm_f
integer, parameter :: maxits = 3

vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar*vkar*ric)

! Modify buoyancy forcing with solar radiation
if (incradbf>0) then
  do ii=1,wlev
    dumbf=d_b0-grav*sum(d_alpha(:,1:ii)*d_rad(:,1:ii),2) ! -ve sign is to account for sign of d_rad
    d_depth=depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
else
  dumbf=d_b0
  do ii=1,wlev
    d_depth=depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
end if

! Estimate surface layer values
usf=0.
vsf=0.
rsf=0.
dsf=0.
isf=wlev-1
do iqw=1,wfull
  dsf=0.
  do ii=1,wlev-1
    aa=dz(iqw,ii)*d_zcr(iqw)
    bb=max(minsfc-dsf,0.)
    deldz=min(aa,bb)
    usf(iqw)=usf(iqw)+water%u(iqw,ii)*deldz
    vsf(iqw)=vsf(iqw)+water%v(iqw,ii)*deldz
    rsf(iqw)=rsf(iqw)+d_rho(iqw,ii)*deldz
    dsf=dsf+deldz
    if (bb<=0.) exit
  end do
  if (depth(iqw,ii)*d_zcr(iqw)>minsfc) then
    isf(iqw)=ii
  else
    isf(iqw)=ii+1
  end if
  usf(iqw)=usf(iqw)/dsf
  vsf(iqw)=vsf(iqw)/dsf
  rsf(iqw)=rsf(iqw)/dsf
end do

! Calculate local buoyancy
dumbuoy=0.
do iqw=1,wfull
  do ii=isf(iqw),wlev
    dumbuoy(iqw,ii)=grav*(d_rho(iqw,ii)-rsf(iqw))
  end do
end do

! Calculate mixed layer depth from critical Ri
dgwater%mixind=wlev-1
dgwater%mixdepth=depth(:,wlev)*d_zcr
rib=0.
do iqw=1,wfull
  do ii=isf(iqw),wlev
    jj=min(ii+1,wlev)
    vtsq=depth(iqw,ii)*d_zcr(iqw)*ws(iqw,ii)*sqrt(0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.))*vtc
    dvsq=(usf(iqw)-water%u(iqw,ii))**2+(vsf(iqw)-water%v(iqw,ii))**2
    rib(iqw,ii)=(depth(iqw,ii)*d_zcr(iqw)-minsfc)*dumbuoy(iqw,ii)/(max(dvsq+vtsq,1.E-20)*d_rho(iqw,ii))
    if (rib(iqw,ii)>ric) then
      jj=max(ii-1,1)
      dgwater%mixind(iqw)=jj
      xp=min(max((ric-rib(iqw,jj))/max(rib(iqw,ii)-rib(iqw,jj),1.E-20),0.),1.)
      dgwater%mixdepth(iqw)=((1.-xp)*depth(iqw,jj)+xp*depth(iqw,ii))*d_zcr(iqw)
      exit
    end if
  end do 
end do

! Refine mixed-layer-depth calculation by improving vertical profile of buoyancy
if (mixmeth==1) then
  do iqw=1,wfull
    ii=dgwater%mixind(iqw)+1
    jj=min(ii+1,wlev)
    oldxp=0.
    oldtrib=rib(iqw,ii-1)
    xp=(dgwater%mixdepth(iqw)-depth(iqw,ii-1)*d_zcr(iqw))/((depth(iqw,ii)-depth(iqw,ii-1))*d_zcr(iqw))
    do kk=1,maxits
      if (xp<0.5) then
        tnsq=(1.-2.*xp)*0.5*max(d_nsq(iqw,ii-1)+d_nsq(iqw,ii),0.)+(2.*xp)*max(d_nsq(iqw,ii),0.)
      else
        tnsq=(2.-2.*xp)*max(d_nsq(iqw,ii),0.)+(2.*xp-1.)*0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.)
      end if
      tws=(1.-xp)*ws(iqw,ii-1)+xp*ws(iqw,ii)
      twu=(1.-xp)*water%u(iqw,ii-1)+xp*water%u(iqw,ii)
      twv=(1.-xp)*water%v(iqw,ii-1)+xp*water%v(iqw,ii)
      tdepth=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
      tbuoy=(1.-xp)*dumbuoy(iqw,ii-1)+xp*dumbuoy(iqw,ii)
      trho=(1.-xp)*d_rho(iqw,ii-1)+xp*d_rho(iqw,ii)
      vtsq=tdepth*tws*sqrt(tnsq)*vtc
      dvsq=(usf(iqw)-twu)**2+(vsf(iqw)-twv)**2
      trib=(tdepth-minsfc)*tbuoy/(max(dvsq+vtsq,1.E-20)*trho)
      if (abs(trib-oldtrib)<1.E-5) exit
      newxp=xp-(trib-ric)*(xp-oldxp)/(trib-oldtrib) ! i.e., (trib-ric-oldtrib+ric)
      oldtrib=trib
      oldxp=xp
      xp=newxp
      xp=min(max(xp,0.),1.)
    end do
    dgwater%mixdepth(iqw)=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
  end do
end if

! calculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

! impose limits for stable conditions
where(dgwater%bf>1.E-10.and.abs(atm_f)>1.E-10)
  l=0.7*d_ustar/abs(atm_f)
elsewhere (dgwater%bf>1.E-10)
  l=d_ustar*d_ustar*d_ustar/(vkar*dgwater%bf)
elsewhere
  l=depth(:,wlev)*d_zcr
end where
dgwater%mixdepth=min(dgwater%mixdepth,l)
dgwater%mixdepth=max(dgwater%mixdepth,depth(:,1)*d_zcr)
dgwater%mixdepth=min(dgwater%mixdepth,depth(:,wlev)*d_zcr)

! recalculate index for mixdepth
dgwater%mixind=wlev-1
do iqw=1,wfull
  do ii=2,wlev
    if (depth(iqw,ii)*d_zcr(iqw)>dgwater%mixdepth(iqw)) then
      dgwater%mixind(iqw)=ii-1
      exit
    end if
  end do
end do

! recalculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

return
end subroutine getmixdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate bouyancy forcing

subroutine getbf(d_rad,d_alpha,d_b0)

implicit none

integer iqw
real, dimension(wfull,wlev), intent(in) :: d_rad,d_alpha
real, dimension(wfull), intent(in) :: d_b0

if (incradbf>0) then
  do iqw=1,wfull
    ! -ve sign is to account for sign of d_rad
    dgwater%bf(iqw)=d_b0(iqw)-grav*sum(d_alpha(iqw,1:dgwater%mixind(iqw))*d_rad(iqw,1:dgwater%mixind(iqw))) 
  end do
else
  dgwater%bf(:)=d_b0(:)
end if

return
end subroutine getbf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the stability functions

subroutine getwx(wm,ws,dep,bf,d_ustar,mixdp)

implicit none

real, dimension(wfull), intent(out) :: wm,ws
real, dimension(wfull), intent(in) :: bf,mixdp,dep
real, dimension(wfull) :: zeta,sig,invl,uuu
real, dimension(wfull), intent(in) :: d_ustar
real, parameter :: zetam=-0.2
real, parameter :: zetas=-1.0
real, parameter :: am=1.26
real, parameter :: cm=8.38
real, parameter :: as=-28.86
real, parameter :: cs=98.96

sig=dep/mixdp                       ! stable
where (bf<=0..and.sig>epsilon) ! unstable
  sig=epsilon
end where
uuu=d_ustar**3
invl=vkar*bf      ! invl = ustar**3/L or L=ustar**3/(vkar*bf)
zeta=sig*mixdp*invl
zeta=min(zeta,uuu) ! MJT suggestion

where (zeta>0.)
  wm=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta>zetam*uuu)
  wm=vkar*(d_ustar*uuu-16.*d_ustar*zeta)**(1./4.)
elsewhere
  wm=vkar*(am*uuu-cm*zeta)**(1./3.)
end where

where (zeta>0.)
  ws=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta>zetas*uuu)
  ws=vkar*(d_ustar*d_ustar-16.*zeta/d_ustar)**(1./2.)
elsewhere
  ws=vkar*(as*uuu-cs*zeta)**(1./3.)
end where

wm=max(wm,1.E-20)
ws=max(ws,1.E-20)

return
end subroutine getwx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rho from equation of state
! From GFDL (MOM3)

subroutine getrho(atm_ps,d_rho,d_alpha,d_beta,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: rho0,pxtr
real, dimension(wfull,wlev) :: d_dz
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_alpha,d_beta
real, dimension(wfull), intent(in) :: atm_ps
real, dimension(wfull), intent(inout) :: d_zcr

pxtr=atm_ps+grav*ice%fracice*(ice%thick*rhoic+ice%snowd*rhosn)
do ii=1,wlev
  d_dz(:,ii)=dz(:,ii)*d_zcr
end do
call calcdensity(d_rho,d_alpha,d_beta,rho0,water%temp,water%sal,d_dz,pxtr)

return
end subroutine getrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate surface density (used for sea-ice melting)

!subroutine getrho1(salin,atm_ps,d_rho,d_zcr)
!
!implicit none
!
!integer ii
!real, dimension(wfull) :: rho0,pxtr
!real, dimension(wfull,1) :: d_dz,d_alpha,d_beta,d_isal
!real, dimension(wfull,1), intent(inout) :: d_rho
!real, dimension(wfull), intent(in) :: atm_ps,salin
!real, dimension(wfull), intent(inout) :: d_zcr
!
!pxtr=atm_ps+grav*ice%fracice*(ice%thick*rhoic+ice%snowd*rhosn)
!d_dz(:,1)=dz(:,1)*d_zcr
!d_isal(:,1)=salin
!call calcdensity(d_rho(:,1:1),d_alpha(:,1:1),d_beta(:,1:1),rho0,water%temp(:,1:1),d_isal(:,1:1),d_dz(:,1:1),pxtr)
!
!return
!end subroutine getrho1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate water boundary conditions

subroutine getwflux(atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow,d_rho,d_nsq,d_rad,d_alpha, &
                    d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: visalb,niralb,netvis,netnir
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr
real, dimension(wfull), intent(in) :: atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow

! buoyancy frequency (calculated at half levels)
do ii=2,wlev
  d_nsq(:,ii)=-2.*grav/(d_rho(:,ii-1)+d_rho(:,ii))*(d_rho(:,ii-1)-d_rho(:,ii))/(dz_hl(:,ii)*d_zcr)
end do
d_nsq(:,1)=2.*d_nsq(:,2)-d_nsq(:,3) ! not used

! shortwave
! use -ve as depth is down
visalb=dgwater%visdiralb*atm_fbvis+dgwater%visdifalb*(1.-atm_fbvis)
niralb=dgwater%nirdiralb*atm_fbnir+dgwater%nirdifalb*(1.-atm_fbnir)
netvis=(1.-visalb)*atm_vnratio
netnir=(1.-niralb)*(1.-atm_vnratio)
d_rad(:,1)=netvis*(exp(-depth_hl(:,2)*d_zcr/mu_1)-1.) &
          +netnir*(exp(-depth_hl(:,2)*d_zcr/mu_2)-1.)
do ii=2,wlev-1 
  d_rad(:,ii)=netvis*(exp(-depth_hl(:,ii+1)*d_zcr/mu_1)-exp(-depth_hl(:,ii)*d_zcr/mu_1)) &
             +netnir*(exp(-depth_hl(:,ii+1)*d_zcr/mu_2)-exp(-depth_hl(:,ii)*d_zcr/mu_2))
end do
d_rad(:,wlev)=-netvis*exp(-depth_hl(:,wlev)*d_zcr/mu_1) &
              -netnir*exp(-depth_hl(:,wlev)*d_zcr/mu_2) ! remainder
do ii=1,wlev
  d_rad(:,ii)=d_rad(:,ii)*(1.-ice%fracice)*atm_sg/(cp0*rhowt)
end do

! Boundary conditions
! MJT notes - use rhowt reference density for Boussinesq fluid approximation
d_wu0=-(1.-ice%fracice)*dgwater%taux/rhowt
d_wv0=-(1.-ice%fracice)*dgwater%tauy/rhowt
d_wt0=-(1.-ice%fracice)*(-dgwater%fg-dgwater%eg+atm_rg-sbconst*(water%temp(:,1)+wrtemp)**4)/(rhowt*cp0)
d_wt0=d_wt0+(1.-ice%fracice)*lf*atm_snd/(rhowt*cp0) ! melting snow
d_ws0=(1.-ice%fracice)*(atm_rnd+atm_snd-dgwater%eg/lv)*water%sal(:,1)/rhowt
d_ws0=d_ws0+atm_inflow*water%sal(:,1)/rhowt ! inflow under ice

d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-6)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign to compensate for sign of wt0 and ws0
                                                  ! This is the opposite sign used by Large for wb0, but
                                                  ! is same sign as Large used for Bf (+ve stable, -ve unstable)

return
end subroutine getwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate density

subroutine calcdensity(d_rho,d_alpha,d_beta,rho0,tt,ss,ddz,pxtr)

implicit none

integer wlx,ii
!integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(:,:), intent(in) :: tt
real, dimension(size(tt,1),size(tt,2)), intent(in) :: ss,ddz
real, dimension(size(tt,1),size(tt,2)), intent(out) :: d_rho,d_alpha,d_beta
real, dimension(size(tt,1)), intent(in) :: pxtr
real, dimension(size(tt,1)), intent(out) :: rho0
real, dimension(size(tt,1)) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32,ptot
real, dimension(size(tt,1)) :: drho0dt,drho0ds,dskdt,dskds,sk,sks
real, dimension(size(tt,1)) :: drhodt,drhods,rs0
real, parameter :: density = 1035.

wlx =size(tt,2)
d_rho=density

t = min(max(tt(:,1)-273.16+wrtemp,-2.2),100.)
s = min(max(ss(:,1),0.),maxsal) ! limit max salinity for equation of state
t2 = t*t
t3 = t2*t
t4 = t3*t
t5 = t4*t
s2 = s*s
s3 = s2*s
s32 = sqrt(s3)

rs0 = 999.842594 + 6.793952e-2*t(:)                                     &
       - 9.095290e-3*t2(:) + 1.001685e-4*t3(:)                          &
       - 1.120083e-6*t4(:) + 6.536332e-9*t5(:) ! density for sal=0.
rho0 = rs0+ s(:)*(0.824493 - 4.0899e-3*t(:)                             &
       + 7.6438e-5*t2(:)                                                &
       - 8.2467e-7*t3(:) + 5.3875e-9*t4(:))                             &
       + s32(:)*(-5.72466e-3 + 1.0227e-4*t(:)                           &
       - 1.6546e-6*t2(:)) + 4.8314e-4*s2(:)     ! + sal terms    
drho0dt=6.793952e-2                                                     &
       - 2.*9.095290e-3*t(:) + 3.*1.001685e-4*t2(:)                     &
       - 4.*1.120083e-6*t3(:) + 5.*6.536332e-9*t4(:)                    &
       + s(:)*( - 4.0899e-3 + 2.*7.6438e-5*t(:)                         &
       - 3.*8.2467e-7*t2(:) + 4.*5.3875e-9*t3(:))                       &
       + s32(:)*(1.0227e-4 - 2.*1.6546e-6*t(:))
drho0ds= (0.824493 - 4.0899e-3*t(:) + 7.6438e-5*t2(:)                   &
       - 8.2467e-7*t3(:) + 5.3875e-9*t4(:))                             &
       + 1.5*sqrt(s(:))*(-5.72466e-3 + 1.0227e-4*t(:)                   &
       - 1.6546e-6*t2(:)) + 2.*4.8314e-4*s(:)

!do i=1,nits
  ptot=pxtr*1.E-5
  do ii=1,wlx
    t = min(max(tt(:,ii)-273.16,-2.2),100.)
    s = min(max(ss(:,ii),0.),maxsal)
    p1   = ptot+grav*d_rho(:,ii)*0.5*ddz(:,ii)*1.E-5 ! hydrostatic approximation
    ptot = ptot+grav*d_rho(:,ii)*ddz(:,ii)*1.E-5
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    s2 = s*s
    s3 = s2*s
    p2 = p1*p1
    s32 = sqrt(s3)
    
    sks = 1.965933e4 + 1.444304e2*t(:) - 1.706103*t2(:)                 &
                + 9.648704e-3*t3(:)  - 4.190253e-5*t4(:)                &
                + p1(:)*(3.186519 + 2.212276e-2*t(:)                    &
                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:))                &
                + p2(:)*(2.102898e-4 - 1.202016e-5*t(:)                 &
                + 1.394680e-7*t2(:)) ! sal=0.
    sk=sks+ s(:)*(52.84855 - 3.101089e-1*t(:)                           &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:))                 &
                + s32(:)*(3.886640e-1 + 9.085835e-3*t(:)                &
                - 4.619924e-4*t2(:))                                    &
                + p1(:)*s(:)*(6.704388e-3  -1.847318e-4*t(:)            &
                + 2.059331e-7*t2(:)) + 1.480266e-4*p1(:)*s32(:)         &
                +p2(:)*s(:)*(-2.040237e-6 &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:)) ! + sal terms             
    dskdt= 1.444304e2 - 2.*1.706103*t(:)                                &
                + 3.*9.648704e-3*t2(:)  - 4.*4.190253e-5*t3(:)          &
                + s(:)*( - 3.101089e-1                                  &
                + 2.*6.283263e-3*t(:) -3.*5.084188e-5*t2(:))            &
                + s32(:)*(9.085835e-3                                   &
                - 2.*4.619924e-4*t(:))                                  &
                + p1(:)*(2.212276e-2                                    &
                - 2.*2.984642e-4*t(:) + 3.*1.956415e-6*t2(:))           &
                + p1(:)*s(:)*(-1.847318e-4                              &
                + 2.*2.059331e-7*t(:))                                  &
                + p2(:)*(- 1.202016e-5                                  &
                + 2.*1.394680e-7*t(:)) +p2(:)*s(:)*(                    &
                + 6.128773e-8 + 2.*6.207323e-10*t(:))
    dskds=(52.84855 - 3.101089e-1*t(:)                                  &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:))                 &
                + 1.5*sqrt(s(:))*(3.886640e-1 + 9.085835e-3*t(:)        &
                - 4.619924e-4*t2(:))                                    &
                + p1(:)*(6.704388e-3  -1.847318e-4*t(:)                 &
                + 2.059331e-7*t2(:)) + 1.5*1.480266e-4*p1(:)*sqrt(s(:)) &
                +p2(:)*(-2.040237e-6                                    &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
!    dskdp=(3.186519 + 2.212276e-2*t(:)                                 &
!                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:))               &
!                + 2.*p1*(2.102898e-4 - 1.202016e-5*t(:)                &
!                + 1.394680e-7*t2(:))                                   &
!                + s(:)*(6.704388e-3  -1.847318e-4*t(:)                 &
!                + 2.059331e-7*t2(:))                                   &
!                + 1.480266e-4*s32(:)                                   &
!                + 2.*p1(:)*s(:)*(-2.040237e-6                          &
!                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
       
    d_rho(:,ii)=rho0/(1.-p1/sk)
  
    drhodt=drho0dt/(1.-p1/sk)-rho0*p1*dskdt/((sk-p1)**2) ! neglected dp1drho*drhodt terms
    drhods=drho0ds/(1.-p1/sk)-rho0*p1*dskds/((sk-p1)**2) ! neglected dp1drho*drhods terms
    
    d_alpha(:,ii)=-drhodt              ! Large et al (1993) convention
    d_beta(:,ii)=drhods                ! Large et al (1993) convention
    
  end do

!end do

return
end subroutine calcdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate melting/freezing point
subroutine calcmelt(d_timelt,d_zcr)

implicit none

integer iqw,ii
real, dimension(wfull), intent(out) :: d_timelt
real, dimension(wfull), intent(in) :: d_zcr
real, dimension(wfull) :: ssf
real dsf,aa,bb,deldz

! Integrate over minsfc to estimate near surface salinity
ssf=0.
do iqw=1,wfull
  dsf=0.
  do ii=1,wlev-1
    aa=dz(iqw,ii)*d_zcr(iqw)
    bb=minsfc-dsf
    deldz=min(aa,bb)
    ssf(iqw)=ssf(iqw)+water%sal(iqw,ii)*deldz
    dsf=dsf+deldz
    if (aa>=bb) exit
  end do
  ssf(iqw)=ssf(iqw)/dsf
end do
d_timelt=273.16-0.054*min(max(ssf,0.),maxsal) ! ice melting temperature from CICE

return
end subroutine calcmelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert temp to theta based on MOM3 routine
!subroutine gettheta(theta,tt,ss,pxtr,ddz)
!
!implicit none
!
!integer ii
!real, dimension(wfull,wlev), intent(out) :: theta
!real, dimension(wfull,wlev), intent(in) :: tt,ss,ddz
!real, dimension(wfull), intent(in) :: pxtr
!real, dimension(wfull) :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11
!real, dimension(wfull) :: t,t2,t3,s,s2,p,p2,potmp,ptot
!real, dimension(wfull,wlev) :: d_rho
!real, parameter :: density = 1035.
!
!d_rho=density
!
!ptot=pxtr*1.E-5
!do ii=1,wlev
!  t = max(tt(:,ii)-273.16,-2.)
!  s = max(ss(:,ii),0.)
!  p   = ptot+grav*d_rho(:,ii)*0.5*ddz(:,ii)*1.E-5 ! hydrostatic approximation
!  ptot = ptot+grav*d_rho(:,ii)*ddz(:,ii)*1.E-5
!    
!  b1    = -1.60e-5*p
!  b2    = 1.014e-5*p*t
!  t2    = t*t
!  t3    = t2*t
!  b3    = -1.27e-7*p*t2
!  b4    = 2.7e-9*p*t3
!  b5    = 1.322e-6*p*s
!  b6    = -2.62e-8*p*s*t
!  s2    = s*s
!  p2    = p*p
!  b7    = 4.1e-9*p*s2
!  b8    = 9.14e-9*p2
!  b9    = -2.77e-10*p2*t
!  b10   = 9.5e-13*p2*t2
!  b11   = -1.557e-13*p2*p
!  potmp = b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11
!  theta(:,ii) = t-potmp+273.16
!end do
!
!return
!end subroutine gettheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes between MLO and atmosphere

subroutine fluxcalc(dt,atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins, &
                    d_neta,atm_oldu,atm_oldv,diag)

implicit none

integer, intent(in) :: diag
integer it
real, intent(in) :: dt
real, dimension(wfull), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
real, dimension(wfull), intent(inout) :: d_neta,atm_oldu,atm_oldv
real, dimension(wfull) :: qsat,dqdt,ri,rho,srcp
real, dimension(wfull) :: fm,fh,fq,con,consea,afroot,af,daf
real, dimension(wfull) :: den,dfm,dden,dcon,sig,factch,root
real, dimension(wfull) :: aft,afq,atu,atv,dcs,facqch
real, dimension(wfull) :: vmagn,egmax,d_wavail,dumwatertemp
real ztv
! momentum flux parameters
real, parameter :: charnck=0.018
real, parameter :: chn10=0.00125
! ..... high wind speed - rough sea
real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
real, parameter :: zcoq1 = 0.0
! ..... low wind speed - smooth sea
real, parameter :: gnu   = 1.5e-5
real, parameter :: zcom2 = 0.11
real, parameter :: zcoh2 = 0.40
real, parameter :: zcoq2 = 0.62

dumwatertemp=max(water%temp(:,1)+wrtemp,271.)
sig=exp(-grav*atm_zmins/(rdry*atm_temp))
srcp=sig**(rdry/cpair)
atu=atm_u-fluxwgt*water%u(:,1)-(1.-fluxwgt)*atm_oldu
atv=atm_v-fluxwgt*water%v(:,1)-(1.-fluxwgt)*atm_oldv
vmagn=sqrt(max(atu*atu+atv*atv,1.e-4))
rho=atm_ps/(rdry*dumwatertemp)
ri=min(grav*(atm_zmin*atm_zmin/atm_zmins)*(1.-dumwatertemp*srcp/atm_temp)/vmagn**2,rimax)

call getqsat(qsat,dqdt,dumwatertemp,atm_ps)
if (zomode==0) then ! CSIRO9
  qsat=0.98*qsat ! with Zeng 1998 for sea water
end if

select case(zomode)
  case(0,1) ! Charnock
    consea=vmagn*charnck/grav
    dgwater%zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(atm_zmin/dgwater%zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*dgwater%zo)
      where (ri>=0.) ! stable water points                                     
        fm=1./(1.+bprm*ri)**2
        con=consea*fm*vmagn
        dgwater%zo=dgwater%zo-(dgwater%zo-con*af)/(1.-con*daf)
      elsewhere     ! unstable water points
        den=1.+af*cms*2.*bprm*sqrt(-ri*atm_zmin/dgwater%zo)
        fm=1.-2.*bprm*ri/den
        con=consea*fm*vmagn
        dden=daf*cms*2.*bprm*sqrt(-ri*atm_zmin/dgwater%zo)+af*cms*bprm*sqrt(-ri)*atm_zmin &
            /(sqrt(atm_zmin/dgwater%zo)*dgwater%zo*dgwater%zo)
        dfm=2.*bprm*ri*dden/(den*den)
        dcon=consea*dfm*vmagn
        dgwater%zo=dgwater%zo-(dgwater%zo-con*af)/(1.-dcon*af-con*daf)
      end where
      dgwater%zo=min(max(dgwater%zo,1.5e-5),6.)
    enddo    ! it=1,4
  case(2) ! Beljaars
    dgwater%zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(atm_zmin/dgwater%zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*dgwater%zo)
      where (ri>=0.) ! stable water points
        fm=1./(1.+bprm*ri)**2
        consea=zcom1*vmagn*vmagn*af*fm/grav+zcom2*gnu/(vmagn*sqrt(fm*af))
        dcs=(zcom1*vmagn*vmagn/grav-0.5*zcom2*gnu/(vmagn*sqrt(fm*af)*fm*af))*(fm*daf)
      elsewhere     ! unstable water points
        con=cms*2.*bprm*sqrt(-ri*atm_zmin/dgwater%zo)
        den=1.+af*con
        fm=1.-2.*bprm*ri/den
        dfm=2.*bprm*ri*(con*daf+af*cms*bprm*sqrt(-ri)*atm_zmin/(sqrt(atm_zmin/dgwater%zo)*dgwater%zo*dgwater%zo))/(den*den)
        consea=zcom1*vmagn*vmagn*af*fm/grav+zcom2*gnu/(vmagn*sqrt(fm*af))
        dcs=(zcom1*vmagn*vmagn/grav-0.5*zcom2*gnu/(vmagn*sqrt(fm*af)*fm*af))*(fm*daf+dfm*af)
      end where
      dgwater%zo=dgwater%zo-(dgwater%zo-consea)/(1.-dcs)      
      dgwater%zo=min(max(dgwater%zo,1.5e-5),6.)
    enddo    ! it=1,4
end select
afroot=vkar/log(atm_zmin/dgwater%zo)
af=afroot*afroot

select case(zomode)
  case(0) ! Charnock CSIRO9
    ztv=exp(vkar/sqrt(chn10))/10.
    aft=(vkar/log(atm_zmins*ztv))**2
    afq=aft
    dgwater%zoh=1./ztv
    dgwater%zoq=dgwater%zoh
    factch=sqrt(dgwater%zo/dgwater%zoh)
    facqch=factch
  case(1) ! Charnock zot=zom
    dgwater%zoh=dgwater%zo
    dgwater%zoq=dgwater%zo
    aft=vkar*vkar/(log(atm_zmins/dgwater%zo)*log(atm_zmins/dgwater%zo))
    afq=aft
    factch=1.
    facqch=1.
  case(2) ! Beljaars
    dgwater%zoh=max(zcoh1+zcoh2*gnu/(vmagn*sqrt(fm*af)),1.5E-7)
    dgwater%zoq=max(zcoq1+zcoq2*gnu/(vmagn*sqrt(fm*af)),1.5E-7)
    factch=sqrt(dgwater%zo/dgwater%zoh)
    facqch=sqrt(dgwater%zo/dgwater%zoq)
    aft=vkar*vkar/(log(atm_zmins/dgwater%zo)*log(atm_zmins/dgwater%zoh))
    afq=vkar*vkar/(log(atm_zmins/dgwater%zo)*log(atm_zmins/dgwater%zoq))
end select

where (ri>=0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
  fq=fm
elsewhere        ! ri is -ve
  root=sqrt(-ri*atm_zmin/dgwater%zo)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  root=sqrt(-ri*atm_zmins/dgwater%zo)
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*facqch*afq*root
  fq=1.-2.*bprm*ri/den
end where

! update drag
dgwater%cd=af*fm
dgwater%cdh=aft*fh
dgwater%cdq=afq*fq

! turn off lake evaporation when minimum depth is reached
! fg should be replaced with bare ground value
d_wavail=max(depth_hl(:,wlev+1)+d_neta-minwater,0.)
egmax=1000.*lv*d_wavail/(dt*max(1.-ice%fracice,0.01))

! explicit estimate of fluxes
! (replace with implicit scheme if water becomes too shallow)
dgwater%fg=rho*dgwater%cdh*cpair*vmagn*(dumwatertemp-atm_temp/srcp)
dgwater%fg=min(max(dgwater%fg,-3000.),3000.)
dgwater%eg=min(rho*dgwater%cdq*lv*vmagn*(qsat-atm_qg),egmax)
dgwater%eg=min(max(dgwater%eg,-3000.),3000.)
dgwater%taux=rho*dgwater%cd*vmagn*atu
dgwater%tauy=rho*dgwater%cd*vmagn*atv

! update free surface after evaporation
d_neta=d_neta-0.001*dt*(1.-ice%fracice)*dgwater%eg/lv

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from CCAM)
! following version of getqsat copes better with T < 0C over ice
subroutine getqsat(qsat,dqdt,temp,ps)

implicit none

real, dimension(wfull), intent(in) :: temp,ps
real, dimension(wfull), intent(out) :: qsat,dqdt
real, dimension(0:220), save :: table
real, dimension(wfull) :: esatf,tdiff,dedt,rx
integer, dimension(wfull) :: ix
logical, save :: first=.true.

if (first) then
  table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
  table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
  table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
  table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
  table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
  table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
  table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
  table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
  table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
  table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
  table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
  table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
  table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
  table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
  table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
  table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
  table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
  table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
  table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
  table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
  table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
  table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
  table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
  table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
  table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
  table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
  table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
  table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
  table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
  table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
  table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
  table(219:220)=(/ 29845.0, 31169.0 /)
  first=.false.
end if

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

! method #1
dedt=table(ix+1)-table(ix) ! divide by 1C
dqdt=qsat*dedt*ps/(esatf*max(ps-esatf,0.1))
! method #2
!dqdt=qsat*lv/rvap*qsat*ps/(temp*temp*(ps-0.378*esatf))

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Pack sea ice for calcuation

subroutine mloice(dt,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt, &
                  d_ustar,d_nk,d_neta,d_imass,diag)

implicit none

integer, intent(in) :: diag
integer ii
integer, dimension(wfull), intent(inout) :: d_nk
real, intent(in) :: dt
real, dimension(wfull,wlev), intent(in) :: d_alpha, d_beta
real, dimension(wfull), intent(inout) :: d_b0, d_wu0, d_wv0, d_wt0, d_ws0, d_ftop, d_tb, d_fb, d_timelt
real, dimension(wfull), intent(inout) :: d_ustar, d_neta, d_imass
real, dimension(wfull) :: d_salflxf, d_salflxs, d_wavail, d_avewtemp
real, dimension(wfull) :: deld, xxx, newthick

d_salflxf=0.                                        ! fresh water flux
d_salflxs=0.                                        ! salt water flux
d_wavail=max(depth_hl(:,wlev+1)+d_neta-minwater,0.) ! water avaliable for freezing

! Calculate average temperature of water column for energy conservation
d_avewtemp=0.
do ii=1,wlev
  d_avewtemp=d_avewtemp+water%temp(:,ii)*dz(:,ii)
end do
d_avewtemp=d_avewtemp/depth_hl(:,wlev+1)
d_avewtemp=d_avewtemp+wrtemp

! update ice prognostic variables
call seaicecalc(dt,d_ftop,d_tb,d_fb,d_timelt,d_salflxf,d_salflxs,d_nk,d_wavail,d_avewtemp,diag)

! Ice depth limitation for poor initial conditions
xxx=max(ice%thick-icemax,0.)
where ( xxx>0.001 )
  newthick=ice%thick-xxx    
  ice%temp(:,1)=(ice%temp(:,1)*cpi*ice%thick-cp0*rhoic*d_avewtemp*xxx)/(cpi*newthick)
  ice%temp(:,2)=(ice%temp(:,2)*cpi*ice%thick-cp0*rhoic*d_avewtemp*xxx)/(cpi*newthick)
  ice%thick=newthick
  d_salflxs=d_salflxs-rhoic*xxx/dt ! saltwater leaving ocean to ice
end where

! update ice velocities due to stress terms
ice%u=ice%u+dt*(dgice%tauxica-dgice%tauxicw)/d_imass
ice%v=ice%v+dt*(dgice%tauyica-dgice%tauyicw)/d_imass

! Remove excessive salinity from ice
where ( ice%thick>=icemin )
  deld=ice%thick*(1.-maxicesal/max(min(ice%sal,2.*maxicesal),maxicesal))
  ice%sal=ice%sal*(1.-deld/ice%thick)
elsewhere
  deld=0.
end where
! no change in ice thickness as fresh water and salt water fluxes are in balance
d_salflxf=d_salflxf+deld*rhoic/dt ! freshwater leaving ocean to ice
d_salflxs=d_salflxs-deld*rhoic/dt ! saltwater leaving ocean to ice

! update water boundary conditions
! MJT notes - use rhowt reference density for Boussinesq fluid approximation
d_wu0=d_wu0-ice%fracice*dgice%tauxicw/rhowt
d_wv0=d_wv0-ice%fracice*dgice%tauyicw/rhowt
d_wt0=d_wt0+ice%fracice*d_fb/(rhowt*cp0)
d_ws0=d_ws0-ice%fracice*(d_salflxf*water%sal(:,1)/rhowt+d_salflxs*(water%sal(:,1)-ice%sal)/rhowt)
d_ustar=sqrt(sqrt(max(d_wu0*d_wu0+d_wv0*d_wv0,1.E-24)))
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign is to account for sign of wt0 and ws0

! update free surface height with water flux from ice
d_neta=d_neta-dt*ice%fracice*(d_salflxf+d_salflxs)/rhowt

return
end subroutine mloice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice(d_timelt,d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer iqw, ii, maxlevel
real aa, bb, dsf, deldz, delt, dels, neutralthick
real, dimension(wfull,wlev) :: sdic
real, dimension(wfull) :: newdic, newicesal, newtn, cdic
real, dimension(wfull), intent(inout) :: d_timelt, d_zcr
real, dimension(wfull) :: maxnewice, d_wavail
real, dimension(wfull) :: newthick
real, dimension(wfull) :: avesal, avetemp, newicetemp
logical, dimension(wfull) :: lnewice, lremove

! limits on ice formation
d_wavail=max(depth_hl(:,wlev+1)+water%eta-minwater,0.)
where ( ice%fracice<0.999 )
  maxnewice=d_wavail*rhowt/rhoic/(1.-ice%fracice)
elsewhere
  maxnewice=0.
end where

! search for water temperatures that are below freezing
sdic=0.
maxlevel=1
do iqw=1,wfull
  dsf=0.
  if ( ice%fracice(iqw)<0.999 ) then
    do ii=1,wlev-1
      aa=dz(iqw,ii)*d_zcr(iqw)
      bb=max(minsfc-dsf,0.)
      deldz=min(aa,bb)
      sdic(iqw,ii)=max(d_timelt(iqw)-water%temp(iqw,ii)-wrtemp,0.)*cp0*rhowt*deldz/qice/(1.-ice%fracice(iqw))
      dsf=dsf+deldz
      maxlevel=max(maxlevel,ii)
      if (bb<=0.) exit
    end do
  end if
end do
newdic=sum(sdic,2)
newdic=min(newdic,maxnewice)
neutralthick=gammi/(cp0*rhoic) ! should be around 0.02 m or 2*icemin
newdic=min(newdic,1.01*neutralthick)
lnewice=newdic>0.99*neutralthick
where ( .not.lnewice )
  newdic=0.
end where

! calculate temperature of water to be converted to ice
avetemp=0.
avesal=0.
do ii=1,wlev
  avetemp=avetemp+water%temp(:,ii)*dz(:,ii)
  avesal=avesal+water%sal(:,ii)*dz(:,ii)
end do
avetemp=avetemp/depth_hl(:,wlev+1)
avetemp=avetemp+wrtemp
avesal=avesal/depth_hl(:,wlev+1)
! balance the energy leaving the water column with that in surface ice layer
newicetemp=avetemp*cp0*rhoic*newdic/gammi
newicesal=avesal

water%eta=water%eta-newdic*(1.-ice%fracice)*rhoic/rhowt
d_zcr=max(1.+water%eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))

! Adjust temperature in water column to balance ice formation
cdic=0.
do ii=1,maxlevel
  sdic(:,ii)=max(min(sdic(:,ii),newdic-cdic),0.)
  cdic=cdic+sdic(:,ii)  
  where ( lnewice )
    newtn=water%temp(:,ii)+qice*sdic(:,ii)/(cp0*rhowt*dz(:,ii)*d_zcr)
    water%temp(:,ii)=water%temp(:,ii)*ice%fracice+newtn*(1.-ice%fracice)
  end where
end do

! form new sea-ice
do iqw=1,wfull
  if ( lnewice(iqw) ) then
    newthick(iqw) =ice%thick(iqw)*ice%fracice(iqw)+newdic(iqw)*(1.-ice%fracice(iqw))
    ice%tsurf(iqw)=ice%tsurf(iqw)*ice%fracice(iqw)+newicetemp(iqw)*(1.-ice%fracice(iqw))
    !ice%temp(iqw,0)*newsnowd(iqw)=ice%temp(iqw,0)*ice%fracice(iqw)*ice%snowd(iqw)   ! unchanged
    if ( newthick(iqw)>himin ) then
      ! ice thickness is decreasing after combining with new thin ice
      ice%temp(iqw,1)=(ice%temp(iqw,1)*ice%thick(iqw)*ice%fracice(iqw)+newicetemp(iqw)*newdic(iqw)*(1.-ice%fracice(iqw))) &
                     /newthick(iqw)
      ice%temp(iqw,2)=(ice%temp(iqw,2)*ice%thick(iqw)*ice%fracice(iqw)+newicetemp(iqw)*newdic(iqw)*(1.-ice%fracice(iqw))) &
                     /newthick(iqw)
    else
      ice%tsurf(iqw)=(ice%tsurf(iqw)*gammi+0.5*(ice%temp(iqw,1)+ice%temp(iqw,2))*cpi*ice%thick(iqw)*ice%fracice(iqw)      &
                     +newicetemp(iqw)*cpi*newdic(iqw)*(1.-ice%fracice(iqw)))/(gammi+cpi*newthick(iqw))
      ice%temp(iqw,1)=ice%tsurf(iqw)
      ice%temp(iqw,2)=ice%tsurf(iqw)
    end if
    ice%store(iqw)=ice%store(iqw)*ice%fracice(iqw)
    ice%sal(iqw)=(ice%sal(iqw)*ice%fracice(iqw)*ice%thick(iqw)+newicesal(iqw)*(1.-ice%fracice(iqw))*newdic(iqw))          &
                 /newthick(iqw)
    ice%thick(iqw)=newthick(iqw)
    ice%snowd(iqw)=ice%snowd(iqw)*ice%fracice(iqw)
    ice%fracice(iqw)=1. ! this will be adjusted for fracbreak below  
  end if
end do

! 1D model of ice break-up
! MJT notes - ice%thick=icemin implies zero energy stored in t1 and t2
!newthick=1.01*icemin
!newfracice=ice%fracice*ice%thick/newthick
!where ( ice%thick<newthick .and. newfracice>fracbreak )
!  ice%tsurf    =ice%tsurf*ice%fracice/newfracice
!  !ice%temp(:,0)*newsnowd=ice%temp(:,0)*ice%fracice*ice%snowd
!  ice%temp(:,1)=ice%tsurf
!  ice%temp(:,2)=ice%tsurf
!  ice%snowd    =ice%snowd*ice%fracice/newfracice
!  ice%thick    =newthick
!  ice%fracice  =newfracice
!end where
!if ( onedice==1 ) then
!  where ( ice%fracice<1. .and. ice%thick>1.01*icemin )
!     worka=min(ice%thick/(1.01*icemin),1./max(ice%fracice,fracbreak))
!     worka=max(worka,1.)
!     newfracice   =ice%fracice*worka
!     newthick     =ice%thick/worka
!     ice%tsurf    =ice%tsurf*ice%fracice/newfracice
!     !ice%temp(:,0)*newsnowd=ice%temp(:,0)*ice%fracice*ice%snowd
!     ! assume that ice%thick is < 0.5*himin so that the following terms are zero
!     !ice%temp(:,1)*newthick=ice%temp(:,1)*ice%thick
!     !ice%temp(:,2)*newthick=ice%temp(:,2)*ice%thick
!     ice%fracice  =newfracice
!     ice%thick    =newthick
!     ice%snowd    =ice%snowd*ice%fracice/newfracice
!  end where
!end if

! calculate temperature of water to be converted to ice
avetemp=0.
avesal=0.
do ii=1,wlev
  avetemp=avetemp+water%temp(:,ii)*dz(:,ii)
  avesal=avesal+water%sal(:,ii)*dz(:,ii)
end do
avetemp=avetemp/depth_hl(:,wlev+1)
avetemp=avetemp+wrtemp
avesal=avesal/depth_hl(:,wlev+1)

! removal
lremove = ice%thick<=icemin .and. ice%fracice>0.
where ( lremove )
  newdic=(ice%thick*rhoic+ice%snowd*rhosn)*ice%fracice/rhowt
elsewhere
  newdic=0.
end where
water%eta=water%eta+newdic
d_zcr=max(1.+water%eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))

do iqw=1,wfull
  if ( lremove(iqw) ) then

    ! update average temperature and salinity
    dsf=minsfc
    dels=ice%sal(iqw)*ice%thick(iqw)*rhoic/rhowt*ice%fracice(iqw)/dsf
    dels=dels-avesal(iqw)*newdic(iqw)/dsf
    delt=ice%fracice(iqw)*gammi*ice%tsurf(iqw)/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*cps*ice%snowd(iqw)*ice%temp(iqw,0)/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*0.5*cpi*ice%thick(iqw)*ice%temp(iqw,1)/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*0.5*cpi*ice%thick(iqw)*ice%temp(iqw,2)/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*ice%thick(iqw)*qice/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*ice%snowd(iqw)*qsnow/(cp0*rhowt*dsf)
    delt=delt+ice%fracice(iqw)*ice%store(iqw)/(cp0*rhowt*dsf)
    delt=delt-cp0*rhowt*avetemp(iqw)*newdic(iqw)/(cp0*rhowt*dsf)
    
    ! adjust temperature and salinity in water column
    dsf=0.
    do ii=1,wlev
      aa=dz(iqw,ii)*d_zcr(iqw)
      bb=max(minsfc-dsf,0.)
      deldz=min(aa,bb)
      water%temp(iqw,ii)=water%temp(iqw,ii)+delt*deldz/aa
      water%sal(iqw,ii) =water%sal(iqw,ii) +dels*deldz/aa
      water%sal(iqw,ii)=max(water%sal(iqw,ii),0.)
      dsf=dsf+deldz
      if ( bb<=0. ) exit
    end do

    ! remove ice
    ice%fracice(iqw)=0.
    ice%thick(iqw)=0.
    ice%snowd(iqw)=0.
    ice%store(iqw)=0.
    ice%sal(iqw)=0.
    ice%tsurf(iqw)=273.16
    ice%temp(iqw,:)=273.16
  end if
end do

return
end subroutine mlonewice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(dt,d_ftop,d_tb,d_fb,d_timelt,d_salflxf,d_salflxs,d_nk,    &
                      d_wavail,d_avewtemp,diag)

implicit none

integer, intent(in) :: diag
integer pc
integer, dimension(5) :: nc
integer, dimension(wfull), intent(inout) :: d_nk
integer, dimension(wfull) :: dt_nk
real, intent(in) :: dt
real, dimension(wfull), intent(inout) :: d_ftop,d_tb,d_fb,d_timelt,d_salflxf,d_salflxs
real, dimension(wfull), intent(inout) :: d_wavail,d_avewtemp
real, dimension(wfull) :: it_tn0,it_tn1,it_tn2
real, dimension(wfull) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(wfull) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(wfull) :: dt_avewtemp
real, dimension(wfull) :: pt_egice
logical, dimension(wfull,5) :: pqpack

! Pack different ice configurations
pqpack(:,1)=( d_nk==2 .and. ice%snowd>0.05 .and. ice%thick>icemin )   ! thick snow + 2 ice layers
nc(1)=count(pqpack(:,1))
pqpack(:,2)=( d_nk==1 .and. ice%snowd>0.05 .and. ice%thick>icemin )   ! thick snow + 1 ice layer
nc(2)=count(pqpack(:,2))
pqpack(:,3)=( d_nk==2 .and. ice%snowd<=0.05 .and. ice%thick>icemin )  ! 2 ice layers (+ snow?)
nc(3)=count(pqpack(:,3))
pqpack(:,4)=( d_nk==1 .and. ice%snowd<=0.05 .and. ice%thick>icemin )  ! 1 ice layer (+ snow?)
nc(4)=count(pqpack(:,4))
pqpack(:,5)=( d_nk==0 .and. ice%thick>icemin )                        ! thin ice (+ snow?)
nc(5)=count(pqpack(:,5))

! Update ice and snow temperatures
do pc=1,5
  if ( nc(pc)>0 ) then
    it_tsurf(1:nc(pc))  =pack(ice%tsurf,pqpack(:,pc))
    it_tn0(1:nc(pc))    =pack(ice%temp(:,0),pqpack(:,pc))
    it_tn1(1:nc(pc))    =pack(ice%temp(:,1),pqpack(:,pc))
    it_tn2(1:nc(pc))    =pack(ice%temp(:,2),pqpack(:,pc))
    pt_egice(1:nc(pc))  =pack(dgice%eg,pqpack(:,pc))
    it_dic(1:nc(pc))    =pack(ice%thick,pqpack(:,pc))
    it_dsn(1:nc(pc))    =pack(ice%snowd,pqpack(:,pc))
    it_sto(1:nc(pc))    =pack(ice%store,pqpack(:,pc))
    dt_ftop(1:nc(pc))   =pack(d_ftop,pqpack(:,pc))
    dt_tb(1:nc(pc))     =pack(d_tb,pqpack(:,pc))
    dt_fb(1:nc(pc))     =pack(d_fb,pqpack(:,pc))
    dt_timelt(1:nc(pc)) =pack(d_timelt,pqpack(:,pc))
    dt_salflxf(1:nc(pc))=pack(d_salflxf,pqpack(:,pc))
    dt_salflxs(1:nc(pc))=pack(d_salflxs,pqpack(:,pc))
    dt_wavail(1:nc(pc)) =pack(d_wavail,pqpack(:,pc))
    dt_nk(1:nc(pc))     =pack(d_nk,pqpack(:,pc))
    dt_avewtemp(1:nc(pc))=pack(d_avewtemp,pqpack(:,pc))
    select case(pc)
      case(1)
        call icetemps1s2i(nc(pc),dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,            &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,         &
                 dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)
      case(2)
        call icetemps1s1i(nc(pc),dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,            &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,         &
                 dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)
      case(3)
        call icetempi2i(nc(pc),dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,              &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,         &
                 dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)
      case(4)
        call icetempi1i(nc(pc),dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,              &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,         &
                 dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)
      case(5)
        call icetemps(nc(pc),dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,                &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,         &
                 dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)
    end select
    ice%tsurf    =unpack(it_tsurf(1:nc(pc)),pqpack(:,pc),ice%tsurf)
    ice%temp(:,0)=unpack(it_tn0(1:nc(pc)),pqpack(:,pc),ice%temp(:,0))
    ice%temp(:,1)=unpack(it_tn1(1:nc(pc)),pqpack(:,pc),ice%temp(:,1))
    ice%temp(:,2)=unpack(it_tn2(1:nc(pc)),pqpack(:,pc),ice%temp(:,2))
    ice%thick    =unpack(it_dic(1:nc(pc)),pqpack(:,pc),ice%thick)
    ice%snowd    =unpack(it_dsn(1:nc(pc)),pqpack(:,pc),ice%snowd)
    ice%store    =unpack(it_sto(1:nc(pc)),pqpack(:,pc),ice%store)
    d_salflxf    =unpack(dt_salflxf(1:nc(pc)),pqpack(:,pc),d_salflxf)
    d_salflxs    =unpack(dt_salflxs(1:nc(pc)),pqpack(:,pc),d_salflxs)
    d_nk         =unpack(dt_nk(1:nc(pc)),pqpack(:,pc),d_nk)
  end if
end do

return
end subroutine seaicecalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and snow (from Mk3 seaice.f)

! Solve for snow surface temperature (implicit time scheme)
!
! Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ts = Snow temp for thin layer near surface................. S
! -------------------------------------------------------------- N
!                                                         fs(up) O
! ...t0 = Mid point temp of snow layer.......................... W
!
! ## ti = Temp at snow/ice interface #############################
!                                                         f0(up) I
! ...t(1) = Mid point temp of first ice layer................... C
!                                                                E
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetemps1s2i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real htdown
real, dimension(nc) :: con,rhin,rhsn,conb,qmax,fl,qneed,xxx,excess,conc
real, dimension(nc) :: subl,snmelt,dhb,isubl,ssubl,smax
real, dimension(nc) :: simelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nc), intent(inout) :: dt_avewtemp
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice
real(kind=8), dimension(nc,2:4) :: aa
real(kind=8), dimension(nc,4) :: bb,dd
real(kind=8), dimension(nc,3) :: cc
real, dimension(nc,4) :: ans

! Thickness of each layer
rhsn=1./it_dsn
rhin=2./(it_dic-gammi/cpi)
con =condsnw/it_dsn
conb=1./(it_dsn/condsnw+0.5*max(it_dic,icemin)/condice)
conc=2.*condice/max(it_dic,icemin)

! no need to map from generic ice pack into 2 layer + snow

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gammi
cc(:,1)=-dt*2.*con/gammi
dd(:,1)=it_tsurf+dt*(dt_ftop-pt_egice*lf/lv)/gammi
aa(:,2)=-dt*2.*con*rhsn/cps
bb(:,2)=1.+dt*2.*con*rhsn/cps+dt*2.*conb*rhsn/cps
cc(:,2)=-dt*2.*conb*rhsn/cps
dd(:,2)=it_tn0
aa(:,3)=-dt*2.*conb*rhin/cpi
bb(:,3)=1.+dt*2.*conb*rhin/cpi+dt*conc*rhin/cpi
cc(:,3)=-dt*conc*rhin/cpi
dd(:,3)=it_tn1
aa(:,4)=-dt*conc*rhin/cpi
bb(:,4)=1.+dt*conc*rhin/cpi+dt*2.*conc*rhin/cpi
dd(:,4)=it_tn2+dt*2.*conc*dt_tb*rhin/cpi
call thomas(ans,aa,bb,cc,dd)
it_tsurf(:)=ans(:,1)
it_tn0(:)  =ans(:,2)
it_tn1(:)  =ans(:,3)
it_tn2(:)  =ans(:,4)
fl=2.*conc*(dt_tb-it_tn2)
dhb=dt*(fl-dt_fb)/qice                      ! Excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+(flnew-fl)/(cpi*it_dic) ! Modify temperature if limit is reached
it_tn2=it_tn2+(flnew-fl)/(cpi*it_dic) ! Modify temperature if limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflxs=dt_salflxs+dhb*rhoic/dt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+0.5*(it_tn1+it_tn2)*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*lf/(lv*qsnow)        ! net snow+ice sublimation
ssubl=min(subl,it_dsn)                ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic        ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhosn/dt
dt_salflxs=dt_salflxs-isubl*rhosn/dt
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=it_dic-smax
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+0.5*(it_tn1+it_tn2)*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! use stored heat in brine pockets to keep temperature at -0.1 until heat is used up
qneed=(dt_timelt-it_tn1)*0.5*cpi*it_dic ! J/m**2
qneed=max(min(qneed,it_sto),0.)
where ( it_dic>himin )
  it_tn1=it_tn1+2.*qneed/(cpi*it_dic)
  it_sto=it_sto-qneed
end where

! the following are snow to ice processes
xxx=it_dic+it_dsn-(rhosn*it_dsn+rhoic*it_dic)/rhowt ! white ice formation
excess=max(it_dsn-xxx,0.)*rhosn/rhowt               ! white ice formation
excess=excess+max(it_dsn-excess*rhowt/rhosn-0.2,0.)*rhosn/rhowt  ! Snow depth limitation and conversion to ice
it_dsn=it_dsn-excess*rhowt/rhosn
it_dic=it_dic+excess*rhowt/rhoic
dt_salflxf=dt_salflxf-excess*rhowt/dt
dt_salflxs=dt_salflxs+excess*rhowt/dt

! Snow melt
snmelt=max(it_tsurf-273.16,0.)*gammi/qsnow
snmelt=min(snmelt,it_dsn)
it_tsurf=it_tsurf-snmelt*qsnow/gammi
dt_salflxf=dt_salflxf-snmelt*rhosn/dt        ! melt fresh water snow (no salt when melting snow)
it_dsn=it_dsn-snmelt
where ( it_dsn>0.05 )
  it_tn0=(it_tn0*cps*(it_dsn+snmelt)-cp0*rhosn*dt_avewtemp*snmelt)/(cps*it_dsn)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn0*cps*(it_dsn+snmelt)-cp0*rhosn*dt_avewtemp*snmelt) &
          /(gammi+cps*it_dsn)
  it_tn0=it_tsurf
end where

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gammi/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gammi
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=it_dic-simelt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+0.5*(it_tn1+it_tn2)*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! test whether to change number of layers
htdown=2.*himin
where ( it_dic<=htdown )
  ! code to decrease number of layers
  dt_nk=1
  ! merge ice layers into one
  it_tn1=0.5*(it_tn1+it_tn2)
  it_tn2=it_tn1
end where

return
end subroutine icetemps1s2i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and snow (from Mk3 seaice.f)

! Solve for snow surface temperature (implicit time scheme)
!
! Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ts = Snow temp for thin layer near surface................. S
! -------------------------------------------------------------- N
!                                                         fs(up) O
! ...t0 = Mid point temp of snow layer.......................... W
!
! ## ti = Temp at snow/ice interface #############################
!                                                         f0(up) I
! ...t(1) = Mid point temp of first ice layer................... C
!                                                                E
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetemps1s1i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real htup,htdown
real, dimension(nc) :: con,rhin,rhsn,conb,qmax,sbrine,fl,xxx,excess,conc
real, dimension(nc) :: subl,snmelt,dhb,isubl,ssubl,smax
real, dimension(nc) :: simelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nc), intent(inout) :: dt_avewtemp
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice
real(kind=8), dimension(nc,2:3) :: aa
real(kind=8), dimension(nc,3) :: bb,dd
real(kind=8), dimension(nc,2) :: cc
real, dimension(nc,3) :: ans

! Thickness of each layer
rhsn=1./it_dsn
rhin=1./(it_dic-gammi/cpi)
con =condsnw/it_dsn
conb=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)
conc=condice/max(it_dic,icemin)

! map from generic ice pack to 1 layer + snow
!it_tsurf=it_tsurf
!it_tn0=it_tn0
it_tn1=0.5*(it_tn1+it_tn2)

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gammi
cc(:,1)=-dt*2.*con/gammi
dd(:,1)=it_tsurf+dt*(dt_ftop-pt_egice*lf/lv)/gammi
aa(:,2)=-dt*2.*con*rhsn/cps
bb(:,2)=1.+dt*2.*con*rhsn/cps+dt*2.*conb*rhsn/cps
cc(:,2)=-dt*2.*conb*rhsn/cps
dd(:,2)=it_tn0
aa(:,3)=-dt*2.*conb*rhin/cpi
bb(:,3)=1.+dt*2.*conb*rhin/cpi+dt*2.*conc*rhin/cpi
dd(:,3)=it_tn1+dt*2.*conc*dt_tb*rhin/cpi
call thomas(ans,aa,bb,cc,dd)
it_tsurf(:)=ans(:,1)
it_tn0(:)  =ans(:,2)
it_tn1(:)  =ans(:,3)
fl=2.*condice*(dt_tb-it_tn1)
dhb=dt*(fl-dt_fb)/qice                      ! Excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+(flnew-fl)/(cpi*it_dic-gammi) ! modify temperature if limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflxs=dt_salflxs+dhb*rhoic/dt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
end where

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*lf/(lv*qsnow)        ! net snow+ice sublimation
ssubl=min(subl,it_dsn)                ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic        ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhosn/dt
dt_salflxs=dt_salflxs-isubl*rhosn/dt
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=it_dic-smax
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
end where

! remove brine store for single layer
sbrine=min(it_sto/qice,it_dic)
it_sto=it_sto-sbrine*qice
dt_salflxs=dt_salflxs-sbrine*rhoic/dt
it_dic=it_dic-sbrine
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+sbrine)-cp0*rhoic*dt_avewtemp*sbrine)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn1*cpi*(it_dic+sbrine)-cp0*rhoic*dt_avewtemp*sbrine) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
end where

! the following are snow to ice processes
xxx=it_dic+it_dsn-(rhosn*it_dsn+rhoic*it_dic)/rhowt ! white ice formation
excess=max(it_dsn-xxx,0.)*rhosn/rhowt               ! white ice formation
excess=excess+max(it_dsn-excess*rhowt/rhosn-0.2,0.)*rhosn/rhowt  ! Snow depth limitation and conversion to ice
it_dsn=it_dsn-excess*rhowt/rhosn
it_dic=it_dic+excess*rhowt/rhoic
dt_salflxf=dt_salflxf-excess*rhowt/dt
dt_salflxs=dt_salflxs+excess*rhowt/dt

! Snow melt
snmelt=max(it_tsurf-273.16,0.)*gammi/qsnow
snmelt=min(snmelt,it_dsn)
it_tsurf=it_tsurf-snmelt*qsnow/gammi
dt_salflxf=dt_salflxf-snmelt*rhosn/dt        ! melt fresh water snow (no salt when melting snow)
it_dsn=it_dsn-snmelt
where ( it_dsn>0.05 )
  it_tn0=(it_tn0*cps*(it_dsn+snmelt)-cp0*rhosn*dt_avewtemp*snmelt)/(cps*it_dsn)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn0*cps*(it_dsn+snmelt)-cp0*rhosn*dt_avewtemp*snmelt) &
          /(gammi+cps*it_dsn)
  it_tn0=it_tsurf
end where

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gammi/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gammi
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=it_dic-simelt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gammi+it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt) &
          /(gammi+cpi*it_dic)
  it_tn1=it_tsurf
end where

! test whether to change number of layers
htup=2.*himin
htdown=himin
where ( it_dic>htup )
  ! code to increase number of layers
  dt_nk=2
  !it_tn2=it_tn1
  ! snow depth has not changed, so split ice layer into two
elsewhere ( it_dic<htdown )
  ! code to decrease number of layers
  dt_nk=0
  it_tsurf=(it_tsurf*gammi+it_tn0*cps*it_dsn+it_tn1*cpi*it_dic)/(gammi+cps*it_dsn+cpi*it_dic)
  it_tn0=it_tsurf
  it_tn1=it_tsurf
end where

it_tn2=it_tn1

return
end subroutine icetemps1s1i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and no snow

! Solve for surface ice temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ti = Ice temp for thin layer near surface..................
! -------------------------------------------------------------- I
!                                                         f0(up) C
! ...t(1) = Mid point temp of first ice layer................... E
!                                                                 
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetempi2i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer, dimension(nc), intent(inout) :: dt_nk
real, intent(in) :: dt
real htdown
real, dimension(nc) :: rhin,qmax,qneed,fl,con,gamm,ssubl,isubl,conb
real, dimension(nc) :: subl,simelt,smax,dhb,snmelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nc), intent(inout) :: dt_avewtemp
real, dimension(nc), intent(in) :: pt_egice
real(kind=8), dimension(nc,2:3) :: aa
real(kind=8), dimension(nc,1:3) :: bb,dd
real(kind=8), dimension(nc,1:2) :: cc
real, dimension(nc,3) :: ans
logical, dimension(nc) :: ltest

con =1./(it_dsn/condsnw+0.5*max(it_dic,icemin)/condice)
conb=2.*condice/max(it_dic,icemin)
gamm=cps*it_dsn+gammi  ! for energy conservation

! map from generic ice pack to 2 layer
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0)/gamm
!it_tn1=it_tn1
!it_tn2=it_tn2

! Thickness of each layer
rhin=2./(it_dic-gammi/cpi)

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gamm
cc(:,1)=-dt*2.*con/gamm
dd(:,1)=it_tsurf+dt*(dt_ftop-pt_egice*lf/lv)/gamm
aa(:,2)=-dt*2.*con*rhin/cpi
bb(:,2)=1.+dt*2.*con*rhin/cpi+dt*conb*rhin/cpi
cc(:,2)=-dt*conb*rhin/cpi
dd(:,2)=it_tn1
aa(:,3)=-dt*conb*rhin/cpi
bb(:,3)=1.+dt*conb*rhin/cpi+dt*2.*conb*rhin/cpi
dd(:,3)=it_tn2+dt*2.*conb*dt_tb*rhin/cpi
call thomas(ans,aa,bb,cc,dd)
it_tsurf(:)=ans(:,1)
it_tn1(:)  =ans(:,2)
it_tn2(:)  =ans(:,3)
fl=2.*conb*(dt_tb-it_tn2)
dhb=dt*(fl-dt_fb)/qice                   ! first guess of excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+2.*(flnew-fl)/(cpi*it_dic) ! modify temperature if limit is reached
it_tn2=it_tn2+2.*(flnew-fl)/(cpi*it_dic) ! modify temperature if limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflxs=dt_salflxs+dhb*rhoic/dt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+0.5*(it_tn1+it_tn2)*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhoic/dt
dt_salflxs=dt_salflxs-isubl*rhoic/dt
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl
gamm=gammi+cps*it_dsn

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=it_dic-smax
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+0.5*(it_tn1+it_tn2)*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! update brine store for middle layer
qneed=(dt_timelt-it_tn1)*0.5*cpi*it_dic ! J/m**2
qneed=max(min(qneed,it_sto),0.)
where ( it_dic>himin )
  it_tn1=it_tn1+2.*qneed/(cpi*it_dic)
  it_sto=it_sto-qneed
end where

! Snow melt
ltest = it_tsurf>273.16+0.01.and.it_dsn>icemin
do while ( any(ltest) )
  where ( ltest )
    snmelt=max(it_tsurf-273.16,0.)*gamm/(qsnow+cp0*rhosn*dt_avewtemp-cps*273.16)
    !snmelt=max(it_tsurf-273.16,0.)*gamm/qsnow
    snmelt=min(snmelt,it_dsn)
    it_tsurf=it_tsurf-snmelt*qsnow/gamm
    dt_salflxf=dt_salflxf-snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)
    it_dsn=it_dsn-snmelt
    gamm=gammi+cps*it_dsn
    it_tsurf=(it_tsurf*(gammi+cps*(it_dsn+snmelt))-cp0*rhosn*dt_avewtemp*snmelt)/(gammi+cps*it_dsn)
  end where
  ltest = it_tsurf>273.16+0.01.and.it_dsn>icemin
end do

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gamm
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=it_dic-simelt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
  it_tn2=(it_tn2*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+0.5*(it_tn1+it_tn2)*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
  it_tn2=it_tsurf
end where

! test whether to change number of layers
htdown=2.*himin
where ( it_dic<htdown )
  ! code to decrease number of layers
  dt_nk=1
  ! merge ice layers into one
  it_tn1=0.5*(it_tn1+it_tn2)
  it_tn2=it_tn1
end where

it_tn0=it_tsurf

return
end subroutine icetempi2i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and no snow

! Solve for surface ice temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ti = Ice temp for thin layer near surface..................
! -------------------------------------------------------------- I
!                                                         f0(up) C
! ...t(1) = Mid point temp of first ice layer................... E
!                                                                 
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetempi1i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer, dimension(nc), intent(inout) :: dt_nk
real, intent(in) :: dt
real htup,htdown
real, dimension(nc) :: rhin,qmax,sbrine,fl,con,gamm,ssubl,isubl,conb
real, dimension(nc) :: subl,simelt,smax,dhb,snmelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nc), intent(inout) :: dt_avewtemp
real, dimension(nc), intent(in) :: pt_egice
real(kind=8), dimension(nc,2:2) :: aa
real(kind=8), dimension(nc,1:2) :: bb,dd
real(kind=8), dimension(nc,1:1) :: cc
real, dimension(nc,2) :: ans
logical, dimension(nc) :: ltest

con =1./(it_dsn/condsnw+max(it_dic,icemin)/condice)        ! snow/ice conduction
conb=condice/max(it_dic,icemin)                            ! ice conduction
gamm=gammi+cps*it_dsn                                      ! for energy conservation

! map from generic ice pack to 1 layer
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0)/gamm
it_tn1=0.5*(it_tn1+it_tn2)

! Thickness of layer
rhin=1./(it_dic-gammi/cpi)

! Solve implicit ice temperature matrix
! f0=2.*con*(newt1-newtsurf) is the flux between tsurf and t1
! fl=2.*conb*(dt_tb-newt1) is the flux between t1 and the bottom
bb(:,1)=1.+dt*2.*con/gamm
cc(:,1)=-dt*2.*con/gamm
dd(:,1)=it_tsurf+dt*(dt_ftop-pt_egice*lf/lv)/gamm   
aa(:,2)=-dt*2.*con*rhin/cpi      
bb(:,2)=1.+dt*2.*con*rhin/cpi+dt*2.*conb*rhin/cpi
dd(:,2)=it_tn1+dt*2.*conb*dt_tb*rhin/cpi  
call thomas(ans,aa,bb,cc,dd)
it_tsurf(:)=ans(:,1)
it_tn1(:)  =ans(:,2)
fl=2.*conb*(dt_tb-it_tn1)                   ! flux between t1 and bottom
dhb=dt*(fl-dt_fb)/qice                      ! first guess of excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
flnew=dt_fb+dhb*qice/dt                     ! final excess flux from below
it_tn1=it_tn1+(flnew-fl)/(cpi*it_dic)       ! does nothing unless a limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflxs=dt_salflxs+dhb*rhoic/dt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+it_tn1*cpi*(it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
end where

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhoic/dt
dt_salflxs=dt_salflxs-isubl*rhoic/dt
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl
gamm=gammi+cps*it_dsn

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=it_dic-smax
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+it_tn1*cpi*(it_dic+smax)-cp0*rhoic*dt_avewtemp*smax) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
end where

! remove brine store for single layer
sbrine=min(it_sto/qice,it_dic)
it_sto=it_sto-sbrine*qice
dt_salflxs=dt_salflxs-sbrine*rhoic/dt
it_dic=it_dic-sbrine
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+sbrine)-cp0*rhoic*dt_avewtemp*sbrine)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+it_tn1*cpi*(it_dic+sbrine)-cp0*rhoic*dt_avewtemp*sbrine) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
end where

! Snow melt
ltest = it_tsurf>273.16+0.1.and.it_dsn>icemin
do while ( any(ltest) )
  where ( ltest )
    snmelt=max(it_tsurf-273.16,0.)*gamm/(qsnow+cp0*rhosn*dt_avewtemp-cps*273.16)
    !snmelt=max(it_tsurf-273.16,0.)*gamm/qsnow
    snmelt=min(snmelt,it_dsn)
    it_tsurf=it_tsurf-snmelt*qsnow/gamm
    dt_salflxf=dt_salflxf-snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)
    it_dsn=it_dsn-snmelt
    gamm=gammi+cps*it_dsn
    it_tsurf=(it_tsurf*(gammi+cps*(it_dsn+snmelt))-cp0*rhosn*dt_avewtemp*snmelt)/(gammi+cps*it_dsn)
  end where
  ltest = it_tsurf>273.16+0.1.and.it_dsn>icemin
end do

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gamm
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=it_dic-simelt
where ( it_dic>himin )
  it_tn1=(it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt)/(cpi*it_dic)
elsewhere
  it_tsurf=(it_tsurf*gamm+it_tn1*cpi*(it_dic+simelt)-cp0*rhoic*dt_avewtemp*simelt) &
          /(gamm+cpi*it_dic)
  it_tn1=it_tsurf
end where
  
! test whether to change number of layers
htup=2.*himin
htdown=himin
where ( it_dic>htup )
  dt_nk=2
  !it_tn2=it_tn1
elsewhere ( it_dic<htdown )
  ! code to decrease number of layers
  dt_nk=0
  it_tsurf=(it_tsurf*(cps*it_dsn+gammi)+it_tn1*cpi*it_dic)/(gammi+cps*it_dsn+cpi*it_dic)
  it_tn1=it_tsurf
end where

it_tn0=it_tsurf
it_tn2=it_tn1

return
end subroutine icetempi1i
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for one layer and snow

! Solve for snow surface temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)#S
!                                                                N
! ...ts = Snow temp for thin layer near surface................. O
!                                                                W
! --------------------------------------------------------fs(up)
!                                                                I
! ...ti = Temp at snow/ice interface............................ C
!                                                                E
! ## tb ##########################################################

subroutine icetemps(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,dt_avewtemp,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,f0,tnew
real, dimension(nc) :: subl,snmelt,smax,dhb,isubl,ssubl,gamm,simelt
real, dimension(nc) :: aa,bb,cc
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nc), intent(inout) :: dt_avewtemp
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice
logical, dimension(nc) :: ltest

con=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)         ! conductivity
gamm=cps*it_dsn+gammi+cpi*it_dic                           ! heat capacity

! map from generic ice pack to thin ice
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0+0.5*cpi*it_dic*(it_tn1+it_tn2))/gamm

! Update tsurf based on fluxes from above and below
! MJT notes - we need an implicit form of temperature (tnew) to account for
! fluxes and changes in ice thickness
!tnew=(it_tsurf*gamm/dt+dt_ftop-pt_egice*lf/lv+con*dt_tb)/(gamm/dt+con)    ! predictor temperature for flux calculation
where ( cpi*it_dic>gammi )
  aa=-cpi*con/qice
  bb=gamm/dt+con*(1.+cp0*rhoic*dt_avewtemp/qice)+cpi*(con*dt_tb-dt_fb)/qice
  cc=-it_tsurf*gamm/dt-dt_ftop+pt_egice*lf/lv-con*dt_tb-cp0*rhoic*dt_avewtemp*(con*dt_tb-dt_fb)/qice
  tnew=-2.*cc/(bb+sqrt(bb*bb-4.*aa*cc))
elsewhere
  tnew=(it_tsurf*gamm/dt+dt_ftop-pt_egice*lf/lv+con*dt_tb+cp0*rhoic*dt_avewtemp*(con*dt_tb-dt_fb)/qice) &
      /(gamm/dt+con*(1.+cp0*rhoic*dt_avewtemp/qice))
end where
f0=con*(dt_tb-tnew)                                                       ! first guess of flux from below
dhb=dt*(f0-dt_fb)/qice                                                    ! excess flux converted to change in ice thickness
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
f0=dhb*qice/dt+dt_fb                                                      ! final flux from below
it_tsurf=it_tsurf+dt*(dt_ftop-pt_egice*lf/lv+f0)/gamm                     ! update ice/snow temperature

! Bottom ablation or accretion to balance energy budget
dt_salflxs=dt_salflxs+dhb*rhoic/dt
it_dic=it_dic+dhb
it_tsurf=(it_tsurf*(gammi+cps*it_dsn+cpi*it_dic-dhb)+cp0*rhoic*dt_avewtemp*dhb)  &
        /(gammi+cps*it_dsn+cpi*it_dic)

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhoic/dt
dt_salflxs=dt_salflxs-isubl*rhoic/dt
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! Limit maximum energy stored in brine pockets
smax=min(it_sto/qice,it_dic)
it_sto=it_sto-smax*qice
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=it_dic-smax
it_tsurf=(it_tsurf*(gammi+cps*it_dsn+cpi*it_dic+smax)-cp0*rhoic*dt_avewtemp*smax)       &
        /(gammi+cps*it_dsn+cpi*it_dic)

! Snow melt
ltest = it_tsurf>273.16+0.01.and.it_dsn>icemin
do while ( any(ltest) )
  where ( ltest )
    gamm=gammi+cps*it_dsn+cpi*it_dic
    snmelt=max(it_tsurf-273.16,0.)*gamm/(qsnow+cp0*rhosn*dt_avewtemp-cps*273.16)
    !snmelt=max(it_tsurf-273.16,0.)*gamm/qsnow
    snmelt=min(snmelt,it_dsn)
    it_tsurf=it_tsurf-snmelt*qsnow/gamm
    dt_salflxf=dt_salflxf-snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)
    it_dsn=it_dsn-snmelt
    it_tsurf=(it_tsurf*gamm-cp0*rhosn*dt_avewtemp*snmelt)/(gammi+cps*it_dsn+cpi*it_dic)
  end where
  ltest = it_tsurf>273.16+0.01.and.it_dsn>icemin
end do

! Ice melt
ltest = it_tsurf>dt_timelt+0.01.and.it_dic>icemin
do while ( any(ltest) )
  where ( ltest )
    gamm=gammi+cps*it_dsn+cpi*it_dic
    simelt=max(it_tsurf-dt_timelt,0.)*gamm/(qice+cp0*rhoic*dt_avewtemp-cpi*dt_timelt)
    !simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
    simelt=min(simelt,it_dic)
    it_tsurf=it_tsurf-simelt*qice/gamm
    dt_salflxs=dt_salflxs-simelt*rhoic/dt
    it_dic=it_dic-simelt
    it_tsurf=(it_tsurf*gamm-cp0*rhoic*dt_avewtemp*simelt)/(gammi+cps*it_dsn+cpi*it_dic)
  end where
  ltest = it_tsurf>dt_timelt+0.01.and.it_dic>icemin
end do

! Recalculate thickness index
where ( it_dic>=himin )
  dt_nk=1
end where

! Update unassigned levels
it_tn0=it_tsurf
it_tn1=it_tsurf
it_tn2=it_tsurf

return
end subroutine icetemps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solves matrix for temperature conduction
                     
subroutine thomas(outo,aai,bbi,cci,ddi)

implicit none

integer ii, nlev
real, dimension(:,:), intent(out) :: outo
real(kind=8), dimension(:,2:), intent(in) :: aai
real(kind=8), dimension(:,:), intent(in) :: bbi,ddi
real(kind=8), dimension(:,:), intent(in) :: cci
real(kind=8), dimension(size(outo,1),size(outo,2)) :: cc, dd, ans
real(kind=8), dimension(size(outo,1)) :: n

nlev=size(outo,2)
cc(:,1)=cci(:,1)/bbi(:,1)
dd(:,1)=ddi(:,1)/bbi(:,1)

do ii=2,nlev-1
  n=bbi(:,ii)-cc(:,ii-1)*aai(:,ii)
  cc(:,ii)=cci(:,ii)/n
  dd(:,ii)=(ddi(:,ii)-dd(:,ii-1)*aai(:,ii))/n
end do
n=bbi(:,nlev)-cc(:,nlev-1)*aai(:,nlev)
dd(:,nlev)=(ddi(:,nlev)-dd(:,nlev-1)*aai(:,nlev))/n
ans(:,nlev)=dd(:,nlev)
do ii=nlev-1,1,-1
  ans(:,ii)=dd(:,ii)-cc(:,ii)*ans(:,ii+1)
end do
outo=real(ans)

return
end subroutine thomas
                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dt,atm_sg,atm_rg,atm_rnd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v,atm_temp,atm_qg, &
                   atm_ps,atm_zmin,atm_zmins,d_ftop,d_tb,d_fb,d_timelt,d_nk,                             &
                   d_ndsn,d_ndic,d_nsto,d_delstore,d_imass,atm_oldu,atm_oldv,diag)

implicit none

integer, intent(in) :: diag
integer itr
real, dimension(wfull), intent(in) :: atm_sg,atm_rg,atm_rnd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v
real, dimension(wfull), intent(in) :: atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
real, dimension(wfull), intent(inout) :: d_ftop,d_tb,d_fb,d_timelt,d_ndsn,d_ndic
real, dimension(wfull), intent(inout) :: d_nsto,d_delstore,d_imass,atm_oldu,atm_oldv
integer, dimension(wfull), intent(inout) :: d_nk
real, intent(in) :: dt
real, dimension(wfull) :: qsat,dqdt,ri,rho,srcp,tnew,qsatnew,gamm,bot
real, dimension(wfull) :: fm,fh,fq,af,aft,afq
real, dimension(wfull) :: den,sig,root
real, dimension(wfull) :: alb,qmax,eye
real, dimension(wfull) :: uu,vv,du,dv,vmagn,icemagn
real, dimension(wfull) :: ustar,g,h,dgu,dgv,dhu,dhv,det
real, dimension(wfull) :: newiu,newiv,dtsurf
real factch

! Prevent unrealistic fluxes due to poor input surface temperature
dtsurf=min(ice%tsurf,273.2)
uu=atm_u-ice%u
vv=atm_v-ice%v
vmagn=sqrt(max(uu*uu+vv*vv,1.E-4))
sig=exp(-grav*atm_zmins/(rdry*atm_temp))
srcp=sig**(rdry/cpair)
rho=atm_ps/(rdry*dtsurf)

dgice%zo=zoseaice ! Mk3.6 (0.01m), CCAM sflux (0.001m), CICE (0.0005m)
factch=factchseaice ! following CSIRO9, CICE (1.) or CCAM sflux (2.72)
dgice%zoh=dgice%zo/(factch*factch)
dgice%zoq=dgice%zoh
af=vkar*vkar/(log(atm_zmin/dgice%zo)*log(atm_zmin/dgice%zo))
aft=vkar*vkar/(log(atm_zmin/dgice%zo)*log(atm_zmin/dgice%zoh))
afq=aft

! number of (thick) ice layers
d_nk=min(int(d_ndic/himin),2)

! radiation
alb=     atm_vnratio*(dgice%visdiralb*atm_fbvis+dgice%visdifalb*(1.-atm_fbvis))+ &
    (1.-atm_vnratio)*(dgice%visdifalb*atm_fbnir+dgice%visdifalb*(1.-atm_fbnir))
qmax=max(qice*0.5*(d_ndic-himin),1.E-10)
eye=0.35*max(1.-d_ndsn/icemin,0.)*max(1.-d_nsto/qmax,0.)*max(min(d_ndic/himin-1.,1.),0.)
where ( d_ndic>icemin )
  d_delstore=atm_sg*(1.-alb)*eye
  d_nsto=d_nsto+dt*d_delstore
elsewhere
  d_delstore=0.
end where

where ( d_nk>0 .and. d_ndsn>0.05 )
  gamm=gammi
elsewhere
  gamm=cps*d_ndsn+gammi
end where

! water temperature at bottom of ice
d_tb=water%temp(:,1)+wrtemp

! Explicit estimate of fluxes
call getqsat(qsat,dqdt,dtsurf,atm_ps)
ri=min(grav*(atm_zmin**2/atm_zmins)*(1.-dtsurf*srcp/atm_temp)/vmagn**2,rimax)
where (ri>=0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
  fq=fh
elsewhere        ! ri is -ve
  root=sqrt(-ri*atm_zmin/dgice%zo)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  root=sqrt(-ri*atm_zmins/dgice%zo)
  den=1.+chs*2.*bprm*aft*factch*root
  fh=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*afq*factch*root
  fq=1.-2.*bprm*ri/den
end where
! egice is for evaporating (lv).  Melting is included with lf.
dgice%wetfrac=max(1.+.008*min(dtsurf-273.16,0.),0.)
dgice%cd =af*fm
dgice%cdh=aft*fh
dgice%cdq=afq*fq
dgice%fg=rho*dgice%cdh*cpair*vmagn*(dtsurf-atm_temp/srcp)
dgice%fg=min(max(dgice%fg,-3000.),3000.)
dgice%eg=dgice%wetfrac*rho*dgice%cdq*lv*vmagn*(qsat-atm_qg)
dgice%eg=min(dgice%eg,d_ndic*qice/(lf*dt))
dgice%eg=min(max(dgice%eg,-3000.),3000.)

! MJT notes - use dtsurf for outgoing longwave for consistency with radiation code
! energy conservation is violated if initial conditions are poor
d_ftop=-dgice%fg-dgice%eg+atm_rg-emisice*sbconst*dtsurf**4+atm_sg*(1.-alb)*(1.-eye) ! first guess
d_ftop=d_ftop+lf*atm_rnd ! converting any rain to snowfall over ice
bot=rho*vmagn*(dgice%cdh*cpair+dgice%cdq*dgice%wetfrac*dqdt*lv)

! iterative method to estimate ice velocity after stress from wind and currents are applied
newiu=ice%u
newiv=ice%v
do itr=1,10 ! max iterations
  uu=atm_u-newiu
  vv=atm_v-newiv
  du=fluxwgt*water%u(:,1)+(1.-fluxwgt)*atm_oldu-newiu
  dv=fluxwgt*water%v(:,1)+(1.-fluxwgt)*atm_oldv-newiv
  vmagn=sqrt(max(uu*uu+vv*vv,1.E-4))
  icemagn=sqrt(max(du*du+dv*dv,1.E-8))
  ! MJT notes - use rhowt reference density for Boussinesq fluid approximation
  g=ice%u-newiu+dt*(rho*dgice%cd*vmagn*uu+rhowt*0.00536*icemagn*du)/d_imass
  h=ice%v-newiv+dt*(rho*dgice%cd*vmagn*vv+rhowt*0.00536*icemagn*dv)/d_imass
  dgu=-1.-dt*(rho*dgice%cd*vmagn*(1.+(uu/vmagn)**2)+rhowt*0.00536*icemagn*(1.+(du/icemagn)**2))/d_imass
  dhu=-dt*(rho*dgice%cd*uu*vv/vmagn+rhowt*0.00536*du*dv/icemagn)/d_imass
  dgv=dhu
  dhv=-1.-dt*(rho*dgice%cd*vmagn*(1.+(vv/vmagn)**2)+rhowt*0.00536*icemagn*(1.+(dv/icemagn)**2))/d_imass
  ! Min det is around 1.
  det=dgu*dhv-dgv*dhu
  newiu=newiu-0.9*( g*dhv-h*dgv)/det
  newiv=newiv-0.9*(-g*dhu+h*dgu)/det
end do

! momentum transfer
uu=atm_u-newiu
vv=atm_v-newiv
du=fluxwgt*water%u(:,1)+(1.-fluxwgt)*atm_oldu-newiu
dv=fluxwgt*water%v(:,1)+(1.-fluxwgt)*atm_oldv-newiv
vmagn=sqrt(max(uu*uu+vv*vv,1.E-4))
icemagn=sqrt(max(du*du+dv*dv,1.E-8))
dgice%tauxica=rho*dgice%cd*vmagn*uu
dgice%tauyica=rho*dgice%cd*vmagn*vv
! MJT notes - use rhowt reference density for Boussinesq fluid approximation
dgice%tauxicw=-rhowt*0.00536*icemagn*du
dgice%tauyicw=-rhowt*0.00536*icemagn*dv
!dgice%tauxicw=(ice%u-newiu)*d_imass/dt+dgice%tauxica
!dgice%tauyicw=(ice%v-newiv)*d_imass/dt+dgice%tauyica
ustar=sqrt(sqrt(max(dgice%tauxicw*dgice%tauxicw+dgice%tauyicw*dgice%tauyicw,0.))/rhowt)
ustar=max(ustar,5.E-4)
d_fb=cp0*rhowt*0.006*ustar*(d_tb-d_timelt)
d_fb=min(max(d_fb,-1000.),1000.)  

! Re-calculate fluxes to prevent overshoot (predictor-corrector)
! MJT notes - use dtsurf for outgoing longwave for consistency with radiation code
tnew=min(dtsurf+d_ftop/(gamm/dt+bot),273.2)
tnew=0.5*(tnew+dtsurf)
call getqsat(qsatnew,dqdt,tnew,atm_ps)
dgice%fg=rho*dgice%cdh*cpair*vmagn*(tnew-atm_temp/srcp)
dgice%fg=min(max(dgice%fg,-3000.),3000.)
dgice%eg=dgice%wetfrac*rho*dgice%cdq*lv*vmagn*(qsatnew-atm_qg)
dgice%eg=min(dgice%eg,d_ndic*qice*lv/(lf*dt))
dgice%eg=min(max(dgice%eg,-3000.),3000.)
d_ftop=-dgice%fg-dgice%eg+atm_rg-emisice*sbconst*dtsurf**4+atm_sg*(1.-alb)*(1.-eye)
! Add flux of heat due to converting any rain to snowfall over ice
d_ftop=d_ftop+lf*atm_rnd ! rain (mm/sec) to W/m**2

return
end subroutine iceflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine screen diagnostics

subroutine scrncalc(atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,diag)

implicit none

integer, intent(in) :: diag
real, dimension(wfull), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
real, dimension(wfull) :: tscrn,qgscrn,uscrn,u10,dumtemp
real, dimension(wfull) :: smixr,qsat,dqdt,atu,atv,dmag

! water
dumtemp=water%temp(:,1)+wrtemp
call getqsat(qsat,dqdt,dumtemp,atm_ps)
if (zomode==0) then
  smixr=0.98*qsat
else
  smixr=qsat
end if
atu=atm_u-water%u(:,1)
atv=atm_v-water%v(:,1)
call scrntile(dgscrn%temp,dgscrn%qg,uscrn,u10,dgwater%zo,dgwater%zoh,dgwater%zoq,dumtemp, &
    smixr,atu,atv,atm_temp,atm_qg,atm_zmin,atm_zmins,diag)
dmag=sqrt(max(atu*atu+atv*atv,1.E-4))
atu=(atm_u-water%u(:,1))*uscrn/dmag+water%u(:,1)
atv=(atm_v-water%v(:,1))*uscrn/dmag+water%v(:,1)
dgscrn%u2=sqrt(atu*atu+atv*atv)
atu=(atm_u-water%u(:,1))*u10/dmag+water%u(:,1)
atv=(atm_v-water%v(:,1))*u10/dmag+water%v(:,1)
dgscrn%u10=sqrt(atu*atu+atv*atv)

! ice
call getqsat(qsat,dqdt,ice%tsurf,atm_ps)
smixr=dgice%wetfrac*qsat+(1.-dgice%wetfrac)*min(qsat,atm_qg)
atu=atm_u-ice%u
atv=atm_v-ice%v
call scrntile(tscrn,qgscrn,uscrn,u10,dgice%zo,dgice%zoh,dgice%zoq,ice%tsurf, &
    smixr,atu,atv,atm_temp,atm_qg,atm_zmin,atm_zmins,diag)
dgscrn%temp=(1.-ice%fracice)*dgscrn%temp+ice%fracice*tscrn
dgscrn%qg=(1.-ice%fracice)*dgscrn%qg+ice%fracice*qgscrn
dmag=sqrt(max(atu*atu+atv*atv,1.E-4))
atu=(atm_u-ice%u)*uscrn/dmag+ice%u
atv=(atm_v-ice%v)*uscrn/dmag+ice%v
dgscrn%u2=(1.-ice%fracice)*dgscrn%u2+ice%fracice*sqrt(atu*atu+atv*atv)
atu=(atm_u-ice%u)*u10/dmag+ice%u
atv=(atm_v-ice%v)*u10/dmag+ice%v
dgscrn%u10=(1.-ice%fracice)*dgscrn%u10+ice%fracice*sqrt(atu*atu+atv*atv)

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! screen diagnostic for individual tile

subroutine scrntile(tscrn,qgscrn,uscrn,u10,zo,zoh,zoq,stemp,smixr,atm_u,atm_v,atm_temp,atm_qg,atm_zmin,atm_zmins,diag)
      
implicit none

integer, intent(in) :: diag
integer ic
real, dimension(wfull), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_zmin,atm_zmins
real, dimension(wfull), intent(in) :: zo,zoh,zoq,stemp,smixr
real, dimension(wfull), intent(out) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: umag,sig
real, dimension(wfull) :: lzom,lzoh,lzoq,thetav,sthetav
real, dimension(wfull) :: thetavstar,z_on_l,zs_on_l,z0_on_l,z0s_on_l,zt_on_l,zq_on_l
real, dimension(wfull) :: pm0,ph0,pq0,pm1,ph1,integralm,integralh,integralq
real, dimension(wfull) :: ustar,qstar,z10_on_l
real, dimension(wfull) :: neutrals,neutral,neutral10,pm10
real, dimension(wfull) :: integralm10,tstar,scrp
integer, parameter ::  nc     = 5
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

umag=max(sqrt(atm_u*atm_u+atm_v*atm_v),0.01)
sig=exp(-grav*atm_zmins/(rdry*atm_temp))
scrp=sig**(rdry/cpair)
thetav=atm_temp*(1.+0.61*atm_qg)/scrp
sthetav=stemp*(1.+0.61*smixr)

! Roughness length for heat
lzom=log(atm_zmin/zo)
lzoh=log(atm_zmins/zoh)
lzoq=log(atm_zmins/zoq)

! Dyer and Hicks approach 
thetavstar=vkar*(thetav-sthetav)/lzoh
ustar=vkar*umag/lzom
do ic=1,nc
  z_on_l  = atm_zmin*vkar*grav*thetavstar/(thetav*ustar*ustar)
  z_on_l  = min(z_on_l,10.)
  zs_on_l = z_on_l*atm_zmins/atm_zmin  
  z0_on_l = z_on_l*zo/atm_zmin
  zt_on_l = z_on_l*zoh/atm_zmin
  zq_on_l = z_on_l*zoq/atm_zmin
  where (z_on_l<0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pq0     = (1.-16.*zq_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*zs_on_l)**(-0.5) !=pq1
    integralm = lzom-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = lzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
    integralq = lzoq-2.*log((1.+1./ph1)/(1.+1./pq0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
           +b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
           +b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5 &
           +b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l) &
           +b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*zs_on_l)**1.5 &
           +b_1*(zs_on_l-(c_1/d_1))*exp(-d_1*zs_on_l) &
           +b_1*c_1/d_1-1.) !=pq1
    pq0 = -((1.+(2./3.)*a_1*zq_on_l)**1.5 &
           +b_1*(zq_on_l-(c_1/d_1))*exp(-d_1*zq_on_l) &
           +b_1*c_1/d_1-1.)
    integralm = lzom-(pm1-pm0)
    integralh = lzoh-(ph1-ph0)
    integralq = lzoq-(ph1-pq0)
  endwhere
  integralm = max(integralm,1.e-10)
  integralh = max(integralh,1.e-10)
  integralq = max(integralq,1.e-10)
  thetavstar = vkar*(thetav-sthetav)/integralh
  ustar      = vkar*umag/integralm
end do
tstar = vkar*(atm_temp-stemp)/integralh
qstar = vkar*(atm_qg-smixr)/integralq
      
! estimate screen diagnostics
z0s_on_l  = z0*zs_on_l/atm_zmins
z0_on_l   = z0*z_on_l/atm_zmin
z10_on_l  = z10*z_on_l/atm_zmin
z0s_on_l  = min(z0s_on_l,10.)
z0_on_l   = min(z0_on_l,10.)
z10_on_l  = min(z10_on_l,10.)
neutrals  = log(atm_zmins/z0)
neutral   = log(atm_zmin/z0)
neutral10 = log(atm_zmin/z10)
where (z_on_l<0.)
  ph0     = (1.-16.*z0s_on_l)**(-0.50)
  ph1     = (1.-16.*zs_on_l)**(-0.50)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm10    = (1.-16.*z10_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralh = neutrals-2.*log((1.+1./ph1)/(1.+1./ph0))
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                 -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                 +2.*(atan(1./pm1)-atan(1./pm0))
  integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                -log((1.+1./pm1**2)/(1.+1./pm10**2)) &
                +2.*(atan(1./pm1)-atan(1./pm10))     
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  ph0  = -((1.+(2./3.)*a_1*z0s_on_l)**1.5 &
         +b_1*(z0s_on_l-(c_1/d_1)) &
         *exp(-d_1*z0s_on_l)+b_1*c_1/d_1-1.)
  ph1  = -((1.+(2./3.)*a_1*zs_on_l)**1.5 &
         +b_1*(zs_on_l-(c_1/d_1)) &
         *exp(-d_1*zs_on_l)+b_1*c_1/d_1-1.)
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
         +b_1*c_1/d_1)
  pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1)) &
         *exp(-d_1*z10_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
         +b_1*c_1/d_1)
  integralh = neutrals-(ph1-ph0)
  integralm = neutral-(pm1-pm0)
  integralm10 = neutral10-(pm1-pm10)
endwhere
integralh = max(integralh,1.e-10)
integralm = max(integralm,1.e-10)
integralm10 = max(integralm10,1.e-10)
tscrn  = atm_temp-tstar*integralh/vkar
qgscrn = atm_qg-qstar*integralh/vkar
qgscrn = max(qgscrn,1.E-4)

uscrn=max(umag-ustar*integralm/vkar,0.)
u10  =max(umag-ustar*integralm10/vkar,0.)

return
end subroutine scrntile

end module mlo
