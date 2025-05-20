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
    
! This is a 1D, mixed layer ocean model based on Large, et al (1994), for ensemble regional climate simulations.
! In CCAM this module interfaces withmlodynamics.f90 for river routing, diffusion and advection routines.

! This version supports a thermodynamic model of sea ice based on O'Farrell from Mk3.5.  We have included a
! free surface so that lakes can change depth, etc.

! This version can assimilate SSTs from GCMs, using a convolution based digital filter (see nesting.f90),
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

#ifdef CCAM
use newmpar_m, only : imax, ntiles
#endif

implicit none

private
public mloeval_work, mlo_calc_k, vgrid, calcdensity, calcmelt, mlocheck
public mlonewice_work, mlo_ema_uvw, mlo_ema_ts, mlo_ema_reset

public wlev,zomode,wrtemp,wrtrho,mxd,mindep,minwater,zoseaice,factchseaice,otaumode,mlosigma
public oclosure,pdl,pdu,usepice,minicemass,cdbot,cp0,ominl,omaxl,mlo_adjeta,mlo_limitsal
public mlo_timeave_length,kemaxdt,mlo_step,mlo_uvcoupl,fluxwgt,delwater
public alphavis_seaice, alphanir_seaice
public alphavis_seasnw, alphanir_seasnw
public omink,omineps
public k_mode,eps_mode,limitL,fixedce3
public nops,nopb,fixedstabfunc

public minsfc, minsal, maxsal, icemax, emisice, incradgam, gammi, cps, cpi
public lv, ls

public turbdata, depthdata, waterdata, icedata, dgwaterdata, dgicedata, dgscrndata


type turbdata
  real, dimension(:,:), allocatable :: k            ! turbulent kinnetic energy
  real, dimension(:,:), allocatable :: eps          ! turbulent dissipation rate
end type

type depthdata
  real, dimension(:,:), allocatable :: depth, depth_hl  ! Column depth (m)
  real, dimension(:,:), allocatable :: dz, dz_hl        ! Column thickness (m)
  real, dimension(:), allocatable :: f                  ! Coriolis
  integer, dimension(:), allocatable :: ibot            ! index of bottom layer  
  logical :: data_allocated
end type depthdata

type waterdata
  real, dimension(:,:), allocatable :: temp         ! water layer temperature with respect to wrtemp (K)
  real, dimension(:,:), allocatable :: sal          ! water layer salinity (PSU)
  real, dimension(:,:), allocatable :: u            ! water u-current (m/s)
  real, dimension(:,:), allocatable :: v            ! water v-current (m/s)
  real, dimension(:,:), allocatable :: w            ! water v-current (m/s)
  real, dimension(:), allocatable :: eta            ! water free surface height (m)
  real, dimension(:), allocatable :: ubot           ! water u-current at bottom from previous time-step (m/s)
  real, dimension(:), allocatable :: vbot           ! water v-current at bottom from previous time-step (m/s)  
  real, dimension(:), allocatable :: utop           ! water u-current at top from previous time-step (m/s)
  real, dimension(:), allocatable :: vtop           ! water v-current at top from previous time-step (m/s)  
  real, dimension(:,:), allocatable :: u_ema        ! exponential average of u
  real, dimension(:,:), allocatable :: v_ema        ! exponential average of v
  real, dimension(:,:), allocatable :: w_ema        ! exponential average of w 
  real, dimension(:,:), allocatable :: temp_ema     ! exponential average of temp
  real, dimension(:,:), allocatable :: sal_ema      ! exponential average of sal
  real, dimension(:,:), allocatable :: dudz         ! for shear calculation
  real, dimension(:,:), allocatable :: dvdz         ! for shear calculation
  real, dimension(:,:), allocatable :: dwdx         ! for shear calculation
  real, dimension(:,:), allocatable :: dwdy         ! for shear calculation  
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
  real, dimension(:), allocatable :: cd             ! water drag coeff for momentum at surface
  real, dimension(:), allocatable :: cdh            ! water drag coeff for heat at surface
  real, dimension(:), allocatable :: cdq            ! water drag coeff for moisture at surface
  real, dimension(:), allocatable :: cd_bot         ! water drag coeff for mementum at bottom
  real, dimension(:), allocatable :: umod           ! relative speed between air and surface
  real, dimension(:), allocatable :: fg             ! water sensible heat flux (W/m2)
  real, dimension(:), allocatable :: eg             ! water latent heat flux (W/m2)
  real, dimension(:), allocatable :: taux           ! wind/water u-component stress (N/m2)
  real, dimension(:), allocatable :: tauy           ! wind/water v-component stress (N/m2)
  real, dimension(:,:), allocatable :: rho          ! water density (kg/m3)
  real, dimension(:,:), allocatable :: alpha        ! gradient of density with respect to temperature
  real, dimension(:,:), allocatable :: beta         ! gradient of density with respect to salinity
  real, dimension(:,:), allocatable :: rad          ! heating profile due to solar radiation
  real, dimension(:), allocatable :: wt0            ! flux for surface theta
  real, dimension(:), allocatable :: wt0_rad        ! radiation component of wt0
  real, dimension(:), allocatable :: wt0_melt       ! melt component of wt0
  real, dimension(:), allocatable :: wt0_eg         ! latent heat flux component of wt0
  real, dimension(:), allocatable :: wt0_fb         ! ice component of wt0
  real, dimension(:), allocatable :: ws0            ! flux for surface salinity
  real, dimension(:), allocatable :: ws0_subsurf    ! flux for sub-surface salinity 
  real, dimension(:), allocatable :: wu0            ! flux for surface u
  real, dimension(:), allocatable :: wv0            ! flux for surface v
  !real, dimension(:), allocatable :: deleng        ! Change in energy stored
end type dgwaterdata

type dgicedata
  real, dimension(:), allocatable :: wetfrac        ! ice area cover of liquid water
  real, dimension(:), allocatable :: visdiralb      ! ice VIS direct albedo
  real, dimension(:), allocatable :: visdifalb      ! ice VIS difuse albedo
  real, dimension(:), allocatable :: nirdiralb      ! ice NIR direct albedo
  real, dimension(:), allocatable :: nirdifalb      ! ice NIR direct albedo
  real, dimension(:), allocatable :: cd             ! ice drag coeff for momentum at surface
  real, dimension(:), allocatable :: cdh            ! ice drag coeff for heat at surface
  real, dimension(:), allocatable :: cdq            ! ice drag coeff for moisture at surface
  real, dimension(:), allocatable :: cd_bot         ! ice drag coeff for momentum at bottom
  real, dimension(:), allocatable :: umod           ! relative speed between air and surface
  real, dimension(:), allocatable :: fg             ! ice sensible heat flux (W/m2)
  real, dimension(:), allocatable :: eg             ! ice latent heat flux (W/m2)
  real, dimension(:), allocatable :: tauxica        ! water/ice u-component stress (N/m2)
  real, dimension(:), allocatable :: tauyica        ! water/ice v-component stress (N/m2)
  real, dimension(:), allocatable :: tauxicw        ! water/ice u-component stress (N/m2)
  real, dimension(:), allocatable :: tauyicw        ! water/ice v-component stress (N/m2)
  real, dimension(:), allocatable :: imass          ! mass of ice (kg)
  !real(kind=8), dimension(:), allocatable :: deleng         ! Change in energy stored
end type dgicedata

type dgscrndata
  real, dimension(:), allocatable :: temp           ! 2m air temperature (K)    
  real, dimension(:), allocatable :: qg             ! 2m water vapor mixing ratio (kg/kg)
  real, dimension(:), allocatable :: u2             ! 2m wind speed (m/s)
  real, dimension(:), allocatable :: u10            ! 10m wind speed (m/s)
end type dgscrndata

! model
integer, save :: wlev         = 20        ! Number of water layers
integer, save :: zomode       = 2         ! roughness calculation (0=Charnock (CSIRO9), 1=Charnock (zot=zom), 2=Beljaars, 3=Moon)
integer, save :: otaumode     = 0         ! momentum coupling (0=Explicit, 1=Implicit)
integer, save :: mlosigma     = 6         ! vertical levels (4=zstar-cubic, 5=zstar-quad, 6=zstar-gotm, 7=zstar-linear)
integer, save :: oclosure     = 0         ! 0=kpp, 1=k-eps
integer, save :: usepice      = 0         ! include ice in surface pressure (0=without ice, 1=with ice)
integer, save :: mlo_adjeta   = 1         ! allow adjustment to surface outside dynamics (0=off, 1=all)
integer, save :: mlo_limitsal = 0         ! limit salinity to maxsal when loading data ( 0=off, 1=limit )
integer, save :: mlo_step     = 0         ! ocean floor (0=full-step, 1=pseduo partial-step, 2=partial-step)
integer, save :: mlo_uvcoupl  = 1         ! wind coupling (0=off, 1=on)
real, save :: pdu    = 2.7                ! zoom factor near the surface for mlosigma==gotm
real, save :: pdl    = 0.0                ! zoom factor near the bottom for mlosigma==gotm
real, save :: kemaxdt = 120.              ! max time-step for k-e coupling

! kpp parameters
integer, parameter :: incradbf  = 1       ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 1       ! include shortwave in non-local term

! k-eps parameters
integer, save :: k_mode        = 2        ! 0=fully explicit k, 1=implicit k, 2=implicit k & pb
integer, save :: eps_mode      = 2        ! 0=fully explicit eps, 1=implicit eps, 2=implicit eps & pb
integer, save :: limitL        = 1        ! 0=no length scale limit, 1=limit length scale
integer, save :: fixedce3      = 0        ! 0=dynamic ce3, 1=fixed ce3
integer, save :: nops          = 0        ! 0=calculate shear production, 1=no shear production
integer, save :: nopb          = 0        ! 0=calculate buoyancy production, 1=no buoyancy production
integer, save :: fixedstabfunc = 0        ! 0=dynamic stability functions, 1=fixed stability functions
real, save :: omink   = 1.e-8             ! minimum k
real, save :: omineps = 1.e-11            ! minimum eps
real, save :: ominl  = 1.e-2              ! minimum L
real, save :: omaxl  = 1.e3               ! maximum L

! model parameters
real, save :: mxd      = 5002.18          ! max depth (m)
real, save :: mindep   = 1.               ! thickness of first layer (m)
real, save :: minwater = 10.              ! absolute minimum water depth (m)
real, save :: delwater = 100.             ! allowed anomoly depth (m)
real, parameter :: ric     = 0.3          ! critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1          ! ratio of surface layer and mixed layer thickness
real, parameter :: minsfc  = 1.           ! minimum thickness to average surface layer properties (m)
real, parameter :: minsal  = 28.          ! minimum non-zero salinity for error checking (PSU)
real, parameter :: maxsal  = 60.          ! maximum salinity used in density and melting point calculations (PSU)
real, parameter :: mu_1    = 23.          ! VIS depth (m) - Type I
real, parameter :: mu_2    = 0.35         ! NIR depth (m) - Type I
real, save :: fluxwgt = 0.7               ! time filter for flux calculations
real, parameter :: wrtemp  = 290.         ! water reference temperature (K)
real, parameter :: wrtrho  = 1035.        ! water reference density (kg m^-3)
! physical parameters
real, parameter :: vkar    = 0.4          ! von Karman constant
real, parameter :: lv      = 2.501e6      ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf      = 3.337e5      ! Latent heat of fusion (J kg^-1)
real, parameter :: ls      = lv + lf      ! Latent heat of sublimation (J kg^-1)
real, parameter :: grav    = 9.8          ! graviational constant (m s^-2)
real, parameter :: sbconst = 5.67e-8      ! Stefan-Boltzmann constant
real, parameter :: cdbot   = 2.4e-3       ! bottom drag coefficent
real, parameter :: utide   = 0.05         ! Tide velocity for bottom drag (m/s)
real, parameter :: cpair   = 1004.64      ! Specific heat of dry air at const P
real, parameter :: rdry    = 287.04       ! Specific gas const for dry air
real, parameter :: rvap    = 461.5        ! Gas constant for water vapor
! water parameters
real, parameter :: cp0   = 3990.          ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: salwt = 34.72          ! reference water salinity (PSU)
! ice parameters
real, save      :: alphavis_seaice = 0.85 ! visible seaice albedo
real, save      :: alphanir_seaice = 0.45 ! near-IR seaice albedo
real, save      :: alphavis_seasnw = 0.95 ! visible seaice albedo
real, save      :: alphanir_seasnw = 0.65 ! near-IR seaice albedo
real, save      :: zoseaice     = 0.0005  ! roughnes length for sea-ice (m)
real, save      :: factchseaice = 1.      ! =sqrt(zo/zoh) for sea-ice
real, parameter :: himin        = 0.1     ! minimum ice thickness for multiple layers (m)
real, parameter :: icemin       = 0.01    ! minimum ice thickness (m)
real, parameter :: icemax       = 6.      ! maximum ice thickness (m)
real, parameter :: rhoic        = 900.    ! ice density (kg/m3)
real, parameter :: rhosn        = 330.    ! snow density (kg/m3)
real, parameter :: qice    = lf*rhoic     ! latent heat of fusion for ice (J m^-3)
real, parameter :: qsnow   = lf*rhosn     ! latent heat of fusion for snow (J m^-3)
real, parameter :: cpi     = 1.8837e6     ! specific heat ice  (J/m**3/K)
real, parameter :: cps     = 6.9069e5     ! specific heat snow (J/m**3/K)
real, parameter :: condice = 2.03439      ! conductivity ice
real, parameter :: condsnw = 0.30976      ! conductivity snow
real, parameter :: gammi = 0.5*cpi*himin  ! specific heat*depth (for ice/snow) (J m^-2 K^-1)
real, parameter :: emisice   = 1.         ! emissivity of ice (0.95?)
real, parameter :: maxicesal = 4.         ! maximum salinity for sea-ice (PSU)
real, parameter :: minicemass = 100.      ! minimum ice mass
! stability function parameters
real, parameter :: bprm=5.                ! 4.7 in rams
real, parameter :: chs=2.6                ! 5.3 in rams
real, parameter :: cms=5.                 ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm

! Time averaging
real, save :: mlo_timeave_length = 0. ! Time period for averaging source terms (Ps, Pb, Pt) in seconds
                                      ! 0 indicates alpha=1.


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define vertical grid

subroutine vgrid(wlin,depin,depthout,depth_hlout)

implicit none

integer, intent(in) :: wlin
integer ii
real, intent(in) :: depin
real, dimension(wlin), intent(out) :: depthout
real, dimension(wlin+1), intent(out) :: depth_hlout
real x, al, bt

x = real(wlin)
select case(mlosigma)
    
  case(4) ! Adcroft and Campin 2003 - cubic
    al = (mindep*x-mxd)/(x-x**3)     ! z* levels
    bt = (mindep*x**3-mxd)/(x**3-x)  
    do ii = 1,wlin+1
      x = real(ii-1)
      depth_hlout(ii) = al*x**3 + bt*x ! ii is for half level ii-0.5
    end do
    
  case(5) ! Adcroft and Campin 2003 - quadratic  
    al = (mxd-mindep*x)/(x**2-x)      ! z* levels
    bt = (mxd-mindep*x*x)/(x-x**2)
    do ii = 1,wlin+1
      x = real(ii-1)
      depth_hlout(ii) = al*x**2 + bt*x ! ii is for half leel ii-0.5
    end do
    
  case(6) ! Adcroft and Campin 2003 - gotm dynamic
    do ii = 1,wlin+1
      x = real(ii-1)
      depth_hlout(ii) = mxd*(tanh((pdu+pdl)*x/wlin -pdu) + tanh(pdu))/(tanh(pdu)+tanh(pdl))
    end do

  case(7) ! Adcroft and Campin 2003 - linear
    do ii = 1,wlin+1
      x = real(ii-1)
      depth_hlout(ii) = x*mxd/wlin
    end do
    
  case default
    write(6,*) "ERROR: Unknown option mlosigma = ",mlosigma
    stop
    
end select
    
! calculate cell mid-points  
do ii = 1,wlin
  depthout(ii) = 0.5*(depth_hlout(ii)+depth_hlout(ii+1))
end do

select case(mlo_step)
  case(0)
    ! full step version
    do ii = 1,wlin
      if ( depthout(ii)>depin ) then
        depth_hlout(ii+1) = depth_hlout(ii)
      end if
    end do
  
  case(1)  
    ! pseduo partial step version (cell point at full step)
    if ( depin>1.e-4 ) then      
      do ii = 1,wlin
        depth_hlout(ii+1) = min( depth_hlout(ii+1), max(depin,depthout(1)+0.1) )
        if ( depthout(ii)>depin .and. ii>1 ) then
          ! avoids thin layers by extending the previous layer  
          depth_hlout(ii) = depth_hlout(ii+1)
        end if
      end do
    else
      depth_hlout(:) = 0.
    end if    

  case(2)
    ! partial step version (cell point at midpoint)
    if ( depin>1.e-4 ) then      
      do ii = 1,wlin
        depth_hlout(ii+1) = min( depth_hlout(ii+1), max(depin,depthout(1)+0.1) )
        ! avoid thin layers by extending the previous layer
        if ( depth_hlout(ii+1)-depth_hlout(ii)<mindep .and. &
             depth_hlout(ii+1)-depth_hlout(ii)>1.e-4  .and. &
             ii>1 ) then
          depth_hlout(ii) = depth_hlout(ii+1)
        end if
      end do
      do ii = 1,wlev
        if ( depth_hlout(ii+1)-depth_hlout(ii)>1.e-4 ) then
          depthout(ii) = 0.5*(depth_hlout(ii)+depth_hlout(ii+1))
        end if  
      end do      
    else
      depth_hlout(:) = 0.
    end if

  case default
    write(6,*) "ERROR: Unknown option mlo_step = ",mlo_step
    stop
    
end select
    
return
end subroutine vgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water and ice

subroutine mloeval_work(dt,atm_zmin,atm_zmins,atm_sg,atm_rg,atm_rnd,atm_snd,atm_u,atm_v,   &
                   atm_temp,atm_qg,atm_ps,                                                 &
                   atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow,diag,                        &
                   calcprog,depth,dgice,dgscrn,dgwater,ice,water,turb)

implicit none

integer, intent(in) :: diag, calcprog
integer iqw, ii
real, intent(in) :: dt
real, dimension(imax), intent(in) :: atm_sg, atm_rg, atm_rnd, atm_snd, atm_u, atm_v
real, dimension(imax), intent(in) :: atm_temp, atm_qg, atm_ps
real, dimension(imax), intent(in) :: atm_vnratio, atm_fbvis, atm_fbnir, atm_inflow, atm_zmin, atm_zmins
type(dgicedata), intent(inout) :: dgice
type(dgscrndata), intent(inout) :: dgscrn
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(inout) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb
real, dimension(imax) :: d_ftop,d_tb,d_zcr
real, dimension(imax) :: d_fb,d_timelt,d_neta,d_ndsn
real, dimension(imax) :: d_ndic,d_nsto,d_delstore
integer, dimension(imax) :: d_nk

if (diag>=1) write(6,*) "Evaluate MLO"

call mlocheck("MLO-start",water_temp=water%temp,water_sal=water%sal,water_u=water%u, &
              water_v=water%v,ice_tsurf=ice%tsurf,ice_temp=ice%temp)

! Set default values for invalid points
do ii = 1,wlev
  where ( depth%dz(:,ii)<1.e-4 )
    water%temp(:,ii) = 288.-wrtemp
    water%sal(:,ii) = 34.72
    water%u(:,ii) = 0.
    water%v(:,ii) = 0.
    turb%k(:,ii)   = omink
    turb%eps(:,ii) = omineps
  end where
end do

! impose limits for rounding errors
if ( mlo_limitsal==1 ) then
  water%sal = min( water%sal, maxsal )
end if  

! adjust levels for free surface
where ( depth%dz(:,1)>1.e-4 )
  d_zcr = max(1.+max(water%eta,-delwater)/depth%depth_hl(:,wlev+1),minwater/depth%depth_hl(:,wlev+1))
elsewhere
  d_zcr = 1.
end where

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
call calcmelt(d_timelt,water)

! equation of state
call getrho(atm_ps,depth,ice,dgwater,water)

! ice mass per unit area
! MJT notes - a limit of minicemass=10 can cause reproducibility issues with
! single precision and multple processes
dgice%imass = max(rhoic*ice%thick+rhosn*ice%snowd, minicemass) 

! split adjustment of free surface and ice thickness to ensure conservation
d_ndsn=ice%snowd
d_ndic=ice%thick
d_nsto=ice%store
d_neta=water%eta

! water fluxes
call fluxcalc(dt,atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,d_neta,     &
              diag,depth,dgwater,ice,water)

! boundary conditions
call getwflux(dt,atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow,  &
              d_zcr,d_neta,depth,dgwater,ice,water)

! ice fluxes
call iceflux(dt,atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v,atm_temp,atm_qg,   &
             atm_ps,atm_zmin,atm_zmins,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_ndsn,d_ndic,d_nsto,d_neta,           &
             d_delstore,diag,dgwater,dgice,ice,water,depth)

! solve for mixed layer depth (calculated at full levels)
call getmixdepth(d_zcr,depth,dgwater,water)

if ( calcprog==0 .or. calcprog==2 .or. calcprog==3 ) then

  ! update ice
  ice%thick = d_ndic
  ice%snowd = d_ndsn
  ice%store = d_nsto
  call mloice(dt,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_neta,diag, &
              depth,dgice,ice,dgwater,water)
  d_ndic = ice%thick
  d_ndsn = ice%snowd
  d_nsto = ice%store
  
end if
if ( calcprog==0 .or. calcprog==2 ) then

  ! update ice velocities due to stress terms
  ice%u = ice%u + dt*(dgice%tauxica-dgice%tauxicw)/dgice%imass
  ice%v = ice%v + dt*(dgice%tauyica-dgice%tauyicw)/dgice%imass
    
end if
#ifndef CCAM
if ( calcprog==0 .or. calcprog==2 .or. calcprog==3 ) then

  ! create or destroy ice
  ! MJT notes - this is done after the flux calculations to agree with the albedo passed to the radiation
  call mlonewice_work(depth,ice,water)

  ! Update exponential time weighted average
  call mlo_ema_uvw(dt,water,depth)
  
end if
#endif
if ( calcprog==0 ) then
  
  ! update water
  call mlocalc(dt,atm_u,atm_v,atm_ps,d_zcr,diag,         &
               depth,dgice,dgwater,ice,water,turb)

end if
if ( calcprog==0 .or. calcprog==2 .or. calcprog==3 ) then

  if ( mlo_adjeta>0 ) then
    ! adjust surface height
    water%eta = d_neta
  end if

end if

! screen diagnostics
call scrncalc(atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,diag, &
              dgice,dgscrn,dgwater,ice,water)

call mlocheck("MLO-end",water_temp=water%temp,water_sal=water%sal,water_u=water%u, &
              water_v=water%v,ice_tsurf=ice%tsurf,ice_temp=ice%temp,               &
              ice_thick=ice%thick)


! energy conservation check
!d_zcr=max(1.+max(water%eta,-delwater)/depth%depth_hl(:,wlev+1),minwater/depth%depth_hl(:,wlev+1))
!dgwater%deleng=0._8
!do ii=1,wlev
!  dgwater%deleng=dgwater%deleng+real((water%temp(:,ii)*d_zcr-oldwatertemp(:,ii)*oldzcr)*depth%dz(:,ii),8)
!end do
!dgwater%deleng=dgwater%deleng*real(wrtrho*cp0/dt,8)
!dgice%deleng=real(ice%store-oldicestore,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%tsurf-oldicefrac*oldicetemp(:,0))*gammi,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,0)*ice%snowd-oldicefrac*oldicetemp(:,1)*oldicesnowd)*cps,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,1)*ice%thick-oldicefrac*oldicetemp(:,2)*oldicethick)*0.5*cpi,8)
!dgice%deleng=dgice%deleng+real((ice%fracice*ice%temp(:,2)*ice%thick-oldicefrac*oldicetemp(:,3)*oldicethick)*0.5*cpi,8)
!dgice%deleng=dgice%deleng/real(dt,8)

return
end subroutine mloeval_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,atm_u,atm_v,atm_ps,d_zcr,diag,         &
                   depth,dgice,dgwater,ice,water,turb)

implicit none

integer, intent(in) :: diag
integer ii, iqw
real, intent(in) :: dt
real, dimension(imax,wlev) :: gammas, rhs, km_hl, ks_hl
real, dimension(imax,2:wlev) :: aa
real, dimension(imax,wlev) :: bb, dd
real, dimension(imax,1:wlev-1) :: cc
real, dimension(imax) :: dumt0, fulldepth, targetdepth, deltaz
real, dimension(imax) :: vmagn, rho, atu, atv
real, dimension(imax), intent(in) :: atm_u, atm_v, atm_ps
real, dimension(imax), intent(inout) :: d_zcr
type(dgicedata), intent(in) :: dgice
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(in) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb

if ( diag>=1 ) write(6,*) "Calculate ocean mixing"

rhs = 0.
aa = 0.
bb = 1.
cc = 0.
dd = 0.

call mlo_calc_k(km_hl,ks_hl,gammas,dt,d_zcr,depth,dgwater,water,turb)

! Counter-gradient term for scalars (rhs)
! +ve sign for rhs terms since z +ve is down
where ( depth%dz(:,1)>1.e-4 )
  rhs(:,1) = ks_hl(:,2)*gammas(:,2)/(depth%dz(:,1)*d_zcr)
end where
do ii = 2,wlev-1
  where ( depth%dz(:,ii)>1.e-4 )  
    rhs(:,ii) = (ks_hl(:,ii+1)*gammas(:,ii+1)-ks_hl(:,ii)*gammas(:,ii))/(depth%dz(:,ii)*d_zcr)
  end where  
end do
where ( depth%dz(:,wlev)>1.e-4 )
  rhs(:,wlev) = -ks_hl(:,wlev)*gammas(:,wlev)/(depth%dz(:,wlev)*d_zcr)
end where


! Diffusion term for scalars (aa,bb,cc)
where ( depth%dz(:,2)*depth%dz(:,1)>1.e-4 )
  cc(:,1) = -dt*ks_hl(:,2)/(depth%dz_hl(:,2)*depth%dz(:,1)*d_zcr**2)
end where  
bb(:,1) = 1. - cc(:,1)
do ii = 2,wlev-1
  where ( depth%dz(:,ii-1)*depth%dz(:,ii)>1.e-4 )  
    aa(:,ii) = -dt*ks_hl(:,ii)/(depth%dz_hl(:,ii)*depth%dz(:,ii)*d_zcr**2)
  end where
  where ( depth%dz(:,ii+1)*depth%dz(:,ii)>1.e-4 )
    cc(:,ii) = -dt*ks_hl(:,ii+1)/(depth%dz_hl(:,ii+1)*depth%dz(:,ii)*d_zcr**2)
  end where  
  bb(:,ii) = 1. - aa(:,ii) - cc(:,ii)
end do
where ( depth%dz(:,wlev-1)*depth%dz(:,wlev)>1.e-4 )
  aa(:,wlev) = -dt*ks_hl(:,wlev)/(depth%dz_hl(:,wlev)*depth%dz(:,wlev)*d_zcr**2)
end where  
bb(:,wlev) = 1. - aa(:,wlev)


! POTENTIAL TEMPERATURE
if ( incradgam>0 ) then
  ! include radiation in counter-gradient term
  do iqw = 1,imax
    dumt0(iqw) = dgwater%wt0(iqw) + sum(dgwater%rad(iqw,1:dgwater%mixind(iqw)))
  end do
else
  dumt0 = dgwater%wt0
end if
do ii = 1,wlev
  dd(:,ii) = water%temp(:,ii) + dt*rhs(:,ii)*dumt0
  where ( depth%dz(:,ii)>1.e-4 )
    dd(:,ii) = dd(:,ii) - dt*dgwater%rad(:,ii)/(depth%dz(:,ii)*d_zcr)
  end where  
end do
where ( depth%dz(:,1)>=1.e-4 )
  dd(:,1) = dd(:,1) - dt*dgwater%wt0/(depth%dz(:,1)*d_zcr)
end where
call thomas(water%temp,aa,bb,cc,dd)


! SALINITY
fulldepth = 0.
targetdepth = min( minsfc, sum( depth%dz, dim=2 )*d_zcr )
do ii = 1,wlev
  deltaz = max( min( depth%dz(:,ii)*d_zcr, targetdepth-fulldepth ), 0.)
  fulldepth = fulldepth + deltaz
  !dd(:,ii) = water%sal(:,ii) - dt*dgwater%ws0_subsurf*deltaz/max(depth%dz(:,ii)*d_zcr*targetdepth,1.e-4)
  !dd(:,ii) = dd(:,ii) + dt*rhs(:,ii)*dgwater%ws0  
  dd(:,ii) = water%sal(:,ii)
  dd(:,ii) = dd(:,ii) + dt*rhs(:,ii)*dgwater%ws0*water%sal(:,1) 
  bb(:,ii) = bb(:,ii) + dt*dgwater%ws0_subsurf*deltaz/max(depth%dz(:,ii)*d_zcr*targetdepth,1.e-4)
end do
where ( depth%dz(:,1)>=1.e-4 )
  !dd(:,1) = dd(:,1) - dt*dgwater%ws0/(depth%dz(:,1)*d_zcr)
  bb(:,1) = bb(:,1) + dt*dgwater%ws0/(depth%dz(:,1)*d_zcr)
end where
call thomas(water%sal,aa,bb,cc,dd)
water%sal = max(0.,water%sal)


! Diffusion term for momentum (aa,bb,cc)
bb(:,1) = 1. - cc(:,1)
do ii = 2,wlev-1
  bb(:,ii) = 1. - aa(:,ii) - cc(:,ii)
end do
bb(:,wlev) = 1. - aa(:,wlev)
where ( depth%dz(:,2)*depth%dz(:,1)>1.e-4 )
  cc(:,1) = -dt*km_hl(:,2)/(depth%dz_hl(:,2)*depth%dz(:,1)*d_zcr**2)
end where
if ( otaumode==1 ) then
  atu = atm_u - fluxwgt*water%u(:,1) - (1.-fluxwgt)*water%utop    ! implicit
  atv = atm_v - fluxwgt*water%v(:,1) - (1.-fluxwgt)*water%vtop    ! implicit
  vmagn = sqrt(max(atu**2+atv**2,1.e-4))                          ! implicit
  rho = atm_ps/(rdry*max(water%temp(:,1)+wrtemp,271.))            ! implicit
  bb(:,1) = 1. - cc(:,1)
  where ( depth%dz(:,1)>=1.e-4 )
    bb(:,1) = bb(:,1) + dt*(1.-ice%fracice)*rho*dgwater%cd        &
                        /(wrtrho*depth%dz(:,1)*d_zcr)             ! implicit
  end where  
else
  bb(:,1) = 1. - cc(:,1)                                          ! explicit  
end if
do ii = 2,wlev-1
  where ( depth%dz(:,ii-1)*depth%dz(:,ii)>1.e-4 )  
    aa(:,ii) = -dt*km_hl(:,ii)/(depth%dz_hl(:,ii)*depth%dz(:,ii)*d_zcr**2)
  end where
  where ( depth%dz(:,ii+1)*depth%dz(:,ii)>1.e-4 )
    cc(:,ii) = -dt*km_hl(:,ii+1)/(depth%dz_hl(:,ii+1)*depth%dz(:,ii)*d_zcr**2)
  end where  
  bb(:,ii) = 1. - aa(:,ii) - cc(:,ii)
end do
where ( depth%dz(:,wlev-1)*depth%dz(:,wlev)>1.e-4 )
  aa(:,wlev) = -dt*km_hl(:,wlev)/(depth%dz_hl(:,wlev)*depth%dz(:,wlev)*d_zcr**2)
end where
bb(:,wlev) = 1. - aa(:,wlev)
! bottom drag
do iqw = 1,imax
  ii = depth%ibot(iqw)  
  if ( depth%dz(iqw,ii)>=1.e-4 ) then
    bb(iqw,ii) = bb(iqw,ii) + dt*dgwater%cd_bot(iqw)/(depth%dz(iqw,ii)*d_zcr(iqw))
  end if
end do


! U diffusion term
do ii = 1,wlev
  dd(:,ii) = water%u(:,ii)
end do
if ( otaumode==1 ) then
  where ( depth%dz(:,1)>=1.e-4 )
    dd(:,1) = dd(:,1) + dt*((1.-ice%fracice)*rho*dgwater%cd*atm_u               &
                        +ice%fracice*dgice%tauxicw)/(wrtrho*depth%dz(:,1)*d_zcr)   ! implicit
  end where
else
  where ( depth%dz(:,1)>=1.e-4 )
    dd(:,1) = dd(:,1) - dt*dgwater%wu0/(depth%dz(:,1)*d_zcr)                       ! explicit
  end where
end if
call thomas(water%u,aa,bb,cc,dd)
if ( otaumode==1 ) then
  dgwater%wu0 = -((1.-ice%fracice)*rho*dgwater%cd*(atm_u-water%u(:,1))  &
          +ice%fracice*dgice%tauxicw)/wrtrho                                     ! implicit
end if

! V diffusion term
do ii = 1,wlev
  dd(:,ii) = water%v(:,ii)
end do
if ( otaumode==1 ) then
  where ( depth%dz(:,1)>=1.e-4 )
    dd(:,1) = dd(:,1) + dt*((1.-ice%fracice)*rho*dgwater%cd*atm_v               &
                        +ice%fracice*dgice%tauyicw)/(wrtrho*depth%dz(:,1)*d_zcr)   ! implicit
  end where
else
  where ( depth%dz(:,1)>=1.e-4 )
    dd(:,1) = dd(:,1) - dt*dgwater%wv0/(depth%dz(:,1)*d_zcr)                      ! explicit
  end where
end if
call thomas(water%v,aa,bb,cc,dd)
if ( otaumode==1 ) then
  dgwater%wv0 = -((1.-ice%fracice)*rho*dgwater%cd*(atm_v-water%v(:,1))  &
          +ice%fracice*dgice%tauyicw)/wrtrho                                     ! implicit
end if


! --- Turn off coriolis terms as this is processed in mlodynamics.f90 ---
!xp = 1. + (0.5*dt*atm_f)**2
!xm = 1. - (0.5*dt*atm_f)**2
!do ii = 1,wlev
!  newa = (water%u(:,ii)*xm+water%v(:,ii)*dt*atm_f)/xp
!  newb = (water%v(:,ii)*xm-water%u(:,ii)*dt*atm_f)/xp
!  water%u(:,ii) = newa
!  water%v(:,ii) = newb
!end do

call mlocheck("MLO-mixing",water_temp=water%temp,water_sal=water%sal,water_u=water%u, &
              water_v=water%v)  
  
return
end subroutine mlocalc

subroutine mlo_calc_k(km_hl,ks_hl,gammas,dt,d_zcr,depth,dgwater,water,turb)                   
                   
implicit none

integer iqw, ii
real, intent(in) :: dt
real, dimension(imax,wlev), intent(inout) :: km_hl, ks_hl, gammas
real, dimension(imax), intent(in) :: d_zcr
type(dgwaterdata), intent(inout) :: dgwater
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb

km_hl = 0.
ks_hl = 0.
gammas = 0.

! solve for stability functions and non-local term (calculated at half levels)
select case(oclosure)
  case(1)
    ! k-e  
    call keps(km_hl,ks_hl,depth,dgwater,water,turb,dt)
  case default
    ! kpp
    call getstab(km_hl,ks_hl,gammas,d_zcr,depth,dgwater,water)
end select

! store currents for next time-step  
do iqw = 1,imax
  water%utop(iqw) = water%u(iqw,1)  
  water%vtop(iqw) = water%v(iqw,1)
  ii =depth%ibot(iqw)
  water%ubot(iqw) = water%u(iqw,ii)
  water%vbot(iqw) = water%v(iqw,ii)
end do  
! Assume ocean mixing occurs after this routine is called
  
return
end subroutine mlo_calc_k
                   
subroutine keps(km_out,ks_out,depth,dgwater,water,turb,dt)

implicit none

real, dimension(imax,wlev), intent(inout) :: km_out, ks_out
type(dgwaterdata), intent(in) :: dgwater
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb
real, intent(in) :: dt

real, dimension(imax,wlev) :: k    !kinetic energy
real, dimension(imax,wlev) :: eps  !dissipation rate
real, dimension(imax,wlev) :: aa  !lower diagnonal
real, dimension(imax,wlev) :: bb  !diagnoal
real, dimension(imax,wlev) :: cc  !upper diagonal
real, dimension(imax,wlev) :: dd  !right hand side
real, dimension(imax,wlev) :: km  !km on full level
real, dimension(imax,wlev) :: ks  !ks on full level
real, dimension(imax,wlev) :: ps  !shear production
real, dimension(imax,wlev) :: pb  !buoyancy production
real, dimension(imax,wlev) :: n2  !Brunt-Vaisala frequency
real, dimension(imax,wlev) :: ce3 !eps Buoyancy coefficient
real, dimension(imax,wlev) :: L   !Length scale
real, dimension(imax,wlev) :: cu  !Galperin stability function
real, dimension(imax,wlev) :: cud !Galperin stability function
real, dimension(imax,wlev) :: alpha !non-dimensional buoyancy parameter

real, dimension(imax,wlev) :: fdepth_hl !fraction depth
real, dimension(imax,wlev) :: d_rho_hl  !d_rho at half level
real, dimension(imax,wlev) :: km_hl     !km on half level
real, dimension(imax,wlev) :: ks_hl     !ks on half level

real, dimension(imax) :: d_ustar
real, dimension(imax) :: pxtr_ema, shear
real, dimension(imax,wlev) :: rho_ema, alpha_ema, beta_ema

real :: dtt, minL
real :: umag, zrough, uoave, voave

integer :: ii,step,iqw,nsteps

real, parameter :: ce1 = 1.44        !eps production coefficient
real, parameter :: ce2 = 1.92        !eps sink coefficient
real, parameter :: ce3stable = -0.4  !eps buoyancy coefficient for stable stratification
real, parameter :: ce3unstable = 1.0 !eps buoyancy coefficient for unstable stratification
real, parameter :: cu0 = 0.5562     
real, parameter :: sigmaeps = 1.08   !eps Schmidt number

d_ustar = max(sqrt(sqrt(dgwater%wu0**2+dgwater%wv0**2)),1.E-6)

fdepth_hl(:,2:wlev) = (depth%depth_hl(:,2:wlev)-depth%depth(:,1:wlev-1))/max(depth%depth(:,2:wlev)-depth%depth(:,1:wlev-1),1.e-8)

! calculate rho_ema from temp_ema and sal_ema
call mlo_ema_ts(dt,water,depth)
pxtr_ema(:) = 0. ! neglect surface pressure as rho is only used for n2 calculation
call calcdensity(rho_ema,alpha_ema,beta_ema,water%temp_ema,water%sal_ema,depth%dz,pxtr_ema)
call interpolate_hl(rho_ema,fdepth_hl,d_rho_hl)

!n2 (full levels)
n2 = 0.
do ii = 2,wlev-1
  where ( depth%dz(:,ii)>1.e-4 )
    n2(:,ii) = -grav/wrtrho*(d_rho_hl(:,ii)-d_rho_hl(:,ii+1))/depth%dz(:,ii)
  end where  
end do
where ( depth%dz(:,1)>1.e-4 )
  ! MJT suggestion
  n2(:,1) = -grav/wrtrho*(rho_ema(:,1)-d_rho_hl(:,2))/(depth%depth_hl(:,2)-depth%depth(:,1))
end where
where ( depth%dz(:,wlev)>1.e-4 )
  ! MJT suggestion
  n2(:,wlev) = -grav/wrtrho*(d_rho_hl(:,wlev-1)-rho_ema(:,wlev))/(depth%depth(:,wlev)-depth%depth_hl(:,wlev-1))
end where

!update arrays based on current state
do ii = 1,wlev
  where ( depth%dz(:,ii)>1.e-4 )  
    k(:,ii) = turb%k(:,ii)
    eps(:,ii) = turb%eps(:,ii)
  elsewhere
    k(:,ii) = omink
    eps(:,ii) = omineps
  end where
end do

!boundary conditions
k(:,1) = (d_ustar(:)/cu0)**2
eps(:,1) = min((cu0)**3*k(:,1)**1.5,1.e9)/(vkar*(0.5*max(depth%dz(:,1),1.e-4)+dgwater%zo(:)))
do iqw = 1,imax
  ii = depth%ibot(iqw)
  uoave = fluxwgt*water%u(iqw,ii) + (1.-fluxwgt)*water%ubot(iqw)
  voave = fluxwgt*water%v(iqw,ii) + (1.-fluxwgt)*water%vbot(iqw)
  umag = sqrt(max(uoave**2+voave**2,1.e-4))
  umag = max( umag, utide )
  zrough = 0.5*max(depth%dz(iqw,ii),1.e-4)/exp(vkar/sqrt(cdbot))
  if ( ii==1 ) then
    k(iqw,1) = 0.5*( k(iqw,1) + cdbot*(umag/cu0)**2 )
    eps(iqw,1) = 0.5*( eps(iqw,1) + min((cu0)**3*max(k(iqw,1),1.e-8)**1.5,1.e9) &
        /(vkar*(0.5*max(depth%dz(iqw,1),1.e-4)+zrough)) )
  else
    k(iqw,ii) = cdbot*(umag/cu0)**2
    eps(iqw,ii) = min((cu0)**3*max(k(iqw,ii),1.e-8)**1.5,1.e9)/(vkar*(0.5*max(depth%dz(iqw,ii),1.e-4)+zrough))
  end if
end do  
k = max( k, omink )
eps = max( eps, omineps )

!limit length scale
L = cu0**3*k**1.5/eps
if ( limitL==1 ) then
  minL = cu0**3*omink**1.5/omineps
  do ii = 2,wlev-1  
    L(:,ii) = max( L(:,ii), minL )
    where ( n2(:,ii) > 0. )
      L(:,ii) = min( L(:,ii), sqrt(0.56*k(:,ii)/n2(:,ii)) )
    end where
  end do
end if
L(:,:) = max( min( L(:,:), omaxl), ominl )

!stability functions
if ( fixedstabfunc==1 ) then
  alpha = 0.
  cu = cu0
  cud = 0.6985
else
  alpha = L**2*n2/k
  alpha = max( min( alpha, 0.56), -0.0466 ) ! SHOC (before eq 6.7.5)
  cu = (cu0 + 2.182*alpha)/(1. + 20.4*alpha + 53.12*alpha**2)
  cud = 0.6985/(1. + 17.34*alpha)
end if

km = max( cu*sqrt(k)*L, 1.e-6 )
ks = max( cud*sqrt(k)*L, 1.e-6 )

!shear production
if ( nops==0 ) then
  do ii=2,wlev-1
    shear = (water%dudz(:,ii)+water%dwdx(:,ii))**2 &
          + (water%dvdz(:,ii)+water%dwdy(:,ii))**2
    ps(:,ii) = km(:,ii)*shear
  end do
else
  do ii=2,wlev-1
    ps(:,ii) = 0.
  end do
end if

!km & ks at half levels
call interpolate_hl(km,fdepth_hl,km_hl)
call interpolate_hl(ks,fdepth_hl,ks_hl)

!coupling loop
nsteps = int(dt/(kemaxdt+0.01)) + 1
dtt = dt/real(nsteps)
do step = 1,nsteps

  !buoyancy production
  if ( nopb==0 ) then
    do ii = 2,wlev-1
      pb(:,ii) = -ks(:,ii)*n2(:,ii)
    end do
  else
    do ii=2,wlev-1
      pb(:,ii) = 0.
    end do
  end if

  !calculate ce3
  if ( fixedce3==1 ) then
    ce3 = ce3stable
  else if ( fixedce3==0 ) then
    do ii=2,wlev-1
      where( pb(:,ii) < 0. )
        ce3(:,ii) = ce3unstable
      else where
        ce3(:,ii) = ce3stable
      endwhere
    end do
  end if

  !solve eps
  !setup diagonals
  aa = 0.
  bb = 1.
  cc = 0.
  dd = eps
  do ii=2,wlev-1
    where ( depth%dz(:,ii)*depth%dz(:,ii-1  )>1.e-4 .and. depth%ibot(:)>ii )  
      aa(:,ii) = -dtt*km_hl(:,ii  )/(depth%dz(:,ii)*depth%dz_hl(:,ii  )*sigmaeps)
    end where
    where ( depth%dz(:,ii)*depth%dz(:,ii+1)>1.e-4 .and. depth%ibot(:)>ii )
      cc(:,ii) = -dtt*km_hl(:,ii+1)/(depth%dz(:,ii)*depth%dz_hl(:,ii+1)*sigmaeps)
    end where  
    bb(:,ii) = 1. - aa(:,ii) - cc(:,ii)
  end do
  if ( eps_mode==0 ) then !explicit eps
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )  
        dd(:,ii) = dd(:,ii) + dtt*eps(:,ii)/k(:,ii)*(ce1*ps(:,ii) + ce3(:,ii)*pb(:,ii) - ce2*eps(:,ii))
      end where
    end do  
  else if ( eps_mode==1 ) then !quasi implicit for eps, Patanker (1980)
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )  
        bb(:,ii) = bb(:,ii) + dtt*ce2*eps(:,ii)/k(:,ii)
        dd(:,ii) = dd(:,ii) + dtt*eps(:,ii)/k(:,ii)*(ce1*ps(:,ii) + ce3(:,ii)*pb(:,ii))
      end where
    end do  
  else if ( eps_mode==2 ) then !quasi implicit for eps & pb, Patanker (1980) & Burchard et al sect 4 (1998)
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. (ce1*ps(:,ii)+ce3(:,ii)*pb(:,ii))>0. .and. depth%ibot(:)>ii )
        bb(:,ii) = bb(:,ii) + dtt*ce2*eps(:,ii)/k(:,ii)
        dd(:,ii) = dd(:,ii) + dtt*eps(:,ii)/k(:,ii)*(ce1*ps(:,ii) + ce3(:,ii)*pb(:,ii))
      elsewhere ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )
        bb(:,ii) = bb(:,ii) + dtt*eps(:,ii)/k(:,ii)*(ce2 - ce3(:,ii)*pb(:,ii)/eps(:,ii))
        dd(:,ii) = dd(:,ii) + dtt*eps(:,ii)/k(:,ii)*(ce1*ps(:,ii))
      end where
    end do  
  end if
  dd(:,2     ) = dd(:,2     ) - aa(:,2     )*eps(:,1)
  dd(:,wlev-1) = dd(:,wlev-1) - cc(:,wlev-1)*eps(:,wlev)

  !solve using thomas algorithm
  call thomas(eps(:,2:wlev-1),aa(:,3:wlev-1),bb(:,2:wlev-1),cc(:,2:wlev-2),dd(:,2:wlev-1))

  !solve k
  !setup diagonals
  aa = 0.
  bb = 1.
  cc = 0.
  dd = k
  do ii = 2,wlev-1
    where ( depth%dz(:,ii)*depth%dz(:,ii-1  )>1.e-4 .and. depth%ibot(:)>ii )  
      aa(:,ii) = -dtt*km_hl(:,ii  )/(depth%dz(:,ii)*depth%dz_hl(:,ii  ))
    end where
    where ( depth%dz(:,ii)*depth%dz(:,ii+1)>1.e-4 .and. depth%ibot(:)>ii )
      cc(:,ii) = -dtt*km_hl(:,ii+1)/(depth%dz(:,ii)*depth%dz_hl(:,ii+1))
    end where
    bb(:,ii) = 1. - aa(:,ii) - cc(:,ii)
  end do
  if ( k_mode==0 ) then !explicit eps
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )  
        dd(:,ii) = dd(:,ii) + dtt*(ps(:,ii) + pb(:,ii) - eps(:,ii))
      end where
    end do    
  else if ( k_mode==1 ) then !quasi impliciit for eps, Patanker (1980)
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )  
        bb(:,ii) = bb(:,ii) + dtt*eps(:,ii)/k(:,ii)
        dd(:,ii) = dd(:,ii) + dtt*(ps(:,ii) + pb(:,ii))
      end where    
    end do    
  else if ( k_mode==2 ) then !quasi implicit for eps & pb, Patanker (1980) & Burchard et al sect 4 (1998)
    do ii = 2,wlev-1
      where ( depth%dz(:,ii)>1.e-4 .and. (ps(:,ii)+pb(:,ii))>0. .and. depth%ibot(:)>ii )  
        bb(:,ii) = bb(:,ii) + dtt*eps(:,ii)/k(:,ii)
        dd(:,ii) = dd(:,ii) + dtt*(ps(:,ii) + pb(:,ii))
      elsewhere ( depth%dz(:,ii)>1.e-4 .and. depth%ibot(:)>ii )
        bb(:,ii) = bb(:,ii) + dtt/k(:,ii)*(eps(:,ii) - pb(:,ii))
        dd(:,ii) = dd(:,ii) + dtt*ps(:,ii)
      end where     
    end do
  end if
  dd(:,2     ) = dd(:,2     ) - aa(:,2     )*k(:,1)
  dd(:,wlev-1) = dd(:,wlev-1) - cc(:,wlev-1)*k(:,wlev)

  !solve using thomas algorithm
  call thomas(k(:,2:wlev-1),aa(:,3:wlev-1),bb(:,2:wlev-1),cc(:,2:wlev-2),dd(:,2:wlev-1))
 
  
  !limit k & eps
  k = max( k, omink )
  eps = max( eps, omineps )

  !limit length scale
  L = cu0**3*k**1.5/eps
  if ( limitL==1 ) then
    minL = cu0**3*omink**1.5/omineps
    do ii = 2,wlev-1
      L(:,ii) = max( L(:,ii), minL )
      where ( n2(:,ii) > 0. )
        L(:,ii) = min( L(:,ii), sqrt(0.56*k(:,ii)/n2(:,ii)) )
      end where
    end do
  end if
  L(:,:) = max( min( L(:,:), omaxl), ominl )

  !stability functions
  if ( fixedstabfunc==1 ) then
    alpha = 0.
    cu = cu0
    cud = 0.6985
  else
    alpha = L**2*n2/k
    alpha = max( min( alpha, 0.56), -0.0466 ) ! SHOC (before eq 6.7.5)
    cu = (cu0 + 2.182*alpha)/(1. + 20.4*alpha + 53.12*alpha**2)
    cud = 0.6985/(1. + 17.34*alpha)
  end if

  km = max( cu*sqrt(k)*L, 1.e-6 )
  ks = max( cud*sqrt(k)*L, 1.e-6 )

  !km & ks at half levels
  call interpolate_hl(km,fdepth_hl,km_hl)
  call interpolate_hl(ks,fdepth_hl,ks_hl)

end do

!update the output variables (internal variables are double precision)
turb%k = k
turb%eps = eps
km_out = km_hl
ks_out = ks_hl

return
end subroutine keps

pure subroutine interpolate_hl(u,fdepth_hl,u_hl)

implicit none

integer k
real, dimension(:,:), intent(in) :: u
real, dimension(:,:), intent(in) :: fdepth_hl
real, dimension(:,:), intent(out) :: u_hl

u_hl(:,1) = 0.
do k = 2,wlev
  u_hl(:,k) = u(:,k-1) + fdepth_hl(:,k)*(u(:,k)-u(:,k-1))
end do

return
end subroutine interpolate_hl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for stability functions
! GFDL use a look-up table to speed-up the code...

pure subroutine getstab(km,ks,gammas,d_zcr,depth,dgwater,water)

implicit none

integer ii,jj,iqw
integer, dimension(imax) :: mixind_hl
real, dimension(imax,wlev), intent(out) :: km,ks,gammas
real, dimension(imax,wlev) :: num,nus,wm,ws,ri,d_nsq
real, dimension(imax) :: sigma,d_depth_hl,d_depth_hlp1
real, dimension(imax) :: a2m,a3m,a2s,a3s
real, dimension(imax) :: numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real, dimension(imax) :: nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real, dimension(imax) :: d_ustar
real xp,cg
real, dimension(imax), intent(in) :: d_zcr
type(dgwaterdata), intent(in) :: dgwater
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.0E-4
real, parameter :: nusw = 0.1E-4

d_ustar = max(sqrt(sqrt(dgwater%wu0**2+dgwater%wv0**2)),1.E-6)
d_nsq = 0.
do ii = 2,wlev
  where ( depth%dz(:,ii-1)*depth%dz(:,ii)>1.e-4 )  
    d_nsq(:,ii) = -grav/wrtrho*(dgwater%rho(:,ii-1)-dgwater%rho(:,ii))/(depth%dz_hl(:,ii)*d_zcr)
  end where  
end do

! stability ---------------------------------------------------------
wm=0.
ws=0.
do ii = 2,wlev
  d_depth_hl=depth%depth_hl(:,ii)*d_zcr
  call getwx(wm(:,ii),ws(:,ii),d_depth_hl,dgwater%bf,d_ustar,dgwater%mixdepth)
end do
wm(:,1)=wm(:,2) ! to avoid problems calculating shallow mixed layer
ws(:,1)=ws(:,2)
!--------------------------------------------------------------------

! Calculate Ri
ri=0. ! ri(:,1) is not used
do ii = 2,wlev
  ri(:,ii)=d_nsq(:,ii)*(depth%dz_hl(:,ii)*d_zcr)**2  & 
           /max((water%u(:,ii-1)-water%u(:,ii))**2         &
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

! calculate G profile -----------------------------------------------
! -ve as z is down
dnumhdz = 0.
dnushdz = 0.
do iqw=1,imax
  mixind_hl(iqw)=dgwater%mixind(iqw)
  jj=min(mixind_hl(iqw)+1,wlev-1)
  d_depth_hl(iqw)=depth%depth_hl(iqw,jj)*d_zcr(iqw)
  if (dgwater%mixdepth(iqw)>d_depth_hl(iqw)) mixind_hl(iqw)=jj
  d_depth_hl(iqw)=depth%depth_hl(iqw,mixind_hl(iqw))*d_zcr(iqw)
  d_depth_hlp1(iqw)=depth%depth_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw)
  xp=(dgwater%mixdepth(iqw)-d_depth_hl(iqw))/max(d_depth_hlp1(iqw)-d_depth_hl(iqw),1.e-8)
  xp=max(0.,min(1.,xp))
  numh(iqw)=(1.-xp)*num(iqw,mixind_hl(iqw))+xp*num(iqw,mixind_hl(iqw)+1)
  nush(iqw)=(1.-xp)*nus(iqw,mixind_hl(iqw))+xp*nus(iqw,mixind_hl(iqw)+1)
  wm1(iqw)=(1.-xp)*wm(iqw,mixind_hl(iqw))+xp*wm(iqw,mixind_hl(iqw)+1)
  ws1(iqw)=(1.-xp)*ws(iqw,mixind_hl(iqw))+xp*ws(iqw,mixind_hl(iqw)+1)
  if ( depth%dz_hl(iqw,mixind_hl(iqw)+1)>1.e-4 ) then
    dnumhdz(iqw)=min(num(iqw,mixind_hl(iqw)+1)-num(iqw,mixind_hl(iqw)),0.) &
        /(depth%dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
    dnushdz(iqw)=min(nus(iqw,mixind_hl(iqw)+1)-nus(iqw,mixind_hl(iqw)),0.) &
      /(depth%dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
  end if  
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
  where (ii<=mixind_hl .and. depth%dz(:,ii)>1.e-4 )
    sigma=depth%depth_hl(:,ii)*d_zcr/dgwater%mixdepth
    km(:,ii)=max(dgwater%mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma)),num(:,ii))
    ks(:,ii)=max(dgwater%mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2s+a3s*sigma)),nus(:,ii))
  elsewhere ( depth%dz(:,ii)>1.e-4 )
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
  where (dgwater%bf<0..and.ii<=mixind_hl.and.depth%dz(:,ii)>1.e-4) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*dgwater%mixdepth,1.E-20)
  end where
end do

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixed layer depth

subroutine getmixdepth(d_zcr,depth,dgwater,water)

implicit none

integer ii,jj,iqw
real vtc,dvsq,vtsq,xp
real, dimension(imax,wlev) :: ws,wm,dumbuoy,rib,d_nsq
real, dimension(imax) :: dumbf,l,d_depth,usf,vsf,rsf,d_b0,d_ustar
real, dimension(imax), intent(in) :: d_zcr
type(dgwaterdata), intent(inout) :: dgwater
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth

d_ustar = max(sqrt(sqrt(dgwater%wu0**2+dgwater%wv0**2)),1.E-6)
vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar**2*ric)
d_b0=-grav*(dgwater%alpha(:,1)*dgwater%wt0-dgwater%beta(:,1)*dgwater%ws0*water%sal(:,1))
d_nsq = 0.
do ii = 2,wlev
  where ( depth%dz(:,ii-1)*depth%dz(:,ii)>1.e-4 )  
    d_nsq(:,ii) = -grav/wrtrho*(dgwater%rho(:,ii-1)-dgwater%rho(:,ii))/(depth%dz_hl(:,ii)*d_zcr)
  end where  
end do

! Modify buoyancy forcing with solar radiation
if (incradbf>0) then
  do ii=1,wlev
    dumbf=d_b0-grav*sum(dgwater%alpha(:,1:ii)*dgwater%rad(:,1:ii),2) ! -ve sign is to account for sign of dgwater%rad
    d_depth=depth%depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
else
  dumbf=d_b0
  do ii=1,wlev
    d_depth=depth%depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
end if

! Estimate surface layer values
usf=water%u(:,1)
vsf=water%v(:,1)
rsf=dgwater%rho(:,1)

! Calculate local buoyancy
dumbuoy=0.
do ii=1,wlev
  dumbuoy(:,ii)=grav*(dgwater%rho(:,ii)-rsf(:))
end do

! Calculate mixed layer depth from critical Ri
dgwater%mixind=wlev-1
dgwater%mixdepth=depth%depth(:,wlev)*d_zcr
rib=0.
do iqw=1,imax
  do ii=1,wlev
    jj=min(ii+1,wlev)
    vtsq=depth%depth(iqw,ii)*d_zcr(iqw)*ws(iqw,ii)*sqrt(0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.))*vtc
    dvsq=(usf(iqw)-water%u(iqw,ii))**2+(vsf(iqw)-water%v(iqw,ii))**2
    rib(iqw,ii)=(depth%depth(iqw,ii)*d_zcr(iqw)-depth%depth(iqw,1))*dumbuoy(iqw,ii) &
        /(max(dvsq+vtsq,1.E-20)*(dgwater%rho(iqw,ii)+wrtrho))
    if (rib(iqw,ii)>ric) then
      jj=max(ii-1,1)
      dgwater%mixind(iqw)=jj
      xp=min(max((ric-rib(iqw,jj))/max(rib(iqw,ii)-rib(iqw,jj),1.E-20),0.),1.)
      dgwater%mixdepth(iqw) = ((1.-xp)*depth%depth(iqw,jj)+xp*depth%depth(iqw,ii))*d_zcr(iqw)
      exit
    end if
  end do 
end do

! calculate buoyancy forcing
call getbf(dgwater,water)

! impose limits for stable conditions
where(dgwater%bf>1.E-10.and.abs(depth%f)>1.E-10)
  l=0.7*d_ustar/abs(depth%f)
elsewhere (dgwater%bf>1.E-10)
  l=d_ustar**3/(vkar*dgwater%bf)
elsewhere
  l=depth%depth(:,wlev)*d_zcr
end where
dgwater%mixdepth=min(dgwater%mixdepth,l)
dgwater%mixdepth=max(dgwater%mixdepth,depth%depth(:,1)*d_zcr)
dgwater%mixdepth=min(dgwater%mixdepth,depth%depth(:,wlev)*d_zcr)

! recalculate index for mixdepth
dgwater%mixind=wlev-1
do iqw=1,imax
  do ii=2,wlev
    if (depth%depth(iqw,ii)*d_zcr(iqw)>dgwater%mixdepth(iqw).or.depth%dz(iqw,ii)<=1.e-4) then
      jj = ii - 1  
      dgwater%mixind(iqw) = jj
      xp=min(max((ric-rib(iqw,jj))/max(rib(iqw,ii)-rib(iqw,jj),1.E-10),0.),1.)
      dgwater%mixdepth(iqw) = ((1.-xp)*depth%depth(iqw,jj)+xp*depth%depth(iqw,ii))*d_zcr(iqw)
      exit
    end if
  end do
end do

! recalculate buoyancy forcing
call getbf(dgwater,water)

return
end subroutine getmixdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate bouyancy forcing

pure subroutine getbf(dgwater,water)

implicit none

integer iqw
real, dimension(imax) :: d_b0
type(dgwaterdata), intent(inout) :: dgwater
type(waterdata), intent(in) :: water

d_b0=-grav*(dgwater%alpha(:,1)*dgwater%wt0-dgwater%beta(:,1)*dgwater%ws0*water%sal(:,1))

if (incradbf>0) then
  do iqw=1,imax
    ! -ve sign is to account for sign of dgwater%rad
    dgwater%bf(iqw)=d_b0(iqw)-grav*sum(dgwater%alpha(iqw,1:dgwater%mixind(iqw)) &
        *dgwater%rad(iqw,1:dgwater%mixind(iqw))) 
  end do
else
  dgwater%bf(:)=d_b0(:)
end if

return
end subroutine getbf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the stability functions

pure subroutine getwx(wm,ws,dep,bf,d_ustar,mixdp)

implicit none

real, dimension(imax), intent(out) :: wm,ws
real, dimension(imax), intent(in) :: bf,mixdp,dep
real, dimension(imax) :: zeta,sig,invl,uuu
real, dimension(imax), intent(in) :: d_ustar
real, parameter :: zetam=-0.2
real, parameter :: zetas=-1.0
real, parameter :: am=1.26
real, parameter :: cm=8.38
real, parameter :: as=-28.86
real, parameter :: cs=98.96

sig=dep/max(mixdp,1.e-8)       ! stable
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
  ws=vkar*(d_ustar**2-16.*zeta/d_ustar)**(1./2.)
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

subroutine getrho(atm_ps,depth,ice,dgwater,water)

implicit none

real, dimension(imax) :: pxtr
real, dimension(imax), intent(in) :: atm_ps
type(icedata), intent(in) :: ice
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth
type(dgwaterdata), intent(inout) :: dgwater

pxtr = atm_ps
if ( usepice==1 ) then
  pxtr = pxtr + grav*ice%fracice*(ice%thick*rhoic+ice%snowd*rhosn)
end if
! neglect eta adjustment
call calcdensity(dgwater%rho,dgwater%alpha,dgwater%beta,water%temp,water%sal,depth%dz,pxtr)

return
end subroutine getrho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate water boundary conditions

subroutine getwflux(dt,atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis, &
                    atm_fbnir,atm_inflow,d_zcr,d_neta,                      &
                    depth,dgwater,ice,water)

implicit none

integer ii
real, dimension(imax) :: visalb,niralb,netvis,netnir
real, dimension(imax), intent(inout) :: d_zcr, d_neta
real, dimension(imax), intent(in) :: atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_inflow
real, intent(in) :: dt
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(in) :: ice
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth

! shortwave
! use -ve as depth is down
visalb=dgwater%visdiralb*atm_fbvis+dgwater%visdifalb*(1.-atm_fbvis)
niralb=dgwater%nirdiralb*atm_fbnir+dgwater%nirdifalb*(1.-atm_fbnir)
netvis=(1.-visalb)*atm_vnratio
netnir=(1.-niralb)*(1.-atm_vnratio)
where ( depth%depth_hl(:,2)<depth%depth_hl(:,wlev+1) )
  dgwater%rad(:,1)=netvis*(exp(-depth%depth_hl(:,2)*d_zcr/mu_1)-1.) &
                  +netnir*(exp(-depth%depth_hl(:,2)*d_zcr/mu_2)-1.)
elsewhere
  dgwater%rad(:,1)=-netvis-netnir ! remainder
end where
do ii=2,wlev-1 
  where ( depth%depth_hl(:,ii+1)<depth%depth_hl(:,wlev+1) )
    dgwater%rad(:,ii)=netvis*(exp(max(-depth%depth_hl(:,ii+1)*d_zcr/mu_1,-50.)) &
                             -exp(max(-depth%depth_hl(:,ii)*d_zcr/mu_1,-50.)))  &
                     +netnir*(exp(max(-depth%depth_hl(:,ii+1)*d_zcr/mu_2,-50.)) &
                             -exp(max(-depth%depth_hl(:,ii)*d_zcr/mu_2,-50.)))
  elsewhere ( depth%dz(:,ii)>1.e-4 )
    dgwater%rad(:,ii)=-netvis*exp(max(-depth%depth_hl(:,ii)*d_zcr/mu_1,-50.)) &
                      -netnir*exp(max(-depth%depth_hl(:,ii)*d_zcr/mu_2,-50.)) ! remainder
  elsewhere
    dgwater%rad(:,ii)=0.  
  end where 
end do
where ( depth%dz(:,wlev)>1.e-4 )
  dgwater%rad(:,wlev)=-netvis*exp(max(-depth%depth_hl(:,wlev)*d_zcr/mu_1,-50.)) &
                      -netnir*exp(max(-depth%depth_hl(:,wlev)*d_zcr/mu_2,-50.)) ! remainder
elsewhere
  dgwater%rad(:,wlev)=0.
end where
do ii=1,wlev
  dgwater%rad(:,ii)=dgwater%rad(:,ii)*(1.-ice%fracice)*atm_sg/(cp0*wrtrho)
end do

! Boundary conditions
! MJT notes - use wrtrho reference density for Boussinesq fluid approximation
dgwater%wu0 = -(1.-ice%fracice)*dgwater%taux/wrtrho
dgwater%wv0 = -(1.-ice%fracice)*dgwater%tauy/wrtrho
dgwater%wt0_eg = -(1.-ice%fracice)*(-dgwater%eg)/(wrtrho*cp0)
dgwater%wt0_rad = -(1.-ice%fracice)*(atm_rg-sbconst*(water%temp(:,1)+wrtemp)**4)/(wrtrho*cp0)
dgwater%wt0_melt = (1.-ice%fracice)*lf*atm_snd/(wrtrho*cp0) ! melting snow
dgwater%wt0 = -(1.-ice%fracice)*(-dgwater%fg)/(wrtrho*cp0) + dgwater%wt0_eg + dgwater%wt0_rad + dgwater%wt0_melt
!dgwater%ws0 = (1.-ice%fracice)*(atm_rnd+atm_snd-dgwater%eg/lv)*water%sal(:,1)/wrtrho
!dgwater%ws0_subsurf = atm_inflow*water%sal(:,1)/wrtrho ! inflow under ice
dgwater%ws0 = (1.-ice%fracice)*(atm_rnd+atm_snd)/wrtrho
! MJT - patch for situations where there is no fresh water inflow
where ( water%sal(:,1)<maxsal ) 
  dgwater%ws0 = dgwater%ws0 - (1.-ice%fracice)*(dgwater%eg/lv)/wrtrho
end where  
dgwater%ws0_subsurf = atm_inflow/wrtrho ! inflow under ice

if ( mlo_adjeta>0 ) then
  d_neta = d_neta + dt*(atm_inflow+(1.-ice%fracice)*(atm_rnd+atm_snd))/wrtrho
end if

return
end subroutine getwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate density

!     compute unesco insitu density (rho in units of gm/cm^3) from
!     salinity (ss in psu), potential temperature (tt in deg C)
!     and pressure (p in bars)
!
!     Reference: Jackett and McDougall, Minimal Adjustment of 
!     Hydrographic Profiles to Achieve Static Stablilty, Journal of
!     Atmospheric and Oceanic Technology, Vol 12, 381-389,  April 1995
                    
subroutine calcdensity(d_rho,d_alpha,d_beta,tt,ss,ddz,pxtr)

implicit none

integer ifx,wlx,ii,iqw
real, dimension(:,:), intent(in) :: tt ! potential temperature
real, dimension(:,:), intent(in) :: ss,ddz
real, dimension(:,:), intent(out) :: d_rho,d_alpha,d_beta
real, dimension(:), intent(in) :: pxtr
real, dimension(size(tt,1)) :: ptot
real t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32,rho0
real drho0dt,drho0ds,dskdt,dskds,sk,sks
real drhodt,drhods,rs0

wlx = size(tt,2)
ifx = size(tt,1)

ptot = pxtr*1.E-5 ! convert Pa to bars

do ii = 1,wlx
  do iqw = 1,ifx
    t = min(max(tt(iqw,ii)+(wrtemp-273.16),-2.2),100.)
    s = min(max(ss(iqw,ii),0.),maxsal)
    p1   = ptot(iqw) + grav*wrtrho*0.5*ddz(iqw,ii)*1.E-5 ! hydrostatic approximation
    ptot(iqw) = ptot(iqw) + grav*wrtrho*ddz(iqw,ii)*1.E-5
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    s2 = s*s
    s3 = s2*s
    p2 = p1*p1
    s32 = sqrt(s3)

    rs0 = (999.842594 - wrtrho) + 6.793952e-2*t                          &
           - 9.095290e-3*t2 + 1.001685e-4*t3                             &
           - 1.120083e-6*t4 + 6.536332e-9*t5 ! density for sal=0.
    rho0 = rs0+ s*(0.824493 - 4.0899e-3*t                                &
           + 7.6438e-5*t2                                                &
           - 8.2467e-7*t3 + 5.3875e-9*t4)                                &
           + s32*(-5.72466e-3 + 1.0227e-4*t                              &
           - 1.6546e-6*t2) + 4.8314e-4*s2    ! + sal terms    
    drho0dt=6.793952e-2                                                  &
           - 2.*9.095290e-3*t + 3.*1.001685e-4*t2                        &
           - 4.*1.120083e-6*t3 + 5.*6.536332e-9*t4                       &
           + s*( -4.0899e-3 + 2.*7.6438e-5*t                             &
           - 3.*8.2467e-7*t2 + 4.*5.3875e-9*t3)                          &
           + s32*(1.0227e-4 - 2.*1.6546e-6*t)
    drho0ds= (0.824493 - 4.0899e-3*t + 7.6438e-5*t2                      &
           - 8.2467e-7*t3 + 5.3875e-9*t4)                                &
           + 1.5*sqrt(s)*(-5.72466e-3 + 1.0227e-4*t                      &
           - 1.6546e-6*t2) + 2.*4.8314e-4*s
    
    sks = 1.965933e4 + 1.444304e2*t - 1.706103*t2                        &
                + 9.648704e-3*t3  - 4.190253e-5*t4                       &
                + p1*(3.186519 + 2.212276e-2*t                           &
                - 2.984642e-4*t2 + 1.956415e-6*t3)                       &
                + p2*(2.102898e-4 - 1.202016e-5*t                        &
                + 1.394680e-7*t2) ! sal=0.
    sk=sks+ s*(52.84855 - 3.101089e-1*t                                  &
                + 6.283263e-3*t2 -5.084188e-5*t3)                        &
                + s32*(3.886640e-1 + 9.085835e-3*t                       &
                - 4.619924e-4*t2)                                        &
                + p1*s*(6.704388e-3  -1.847318e-4*t                      &
                + 2.059331e-7*t2) + 1.480266e-4*p1*s32                   &
                +p2*s*(-2.040237e-6 &
                + 6.128773e-8*t + 6.207323e-10*t2) ! + sal terms             
    dskdt= 1.444304e2 - 2.*1.706103*t                                   &
                + 3.*9.648704e-3*t2  - 4.*4.190253e-5*t3                &
                + s*( - 3.101089e-1                                     &
                + 2.*6.283263e-3*t -3.*5.084188e-5*t2)                  &
                + s32*(9.085835e-3                                      &
                - 2.*4.619924e-4*t)                                     &
                + p1*(2.212276e-2                                       &
                - 2.*2.984642e-4*t + 3.*1.956415e-6*t2)                 &
                + p1*s*(-1.847318e-4                                    &
                + 2.*2.059331e-7*t)                                     &
                + p2*(- 1.202016e-5                                     &
                + 2.*1.394680e-7*t) +p2*s*(                             &
                + 6.128773e-8 + 2.*6.207323e-10*t)
    dskds=(52.84855 - 3.101089e-1*t                                     &
                + 6.283263e-3*t2 -5.084188e-5*t3)                       &
                + 1.5*sqrt(s)*(3.886640e-1 + 9.085835e-3*t              &
                - 4.619924e-4*t2)                                       &
                + p1*(6.704388e-3  -1.847318e-4*t                       &
                + 2.059331e-7*t2) + 1.5*1.480266e-4*p1*sqrt(s)          &
                +p2*(-2.040237e-6                                       &
                + 6.128773e-8*t + 6.207323e-10*t2)
       
    d_rho(iqw,ii)=(rho0*sk + wrtrho*p1)/max(sk-p1,1.)
  
    drhodt=drho0dt*sk/max(sk-p1,1.)-(rho0+wrtrho)*p1*dskdt/(max(sk-p1,1.)**2)
    drhods=drho0ds*sk/max(sk-p1,1.)-(rho0+wrtrho)*p1*dskds/(max(sk-p1,1.)**2)
  
    d_alpha(iqw,ii)=-drhodt              ! Large et al (1993) convention
    d_beta(iqw,ii)=drhods                ! Large et al (1993) convention

  end do  
end do

return
end subroutine calcdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate melting/freezing point
pure subroutine calcmelt(d_timelt,water)

implicit none

real, dimension(imax), intent(out) :: d_timelt
type(waterdata), intent(in) :: water

d_timelt=273.16-0.054*min(max(water%sal(:,1),0.),maxsal) ! ice melting temperature from CICE

return
end subroutine calcmelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert temp to theta based on MOM3 routine
!subroutine gettheta(theta,tt,ss,pxtr,ddz)
!
!implicit none
!
!integer ii
!real, dimension(imax,wlev), intent(out) :: theta
!real, dimension(imax,wlev), intent(in) :: tt,ss,ddz
!real, dimension(imax), intent(in) :: pxtr
!real, dimension(imax) :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11
!real, dimension(imax) :: t,t2,t3,s,s2,p,p2,potmp,ptot
!real, dimension(imax,wlev) :: d_rho
!real, parameter :: density = 1035.
!
!d_rho=density - wrtrho
!
!ptot=pxtr*1.E-5
!do ii=1,wlev
!  t = max(tt(:,ii)+(wrtemp-273.16),-2.)
!  s = max(ss(:,ii),0.)
!  p   = ptot+grav*(d_rho(:,ii)+wrtrho)*0.5*ddz(:,ii)*1.E-5 ! hydrostatic approximation
!  ptot = ptot+grav*(d_rho(:,ii)+wrtrho)*ddz(:,ii)*1.E-5
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
                    d_neta,diag,depth,dgwater,ice,water)

implicit none

integer, intent(in) :: diag
integer it, iqw, ii
real, intent(in) :: dt
real, dimension(imax), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
real, dimension(imax), intent(inout) :: d_neta
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(in) :: ice
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth
real, dimension(imax) :: qsat,dqdt,ri,rho,srcp
real, dimension(imax) :: af, aft, afq, factch, facqch
real, dimension(imax) :: atu, atv
real, dimension(imax) :: vmagn,egmax,d_wavail,dumwatertemp
real ztv, umag, uoave, voave
real consea, dumazmin, afroot, daf, fm, con, dcon, dden
real dfm, den, dcs, sig, root, dumazmins, fh, fq
real est_ustar, est_ustar2, est_dustar2dzo, est_dustardzo
real u10_neutral, est_zo, est_grad
! momentum flux parameters
real, parameter :: charnck = 0.018
real, parameter :: chn10   = 0.00125
! ..... high wind speed - rough sea
real, parameter :: zcom1 = 0.018     ! Charnock's constant
real, parameter :: zcoh1 = 0.0       ! Beljaars 1995 values ( https://doi.org/10.1002/qj.49712152203 )
real, parameter :: zcoq1 = 0.0
! ..... low wind speed - smooth sea
real, parameter :: gnu   = 1.5e-5
real, parameter :: zcom2 = 0.11
real, parameter :: zcoh2 = 0.40
real, parameter :: zcoq2 = 0.62

if (diag>=1.and.ntiles==1) write(6,*) "Calculate ocean fluxes"

! for gfortran
fm = 0.

dumwatertemp=max(water%temp(:,1)+wrtemp,271.)
do iqw = 1,imax
  sig = exp(-grav*max(atm_zmins(iqw),3.)/(rdry*atm_temp(iqw)))
  srcp(iqw) = sig**(rdry/cpair)
end do  
atu=atm_u-fluxwgt*water%u(:,1)-(1.-fluxwgt)*water%utop
atv=atm_v-fluxwgt*water%v(:,1)-(1.-fluxwgt)*water%vtop
vmagn=sqrt(max(atu**2+atv**2,1.e-4))
rho=atm_ps/(rdry*dumwatertemp)
ri=min(grav*(max(atm_zmin,3.)*max(atm_zmin,3.)/max(atm_zmins,3.))*(1.-dumwatertemp*srcp/atm_temp)/vmagn**2,rimax)

call getqsat(qsat,dqdt,dumwatertemp,atm_ps)
if (zomode==0) qsat=0.98*qsat ! with Zeng 1998 for sea water

select case(zomode)
  case(0,1) ! Charnock
    do iqw = 1,imax
      consea=vmagn(iqw)*charnck/grav
      dgwater%zo(iqw)=0.001    ! first guess
      do it=1,4
        dumazmin=max(atm_zmin(iqw),dgwater%zo(iqw)+0.2)
        afroot=vkar/log(dumazmin/dgwater%zo(iqw))
        af(iqw)=afroot**2
        daf=2.*af(iqw)*afroot/(vkar*dgwater%zo(iqw))
        if ( ri(iqw)>=0. ) then ! stable water points                                     
          fm=1./(1.+bprm*ri(iqw))**2
          con=consea*fm*vmagn(iqw)
          dcon=0.
        else                    ! unstable water points
          den=1.+af(iqw)*cms*2.*bprm*sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))
          fm=1.-2.*bprm*ri(iqw)/den
          con=consea*fm*vmagn(iqw)
          dden=daf*cms*2.*bprm*sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))+af(iqw)*cms*bprm*sqrt(-ri(iqw))*dumazmin &
              /(sqrt(dumazmin/dgwater%zo(iqw))*dgwater%zo(iqw)**2)
          dfm=2.*bprm*ri(iqw)*dden/(den**2)
          dcon=consea*dfm*vmagn(iqw)
        end if
        dgwater%zo(iqw)=dgwater%zo(iqw)-(dgwater%zo(iqw)-con*af(iqw)) &
                                       /(1.-dcon*af(iqw)-con*daf)
        dgwater%zo(iqw)=min(max(dgwater%zo(iqw),1.5e-5),6.)
      end do  ! iqw
    end do    ! it=1,4
  case(2) ! Beljaars
    do iqw = 1,imax  
      dgwater%zo(iqw)=0.001    ! first guess
      do it=1,4
        dumazmin=max(atm_zmin(iqw),dgwater%zo(iqw)+0.2)
        afroot=vkar/log(dumazmin/dgwater%zo(iqw))
        af(iqw)=afroot**2
        daf=2.*af(iqw)*afroot/(vkar*dgwater%zo(iqw))
        if ( ri(iqw)>=0. ) then ! stable water points
          fm=1./(1.+bprm*ri(iqw))**2
          consea=zcom1*vmagn(iqw)**2*af(iqw)*fm/grav+zcom2*gnu/(vmagn(iqw)*sqrt(fm*af(iqw)))
          dcs=(zcom1*vmagn(iqw)**2/grav-0.5*zcom2*gnu/(vmagn(iqw)*sqrt(fm*af(iqw))*fm*af(iqw)))*(fm*daf)
        else              ! unstable water points
          con=cms*2.*bprm*sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))
          den=1.+af(iqw)*con
          fm=1.-2.*bprm*ri(iqw)/den
          dfm=2.*bprm*ri(iqw)*(con*daf+af(iqw)*cms*bprm*sqrt(-ri(iqw))*dumazmin &
              /(sqrt(dumazmin/dgwater%zo(iqw))*dgwater%zo(iqw)**2))/(den**2)
          consea=zcom1*vmagn(iqw)**2*af(iqw)*fm/grav+zcom2*gnu/(vmagn(iqw)*sqrt(fm*af(iqw)))
          dcs=(zcom1*vmagn(iqw)**2/grav-0.5*zcom2*gnu/(vmagn(iqw)*sqrt(fm*af(iqw))*fm*af(iqw)))*(fm*daf+dfm*af(iqw))
        end if
        dgwater%zo(iqw)=dgwater%zo(iqw)-(dgwater%zo(iqw)-consea)/(1.-dcs)      
        dgwater%zo(iqw)=min(max(dgwater%zo(iqw),1.5e-5),6.)
      end do  ! iqw
    end do    ! it=1,4
  case(3) ! Moon
    do iqw = 1,imax
      dgwater%zo(iqw)=0.001    ! first guess
      do it=1,4
        dumazmin=max(atm_zmin(iqw),dgwater%zo(iqw)+0.2)
        afroot=vkar/log(dumazmin/dgwater%zo(iqw))
        af(iqw)=afroot**2
        daf=2.*af(iqw)*afroot/(vkar*dgwater%zo(iqw))
        if ( ri(iqw)>=0. ) then ! stable water points                                     
          fm=1./(1.+bprm*ri(iqw))**2
          est_ustar2 = vmagn(iqw)*fm*af(iqw)
          est_dustar2dzo = vmagn(iqw)*fm*daf 
          est_ustar = sqrt(est_ustar2)
          est_dustardzo = sqrt(vmagn(iqw)*fm)*af(iqw)/(vkar*dgwater%zo(iqw))
          u10_neutral = est_ustar*log(10./dgwater%zo(iqw))/vkar
          if ( u10_neutral<=12.5 ) then
            est_zo = (0.018/grav)*est_ustar2
            est_grad = (0.018/grav)*est_dustar2dzo
          else
            est_zo = 0.000085*(-0.56*est_ustar2+20.255*est_ustar+2.458)-0.00058
            est_grad = 0.000085*(-0.56*est_dustar2dzo+20.255*est_dustardzo)
          end if
        else                    ! unstable water points
          den=1.+af(iqw)*cms*2.*bprm*sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))
          fm=1.-2.*bprm*ri(iqw)/den
          dden=daf*cms*2.*bprm*sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))+af(iqw)*cms*bprm*sqrt(-ri(iqw))*dumazmin &
              /(sqrt(dumazmin/dgwater%zo(iqw))*dgwater%zo(iqw)**2)
          dfm=2.*bprm*ri(iqw)*dden/(den**2)
          est_ustar2 = vmagn(iqw)*fm*af(iqw)
          est_dustar2dzo = vmagn(iqw)*(fm*daf+dfm*af(iqw))
          est_ustar = sqrt(est_ustar2)
          est_dustardzo = sqrt(vmagn(iqw)*fm)*af(iqw)/(vkar*dgwater%zo(iqw)) + &
                          sqrt(vmagn(iqw)/fm)*0.5*dfm*afroot
          u10_neutral = est_ustar*log(10./dgwater%zo(iqw))/vkar
          if ( u10_neutral<=12.5 ) then
            est_zo = (0.018/grav)*est_ustar2
            est_grad = (0.018/grav)*est_dustar2dzo
          else
            est_zo = 0.000085*(-0.56*est_ustar2+20.255*est_ustar+2.458)-0.00058
            est_grad = 0.000085*(-0.56*est_dustar2dzo+20.255*est_dustardzo)
          end if          
        end if
        dgwater%zo(iqw)=dgwater%zo(iqw)-(dgwater%zo(iqw)-est_zo)/(1.-est_grad)
        dgwater%zo(iqw)=min(max(dgwater%zo(iqw),1.5e-5),6.)
      end do  ! iqw
    end do    ! it=1,4
end select
af=(vkar/log(max(atm_zmin,dgwater%zo+0.2)/dgwater%zo))**2

select case(zomode)
  case(0,3) ! Charnock CSIRO9 & Moon
    ztv=exp(vkar/sqrt(chn10))/10.
    aft=(vkar/log(max(atm_zmins*ztv,1.)))**2
    afq=aft
    dgwater%zoh=1./ztv
    dgwater%zoq=dgwater%zoh
    factch=sqrt(dgwater%zo/dgwater%zoh)
    facqch=factch
  case(1) ! Charnock zot=zom
    dgwater%zoh=dgwater%zo
    dgwater%zoq=dgwater%zo
    aft=(vkar/(log(max(atm_zmins,dgwater%zo+0.2)/dgwater%zo)))**2
    afq=aft
    factch=1.
    facqch=1.
  case(2) ! Beljaars
    dgwater%zoh=max(zcoh1+zcoh2*gnu/(vmagn*sqrt(fm*af)),1.5E-7)
    dgwater%zoq=max(zcoq1+zcoq2*gnu/(vmagn*sqrt(fm*af)),1.5E-7)
    do iqw = 1,imax
      dumazmins=max(atm_zmins(iqw),dgwater%zo(iqw)+0.2,dgwater%zoh(iqw)+0.2,dgwater%zoq(iqw)+0.2)
      aft(iqw)=vkar**2/(log(dumazmins/dgwater%zo(iqw))*log(dumazmins/dgwater%zoh(iqw)))
      afq(iqw)=vkar**2/(log(dumazmins/dgwater%zo(iqw))*log(dumazmins/dgwater%zoq(iqw)))
    end do  
    factch=sqrt(dgwater%zo/dgwater%zoh)
    facqch=sqrt(dgwater%zo/dgwater%zoq)
end select

! update drag
do iqw = 1,imax
  ! top  
  if (ri(iqw)>=0.) then
    fm=1./(1.+bprm*ri(iqw))**2  ! no zo contrib for stable
    fh=fm
    fq=fm
  else  ! ri is -ve
    dumazmin=max(atm_zmin(iqw),dgwater%zo(iqw)+0.2)
    root=sqrt(-ri(iqw)*dumazmin/dgwater%zo(iqw))
    den=1.+cms*2.*bprm*af(iqw)*root
    fm=1.-2.*bprm*ri(iqw)/den
    dumazmins=max(atm_zmins(iqw),dgwater%zo(iqw)+0.2)
    root=sqrt(-ri(iqw)*dumazmins/dgwater%zo(iqw))
    den=1.+chs*2.*bprm*factch(iqw)*aft(iqw)*root
    fh=1.-2.*bprm*ri(iqw)/den
    den=1.+chs*2.*bprm*facqch(iqw)*afq(iqw)*root
    fq=1.-2.*bprm*ri(iqw)/den
  end if
  dgwater%cd(iqw)=af(iqw)*fm*vmagn(iqw)
  dgwater%cdh(iqw)=aft(iqw)*fh*vmagn(iqw)
  dgwater%cdq(iqw)=afq(iqw)*fq*vmagn(iqw)
  dgwater%umod(iqw)=vmagn(iqw)

  ! bottom
  ii = depth%ibot(iqw)  
  uoave = fluxwgt*water%u(iqw,ii) + (1.-fluxwgt)*water%ubot(iqw)
  voave = fluxwgt*water%v(iqw,ii) + (1.-fluxwgt)*water%vbot(iqw)
  umag = sqrt(max(uoave**2+voave**2,1.e-4))
  umag = max( umag, utide )
  dgwater%cd_bot(iqw) = cdbot*umag
end do

! disable coupling with winds
if ( mlo_uvcoupl==0 ) then
  do iqw = 1,imax
    dgwater%cd(iqw) = 0.
  end do
end if

! turn off lake evaporation when minimum depth is reached
! fg should be replaced with bare ground value
d_wavail=max(depth%depth_hl(:,wlev+1)+d_neta-minwater,0.)
d_wavail = min( d_wavail, max(delwater+d_neta,0.) )
egmax=1000.*lv*d_wavail/(dt*max(1.-ice%fracice,0.01))

! explicit estimate of fluxes
! (replace with implicit scheme if water becomes too shallow)
dgwater%fg=rho*dgwater%cdh*cpair*(dumwatertemp-atm_temp/srcp)
dgwater%fg=min(max(dgwater%fg,-3000.),3000.)
dgwater%eg=min(rho*dgwater%cdq*lv*(qsat-atm_qg),egmax)
dgwater%eg=min(max(dgwater%eg,-3000.),3000.)
dgwater%taux=rho*dgwater%cd*atu
dgwater%tauy=rho*dgwater%cd*atv

if ( mlo_adjeta>0 ) then
  ! update free surface after evaporation
  d_neta=d_neta-0.001*dt*(1.-ice%fracice)*dgwater%eg/lv
end if

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from CCAM)
! following version of getqsat copes better with T < 0C over ice
pure subroutine getqsat(qsat,dqdt,temp,ps)

implicit none

real, dimension(imax), intent(in) :: temp,ps
real, dimension(imax), intent(out) :: qsat,dqdt
real, dimension(imax) :: esatf,tdiff,dedt,rx
integer, dimension(imax) :: ix

real, dimension(0:220), parameter :: table =                                          &
  (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                               & !-146C
     6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                                            & !-141C
     36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                                         & !-136C
     0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,                 & !-131C
     0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,                 & !-126C
     0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,                  & !-121C
     0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,                      & !-116C
     0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,                        & !-111C
     0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,                           & !-106C
     0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,                            & !-101C 
     0.001403, 0.001719, 0.002101, 0.002561, 0.003117, 0.003784,                      & !-95C
     0.004584, 0.005542, 0.006685, 0.008049, 0.009672, 0.01160, 0.01388, 0.01658,     & !-87C
     0.01977, 0.02353, 0.02796, 0.03316, 0.03925, 0.04638, 0.05472, 0.06444, 0.07577, & !-78C
     0.08894, 0.1042, 0.1220, 0.1425, 0.1662, 0.1936, 0.2252, 0.2615, 0.3032,         & !-69C
     0.3511, 0.4060, 0.4688, 0.5406, 0.6225, 0.7159, 0.8223, 0.9432, 1.080,           & !-60C
     1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,                   & !-51C
     3.935, 4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,                          & !-43C
     10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,                   & !-34C
     27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85,                & !-24C 
     77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67,                     & !-16C
     171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78,                 & !-8C
     353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78,                  & !0C
     656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,                  & !8C
     1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,                  & !16C
     1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,                  & !24C
     3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,                  & !32C
     5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,                  & !40C
     7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,                        & !47C
     11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,                   & !54C
     15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,                   & !61C
     21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,                   & !68C
     29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix) + rx*table(ix+1)
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

subroutine mloice(dt,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_neta,diag,  &
                  depth,dgice,ice,dgwater,water)

implicit none

integer, intent(in) :: diag
integer, dimension(imax), intent(inout) :: d_nk
real, intent(in) :: dt
real, dimension(imax), intent(inout) :: d_ftop, d_tb, d_fb
real, dimension(imax), intent(inout) :: d_neta
real, dimension(imax), intent(in) :: d_timelt
type(dgicedata), intent(in) :: dgice
type(icedata), intent(inout) :: ice
type(waterdata), intent(in) :: water
type(depthdata), intent(in) :: depth
type(dgwaterdata), intent(inout) :: dgwater
real, dimension(imax) :: d_salflx, d_wavail

if (diag>=1) write(6,*) "Update ice thermodynamic model"

d_salflx=0.                                               ! fresh water flux
d_wavail=max(depth%depth_hl(:,wlev+1)+d_neta-minwater,0.) ! water avaliable for freezing
d_wavail=min( d_wavail, max(delwater+d_neta,0.)  )

! Ice depth limitation for poor initial conditions
ice%thick=min(ice%thick, icemax)  
! no temperature or salinity conservation

! update ice prognostic variables
call seaicecalc(dt,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_nk,d_wavail,diag, &
                dgice,ice)

! Ice depth limitation for poor initial conditions
ice%thick=min(ice%thick, icemax)  
! no temperature or salinity conservation

!dgwater%ws0_subsurf = dgwater%ws0_subsurf - ice%fracice*d_salflx*water%sal(:,1)/wrtrho
dgwater%ws0_subsurf = dgwater%ws0_subsurf - ice%fracice*d_salflx/wrtrho

if ( mlo_adjeta>0 ) then
  d_neta = d_neta - dt*ice%fracice*d_salflx/wrtrho
end if

return
end subroutine mloice

subroutine mlonewice_work(depth,ice,water)

implicit none

integer iqw, ii
real aa, bb, deldz
real, dimension(imax,wlev) :: sdic
real, dimension(imax) :: newdic, cdic, dsf, delt
real, dimension(imax) :: d_zcr, d_neta
real, dimension(imax) :: d_timelt
type(icedata), intent(inout) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
real, dimension(imax) :: maxnewice, d_wavail
logical, dimension(imax) :: lremove

call calcmelt(d_timelt,water)
d_zcr = max(1.+water%eta/max(depth%depth_hl(:,wlev+1),1.e-4),minwater/max(depth%depth_hl(:,wlev+1),1.e-4))  
d_neta = water%eta

! limits on ice formation due to water avaliability
d_wavail=max(depth%depth_hl(:,wlev+1)+d_neta-minwater,0.)
d_wavail=min( d_wavail, max(delwater+d_neta,0.) )
where ( ice%fracice<0.999 )
  maxnewice=d_wavail*wrtrho/rhoic/(1.-ice%fracice)
elsewhere
  maxnewice=0.
end where

! search for water temperatures that are below freezing
sdic = 0.
newdic = 0.
dsf = 0.
do ii = 1,wlev-1
  do iqw = 1,imax
    if ( maxnewice(iqw)>0. ) then
      aa = depth%dz(iqw,ii)*d_zcr(iqw)
      bb = max(minsfc-dsf(iqw), 0.)
      deldz = min(aa, bb)
      ! Energy = sdic*qice*(1-fracice) = del_temp*cp0*wrtrho*deldz
      sdic(iqw,ii) = max(d_timelt(iqw)-water%temp(iqw,ii)-wrtemp,0.)*cp0*wrtrho*deldz &
                   /qice/(1.-ice%fracice(iqw))
      dsf(iqw) = dsf(iqw) + deldz
    end if
  end do
  newdic(:) = newdic(:) + sdic(:,ii)
end do
newdic=min(newdic,maxnewice)
where ( newdic<=icemin )
  newdic = 0.
end where

d_neta = d_neta - newdic*(1.-ice%fracice)*rhoic/wrtrho

! Adjust temperature in water column to balance the energy cost of ice formation
! Energy = qice*newdic = del_temp*c0*wrtrho*dz*d_zcr
cdic = 0.
do ii = 1,wlev
  sdic(:,ii)=max(min(sdic(:,ii),newdic-cdic),0.)
  cdic=cdic+sdic(:,ii)  
  where ( newdic>icemin .and. depth%dz(:,ii)>1.e-4 )
    water%temp(:,ii)=water%temp(:,ii)+(1.-ice%fracice)*qice*sdic(:,ii) &
        /(cp0*wrtrho*max(depth%dz(:,ii)*d_zcr,1.e-4))
    ! MJT notes - remove salt flux between ice and water for now
    water%sal(:,ii) =water%sal(:,ii)*(1.+(1.-ice%fracice)*sdic(:,ii)*rhoic &
                      /(wrtrho*max(depth%dz(:,ii)*d_zcr,1.e-4)))
  end where
end do

! form new sea-ice
where ( newdic>icemin )
  ice%thick(:)=ice%thick(:)*ice%fracice(:)+newdic(:)*(1.-ice%fracice(:))
  ice%tsurf(:)=ice%tsurf(:)*ice%fracice(:)+d_timelt(:)*(1.-ice%fracice(:))
  ice%store(:)=ice%store(:)*ice%fracice(:)
  ice%snowd(:)=ice%snowd(:)*ice%fracice(:)
  ice%fracice(:)=1.
end where

! Ice depth limitation for poor initial conditions
ice%thick=ice%thick-max(ice%thick-icemax,0.)  
! no temperature or salinity conservation

call mlocheck("MLO-newice",water_temp=water%temp,water_sal=water%sal,ice_tsurf=ice%tsurf, &
              ice_temp=ice%temp,ice_thick=ice%thick)

! removal
lremove = ice%thick<=icemin .and. ice%fracice>0.
where ( lremove )
  newdic=(ice%thick*rhoic+ice%snowd*rhosn)/wrtrho
elsewhere
  newdic=0.
end where
d_neta=d_neta+ice%fracice*newdic*rhoic/wrtrho

! update average temperature and salinity
where ( lremove )
  delt=ice%fracice(:)*ice%thick(:)*qice
  delt=delt+ice%fracice(:)*ice%snowd(:)*qsnow
  delt=delt-ice%fracice(:)*ice%store(:)
  delt=delt-ice%fracice(:)*gammi*(ice%tsurf(:)-d_timelt(:)) ! change from when ice formed
end where

dsf = 0.
do ii = 1,wlev
  do iqw = 1,imax
    if ( lremove(iqw) .and. depth%dz(iqw,ii)>1.e-4 ) then
      ! adjust temperature and salinity in water column
      aa = depth%dz(iqw,ii)*d_zcr(iqw)
      bb = max(minsfc-dsf(iqw),0.)
      deldz = min(aa,bb)
      sdic(iqw,ii) = newdic(iqw)*deldz/minsfc
      water%temp(iqw,ii) = water%temp(iqw,ii)-delt(iqw)*(deldz/minsfc) &
          /(cp0*wrtrho*depth%dz(iqw,ii)*d_zcr(iqw))
      ! MJT notes - remove salt flux between ice and water for now
      water%sal(iqw,ii) = water%sal(iqw,ii)/(1.+ice%fracice(iqw)*sdic(iqw,ii)*rhoic &
                          /(wrtrho*depth%dz(iqw,ii)*d_zcr(iqw)))
      dsf(iqw) = dsf(iqw) + deldz
    end if
  end do
end do  

! remove ice
where ( lremove )
  ice%fracice(:)=0.
  ice%thick(:)=0.
  ice%snowd(:)=0.
  ice%store(:)=0.
  ice%tsurf(:)=273.16
  ice%temp(:,0)=273.16
  ice%temp(:,1)=273.16
  ice%temp(:,2)=273.16
end where

if ( mlo_adjeta>0 ) then
  ! MJT notes - remove salt flux between ice and water for now
  water%eta = d_neta
end if

call mlocheck("MLO-icemelt",water_temp=water%temp,water_sal=water%sal,ice_tsurf=ice%tsurf, &
              ice_temp=ice%temp,ice_thick=ice%thick)

return
end subroutine mlonewice_work
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(dt,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_nk,    &
                      d_wavail,diag,                                 &
                      dgice,ice)

implicit none

integer, intent(in) :: diag
integer pc, nc
integer, dimension(imax), intent(inout) :: d_nk
integer, dimension(imax) :: dt_nk
real, intent(in) :: dt
real, dimension(imax), intent(inout) :: d_ftop,d_tb,d_fb,d_salflx
real, dimension(imax), intent(inout) :: d_wavail
real, dimension(imax), intent(in) :: d_timelt
type(dgicedata), intent(in) :: dgice
type(icedata), intent(inout) :: ice
real, dimension(imax) :: it_tn0,it_tn1,it_tn2
real, dimension(imax) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(imax) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
real, dimension(imax) :: pt_egice
logical, dimension(imax,5) :: pqpack

if (diag>=1) write(6,*) "Pack ice data"

! Pack different ice configurations
pqpack(:,1)=( d_nk==2 .and. ice%snowd>0.05 .and. ice%thick>icemin )   ! thick snow + 2 ice layers
pqpack(:,2)=( d_nk==1 .and. ice%snowd>0.05 .and. ice%thick>icemin )   ! thick snow + 1 ice layer
pqpack(:,3)=( d_nk==2 .and. ice%snowd<=0.05 .and. ice%thick>icemin )  ! 2 ice layers (+ snow?)
pqpack(:,4)=( d_nk==1 .and. ice%snowd<=0.05 .and. ice%thick>icemin )  ! 1 ice layer (+ snow?)
pqpack(:,5)=( d_nk==0 .and. ice%thick>icemin )                        ! thin ice (+ snow?)

! Update ice and snow temperatures
do pc=1,5
  nc=count(pqpack(:,pc))  
  if ( nc>0 ) then
    it_tsurf(1:nc)  =pack(ice%tsurf,pqpack(:,pc))
    it_tn0(1:nc)    =pack(ice%temp(:,0),pqpack(:,pc))
    it_tn1(1:nc)    =pack(ice%temp(:,1),pqpack(:,pc))
    it_tn2(1:nc)    =pack(ice%temp(:,2),pqpack(:,pc))
    pt_egice(1:nc)  =pack(dgice%eg,pqpack(:,pc))
    it_dic(1:nc)    =pack(ice%thick,pqpack(:,pc))
    it_dsn(1:nc)    =pack(ice%snowd,pqpack(:,pc))
    it_sto(1:nc)    =pack(ice%store,pqpack(:,pc))
    dt_ftop(1:nc)   =pack(d_ftop,pqpack(:,pc))
    dt_tb(1:nc)     =pack(d_tb,pqpack(:,pc))
    dt_fb(1:nc)     =pack(d_fb,pqpack(:,pc))
    dt_timelt(1:nc) =pack(d_timelt,pqpack(:,pc))
    dt_salflx(1:nc) =pack(d_salflx,pqpack(:,pc))
    dt_wavail(1:nc) =pack(d_wavail,pqpack(:,pc))
    dt_nk(1:nc)     =pack(d_nk,pqpack(:,pc))
    select case(pc)
      case(1)
        call icetemps1s2i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,            &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,          &
                 dt_nk,dt_wavail,pt_egice,diag)
      case(2)
        call icetemps1s1i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,            &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,          &
                 dt_nk,dt_wavail,pt_egice,diag)
      case(3)
        call icetempi2i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,              &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,          &
                 dt_nk,dt_wavail,pt_egice,diag)
      case(4)
        call icetempi1i(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,              &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,          &
                 dt_nk,dt_wavail,pt_egice,diag)
      case(5)
        call icetemps(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,                &
                 it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,          &
                 dt_nk,dt_wavail,pt_egice,diag)
    end select
    ice%tsurf    =unpack(it_tsurf(1:nc),pqpack(:,pc),ice%tsurf)
    ice%temp(:,0)=unpack(it_tn0(1:nc),pqpack(:,pc),ice%temp(:,0))
    ice%temp(:,1)=unpack(it_tn1(1:nc),pqpack(:,pc),ice%temp(:,1))
    ice%temp(:,2)=unpack(it_tn2(1:nc),pqpack(:,pc),ice%temp(:,2))
    ice%thick    =unpack(it_dic(1:nc),pqpack(:,pc),ice%thick)
    ice%snowd    =unpack(it_dsn(1:nc),pqpack(:,pc),ice%snowd)
    ice%store    =unpack(it_sto(1:nc),pqpack(:,pc),ice%store)
    d_salflx     =unpack(dt_salflx(1:nc),pqpack(:,pc),d_salflx)
    d_nk         =unpack(dt_nk(1:nc),pqpack(:,pc),d_nk)
    select case(pc)
      case(1)
        call mlocheck("MLO-icetemps1s2i",ice_tsurf=ice%tsurf,ice_temp=ice%temp, &
                      ice_thick=ice%thick)
      case(2)
        call mlocheck("MLO-icetemps1s1i",ice_tsurf=ice%tsurf,ice_temp=ice%temp, &
                      ice_thick=ice%thick)
      case(3)
        call mlocheck("MLO-icetempi2i",ice_tsurf=ice%tsurf,ice_temp=ice%temp,   &
                      ice_thick=ice%thick)
      case(4)
        call mlocheck("MLO-icetempi1i",ice_tsurf=ice%tsurf,ice_temp=ice%temp,   &
                      ice_thick=ice%thick)  
      case(5)
        call mlocheck("MLO-icetemps",ice_tsurf=ice%tsurf,ice_temp=ice%temp,     &
                      ice_thick=ice%thick)
    end select
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
                     dt_timelt,dt_salflx,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real htdown
real, dimension(nc) :: con,rhin,rhsn,conb,qmax,fl,qneed,xxx,excess,conc
real, dimension(nc) :: subl,snmelt,dhb,isubl,ssubl,smax
real, dimension(nc) :: simelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice
real, dimension(nc,2:4) :: aa
real, dimension(nc,4) :: bb,dd
real, dimension(nc,3) :: cc
real, dimension(nc,4) :: ans

if (diag>=1) write(6,*) "Two ice layers + snow"

! Thickness of each layer
rhsn=1./max(it_dsn,0.05)
rhin=2./(max(it_dic,himin)-gammi/cpi)
con =condsnw/max(it_dsn,0.05)
conb=1./(it_dsn/condsnw+0.5*max(it_dic,himin)/condice)
conc=2.*condice/max(it_dic,himin)

! no need to map from generic ice pack into 2 layer + snow

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gammi
cc(:,1)=-dt*2.*con/gammi
dd(:,1)=it_tsurf+dt*dt_ftop/gammi
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
! fix up excess flux
fl=2.*conc*(dt_tb-it_tn2)
dhb=dt*(fl-dt_fb)/qice                ! Excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*wrtrho/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+(flnew-fl)/(cpi*max(it_dic,himin)) ! Modify temperature if limit is reached
it_tn2=it_tn2+(flnew-fl)/(cpi*max(it_dic,himin)) ! Modify temperature if limit is reached

! Bottom ablation or accretion
it_dic = it_dic + dhb
dt_salflx = dt_salflx + dhb*rhoic/dt

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*ls/(lv*qsnow)        ! net snow+ice sublimation
ssubl=min(subl,it_dsn)                ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic        ! ice only component of sublimation
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! use stored heat in brine pockets to keep temperature at -0.1 until heat is used up
qneed=(dt_timelt-it_tn1)*0.5*cpi*it_dic ! J/m**2
qneed=max(min(qneed,it_sto),0.)
where ( it_dic>himin )
  it_tn1=it_tn1+2.*qneed/(cpi*it_dic)
  it_sto=it_sto-qneed
end where

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
it_dic=it_dic-smax
dt_salflx=dt_salflx-smax*rhoic/dt

! the following are snow to ice processes
xxx=it_dic+it_dsn-(rhosn*it_dsn+rhoic*it_dic)/wrtrho ! white ice formation
excess=max(it_dsn-xxx,0.)*rhosn/wrtrho               ! white ice formation
excess=min( excess, it_dsn*rhosn/wrtrho )
it_dsn = it_dsn - excess*wrtrho/rhosn
it_dic = it_dic + excess*wrtrho/rhoic

! Snow melt
snmelt = max(it_tsurf-273.16,0.)*gammi/qsnow
snmelt = min(snmelt,it_dsn)
it_dsn = it_dsn - snmelt
it_tsurf = it_tsurf - snmelt*qsnow/gammi
dt_salflx = dt_salflx - snmelt*rhosn/dt   ! melt fresh water snow (no salt when melting snow)

! Ice melt
simelt = max(it_tsurf-dt_timelt,0.)*gammi/qice
simelt = min(simelt,it_dic)
it_dic = it_dic - simelt
it_tsurf = it_tsurf - simelt*qice/gammi
dt_salflx = dt_salflx - simelt*rhoic/dt

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
! Calculate ice temperature for one layer and snow (from Mk3 seaice.f)

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
                     dt_timelt,dt_salflx,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real htup,htdown
real, dimension(nc) :: con,rhin,rhsn,conb,qmax,sbrine,fl,xxx,excess,conc
real, dimension(nc) :: subl,snmelt,dhb,isubl,ssubl,smax
real, dimension(nc) :: simelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice
real, dimension(nc,2:3) :: aa
real, dimension(nc,3) :: bb,dd
real, dimension(nc,2) :: cc
real, dimension(nc,3) :: ans

if (diag>=1) write(6,*) "One ice layer + snow"

! Thickness of each layer
rhsn=1./max(it_dsn,0.05)
rhin=1./(max(it_dic,himin)-gammi/cpi)
con =condsnw/max(it_dsn,0.05)
conb=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)
conc=condice/max(it_dic,icemin)

! map from generic ice pack to 1 layer + snow
!it_tsurf=it_tsurf
!it_tn0=it_tn0
it_tn1=0.5*(it_tn1+it_tn2)

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gammi
cc(:,1)=-dt*2.*con/gammi
dd(:,1)=it_tsurf+dt*dt_ftop/gammi
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
! fix excess flux
fl=2.*conc*(dt_tb-it_tn1)
dhb=dt*(fl-dt_fb)/qice                      ! Excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*wrtrho/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+(flnew-fl)/(cpi*max(it_dic,himin)-gammi) ! modify temperature if limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflx=dt_salflx+dhb*rhoic/dt

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*ls/(lv*qsnow)        ! net snow+ice sublimation
ssubl=min(subl,it_dsn)                ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic        ! ice only component of sublimation
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! remove brine store for single layer
sbrine=min(it_sto/qice,it_dic)
it_sto=it_sto-sbrine*qice
it_dic=it_dic-sbrine
dt_salflx=dt_salflx-sbrine*rhoic/dt

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
it_dic=it_dic-smax
dt_salflx=dt_salflx-smax*rhoic/dt

! the following are snow to ice processes
xxx=it_dic+it_dsn-(rhosn*it_dsn+rhoic*it_dic)/wrtrho ! white ice formation
excess=max(it_dsn-xxx,0.)*rhosn/wrtrho               ! white ice formation
excess=min( excess, it_dsn*rhosn/wrtrho )
it_dsn=it_dsn-excess*wrtrho/rhosn
it_dic=it_dic+excess*wrtrho/rhoic

! Snow melt
snmelt=max(it_tsurf-273.16,0.)*gammi/qsnow
snmelt=min(snmelt,it_dsn)
it_dsn=it_dsn-snmelt
it_tsurf=it_tsurf-snmelt*qsnow/gammi
dt_salflx=dt_salflx-snmelt*rhosn/dt        ! melt fresh water snow (no salt when melting snow)

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gammi/qice
simelt=min(simelt,it_dic)
it_dic=it_dic-simelt
it_tsurf=it_tsurf-simelt*qice/gammi
dt_salflx=dt_salflx-simelt*rhoic/dt

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

it_tn2 = it_tn1

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
                     dt_timelt,dt_salflx,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer, dimension(nc), intent(inout) :: dt_nk
real, intent(in) :: dt
real htdown
real, dimension(nc) :: rhin,qmax,qneed,fl,con,gamm,ssubl,isubl,conb
real, dimension(nc) :: subl,simelt,smax,dhb,snmelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
real, dimension(nc), intent(in) :: pt_egice
real, dimension(nc,2:3) :: aa
real, dimension(nc,1:3) :: bb,dd
real, dimension(nc,1:2) :: cc
real, dimension(nc,3) :: ans

if (diag>=1) write(6,*) "Two ice layers + without snow"

con =1./(it_dsn/condsnw+0.5*max(it_dic,himin)/condice)
conb=2.*condice/max(it_dic,himin)
gamm=cps*it_dsn+gammi  ! for energy conservation

! map from generic ice pack to 2 layer
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0)/gamm
!it_tn1=it_tn1
!it_tn2=it_tn2

! Thickness of each layer
rhin=2./(max(it_dic,himin)-gammi/cpi)

! Solve implicit ice temperature matrix
bb(:,1)=1.+dt*2.*con/gamm
cc(:,1)=-dt*2.*con/gamm
dd(:,1)=it_tsurf+dt*dt_ftop/gamm
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
! fix excess flux
fl=2.*conb*(dt_tb-it_tn2)
dhb=dt*(fl-dt_fb)/qice                   ! first guess of excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*wrtrho/rhoic)
flnew=dt_fb+dhb*qice/dt
it_tn1=it_tn1+2.*(flnew-fl)/(cpi*max(it_dic,himin)) ! modify temperature if limit is reached
it_tn2=it_tn2+2.*(flnew-fl)/(cpi*max(it_dic,himin)) ! modify temperature if limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflx=dt_salflx+dhb*rhoic/dt

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*ls/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl
gamm=gammi+cps*it_dsn

! update brine store for middle layer
qneed=(dt_timelt-it_tn1)*0.5*cpi*it_dic ! J/m**2
qneed=max(min(qneed,it_sto),0.)
where ( it_dic>himin )
  it_tn1=it_tn1+2.*qneed/(cpi*it_dic)
  it_sto=it_sto-qneed
end where

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=it_sto-smax*qice
it_dic=it_dic-smax
dt_salflx=dt_salflx-smax*rhoic/dt

! Snow melt
snmelt=max(it_tsurf-273.16,0.)*gamm/qsnow
snmelt=min(snmelt,it_dsn)
it_dsn=it_dsn-snmelt
it_tsurf = it_tsurf - snmelt*qsnow/gamm
dt_salflx=dt_salflx-snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)
gamm=gammi+cps*it_dsn

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
simelt=min(simelt,it_dic)
it_dic=it_dic-simelt
it_tsurf = it_tsurf - simelt*qice/gamm
dt_salflx=dt_salflx-simelt*rhoic/dt

! test whether to change number of layers
htdown=2.*himin
where ( it_dic<htdown )
  ! code to decrease number of layers
  dt_nk=1
  ! merge ice layers into one
  it_tn1=0.5*(it_tn1+it_tn2)
  it_tn2=it_tn1
end where

it_tn0 = it_tsurf

return
end subroutine icetempi2i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for one layer and no snow

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
                     dt_timelt,dt_salflx,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer, dimension(nc), intent(inout) :: dt_nk
real, intent(in) :: dt
real htup,htdown
real, dimension(nc) :: rhin,sbrine,fl,con,gamm,ssubl,isubl,conb
real, dimension(nc) :: subl,simelt,dhb,snmelt,flnew
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
real, dimension(nc), intent(in) :: pt_egice
real, dimension(nc,2:2) :: aa
real, dimension(nc,1:2) :: bb,dd
real, dimension(nc,1:1) :: cc
real, dimension(nc,2) :: ans

if (diag>=1) write(6,*) "One ice layer + without snow"

con =1./(it_dsn/condsnw+max(it_dic,icemin)/condice)        ! snow/ice conduction
conb=condice/max(it_dic,icemin)                            ! ice conduction
gamm=gammi+cps*it_dsn                                      ! for energy conservation

! map from generic ice pack to 1 layer
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0)/gamm
it_tn1=0.5*(it_tn1+it_tn2)

! Thickness of layer
rhin=1./(max(it_dic,himin)-gammi/cpi)

! Solve implicit ice temperature matrix
! f0=2.*con*(newt1-newtsurf) is the flux between tsurf and t1
! fl=2.*conb*(dt_tb-newt1) is the flux between t1 and the bottom
bb(:,1)=1.+dt*2.*con/gamm
cc(:,1)=-dt*2.*con/gamm
dd(:,1)=it_tsurf+dt*dt_ftop/gamm   
aa(:,2)=-dt*2.*con*rhin/cpi      
bb(:,2)=1.+dt*2.*con*rhin/cpi+dt*2.*conb*rhin/cpi
dd(:,2)=it_tn1+dt*2.*conb*dt_tb*rhin/cpi  
call thomas(ans,aa,bb,cc,dd)
it_tsurf(:)=ans(:,1)
it_tn1(:)  =ans(:,2)
! fix excessive flux
fl=2.*conb*(dt_tb-it_tn1)                   ! flux between t1 and bottom
dhb=dt*(fl-dt_fb)/qice                      ! first guess of excess flux between water and ice layer
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*wrtrho/rhoic)
flnew=dt_fb+dhb*qice/dt                     ! final excess flux from below
it_tn1=it_tn1+(flnew-fl)/(cpi*max(it_dic,himin))  ! does nothing unless a limit is reached

! Bottom ablation or accretion
it_dic=it_dic+dhb
dt_salflx=dt_salflx+dhb*rhoic/dt

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*ls/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl
gamm=gammi+cps*it_dsn

! remove brine store for single layer
sbrine=min(it_sto/qice,it_dic)
it_sto=it_sto-sbrine*qice
it_dic=it_dic-sbrine
dt_salflx=dt_salflx-sbrine*rhoic/dt

! Snow melt
snmelt=max(it_tsurf-273.16,0.)*gamm/qsnow
snmelt=min(snmelt,it_dsn)
it_dsn=it_dsn-snmelt
it_tsurf=it_tsurf-snmelt*qsnow/gamm
dt_salflx=dt_salflx-snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)
gamm=gammi+cps*it_dsn

! Ice melt
simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
simelt=min(simelt,it_dic)
it_dic=it_dic-simelt
it_tsurf=it_tsurf-simelt*qice/gamm
dt_salflx=dt_salflx-simelt*rhoic/dt
  
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
                     dt_timelt,dt_salflx,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,f0,tnew
real, dimension(nc) :: subl,snmelt,sbrine,dhb,isubl,ssubl,gamm,simelt
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

if (diag>=1) write(6,*) "Combined (thin) ice and snow layer"

con=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)         ! conductivity
gamm=cps*it_dsn+gammi+cpi*it_dic                           ! heat capacity

! map from generic ice pack to thin ice
it_tsurf=(gammi*it_tsurf+cps*it_dsn*it_tn0+0.5*cpi*it_dic*(it_tn1+it_tn2))/gamm

! Update tsurf based on fluxes from above and below
! MJT notes - we need an implicit form of temperature (tnew) to account for
! fluxes
tnew=(it_tsurf*gamm/dt+dt_ftop+con*dt_tb)/(gamm/dt+con)    ! predictor temperature for flux calculation
! fix excessive flux
f0=con*(dt_tb-tnew)                                                       ! first guess of flux from below
dhb=dt*(f0-dt_fb)/qice                                                    ! excess flux converted to change in ice thickness
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*wrtrho/rhoic)
f0=dhb*qice/dt+dt_fb                                                      ! final flux from below
it_tsurf=it_tsurf+dt*(dt_ftop+f0)/gamm                     ! update ice/snow temperature

! Bottom ablation or accretion to balance energy budget
it_dic=it_dic+dhb
dt_salflx=dt_salflx+dhb*rhoic/dt

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*ls/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
it_dsn=it_dsn-ssubl
it_dic=it_dic-isubl

! remove brine store for single layer
sbrine=min(it_sto/qice,it_dic)
it_sto=it_sto-sbrine*qice
it_dic=it_dic-sbrine
dt_salflx=dt_salflx-sbrine*rhoic/dt

! Snow melt
gamm = gammi + cps*it_dsn + cpi*it_dic
snmelt = max(it_tsurf-273.16,0.)*gamm/qsnow
snmelt = min(snmelt,it_dsn)
it_dsn = it_dsn - snmelt
it_tsurf = it_tsurf - snmelt*qsnow/gamm
dt_salflx = dt_salflx - snmelt*rhosn/dt ! melt fresh water snow (no salt when melting snow)

! Ice melt
gamm=gammi+cps*it_dsn+cpi*it_dic
simelt=max(it_tsurf-dt_timelt,0.)*gamm/qice
simelt=min(simelt,it_dic)
it_dic=it_dic-simelt
it_tsurf=it_tsurf-simelt*qice/gamm
dt_salflx=dt_salflx-simelt*rhoic/dt

! Recalculate thickness index
where ( it_dic>=himin )
  dt_nk=1
end where

it_tn0=it_tsurf
it_tn1=it_tsurf
it_tn2=it_tsurf

return
end subroutine icetemps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solves matrix for temperature conduction
                     
pure subroutine thomas(outo,aai,bbi,cci,ddi)

implicit none

integer ii, nlev
real, dimension(:,:), intent(out) :: outo
real, dimension(:,2:), intent(in) :: aai
real, dimension(:,:), intent(in) :: bbi,ddi
real, dimension(:,:), intent(in) :: cci
real, dimension(size(outo,1),size(outo,2)) :: cc, dd
real, dimension(size(outo,1)) :: n

nlev = size(outo,2)
cc(:,1) = cci(:,1)/bbi(:,1)
dd(:,1) = ddi(:,1)/bbi(:,1)

do ii = 2,nlev-1
  n = bbi(:,ii)-cc(:,ii-1)*aai(:,ii)
  cc(:,ii) = cci(:,ii)/n
  dd(:,ii) = (ddi(:,ii)-dd(:,ii-1)*aai(:,ii))/n
end do
n = bbi(:,nlev)-cc(:,nlev-1)*aai(:,nlev)
dd(:,nlev) = (ddi(:,nlev)-dd(:,nlev-1)*aai(:,nlev))/n
outo(:,nlev) = dd(:,nlev)
do ii = nlev-1,1,-1
  outo(:,ii) = dd(:,ii)-cc(:,ii)*outo(:,ii+1)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dt,atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v,atm_temp,atm_qg, &
                   atm_ps,atm_zmin,atm_zmins,d_ftop,d_tb,d_fb,d_timelt,d_nk,                                     &
                   d_ndsn,d_ndic,d_nsto,d_neta,d_delstore,diag,                                                  &
                   dgwater,dgice,ice,water,depth)

implicit none

integer, intent(in) :: diag
integer itr
real, dimension(imax), intent(in) :: atm_sg,atm_rg,atm_rnd,atm_snd,atm_vnratio,atm_fbvis,atm_fbnir,atm_u,atm_v
real, dimension(imax), intent(in) :: atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
real, dimension(imax), intent(inout) :: d_ftop,d_tb,d_fb,d_ndsn,d_ndic
real, dimension(imax), intent(inout) :: d_nsto,d_delstore,d_neta
real, dimension(imax), intent(in) :: d_timelt
integer, dimension(imax), intent(inout) :: d_nk
type(dgwaterdata), intent(inout) :: dgwater
type(dgicedata), intent(inout) :: dgice
type(icedata), intent(inout) :: ice
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth
real, intent(in) :: dt
real, dimension(imax) :: qsat,dqdt,ri,rho,srcp,tnew,qsatnew,gamm,bot
real, dimension(imax) :: fm,fh,fq,af,aft,afq
real, dimension(imax) :: den,sig,root
real, dimension(imax) :: alb,eye
real, dimension(imax) :: uu,vv,du,dv,vmagn,icemagn
real, dimension(imax) :: ustar,g,h,dgu,dgv,dhu,dhv,det
real, dimension(imax) :: newiu,newiv,dtsurf
real, dimension(imax) :: dumazmin, dumazmins
real, dimension(imax) :: wavail, f0
real zohseaice, zoqseaice

if (diag>=1) write(6,*) "Calculate ice fluxes"

! Prevent unrealistic fluxes due to poor input surface temperature
dtsurf=min(ice%tsurf,273.2)
uu=atm_u-ice%u
vv=atm_v-ice%v
vmagn=sqrt(max(uu*uu+vv*vv,1.E-4))
sig=exp(-grav*atm_zmins/(rdry*atm_temp))
srcp=sig**(rdry/cpair)
rho=atm_ps/(rdry*dtsurf)

zohseaice=zoseaice/(factchseaice**2)
zoqseaice=zoseaice/(factchseaice**2)
dumazmin=max(atm_zmin,zoseaice+0.2)
dumazmins=max(atm_zmin,zoseaice+0.2)
af=vkar*vkar/(log(dumazmin/zoseaice)*log(dumazmin/zoseaice))
aft=vkar*vkar/(log(dumazmins/zoseaice)*log(dumazmins/zohseaice))
afq=aft

! number of (thick) ice layers
d_nk=min(int(d_ndic/himin),2)

! radiation
alb=     atm_vnratio*(dgice%visdiralb*atm_fbvis+dgice%visdifalb*(1.-atm_fbvis))+ &
    (1.-atm_vnratio)*(dgice%visdifalb*atm_fbnir+dgice%visdifalb*(1.-atm_fbnir))
! allow penetrating radiation to be stored up to some maximum value
where ( d_ndic>icemin )
  eye=0.35  
  d_delstore=atm_sg*(1.-alb)*eye
  d_nsto=d_nsto+dt*d_delstore
elsewhere
  eye=0.  
  d_delstore=0.
end where

where ( d_nk>0 .and. d_ndsn>0.05 )
  gamm=gammi
elsewhere
  gamm=cps*d_ndsn+gammi
end where

! water temperature at bottom of ice
d_tb=max(water%temp(:,1)+wrtemp,0.)

! Explicit estimate of fluxes
call getqsat(qsat,dqdt,dtsurf,atm_ps)
ri=min(grav*(atm_zmin**2/atm_zmins)*(1.-dtsurf*srcp/atm_temp)/vmagn**2,rimax)
where (ri>=0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
  fq=fh
elsewhere        ! ri is -ve
  dumazmin=max(atm_zmin,zoseaice+0.2)
  root=sqrt(-ri*dumazmin/zoseaice)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  dumazmins=max(atm_zmins,zoseaice+0.2)
  root=sqrt(-ri*dumazmins/zoseaice)
  den=1.+chs*2.*bprm*aft*factchseaice*root
  fh=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*afq*factchseaice*root
  fq=1.-2.*bprm*ri/den
end where
! egice is for evaporating (lv).  Melting is included with lf.
dgice%wetfrac=max(1.+.008*min(dtsurf-273.16,0.),0.)
dgice%cd =af*fm*vmagn
dgice%cdh=aft*fh*vmagn
dgice%cdq=afq*fq*vmagn
dgice%umod=vmagn
dgice%fg=rho*dgice%cdh*cpair*(dtsurf-atm_temp/srcp)
dgice%fg=min(max(dgice%fg,-3000.),3000.)
dgice%eg=dgice%wetfrac*rho*dgice%cdq*lv*(qsat-atm_qg)
dgice%eg=min(dgice%eg,d_ndic*qice/(ls*dt))
dgice%eg=min(max(dgice%eg,-3000.),3000.)

du=fluxwgt*water%u(:,1)+(1.-fluxwgt)*water%utop-ice%u
dv=fluxwgt*water%v(:,1)+(1.-fluxwgt)*water%vtop-ice%v
icemagn=sqrt(max(du**2+dv**2,1.E-8))
dgice%cd_bot = 0.00536*icemagn

! MJT notes - use dtsurf for outgoing longwave for consistency with radiation code
! energy conservation is violated if initial conditions are poor
d_ftop=-dgice%fg-dgice%eg*lf/lv+atm_rg-emisice*sbconst*dtsurf**4+atm_sg*(1.-alb)*(1.-eye) ! first guess
d_ftop=d_ftop+lf*atm_rnd ! converting any rain to snowfall over ice
bot=rho*(dgice%cdh*cpair+dgice%cdq*dgice%wetfrac*dqdt*ls)

! iterative method to estimate ice velocity after stress from wind and currents are applied
newiu=ice%u
newiv=ice%v
do itr=1,10 ! maximum number of iterations
  uu=atm_u-newiu
  vv=atm_v-newiv
  du=fluxwgt*water%u(:,1)+(1.-fluxwgt)*water%utop-newiu
  dv=fluxwgt*water%v(:,1)+(1.-fluxwgt)*water%vtop-newiv
  vmagn=sqrt(max(uu**2+vv**2,1.E-4))
  icemagn=sqrt(max(du**2+dv**2,1.E-8))
  dgice%cd = af*fm*vmagn
  dgice%cd_bot = 0.00536*icemagn
  ! MJT notes - use wrtrho reference density for Boussinesq fluid approximation
  g=ice%u-newiu+dt*(rho*dgice%cd*uu+wrtrho*dgice%cd_bot*du)/dgice%imass
  h=ice%v-newiv+dt*(rho*dgice%cd*vv+wrtrho*dgice%cd_bot*dv)/dgice%imass
  dgu=-1.-dt*(rho*dgice%cd*(1.+(uu/vmagn)**2)+wrtrho*dgice%cd_bot*(1.+(du/icemagn)**2))/dgice%imass
  dhu=-dt*(rho*dgice%cd*uu*vv/(vmagn**2)+wrtrho*dgice%cd_bot*du*dv/(icemagn**2))/dgice%imass
  dgv=dhu
  dhv=-1.-dt*(rho*dgice%cd*(1.+(vv/vmagn)**2)+wrtrho*dgice%cd_bot*(1.+(dv/icemagn)**2))/dgice%imass
  ! Min det is around 1.
  det=dgu*dhv-dgv*dhu
  newiu=newiu-0.9*( g*dhv-h*dgv)/det ! 0.9 is to help solution converge by underestimating the change each iteration.
  newiv=newiv-0.9*(-g*dhu+h*dgu)/det
end do

! momentum transfer
uu=atm_u-newiu
vv=atm_v-newiv
du=fluxwgt*water%u(:,1)+(1.-fluxwgt)*water%utop-newiu
dv=fluxwgt*water%v(:,1)+(1.-fluxwgt)*water%vtop-newiv
vmagn=sqrt(max(uu**2+vv**2,1.E-4))
icemagn=sqrt(max(du**2+dv**2,1.E-4))
dgice%cd = af*fm*vmagn
dgice%cd_bot = 0.00536*icemagn
if ( mlo_uvcoupl==0 ) then
  dgice%cd_bot = 0.  
end if    
dgice%tauxica=rho*dgice%cd*uu
dgice%tauyica=rho*dgice%cd*vv
! MJT notes - use wrtrho reference density for Boussinesq fluid approximation
dgice%tauxicw=-wrtrho*dgice%cd_bot*du
dgice%tauyicw=-wrtrho*dgice%cd_bot*dv
ustar=sqrt(sqrt(max(dgice%tauxicw**2+dgice%tauyicw**2,0.))/wrtrho)
ustar=max(ustar,5.E-4)

! bottom flux
d_fb=cp0*wrtrho*0.006*ustar*(d_tb-d_timelt)
d_fb=min(max(d_fb,-1000.),1000.)

wavail=max(depth%depth_hl(:,wlev+1)+water%eta-minwater,0.)
wavail=min( wavail, max(delwater+water%eta,0.) )
f0=wavail*wrtrho/rhoic*qice/dt
d_fb=max(d_fb,-f0)

! Re-calculate fluxes to prevent overshoot (predictor-corrector)
! MJT notes - use dtsurf for outgoing longwave for consistency with radiation code
tnew=min(dtsurf+d_ftop/(gamm/dt+bot),273.2)
tnew=0.5*(tnew+dtsurf)
call getqsat(qsatnew,dqdt,tnew,atm_ps)
dgice%fg=rho*dgice%cdh*cpair*(tnew-atm_temp/srcp)
dgice%fg=min(max(dgice%fg,-3000.),3000.)
dgice%eg=dgice%wetfrac*rho*dgice%cdq*lv*(qsatnew-atm_qg)
dgice%eg=min(dgice%eg,d_ndic*qice/(ls*dt))
dgice%eg=min(max(dgice%eg,-3000.),3000.)
d_ftop=-dgice%fg-dgice%eg*lf/lv+atm_rg-emisice*sbconst*dtsurf**4+atm_sg*(1.-alb)*(1.-eye)
! Add flux of heat due to converting any rain to snowfall over ice
d_ftop=d_ftop+lf*atm_rnd ! rain (mm/sec) to W/m**2

! update water boundary conditions
! MJT notes - use wrtrho reference density for Boussinesq fluid approximation
dgwater%wu0 = dgwater%wu0 - ice%fracice*dgice%tauxicw/wrtrho
dgwater%wv0 = dgwater%wv0 - ice%fracice*dgice%tauyicw/wrtrho
dgwater%wt0_fb = ice%fracice*d_fb/(wrtrho*cp0)
dgwater%wt0 = dgwater%wt0 + dgwater%wt0_fb

! update snow depth
d_ndsn = d_ndsn + dt*(atm_rnd+atm_snd)/rhosn

return
end subroutine iceflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine screen diagnostics

subroutine scrncalc(atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins,diag, &
                    dgice,dgscrn,dgwater,ice,water)

implicit none

integer, intent(in) :: diag
integer iqw
real, dimension(imax), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_ps,atm_zmin,atm_zmins
type(dgicedata), intent(in) :: dgice
type(dgscrndata), intent(inout) :: dgscrn
type(dgwaterdata), intent(in) :: dgwater
type(icedata), intent(in) :: ice
type(waterdata), intent(in) :: water
real, dimension(imax) :: tscrn,qgscrn,uscrn,u10,dumtemp
real, dimension(imax) :: smixr,qsat,dqdt,atu_v,atv_v
real, dimension(imax) :: zohseaice,zoqseaice,zoxseaice
real dmag, atu, atv

if (diag>=1) write(6,*) "Calculate 2m diagnostics"

! water
dumtemp=max(water%temp(:,1)+wrtemp,0.)
call getqsat(qsat,dqdt,dumtemp,atm_ps)
if (zomode==0) then
  smixr=0.98*qsat
else
  smixr=qsat
end if
atu_v=atm_u-water%u(:,1)
atv_v=atm_v-water%v(:,1)
call scrntile(dgscrn%temp,dgscrn%qg,uscrn,u10,dgwater%zo,dgwater%zoh,dgwater%zoq,dumtemp, &
    smixr,atu_v,atv_v,atm_temp,atm_qg,atm_zmin,atm_zmins,diag)
do iqw = 1,imax
  dmag=sqrt(max(atu_v(iqw)**2+atv_v(iqw)**2,1.E-4))
  atu=(atm_u(iqw)-water%u(iqw,1))*uscrn(iqw)/dmag+water%u(iqw,1)
  atv=(atm_v(iqw)-water%v(iqw,1))*uscrn(iqw)/dmag+water%v(iqw,1)
  dgscrn%u2(iqw)=sqrt(atu**2+atv**2)
  atu=(atm_u(iqw)-water%u(iqw,1))*u10(iqw)/dmag+water%u(iqw,1)
  atv=(atm_v(iqw)-water%v(iqw,1))*u10(iqw)/dmag+water%v(iqw,1)
  dgscrn%u10(iqw)=sqrt(atu**2+atv**2)
end do  

! ice
call getqsat(qsat,dqdt,ice%tsurf,atm_ps)
smixr=dgice%wetfrac*qsat+(1.-dgice%wetfrac)*min(qsat,atm_qg)
atu_v=atm_u-ice%u
atv_v=atm_v-ice%v
zoxseaice(:)=zoseaice
zohseaice(:)=zoseaice/(factchseaice*factchseaice)
zoqseaice(:)=zohseaice(:)
call scrntile(tscrn,qgscrn,uscrn,u10,zoxseaice,zohseaice,zoqseaice,ice%tsurf, &
    smixr,atu_v,atv_v,atm_temp,atm_qg,atm_zmin,atm_zmins,diag)
do iqw = 1,imax
  dgscrn%temp(iqw)=(1.-ice%fracice(iqw))*dgscrn%temp(iqw)+ice%fracice(iqw)*tscrn(iqw)
  dgscrn%qg(iqw)=(1.-ice%fracice(iqw))*dgscrn%qg(iqw)+ice%fracice(iqw)*qgscrn(iqw)
  dmag=sqrt(max(atu_v(iqw)**2+atv_v(iqw)**2,1.E-4))
  atu=(atm_u(iqw)-ice%u(iqw))*uscrn(iqw)/dmag+ice%u(iqw)
  atv=(atm_v(iqw)-ice%v(iqw))*uscrn(iqw)/dmag+ice%v(iqw)
  dgscrn%u2=(1.-ice%fracice(iqw))*dgscrn%u2(iqw)+ice%fracice(iqw)*sqrt(atu**2+atv**2)
  atu=(atm_u(iqw)-ice%u(iqw))*u10(iqw)/dmag+ice%u(iqw)
  atv=(atm_v(iqw)-ice%v(iqw))*u10(iqw)/dmag+ice%v(iqw)
  dgscrn%u10=(1.-ice%fracice(iqw))*dgscrn%u10(iqw)+ice%fracice(iqw)*sqrt(atu**2+atv**2)
end do

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! screen diagnostic for individual tile

pure subroutine scrntile(tscrn,qgscrn,uscrn,u10,zo,zoh,zoq,stemp,smixr,atm_u,atm_v, &
                         atm_temp,atm_qg,atm_zmin,atm_zmins,diag)

implicit none

integer, intent(in) :: diag
integer ic
real, dimension(imax), intent(in) :: atm_u,atm_v,atm_temp,atm_qg,atm_zmin,atm_zmins
real, dimension(imax), intent(in) :: zo,zoh,zoq,stemp,smixr
real, dimension(imax), intent(out) :: tscrn,qgscrn,uscrn,u10
real, dimension(imax) :: umag,sig
real, dimension(imax) :: lzom,lzoh,lzoq,thetav,sthetav
real, dimension(imax) :: thetavstar,z_on_l,zs_on_l,z0_on_l,z0s_on_l,zt_on_l,zq_on_l
real, dimension(imax) :: pm0,ph0,pq0,pm1,ph1,integralm,integralh,integralq
real, dimension(imax) :: ustar,qstar,z10_on_l
real, dimension(imax) :: neutrals,neutral,neutral10,pm10
real, dimension(imax) :: integralm10,tstar,scrp,dumzmin,dumzmins
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
dumzmin=max(atm_zmin,zo+0.2)
dumzmins=max(atm_zmins,zo+0.2)
lzom=log(dumzmin/zo)
lzoh=log(dumzmins/zoh)
lzoq=log(dumzmins/zoq)

! Dyer and Hicks approach 
thetavstar=vkar*(thetav-sthetav)/lzoh
ustar=vkar*umag/lzom
do ic=1,nc
  z_on_l  = dumzmin*vkar*grav*thetavstar/(thetav*ustar**2)
  z_on_l  = min(z_on_l,10.)
  zs_on_l = z_on_l*dumzmins/dumzmin  
  z0_on_l = z_on_l*zo/dumzmin
  zt_on_l = z_on_l*zoh/dumzmin
  zq_on_l = z_on_l*zoq/dumzmin
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
z0s_on_l  = z0*zs_on_l/dumzmins
z0_on_l   = z0*z_on_l/dumzmin
z10_on_l  = z10*z_on_l/dumzmin
z0s_on_l  = min(z0s_on_l,10.)
z0_on_l   = min(z0_on_l,10.)
z10_on_l  = min(z10_on_l,10.)
neutrals  = log(dumzmins/z0)
neutral   = log(dumzmin/z0)
neutral10 = log(dumzmin/z10)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update exponential weighted moving average

subroutine mlo_ema_uvw(dt,water,depth)

implicit none

integer ii, iqw
real, intent(in) :: dt
real alpha, nstep
real, dimension(imax,wlev) :: u_hl, v_hl
real, dimension(imax,wlev) :: fdepth_hl
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth

nstep = max( mlo_timeave_length/dt, 1. ) ! this is a real value
alpha = 2./(nstep + 1.) 

do ii = 1,wlev
  do iqw = 1,imax
    water%u_ema(iqw,ii) = alpha*water%u(iqw,ii) + (1.-alpha)*water%u_ema(iqw,ii)
    water%v_ema(iqw,ii) = alpha*water%v(iqw,ii) + (1.-alpha)*water%v_ema(iqw,ii)
    water%w_ema(iqw,ii) = alpha*water%w(iqw,ii) + (1.-alpha)*water%w_ema(iqw,ii)
  end do
end do
! vertical contribution to shear
fdepth_hl(:,2:wlev) = (depth%depth_hl(:,2:wlev)-depth%depth(:,1:wlev-1)) &
  /max(depth%depth(:,2:wlev)-depth%depth(:,1:wlev-1),1.e-8)
call interpolate_hl(water%u_ema,fdepth_hl,u_hl)    
call interpolate_hl(water%v_ema,fdepth_hl,v_hl)
do ii = 2,wlev-1
  do iqw = 1,imax  
    if ( depth%dz(iqw,ii)>1.e-4 ) then  
      water%dudz(iqw,ii) = (u_hl(iqw,ii+1)-u_hl(iqw,ii))/depth%dz(iqw,ii)
      water%dvdz(iqw,ii) = (v_hl(iqw,ii+1)-v_hl(iqw,ii))/depth%dz(iqw,ii)
    else
      water%dudz(iqw,ii) = 0.
      water%dvdz(iqw,ii) = 0.
    end if  
  end do    
end do    
  
return
end subroutine mlo_ema_uvw

subroutine mlo_ema_ts(dt,water,depth)

implicit none

integer ii, iqw
real, intent(in) :: dt
real alpha, nstep
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth

nstep = max( mlo_timeave_length/dt, 1. ) ! this is a real value
alpha = 2./(nstep + 1.) 

do ii = 1,wlev
  do iqw = 1,imax
    water%temp_ema(iqw,ii) = alpha*water%temp(iqw,ii) + (1.-alpha)*water%temp_ema(iqw,ii)
    water%sal_ema(iqw,ii) = alpha*water%sal(iqw,ii) + (1.-alpha)*water%sal_ema(iqw,ii)
  end do
end do
  
return
end subroutine mlo_ema_ts

subroutine mlo_ema_reset(dt,water,depth)

implicit none

integer ii, iqw
real, intent(in) :: dt
real alpha, nstep
type(waterdata), intent(inout) :: water
type(depthdata), intent(in) :: depth

nstep = max( mlo_timeave_length/dt, 1. ) ! this is a real value
alpha = 2./(nstep + 1.) 

do ii = 1,wlev
  do iqw = 1,imax
    water%u_ema(iqw,ii) = water%u(iqw,ii)
    water%v_ema(iqw,ii) = water%v(iqw,ii)  
    water%w_ema(iqw,ii) = 0.
    water%temp_ema(iqw,ii) = water%temp(iqw,ii)
    water%sal_ema(iqw,ii) = water%sal(iqw,ii)
  end do
end do
  
return
end subroutine mlo_ema_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal error checking routines

subroutine mlocheck(message,water_temp,water_sal,water_u,water_v,ice_tsurf,ice_temp, &
                    ice_u,ice_v,ice_thick)

implicit none

character(len=*), intent(in) :: message
real, dimension(:,:), intent(in), optional :: water_temp, water_sal, water_u, water_v
real, dimension(:), intent(in), optional :: ice_tsurf
real, dimension(:,:), intent(in), optional :: ice_temp
real, dimension(:), intent(in), optional :: ice_u, ice_v
real, dimension(:), intent(in), optional :: ice_thick

if ( present( water_temp ) ) then
  if ( any( water_temp+wrtemp<100. .or. water_temp+wrtemp>400. ) ) then
    write(6,*) "ERROR: water_temp is out-of-range in ",trim(message)
    write(6,*) "minval,maxval ",minval(water_temp+wrtemp),maxval(water_temp+wrtemp)
    write(6,*) "minloc,maxloc ",minloc(water_temp+wrtemp),maxloc(water_temp+wrtemp)
    stop -1
  end if
end if 

if ( present( water_sal ) ) then
  if ( any( water_sal<-1.e-4 .or. water_sal>120. ) ) then
    write(6,*) "ERROR: water_sal is out-of-range in ",trim(message)
    write(6,*) "minval,maxval ",minval(water_sal),maxval(water_sal)
    write(6,*) "minloc,maxloc ",minloc(water_sal),maxloc(water_sal)
    stop -1
  end if
end if

if ( present( water_u ) .and. present( water_v ) ) then
  if ( any(abs(water_u)>40.) .or. any(abs(water_v)>40.) ) then
    write(6,*) "ERROR: current out-of-range in ",trim(message)
    write(6,*) "u ",minval(water_u),maxval(water_u)
    write(6,*) "  ",minloc(water_u),maxloc(water_u)
    write(6,*) "v ",minval(water_v),maxval(water_v)
    write(6,*) "  ",minloc(water_v),maxloc(water_v)
    stop -1
  end if
end if  

if ( present( ice_tsurf ) ) then
  if ( any( ice_tsurf<100. .or. ice_tsurf>400. ) ) then
    write(6,*) "ERROR: ice_tsurf is out-of-range in ",trim(message)
    write(6,*) "minval,maxval ",minval(ice_tsurf),maxval(ice_tsurf)
    write(6,*) "minloc,maxloc ",minloc(ice_tsurf),maxloc(ice_tsurf)
    stop -1
  end if
end if  

if ( present( ice_temp ) ) then
  if ( any( ice_temp<100. .or. ice_temp>400. ) ) then
    write(6,*) "ERROR: ice_temp is out-of-range in ",trim(message)
    write(6,*) "minval,maxval ",minval(ice_temp),maxval(ice_temp)
    write(6,*) "minloc,maxloc ",minloc(ice_temp),maxloc(ice_temp)
    stop -1
  end if
end if  

if ( present( ice_u ) .and. present( ice_v ) ) then
  if ( any(abs(ice_u)>40.) .or. any(abs(ice_v)>40.) ) then
    write(6,*) "ERROR: ice velocity is out-of-range in ",trim(message)
    write(6,*) "u ",minval(ice_u),maxval(ice_u)
    write(6,*) "v ",minval(ice_v),maxval(ice_v)
    stop -1
  end if
end if  

!if ( present( ice_thick ) ) then
!  if ( any(ice_thick>icemax) ) then
!    write(6,*) "ERROR: ice thickness is out-of-range in ",trim(message)
!    write(6,*) "maxval ",maxval(ice_thick)
!    write(6,*) "maxloc ",maxloc(ice_thick)
!    stop -1
!  end if
!end if

return
end subroutine mlocheck

end module mlo
