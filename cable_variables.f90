!!$ cable_variables.f90
!!$
!!$ This file declares all non-local variables for CABLE, 
!!$ CSIRO land surface model
!!$
!!$ Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
!!$ Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!!$
!!$ Fortran-95 coding by Harvey Davies, Gab Abramowitz and Martin Dix
!!$ bugs to gabsun@gmail.com.

!=========================================================================
MODULE define_dimensions
  INTEGER 	        :: mp      ! # land points
  INTEGER, PARAMETER	:: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER	:: nrb = 3 ! # radiation bands
  INTEGER, PARAMETER	:: ms = 6  ! # soil layers
  INTEGER, PARAMETER	:: ncp = 3 ! # vegetation carbon stores
  INTEGER, PARAMETER	:: ncs = 2 ! # soil carbon stores
  ! i_d is default kind for representing integer values.
  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)
END MODULE define_dimensions
!=========================================================================
MODULE define_types
  ! Contains all variables which are not subroutine-internal
  USE define_dimensions
  ! Energy and water balance variables:
  TYPE balances_type 
     REAL(r_1), DIMENSION(:), POINTER :: drybal ! energy balance for dry canopy
     REAL(r_1), DIMENSION(:), POINTER :: ebal   ! energy balance per time step (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER :: ebal_tot ! cumulative energy balance (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER :: evap_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: osnowd0  ! snow depth, first time step
     REAL(r_1), DIMENSION(:), POINTER :: precip_tot ! cumulative precipitation (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnoff_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbal   ! water balance per time step (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbal_tot ! cumulative water balance (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbtot0 ! total soil water (mm), first time step
     REAL(r_1), DIMENSION(:), POINTER :: wetbal ! energy balance for wet canopy
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type 
     REAL(r_1), DIMENSION(:), POINTER :: albsoil ! soil reflectance
     REAL(r_1), DIMENSION(:), POINTER :: bch  ! parameter b in Campbell equation
     REAL(r_1), DIMENSION(:), POINTER :: c3   ! c3 drainage coeff (fraction)
     REAL(r_1), DIMENSION(:), POINTER :: clay ! fraction of soil which is clay
     REAL(r_1), DIMENSION(:), POINTER :: cnsd ! thermal conductivity of dry soil [W/m/K]
     REAL(r_1), DIMENSION(:), POINTER :: css  ! soil specific heat capacity [kJ/kg/K]
     REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL(r_1), DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER(i_d), DIMENSION(:), POINTER :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
     INTEGER(i_d), DIMENSION(:), POINTER :: ibp2   ! parameter two in K vis suction (function of pbch)
     INTEGER(i_d), DIMENSION(:), POINTER :: isoilm ! integer soil type
     REAL(r_1), DIMENSION(:), POINTER :: rhosoil ! soil density [kg/m3]
     REAL(r_1), DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL(r_1), DIMENSION(:), POINTER :: sand  ! fraction of soil which is sand
     REAL(r_1), DIMENSION(:), POINTER :: sfc   ! vol H2O @ field capacity
     REAL(r_1), DIMENSION(:), POINTER :: silt  ! fraction of soil which is silt
     REAL(r_1), DIMENSION(:), POINTER :: ssat  ! vol H2O @ saturation
     REAL(r_1), DIMENSION(:), POINTER :: sucs  ! suction at saturation (m)
     REAL(r_1), DIMENSION(:), POINTER :: swilt ! vol H2O @ wilting
     REAL(r_1), DIMENSION(ms) :: zse   ! thickness of each soil layer (1=top) in m
     REAL(r_1), DIMENSION(ms+1) :: zshh ! distance between consecutive layer midpoints (m)
  END TYPE soil_parameter_type
  ! Soil and snow variables:
  TYPE soil_snow_type 
     REAL(r_1), DIMENSION(:,:), POINTER :: albsoilsn ! soil + snow reflectance
     REAL(r_1), DIMENSION(:), POINTER :: cls     ! factor for latent heat
     REAL(r_1), DIMENSION(:), POINTER :: dfn_dtg ! d(canopy%fns)/d(ssoil%tgg)
     REAL(r_1), DIMENSION(:), POINTER :: dfh_dtg ! d(canopy%fhs)/d(ssoil%tgg)
     REAL(r_1), DIMENSION(:), POINTER :: dfe_ddq ! d(canopy%fes)/d(dq)
     REAL(r_1), DIMENSION(:), POINTER :: ddq_dtg ! d(dq)/d(ssoil%tgg)
     REAL(r_2), DIMENSION(:,:), POINTER :: gammzz ! heat capacity for each soil layer
     INTEGER(i_d), DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow
     REAL(r_1), DIMENSION(:), POINTER :: osnowd  ! snow depth from previous time step
     REAL(r_1), DIMENSION(:), POINTER :: potev   ! potential evapotranspiration
     REAL(r_1), DIMENSION(:), POINTER :: runoff  ! total runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof1   ! surface runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof2   ! deep drainage (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rtsoil  ! turbulent resistance for soil
     REAL(r_1), DIMENSION(:,:), POINTER :: sdepth ! snow depth
     REAL(r_1), DIMENSION(:,:), POINTER  :: smass  ! snow mass
     REAL(r_1), DIMENSION(:), POINTER :: snage   ! snow age
     REAL(r_1), DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
     REAL(r_1), DIMENSION(:,:), POINTER  :: ssdn ! snow densities
     REAL(r_1), DIMENSION(:), POINTER :: ssdnn   ! average snow density
     REAL(r_1), DIMENSION(:,:), POINTER :: tgg   ! soil temperature in K
     REAL(r_1), DIMENSION(:,:), POINTER  :: tggsn ! snow temperature in K
     REAL(r_1), DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(:,:), POINTER :: wb    ! volumetric soil moisture (solid+liq)
     REAL(r_1), DIMENSION(:,:), POINTER :: wbfice !
     REAL(r_2), DIMENSION(:,:), POINTER :: wbice  ! soil ice
     REAL(r_2), DIMENSION(:,:), POINTER :: wblf !
     REAL(r_1), DIMENSION(:), POINTER :: wbtot   ! total soil water (mm)
  END TYPE soil_snow_type
  ! Vegetation parameters:
  TYPE veg_parameter_type
     INTEGER(i_d),DIMENSION(:), POINTER :: iveg ! vegetation type
     INTEGER(i_d),DIMENSION(:), POINTER :: meth ! method for calculation of canopy fluxes and temp.
     REAL(r_1), DIMENSION(:), POINTER :: vlai   ! leaf area index
     REAL(r_1), DIMENSION(:), POINTER :: vlaimax ! ???
     REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of resistances
     REAL(r_1), DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
     REAL(r_1), DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
     REAL(r_1), DIMENSION(:), POINTER :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
     REAL(r_1), DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: vbeta  ! 
     REAL(r_1), DIMENSION(:), POINTER :: hc	! roughness height of canopy (veg - snow)
     REAL(r_1), DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless)
     REAL(r_1), DIMENSION(:), POINTER :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
     REAL(r_1), DIMENSION(:), POINTER :: dleaf  ! chararacteristc legnth of leaf (m)
     REAL(r_1), DIMENSION(:), POINTER :: rp20   ! plant respiration coefficient at 20 C
     REAL(r_1), DIMENSION(:), POINTER :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
  END TYPE veg_parameter_type
  ! Canopy/vegetation variables:
  TYPE canopy_type
     REAL(r_1), DIMENSION(:), POINTER :: cansto ! canopy water storage (mm)
     REAL(r_1), DIMENSION(:), POINTER :: delwc ! change in canopy water store (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: dewmm ! dewfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: fe    ! total latent heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fh    ! total sensible heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fpn   ! plant photosynthesis (g C s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frp   ! plant respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frpw  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: frpr  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: frs   ! soil respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: fnee  ! net carbon flux (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frday ! daytime leaf resp
     REAL(r_1), DIMENSION(:), POINTER :: fnv   ! net rad. avail. to canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fev   ! latent hf from canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fevc  ! dry canopy transpiration (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fevw  ! lat heat fl wet canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fhv   ! sens heatfl from canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fhvw  ! sens heatfl from wet canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fns   ! net rad avail to soil (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fes   ! latent heatfl from soil (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fhs   ! sensible heat flux from soil
     REAL(r_1), DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
     REAL(r_1), DIMENSION(:), POINTER :: ga    ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: ghflux  ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: sghflux ! ground heat flux (W/m2) ???
     REAL(r_2), DIMENSION(:), POINTER :: dgdtg ! derivative of gflux wrt soil temp
     REAL(r_1), DIMENSION(:), POINTER :: through ! canopy throughfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: precis! throughfall to soil, after snow (mm)
     REAL(r_1), DIMENSION(:), POINTER :: rnet  ! net radiation absorbed by surface (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: spill ! can.storage excess after dewfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: wcint ! canopy rainfall interception (mm)
     REAL(r_1), DIMENSION(:), POINTER :: us    ! friction velocity
     REAL(r_1), DIMENSION(:), POINTER :: tscrn ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: qscrn ! specific humudity at screen height (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: uscrn ! wind speed at screen height (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: cduv  ! drag coefficient for momentum
  END TYPE canopy_type
  ! Radiation variables:
  TYPE radiation_type
     REAL(r_1), DIMENSION(:,:), POINTER :: albedo ! canopy+soil albedo
     REAL(r_1), DIMENSION(:), POINTER     :: extkb  ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER     :: extkd2 ! diffuse 2D radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
     REAL(r_1), DIMENSION(:), POINTER ::  extkn ! extinction coef for vertical nitrogen profile in canopy(-)
     REAL(r_1), DIMENSION(:), POINTER :: flws   ! soil long-wave radiation
     REAL(r_1), DIMENSION(:,:), POINTER  :: fvlai  ! leaf area index of big leaf
     REAL(r_1), DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
     REAL(r_1), DIMENSION(:), POINTER :: latitude  ! latitude
     REAL(r_1), DIMENSION(:), POINTER :: lwabv ! long wave absorbed by vegetation
     REAL(r_1), DIMENSION(:,:,:), POINTER :: qcan ! absorbed radiation for canopy (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER     :: qssabs ! absorbed short-wave radiation for soil
     REAL(r_1), DIMENSION(:,:), POINTER ::  rhocdf ! canopy diffuse reflectance (-)
     REAL(r_1), DIMENSION(:,:), POINTER  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
     REAL(r_1), DIMENSION(:,:), POINTER  :: scalex ! scaling PARAMETER for big leaf
     REAL(r_1), DIMENSION(:), POINTER     :: transd ! fraction SW diffuse transmitted through canopy
     REAL(r_1), DIMENSION(:), POINTER ::  trad  !  radiative temperature (soil and veg)
  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL(r_1), DIMENSION(:), POINTER	:: coexp ! Extinction coefficient for wind profile in canopy
     REAL(r_1), DIMENSION(:), POINTER	:: disp  ! zero-plane displacement
     REAL(r_1), DIMENSION(:), POINTER	:: hruff ! canopy height above snow level
     REAL(r_1), DIMENSION(:), POINTER	:: rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usa ! resistance from disp to hruf
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usb ! resistance from hruf to zruffs (or zref if zref<zruffs)
     REAL(r_1), DIMENSION(:), POINTER	:: rt1 ! 1/aerodynamic conductance
     REAL(r_1), DIMENSION(:), POINTER	:: term2, term3, term5, term6 ! for aerodynamic resistance calc.
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL(r_1), DIMENSION(:), POINTER	:: usuh ! Friction velocity/windspeed at canopy height
     REAL(r_1), DIMENSION(:), POINTER	:: za   ! level of lowest atmospheric model layer
     REAL(r_1), DIMENSION(:), POINTER	:: z0m  ! roughness length
     REAL(r_1), DIMENSION(:), POINTER	:: zref ! Reference height for met forcing
     REAL(r_1), DIMENSION(:), POINTER	:: zruffs ! SCALAR Roughness sublayer depth (ground=origin)
     REAL(r_1), DIMENSION(:), POINTER	:: z0soilsn ! roughness length of bare soil surface
     REAL(r_1), DIMENSION(:), POINTER	:: z0soil ! roughness length of bare soil surface
  END TYPE roughness_type
  ! Air variables:
  TYPE air_type
     REAL(r_1), DIMENSION(:), POINTER :: rho  ! dry air density (kg m-3)
     REAL(r_1), DIMENSION(:), POINTER :: volm ! molar volume (m3 mol-1)
     REAL(r_1), DIMENSION(:), POINTER :: rlam ! latent heat for water (j/kg)
     REAL(r_1), DIMENSION(:), POINTER :: qsat ! saturation specific humidity
     REAL(r_1), DIMENSION(:), POINTER :: epsi ! d(qsat)/dT ((kg/kg)/K)
     REAL(r_1), DIMENSION(:), POINTER :: visc ! air kinematic viscosity (m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: psyc ! psychrometric constant
     REAL(r_1), DIMENSION(:), POINTER :: dsatdk ! d(es)/dT (mb/K)
     REAL(r_1), DIMENSION(:), POINTER :: cmolar ! conv. from m/s to mol/m2/s
  END TYPE air_type
  ! Meterological data:
  TYPE met_type
     REAL(r_1), DIMENSION(:), POINTER :: ca   ! CO2 concentration (mol/mol)
     INTEGER(i_d), DIMENSION(:), POINTER :: year ! local time year AD 
     INTEGER(i_d), DIMENSION(:), POINTER :: moy  ! local time month of year 
     REAL(r_1), DIMENSION(:), POINTER :: doy  ! local time day of year = days since 0 hr 1st Jan 
     REAL(r_1), DIMENSION(:), POINTER :: hod  ! local hour of day
     REAL(r_1), DIMENSION(:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: precips ! solid preipitation only (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: tc	 ! surface air temperature (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tk	 ! surface air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvair   ! within canopy air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvrad   ! radiative vegetation temperature (K)
     REAL(r_1), DIMENSION(:), POINTER :: pmb     ! surface air pressure (mbar)
     REAL(r_1), DIMENSION(:), POINTER :: ua	 ! surface wind speed (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: qv	 ! surface specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: qvair   ! within canopy specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: da	 ! water vap pressure deficit at ref height (Pa)
     REAL(r_1), DIMENSION(:), POINTER :: dva	 ! in canopy water vap pressure deficit (Pa)
     REAL(r_1), DIMENSION(:), POINTER :: coszen  ! cos(zenith angle of sun)
  END TYPE met_type
  ! Cumulative flux variables:
  TYPE sum_flux_type
     REAL(r_1), DIMENSION(:), POINTER	:: sumpn  ! sum of canopy photosynthesis (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: sumrp  ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: sumrpw ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: sumrpr ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: sumrs  ! sum of soil respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: sumrd  ! sum of daytime respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER	:: dsumpn ! daily sumpn
     REAL(r_1), DIMENSION(:), POINTER	:: dsumrp ! daily sumrp
     REAL(r_1), DIMENSION(:), POINTER	:: dsumrs ! daily sumrs
     REAL(r_1), DIMENSION(:), POINTER	:: dsumrd ! daily sumrd
     REAL(r_1), DIMENSION(:), POINTER	:: sumxrp ! sum plant resp. modifier
     REAL(r_1), DIMENSION(:), POINTER	:: sumxrs ! sum soil resp. modifier
  END TYPE sum_flux_type
  TYPE bgc_pool_type
     REAL(r_1), DIMENSION(:,:), POINTER :: cplant ! plant carbon (g C/m2))
     REAL(r_1), DIMENSION(:,:), POINTER :: csoil  ! soil carbon (g C/m2)
     REAL(r_1), DIMENSION(ncp)	:: ratecp ! plant carbon rate constant (1/year)
     REAL(r_1), DIMENSION(ncs)	:: ratecs ! soil carbon rate constant (1/year)
  END TYPE bgc_pool_type

  ! Functions for allocating these types
  ! All overloaded so code only needs to call alloc_cbm_var
  ! Alloc routines could all initialise to NaN or zero for debugging?
  ! Don't need the mp argument here as it's a module variable.
  PUBLIC :: alloc_cbm_var
  PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type
  INTERFACE alloc_cbm_var
     MODULE PROCEDURE alloc_balances_type,            &
          alloc_soil_parameter_type,      &
          alloc_soil_snow_type,           &
          alloc_veg_parameter_type,       &
          alloc_canopy_type,              &
          alloc_radiation_type,           &
          alloc_roughness_type,           &
          alloc_air_type,                 &
          alloc_met_type,                 &
          alloc_sum_flux_type,            &
          alloc_bgc_pool_type
  END INTERFACE
  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type,            &
          dealloc_soil_parameter_type,      &
          dealloc_soil_snow_type,           &
          dealloc_veg_parameter_type,       &
          dealloc_canopy_type,              &
          dealloc_radiation_type,           &
          dealloc_roughness_type,           &
          dealloc_air_type,                 &
          dealloc_met_type,                 &
          dealloc_sum_flux_type,            &
          dealloc_bgc_pool_type
  END INTERFACE
CONTAINS
  
  SUBROUTINE alloc_balances_type(var, mp)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % drybal(mp) )
    ALLOCATE ( var % ebal(mp) )
    ALLOCATE ( var % ebal_tot(mp) )
    ALLOCATE ( var % evap_tot(mp) )
    ALLOCATE ( var % osnowd0(mp) )
    ALLOCATE ( var % precip_tot(mp) )
    ALLOCATE ( var % rnoff_tot(mp) )
    ALLOCATE ( var % wbal(mp) )
    ALLOCATE ( var % wbal_tot(mp) )
    ALLOCATE ( var % wbtot0(mp) )
    ALLOCATE ( var % wetbal(mp) )
  END SUBROUTINE alloc_balances_type

  SUBROUTINE alloc_soil_parameter_type(var, mp)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % albsoil(mp) )
    ALLOCATE ( var % bch(mp) )
    ALLOCATE ( var % c3(mp) )
    ALLOCATE ( var % clay(mp) )
    ALLOCATE ( var % cnsd(mp) )
    ALLOCATE ( var % css(mp) )
    ALLOCATE ( var % froot(mp,ms) )
    ALLOCATE ( var % hsbh(mp) )
    ALLOCATE ( var % hyds(mp) )
    ALLOCATE ( var % i2bp3(mp) )
    ALLOCATE ( var % ibp2(mp) )
    ALLOCATE ( var % isoilm(mp) )
    ALLOCATE ( var % rhosoil(mp) )
    ALLOCATE ( var % rs20(mp) )
    ALLOCATE ( var % sand(mp) )
    ALLOCATE ( var % sfc(mp) )
    ALLOCATE ( var % silt(mp) )
    ALLOCATE ( var % ssat(mp) )
    ALLOCATE ( var % sucs(mp) )
    ALLOCATE ( var % swilt(mp) )
  END SUBROUTINE alloc_soil_parameter_type
 
  SUBROUTINE alloc_soil_snow_type(var, mp)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % albsoilsn(mp,nrb) )
    ALLOCATE ( var % cls(mp) )
    ALLOCATE ( var % dfn_dtg(mp) )
    ALLOCATE ( var % dfh_dtg(mp) )
    ALLOCATE ( var % dfe_ddq(mp) )
    ALLOCATE ( var % ddq_dtg(mp) )
    ALLOCATE ( var % gammzz(mp,ms) )
    ALLOCATE ( var % isflag(mp) )
    ALLOCATE ( var % osnowd(mp) )
    ALLOCATE ( var % potev(mp) )
    ALLOCATE ( var % runoff(mp) )
    ALLOCATE ( var % rnof1(mp) )
    ALLOCATE ( var % rnof2(mp) )
    ALLOCATE ( var % rtsoil(mp) )
    ALLOCATE ( var % sdepth(mp,3) )
    ALLOCATE ( var % smass(mp,3) )
    ALLOCATE ( var % snage(mp) )
    ALLOCATE ( var % snowd(mp) )
    ALLOCATE ( var % ssdn(mp,3) )
    ALLOCATE ( var % ssdnn(mp) )
    ALLOCATE ( var % tgg(mp,ms) )
    ALLOCATE ( var % tggsn(mp,3) )
    ALLOCATE ( var % tss(mp) )
    ALLOCATE ( var % wb(mp,ms) )
    ALLOCATE ( var % wbfice(mp,ms) )
    ALLOCATE ( var % wbice(mp,ms) )
    ALLOCATE ( var % wblf(mp,ms) )
    ALLOCATE ( var % wbtot(mp) )
  END SUBROUTINE alloc_soil_snow_type
   
  SUBROUTINE alloc_veg_parameter_type(var, mp)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % iveg(mp) )
    ALLOCATE ( var % meth(mp) )
    ALLOCATE ( var % vlai(mp) )
    ALLOCATE ( var % vlaimax(mp) )
    ALLOCATE ( var % vlaiw(mp) )
    ALLOCATE ( var % fwet(mp) )
    ALLOCATE ( var % canst1(mp) )
    ALLOCATE ( var % ejmax(mp) )
    ALLOCATE ( var % frac4(mp) )
    ALLOCATE ( var % tminvj(mp) )
    ALLOCATE ( var % tmaxvj(mp) )
    ALLOCATE ( var % vbeta(mp) )
    ALLOCATE ( var % hc(mp) )
    ALLOCATE ( var % shelrb(mp) )
    ALLOCATE ( var % vcmax(mp) )
    ALLOCATE ( var % xfang(mp) )
    ALLOCATE ( var % dleaf(mp) )
    ALLOCATE ( var % rp20(mp) )
    ALLOCATE ( var % rpcoef(mp) )
  END SUBROUTINE alloc_veg_parameter_type
   
  SUBROUTINE alloc_canopy_type(var, mp)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % cansto(mp) )
    ALLOCATE ( var % delwc(mp) )
    ALLOCATE ( var % dewmm(mp) )
    ALLOCATE ( var % fe(mp) )
    ALLOCATE ( var % fh(mp) )
    ALLOCATE ( var % fpn(mp) )
    ALLOCATE ( var % frp(mp) )
    ALLOCATE ( var % frpw(mp) )
    ALLOCATE ( var % frpr(mp) )
    ALLOCATE ( var % frs(mp) )
    ALLOCATE ( var % fnee(mp) )
    ALLOCATE ( var % frday(mp) )
    ALLOCATE ( var % fnv(mp) )
    ALLOCATE ( var % fev(mp) )
    ALLOCATE ( var % fevc(mp) )
    ALLOCATE ( var % fevw(mp) )
    ALLOCATE ( var % fhv(mp) )
    ALLOCATE ( var % fhvw(mp) )
    ALLOCATE ( var % fns(mp) )
    ALLOCATE ( var % fes(mp) )
    ALLOCATE ( var % fhs(mp) )
    ALLOCATE ( var % tv(mp) )
    ALLOCATE ( var % ga(mp) )
    ALLOCATE ( var % ghflux(mp) )
    ALLOCATE ( var % sghflux(mp) )
    ALLOCATE ( var % dgdtg(mp) )
    ALLOCATE ( var % through(mp) )
    ALLOCATE ( var % precis(mp) )
    ALLOCATE ( var % rnet(mp) )
    ALLOCATE ( var % spill(mp) )
    ALLOCATE ( var % wcint(mp) )
    ALLOCATE ( var % us(mp) )
    ALLOCATE ( var % tscrn(mp) )
    ALLOCATE ( var % qscrn(mp) )
    ALLOCATE ( var % uscrn(mp) )
    ALLOCATE ( var % cduv(mp) )
  END SUBROUTINE alloc_canopy_type
   
  SUBROUTINE alloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % albedo(mp,nrb) )
    ALLOCATE ( var % extkb(mp) )
    ALLOCATE ( var % extkd2(mp) )
    ALLOCATE ( var % extkd(mp) )
    ALLOCATE ( var % extkn(mp) )
    ALLOCATE ( var % flws(mp) )
    ALLOCATE ( var % fvlai(mp,mf) )
    ALLOCATE ( var % gradis(mp,mf) )
    ALLOCATE ( var % latitude(mp) )
    ALLOCATE ( var % lwabv(mp) )
    ALLOCATE ( var % qcan(mp,mf,nrb) )
    ALLOCATE ( var % qssabs(mp) )
    ALLOCATE ( var % rhocdf(mp,nrb) )
    ALLOCATE ( var % rniso(mp,mf) )
    ALLOCATE ( var % scalex(mp,mf) )
    ALLOCATE ( var % transd(mp) )
    ALLOCATE ( var % trad(mp) )
  END SUBROUTINE alloc_radiation_type
   
  SUBROUTINE alloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % coexp(mp) )
    ALLOCATE ( var % disp(mp) )
    ALLOCATE ( var % hruff(mp) )
    ALLOCATE ( var % rt0us(mp) )
    ALLOCATE ( var % rt1usa(mp) )
    ALLOCATE ( var % rt1usb(mp) )
    ALLOCATE ( var % rt1(mp) )
    ALLOCATE ( var % term2(mp) )
    ALLOCATE ( var % term3(mp) )
    ALLOCATE ( var % term5(mp) )
    ALLOCATE ( var % term6(mp) )
    ALLOCATE ( var % usuh(mp) )
    ALLOCATE ( var % za(mp) )
    ALLOCATE ( var % z0m(mp) )
    ALLOCATE ( var % zref(mp) )
    ALLOCATE ( var % zruffs(mp) )
    ALLOCATE ( var % z0soilsn(mp) )
    ALLOCATE ( var % z0soil(mp) )
  END SUBROUTINE alloc_roughness_type
   
  SUBROUTINE alloc_air_type(var, mp)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % rho(mp) )
    ALLOCATE ( var % volm(mp) )
    ALLOCATE ( var % rlam(mp) )
    ALLOCATE ( var % qsat(mp) )
    ALLOCATE ( var % epsi(mp) )
    ALLOCATE ( var % visc(mp) )
    ALLOCATE ( var % psyc(mp) )
    ALLOCATE ( var % dsatdk(mp) )
    ALLOCATE ( var % cmolar(mp) )
  END SUBROUTINE alloc_air_type
   
  SUBROUTINE alloc_met_type(var, mp)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % ca(mp) )
    ALLOCATE ( var % year(mp) )
    ALLOCATE ( var % moy(mp) )
    ALLOCATE ( var % doy(mp) )
    ALLOCATE ( var % hod(mp) )
    ALLOCATE ( var % fsd(mp) )
    ALLOCATE ( var % fld(mp) )
    ALLOCATE ( var % precip(mp) )
    ALLOCATE ( var % precips(mp) )
    ALLOCATE ( var % tc(mp) )
    ALLOCATE ( var % tk(mp) )
    ALLOCATE ( var % tvair(mp) )
    ALLOCATE ( var % tvrad(mp) )
    ALLOCATE ( var % pmb(mp) )
    ALLOCATE ( var % ua(mp) )
    ALLOCATE ( var % qv(mp) )
    ALLOCATE ( var % qvair(mp) )
    ALLOCATE ( var % da(mp) )
    ALLOCATE ( var % dva(mp) )
    ALLOCATE ( var % coszen(mp) )
  END SUBROUTINE alloc_met_type
   
  SUBROUTINE alloc_sum_flux_type(var, mp)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % sumpn(mp) )
    ALLOCATE ( var % sumrp(mp) )
    ALLOCATE ( var % sumrpw(mp) )
    ALLOCATE ( var % sumrpr(mp) )
    ALLOCATE ( var % sumrs(mp) )
    ALLOCATE ( var % sumrd(mp) )
    ALLOCATE ( var % dsumpn(mp) )
    ALLOCATE ( var % dsumrp(mp) )
    ALLOCATE ( var % dsumrs(mp) )
    ALLOCATE ( var % dsumrd(mp) )
    ALLOCATE ( var % sumxrp(mp) )
    ALLOCATE ( var % sumxrs(mp) )
  END SUBROUTINE alloc_sum_flux_type

  SUBROUTINE alloc_bgc_pool_type(var, mp)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % cplant(mp,ncp) )
    ALLOCATE ( var % csoil(mp,ncs) )
  END SUBROUTINE alloc_bgc_pool_type

  ! Begin deallocation routines:
   SUBROUTINE dealloc_balances_type(var, mp)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % drybal )
    DEALLOCATE ( var % ebal )
    DEALLOCATE ( var % ebal_tot )
    DEALLOCATE ( var % evap_tot )
    DEALLOCATE ( var % osnowd0 )
    DEALLOCATE ( var % precip_tot )
    DEALLOCATE ( var % rnoff_tot )
    DEALLOCATE ( var % wbal )
    DEALLOCATE ( var % wbal_tot )
    DEALLOCATE ( var % wbtot0 )
    DEALLOCATE ( var % wetbal )
  END SUBROUTINE dealloc_balances_type

  SUBROUTINE dealloc_soil_parameter_type(var, mp)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % albsoil )
    DEALLOCATE ( var % bch )
    DEALLOCATE ( var % c3 )
    DEALLOCATE ( var % clay )
    DEALLOCATE ( var % cnsd )
    DEALLOCATE ( var % css )
    DEALLOCATE ( var % froot )
    DEALLOCATE ( var % hsbh )
    DEALLOCATE ( var % hyds )
    DEALLOCATE ( var % i2bp3 )
    DEALLOCATE ( var % ibp2 )
    DEALLOCATE ( var % isoilm )
    DEALLOCATE ( var % rhosoil )
    DEALLOCATE ( var % rs20 )
    DEALLOCATE ( var % sand )
    DEALLOCATE ( var % sfc )
    DEALLOCATE ( var % silt )
    DEALLOCATE ( var % ssat )
    DEALLOCATE ( var % sucs )
    DEALLOCATE ( var % swilt )
  END SUBROUTINE dealloc_soil_parameter_type
 
  SUBROUTINE dealloc_soil_snow_type(var, mp)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % albsoilsn )
    DEALLOCATE ( var % cls )
    DEALLOCATE ( var % dfn_dtg )
    DEALLOCATE ( var % dfh_dtg )
    DEALLOCATE ( var % dfe_ddq )
    DEALLOCATE ( var % ddq_dtg )
    DEALLOCATE ( var % gammzz )
    DEALLOCATE ( var % isflag )
    DEALLOCATE ( var % osnowd )
    DEALLOCATE ( var % potev )
    DEALLOCATE ( var % runoff )
    DEALLOCATE ( var % rnof1 )
    DEALLOCATE ( var % rnof2 )
    DEALLOCATE ( var % rtsoil )
    DEALLOCATE ( var % sdepth )
    DEALLOCATE ( var % smass )
    DEALLOCATE ( var % snage )
    DEALLOCATE ( var % snowd )
    DEALLOCATE ( var % ssdn )
    DEALLOCATE ( var % ssdnn )
    DEALLOCATE ( var % tgg )
    DEALLOCATE ( var % tggsn )
    DEALLOCATE ( var % tss )
    DEALLOCATE ( var % wb )
    DEALLOCATE ( var % wbfice )
    DEALLOCATE ( var % wbice )
    DEALLOCATE ( var % wblf )
    DEALLOCATE ( var % wbtot )
  END SUBROUTINE dealloc_soil_snow_type
   
  SUBROUTINE dealloc_veg_parameter_type(var, mp)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % iveg )
    DEALLOCATE ( var % meth )
    DEALLOCATE ( var % vlai )
    DEALLOCATE ( var % vlaimax )
    DEALLOCATE ( var % vlaiw )
    DEALLOCATE ( var % fwet )
    DEALLOCATE ( var % canst1 )
    DEALLOCATE ( var % ejmax )
    DEALLOCATE ( var % frac4 )
    DEALLOCATE ( var % tminvj )
    DEALLOCATE ( var % tmaxvj )
    DEALLOCATE ( var % vbeta )
    DEALLOCATE ( var % hc )
    DEALLOCATE ( var % shelrb )
    DEALLOCATE ( var % vcmax )
    DEALLOCATE ( var % xfang )
    DEALLOCATE ( var % dleaf )
    DEALLOCATE ( var % rp20 )
    DEALLOCATE ( var % rpcoef )
  END SUBROUTINE dealloc_veg_parameter_type
   
  SUBROUTINE dealloc_canopy_type(var, mp)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % cansto )
    DEALLOCATE ( var % delwc )
    DEALLOCATE ( var % dewmm )
    DEALLOCATE ( var % fe )
    DEALLOCATE ( var % fh )
    DEALLOCATE ( var % fpn )
    DEALLOCATE ( var % frp )
    DEALLOCATE ( var % frpw )
    DEALLOCATE ( var % frpr )
    DEALLOCATE ( var % frs )
    DEALLOCATE ( var % fnee )
    DEALLOCATE ( var % frday )
    DEALLOCATE ( var % fnv )
    DEALLOCATE ( var % fev )
    DEALLOCATE ( var % fevc )
    DEALLOCATE ( var % fevw )
    DEALLOCATE ( var % fhv )
    DEALLOCATE ( var % fhvw )
    DEALLOCATE ( var % fns )
    DEALLOCATE ( var % fes )
    DEALLOCATE ( var % fhs )
    DEALLOCATE ( var % tv )
    DEALLOCATE ( var % ga )
    DEALLOCATE ( var % ghflux )
    DEALLOCATE ( var % sghflux )
    DEALLOCATE ( var % dgdtg )
    DEALLOCATE ( var % through )
    DEALLOCATE ( var % precis )
    DEALLOCATE ( var % rnet )
    DEALLOCATE ( var % spill )
    DEALLOCATE ( var % wcint )
    DEALLOCATE ( var % us )
    DEALLOCATE ( var % tscrn )
    DEALLOCATE ( var % qscrn )
    DEALLOCATE ( var % uscrn )
    DEALLOCATE ( var % cduv )
  END SUBROUTINE dealloc_canopy_type
   
  SUBROUTINE dealloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % albedo )
    DEALLOCATE ( var % extkb )
    DEALLOCATE ( var % extkd2 )
    DEALLOCATE ( var % extkd )
    DEALLOCATE ( var % extkn )
    DEALLOCATE ( var % flws )
    DEALLOCATE ( var % fvlai )
    DEALLOCATE ( var % gradis )
    DEALLOCATE ( var % latitude )
    DEALLOCATE ( var % lwabv )
    DEALLOCATE ( var % qcan )
    DEALLOCATE ( var % qssabs )
    DEALLOCATE ( var % rhocdf )
    DEALLOCATE ( var % rniso )
    DEALLOCATE ( var % scalex )
    DEALLOCATE ( var % transd )
    DEALLOCATE ( var % trad )
  END SUBROUTINE dealloc_radiation_type
   
  SUBROUTINE dealloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % coexp )
    DEALLOCATE ( var % disp )
    DEALLOCATE ( var % hruff )
    DEALLOCATE ( var % rt0us )
    DEALLOCATE ( var % rt1usa )
    DEALLOCATE ( var % rt1usb )
    DEALLOCATE ( var % rt1 )
    DEALLOCATE ( var % term2 )
    DEALLOCATE ( var % term3 )
    DEALLOCATE ( var % term5 )
    DEALLOCATE ( var % term6 )
    DEALLOCATE ( var % usuh )
    DEALLOCATE ( var % za )
    DEALLOCATE ( var % z0m )
    DEALLOCATE ( var % zref )
    DEALLOCATE ( var % zruffs )
    DEALLOCATE ( var % z0soilsn )
    DEALLOCATE ( var % z0soil )
  END SUBROUTINE dealloc_roughness_type
   
  SUBROUTINE dealloc_air_type(var, mp)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % rho )
    DEALLOCATE ( var % volm )
    DEALLOCATE ( var % rlam )
    DEALLOCATE ( var % qsat )
    DEALLOCATE ( var % epsi )
    DEALLOCATE ( var % visc )
    DEALLOCATE ( var % psyc )
    DEALLOCATE ( var % dsatdk )
    DEALLOCATE ( var % cmolar )
  END SUBROUTINE dealloc_air_type
   
  SUBROUTINE dealloc_met_type(var, mp)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % ca )
    DEALLOCATE ( var % year )
    DEALLOCATE ( var % moy )
    DEALLOCATE ( var % doy )
    DEALLOCATE ( var % hod )
    DEALLOCATE ( var % fsd )
    DEALLOCATE ( var % fld )
    DEALLOCATE ( var % precip )
    DEALLOCATE ( var % precips )
    DEALLOCATE ( var % tc )
    DEALLOCATE ( var % tk )
    DEALLOCATE ( var % tvair )
    DEALLOCATE ( var % tvrad )
    DEALLOCATE ( var % pmb )
    DEALLOCATE ( var % ua )
    DEALLOCATE ( var % qv )
    DEALLOCATE ( var % qvair )
    DEALLOCATE ( var % da )
    DEALLOCATE ( var % dva )
    DEALLOCATE ( var % coszen )
  END SUBROUTINE dealloc_met_type
   
  SUBROUTINE dealloc_sum_flux_type(var, mp)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % sumpn )
    DEALLOCATE ( var % sumrp )
    DEALLOCATE ( var % sumrpw )
    DEALLOCATE ( var % sumrpr )
    DEALLOCATE ( var % sumrs )
    DEALLOCATE ( var % sumrd )
    DEALLOCATE ( var % dsumpn )
    DEALLOCATE ( var % dsumrp )
    DEALLOCATE ( var % dsumrs )
    DEALLOCATE ( var % dsumrd )
    DEALLOCATE ( var % sumxrp )
    DEALLOCATE ( var % sumxrs )
  END SUBROUTINE dealloc_sum_flux_type

  SUBROUTINE dealloc_bgc_pool_type(var, mp)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % cplant )
    DEALLOCATE ( var % csoil )
  END SUBROUTINE dealloc_bgc_pool_type

END MODULE define_types
!=========================================================================
MODULE math_constants
  USE define_dimensions, ONLY : i_d, r_1
  REAL(r_1), PARAMETER	:: pi = 3.1415927
  REAL(r_1), PARAMETER	:: pi180 = pi / 180.0 ! radians / degree
  REAL(r_1), PARAMETER	:: two_pi = 2.0 * pi
END MODULE math_constants
!=========================================================================
MODULE physical_constants
  USE define_dimensions, ONLY : i_d, r_1
  REAL(r_1), PARAMETER :: capp   = 1004.64  ! air spec. heat (J/kg/K)
  REAL(r_1), PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
  REAL(r_1), PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
  REAL(r_1), PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
  REAL(r_1), PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
  REAL(r_1), PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water	(kg/mol)
  REAL(r_1), PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
  REAL(r_1), PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
  ! Teten coefficients
  REAL(r_1), PARAMETER :: tetena=6.106  ! ??? refs?
  REAL(r_1), PARAMETER :: tetenb=17.27
  REAL(r_1), PARAMETER :: tetenc=237.3
  ! Aerodynamic parameters, diffusivities, water density:
  REAL(r_1), PARAMETER :: vonk   = 0.40 ! von Karman constant
  REAL(r_1), PARAMETER :: a33    = 1.25 ! inertial sublayer sw/us
  REAL(r_1), PARAMETER :: csw    = 0.50 ! canopy sw decay (Weil theory)
  REAL(r_1), PARAMETER :: ctl    = 0.40 ! Wagga wheat (RDD 1992, Challenges)
  REAL(r_1), PARAMETER :: apol   = 0.70 ! Polhausen coeff: single-sided plate
  REAL(r_1), PARAMETER :: prandt = 0.71 ! Prandtl number: visc/diffh
  REAL(r_1), PARAMETER :: schmid = 0.60 ! Schmidt number: visc/diffw
  REAL(r_1), PARAMETER :: diffwc = 1.60 ! diffw/diffc = H2O/CO2 diffusivity
  REAL(r_1), PARAMETER :: rhow   = 1000.0 ! liquid water density   [kg/m3]
  REAL(r_1), PARAMETER :: emleaf = 1.0  ! leaf emissivity
  REAL(r_1), PARAMETER :: emsoil = 1.0  ! soil emissivity
  REAL(r_1), PARAMETER :: cr = 0.3	! element drag coefficient
  REAL(r_1), PARAMETER :: cs = 0.003    ! substrate drag coefficient
  REAL(r_1), PARAMETER :: beta  = cr/cs ! ratio cr/cs
  REAL(r_1), PARAMETER :: ccd   = 15.0  ! constant in d/h equation
  REAL(r_1), PARAMETER :: ccw   = 2.0   ! ccw=(zw-d)/(h-d)
  REAL(r_1), PARAMETER :: usuhm = 0.3   ! (max of us/uh)
  ! Turbulence parameters:
  INTEGER(i_d), PARAMETER :: niter = 4  ! number of iterations for za/L
  REAL(r_1), PARAMETER :: zetmul = 0.4  ! if niter=2, final zeta=zetmul*zetar(2)
  REAL(r_1), PARAMETER :: zeta0  = 0.0  ! initial value of za/L
  REAL(r_1), PARAMETER :: zetneg = -10. ! negative limit on za/L when niter>=3
  REAL(r_1), PARAMETER :: zetpos = 0.5  ! positive limit on za/L when niter>=3
  REAL(r_1), PARAMETER :: zdlin  = 1.0  ! height frac of d below which TL linear
  REAL(r_1), PARAMETER :: umin   = 1.
END MODULE physical_constants
!=========================================================================
MODULE other_constants
  USE define_dimensions, ONLY : i_d, r_1, nrb
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: refl = (/ 0.1, 0.425, 0.05 /) ! leaf reflectance
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: taul = (/ 0.1, 0.425, 0.05/)  ! leaf transmittance
  INTEGER(i_d), PARAMETER :: istemp = 4 ! soil temp:	 1,2,3,4 = FR,kf,mrr,mrrkf
  INTEGER(i_d), PARAMETER :: ismois = 2 ! soil moist:  1,2,3	 = MP84,NP89,Richards
  INTEGER(i_d), PARAMETER :: isinf  = 2 ! soil infilt: 1,2	 = MP84, FC96
  INTEGER(i_d), PARAMETER :: isevap = 2 ! soil evap: 1,2,3 = alfa,beta,threshold
  INTEGER(i_d), PARAMETER :: itherm = 1 ! VW or KGK algorithm for hconds,rkapps
  INTEGER(i_d), PARAMETER :: irktem = 5 ! RK steps in soil temp schemes
  INTEGER(i_d), PARAMETER :: irkmoi = 5 ! RK steps in soil moisture schemes
  ! soil water parameters:
  REAL(r_1), PARAMETER :: etarct = 0.7  ! rel soil moisture for finding zst1,zst2
  REAL(r_1), PARAMETER :: dbde = 1.3333 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
  INTEGER(i_d), PARAMETER :: istsw = 1  !
  INTEGER(i_d), PARAMETER :: iresp=0	! unscaled (iresp=0) or scaled (iresp=1) respiration
END MODULE other_constants
!=========================================================================
MODULE photosynthetic_constants
  USE define_dimensions, ONLY : i_d, r_1
  INTEGER(i_d), PARAMETER :: maxiter=20 ! max # interations for leaf temperature
  REAL(r_1), PARAMETER :: a1c3 = 9.0
  REAL(r_1), PARAMETER :: a1c4 = 4.0
  REAL(r_1), PARAMETER :: alpha3 = 0.200
  REAL(r_1), PARAMETER :: alpha4  = 0.05
  REAL(r_1), PARAMETER :: cfrd3  = 0.015
  REAL(r_1), PARAMETER :: cfrd4  = 0.025
  REAL(r_1), PARAMETER :: conkc0 = 302.e-6	!mol mol^-1
  REAL(r_1), PARAMETER :: conko0 = 256.e-3	!mol mol^-1
  REAL(r_1), PARAMETER :: convx3 = 1.0E-2
  REAL(r_1), PARAMETER :: convx4 = 0.8
  REAL(r_1), PARAMETER :: d0c3 = 1500.0
  REAL(r_1), PARAMETER :: d0c4 = 1500.0
  REAL(r_1), PARAMETER :: ekc = 59430.0	!J mol^-1
  REAL(r_1), PARAMETER :: eko = 36000.0	!J mol^-1
  REAL(r_1), PARAMETER :: gam0 = 28.0E-6	!mol mol^-1 @ 20C = 36.9 @ 25C
  REAL(r_1), PARAMETER :: gam1 = 0.0509
  REAL(r_1), PARAMETER :: gam2 = 0.0010
  REAL(r_1), PARAMETER :: gsw03  = 0.01
  REAL(r_1), PARAMETER :: gsw04  = 0.04
  REAL(r_1), PARAMETER :: rgbwc  = 1.32
  REAL(r_1), PARAMETER :: rgswc  = 1.57
  REAL(r_1), PARAMETER :: tmaxj  = 45.0
  REAL(r_1), PARAMETER :: tmaxv  = 45.0
  REAL(r_1), PARAMETER :: tminj  = -5.0
  REAL(r_1), PARAMETER :: tminv  = -5.0
  REAL(r_1), PARAMETER :: toptj  = 20.0
  REAL(r_1), PARAMETER :: toptv  = 20.0
  REAL(r_1), PARAMETER :: trefk= 298.2	!reference temperature K
END MODULE photosynthetic_constants
