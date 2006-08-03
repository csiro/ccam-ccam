!!$ cable_variables.f90
!!$
!!$ This file declares all non-local variables for CABLE, 
!!$ CSIRO land surface model
!!$
!!$ Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
!!$ Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!!$
!!$ Fortran-95 coding by Harvey Davies and Gab Abramowitz,
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
     REAL(r_1), DIMENSION(:), POINTER :: cnsd ! specific heat capapcity of soil (unit??)
     REAL(r_1), DIMENSION(:), POINTER :: css  ! soil heat capacity [kJ/kg/C]
     REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL(r_1), DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER(i_d), DIMENSION(:), POINTER :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
     INTEGER(i_d), DIMENSION(:), POINTER :: ibp2   !parameter two in K vis suction (function of pbch)
     INTEGER(i_d), DIMENSION(:), POINTER :: isoilm ! soil type
     REAL(r_1), DIMENSION(:), POINTER :: rhosoil ! soil density [kg/m3]
     REAL(r_1), DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL(r_1), DIMENSION(:), POINTER :: sand  ! fraction of soil which is sand
     REAL(r_1), DIMENSION(:), POINTER :: sfc   ! vol H2O @ field capacity
     REAL(r_1), DIMENSION(:), POINTER :: silt  ! fraction of soil which is silt
     REAL(r_1), DIMENSION(:), POINTER :: ssat  ! vol H2O @ saturation
     REAL(r_1), DIMENSION(:), POINTER :: sucs  ! suction at saturation (m)
     REAL(r_1), DIMENSION(:), POINTER :: swilt ! vol H2O @ wilting
     REAL(r_1), DIMENSION(ms) :: zse   ! thickness of each soil layer (1=top) in m
     REAL(r_1), DIMENSION(ms+1) :: zshh ! depth from surface of mid-point of a layer (m)
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
     REAL(r_1), DIMENSION(:,:), POINTER :: sdepth ! 
     REAL(r_1), DIMENSION(:,:), POINTER  :: smass  ! snow mass
     REAL(r_1), DIMENSION(:), POINTER :: snage   ! snow age
     REAL(r_1), DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
     REAL(r_1), DIMENSION(:,:), POINTER  :: ssdn ! snow densities
     REAL(r_1), DIMENSION(:), POINTER :: ssdnn   ! average snow density
     REAL(r_1), DIMENSION(:,:), POINTER :: tgg  ! soil temperature in K
     REAL(r_1), DIMENSION(:,:), POINTER  :: tggsn  ! snow temperature in K
     REAL(r_1), DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(:,:), POINTER :: wb   ! volumetric soil moisture (solid+liq)
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
     REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of reistances
     REAL(r_1), DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
     REAL(r_1), DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
     REAL(r_1), DIMENSION(:), POINTER :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
     REAL(r_1), DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: vbeta  ! ???
     REAL(r_1), DIMENSION(:), POINTER :: hc	! roughness height of canopy (veg - snow)
     REAL(r_1), DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless)
     REAL(r_1), DIMENSION(:), POINTER :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
     REAL(r_1), DIMENSION(:), POINTER :: dleaf ! chararacteristc legnth of leaf (m)
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
     REAL(r_1), DIMENSION(:), POINTER :: fhs   !sensible heat fl from soil
     REAL(r_1), DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
     REAL(r_1), DIMENSION(:), POINTER :: ga    ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: ghflux  ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: sghflux ! ground heat flux (W/m2) ???
     REAL(r_2), DIMENSION(:), POINTER :: dgdtg ! derivative of gflux wrt soil temp
     REAL(r_1), DIMENSION(:), POINTER :: through ! canopy throughfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: precis! throughfall to soil, after snow (mm)
     REAL(r_1), DIMENSION(:), POINTER :: rnet  ! ???
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
     REAL(r_1), DIMENSION(:), POINTER     :: extkd2 ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
     REAL(r_1), DIMENSION(:), POINTER ::  extkn ! vertical leaf nitrogen profile (-)
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
     REAL(r_1), DIMENSION(:), POINTER ::  trad  ! ???
  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL(r_1), DIMENSION(:), POINTER	:: coexp
     REAL(r_1), DIMENSION(:), POINTER	:: disp     ! zero-plane displacement
     REAL(r_1), DIMENSION(:), POINTER	:: rt0us
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usa
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usb
     REAL(r_1), DIMENSION(:), POINTER	:: rt1 ! 1/aerodynamic conductance
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL(r_1), DIMENSION(:), POINTER	:: usuh 
     REAL(r_1), DIMENSION(:), POINTER	:: za
     REAL(r_1), DIMENSION(:), POINTER	:: z0m      ! roughness length
     REAL(r_1), DIMENSION(:), POINTER	:: zref
     REAL(r_1), DIMENSION(:), POINTER	:: zruffs
     REAL(r_1), DIMENSION(:), POINTER	:: z0soilsn
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
     REAL(r_1), DIMENSION(:), POINTER :: ca   ! CO2 concentration
     INTEGER(i_d), DIMENSION(:), POINTER :: year ! local time year AD 
     INTEGER(i_d), DIMENSION(:), POINTER :: moy  ! local time month of year 
     REAL(r_1), DIMENSION(:), POINTER :: doy  ! local time day of year = days since 0 hr 1st Jan 
     REAL(r_1), DIMENSION(:), POINTER :: hod  ! local hour of day
     REAL(r_1), DIMENSION(:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: precips ! solid only (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: tc	 ! surface air temperature (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tk	 ! surface air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvair   ! within canopy air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvrad   ! ???
     REAL(r_1), DIMENSION(:), POINTER :: pmb     ! surface air pressure (mbar)
     REAL(r_1), DIMENSION(:), POINTER :: ua	 ! surface wind speed (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: qv	 ! surface specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: qvair   ! within canopy specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: da	 ! surf water vap pressure deficit (Pa)
     REAL(r_1), DIMENSION(:), POINTER :: dva	 ! surf water vap pressure deficit (Pa) ???
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
  public :: alloc_cbm_var
  private :: alloc_bgc_pool_type
  interface alloc_cbm_var
    module procedure alloc_balances_type,            &
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
  end interface

CONTAINS

  subroutine alloc_balances_type(var, mp)
    type(balances_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % drybal(mp) )
    allocate ( var % ebal(mp) )
    allocate ( var % ebal_tot(mp) )
    allocate ( var % evap_tot(mp) )
    allocate ( var % osnowd0(mp) )
    allocate ( var % precip_tot(mp) )
    allocate ( var % rnoff_tot(mp) )
    allocate ( var % wbal(mp) )
    allocate ( var % wbal_tot(mp) )
    allocate ( var % wbtot0(mp) )
    allocate ( var % wetbal(mp) )
  end subroutine alloc_balances_type

  subroutine alloc_soil_parameter_type(var, mp)
    type(soil_parameter_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % albsoil(mp) )
    allocate ( var % bch(mp) )
    allocate ( var % c3(mp) )
    allocate ( var % clay(mp) )
    allocate ( var % cnsd(mp) )
    allocate ( var % css(mp) )
    allocate ( var % froot(mp,ms) )
    allocate ( var % hsbh(mp) )
    allocate ( var % hyds(mp) )
    allocate ( var % i2bp3(mp) )
    allocate ( var % ibp2(mp) )
    allocate ( var % isoilm(mp) )
    allocate ( var % rhosoil(mp) )
    allocate ( var % rs20(mp) )
    allocate ( var % sand(mp) )
    allocate ( var % sfc(mp) )
    allocate ( var % silt(mp) )
    allocate ( var % ssat(mp) )
    allocate ( var % sucs(mp) )
    allocate ( var % swilt(mp) )
  end subroutine alloc_soil_parameter_type
 
  subroutine alloc_soil_snow_type(var, mp)
    type(soil_snow_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % albsoilsn(mp,nrb) )
    allocate ( var % cls(mp) )
    allocate ( var % dfn_dtg(mp) )
    allocate ( var % dfh_dtg(mp) )
    allocate ( var % dfe_ddq(mp) )
    allocate ( var % ddq_dtg(mp) )
    allocate ( var % gammzz(mp,ms) )
    allocate ( var % isflag(mp) )
    allocate ( var % osnowd(mp) )
    allocate ( var % potev(mp) )
    allocate ( var % runoff(mp) )
    allocate ( var % rnof1(mp) )
    allocate ( var % rnof2(mp) )
    allocate ( var % rtsoil(mp) )
    allocate ( var % sdepth(mp,3) )
    allocate ( var % smass(mp,3) )
    allocate ( var % snage(mp) )
    allocate ( var % snowd(mp) )
    allocate ( var % ssdn(mp,3) )
    allocate ( var % ssdnn(mp) )
    allocate ( var % tgg(mp,ms) )
    allocate ( var % tggsn(mp,3) )
    allocate ( var % tss(mp) )
    allocate ( var % wb(mp,ms) )
    allocate ( var % wbfice(mp,ms) )
    allocate ( var % wbice(mp,ms) )
    allocate ( var % wblf(mp,ms) )
    allocate ( var % wbtot(mp) )
  end subroutine alloc_soil_snow_type
   
  subroutine alloc_veg_parameter_type(var, mp)
    type(veg_parameter_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % iveg(mp) )
    allocate ( var % meth(mp) )
    allocate ( var % vlai(mp) )
    allocate ( var % vlaimax(mp) )
    allocate ( var % vlaiw(mp) )
    allocate ( var % fwet(mp) )
    allocate ( var % canst1(mp) )
    allocate ( var % ejmax(mp) )
    allocate ( var % frac4(mp) )
    allocate ( var % tminvj(mp) )
    allocate ( var % tmaxvj(mp) )
    allocate ( var % vbeta(mp) )
    allocate ( var % hc(mp) )
    allocate ( var % shelrb(mp) )
    allocate ( var % vcmax(mp) )
    allocate ( var % xfang(mp) )
    allocate ( var % dleaf(mp) )
    allocate ( var % rp20(mp) )
    allocate ( var % rpcoef(mp) )
  end subroutine alloc_veg_parameter_type
   
  subroutine alloc_canopy_type(var, mp)
    type(canopy_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % cansto(mp) )
    allocate ( var % delwc(mp) )
    allocate ( var % dewmm(mp) )
    allocate ( var % fe(mp) )
    allocate ( var % fh(mp) )
    allocate ( var % fpn(mp) )
    allocate ( var % frp(mp) )
    allocate ( var % frpw(mp) )
    allocate ( var % frpr(mp) )
    allocate ( var % frs(mp) )
    allocate ( var % fnee(mp) )
    allocate ( var % frday(mp) )
    allocate ( var % fnv(mp) )
    allocate ( var % fev(mp) )
    allocate ( var % fevc(mp) )
    allocate ( var % fevw(mp) )
    allocate ( var % fhv(mp) )
    allocate ( var % fhvw(mp) )
    allocate ( var % fns(mp) )
    allocate ( var % fes(mp) )
    allocate ( var % fhs(mp) )
    allocate ( var % tv(mp) )
    allocate ( var % ga(mp) )
    allocate ( var % ghflux(mp) )
    allocate ( var % sghflux(mp) )
    allocate ( var % dgdtg(mp) )
    allocate ( var % through(mp) )
    allocate ( var % precis(mp) )
    allocate ( var % rnet(mp) )
    allocate ( var % spill(mp) )
    allocate ( var % wcint(mp) )
    allocate ( var % us(mp) )
    allocate ( var % tscrn(mp) )
    allocate ( var % qscrn(mp) )
    allocate ( var % uscrn(mp) )
    allocate ( var % cduv(mp) )
  end subroutine alloc_canopy_type
   
  subroutine alloc_radiation_type(var, mp)
    type(radiation_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % albedo(mp,nrb) )
    allocate ( var % extkb(mp) )
    allocate ( var % extkd2(mp) )
    allocate ( var % extkd(mp) )
    allocate ( var % extkn(mp) )
    allocate ( var % flws(mp) )
    allocate ( var % fvlai(mp,mf) )
    allocate ( var % gradis(mp,mf) )
    allocate ( var % latitude(mp) )
    allocate ( var % lwabv(mp) )
    allocate ( var % qcan(mp,mf,nrb) )
    allocate ( var % qssabs(mp) )
    allocate ( var % rhocdf(mp,nrb) )
    allocate ( var % rniso(mp,mf) )
    allocate ( var % scalex(mp,mf) )
    allocate ( var % transd(mp) )
    allocate ( var % trad(mp) )
  end subroutine alloc_radiation_type
   
  subroutine alloc_roughness_type(var, mp)
    type(roughness_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % coexp(mp) )
    allocate ( var % disp(mp) )
    allocate ( var % rt0us(mp) )
    allocate ( var % rt1usa(mp) )
    allocate ( var % rt1usb(mp) )
    allocate ( var % rt1(mp) )
    allocate ( var % rt1(mp) )
    allocate ( var % usuh(mp) )
    allocate ( var % za(mp) )
    allocate ( var % z0m(mp) )
    allocate ( var % zref(mp) )
    allocate ( var % zruffs(mp) )
    allocate ( var % z0soilsn(mp) )
  end subroutine alloc_roughness_type
   
  subroutine alloc_air_type(var, mp)
    type(air_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % rho(mp) )
    allocate ( var % volm(mp) )
    allocate ( var % rlam(mp) )
    allocate ( var % qsat(mp) )
    allocate ( var % epsi(mp) )
    allocate ( var % visc(mp) )
    allocate ( var % psyc(mp) )
    allocate ( var % dsatdk(mp) )
    allocate ( var % cmolar(mp) )
  end subroutine alloc_air_type
   
  subroutine alloc_met_type(var, mp)
    type(met_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % ca(mp) )
    allocate ( var % year(mp) )
    allocate ( var % moy(mp) )
    allocate ( var % doy(mp) )
    allocate ( var % hod(mp) )
    allocate ( var % fsd(mp) )
    allocate ( var % fld(mp) )
    allocate ( var % precip(mp) )
    allocate ( var % precips(mp) )
    allocate ( var % tc(mp) )
    allocate ( var % tk(mp) )
    allocate ( var % tvair(mp) )
    allocate ( var % tvrad(mp) )
    allocate ( var % pmb(mp) )
    allocate ( var % ua(mp) )
    allocate ( var % qv(mp) )
    allocate ( var % qvair(mp) )
    allocate ( var % da(mp) )
    allocate ( var % dva(mp) )
    allocate ( var % coszen(mp) )
  end subroutine alloc_met_type
   
  subroutine alloc_sum_flux_type(var, mp)
    type(sum_flux_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % sumpn(mp) )
    allocate ( var % sumrp(mp) )
    allocate ( var % sumrpw(mp) )
    allocate ( var % sumrpr(mp) )
    allocate ( var % sumrs(mp) )
    allocate ( var % sumrd(mp) )
    allocate ( var % dsumpn(mp) )
    allocate ( var % dsumrp(mp) )
    allocate ( var % dsumrs(mp) )
    allocate ( var % dsumrd(mp) )
    allocate ( var % sumxrp(mp) )
    allocate ( var % sumxrs(mp) )
  end subroutine alloc_sum_flux_type

  subroutine alloc_bgc_pool_type(var, mp)
    type(bgc_pool_type), intent(out) :: var
    integer, intent(in) :: mp
    allocate ( var % cplant(mp,ncp) )
    allocate ( var % csoil(mp,ncs) )
  end subroutine alloc_bgc_pool_type
  
END MODULE define_types
!=========================================================================
MODULE math_constants
  USE define_dimensions, only : i_d, r_1
  REAL(r_1), PARAMETER	:: pi = 3.1415927
  REAL(r_1), PARAMETER	:: pi180 = pi / 180.0 ! radians / degree
  REAL(r_1), PARAMETER	:: two_pi = 2.0 * pi
END MODULE math_constants
!=========================================================================
MODULE physical_constants
  USE define_dimensions, only : i_d, r_1
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
  REAL(r_1), PARAMETER :: z0soil = 1.0e-6 ! roughness length of bare soil surface
END MODULE physical_constants
!=========================================================================
MODULE other_constants
  USE define_dimensions, only : i_d, r_1, nrb
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
  USE define_dimensions, only : i_d, r_1
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
