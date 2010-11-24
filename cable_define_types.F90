MODULE define_types
  ! Contains all variables which are not subroutine-internal
  USE define_dimensions
  ! Energy and water balance variables:
  TYPE balances_type 
!!les     REAL(r_1), DIMENSION(:), POINTER :: cansto0 ! canopy water storage (mm)
     REAL(r_1), DIMENSION(:), POINTER :: drybal ! energy balance for dry canopy
     REAL(r_1), DIMENSION(:), POINTER :: ebal   ! energy balance per time step (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER :: ebal_tot ! cumulative energy balance (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER :: evap_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: evapc_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: evaps_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: osnowd0  ! snow depth, first time step
     REAL(r_1), DIMENSION(:), POINTER :: precip_tot ! cumulative precipitation (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnoff_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof1_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof2_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: snowdc_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbal   ! water balance per time step (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbal_tot ! cumulative water balance (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbal_tot1 ! cumulative water balance (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: wbtot0 ! total soil water (mm), first time step
     REAL(r_1), DIMENSION(:), POINTER :: owbtot ! total soil water (mm), first time step
     REAL(r_1), DIMENSION(:), POINTER :: wetbal ! energy balance for wet canopy
     REAL(r_1), DIMENSION(:), POINTER :: delwc_tot ! energy balance for wet canopy
     REAL(r_1), DIMENSION(:), POINTER :: qasrf_tot ! heat advected to the snow by precip. 
     REAL(r_1), DIMENSION(:), POINTER :: qfsrf_tot ! energy of snowpack phase changes 
     REAL(r_1), DIMENSION(:), POINTER :: qssrf_tot ! energy of snowpack phase changes 
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type 
     REAL(r_1), DIMENSION(:), POINTER :: albsoil ! soil reflectance
     REAL(r_1), DIMENSION(:), POINTER :: bch  ! parameter b in Campbell equation
     REAL(r_1), DIMENSION(:), POINTER :: c3   ! c3 drainage coeff (fraction)
     REAL(r_1), DIMENSION(:), POINTER :: clay ! fraction of soil which is clay
     REAL(r_2), DIMENSION(:), POINTER :: cnsd ! thermal conductivity of dry soil [W/m/K]
!les     REAL(r_1), DIMENSION(:), POINTER :: cnsdJ ! thermal conductivity of dry soil [W/m/K]
     REAL(r_1), DIMENSION(:), POINTER :: css  ! soil specific heat capacity [kJ/kg/K]
!les: move froot to veg param
!     REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL(r_1), DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER(i_d), DIMENSION(:), POINTER :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
     INTEGER(i_d), DIMENSION(:), POINTER :: ibp2   ! parameter two in K vis suction (function of pbch)
     INTEGER(i_d), DIMENSION(:), POINTER :: isoilm ! integer soil type
     REAL(r_2), DIMENSION(:), POINTER :: pwb_min ! working variable (swilt/ssat)**ibp2
 
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
     REAL(r_1), DIMENSION(:), POINTER :: fwtop   ! water flux to the soil
     REAL(r_2), DIMENSION(:,:), POINTER :: gammzz ! heat capacity for each soil layer
     INTEGER(i_d), DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow
     REAL(r_1), DIMENSION(:), POINTER :: potev   ! potential evapotranspiration
     REAL(r_1), DIMENSION(:), POINTER :: qasrf ! heat advected to the snow by precip. 
     REAL(r_1), DIMENSION(:), POINTER :: qfsrf ! energy of snowpack phase changes 
     REAL(r_1), DIMENSION(:), POINTER :: qssrf ! sublimation 
     REAL(r_1), DIMENSION(:), POINTER :: runoff  ! total runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof1   ! surface runoff (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rnof2   ! deep drainage (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: rtsoil  ! turbulent resistance for soil
     REAL(r_1), DIMENSION(:,:), POINTER :: sconds !
     REAL(r_1), DIMENSION(:,:), POINTER :: sdepth ! snow depth
     REAL(r_1), DIMENSION(:,:), POINTER  :: smass  ! snow mass
     REAL(r_1), DIMENSION(:), POINTER :: snage   ! snow age
     REAL(r_1), DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
     REAL(r_1), DIMENSION(:), POINTER :: osnowd  ! snow depth from previous time step
!les: Making wetfac global
     REAL(r_1), DIMENSION(:), POINTER :: wetfac ! surface wetness fact. at current time step
     REAL(r_1), DIMENSION(:), POINTER :: owetfac ! surface wetness fact. at previous time step

     REAL(r_1), DIMENSION(:), POINTER :: smelt   ! snow melt 
!     REAL(r_1), DIMENSION(:,0:3), POINTER :: smelt1 ! snow melt (3 layer snowpack) 
     REAL(r_1), DIMENSION(:,:), POINTER  :: ssdn ! snow densities
     REAL(r_1), DIMENSION(:), POINTER :: ssdnn   ! average snow density
     REAL(r_1), DIMENSION(:,:), POINTER :: tgg   ! soil temperature in K
!!les     REAL(r_1), DIMENSION(:), POINTER :: otgg   ! soil temperature in K
     REAL(r_1), DIMENSION(:,:), POINTER  :: tggsn ! snow temperature in K
     REAL(r_1), DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
     REAL(r_1), DIMENSION(:), POINTER :: otss     ! surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(:,:), POINTER :: wb    ! volumetric soil moisture (solid+liq)
     REAL(r_2), DIMENSION(:,:), POINTER :: wbfice !
     REAL(r_2), DIMENSION(:,:), POINTER :: wbice  ! volumetric soil ice
     REAL(r_2), DIMENSION(:,:), POINTER :: wblf !
     REAL(r_1), DIMENSION(:), POINTER :: wbtot   ! total soil water (mm)
!!les     REAL(r_1), DIMENSION(:), POINTER :: tprecip
!!     REAL(r_1), DIMENSION(:), POINTER :: tevap
!!     REAL(r_1), DIMENSION(:), POINTER :: trnoff
!!     REAL(r_1), DIMENSION(:), POINTER :: totenbal
!!     REAL(r_1), DIMENSION(:), POINTER :: totenbal2
!sxy
    REAL(r_1), DIMENSION(:), POINTER :: evapsn  ! snow evaporation  
  END TYPE soil_snow_type
  ! Vegetation parameters:
  TYPE veg_parameter_type
     INTEGER(i_d),DIMENSION(:), POINTER :: iveg ! vegetation type
     INTEGER(i_d),DIMENSION(:), POINTER :: meth ! method for calculation of canopy fluxes and temp.
     REAL(r_1), DIMENSION(:), POINTER :: vlai   ! leaf area index
     REAL(r_1), DIMENSION(:), POINTER :: vlaimax ! ???
!les: move vlaiw to canopy
!     REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of resistances
     REAL(r_1), DIMENSION(:,:), POINTER :: refl
     REAL(r_1), DIMENSION(:,:), POINTER :: taul
     REAL(r_1), DIMENSION(:), POINTER :: xalbnir ! modifier for albedo in near ir bnad (YP Apr08)

!les: moved to canopy_type
!     REAL(r_1), DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
     REAL(r_1), DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
     REAL(r_1), DIMENSION(:), POINTER :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
     REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: vbeta  ! 
!     REAL(r_1), DIMENSION(:), POINTER :: hc_grdmax ! maximum height of canopy from tiles belonging to the same grid
     REAL(r_1), DIMENSION(:), POINTER :: hc	! roughness height of canopy (veg - snow)
     REAL(r_1), DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless)
!les wai,vegcf,xalbnir
     REAL(r_1), DIMENSION(:), POINTER :: wai    ! wood area index (stem+branches+twigs)
     REAL(r_1), DIMENSION(:), POINTER :: vegcf  ! biome-specific soil respiration rate
     REAL(r_1), DIMENSION(:), POINTER :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
     REAL(r_1), DIMENSION(:), POINTER :: dleaf  ! chararacteristc legnth of leaf (m)
     REAL(r_1), DIMENSION(:), POINTER :: rp20   ! plant respiration coefficient at 20 C
     REAL(r_1), DIMENSION(:), POINTER :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
!les
!    rml 22/10/07
     LOGICAL,   DIMENSION(:), POINTER :: deciduous ! flag used for phenology fix
  END TYPE veg_parameter_type
  ! Canopy/vegetation variables:
  TYPE canopy_type
     REAL(r_1), DIMENSION(:), POINTER :: cansto ! canopy water storage (mm)
     REAL(r_1), DIMENSION(:), POINTER :: oldcansto ! canopy water storage (mm)
     REAL(r_1), DIMENSION(:), POINTER :: delwc ! change in canopy water store (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: dewmm ! dewfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: epot  ! total potential evaporation 
     REAL(r_1), DIMENSION(:), POINTER :: fe    ! total latent heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fh    ! total sensible heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fpn   ! plant photosynthesis (g C s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frp   ! plant respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frpw  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: frpr  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: fnpp  ! npp flux
     REAL(r_1), DIMENSION(:), POINTER :: frs   ! soil respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: fnee  ! net carbon flux (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frday ! daytime leaf resp
     REAL(r_1), DIMENSION(:), POINTER :: fnv   ! net rad. avail. to canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fev   ! latent hf from canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fevc  ! dry canopy transpiration (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fevw  ! lat heat fl wet canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fevw_pot ! potential lat heat from canopy
     REAL(r_1), DIMENSION(:), POINTER :: fhv   ! sens heatfl from canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fhvw  ! sens heatfl from wet canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fns   ! net rad avail to soil (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fes   ! latent heatfl from soil (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: segg  ! latent heatfl from soil mm
     REAL(r_1), DIMENSION(:), POINTER :: fhs   ! sensible heat flux from soil
!
     REAL(r_1), DIMENSION(:), POINTER :: fwet  ! fraction of canopy wet
     REAL(r_1), DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
     REAL(r_1), DIMENSION(:,:), POINTER :: gswx ! ! stom cond for water
     REAL(r_1), DIMENSION(:), POINTER :: gswx_T ! ! stom cond for water
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
     REAL(r_1), DIMENSION(:), POINTER :: ua_10m	 ! screen level or 10m wind speed (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: tscrn ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tscrn1 ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tscrn2 ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tscrn3 ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: qscrn ! specific humudity at screen height (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: uscrn ! wind speed at screen height (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: uscrn1 ! wind speed at screen height (m/s)
!les
     REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of resistances
     REAL(r_1), DIMENSION(:), POINTER :: cduv  ! drag coefficient for momentum
     REAL(r_1), DIMENSION(:), POINTER :: cdtq  ! drag coefficient for heat
     REAL(r_1), DIMENSION(:), POINTER :: wetfac_cs ! 
     REAL(r_1), DIMENSION(:,:), POINTER :: zetar ! stability correction

  END TYPE canopy_type
  ! Radiation variables:
  TYPE radiation_type
     REAL(r_1), DIMENSION(:,:), POINTER :: albedo ! canopy+soil albedo
     REAL(r_1), DIMENSION(:), POINTER :: albedo_T ! canopy+soil albedo for VIS+NIR
     REAL(r_1), DIMENSION(:), POINTER     :: extkb  ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER     :: extkd2 ! diffuse 2D radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
     REAL(r_1), DIMENSION(:), POINTER ::  extkn ! extinction coef for vertical nitrogen profile in canopy(-)
     REAL(r_1), DIMENSION(:), POINTER :: flws   ! soil long-wave radiation
     REAL(r_1), DIMENSION(:,:), POINTER  :: fvlai  ! leaf area index of big leaf
     REAL(r_1), DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
     REAL(r_1), DIMENSION(:), POINTER :: latitude  ! latitude
     REAL(r_1), DIMENSION(:), POINTER :: longitude ! longitude
     REAL(r_1), DIMENSION(:), POINTER :: lwabv ! long wave absorbed by vegetation
     REAL(r_1), DIMENSION(:,:,:), POINTER :: qcan ! absorbed radiation for canopy (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER     :: qssabs ! absorbed short-wave radiation for soil
     REAL(r_1), DIMENSION(:,:), POINTER ::  rhocdf ! canopy diffuse reflectance (-)
     REAL(r_1), DIMENSION(:,:), POINTER  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
     REAL(r_1), DIMENSION(:,:), POINTER  :: scalex ! scaling PARAMETER for big leaf
     REAL(r_1), DIMENSION(:), POINTER     :: transd ! fraction SW diffuse transmitted through canopy
     REAL(r_1), DIMENSION(:), POINTER ::  trad  !  radiative temperature (soil and veg)
!===sxy
     REAL(r_1),DIMENSION(:,:), POINTER  :: reffdf  !effective conopy diffuse reflectance
     REAL(r_1),DIMENSION(:,:), POINTER  :: reffbm  !effective conopy beam reflectance
     REAL(r_1), DIMENSION(:,:), POINTER :: rhocbm  !modified canopy beam reflectance(6.21)
     REAL(r_1), DIMENSION(:,:), POINTER :: extkbm  !modified k beam(6.20)(for leaf scattering)
     REAL(r_1), DIMENSION(:,:), POINTER :: extkdm  !modified k diffuse(6.20)(for leaf scattering)
     REAL(r_1), DIMENSION(:,:), POINTER   :: fbeam   !beam fraction 
     REAL(r_1), DIMENSION(:,:), POINTER   :: cexpkbm ! canopy beam transmittance
     REAL(r_1), DIMENSION(:,:), POINTER   :: cexpkdm ! canopy diffuse transmittance
     REAL(r_1), DIMENSION(:), POINTER   :: transb  ! fraction SW beam tranmitted through canopy
!===sxy

  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL(r_1), DIMENSION(:), POINTER	:: coexp ! Extinction coefficient for wind profile in canopy
     REAL(r_1), DIMENSION(:), POINTER   :: hruff_grmx ! maximum height of canopy from tiles belonging to the same grid
     REAL(r_1), DIMENSION(:), POINTER	:: disp  ! zero-plane displacement
     REAL(r_1), DIMENSION(:), POINTER	:: hruff ! canopy height above snow level
     REAL(r_1), DIMENSION(:), POINTER	:: rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usa ! resistance from disp to hruf
     REAL(r_1), DIMENSION(:), POINTER	:: rt1usb ! resistance from hruf to zruffs (or zref if zref<zruffs)
     REAL(r_1), DIMENSION(:), POINTER	:: rt1 ! 1/aerodynamic conductance
     REAL(r_1), DIMENSION(:), POINTER	:: term2, term3, term5, term6 ! for aerodynamic resistance calc.
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL(r_1), DIMENSION(:), POINTER	:: usuh ! Friction velocity/windspeed at canopy height
     REAL(r_1), DIMENSION(:), POINTER	:: za_uv   ! level of lowest atmospheric model layer
     REAL(r_1), DIMENSION(:), POINTER	:: za_tq   ! level of lowest atmospheric model layer
     REAL(r_1), DIMENSION(:), POINTER	:: z0m  ! roughness length
     REAL(r_1), DIMENSION(:), POINTER	:: zref_uv ! Reference height for met forcing
     REAL(r_1), DIMENSION(:), POINTER	:: zref_tq ! Reference height for met forcing
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
     REAL(r_1), DIMENSION(:,:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: precip_s ! solid preipitation only (mm/dels)
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
!les
!!  TYPE grid_type
!!     REAL(r_1), DIMENSION(:,:), POINTER :: cplant_g ! plant carbon (g C/m2))
!!     REAL(r_1), DIMENSION(:,:), POINTER :: csoil_g  ! soil carbon (g C/m2)
!!     REAL(r_1), DIMENSION(:,:), POINTER :: frac_g  ! tile fractions 
!!     REAL(r_1), DIMENSION(:), POINTER :: fe_g    ! total latent heat (W/m2)
!!     REAL(r_1), DIMENSION(:), POINTER :: fh_g    ! total sensible heat (W/m2)
!!  END TYPE grid_type

!les
  TYPE (air_type)       :: air ! air property variables
  TYPE (bgc_pool_type)  :: bgc ! carbon pool variables
  TYPE (canopy_type)    :: canopy ! vegetation variables
  TYPE (met_type)       :: met ! met input variables
  TYPE (balances_type)  :: bal ! energy and water balance variables
  TYPE (radiation_type) :: rad ! radiation variables
  TYPE (roughness_type) :: rough ! roughness variables
  TYPE (soil_parameter_type) :: soil ! soil parameters
  TYPE (soil_snow_type) :: ssoil ! soil and snow variables
  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables 
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters

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
          alloc_bgc_pool_type!,            &
!!les          alloc_grid_type
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
          dealloc_bgc_pool_type!,            &
!!les          dealloc_grid_type
  END INTERFACE
CONTAINS
  
  SUBROUTINE alloc_balances_type(var, mp)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
!!les    ALLOCATE ( var % cansto0(mp) )
    ALLOCATE ( var % drybal(mp) )
    ALLOCATE ( var % ebal(mp) )
    ALLOCATE ( var % ebal_tot(mp) )
    ALLOCATE ( var % evap_tot(mp) )
    ALLOCATE ( var % evapc_tot(mp) )
    ALLOCATE ( var % evaps_tot(mp) )
    ALLOCATE ( var % delwc_tot(mp) )
    ALLOCATE ( var % qasrf_tot(mp) )
    ALLOCATE ( var % qfsrf_tot(mp) )
    ALLOCATE ( var % qssrf_tot(mp) )
    ALLOCATE ( var % osnowd0(mp) )
    ALLOCATE ( var % precip_tot(mp) )
    ALLOCATE ( var % rnoff_tot(mp) )
    ALLOCATE ( var % rnof1_tot(mp) )
    ALLOCATE ( var % rnof2_tot(mp) )
    ALLOCATE ( var % snowdc_tot(mp) )
    ALLOCATE ( var % wbal(mp) )
    ALLOCATE ( var % wbal_tot(mp) )
    ALLOCATE ( var % wbal_tot1(mp) )
    ALLOCATE ( var % wbtot0(mp) )
    ALLOCATE ( var % owbtot(mp) )
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
!les    ALLOCATE ( var % cnsdJ(mp) )
    ALLOCATE ( var % css(mp) )
!    ALLOCATE ( var % froot(mp,ms) )
    ALLOCATE ( var % hsbh(mp) )
    ALLOCATE ( var % hyds(mp) )
    ALLOCATE ( var % i2bp3(mp) )
    ALLOCATE ( var % ibp2(mp) )
    ALLOCATE ( var % isoilm(mp) )
    ALLOCATE ( var % pwb_min(mp) )
    ALLOCATE ( var % rhosoil(mp) )
    ALLOCATE ( var % rs20(mp) )
    ALLOCATE ( var % sand(mp) )
    ALLOCATE ( var % sfc(mp) )
    ALLOCATE ( var % silt(mp) )
    ALLOCATE ( var % ssat(mp) )
    ALLOCATE ( var % sucs(mp) )
    ALLOCATE ( var % swilt(mp) )
   print *,'testing cable_define_types.F9a v674'
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
!les
    ALLOCATE ( var % wetfac(mp) )
    ALLOCATE ( var % owetfac(mp) )
    ALLOCATE ( var % potev(mp) )
    ALLOCATE ( var % runoff(mp) )
    ALLOCATE ( var % rnof1(mp) )
    ALLOCATE ( var % rnof2(mp) )
    ALLOCATE ( var % rtsoil(mp) )
    ALLOCATE ( var % sdepth(mp,3) )
    ALLOCATE ( var % smass(mp,3) )
    ALLOCATE ( var % sconds(mp,3) )
    ALLOCATE ( var % snage(mp) )
    ALLOCATE ( var % snowd(mp) )
    ALLOCATE ( var % qasrf(mp) )
    ALLOCATE ( var % qfsrf(mp) )
    ALLOCATE ( var % qssrf(mp) )
    ALLOCATE ( var % smelt(mp) )
    ALLOCATE ( var % fwtop(mp) )
    ALLOCATE ( var % ssdn(mp,3) )
    ALLOCATE ( var % ssdnn(mp) )
    ALLOCATE ( var % tgg(mp,ms) )
!!les    ALLOCATE ( var % otgg(mp) )
    ALLOCATE ( var % tggsn(mp,3) )
    ALLOCATE ( var % tss(mp) )
    ALLOCATE ( var % otss(mp) )
    ALLOCATE ( var % wb(mp,ms) )
    ALLOCATE ( var % wbfice(mp,ms) )
    ALLOCATE ( var % wbice(mp,ms) )
    ALLOCATE ( var % wblf(mp,ms) )
    ALLOCATE ( var % wbtot(mp) )
    ALLOCATE ( var % evapsn(mp))
  END SUBROUTINE alloc_soil_snow_type
   
  SUBROUTINE alloc_veg_parameter_type(var, mp)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % iveg(mp) )
    ALLOCATE ( var % meth(mp) )
    ALLOCATE ( var % vlai(mp) )
    ALLOCATE ( var % vlaimax(mp) )
!    ALLOCATE ( var % vlaiw(mp) )
!    ALLOCATE ( var % fwet(mp) )
    ALLOCATE ( var % refl(mp,2) )
    ALLOCATE ( var % taul(mp,2) )
    ALLOCATE ( var % xalbnir(mp) )
    ALLOCATE ( var % canst1(mp) )
    ALLOCATE ( var % ejmax(mp) )
    ALLOCATE ( var % frac4(mp) )
    ALLOCATE ( var % froot(mp,ms) )
!les
    ALLOCATE ( var % wai(mp) )     ! new addition in Oct 2007 (YP)
    ALLOCATE ( var % vegcf(mp) )   ! new addition in Oct 2007 (YP) 
    ALLOCATE ( var % xalbnir(mp) ) ! new addition in Apr 2008 (YP)
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
    ALLOCATE ( var % deciduous(mp) ) ! les: rml addition 22/10/07
  END SUBROUTINE alloc_veg_parameter_type
   
  SUBROUTINE alloc_canopy_type(var, mp)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % cansto(mp) )
    ALLOCATE ( var % oldcansto(mp) )
    ALLOCATE ( var % delwc(mp) )
    ALLOCATE ( var % dewmm(mp) )
    ALLOCATE ( var % fe(mp) )
    ALLOCATE ( var % epot(mp) )
    ALLOCATE ( var % fh(mp) )
    ALLOCATE ( var % fpn(mp) )
    ALLOCATE ( var % frp(mp) )
    ALLOCATE ( var % fnpp(mp) )
    ALLOCATE ( var % frpw(mp) )
    ALLOCATE ( var % frpr(mp) )
    ALLOCATE ( var % frs(mp) )
    ALLOCATE ( var % fnee(mp) )
    ALLOCATE ( var % frday(mp) )
    ALLOCATE ( var % fnv(mp) )
    ALLOCATE ( var % fev(mp) )
    ALLOCATE ( var % fevc(mp) )
    ALLOCATE ( var % fevw(mp) )
    ALLOCATE ( var % fevw_pot(mp) )
    ALLOCATE ( var % fhv(mp) )
    ALLOCATE ( var % fhvw(mp) )
    ALLOCATE ( var % fns(mp) )
    ALLOCATE ( var % fes(mp) )
    ALLOCATE ( var % segg(mp) )
    ALLOCATE ( var % fhs(mp) )
!les
    ALLOCATE ( var % fwet(mp) )
    ALLOCATE ( var % tv(mp) )
    ALLOCATE ( var % gswx(mp,mf) )
    ALLOCATE ( var % gswx_T(mp) )
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
    ALLOCATE ( var % ua_10m(mp) )
    ALLOCATE ( var % tscrn(mp) )
    ALLOCATE ( var % tscrn1(mp) )
    ALLOCATE ( var % tscrn2(mp) )
    ALLOCATE ( var % tscrn3(mp) )
    ALLOCATE ( var % qscrn(mp) )
    ALLOCATE ( var % uscrn(mp) )
    ALLOCATE ( var % uscrn1(mp) )
!les
    ALLOCATE ( var % vlaiw(mp) )
    ALLOCATE ( var % cduv(mp) )
    ALLOCATE ( var % cdtq(mp) )
    ALLOCATE ( var % wetfac_cs(mp) )
    ALLOCATE ( var % zetar(mp,niter) )

  END SUBROUTINE alloc_canopy_type
   
  SUBROUTINE alloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % albedo(mp,nrb) )
    ALLOCATE ( var % albedo_T(mp) )
    ALLOCATE ( var % extkb(mp) )
    ALLOCATE ( var % extkd2(mp) )
    ALLOCATE ( var % extkd(mp) )
    ALLOCATE ( var % extkn(mp) )
    ALLOCATE ( var % flws(mp) )
    ALLOCATE ( var % fvlai(mp,mf) )
    ALLOCATE ( var % gradis(mp,mf) )
    ALLOCATE ( var % latitude(mp) )
    ALLOCATE ( var % longitude(mp) )
    ALLOCATE ( var % lwabv(mp) )
    ALLOCATE ( var % qcan(mp,mf,nrb) )
    ALLOCATE ( var % qssabs(mp) )
    ALLOCATE ( var % rhocdf(mp,nrb) )
    ALLOCATE ( var % rniso(mp,mf) )
    ALLOCATE ( var % scalex(mp,mf) )
    ALLOCATE ( var % transd(mp) )
    ALLOCATE ( var % trad(mp) )
!==sxy
    ALLOCATE ( var % reffdf(mp,nrb))
    ALLOCATE ( var % reffbm(mp,nrb))
    ALLOCATE ( var % rhocbm(mp,nrb))
    ALLOCATE ( var % extkbm(mp,nrb))
    ALLOCATE ( var % extkdm(mp,nrb))
!les: remove hard wiring
    ALLOCATE ( var % fbeam (mp,nrb))  ! needs 3 in second dimension
    ALLOCATE ( var % cexpkbm(mp,nrb)) ! only needs 2 in the second dimension
    ALLOCATE ( var % cexpkdm(mp,nrb)) ! only needs 2 in the second dimension
!    ALLOCATE ( var % fbeam (mp,3))
!    ALLOCATE ( var % cexpkbm(mp,2))
!    ALLOCATE ( var % cexpkdm(mp,2))
    ALLOCATE ( var % transb (mp))
!==sxy
  END SUBROUTINE alloc_radiation_type
   
  SUBROUTINE alloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % hruff_grmx(mp) )
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
    ALLOCATE ( var % za_uv(mp) )
    ALLOCATE ( var % za_tq(mp) )
    ALLOCATE ( var % z0m(mp) )
    ALLOCATE ( var % zref_uv(mp) )
    ALLOCATE ( var % zref_tq(mp) )
    ALLOCATE ( var % zruffs(mp) )
    ALLOCATE ( var % z0soilsn(mp) )
    ALLOCATE ( var % z0soil(mp) )
  END SUBROUTINE alloc_roughness_type
   
  SUBROUTINE alloc_air_type(var, mp)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    print*,'mp==',mp
    ALLOCATE ( var % rho(mp) )
    print*,'alloc_air_type2',mp
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
!les: remove hardwiring
    ALLOCATE ( var % fsd(mp,nrb) )
!    ALLOCATE ( var % fsd(mp,3) )
    ALLOCATE ( var % fld(mp) )
    ALLOCATE ( var % precip(mp) )
    ALLOCATE ( var % precip_s(mp) )
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

!les
!!  SUBROUTINE alloc_grid_type(var, mg)
!!    TYPE(grid_type), INTENT(inout) :: var
!!    INTEGER, INTENT(in) :: mg
!!    ALLOCATE ( var % frac_g(mg,12) )
!!    ALLOCATE ( var % csoil_g(mg,ncs) )
!!    ALLOCATE ( var % fe_g(mg) )
!!    ALLOCATE ( var % fh_g(mg) )
!!  END SUBROUTINE alloc_grid_type

  ! Begin deallocation routines:
   SUBROUTINE dealloc_balances_type(var, mp)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
!!les    DEALLOCATE ( var % cansto0 )
    DEALLOCATE ( var % drybal )
    DEALLOCATE ( var % ebal )
    DEALLOCATE ( var % ebal_tot )
    DEALLOCATE ( var % evap_tot )
    DEALLOCATE ( var % evapc_tot )
    DEALLOCATE ( var % evaps_tot )
    DEALLOCATE ( var % delwc_tot )
    DEALLOCATE ( var % osnowd0 )
    DEALLOCATE ( var % precip_tot )
    DEALLOCATE ( var % rnoff_tot )
    DEALLOCATE ( var % rnof1_tot )
    DEALLOCATE ( var % rnof2_tot )
    DEALLOCATE ( var % snowdc_tot )
    DEALLOCATE ( var % wbal )
    DEALLOCATE ( var % wbal_tot )
    DEALLOCATE ( var % wbal_tot1 )
    DEALLOCATE ( var % wbtot0 )
    DEALLOCATE ( var % owbtot )
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
!les    DEALLOCATE ( var % cnsdJ )
    DEALLOCATE ( var % css )
!    DEALLOCATE ( var % froot )
    DEALLOCATE ( var % hsbh )
    DEALLOCATE ( var % hyds )
    DEALLOCATE ( var % i2bp3 )
    DEALLOCATE ( var % ibp2 )
    DEALLOCATE ( var % isoilm )
    DEALLOCATE ( var % pwb_min )
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
!les
    DEALLOCATE ( var % wetfac )
    DEALLOCATE ( var % owetfac )
    DEALLOCATE ( var % potev )
    DEALLOCATE ( var % runoff )
    DEALLOCATE ( var % rnof1 )
    DEALLOCATE ( var % rnof2 )
    DEALLOCATE ( var % rtsoil )
    DEALLOCATE ( var % sdepth )
    DEALLOCATE ( var % smass )
    DEALLOCATE ( var % sconds )
    DEALLOCATE ( var % snage )
    DEALLOCATE ( var % snowd )
    DEALLOCATE ( var % qasrf )
    DEALLOCATE ( var % qfsrf )
    DEALLOCATE ( var % qssrf )
    DEALLOCATE ( var % smelt )
    DEALLOCATE ( var % fwtop )
    DEALLOCATE ( var % ssdn )
    DEALLOCATE ( var % ssdnn )
    DEALLOCATE ( var % tgg )
!!les    DEALLOCATE ( var % otgg )
    DEALLOCATE ( var % tggsn )
    DEALLOCATE ( var % tss )
    DEALLOCATE ( var % otss )
    DEALLOCATE ( var % wb )
    DEALLOCATE ( var % wbfice )
    DEALLOCATE ( var % wbice )
    DEALLOCATE ( var % wblf )
    DEALLOCATE ( var % wbtot )
    DEALLOCATE ( var % evapsn)
  END SUBROUTINE dealloc_soil_snow_type
   
  SUBROUTINE dealloc_veg_parameter_type(var, mp)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % iveg )
    DEALLOCATE ( var % meth )
    DEALLOCATE ( var % vlai )
    DEALLOCATE ( var % vlaimax )
!    DEALLOCATE ( var % vlaiw )
!    DEALLOCATE ( var % fwet )
    DEALLOCATE ( var % refl )
    DEALLOCATE ( var % taul )
    DEALLOCATE ( var % xalbnir )
    DEALLOCATE ( var % canst1 )
    DEALLOCATE ( var % ejmax )
    DEALLOCATE ( var % frac4 )
    DEALLOCATE ( var % froot )
!les
    DEALLOCATE ( var % wai )     ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % vegcf )   ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % xalbnir ) ! new addition in Apr 2008 (YP)
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
    DEALLOCATE ( var % deciduous ) ! les: rml addition 22/10/07
  END SUBROUTINE dealloc_veg_parameter_type
   
  SUBROUTINE dealloc_canopy_type(var, mp)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % cansto )
    DEALLOCATE ( var % oldcansto )
    DEALLOCATE ( var % delwc )
    DEALLOCATE ( var % dewmm )
    DEALLOCATE ( var % fe )
    DEALLOCATE ( var % epot )
    DEALLOCATE ( var % fh )
    DEALLOCATE ( var % fpn )
    DEALLOCATE ( var % fnpp )
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
    DEALLOCATE ( var % fevw_pot )
    DEALLOCATE ( var % fhv )
    DEALLOCATE ( var % fhvw )
    DEALLOCATE ( var % fns )
    DEALLOCATE ( var % fes )
    DEALLOCATE ( var % fhs )
!les
    DEALLOCATE ( var % fwet )
    DEALLOCATE ( var % tv )
    DEALLOCATE ( var % gswx )
    DEALLOCATE ( var % gswx_T )
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
    DEALLOCATE ( var % ua_10m )
    DEALLOCATE ( var % tscrn )
    DEALLOCATE ( var % tscrn1 )
    DEALLOCATE ( var % tscrn2 )
    DEALLOCATE ( var % tscrn3 )
    DEALLOCATE ( var % qscrn )
    DEALLOCATE ( var % uscrn )
    DEALLOCATE ( var % uscrn1 )
!les
    DEALLOCATE ( var % vlaiw )
    DEALLOCATE ( var % cduv )
    DEALLOCATE ( var % cdtq )
    DEALLOCATE ( var % wetfac_cs )
    DEALLOCATE ( var % zetar )

  END SUBROUTINE dealloc_canopy_type
   
  SUBROUTINE dealloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % albedo )
    DEALLOCATE ( var % albedo_T )
    DEALLOCATE ( var % extkb )
    DEALLOCATE ( var % extkd2 )
    DEALLOCATE ( var % extkd )
    DEALLOCATE ( var % extkn )
    DEALLOCATE ( var % flws )
    DEALLOCATE ( var % fvlai )
    DEALLOCATE ( var % gradis )
    DEALLOCATE ( var % latitude )
    DEALLOCATE ( var % longitude )
    DEALLOCATE ( var % lwabv )
    DEALLOCATE ( var % qcan )
    DEALLOCATE ( var % qssabs )
    DEALLOCATE ( var % rhocdf )
    DEALLOCATE ( var % rniso )
    DEALLOCATE ( var % scalex )
    DEALLOCATE ( var % transd )
    DEALLOCATE ( var % trad )
!===sxy
    DEALLOCATE ( var % reffdf)
    DEALLOCATE ( var % reffbm)
    DEALLOCATE ( var % rhocbm)
    DEALLOCATE ( var % extkbm)
    DEALLOCATE ( var % extkdm)
    DEALLOCATE ( var % fbeam )
    DEALLOCATE ( var % cexpkbm)
    DEALLOCATE ( var % cexpkdm)
    DEALLOCATE ( var % transb )
!===sxy   
  END SUBROUTINE dealloc_radiation_type
   
  SUBROUTINE dealloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % hruff_grmx )
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
    DEALLOCATE ( var % za_uv )
    DEALLOCATE ( var % za_tq )
    DEALLOCATE ( var % z0m )
    DEALLOCATE ( var % zref_uv )
    DEALLOCATE ( var % zref_tq )
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
    DEALLOCATE ( var % precip_s )
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

!les  
!!  SUBROUTINE dealloc_grid_type(var, mg)
!!    TYPE(grid_type), INTENT(inout) :: var
!!    INTEGER, INTENT(in) :: mg
!!    DEALLOCATE ( var % frac_g )
!!    DEALLOCATE ( var % csoil_g )
!!    DEALLOCATE ( var % fe_g )
!!    DEALLOCATE ( var % fh_g )
!!  END SUBROUTINE dealloc_grid_type


END MODULE define_types
