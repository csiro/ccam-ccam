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
  INTEGER, PARAMETER	:: mp = 3864  ! # land points
  INTEGER, PARAMETER	:: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER	:: nrb = 3 ! # radiation bands
      include 'newmpar.h'
!  INTEGER, PARAMETER	:: ms = 6  ! # soil layers
!  INTEGER, PARAMETER	:: kl = 18
!  INTEGER, PARAMETER	:: ncp = 3 ! # vegetation carbon stores
!  INTEGER, PARAMETER	:: ncs = 2 ! # soil carbon stores
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
     REAL(r_1), DIMENSION(mp) :: drybal ! energy balance for dry canopy
     REAL(r_1), DIMENSION(mp) :: ebal   ! energy balance per time step (W/m^2)
     REAL(r_1), DIMENSION(mp) :: ebal_tot ! cumulative energy balance (W/m^2)
     REAL(r_1), DIMENSION(mp) :: evap_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_1), DIMENSION(mp) :: osnowd0  ! snow depth, first time step
     REAL(r_1), DIMENSION(mp) :: precip_tot ! cumulative precipitation (mm/dels)
     REAL(r_1), DIMENSION(mp) :: rnoff_tot  ! cumulative runoff (mm/dels)
     REAL(r_1), DIMENSION(mp) :: wbal   ! water balance per time step (mm/dels)
     REAL(r_1), DIMENSION(mp) :: wbal_tot ! cumulative water balance (mm/dels)
     REAL(r_1), DIMENSION(mp) :: wbtot0 ! total soil water (mm), first time step
     REAL(r_1), DIMENSION(mp) :: wetbal ! energy balance for wet canopy
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type 
     REAL(r_1), DIMENSION(mp) :: albsoil ! soil reflectance
     REAL(r_1), DIMENSION(mp) :: bch  ! parameter b in Campbell equation
     REAL(r_1), DIMENSION(mp) :: c3   ! c3 drainage coeff (fraction)
     REAL(r_1), DIMENSION(mp) :: clay ! fraction of soil which is clay
     REAL(r_1), DIMENSION(mp) :: cnsd ! specific heat capapcity of soil (unit??)
     REAL(r_1), DIMENSION(mp) :: css  ! soil heat capacity [kJ/kg/C]
     REAL(r_1), DIMENSION(mp,ms) :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(mp) :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL(r_1), DIMENSION(mp) :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER(i_d), DIMENSION(mp) :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
     INTEGER(i_d), DIMENSION(mp) :: ibp2   !parameter two in K vis suction (function of pbch)
     INTEGER(i_d), DIMENSION(mp) :: isoilm ! soil type
     REAL(r_1), DIMENSION(mp) :: rhosoil ! soil density [kg/m3]
     REAL(r_1), DIMENSION(mp) :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL(r_1), DIMENSION(mp) :: sand  ! fraction of soil which is sand
     REAL(r_1), DIMENSION(mp) :: sfc   ! vol H2O @ field capacity
     REAL(r_1), DIMENSION(mp) :: silt  ! fraction of soil which is silt
     REAL(r_1), DIMENSION(mp) :: ssat  ! vol H2O @ saturation
     REAL(r_1), DIMENSION(mp) :: sucs  ! suction at saturation (m)
     REAL(r_1), DIMENSION(mp) :: swilt ! vol H2O @ wilting
     REAL(r_1), DIMENSION(ms) :: zse   ! thickness of each soil layer (1=top) in m
     REAL(r_1), DIMENSION(ms+1) :: zshh ! depth from surface of mid-point of a layer (m)
  END TYPE soil_parameter_type
  ! Soil and snow variables:
  TYPE soil_snow_type 
     REAL(r_1), DIMENSION(mp,nrb) :: albsoilsn ! soil + snow reflectance
     REAL(r_1), DIMENSION(mp) :: cls     ! factor for latent heat
     REAL(r_1), DIMENSION(mp) :: dfn_dtg ! d(canopy%fns)/d(ssoil%tgg)
     REAL(r_1), DIMENSION(mp) :: dfh_dtg ! d(canopy%fhs)/d(ssoil%tgg)
     REAL(r_1), DIMENSION(mp) :: dfe_ddq ! d(canopy%fes)/d(dq)
     REAL(r_1), DIMENSION(mp) :: ddq_dtg ! d(dq)/d(ssoil%tgg)
     REAL(r_2), DIMENSION(mp,ms) :: gammzz ! heat capacity for each soil layer
     INTEGER(i_d), DIMENSION(mp) :: isflag ! 0 => no snow 1 => snow
     REAL(r_1), DIMENSION(mp) :: osnowd  ! snow depth from previous time step
     REAL(r_1), DIMENSION(mp) :: potev   ! potential evapotranspiration
     REAL(r_1), DIMENSION(mp) :: runoff  ! total runoff (mm/dels)
     REAL(r_1), DIMENSION(mp) :: rnof1   ! surface runoff (mm/dels)
     REAL(r_1), DIMENSION(mp) :: rnof2   ! deep drainage (mm/dels)
     REAL(r_1), DIMENSION(mp) :: rtsoil  ! turbulent resistance for soil
     REAL(r_1), DIMENSION(mp,3) :: sdepth ! 
     REAL(r_1), DIMENSION(mp,3)  :: smass  ! snow mass
     REAL(r_1), DIMENSION(mp) :: snage   ! snow age
     REAL(r_1), DIMENSION(mp) :: snowd   ! snow depth (liquid water)
     REAL(r_1), DIMENSION(mp,3)  :: ssdn ! snow densities
     REAL(r_1), DIMENSION(mp) :: ssdnn   ! average snow density
     REAL(r_1), DIMENSION(mp,ms) :: tgg  ! soil temperature in K
     REAL(r_1), DIMENSION(mp,3)  :: tggsn  ! snow temperature in K
     REAL(r_1), DIMENSION(mp) :: tss     ! surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(mp,ms) :: wb   ! volumetric soil moisture (solid+liq)
     REAL(r_1), DIMENSION(mp,ms) :: wbfice !
     REAL(r_2), DIMENSION(mp,ms) :: wbice  ! soil ice
     REAL(r_2), DIMENSION(mp,ms) :: wblf !
     REAL(r_1), DIMENSION(mp) :: wbtot   ! total soil water (mm)
  END TYPE soil_snow_type
  ! Vegetation parameters:
  TYPE veg_parameter_type
     INTEGER(i_d),DIMENSION(mp) :: iveg ! vegetation type
     INTEGER(i_d),DIMENSION(mp) :: meth ! method for calculation of canopy fluxes and temp.
     REAL(r_1), DIMENSION(mp) :: vlai   ! leaf area index
     REAL(r_1), DIMENSION(mp) :: vlaimax ! ???
     REAL(r_1), DIMENSION(mp) :: vlaiw  ! lai adjusted for snow depth for calculation of reistances
     REAL(r_1), DIMENSION(mp) :: fwet   ! fraction of canopy wet
     REAL(r_1), DIMENSION(mp) :: canst1 ! max intercepted water by canopy (mm/LAI)
     REAL(r_1), DIMENSION(mp) :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
     REAL(r_1), DIMENSION(mp) :: frac4  ! fraction of c4 plants
     REAL(r_1), DIMENSION(mp) :: tminvj ! min temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(mp) :: tmaxvj ! max temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(mp) :: vbeta  ! ???
     REAL(r_1), DIMENSION(mp) :: hc	! roughness height of canopy (veg - snow)
     REAL(r_1), DIMENSION(mp) :: shelrb ! sheltering factor (dimensionless)
     REAL(r_1), DIMENSION(mp) :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
     REAL(r_1), DIMENSION(mp) :: xfang  ! leaf angle PARAMETER
     REAL(r_1), DIMENSION(mp) :: dleaf ! chararacteristc legnth of leaf (m)
     REAL(r_1), DIMENSION(mp) :: rp20   ! plant respiration coefficient at 20 C
     REAL(r_1), DIMENSION(mp) :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
  END TYPE veg_parameter_type
  ! Canopy/vegetation variables:
  TYPE canopy_type
     REAL(r_1), DIMENSION(mp) :: cansto ! canopy water storage (mm)
     REAL(r_1), DIMENSION(mp) :: delwc ! change in canopy water store (mm/dels)
     REAL(r_1), DIMENSION(mp) :: dewmm ! dewfall (mm)
     REAL(r_1), DIMENSION(mp) :: fe    ! total latent heat (W/m2)
     REAL(r_1), DIMENSION(mp) :: fh    ! total sensible heat (W/m2)
     REAL(r_1), DIMENSION(mp) :: fpn   ! plant photosynthesis (g C s-1)
     REAL(r_1), DIMENSION(mp) :: frp   ! plant respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(mp) :: frpw  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(mp) :: frpr  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(mp) :: frs   ! soil respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(mp) :: fnee  ! net carbon flux (g C m-2 s-1)
     REAL(r_1), DIMENSION(mp) :: frday ! daytime leaf resp
     REAL(r_1), DIMENSION(mp) :: fnv   ! net rad. avail. to canopy (W/m2)
     REAL(r_1), DIMENSION(mp) :: fev   ! latent hf from canopy (W/m2)
     REAL(r_2), DIMENSION(mp) :: fevc  ! dry canopy transpiration (W/m2)
     REAL(r_1), DIMENSION(mp) :: fevw  ! lat heat fl wet canopy (W/m2)
     REAL(r_1), DIMENSION(mp) :: fhv   ! sens heatfl from canopy (W/m2)
     REAL(r_1), DIMENSION(mp) :: fhvw  ! sens heatfl from wet canopy (W/m2)
     REAL(r_1), DIMENSION(mp) :: fns   ! net rad avail to soil (W/m2)
     REAL(r_1), DIMENSION(mp) :: fes   ! latent heatfl from soil (W/m2)
     REAL(r_1), DIMENSION(mp) :: fhs   !sensible heat fl from soil
     REAL(r_1), DIMENSION(mp) :: tv    ! vegetation temp (K)
     REAL(r_1), DIMENSION(mp) :: ga    ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(mp) :: ghflux  ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(mp) :: sghflux ! ground heat flux (W/m2) ???
     REAL(r_2), DIMENSION(mp) :: dgdtg ! derivative of gflux wrt soil temp
     REAL(r_1), DIMENSION(mp) :: through ! canopy throughfall (mm)
     REAL(r_1), DIMENSION(mp) :: precis! throughfall to soil, after snow (mm)
     REAL(r_1), DIMENSION(mp) :: rnet  ! ???
     REAL(r_1), DIMENSION(mp) :: spill ! can.storage excess after dewfall (mm)
     REAL(r_1), DIMENSION(mp) :: wcint ! canopy rainfall interception (mm)
     REAL(r_1), DIMENSION(mp) :: us    ! friction velocity
     REAL(r_1), DIMENSION(mp) :: tscrn ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(mp) :: qscrn ! specific humudity at screen height (g/g)
     REAL(r_1), DIMENSION(mp) :: uscrn ! wind speed at screen height (m/s)
     REAL(r_1), DIMENSION(mp) :: cduv  ! drag coefficient for momentum
  END TYPE canopy_type
  ! Radiation variables:
  TYPE radiation_type
     REAL(r_1), DIMENSION(mp,nrb) :: albedo ! canopy+soil albedo
     REAL(r_1), DIMENSION(mp)     :: extkb  ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(mp)     :: extkd2 ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(mp) ::  extkd ! diffuse radiation extinction coeff (-)
     REAL(r_1), DIMENSION(mp) ::  extkn ! vertical leaf nitrogen profile (-)
     REAL(r_1), DIMENSION(mp) :: flws   ! soil long-wave radiation
     REAL(r_1), DIMENSION(mp,mf)  :: fvlai  ! leaf area index of big leaf
     REAL(r_1), DIMENSION(mp,mf)  :: gradis ! radiative conductance
     REAL(r_1), DIMENSION(mp) :: latitude  ! latitude
     REAL(r_1), DIMENSION(mp) :: lwabv ! long wave absorbed by vegetation
     REAL(r_1), DIMENSION(mp,mf,nrb) :: qcan ! absorbed radiation for canopy (W/m^2)
     REAL(r_1), DIMENSION(mp)     :: qssabs ! absorbed short-wave radiation for soil
     REAL(r_1), DIMENSION(mp,nrb) ::  rhocdf ! canopy diffuse reflectance (-)
     REAL(r_1), DIMENSION(mp,mf)  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
     REAL(r_1), DIMENSION(mp,mf)  :: scalex ! scaling PARAMETER for big leaf
     REAL(r_1), DIMENSION(mp)     :: transd ! fraction SW diffuse transmitted through canopy
     REAL(r_1), DIMENSION(mp) ::  trad  ! ???
  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL(r_1), DIMENSION(mp)	:: coexp
     REAL(r_1), DIMENSION(mp)	:: disp     ! zero-plane displacement
     REAL(r_1), DIMENSION(mp)	:: rt0us
     REAL(r_1), DIMENSION(mp)	:: rt1usa
     REAL(r_1), DIMENSION(mp)	:: rt1usb
     REAL(r_1), DIMENSION(mp)	:: rt1 ! 1/aerodynamic conductance
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL(r_1), DIMENSION(mp)	:: usuh 
     REAL(r_1), DIMENSION(mp)	:: za
     REAL(r_1), DIMENSION(mp)	:: z0m      ! roughness length
     REAL(r_1), DIMENSION(mp)	:: zref
     REAL(r_1), DIMENSION(mp)	:: zruffs
     REAL(r_1), DIMENSION(mp)	:: z0soilsn
  END TYPE roughness_type
  ! Air variables:
  TYPE air_type
     REAL(r_1), DIMENSION(mp) :: rho  ! dry air density (kg m-3)
     REAL(r_1), DIMENSION(mp) :: volm ! molar volume (m3 mol-1)
     REAL(r_1), DIMENSION(mp) :: rlam ! latent heat for water (j/kg)
     REAL(r_1), DIMENSION(mp) :: qsat ! saturation specific humidity
     REAL(r_1), DIMENSION(mp) :: epsi ! d(qsat)/dT ((kg/kg)/K)
     REAL(r_1), DIMENSION(mp) :: visc ! air kinematic viscosity (m2/s)
     REAL(r_1), DIMENSION(mp) :: psyc ! psychrometric constant
     REAL(r_1), DIMENSION(mp) :: dsatdk ! d(es)/dT (mb/K)
     REAL(r_1), DIMENSION(mp) :: cmolar ! conv. from m/s to mol/m2/s
  END TYPE air_type
  ! Meterological data:
  TYPE met_type
     REAL(r_1), DIMENSION(mp) :: ca   ! CO2 concentration
     INTEGER(i_d), DIMENSION(mp) :: year ! local time year AD 
     INTEGER(i_d), DIMENSION(mp) :: moy  ! local time month of year 
     REAL(r_1), DIMENSION(mp) :: doy  ! local time day of year = days since 0 hr 1st Jan 
     REAL(r_1), DIMENSION(mp) :: hod  ! local hour of day
     REAL(r_1), DIMENSION(mp) :: fsd  ! downward short-wave radiation (W/m2)
     REAL(r_1), DIMENSION(mp) :: fld  ! downward long-wave radiation (W/m2)
     REAL(r_1), DIMENSION(mp) :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL(r_1), DIMENSION(mp) :: precips ! solid only (mm/dels)
     REAL(r_1), DIMENSION(mp) :: tc	 ! surface air temperature (oC)
     REAL(r_1), DIMENSION(mp) :: tk	 ! surface air temperature (oK)
     REAL(r_1), DIMENSION(mp) :: tvair   ! within canopy air temperature (oK)
     REAL(r_1), DIMENSION(mp) :: tvrad   ! ???
     REAL(r_1), DIMENSION(mp) :: pmb     ! surface air pressure (mbar)
     REAL(r_1), DIMENSION(mp) :: ua	 ! surface wind speed (m/s)
     REAL(r_1), DIMENSION(mp) :: qv	 ! surface specific humidity (g/g)
     REAL(r_1), DIMENSION(mp) :: qvair   ! within canopy specific humidity (g/g)
     REAL(r_1), DIMENSION(mp) :: da	 ! surf water vap pressure deficit (Pa)
     REAL(r_1), DIMENSION(mp) :: dva	 ! surf water vap pressure deficit (Pa) ???
     REAL(r_1), DIMENSION(mp) :: coszen  ! cos(zenith angle of sun)
  END TYPE met_type
  ! Cumulative flux variables:
  TYPE sum_flux_type
     REAL(r_1), DIMENSION(mp)	:: sumpn  ! sum of canopy photosynthesis (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: sumrp  ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: sumrpw ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: sumrpr ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: sumrs  ! sum of soil respiration (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: sumrd  ! sum of daytime respiration (g C m-2)
     REAL(r_1), DIMENSION(mp)	:: dsumpn ! daily sumpn
     REAL(r_1), DIMENSION(mp)	:: dsumrp ! daily sumrp
     REAL(r_1), DIMENSION(mp)	:: dsumrs ! daily sumrs
     REAL(r_1), DIMENSION(mp)	:: dsumrd ! daily sumrd
     REAL(r_1), DIMENSION(mp)	:: sumxrp ! sum plant resp. modifier
     REAL(r_1), DIMENSION(mp)	:: sumxrs ! sum soil resp. modifier
  END TYPE sum_flux_type
  TYPE bgc_pool_type
     REAL(r_1), DIMENSION(mp,ncp) :: cplant ! plant carbon (g C/m2))
     REAL(r_1), DIMENSION(mp,ncs) :: csoil  ! soil carbon (g C/m2)
     REAL(r_1), DIMENSION(ncp)	:: ratecp ! plant carbon rate constant (1/year)
     REAL(r_1), DIMENSION(ncs)	:: ratecs ! soil carbon rate constant (1/year)
  END TYPE bgc_pool_type
END MODULE define_types
!=========================================================================
MODULE math_constants
  USE define_dimensions
  REAL(r_1), PARAMETER	:: pi = 3.1415927
  REAL(r_1), PARAMETER	:: pi180 = pi / 180.0 ! radians / degree
  REAL(r_1), PARAMETER	:: two_pi = 2.0 * pi
END MODULE math_constants
!=========================================================================
MODULE physical_constants
  USE define_dimensions
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
  USE define_dimensions
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
  USE define_dimensions
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
