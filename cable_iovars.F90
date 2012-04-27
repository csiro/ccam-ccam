! cable_iovars.f90
!
! Contains input/output related variables for CABLE; 
!
! Gab Abramowitz 2009 University of New South Wales, gabsun@gmail.com
!
! Contains one module without subroutines.
!
! Please send bug reports to Bernard.Pak@csiro.au
!
MODULE io_variables
  USE define_dimensions, ONLY: r_1, r_2, i_d, mvtype, mstype
  IMPLICIT NONE
  PUBLIC
  PRIVATE r_1, r_2, i_d, mvtype, mstype
  
  ! ============ Timing variables =====================
  REAL(r_1)    :: shod ! start time hour-of-day
  INTEGER(i_d) :: sdoy,smoy,syear ! start time day-of-year month and year
  CHARACTER(LEN=33) :: timeunits ! timing info read from nc file
  CHARACTER(LEN=3) :: time_coord ! GMT or LOCal time variables
  REAL(r_2),POINTER,DIMENSION(:) :: timevar ! time variable from file
  INTEGER(i_d),DIMENSION(12) :: daysm = &
       (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER(i_d),DIMENSION(12) :: daysml = &
       (/31,29,31,30,31,30,31,31,30,31,30,31/)
  INTEGER(i_d),DIMENSION(12) :: lastday = &
       (/31,59,90,120,151,181,212,243,273,304,334,365/)
  INTEGER(i_d),DIMENSION(12) :: lastdayl = &
       (/31,60,91,121,152,182,213,244,274,305,335,366/)
  LOGICAL :: leaps   ! use leap year timing?
  
  ! ============ Structure variables ===================
  REAL(r_1), POINTER,DIMENSION(:) :: latitude, longitude
  REAL(r_1),POINTER, DIMENSION(:,:) :: lat_all, lon_all ! lat and lon
  CHARACTER(LEN=4) :: metGrid ! Either 'land' or 'mask'
  INTEGER(i_d),POINTER,DIMENSION(:,:) :: mask ! land/sea mask from met file
  INTEGER(i_d),POINTER,DIMENSION(:) :: land_x,land_y ! indicies of land in mask
  INTEGER(i_d) :: xdimsize,ydimsize ! sizes of x and y dimensions
  INTEGER(i_d) :: ngridcells ! number of gridcells in simulation
  TYPE patch_type ! For vegetated surface type
     REAL(r_1) :: frac ! fractional cover of each veg patch
     REAL(r_1) :: latitude
     REAL(r_1) :: longitude
  END TYPE patch_type
  TYPE land_type  
     INTEGER(i_d) :: nap ! number of active (>0%) patches (<=max_vegpatches)
     INTEGER(i_d) :: cstart ! pos of 1st gridcell veg patch in main arrays
     INTEGER(i_d) :: cend ! pos of last gridcell veg patch in main arrays
     INTEGER(i_d) :: ilat ! replacing land_y  ! ??
     INTEGER(i_d) :: ilon ! replacing land_x  ! ??
  END TYPE land_type
  TYPE(land_type),DIMENSION(:),POINTER :: landpt 
  TYPE(patch_type), DIMENSION(:), POINTER :: patch
  INTEGER(i_d) :: max_vegpatches ! The maximum # of patches in any grid cell
  INTEGER(i_d) :: nmetpatches ! size of patch dimension in met file, if exists
 
  ! =============== File details ==========================
  TYPE filenames_type
     CHARACTER(LEN=99) :: met         ! name of file for CABLE input
     CHARACTER(LEN=99) :: out         ! name of file for CABLE output
     CHARACTER(LEN=99) :: log         ! name of file for execution log
     CHARACTER(LEN=99) :: restart_in  ! name of restart file to read
     CHARACTER(LEN=99) :: restart_out ! name of restart file to read
     CHARACTER(LEN=99) :: LAI         ! name of file for default LAI
     CHARACTER(LEN=99) :: type        ! file for default veg/soil type
     CHARACTER(LEN=99) :: veg         ! file for vegetation parameters
     CHARACTER(LEN=99) :: soil        ! name of file for soil parameters
     CHARACTER(LEN=99) :: inits       ! name of file for initialisations
     CHARACTER(LEN=99) :: soilIGBP    ! name of file for IGBP soil map
  END TYPE filenames_type
  TYPE(filenames_type) :: filename
  TYPE gswp_type
     CHARACTER(LEN=99) :: rainf
     CHARACTER(LEN=99) :: snowf
     CHARACTER(LEN=99) :: LWdown
     CHARACTER(LEN=99) :: SWdown
     CHARACTER(LEN=99) :: PSurf
     CHARACTER(LEN=99) :: Qair
     CHARACTER(LEN=99) :: Tair
     CHARACTER(LEN=99) :: wind
  END TYPE gswp_type
  TYPE(gswp_type)      :: gswpfile
  INTEGER(i_d) :: ncciy    ! year number (& switch) for gswp run
  INTEGER(i_d) :: ncid_rin ! input netcdf restart file ID
  INTEGER(i_d) :: logn     ! log file unit number
  LOGICAL :: verbose ! print init and param details of all grid cells?
  LOGICAL :: soilparmnew ! read IGBP new soil map. Q.Zhang @ 12/20/2010
  LOGICAL :: redistrb    ! Turn on/off the hydraulic redistribution
  REAL(r_1) :: wiltParam
  REAL(r_1) :: satuParam

  ! ================ Veg and soil type variables ====================================
  INTEGER(i_d),POINTER :: soiltype_metfile(:,:) ! user defined soil type (from met file)
  INTEGER(i_d),POINTER :: vegtype_metfile(:,:) ! user-def veg type (from met file)
! added mvtype and mstype to define_dimensions (BP sep2010)
!  INTEGER(i_d) :: mvtype ! # of vegetation types in veg classification scheme
!  INTEGER(i_d) :: mstype ! # of soil types in soil classification scheme
  CHARACTER(LEN=70), DIMENSION(:), POINTER :: veg_desc ! decriptions of veg type
  CHARACTER(LEN=70), DIMENSION(:), POINTER :: soil_desc ! decriptns of soil type 
  TYPE parID_type ! model parameter IDs in netcdf file
     INTEGER(i_d) :: bch,latitude,clay,css,rhosoil,hyds,rs20,sand,sfc,silt, &
          ssat,sucs,swilt,froot,zse,canst1,dleaf,meth,za_tq,za_uv, &
          ejmax,frac4,hc,lai,rp20,rpcoef,shelrb, vbeta, xalbnir, &
          vcmax,xfang,ratecp,ratecs,refsbare,isoil,iveg,albsoil,&
          taul,refl,tauw,refw,wai,vegcf,extkn,tminvj,tmaxvj, & 
          veg_class,soil_class,mvtype,mstype,patchfrac
  END TYPE parID_type
  
  ! =============== Logical  variables ============================
  TYPE input_details_type
     LOGICAL :: Wind ! T => 'Wind' is present; F => use vector component wind
     LOGICAL :: LWdown ! T=> downward longwave is present in met file
     LOGICAL :: CO2air ! T=> air CO2 concentration is present in met file
     LOGICAL :: PSurf ! T=> surface air pressure is present in met file
     LOGICAL :: Snowf ! T=> snowfall variable is present in met file
     LOGICAL :: avPrecip! T=> ave rainfall present in met file (use for spinup)
     LOGICAL :: LAI   ! T=> LAI is present in the met file
     LOGICAL :: LAI_T ! T=> LAI is time dependent, for each time step
     LOGICAL :: LAI_M ! T=> LAI is time dependent, for each month
     LOGICAL :: LAI_P ! T=> LAI is patch dependent
     LOGICAL :: parameters ! TRUE if non-default parameters are found
     LOGICAL :: initial ! switched to TRUE when initialisation data are loaded
     LOGICAL :: patch ! T=> met file have a subgrid veg/soil patch dimension
     LOGICAL :: laiPatch ! T=> LAI file have a subgrid veg patch dimension
  END TYPE input_details_type
  TYPE(input_details_type) :: exists
  TYPE output_inclusion_type
    ! Which variables to include in output file, values initialised here
    ! and can be reset by namelist file read in driver:
    ! Groups of output variables:
    LOGICAL :: met = .FALSE. ! input met data
    LOGICAL :: flux = .FALSE.  ! convective, runoff, NEE
    LOGICAL :: radiation = .FALSE. ! net rad, albedo
    LOGICAL :: carbon = .FALSE. ! NEE, GPP, NPP, stores 
    LOGICAL :: soil = .FALSE.  ! soil states
    LOGICAL :: snow = .FALSE.  ! snow states
    LOGICAL :: veg = .FALSE. ! vegetation states
    LOGICAL :: params = .FALSE. ! input parameters used to produce run
    LOGICAL :: balances = .FALSE. ! energy and water balances
    LOGICAL :: restart = .FALSE. ! create restart file?
    LOGICAL :: ensemble = .FALSE. ! are we creating an ensemble run?
    LOGICAL :: patch = .FALSE.  ! should patch-specific info be written to output file?
    ! Should output grid follow met file 'default'; force with 'land' or 'mask':
    CHARACTER(LEN=7) :: grid = 'default' 
    CHARACTER(LEN=7) :: averaging = 'all' ! 'all', 'daily', 'monthly', 'user6'(6hrly)
    INTEGER(i_d) :: interval ! in case of 'user6' above, interval will be 6
    ! variables specified individually:
    LOGICAL :: SWdown = .FALSE.  ! 6 downward short-wave radiation [W/m2]
    LOGICAL :: LWdown = .FALSE.   ! 7 downward long-wave radiation [W/m2]
    LOGICAL :: Rainf = .FALSE.    ! 8 rainfall [kg/m2/s]
    LOGICAL :: Snowf = .FALSE.    ! 9 snowfall [kg/m2/s]
    LOGICAL :: PSurf = .FALSE.    ! 10 surface pressure [Pa]
    LOGICAL :: Tair = .FALSE.   ! 11 surface air temperature [K]
    LOGICAL :: Qair = .FALSE.  ! 12 specific humidity [kg/kg]
    LOGICAL :: CO2air = .FALSE. ! 13 CO2 concentration [ppmv]
    LOGICAL :: Wind = .FALSE.   ! 14 windspeed [m/s]
    LOGICAL :: Wind_N = .FALSE. ! 15 surface wind speed, N component [m/s]
    LOGICAL :: Wind_E = .FALSE. ! 16 surface wind speed, E component [m/s]
    LOGICAL :: LAI = .FALSE.      !
    LOGICAL :: Qh = .FALSE.      ! 17 sensible heat flux [W/m2]
    LOGICAL :: Qle = .FALSE.      ! 18 latent heat flux [W/m2]
    LOGICAL :: Qg = .FALSE.       ! 19 ground heat flux [W/m2]
    LOGICAL :: SWnet = .FALSE.    ! 20 net shortwave [W/m2]
    LOGICAL :: LWnet = .FALSE.    ! 21 net longwave [W/m2]
    LOGICAL :: Evap = .FALSE.     ! 22 total evapotranspiration [kg/m2/s]
    LOGICAL :: Ewater = .FALSE. ! 23 evap. from surface water storage [kg/m2/s]
    LOGICAL :: ESoil = .FALSE.    ! 24 bare soil evaporation [kg/m2/s]
    LOGICAL :: TVeg = .FALSE.     ! 25 vegetation transpiration [kg/m2/s]
    LOGICAL :: ECanop = .FALSE.   ! 26 interception evaporation [kg/m2/s]
    LOGICAL :: PotEvap = .FALSE.  ! 27 potential evapotranspiration [kg/m2/s]
    LOGICAL :: ACond = .FALSE.    ! 28 aerodynamic conductance [m/s]
    LOGICAL :: SoilWet = .FALSE.  ! 29 total soil wetness [-] 
    LOGICAL :: Albedo = .FALSE.   ! 30 albedo [-] 
    LOGICAL :: VegT = .FALSE.    ! 31 vegetation temperature [K]
    LOGICAL :: SoilTemp = .FALSE.  ! 32 av.layer soil temperature [K]
    LOGICAL :: SoilMoist = .FALSE. ! 33 av.layer soil moisture [kg/m2]
    LOGICAL :: Qs = .FALSE.           ! 34 surface runoff [kg/m2/s]
    LOGICAL :: Qsb = .FALSE.          ! 35 subsurface runoff [kg/m2/s]
    LOGICAL :: DelSoilMoist = .FALSE. ! 36 change in soilmoisture (sum layers) [kg/m2]
    LOGICAL :: DelSWE = .FALSE.       ! 37 change in snow water equivalent [kg/m2]
    LOGICAL :: DelIntercept = .FALSE. ! 38 change in interception storage [kg/m2]
    LOGICAL :: SnowT = .FALSE.        ! 39 snow surface temp [K]
    LOGICAL :: BaresoilT = .FALSE.    ! 40 surface bare soil temp [K]
    LOGICAL :: AvgSurfT = .FALSE.     ! 41 Average surface temperature [K]
    LOGICAL :: RadT = .FALSE.         ! 42 Radiative surface temperature [K]
    LOGICAL :: SWE = .FALSE.          ! 43 snow water equivalent [kg/m2]
    LOGICAL :: RootMoist = .FALSE.    ! 44 root zone soil moisture [kg/m2]
    LOGICAL :: CanopInt = .FALSE.     ! 45 total canopy water storage [kg/m2]
    LOGICAL :: NEE  = .FALSE.         ! 46 net ecosystem exchange [umol/m2/s]
    LOGICAL :: NPP  = .FALSE.         ! 47 net primary production of C by veg [umol/m2/s]
    LOGICAL :: GPP = .FALSE.          ! 48 gross primary production C by veg [umol/m2/s]
    LOGICAL :: AutoResp = .FALSE.     ! 49 autotrophic respiration [umol/m2/s]
    LOGICAL :: LeafResp = .FALSE.     ! 51 autotrophic respiration [umol/m2/s]
    LOGICAL :: HeteroResp = .FALSE.   ! 50 heterotrophic respiration [umol/m2/s]
    LOGICAL :: SnowDepth = .FALSE.    ! actual depth of snow in [m]
    ! Non-Alma variables
    LOGICAL :: Rnet = .FALSE.         ! net absorbed radiation [W/m2]
    LOGICAL :: HVeg = .FALSE.         ! sensible heat from vegetation [W/m2]
    LOGICAL :: HSoil = .FALSE.        ! sensible heat from soil [W/m2]
    LOGICAL :: Ebal = .FALSE.         ! cumulative energy balance [W/m2]
    LOGICAL :: Wbal = .FALSE.         ! cumulative water balance [W/m2]
    ! Model parameters
    LOGICAL :: bch = .FALSE.       ! parameter b in Campbell equation 1985
    LOGICAL :: latitude = .FALSE.  ! site latitude
    LOGICAL :: clay = .FALSE.      ! fraction of clay in soil
    LOGICAL :: css = .FALSE.       ! heat capacity of soil minerals [J/kg/C]
    LOGICAL :: rhosoil = .FALSE.   ! soil density [kg/m3]
    LOGICAL :: hyds = .FALSE.      ! hydraulic conductivity @ saturation [m/s], Ksat
    LOGICAL :: rs20 = .FALSE.      ! soil respiration at 20 C [dimensionless], (0.1 - 10), prop to om
    LOGICAL :: sand  = .FALSE.     ! fraction of sand in soil
    LOGICAL :: sfc = .FALSE.       ! vol H2O @ field capacity
    LOGICAL :: silt  = .FALSE.     ! fraction of silt in soil
    LOGICAL :: ssat = .FALSE.      ! vol H2O @ saturation
    LOGICAL :: sucs = .FALSE.      ! suction at saturation [m]
    LOGICAL :: swilt = .FALSE.     ! vol H2O @ wilting
    LOGICAL :: froot = .FALSE.   ! fraction of roots in each soil layer
    LOGICAL :: zse = .FALSE.     ! thickness of each soil layer (1=top) (m)
    LOGICAL :: canst1 = .FALSE.    ! max intercepted water by canopy [mm/LAI] (0.08 - 0.12) {avoid}
    LOGICAL :: dleaf = .FALSE.     ! chararacteristic length of leaf [m], (0.005 - 0.2) pine -> tropical
    LOGICAL :: ejmax  = .FALSE.    ! max pot. electron transport rate top leaf[mol/m2/s](1e-5 - 3e-4) {use}
    LOGICAL :: frac4  = .FALSE.    ! fraction of c4 plants [-]
    LOGICAL :: hc = .FALSE.        ! height of canopy [m]
    LOGICAL :: rp20  = .FALSE.     ! plant respiration coefficient at 20 C [-] 0.1 - 10 (frp 0 - 15e-6 mol/m2/s)
    LOGICAL :: rpcoef  = .FALSE.   ! temperature coef nonleaf plant respiration [1/C] (0.8 - 1.5)
    LOGICAL :: shelrb  = .FALSE.   ! sheltering factor [-] {avoid - insensitive?}
    LOGICAL :: vcmax  = .FALSE.    ! maximum RuBP carboxylation rate top leaf [mol/m2/s](5e-6 - 1.5e-4){use}
    LOGICAL :: xfang  = .FALSE.    ! leaf angle PARAMETER (dimensionless) (v leaf -1.0 horiz 1.0 sphere 0 (-1 - 1))
    LOGICAL :: wai    = .FALSE.    ! wood area index
    LOGICAL :: vegcf  = .FALSE.    ! 
    LOGICAL :: extkn  = .FALSE.    ! 
!    LOGICAL :: rootbeta = .FALSE.  !
    LOGICAL :: ratecp = .FALSE.  ! plant carbon pool rate constant (1/year)
    LOGICAL :: ratecs = .FALSE.  ! soil carbon pool rate constant (1/year)
    LOGICAL :: albsoil = .FALSE.  ! soil reflectance [-]
    LOGICAL :: taul = .FALSE.    ! leaf transmissivity [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
    LOGICAL :: refl = .FALSE.    ! leaf reflectance [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
    LOGICAL :: tminvj = .FALSE.    ! min temperature of the start of photosynthesis(leaf phenology)[-] (-10 - 10)
    LOGICAL :: tmaxvj  = .FALSE.   ! max temperature of the start of photosynthesis(leaf phenology)[-] (-5 - 15)
    LOGICAL :: vbeta = .FALSE.     ! stomatal sensitivity to soil water
    LOGICAL :: xalbnir = .FALSE.   ! modifier for albedo in near ir band
    LOGICAL :: iveg  = .FALSE.     ! vegetation type from global index
    LOGICAL :: patchfrac  = .FALSE.! fractional cover of each veg/soil patch
    LOGICAL :: isoil  = .FALSE.    ! soil type from global index
    LOGICAL :: meth  = .FALSE.     ! method for solving turbulence in canopy scheme
    LOGICAL :: za  = .FALSE.       ! something to do with roughness ????
  END TYPE output_inclusion_type
  TYPE(output_inclusion_type),SAVE :: output ! do these vars get written to output?
  TYPE(output_inclusion_type),SAVE :: patchout ! do we want patch-specific info
  TYPE checks_type
     LOGICAL :: ranges, energy_bal, mass_bal
  END TYPE checks_type
  TYPE(checks_type) :: check ! what types of checks to perform
  
  ! ============== Proxy input variables ================================
  REAL(r_1),POINTER,DIMENSION(:)  :: PrecipScale! precip scaling per site for spinup
  REAL(r_1),POINTER,DIMENSION(:,:)  :: defaultLAI ! in case met file/host model has no LAI
  REAL(r_1) :: fixedCO2 ! CO2 level if CO2air not in met file

  ! For threading:
  !$OMP THREADPRIVATE(landpt,patch)

END MODULE io_variables
