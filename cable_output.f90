!!$ cable_output.f90
!!$
!!$ Output module for CABLE land surface scheme offline netcdf driver; 
!!$
!!$ Gab Abramowitz 2007 University of New South Wales/
!!$ CSIRO Marine and Atmospheric Research; gabsun@gmail.com
!!$
!!$ Subroutines in this file open the output netcdf file for offline CABLE
!!$ and write to it, providing details about variable units and names. 

MODULE output_module
  USE input_module
  IMPLICIT NONE
  INTEGER :: ncid_out ! output data netcdf file ID
  INTEGER :: tvarID ! time variable ID
  TYPE out_varID_type ! output variable IDs in netcdf file
     INTEGER :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair,Qair,Rainf, &
          Snowf,CO2air,Qle,Qh,Qg,NEE,SWnet,LWnet,SoilMoist, &
          SoilTemp,Albedo,Qs,Qsb,Evap,BaresoilT,SWE,SnowT,RadT, &
          VegT,Ebal,Wbal,AutoResp,HeteroResp,GPP,NPP,LAI,ECanop, &
          TVeg,ESoil,CanopInt,SnowDepth,HVeg,HSoil,Rnet
  END TYPE out_varID_type
  TYPE(out_varID_type) :: ovid ! netcdf variable IDs for output variables
  TYPE(parID_type) :: opid ! netcdf variable IDs for output variables
  TYPE dataset_type
     ! Possible forcing variables for CABLE:
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SWdown   ! 6 downward short-wave radiation [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: LWdown   ! 7 downward long-wave radiation [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Rainf    ! 8 rainfall [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Snowf    ! 9 snowfall [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: PSurf    ! 10 surface pressure [Pa]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: Tair   ! 11 surface air temperature [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: Qair   ! 12 specific humidity [kg/kg]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: CO2air ! 13 CO2 concentration [ppmv]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: Wind   ! 14 windspeed [m/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: Wind_N ! 15 surface wind speed, N component [m/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: Wind_E ! 16 surface wind speed, E component [m/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: LAI      !
     ! Possible output variables:
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Qh       ! 17 sensible heat flux [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Qle      ! 18 latent heat flux [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Qg       ! 19 ground heat flux [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SWnet    ! 20 net shortwave [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: LWnet    ! 21 net longwave [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Evap     ! 22 total evapotranspiration [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Ewater ! 23 evap. from surface water storage [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: ESoil    ! 24 bare soil evaporation [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: TVeg     ! 25 vegetation transpiration [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: ECanop   ! 26 interception evaporation [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: PotEvap  ! 27 potential evapotranspiration [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: ACond    ! 28 aerodynamic conductance [m/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SoilWet  ! 29 total soil wetness [-] 
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Albedo   ! 30 albedo [-] 
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: VegT     ! 31 vegetation temperature [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: SoilTemp  ! 32 av.layer soil temperature [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:,:) :: SoilMoist ! 33 av.layer soil moisture [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Qs           ! 34 surface runoff [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Qsb          ! 35 subsurface runoff [kg/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: DelSoilMoist ! 36 change in soilmoisture (sum layers) [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: DelSWE       ! 37 change in snow water equivalent [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: DelIntercept ! 38 change in interception storage [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SnowT        ! 39 snow surface temp [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: BaresoilT    ! 40 surface bare soil temp [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: AvgSurfT     ! 41 Average surface temperature [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: RadT         ! 42 Radiative surface temperature [K]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SWE          ! 43 snow water equivalent [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: RootMoist    ! 44 root zone soil moisture [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: CanopInt     ! 45 total canopy water storage [kg/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: NEE          ! 46 net ecosystem exchange [umol/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: NPP          ! 47 net primary production of C by veg [umol/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: GPP          ! 48 gross primary production C by veg [umol/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: AutoResp     ! 49 autotrophic respiration [umol/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: HeteroResp   ! 50 heterotrophic respiration [umol/m2/s]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: SnowDepth    ! actual depth of snow in [m]
     ! Non-Alma variables
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Rnet         ! net absorbed radiation [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: HVeg         ! sensible heat from vegetation [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: HSoil        ! sensible heat from soil [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Ebal         ! cumulative energy balance [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Wbal         ! cumulative water balance [W/m2]
     ! Model parameters
     REAL(r_1),POINTER,DIMENSION(:,:) :: bch       ! parameter b in Campbell equation 1985
     REAL(r_1),POINTER,DIMENSION(:,:) :: latitude  ! site latitude
     ! Frac of lowest layer soil moisture above field capacity which drains
     REAL(r_1),POINTER,DIMENSION(:,:) :: clay      ! fraction of clay in soil
     REAL(r_1),POINTER,DIMENSION(:,:) :: css       ! heat capacity of soil minerals [J/kg/C]
     REAL(r_1),POINTER,DIMENSION(:,:) :: rhosoil   ! soil density [kg/m3]
     REAL(r_1),POINTER,DIMENSION(:,:) :: hyds      ! hydraulic conductivity @ saturation [m/s], Ksat
     REAL(r_1),POINTER,DIMENSION(:,:) :: rs20      ! soil respiration at 20 C [dimensionless], (0.1 - 10), prop to om
     REAL(r_1),POINTER,DIMENSION(:,:) :: sand      ! fraction of sand in soil
     REAL(r_1),POINTER,DIMENSION(:,:) :: sfc       ! vol H2O @ field capacity
     REAL(r_1),POINTER,DIMENSION(:,:) :: silt      ! fraction of silt in soil
     REAL(r_1),POINTER,DIMENSION(:,:) :: ssat      ! vol H2O @ saturation
     REAL(r_1),POINTER,DIMENSION(:,:) :: sucs      ! suction at saturation [m]
     REAL(r_1),POINTER,DIMENSION(:,:) :: swilt     ! vol H2O @ wilting
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: froot   ! fraction of roots in each soil layer
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: zse     ! thickness of each soil layer (1=top) (m)
     REAL(r_1),POINTER,DIMENSION(:,:) :: canst1    ! max intercepted water by canopy [mm/LAI] (0.08 - 0.12) {avoid}
     REAL(r_1),POINTER,DIMENSION(:,:) :: dleaf     ! chararacteristic length of leaf [m], (0.005 - 0.2) pine -> tropical
     REAL(r_1),POINTER,DIMENSION(:,:) :: ejmax     ! max pot. electron transport rate top leaf[mol/m2/s](1e-5 - 3e-4) {use}
     REAL(r_1),POINTER,DIMENSION(:,:) :: frac4     ! fraction of c4 plants [-]
     REAL(r_1),POINTER,DIMENSION(:,:) :: hc        ! height of canopy [m]
     REAL(r_1),POINTER,DIMENSION(:,:) :: rp20      ! plant respiration coefficient at 20 C [-] 0.1 - 10 (frp 0 - 15e-6 mol/m2/s)
     REAL(r_1),POINTER,DIMENSION(:,:) :: rpcoef    ! temperature coef nonleaf plant respiration [1/C] (0.8 - 1.5)
     REAL(r_1),POINTER,DIMENSION(:,:) :: shelrb    ! sheltering factor [-] {avoid - insensitive?}
     REAL(r_1),POINTER,DIMENSION(:,:) :: vcmax     ! maximum RuBP carboxylation rate top leaf [mol/m2/s](5e-6 - 1.5e-4){use}
     REAL(r_1),POINTER,DIMENSION(:,:) :: xfang     ! leaf angle PARAMETER (dimensionless) (v leaf -1.0 horiz 1.0 sphere 0 (-1 - 1))
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: ratecp  ! plant carbon pool rate constant (1/year)
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: ratecs  ! soil carbon pool rate constant (1/year)
     REAL(r_1),POINTER,DIMENSION(:,:) :: albsoil   ! soil reflectance [-]
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: taul    ! leaf transmissivity [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
     REAL(r_1),POINTER,DIMENSION(:,:,:) :: refl    ! leaf reflectance [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
     REAL(r_1),POINTER,DIMENSION(:,:) :: tminvj    ! min temperature of the start of photosynthesis(leaf phenology)[-] (-10 - 10)
     REAL(r_1),POINTER,DIMENSION(:,:) :: tmaxvj    ! max temperature of the start of photosynthesis(leaf phenology)[-] (-5 - 15)
     REAL(r_1),POINTER,DIMENSION(:,:) :: vbeta     ! stomatal sensitivity to soil water
     REAL(r_1),POINTER,DIMENSION(:,:) :: iveg      ! vegetation type from global index
     REAL(r_1),POINTER,DIMENSION(:,:) :: isoil     ! soil type from global index
     REAL(r_1),POINTER,DIMENSION(:,:) :: meth      ! canopy turbulence parameterisation selection
     REAL(r_1),POINTER,DIMENSION(:,:) :: za        ! reference height - lowest atm model layer
  END TYPE dataset_type
  TYPE (dataset_type) :: out ! CABLE output after units adjustment
  TYPE output_inclusion_type
     ! Which variables to include in output file, initialised here and set in driver:
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
     LOGICAL :: ratecp = .FALSE.  ! plant carbon pool rate constant (1/year)
     LOGICAL :: ratecs = .FALSE.  ! soil carbon pool rate constant (1/year)
     LOGICAL :: albsoil = .FALSE.  ! soil reflectance [-]
     LOGICAL :: taul = .FALSE.    ! leaf transmissivity [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
     LOGICAL :: refl = .FALSE.    ! leaf reflectance [-](V:0.07 - 0.15 NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
     LOGICAL :: tminvj = .FALSE.    ! min temperature of the start of photosynthesis(leaf phenology)[-] (-10 - 10)
     LOGICAL :: tmaxvj  = .FALSE.   ! max temperature of the start of photosynthesis(leaf phenology)[-] (-5 - 15)
     LOGICAL :: vbeta = .FALSE.     ! stomatal sensitivity to soil water
     LOGICAL :: iveg  = .FALSE.     ! vegetation type from global index
     LOGICAL :: isoil  = .FALSE.    ! soil type from global index
     LOGICAL :: meth  = .FALSE.     ! method for solving ??? in canopy scheme
     LOGICAL :: za  = .FALSE.       ! something to do with roughness ????
  END TYPE output_inclusion_type
  TYPE(output_inclusion_type),SAVE :: output

CONTAINS

  SUBROUTINE open_output_file(filename_out,filename_met,dels,soil,veg,bgc,rough)
    CHARACTER(LEN=*), INTENT(IN) :: filename_out 
    CHARACTER(LEN=*), INTENT(IN) :: filename_met 
    REAL(r_1),INTENT(IN) :: dels ! time step size
    TYPE (soil_parameter_type),INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type),INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (bgc_pool_type),INTENT(IN)       :: bgc
    TYPE (roughness_type),INTENT(IN)      :: rough
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
    INTEGER :: xID,yID,zID,radID,soilID,soilcarbID,plantcarbID,tID,landID ! dimension IDs
    INTEGER :: tvarID,latID,lonID ! time,lat,lon variable ID
    INTEGER :: xvID,yvID          ! coordinate variable IDs for GrADS readability
    INTEGER :: i ! do loop counter
    ! Create output file:
    status = NF90_CREATE(filename_out,NF90_CLOBBER,ncid_out)
    IF(status/=NF90_NOERR) CALL nc_abort('Error creating output file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    ! Put the file in define mode:
    status = NF90_REDEF(ncid_out)
    ! Define dimensions:
    status = NF90_DEF_DIM(ncid_out,'x',xdimsize,xID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining x dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'y',ydimsize,yID) 
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining y dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'soil',ms,soilID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining vertical soil dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'rad',nrb,radID) ! number of radiation bands
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining radiation dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'soil_carbon_pools',ncs,soilcarbID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining soil carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'plant_carbon_pools',ncp,plantcarbID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining plant carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_DEF_DIM(ncid_out,'time',NF90_UNLIMITED,tID) 
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    IF(gridType=='mask') THEN ! for land/sea mask type grid:
       ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type:
       status = NF90_DEF_DIM(ncid_out,'z',1,zID) 
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining z dimension in output file. '// &
            '(SUBROUTINE open_output_file)')
    ELSEIF(gridType=='land') THEN ! For land only compression grid:
       status = NF90_DEF_DIM(ncid_out,'land',mp,landID) ! number of land points
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining land dimension in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Define "time" variable and its attributes:
    status=NF90_DEF_VAR(ncid_out,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,tvarID,'units',timeunits)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,tvarID,'coordinate',time_coord)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Define latitude and longitude variable (ALMA):
    status=NF90_DEF_VAR(ncid_out,'latitude',NF90_FLOAT,(/xID,yID/),latID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining latitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,latID,'units','degrees_north')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining latitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    status=NF90_DEF_VAR(ncid_out,'longitude',NF90_FLOAT,(/xID,yID/),lonID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining longitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,lonID,'units','degrees_east')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining longitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Write "cordinate variables" to enable reading by GrADS:
    status=NF90_DEF_VAR(ncid_out,'x',NF90_FLOAT,(/xID/),xvID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining "x" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,xvID,'units','degrees_north')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error writing x coordinate variable (GrADS) units in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,xvID,'comment','x coordinate variable for GrADS compatibility')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error writing x variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
    status=NF90_DEF_VAR(ncid_out,'y',NF90_FLOAT,(/yID/),yvID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining "y" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
     status = NF90_PUT_ATT(ncid_out,yvID,'units','degrees_east')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error writing y coordinate variable (GrADS) units in output file. '// &
         '(SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,yvID,'comment','y coordinate variable for GrADS compatibility')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error writing y variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
   
    !=============DEFINE OUTPUT VARIABLES===================================
    ! SWdown:
    IF(output%met.OR.output%SWdown) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SWdown to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SWdown',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%SWdown) ! Define SWdown
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWdown variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWdown(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SWdown to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SWdown',NF90_FLOAT,(/landID,tID/), &
               ovid%SWdown)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWdown variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWdown(mp,1,1))
       END IF
       ! Define SWdown units:
       status = NF90_PUT_ATT(ncid_out,ovid%SWdown,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWdown,'long_name', &
            'Downward shortwave radiation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! LWdown:
    IF(output%met.OR.output%LWdown) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing LWdown to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'LWdown',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%LWdown) ! Define LWdown
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LWdown variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LWdown(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing LWdown to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'LWdown',NF90_FLOAT,(/landID,tID/), &
               ovid%LWdown)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LWdown variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LWdown(mp,1,1))
       END IF
       ! Define LWdown units:
       status = NF90_PUT_ATT(ncid_out,ovid%LWdown,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWdown,'long_name', &
            'Downward longwave radiation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Tair:
    IF(output%met.OR.output%Tair) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Tair to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Tair',NF90_FLOAT,(/xID,yID,zID,tID/), &
               ovid%Tair) ! Define Tair
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Tair variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Tair(xdimsize,ydimsize,1,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Tair to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Tair',NF90_FLOAT,(/landID,tID/), &
               ovid%Tair)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Tair variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Tair(mp,1,1,1))
       END IF
       ! Define Tair units:
       status = NF90_PUT_ATT(ncid_out,ovid%Tair,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Tair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Tair,'long_name', &
            'Surface air temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Tair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Rainf:
    IF(output%met.OR.output%Rainf) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Rainf to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Rainf',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Rainf) ! Define Rainf
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Rainf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Rainf(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Rainf to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Rainf',NF90_FLOAT,(/landID,tID/), &
               ovid%Rainf)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Rainf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Rainf(mp,1,1))
       END IF
       ! Define Rainf units:
       status = NF90_PUT_ATT(ncid_out,ovid%Rainf,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rainf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Rainf,'long_name', &
            'Rainfall AND snowfall')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rainf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Snowf:
    IF(output%met.OR.output%Snowf) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Snowf to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Snowf',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Snowf) ! Define Snowf
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Snowf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Snowf(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Snowf to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Snowf',NF90_FLOAT,(/landID,tID/), &
               ovid%Snowf)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Snowf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Snowf(mp,1,1))
       END IF
       ! Define Snowf units:
       status = NF90_PUT_ATT(ncid_out,ovid%Snowf,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Snowf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Snowf,'long_name', &
            'Snowfall')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Snowf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qair:
    IF(output%met.OR.output%Qair) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qair to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qair',NF90_FLOAT,(/xID,yID,zID,tID/), &
               ovid%Qair) ! Define Qair
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qair variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qair(xdimsize,ydimsize,1,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qair to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qair',NF90_FLOAT,(/landID,tID/), &
               ovid%Qair)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qair variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qair(mp,1,1,1))
       END IF
       ! Define Qair units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qair,'units','kg/kg')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qair,'long_name', &
            'Surface specific humidity')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Wind:
    IF(output%met.OR.output%Wind) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Wind to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Wind',NF90_FLOAT,(/xID,yID,zID,tID/), &
               ovid%Wind) ! Define Wind
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Wind variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Wind(xdimsize,ydimsize,1,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Wind to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Wind',NF90_FLOAT,(/landID,tID/), &
               ovid%Wind)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Wind variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Wind(mp,1,1,1))
       END IF
       ! Define Wind units:
       status = NF90_PUT_ATT(ncid_out,ovid%Wind,'units','m/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wind variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Wind,'long_name', &
            'Scalar surface wind speed')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wind variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! PSurf:
    IF(output%met.OR.output%PSurf) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing PSurf to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'PSurf',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%PSurf) ! Define PSurf
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining PSurf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%PSurf(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing PSurf to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'PSurf',NF90_FLOAT,(/landID,tID/), &
               ovid%PSurf)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining PSurf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%PSurf(mp,1,1))
       END IF
       ! Define PSurf units:
       status = NF90_PUT_ATT(ncid_out,ovid%PSurf,'units','hPa')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining PSurf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%PSurf,'long_name', &
            'Surface air pressure')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining PSurf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! CO2air:
    IF(output%met.OR.output%CO2air) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing CO2air to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'CO2air',NF90_FLOAT,(/xID,yID,zID,tID/), &
               ovid%CO2air) ! Define CO2air
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining CO2air variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%CO2air(xdimsize,ydimsize,1,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing CO2air to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'CO2air',NF90_FLOAT,(/landID,tID/), &
               ovid%CO2air)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining CO2air variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%CO2air(mp,1,1,1))
       END IF
       ! Define CO2air units:
       status = NF90_PUT_ATT(ncid_out,ovid%CO2air,'units','ppmv')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CO2air variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%CO2air,'long_name', &
            'Surface air CO2 concentration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CO2air variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qle:
    IF(output%flux.OR.output%Qle) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qle to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qle',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Qle) ! Define Qle
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qle variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qle(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qle to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qle',NF90_FLOAT,(/landID,tID/), &
               ovid%Qle)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qle variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qle(mp,1,1))
       END IF
       ! Define Qle units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qle,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qle variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qle,'long_name', &
            'Surface latent heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qle variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qh:
    IF(output%flux.OR.output%Qh) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qh to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qh',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Qh) ! Define Qh
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qh variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qh(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qh to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qh',NF90_FLOAT,(/landID,tID/), &
               ovid%Qh)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qh variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qh(mp,1,1))
       END IF
       ! Define Qh units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qh,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qh variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qh,'long_name', &
            'Surface sensible heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qh variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qg:
    IF(output%flux.OR.output%Qg) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qg to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qg',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Qg) ! Define Qg
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qg(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qg to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qg',NF90_FLOAT,(/landID,tID/), &
               ovid%Qg)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qg(mp,1,1))
       END IF
       ! Define Qg units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qg,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qg,'long_name', &
            'Surface ground heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qs:
    IF(output%flux.OR.output%Qs) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qs to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qs',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Qs) ! Define Qs
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qs(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qs to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qs',NF90_FLOAT,(/landID,tID/), &
               ovid%Qs)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qs(mp,1,1))
       END IF
       ! Define Qs units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qs,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qs,'long_name', &
            'Surface runoff')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Qsb:
    IF(output%flux.OR.output%Qsb) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Qsb to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Qsb',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Qsb) ! Define Qsb
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qsb variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qsb(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Qsb to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Qsb',NF90_FLOAT,(/landID,tID/), &
               ovid%Qsb)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Qsb variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Qsb(mp,1,1))
       END IF
       ! Define Qsb units:
       status = NF90_PUT_ATT(ncid_out,ovid%Qsb,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qsb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qsb,'long_name', &
            'Subsurface runoff')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qsb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Evap:
    IF(output%flux.OR.output%Evap) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Evap to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Evap',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Evap) ! Define Evap
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Evap variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Evap(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Evap to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Evap',NF90_FLOAT,(/landID,tID/), &
               ovid%Evap)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Evap variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Evap(mp,1,1))
       END IF
       ! Define Evap units:
       status = NF90_PUT_ATT(ncid_out,ovid%Evap,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Evap variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Evap,'long_name', &
            'Total evapotranspiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Evap variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ECanop:
    IF(output%flux.OR.output%veg.OR.output%ECanop) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ECanop to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ECanop',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%ECanop) ! Define ECanop
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ECanop variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ECanop(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ECanop to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ECanop',NF90_FLOAT,(/landID,tID/), &
               ovid%ECanop)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ECanop variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ECanop(mp,1,1))
       END IF
       ! Define ECanop units:
       status = NF90_PUT_ATT(ncid_out,ovid%ECanop,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ECanop variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%ECanop,'long_name', &
            'Wet canopy evaporation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ECanop variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! TVeg:
    IF(output%flux.OR.output%veg.OR.output%TVeg) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing TVeg to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'TVeg',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%TVeg) ! Define TVeg
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining TVeg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%TVeg(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing TVeg to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'TVeg',NF90_FLOAT,(/landID,tID/), &
               ovid%TVeg)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining TVeg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%TVeg(mp,1,1))
       END IF
       ! Define TVeg units:
       status = NF90_PUT_ATT(ncid_out,ovid%TVeg,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining TVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%TVeg,'long_name', &
            'Vegetation transpiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining TVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ESoil:
    IF(output%flux.OR.output%soil.OR.output%ESoil) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ESoil to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ESoil',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%ESoil) ! Define ESoil
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ESoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ESoil(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ESoil to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ESoil',NF90_FLOAT,(/landID,tID/), &
               ovid%ESoil)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ESoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ESoil(mp,1,1))
       END IF
       ! Define ESoil units:
       status = NF90_PUT_ATT(ncid_out,ovid%ESoil,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ESoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%ESoil,'long_name', &
            'Evaporation from soil')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ESoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! HVeg:
    IF(output%flux.OR.output%veg.OR.output%HVeg) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing HVeg to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'HVeg',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%HVeg) ! Define HVeg
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HVeg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HVeg(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing HVeg to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'HVeg',NF90_FLOAT,(/landID,tID/), &
               ovid%HVeg)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HVeg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HVeg(mp,1,1))
       END IF
       ! Define HVeg units:
       status = NF90_PUT_ATT(ncid_out,ovid%HVeg,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HVeg,'long_name', &
            'Sensible heat from vegetation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! HSoil:
    IF(output%flux.OR.output%soil.OR.output%HSoil) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing HSoil to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'HSoil',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%HSoil) ! Define HSoil
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HSoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HSoil(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing HSoil to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'HSoil',NF90_FLOAT,(/landID,tID/), &
               ovid%HSoil)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HSoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HSoil(mp,1,1))
       END IF
       ! Define HSoil units:
       status = NF90_PUT_ATT(ncid_out,ovid%HSoil,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HSoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HSoil,'long_name', &
            'Sensible heat from soil')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HSoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! NEE:
    IF(output%flux.OR.output%carbon.OR.output%NEE) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing NEE to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'NEE',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%NEE) ! Define NEE
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NEE variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%NEE(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing NEE to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'NEE',NF90_FLOAT,(/landID,tID/), &
               ovid%NEE)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NEE variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%NEE(mp,1,1))
       END IF
       ! Define NEE units:
       status = NF90_PUT_ATT(ncid_out,ovid%NEE,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NEE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%NEE,'long_name', &
            'Net ecosystem exchange of CO2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NEE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SoilMoist:
    IF(output%soil.OR.output%SoilMoist) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SoilMoist to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SoilMoist',NF90_FLOAT,(/xID,yID,soilID,tID/), &
               ovid%SoilMoist) ! Define SoilMoist
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SoilMoist variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SoilMoist(xdimsize,ydimsize,ms,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SoilMoist to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SoilMoist',NF90_FLOAT,(/landID,soilID,tID/), &
               ovid%SoilMoist)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SoilMoist variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SoilMoist(mp,ms,1,1))
       END IF
       ! Define SoilMoist units:
       status = NF90_PUT_ATT(ncid_out,ovid%SoilMoist,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilMoist variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SoilMoist,'long_name', &
            'Average layer soil moisture')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilMoist variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SoilTemp:
    IF(output%soil.OR.output%SoilTemp) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SoilTemp to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SoilTemp',NF90_FLOAT,(/xID,yID,soilID,tID/), &
               ovid%SoilTemp) ! Define SoilTemp
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SoilTemp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SoilTemp(xdimsize,ydimsize,ms,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SoilTemp to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SoilTemp',NF90_FLOAT,(/landID,soilID,tID/), &
               ovid%SoilTemp)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SoilTemp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SoilTemp(mp,ms,1,1))
       END IF
       ! Define SoilTemp units:
       status = NF90_PUT_ATT(ncid_out,ovid%SoilTemp,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilTemp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SoilTemp,'long_name', &
            'Average layer soil temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilTemp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! BaresoilT:
    IF(output%soil.OR.output%BaresoilT) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing BaresoilT to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'BaresoilT',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%BaresoilT) ! Define BaresoilT
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining BaresoilT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%BaresoilT(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing BaresoilT to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'BaresoilT',NF90_FLOAT,(/landID,tID/), &
               ovid%BaresoilT)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining BaresoilT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%BaresoilT(mp,1,1))
       END IF
       ! Define BaresoilT units:
       status = NF90_PUT_ATT(ncid_out,ovid%BaresoilT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining BaresoilT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%BaresoilT,'long_name', &
            'Bare soil temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining BaresoilT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SWE:
    IF(output%snow.OR.output%SWE) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SWE to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SWE',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%SWE) ! Define SWE
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWE variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWE(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SWE to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SWE',NF90_FLOAT,(/landID,tID/), &
               ovid%SWE)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWE variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWE(mp,1,1))
       END IF
       ! Define SWE units:
       status = NF90_PUT_ATT(ncid_out,ovid%SWE,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWE,'long_name', &
            'Snow water equivalent')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SnowT:
    IF(output%snow.OR.output%SnowT) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SnowT to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SnowT',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%SnowT) ! Define SnowT
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SnowT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SnowT(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SnowT to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SnowT',NF90_FLOAT,(/landID,tID/), &
               ovid%SnowT)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SnowT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SnowT(mp,1,1))
       END IF
       ! Define SnowT units:
       status = NF90_PUT_ATT(ncid_out,ovid%SnowT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowT,'long_name', &
            'Snow surface temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SnowDepth:
    IF(output%snow.OR.output%SnowDepth) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SnowDepth to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SnowDepth',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%SnowDepth) ! Define SnowDepth
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SnowDepth variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SnowDepth(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SnowDepth to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SnowDepth',NF90_FLOAT,(/landID,tID/), &
               ovid%SnowDepth)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SnowDepth variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SnowDepth(mp,1,1))
       END IF
       ! Define SnowDepth units:
       status = NF90_PUT_ATT(ncid_out,ovid%SnowDepth,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowDepth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowDepth,'long_name', &
            'Snow depth')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowDepth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! SWnet:
    IF(output%radiation.OR.output%SWnet) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing SWnet to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'SWnet',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%SWnet) ! Define SWnet
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWnet(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing SWnet to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'SWnet',NF90_FLOAT,(/landID,tID/), &
               ovid%SWnet)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining SWnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%SWnet(mp,1,1))
       END IF
       ! Define SWnet units:
       status = NF90_PUT_ATT(ncid_out,ovid%SWnet,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWnet,'long_name', &
            'Net shortwave radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! LWnet:
    IF(output%radiation.OR.output%LWnet) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing LWnet to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'LWnet',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%LWnet) ! Define LWnet
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LWnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LWnet(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing LWnet to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'LWnet',NF90_FLOAT,(/landID,tID/), &
               ovid%LWnet)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LWnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LWnet(mp,1,1))
       END IF
       ! Define LWnet units:
       status = NF90_PUT_ATT(ncid_out,ovid%LWnet,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWnet,'long_name', &
            'Net longwave radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Rnet:
    IF(output%radiation.OR.output%Rnet) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Rnet to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Rnet',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Rnet) ! Define Rnet
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Rnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Rnet(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Rnet to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Rnet',NF90_FLOAT,(/landID,tID/), &
               ovid%Rnet)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Rnet variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Rnet(mp,1,1))
       END IF
       ! Define Rnet units:
       status = NF90_PUT_ATT(ncid_out,ovid%Rnet,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Rnet,'long_name', &
            'Net radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Albedo:
    IF(output%radiation.OR.output%Albedo) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Albedo to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Albedo',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Albedo) ! Define Albedo
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Albedo variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Albedo(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Albedo to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Albedo',NF90_FLOAT,(/landID,tID/), &
               ovid%Albedo)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Albedo variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Albedo(mp,1,1))
       END IF
       ! Define Albedo units:
       status = NF90_PUT_ATT(ncid_out,ovid%Albedo,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Albedo variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Albedo,'long_name', &
            'Net longwave radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Albedo variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! RadT:
    IF(output%radiation.OR.output%RadT) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing RadT to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'RadT',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%RadT) ! Define RadT
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining RadT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%RadT(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing RadT to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'RadT',NF90_FLOAT,(/landID,tID/), &
               ovid%RadT)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining RadT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%RadT(mp,1,1))
       END IF
       ! Define RadT units:
       status = NF90_PUT_ATT(ncid_out,ovid%RadT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining RadT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%RadT,'long_name', &
            'Radiative surface temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining RadT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! VegT:
    IF(output%veg.OR.output%VegT) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing VegT to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'VegT',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%VegT) ! Define VegT
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining VegT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%VegT(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing VegT to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'VegT',NF90_FLOAT,(/landID,tID/), &
               ovid%VegT)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining VegT variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%VegT(mp,1,1))
       END IF
       ! Define VegT units:
       status = NF90_PUT_ATT(ncid_out,ovid%VegT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining VegT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%VegT,'long_name', &
            'Average vegetation temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining VegT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! CanopInt:
    IF(output%veg.OR.output%CanopInt) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing CanopInt to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'CanopInt',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%CanopInt) ! Define CanopInt
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining CanopInt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%CanopInt(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing CanopInt to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'CanopInt',NF90_FLOAT,(/landID,tID/), &
               ovid%CanopInt)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining CanopInt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%CanopInt(mp,1,1))
       END IF
       ! Define CanopInt units:
       status = NF90_PUT_ATT(ncid_out,ovid%CanopInt,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CanopInt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%CanopInt,'long_name', &
            'Canopy water storage')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CanopInt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! LAI:
    IF(output%veg.OR.output%LAI) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing LAI to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'LAI',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%LAI) ! Define LAI
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LAI variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LAI(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing LAI to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'LAI',NF90_FLOAT,(/landID,tID/), &
               ovid%LAI)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining LAI variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%LAI(mp,1,1))
       END IF
       ! Define LAI units:
       status = NF90_PUT_ATT(ncid_out,ovid%LAI,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LAI variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LAI,'long_name', &
            'Leaf area index')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LAI variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Ebal:
    IF(output%balances.OR.output%Ebal) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Ebal to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Ebal',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Ebal) ! Define Ebal
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Ebal variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Ebal(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Ebal to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Ebal',NF90_FLOAT,(/landID,tID/), &
               ovid%Ebal)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Ebal variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Ebal(mp,1,1))
       END IF
       ! Define Ebal units:
       status = NF90_PUT_ATT(ncid_out,ovid%Ebal,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Ebal variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Ebal,'long_name', &
            'Cumulative energy balance')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Ebal variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! Wbal:
    IF(output%balances.OR.output%Wbal) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing Wbal to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'Wbal',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%Wbal) ! Define Wbal
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Wbal variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Wbal(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing Wbal to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'Wbal',NF90_FLOAT,(/landID,tID/), &
               ovid%Wbal)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining Wbal variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%Wbal(mp,1,1))
       END IF
       ! Define Wbal units:
       status = NF90_PUT_ATT(ncid_out,ovid%Wbal,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wbal variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Wbal,'long_name', &
            'Cumulative water balance')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wbal variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! AutoResp:
    IF(output%carbon.OR.output%AutoResp) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing AutoResp to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'AutoResp',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%AutoResp) ! Define AutoResp
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining AutoResp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%AutoResp(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing AutoResp to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'AutoResp',NF90_FLOAT,(/landID,tID/), &
               ovid%AutoResp)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining AutoResp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%AutoResp(mp,1,1))
       END IF
       ! Define AutoResp units:
       status = NF90_PUT_ATT(ncid_out,ovid%AutoResp,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining AutoResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%AutoResp,'long_name', &
            'Autotrophic respiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining AutoResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! HeteroResp:
    IF(output%carbon.OR.output%HeteroResp) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing HeteroResp to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'HeteroResp',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%HeteroResp) ! Define HeteroResp
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HeteroResp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HeteroResp(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing HeteroResp to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'HeteroResp',NF90_FLOAT,(/landID,tID/), &
               ovid%HeteroResp)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining HeteroResp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%HeteroResp(mp,1,1))
       END IF
       ! Define HeteroResp units:
       status = NF90_PUT_ATT(ncid_out,ovid%HeteroResp,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HeteroResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HeteroResp,'long_name', &
            'Heterotrophic respiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HeteroResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! GPP:
    IF(output%carbon.OR.output%GPP) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing GPP to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'GPP',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%GPP) ! Define GPP
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining GPP variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%GPP(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing GPP to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'GPP',NF90_FLOAT,(/landID,tID/), &
               ovid%GPP)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining GPP variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%GPP(mp,1,1))
       END IF
       ! Define GPP units:
       status = NF90_PUT_ATT(ncid_out,ovid%GPP,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining GPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%GPP,'long_name', &
            'Gross primary production')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining GPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! NPP:
    IF(output%carbon.OR.output%NPP) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing NPP to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'NPP',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%NPP) ! Define NPP
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NPP variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%NPP(xdimsize,ydimsize,1))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing NPP to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'NPP',NF90_FLOAT,(/landID,tID/), &
               ovid%NPP)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NPP variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%NPP(mp,1,1))
       END IF
       ! Define NPP units:
       status = NF90_PUT_ATT(ncid_out,ovid%NPP,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%NPP,'long_name', &
            'Net primary production')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! iveg:
    IF(output%params.OR.output%iveg) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing iveg to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'iveg',NF90_FLOAT,(/xID,yID/), &
               opid%iveg) ! Define iveg
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining iveg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%iveg(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing iveg to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'iveg',NF90_FLOAT,(/landID/), &
               opid%iveg)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining iveg variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%iveg(mp,1))
       END IF
       ! Define iveg units:
       status = NF90_PUT_ATT(ncid_out,opid%iveg,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining iveg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%iveg,'long_name', &
            'Vegetation type')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining iveg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! isoil:
    IF(output%params.OR.output%isoil) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing isoil to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'isoil',NF90_FLOAT,(/xID,yID/), &
               opid%isoil) ! Define isoil
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining isoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%isoil(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing isoil to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'isoil',NF90_FLOAT,(/landID/), &
               opid%isoil)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining isoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%isoil(mp,1))
       END IF
       ! Define isoil units:
       status = NF90_PUT_ATT(ncid_out,opid%isoil,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining isoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%isoil,'long_name', &
            'Soil type')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining isoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! clay:
    IF(output%params.OR.output%clay) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing clay to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'clay',NF90_FLOAT,(/xID,yID/), &
               opid%clay) ! Define clay
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining clay variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%clay(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing clay to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'clay',NF90_FLOAT,(/landID/), &
               opid%clay)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining clay variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%clay(mp,1))
       END IF
       ! Define clay units:
       status = NF90_PUT_ATT(ncid_out,opid%clay,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining clay variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%clay,'long_name', &
            'Fraction of soil which is clay')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining clay variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! sand:
    IF(output%params.OR.output%sand) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing sand to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'sand',NF90_FLOAT,(/xID,yID/), &
               opid%sand) ! Define sand
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sand variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sand(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing sand to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'sand',NF90_FLOAT,(/landID/), &
               opid%sand)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sand variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sand(mp,1))
       END IF
       ! Define sand units:
       status = NF90_PUT_ATT(ncid_out,opid%sand,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sand variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sand,'long_name', &
            'Fraction of soil which is sand')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sand variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! silt:
    IF(output%params.OR.output%silt) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing silt to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'silt',NF90_FLOAT,(/xID,yID/), &
               opid%silt) ! Define silt
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining silt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%silt(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing silt to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'silt',NF90_FLOAT,(/landID/), &
               opid%silt)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining silt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%silt(mp,1))
       END IF
       ! Define silt units:
       status = NF90_PUT_ATT(ncid_out,opid%silt,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining silt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%silt,'long_name', &
            'Fraction of soil which is silt')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining silt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ssat:
    IF(output%params.OR.output%ssat) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ssat to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ssat',NF90_FLOAT,(/xID,yID/), &
               opid%ssat) ! Define ssat
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ssat variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ssat(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ssat to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ssat',NF90_FLOAT,(/landID/), &
               opid%ssat)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ssat variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ssat(mp,1))
       END IF
       ! Define ssat units:
       status = NF90_PUT_ATT(ncid_out,opid%ssat,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ssat variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ssat,'long_name', &
            'Fraction of soil volume which is water @ saturation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ssat variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! sfc:
    IF(output%params.OR.output%sfc) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing sfc to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'sfc',NF90_FLOAT,(/xID,yID/), &
               opid%sfc) ! Define sfc
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sfc variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sfc(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing sfc to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'sfc',NF90_FLOAT,(/landID/), &
               opid%sfc)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sfc variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sfc(mp,1))
       END IF
       ! Define sfc units:
       status = NF90_PUT_ATT(ncid_out,opid%sfc,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sfc variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sfc,'long_name', &
            'Fraction of soil volume which is water @ field capacity')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sfc variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! swilt:
    IF(output%params.OR.output%swilt) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing swilt to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'swilt',NF90_FLOAT,(/xID,yID/), &
               opid%swilt) ! Define swilt
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining swilt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%swilt(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing swilt to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'swilt',NF90_FLOAT,(/landID/), &
               opid%swilt)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining swilt variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%swilt(mp,1))
       END IF
       ! Define swilt units:
       status = NF90_PUT_ATT(ncid_out,opid%swilt,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining swilt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%swilt,'long_name', &
            'Fraction of soil volume which is water @ wilting point')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining swilt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! zse:
    IF(output%params.OR.output%zse) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing zse to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'zse',NF90_FLOAT,(/xID,yID,soilID/), &
               opid%zse) ! Define zse
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining zse variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%zse(xdimsize,ydimsize,ms))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing zse to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'zse',NF90_FLOAT,(/landID/), &
               opid%zse)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining zse variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%zse(mp,ms,1))
       END IF
       ! Define zse units:
       status = NF90_PUT_ATT(ncid_out,opid%zse,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining zse variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%zse,'long_name', &
            'Depth of each soil layer')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining zse variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! froot:
    IF(output%params.OR.output%froot) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing froot to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'froot',NF90_FLOAT,(/xID,yID,soilID/), &
               opid%froot) ! Define froot
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining froot variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%froot(xdimsize,ydimsize,ms))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing froot to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'froot',NF90_FLOAT,(/landID/), &
               opid%froot)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining froot variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%froot(mp,ms,1))
       END IF
       ! Define froot units:
       status = NF90_PUT_ATT(ncid_out,opid%froot,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining froot variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%froot,'long_name', &
            'Fraction of roots in each soil layer')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining froot variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! bch:
    IF(output%params.OR.output%bch) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing bch to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'bch',NF90_FLOAT,(/xID,yID/), &
               opid%bch) ! Define bch
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining bch variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%bch(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing bch to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'bch',NF90_FLOAT,(/landID/), &
               opid%bch)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining bch variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%bch(mp,1))
       END IF
       ! Define bch units:
       status = NF90_PUT_ATT(ncid_out,opid%bch,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining bch variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%bch,'long_name', &
            'Parameter b, Campbell eqn 1985')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining bch variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! hyds:
    IF(output%params.OR.output%hyds) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing hyds to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'hyds',NF90_FLOAT,(/xID,yID/), &
               opid%hyds) ! Define hyds
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining hyds variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%hyds(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing hyds to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'hyds',NF90_FLOAT,(/landID/), &
               opid%hyds)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining hyds variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%hyds(mp,1))
       END IF
       ! Define hyds units:
       status = NF90_PUT_ATT(ncid_out,opid%hyds,'units','m/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hyds variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%hyds,'long_name', &
            'Hydraulic conductivity @ saturation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hyds variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! sucs:
    IF(output%params.OR.output%sucs) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing sucs to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'sucs',NF90_FLOAT,(/xID,yID/), &
               opid%sucs) ! Define sucs
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sucs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sucs(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing sucs to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'sucs',NF90_FLOAT,(/landID/), &
               opid%sucs)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining sucs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%sucs(mp,1))
       END IF
       ! Define sucs units:
       status = NF90_PUT_ATT(ncid_out,opid%sucs,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sucs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sucs,'long_name', &
            'Suction @ saturation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sucs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! css:
    IF(output%params.OR.output%css) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing css to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'css',NF90_FLOAT,(/xID,yID/), &
               opid%css) ! Define css
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining css variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%css(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing css to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'css',NF90_FLOAT,(/landID/), &
               opid%css)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining css variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%css(mp,1))
       END IF
       ! Define css units:
       status = NF90_PUT_ATT(ncid_out,opid%css,'units','J/kg/C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining css variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%css,'long_name', &
            'Heat capacity of soil minerals')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining css variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! rhosoil:
    IF(output%params.OR.output%rhosoil) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing rhosoil to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'rhosoil',NF90_FLOAT,(/xID,yID/), &
               opid%rhosoil) ! Define rhosoil
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rhosoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rhosoil(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing rhosoil to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'rhosoil',NF90_FLOAT,(/landID/), &
               opid%rhosoil)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rhosoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rhosoil(mp,1))
       END IF
       ! Define rhosoil units:
       status = NF90_PUT_ATT(ncid_out,opid%rhosoil,'units','kg/m^3')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rhosoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rhosoil,'long_name', &
            'Density of soil minerals')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rhosoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! rs20:
    IF(output%params.OR.output%rs20) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing rs20 to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'rs20',NF90_FLOAT,(/xID,yID/), &
               opid%rs20) ! Define rs20
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rs20 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rs20(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing rs20 to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'rs20',NF90_FLOAT,(/landID/), &
               opid%rs20)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rs20 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rs20(mp,1))
       END IF
       ! Define rs20 units:
       status = NF90_PUT_ATT(ncid_out,opid%rs20,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rs20 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rs20,'long_name', &
            'Soil respiration coefficient at 20C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rs20 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! albsoil:
    IF(output%params.OR.output%albsoil) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing albsoil to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'albsoil',NF90_FLOAT,(/xID,yID/), &
               opid%albsoil) ! Define albsoil
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining albsoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%albsoil(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing albsoil to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'albsoil',NF90_FLOAT,(/landID/), &
               opid%albsoil)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining albsoil variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%albsoil(mp,1))
       END IF
       ! Define albsoil units:
       status = NF90_PUT_ATT(ncid_out,opid%albsoil,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining albsoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%albsoil,'long_name', &
            'Snow free shortwave soil reflectance fraction')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining albsoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! hc:
    IF(output%params.OR.output%hc) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing hc to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'hc',NF90_FLOAT,(/xID,yID/), &
               opid%hc) ! Define hc
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining hc variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%hc(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing hc to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'hc',NF90_FLOAT,(/landID/), &
               opid%hc)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining hc variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%hc(mp,1))
       END IF
       ! Define hc units:
       status = NF90_PUT_ATT(ncid_out,opid%hc,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hc variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%hc,'long_name', &
            'Height of canopy')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hc variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! canst1:
    IF(output%params.OR.output%canst1) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing canst1 to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'canst1',NF90_FLOAT,(/xID,yID/), &
               opid%canst1) ! Define canst1
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining canst1 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%canst1(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing canst1 to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'canst1',NF90_FLOAT,(/landID/), &
               opid%canst1)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining canst1 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%canst1(mp,1))
       END IF
       ! Define canst1 units:
       status = NF90_PUT_ATT(ncid_out,opid%canst1,'units','mm/LAI')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining canst1 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%canst1,'long_name', &
            'Max water intercepted by canopy')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining canst1 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! dleaf:
    IF(output%params.OR.output%dleaf) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing dleaf to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'dleaf',NF90_FLOAT,(/xID,yID/), &
               opid%dleaf) ! Define dleaf
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining dleaf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%dleaf(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing dleaf to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'dleaf',NF90_FLOAT,(/landID/), &
               opid%dleaf)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining dleaf variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%dleaf(mp,1))
       END IF
       ! Define dleaf units:
       status = NF90_PUT_ATT(ncid_out,opid%dleaf,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining dleaf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%dleaf,'long_name', &
            'Chararacteristic length of leaf')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining dleaf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! frac4:
    IF(output%params.OR.output%frac4) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing frac4 to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'frac4',NF90_FLOAT,(/xID,yID/), &
               opid%frac4) ! Define frac4
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining frac4 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%frac4(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing frac4 to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'frac4',NF90_FLOAT,(/landID/), &
               opid%frac4)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining frac4 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%frac4(mp,1))
       END IF
       ! Define frac4 units:
       status = NF90_PUT_ATT(ncid_out,opid%frac4,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining frac4 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%frac4,'long_name', &
            'Fraction of plants which are C4')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining frac4 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ejmax:
    IF(output%params.OR.output%ejmax) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ejmax to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ejmax',NF90_FLOAT,(/xID,yID/), &
               opid%ejmax) ! Define ejmax
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ejmax variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ejmax(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ejmax to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ejmax',NF90_FLOAT,(/landID/), &
               opid%ejmax)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ejmax variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ejmax(mp,1))
       END IF
       ! Define ejmax units:
       status = NF90_PUT_ATT(ncid_out,opid%ejmax,'units','mol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ejmax variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ejmax,'long_name', &
            'Max potential electron transport rate top leaf')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ejmax variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! vcmax:
    IF(output%params.OR.output%vcmax) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing vcmax to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'vcmax',NF90_FLOAT,(/xID,yID/), &
               opid%vcmax) ! Define vcmax
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining vcmax variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%vcmax(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing vcmax to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'vcmax',NF90_FLOAT,(/landID/), &
               opid%vcmax)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining vcmax variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%vcmax(mp,1))
       END IF
       ! Define vcmax units:
       status = NF90_PUT_ATT(ncid_out,opid%vcmax,'units','mol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vcmax variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%vcmax,'long_name', &
            'Maximum RuBP carboxylation rate top leaf')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vcmax variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! rp20:
    IF(output%params.OR.output%rp20) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing rp20 to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'rp20',NF90_FLOAT,(/xID,yID/), &
               opid%rp20) ! Define rp20
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rp20 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rp20(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing rp20 to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'rp20',NF90_FLOAT,(/landID/), &
               opid%rp20)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rp20 variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rp20(mp,1))
       END IF
       ! Define rp20 units:
       status = NF90_PUT_ATT(ncid_out,opid%rp20,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rp20 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rp20,'long_name', &
            'Plant respiration coefficient at 20C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rp20 variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! rpcoef:
    IF(output%params.OR.output%rpcoef) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing rpcoef to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'rpcoef',NF90_FLOAT,(/xID,yID/), &
               opid%rpcoef) ! Define rpcoef
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rpcoef variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rpcoef(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing rpcoef to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'rpcoef',NF90_FLOAT,(/landID/), &
               opid%rpcoef)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining rpcoef variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%rpcoef(mp,1))
       END IF
       ! Define rpcoef units:
       status = NF90_PUT_ATT(ncid_out,opid%rpcoef,'units','1/C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rpcoef variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rpcoef,'long_name', &
            'Temperature coef nonleaf plant respiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rpcoef variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! shelrb:
    IF(output%params.OR.output%shelrb) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing shelrb to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'shelrb',NF90_FLOAT,(/xID,yID/), &
               opid%shelrb) ! Define shelrb
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining shelrb variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%shelrb(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing shelrb to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'shelrb',NF90_FLOAT,(/landID/), &
               opid%shelrb)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining shelrb variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%shelrb(mp,1))
       END IF
       ! Define shelrb units:
       status = NF90_PUT_ATT(ncid_out,opid%shelrb,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining shelrb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%shelrb,'long_name', &
            'Sheltering factor')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining shelrb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! xfang:
    IF(output%params.OR.output%xfang) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing xfang to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'xfang',NF90_FLOAT,(/xID,yID/), &
               opid%xfang) ! Define xfang
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining xfang variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%xfang(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing xfang to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'xfang',NF90_FLOAT,(/landID/), &
               opid%xfang)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining xfang variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%xfang(mp,1))
       END IF
       ! Define xfang units:
       status = NF90_PUT_ATT(ncid_out,opid%xfang,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining xfang variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%xfang,'long_name', &
            'Leaf angle parameter')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining xfang variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! tminvj:
    IF(output%params.OR.output%tminvj) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing tminvj to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'tminvj',NF90_FLOAT,(/xID,yID/), &
               opid%tminvj) ! Define tminvj
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining tminvj variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%tminvj(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing tminvj to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'tminvj',NF90_FLOAT,(/landID/), &
               opid%tminvj)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining tminvj variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%tminvj(mp,1))
       END IF
       ! Define tminvj units:
       status = NF90_PUT_ATT(ncid_out,opid%tminvj,'units','C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tminvj variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%tminvj,'long_name', &
            'Min temperature for the start of photosynthesis')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tminvj variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! tmaxvj:
    IF(output%params.OR.output%tmaxvj) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing tmaxvj to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'tmaxvj',NF90_FLOAT,(/xID,yID/), &
               opid%tmaxvj) ! Define tmaxvj
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining tmaxvj variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%tmaxvj(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing tmaxvj to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'tmaxvj',NF90_FLOAT,(/landID/), &
               opid%tmaxvj)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining tmaxvj variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%tmaxvj(mp,1))
       END IF
       ! Define tmaxvj units:
       status = NF90_PUT_ATT(ncid_out,opid%tmaxvj,'units','C')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tmaxvj variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%tmaxvj,'long_name', &
            'Max temperature for photosynthesis')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tmaxvj variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! vbeta:
    IF(output%params.OR.output%vbeta) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing vbeta to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'vbeta',NF90_FLOAT,(/xID,yID/), &
               opid%vbeta) ! Define vbeta
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining vbeta variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%vbeta(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing vbeta to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'vbeta',NF90_FLOAT,(/landID/), &
               opid%vbeta)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining vbeta variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%vbeta(mp,1))
       END IF
       ! Define vbeta units:
       status = NF90_PUT_ATT(ncid_out,opid%vbeta,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vbeta variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%vbeta,'long_name', &
            'Stomatal sensitivity to soil water')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vbeta variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ratecp:
    IF(output%params.OR.output%ratecp) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ratecp to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ratecp',NF90_FLOAT,(/xID,yID,plantcarbID/), &
               opid%ratecp) ! Define ratecp
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ratecp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ratecp(xdimsize,ydimsize,ncp))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ratecp to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ratecp',NF90_FLOAT,(/landID/), &
               opid%ratecp)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ratecp variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ratecp(mp,ncp,1))
       END IF
       ! Define ratecp units:
       status = NF90_PUT_ATT(ncid_out,opid%ratecp,'units','1/year')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ratecp,'long_name', &
            'Plant carbon rate constant')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! ratecs:
    IF(output%params.OR.output%ratecs) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing ratecs to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'ratecs',NF90_FLOAT,(/xID,yID,soilcarbID/), &
               opid%ratecs) ! Define ratecs
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ratecs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ratecs(xdimsize,ydimsize,ncs))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing ratecs to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'ratecs',NF90_FLOAT,(/landID,soilcarbID/), &
               opid%ratecs)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining ratecs variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%ratecs(mp,ncs,1))
       END IF
       ! Define ratecs units:
       status = NF90_PUT_ATT(ncid_out,opid%ratecs,'units','1/year')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ratecs,'long_name', &
            'Soil carbon rate constant')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! meth:
    IF(output%params.OR.output%meth) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing meth to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'meth',NF90_FLOAT,(/xID,yID/), &
               opid%meth) ! Define meth
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining meth variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%meth(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing meth to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'meth',NF90_FLOAT,(/landID/), &
               opid%meth)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining meth variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%meth(mp,1))
       END IF
       ! Define meth units:
       status = NF90_PUT_ATT(ncid_out,opid%meth,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining meth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%meth,'long_name', &
            'Canopy turbulence parameterisation choice')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining meth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    ! za:
    IF(output%params.OR.output%za) THEN
       IF(gridType=='mask') THEN
          WRITE(logn,*) 'Writing za to output file using mask grid'
          status=NF90_DEF_VAR(ncid_out,'za',NF90_FLOAT,(/xID,yID/), &
               opid%za) ! Define za
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining za variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%za(xdimsize,ydimsize))
       ELSE IF(gridType=='land') THEN
          WRITE(logn,*) 'Writing za to output file using land grid'
          status=NF90_DEF_VAR(ncid_out,'za',NF90_FLOAT,(/landID/), &
               opid%za)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining za variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%za(mp,1))
       END IF
       ! Define za units:
       status = NF90_PUT_ATT(ncid_out,opid%za,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining za variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%za,'long_name', &
            'Reference height (lowest atm. model layer)')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining za variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
    END IF
    
    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    status = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"Production", &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_out)// ' (SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"Source", &
         'CABLE LSM output file')
     IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_out)// ' (SUBROUTINE open_output_file)')
    status = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"CABLE_input_file", &
         TRIM(filename_met))
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_out)// ' (SUBROUTINE open_output_file)')
    
    ! End netcdf define mode:
    status = NF90_ENDDEF(ncid_out)
    IF(status/=NF90_NOERR) CALL nc_abort('Error creating output file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
   
    ! Write time variable:
    status=NF90_PUT_VAR(ncid_out,tvarID,timevar)
    IF(status/=NF90_NOERR) CALL nc_abort('Error time variable to file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    ! Write latitude and longitude variables:
    status=NF90_PUT_VAR(ncid_out,latID,lat_all)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing latitude variable to file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    status=NF90_PUT_VAR(ncid_out,lonID,lon_all)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing longitude variable to file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    ! Write GrADS coordinate variables
    status=NF90_PUT_VAR(ncid_out,xvID,lon_all(:,1))
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error writing GrADS x coordinate variable to file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    status=NF90_PUT_VAR(ncid_out,yvID,lat_all(1,:))
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error writing GrADS y coordinate variable to file ' &
         //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
    ! Write model parameters if requested:
    IF(output%params) THEN
       IF(gridType=='mask') THEN 
          ! Define global grid parameter variables
          DO i = 1, mp
             out%iveg(land_x(i),land_y(i)) = veg%iveg(i)
             out%isoil(land_x(i),land_y(i)) = soil%isoilm(i)
             out%bch(land_x(i),land_y(i)) = soil%bch(i)
             out%clay(land_x(i),land_y(i)) = soil%clay(i)
             out%sand(land_x(i),land_y(i)) = soil%sand(i)
             out%silt(land_x(i),land_y(i)) = soil%silt(i)
             out%css(land_x(i),land_y(i)) = soil%css(i)
             out%rhosoil(land_x(i),land_y(i))= soil%rhosoil(i)
             out%hyds(land_x(i),land_y(i)) = soil%hyds(i)
             out%sucs(land_x(i),land_y(i)) = soil%sucs(i)
             out%rs20(land_x(i),land_y(i)) = soil%rs20(i)
             out%ssat(land_x(i),land_y(i)) = soil%ssat(i)
             out%sfc(land_x(i),land_y(i)) = soil%sfc(i)
             out%swilt(land_x(i),land_y(i)) = soil%swilt(i)
             out%froot(land_x(i),land_y(i),:) = veg%froot(i,:)
             out%zse(land_x(i),land_y(i),:) = soil%zse
             out%albsoil(land_x(i),land_y(i)) = soil%albsoil(i)
             out%canst1(land_x(i),land_y(i)) = veg%canst1(i)
             out%dleaf(land_x(i),land_y(i)) = veg%dleaf(i)
             out%ejmax(land_x(i),land_y(i)) = veg%ejmax(i)
             out%vcmax(land_x(i),land_y(i)) = veg%vcmax(i)
             out%frac4(land_x(i),land_y(i)) = veg%frac4(i)
             out%hc(land_x(i),land_y(i)) = veg%hc(i)
             out%rp20(land_x(i),land_y(i)) = veg%rp20(i)
             out%rpcoef(land_x(i),land_y(i)) = veg%rpcoef(i)
             out%shelrb(land_x(i),land_y(i)) = veg%shelrb(i)
             out%xfang(land_x(i),land_y(i)) = veg%xfang(i)
             out%tminvj(land_x(i),land_y(i)) = veg%tminvj(i)
             out%tmaxvj(land_x(i),land_y(i)) = veg%tmaxvj(i)
             out%vbeta(land_x(i),land_y(i)) = veg%vbeta(i)
             out%ratecp(land_x(i),land_y(i),:) = bgc%ratecp(:) !no spatial dim at present
             out%ratecs(land_x(i),land_y(i),:) = bgc%ratecs(:) ! no spatial dim at present
             out%meth(land_x(i),land_y(i)) = veg%meth(i)
             out%za(land_x(i),land_y(i)) = rough%za(i)
          END DO
          ! Write parameters to netcdf file:
          status=NF90_PUT_VAR(ncid_out,opid%iveg,out%iveg,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing iveg parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%isoil,out%isoil,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing isoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%bch,out%bch,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing bch parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%clay,out%clay,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing clay parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sand,out%sand,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sand parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%silt,out%silt,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing silt parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%css,out%css,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing css parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rhosoil,out%rhosoil,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rhosoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%hyds,out%hyds,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing hyds parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sucs,out%sucs,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sucs parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rs20,out%rs20,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rs20 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ssat,out%ssat,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ssat parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sfc,out%sfc,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sfc parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%swilt,out%swilt,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing swilt parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%froot,out%froot,start=(/1,1,1/), &
               count=(/xdimsize,ydimsize,ms/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing froot parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%zse,out%zse,start=(/1,1,1/), &
               count=(/xdimsize,ydimsize,ms/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing zse parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%albsoil,out%albsoil,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing albsoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%canst1,out%canst1,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing canst1 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%dleaf,out%dleaf,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing dleaf parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ejmax,out%ejmax,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ejmax parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%vcmax,out%vcmax,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing vcmax parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%frac4,out%frac4,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing frac4 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%hc,out%hc,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing hc parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rp20,out%rp20,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rp20 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rpcoef,out%rpcoef,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rpcoef parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%shelrb,out%shelrb,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing shelrb parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%xfang,out%xfang,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing xfang parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%tminvj,out%tminvj,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing tminvj parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%tmaxvj,out%tmaxvj,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing tmaxvj parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%vbeta,out%vbeta,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing vbeta parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ratecp,out%ratecp,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecp parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ratecs,out%ratecs,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecs parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%meth,out%meth,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing meth parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%za,out%za,start=(/1,1/), &
               count=(/xdimsize,ydimsize/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing za parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')

       ELSE IF(gridType=='land') THEN
          DO i = 1, mp
             out%iveg(i,1) = veg%iveg(i)
             out%isoil(i,1) = soil%isoilm(i)
             out%bch(i,1) = soil%bch(i)
             out%clay(i,1) = soil%clay(i)
             out%sand(i,1) = soil%sand(i)
             out%silt(i,1) = soil%silt(i)
             out%css(i,1) = soil%css(i)
             out%rhosoil(i,1)= soil%rhosoil(i)
             out%hyds(i,1) = soil%hyds(i)
             out%sucs(i,1) = soil%sucs(i)
             out%rs20(i,1) = soil%rs20(i)
             out%ssat(i,1) = soil%ssat(i)
             out%sfc(i,1) = soil%sfc(i)
             out%swilt(i,1) = soil%swilt(i)
             out%froot(i,:,1) = veg%froot(i,:)
             out%zse(i,:,1) = soil%zse
             out%albsoil(i,1) = soil%albsoil(i)
             out%canst1(i,1) = veg%canst1(i)
             out%dleaf(i,1) = veg%dleaf(i)
             out%ejmax(i,1) = veg%ejmax(i)
             out%vcmax(i,1) = veg%vcmax(i)
             out%frac4(i,1) = veg%frac4(i)
             out%hc(i,1) = veg%hc(i)
             out%rp20(i,1) = veg%rp20(i)
             out%rpcoef(i,1) = veg%rpcoef(i)
             out%shelrb(i,1) = veg%shelrb(i)
             out%xfang(i,1) = veg%xfang(i)
             out%tminvj(i,1) = veg%tminvj(i)
             out%tmaxvj(i,1) = veg%tmaxvj(i)
             out%vbeta(i,1) = veg%vbeta(i)
             out%ratecp(i,:,1) = bgc%ratecp(:) !no spatial dim at present
             out%ratecs(i,:,1) = bgc%ratecs(:) ! no spatial dim at present
             out%meth(i,1) = veg%meth(i)
             out%za(i,1) = rough%za(i)
          END DO
          ! Write parameters to netcdf file:
          status=NF90_PUT_VAR(ncid_out,opid%iveg,out%iveg,start=(/1/), &
               count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing iveg parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%isoil,out%isoil(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing isoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%bch,out%bch(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing bch parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%clay,out%clay(:,1),start=(/1/), count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing clay parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sand,out%sand(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sand parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%silt,out%silt(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing silt parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%css,out%css(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing css parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rhosoil,out%rhosoil(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rhosoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%hyds,out%hyds(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing hyds parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sucs,out%sucs(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sucs parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rs20,out%rs20(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rs20 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ssat,out%ssat(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ssat parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%sfc,out%sfc(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing sfc parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%swilt,out%swilt(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing swilt parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%froot,out%froot(:,:,1),start=(/1,1/), &
               count=(/mp,ms/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing froot parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%zse,out%zse(:,:,1),start=(/1,1/), &
               count=(/mp,ms/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing zse parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%albsoil,out%albsoil(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing albsoil parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%canst1,out%canst1(:,1),start=(/1/), count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing canst1 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%dleaf,out%dleaf(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing dleaf parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ejmax,out%ejmax(:,1),start=(/1/), count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ejmax parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%vcmax,out%vcmax(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing vcmax parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%frac4,out%frac4(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing frac4 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%hc,out%hc(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing hc parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rp20,out%rp20(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rp20 parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%rpcoef,out%rpcoef(:,1),start=(/1/), count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing rpcoef parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%shelrb,out%shelrb(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing shelrb parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%xfang,out%xfang(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing xfang parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%tminvj,out%tminvj(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing tminvj parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%tmaxvj,out%tmaxvj(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing tmaxvj parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%vbeta,out%vbeta(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing vbeta parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ratecp,out%ratecp(:,:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecp parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%ratecs,out%ratecs(:,:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecs parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%meth,out%meth(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing meth parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
          status=NF90_PUT_VAR(ncid_out,opid%za,out%za(:,1),start=(/1/),count=(/mp/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing za parameter to file ' &
               //TRIM(filename_out)// '(SUBROUTINE open_output_file)')
       END IF ! grid type, land/sea or compression land only

       ! Deallocate output parameter variables, since they've now been written to file:
       DEALLOCATE(out%bch,out%clay,out%sand,out%silt,out%css,out%rhosoil, &
            out%hyds,out%rs20,out%ssat,out%sfc,out%swilt,out%sucs,out%froot,out%zse, &
            out%canst1,out%dleaf,out%ejmax,out%vcmax,out%frac4,out%hc,out%rp20, &
            out%rpcoef,out%shelrb,out%xfang,out%albsoil,out%tminvj,out%tmaxvj, &
            out%vbeta,out%ratecp,out%ratecs,out%meth,out%za)
    END IF
       
  END SUBROUTINE open_output_file
  !=====================================================================
  SUBROUTINE write_output(ktau,dels,filename_out,met,canopy,ssoil,rad, &
       bal,air,soil,veg)
    ! Writes model outputs; called each timestep.
    INTEGER(i_d), INTENT(IN)          :: ktau ! time step
    REAL(r_1),INTENT(IN)              :: dels ! time step size
    CHARACTER(LEN=*), INTENT(IN) :: filename_out ! output file name
    TYPE(met_type),INTENT(IN) :: met  ! met data
    TYPE (canopy_type),INTENT(IN):: canopy ! canopy variable data
    TYPE (soil_snow_type),INTENT(IN)  :: ssoil ! soil data
    TYPE (soil_parameter_type),INTENT(IN)  :: soil ! soil data
    TYPE (radiation_type),INTENT(IN)  :: rad   ! radiation data
    TYPE (balances_type),INTENT(INOUT):: bal 
    TYPE (air_type),INTENT(IN) 	:: air
    TYPE (veg_parameter_type),INTENT(IN) :: veg
    INTEGER :: i,j ! do loop counter
    
    ! IF asked to check mass/water balance:
    IF(check%mass_bal) CALL mass_balance(ktau,dels,ssoil,soil,canopy, &
         met,air,bal)

    ! IF asked to check energy balance:
    IF(check%energy_bal) CALL energy_balance(ktau,dels,met,rad, &
         canopy,bal,ssoil)
    
    !-------------------WRITE MET DATA-----------------------------
    ! SWdown:
    IF(output%met.OR.output%SWdown) THEN ! If SWdown is requested output
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SWdown(land_x(i),land_y(i),1) = met%fsd(i)
          END DO
          WHERE(mask/=1) out%SWdown(:,:,1) = -9.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SWdown,out%SWdown,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%SWdown(:,1,1) = met%fsd
          status=NF90_PUT_VAR(ncid_out,ovid%SWdown,out%SWdown(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       ! Check writing was successful:
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWdown variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! LWdown:
    IF(output%met.OR.output%LWdown) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%LWdown(land_x(i),land_y(i),1) = met%fld(i)
          END DO
          WHERE(mask/=1) out%LWdown(:,:,1) = 100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%LWdown,out%LWdown,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%LWdown(:,1,1) = met%fld
          status=NF90_PUT_VAR(ncid_out,ovid%LWdown,out%LWdown(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LWdown variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Rainf:
    IF(output%met.OR.output%Rainf) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Rainf(land_x(i),land_y(i),1) = met%precip(i)/dels
          END DO
          WHERE(mask/=1) out%Rainf(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Rainf,out%Rainf,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Rainf(:,1,1) = met%precip/dels
          status=NF90_PUT_VAR(ncid_out,ovid%Rainf,out%Rainf(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Rainf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Snowf:
    IF(output%met.OR.output%Snowf) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Snowf(land_x(i),land_y(i),1) = met%precips(i)/dels
          END DO
          WHERE(mask/=1) out%Snowf(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Snowf,out%Snowf,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Snowf(:,1,1) = met%precips/dels
          status=NF90_PUT_VAR(ncid_out,ovid%Snowf,out%Snowf(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Snowf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! PSurf:
    IF(output%met.OR.output%PSurf) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%PSurf(land_x(i),land_y(i),1) = met%pmb(i)
          END DO
          WHERE(mask/=1) out%PSurf(:,:,1) = 1200.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%PSurf,out%PSurf,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%PSurf(:,1,1) = met%pmb
          status=NF90_PUT_VAR(ncid_out,ovid%PSurf,out%PSurf(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing PSurf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Tair:
    IF(output%met.OR.output%Tair) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Tair(land_x(i),land_y(i),1,1) = met%tk(i)
          END DO
          WHERE(mask/=1) out%Tair(:,:,1,1) = 180.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Tair,out%Tair,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,1,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Tair(:,1,1,1) = met%tk
          status=NF90_PUT_VAR(ncid_out,ovid%Tair,out%Tair(:,:,1,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Tair variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Qair:
    IF(output%met.OR.output%Qair) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qair(land_x(i),land_y(i),1,1) = met%qv(i)
          END DO
          WHERE(mask/=1) out%Qair(:,:,1,1) = 0.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qair,out%Qair,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,1,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qair(:,1,1,1) = met%qv
          status=NF90_PUT_VAR(ncid_out,ovid%Qair,out%Qair(:,:,1,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qair variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Wind:
    IF(output%met.OR.output%Wind) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Wind(land_x(i),land_y(i),1,1) = met%ua(i)
          END DO
          WHERE(mask/=1) out%Wind(:,:,1,1) = -1.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Wind,out%Wind,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,1,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Wind(:,1,1,1) = met%ua
          status=NF90_PUT_VAR(ncid_out,ovid%Wind,out%Wind(:,:,1,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Wind variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! CO2air:
    IF(output%met.OR.output%CO2air) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%CO2air(land_x(i),land_y(i),1,1) = met%ca(i)*1000000.0
          END DO
          WHERE(mask/=1) out%CO2air(:,:,1,1) = 100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%CO2air,out%CO2air,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,1,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%CO2air(:,1,1,1) = met%ca*1000000.0
          status=NF90_PUT_VAR(ncid_out,ovid%CO2air,out%CO2air(:,:,1,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing CO2air variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    !--------------------------WRITE FLUX DATA-------------------------------------
    ! Qle:
    IF(output%flux.OR.output%Qle) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qle(land_x(i),land_y(i),1) = canopy%fe(i) 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Qle(land_x(i),land_y(i),1)<ranges%Qle(1)).OR. &
                     (out%Qle(land_x(i),land_y(i),1)>ranges%Qle(2))) &
                     CALL range_abort('Qle out of specified ranges!', ktau, met, &
                     out%Qle(land_x(i),land_y(i),1),ranges%Qle,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Qle(:,:,1) = -1000.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qle,out%Qle,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qle(:,1,1) = canopy%fe
          status=NF90_PUT_VAR(ncid_out,ovid%Qle,out%Qle(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Qle(i,1,1)<ranges%Qle(1)).OR.(out%Qle(i,1,1)>ranges%Qle(2))) &
                     CALL range_abort('Qle out of specified ranges!', ktau, met, &
                     out%Qle(i,1,1),ranges%Qle,i)
             END DO
          END IF
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qle variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Qh:
    IF(output%flux.OR.output%Qh) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qh(land_x(i),land_y(i),1) = canopy%fh(i)  
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Qh(land_x(i),land_y(i),1)<ranges%Qh(1)).OR. &
                     (out%Qh(land_x(i),land_y(i),1)>ranges%Qh(2))) &
                     CALL range_abort('Qh out of specified ranges!', ktau, met, &
                     out%Qh(land_x(i),land_y(i),1),ranges%Qh,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Qh(:,:,1) = -1000.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qh,out%Qh,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qh(:,1,1) = canopy%fh
          status=NF90_PUT_VAR(ncid_out,ovid%Qh,out%Qh(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Qh(i,1,1)<ranges%Qh(1)).OR.(out%Qh(i,1,1)>ranges%Qh(2))) &
                     CALL range_abort('Qh out of specified ranges!', ktau, met, &
                     out%Qh(i,1,1),ranges%Qh,i)
             END DO
          END IF
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qh variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Qg:
    IF(output%flux.OR.output%Qg) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qg(land_x(i),land_y(i),1) = canopy%ga(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Qg(land_x(i),land_y(i),1)<ranges%Qg(1)).OR. &
                     (out%Qg(land_x(i),land_y(i),1)>ranges%Qg(2))) &
                     CALL range_abort('Qg out of specified ranges!', ktau, met, &
                     out%Qg(land_x(i),land_y(i),1),ranges%Qg,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Qg(:,:,1) = -1000.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qg,out%Qg,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qg(:,1,1) = canopy%ga
          status=NF90_PUT_VAR(ncid_out,ovid%Qg,out%Qg(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Qg(i,1,1)<ranges%Qg(1)).OR.(out%Qg(i,1,1)>ranges%Qg(2))) &
                     CALL range_abort('Qg out of specified ranges!', ktau, met, &
                     out%Qg(i,1,1),ranges%Qg,i)
             END DO
          END IF
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Qs:
    IF(output%flux.OR.output%Qs) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qs(land_x(i),land_y(i),1) = ssoil%rnof1(i)/dels 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Qs(land_x(i),land_y(i),1)<ranges%Qs(1)).OR. &
                     (out%Qs(land_x(i),land_y(i),1)>ranges%Qs(2))) &
                     CALL range_abort('Qs out of specified ranges!', ktau, met, &
                     out%Qs(land_x(i),land_y(i),1),ranges%Qs,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Qs(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qs,out%Qs,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qs(:,1,1) = ssoil%rnof1/dels
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Qs(i,1,1)<ranges%Qs(1)).OR.(out%Qs(i,1,1)>ranges%Qs(2))) &
                     CALL range_abort('Qs out of specified ranges!', ktau, met, &
                     out%Qs(i,1,1),ranges%Qs,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%Qs,out%Qs(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qs variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Qsb:
    IF(output%flux.OR.output%Qsb) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Qsb(land_x(i),land_y(i),1) = ssoil%rnof2(i)/dels 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Qsb(land_x(i),land_y(i),1)<ranges%Qsb(1)).OR. &
                     (out%Qsb(land_x(i),land_y(i),1)>ranges%Qsb(2))) &
                     CALL range_abort('Qsb out of specified ranges!', ktau, met, &
                     out%Qsb(land_x(i),land_y(i),1),ranges%Qsb,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Qsb(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Qsb,out%Qsb,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Qsb(:,1,1) = ssoil%rnof2/dels
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Qsb(i,1,1)<ranges%Qsb(1)).OR.(out%Qsb(i,1,1)>ranges%Qsb(2))) &
                     CALL range_abort('Qsb out of specified ranges!', ktau, met, &
                     out%Qsb(i,1,1),ranges%Qsb,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%Qsb,out%Qsb(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qsb variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Evap:
    IF(output%flux.OR.output%Evap) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Evap(land_x(i),land_y(i),1) = canopy%fe(i)/air%rlam(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Evap(land_x(i),land_y(i),1)<ranges%Evap(1)).OR. &
                     (out%Evap(land_x(i),land_y(i),1)>ranges%Evap(2))) &
                     CALL range_abort('Evap out of specified ranges!', ktau, met, &
                     out%Evap(land_x(i),land_y(i),1),ranges%Evap,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Evap(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Evap,out%Evap,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Evap(:,1,1) = canopy%fe/air%rlam
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Evap(i,1,1)<ranges%Evap(1)).OR.(out%Evap(i,1,1)>ranges%Evap(2))) &
                     CALL range_abort('Evap out of specified ranges!', ktau, met, &
                     out%Evap(i,1,1),ranges%Evap,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%Evap,out%Evap(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Evap variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! ECanop:
    IF(output%flux.OR.output%veg.OR.output%ECanop) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%ECanop(land_x(i),land_y(i),1) = canopy%fevw(i)/air%rlam(i) 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%ECanop(land_x(i),land_y(i),1)<ranges%ECanop(1)).OR. &
                     (out%ECanop(land_x(i),land_y(i),1)>ranges%ECanop(2))) &
                     CALL range_abort('ECanop out of specified ranges!', ktau, met, &
                     out%ECanop(land_x(i),land_y(i),1),ranges%ECanop,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%ECanop(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%ECanop,out%ECanop,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%ECanop(:,1,1) = canopy%fevw/air%rlam
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%ECanop(i,1,1)<ranges%ECanop(1)).OR.(out%ECanop(i,1,1)>ranges%ECanop(2))) &
                     CALL range_abort('ECanop out of specified ranges!', ktau, met, &
                     out%ECanop(i,1,1),ranges%ECanop,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%ECanop,out%ECanop(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing ECanop variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! TVeg:
    IF(output%flux.OR.output%veg.OR.output%TVeg) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%TVeg(land_x(i),land_y(i),1) = REAL(canopy%fevc(i)/air%rlam(i),r_1) 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%TVeg(land_x(i),land_y(i),1)<ranges%TVeg(1)).OR. &
                     (out%TVeg(land_x(i),land_y(i),1)>ranges%TVeg(2))) &
                     CALL range_abort('TVeg out of specified ranges!', ktau, met, &
                     out%TVeg(land_x(i),land_y(i),1),ranges%TVeg,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%TVeg(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%TVeg,out%TVeg,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%TVeg(:,1,1) = REAL(canopy%fevc/air%rlam,r_1)
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%TVeg(i,1,1)<ranges%TVeg(1)).OR.(out%TVeg(i,1,1)>ranges%TVeg(2))) &
                     CALL range_abort('TVeg out of specified ranges!', ktau, met, &
                     out%TVeg(i,1,1),ranges%TVeg,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%TVeg,out%TVeg(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing TVeg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! ESoil:
    IF(output%flux.OR.output%soil.OR.output%ESoil) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%ESoil(land_x(i),land_y(i),1) = canopy%fes(i)/air%rlam(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%ESoil(land_x(i),land_y(i),1)<ranges%ESoil(1)).OR. &
                     (out%ESoil(land_x(i),land_y(i),1)>ranges%ESoil(2))) &
                     CALL range_abort('ESoil out of specified ranges!', ktau, met, &
                     out%ESoil(land_x(i),land_y(i),1),ranges%ESoil,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%ESoil(:,:,1) = -0.01 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%ESoil,out%ESoil,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%ESoil(:,1,1) = canopy%fes/air%rlam
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%ESoil(i,1,1)<ranges%ESoil(1)).OR.(out%ESoil(i,1,1)>ranges%ESoil(2))) &
                     CALL range_abort('ESoil out of specified ranges!', ktau, met, &
                     out%ESoil(i,1,1),ranges%ESoil,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%ESoil,out%ESoil(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing ESoil variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! HVeg:
    IF(output%flux.OR.output%veg.OR.output%HVeg) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%HVeg(land_x(i),land_y(i),1) = canopy%fhv(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%HVeg(land_x(i),land_y(i),1)<ranges%HVeg(1)).OR. &
                     (out%HVeg(land_x(i),land_y(i),1)>ranges%HVeg(2))) &
                     CALL range_abort('HVeg out of specified ranges!', ktau, met, &
                     out%HVeg(land_x(i),land_y(i),1),ranges%HVeg,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%HVeg(:,:,1) = -1000.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%HVeg,out%HVeg,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%HVeg(:,1,1) = canopy%fhv
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%HVeg(i,1,1)<ranges%HVeg(1)).OR.(out%HVeg(i,1,1)>ranges%HVeg(2))) &
                     CALL range_abort('HVeg out of specified ranges!', ktau, met, &
                     out%HVeg(i,1,1),ranges%HVeg,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%HVeg,out%HVeg(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HVeg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! HSoil:
    IF(output%flux.OR.output%soil.OR.output%HSoil) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%HSoil(land_x(i),land_y(i),1) = canopy%fhs(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%HSoil(land_x(i),land_y(i),1)<ranges%HSoil(1)).OR. &
                     (out%HSoil(land_x(i),land_y(i),1)>ranges%HSoil(2))) &
                     CALL range_abort('HSoil out of specified ranges!', ktau, met, &
                     out%HSoil(land_x(i),land_y(i),1),ranges%HSoil,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%HSoil(:,:,1) = -1000.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%HSoil,out%HSoil,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%HSoil(:,1,1) = canopy%fhs
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%HSoil(i,1,1)<ranges%HSoil(1)).OR.(out%HSoil(i,1,1)>ranges%HSoil(2))) &
                     CALL range_abort('HSoil out of specified ranges!', ktau, met, &
                     out%HSoil(i,1,1),ranges%HSoil,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%HSoil,out%HSoil(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HSoil variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! NEE:
    IF(output%flux.OR.output%carbon.OR.output%NEE) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%NEE(land_x(i),land_y(i),1) = canopy%fnee(i)/1.201E-5 ! g/m2/s to umol/m2/s
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%NEE(land_x(i),land_y(i),1)<ranges%NEE(1)).OR. &
                     (out%NEE(land_x(i),land_y(i),1)>ranges%NEE(2))) &
                     CALL range_abort('NEE out of specified ranges!', ktau, met, &
                     out%NEE(land_x(i),land_y(i),1),ranges%NEE,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%NEE(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%NEE,out%NEE,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%NEE(:,1,1) = canopy%fnee/1.201E-5 ! g/m2/s to umol/m2/s 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%NEE(i,1,1)<ranges%NEE(1)).OR.(out%NEE(i,1,1)>ranges%NEE(2))) &
                     CALL range_abort('NEE out of specified ranges!', ktau, met, &
                     out%NEE(i,1,1),ranges%NEE,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%NEE,out%NEE(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing NEE variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    !--------------------------WRITE SOIL DATA-------------------------------------
    ! SoilMoist:
    IF(output%soil.OR.output%SoilMoist) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SoilMoist(land_x(i),land_y(i),:,1) = REAL(ssoil%wb(i,:)*soil%zse*1000.0,r_1)
             IF(check%ranges) THEN  ! Check ranges:
                DO j=1,ms
                   IF((out%SoilMoist(land_x(i),land_y(i),j,1)<ranges%SoilMoist(1)).OR. &
                        (out%SoilMoist(land_x(i),land_y(i),j,1)>ranges%SoilMoist(2))) &
                        CALL range_abort('SoilMoist out of specified ranges!', ktau, met, &
                        out%SoilMoist(land_x(i),land_y(i),j,1),ranges%SoilMoist,land_x(i),land_y(i))
                END DO
             END IF
          END DO
          WHERE(mask/=1) out%SoilMoist(:,:,1,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SoilMoist,out%SoilMoist,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,ms,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          DO i = 1, mp
             out%SoilMoist(i,:,1,1)=REAL(ssoil%wb(i,:)*soil%zse*1000.0,r_1)
             IF(check%ranges) THEN  ! Check ranges:
                DO j=1,ms
                   IF((out%SoilMoist(i,j,1,1)<ranges%SoilMoist(1)).OR. &
                        (out%SoilMoist(i,j,1,1)>ranges%SoilMoist(2))) &
                        CALL range_abort('SoilMoist out of specified ranges!', ktau, met, &
                        out%SoilMoist(i,j,1,1),ranges%SoilMoist,i)
                END DO
             END IF
          END DO
          status=NF90_PUT_VAR(ncid_out,ovid%SoilMoist,out%SoilMoist(:,:,:,1), &
               start=(/1,1,ktau/),count=(/mp,ms,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SoilMoist variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! SoilTemp:
    IF(output%soil.OR.output%SoilTemp) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SoilTemp(land_x(i),land_y(i),:,1) = ssoil%tgg(i,:)
             IF(check%ranges) THEN  ! Check ranges:
                DO j=1,ms
                   IF((out%SoilTemp(land_x(i),land_y(i),j,1)<ranges%SoilTemp(1)).OR. &
                        (out%SoilTemp(land_x(i),land_y(i),j,1)>ranges%SoilTemp(2))) &
                        CALL range_abort('SoilTemp out of specified ranges!', ktau, met, &
                        out%SoilTemp(land_x(i),land_y(i),j,1),ranges%SoilTemp,land_x(i),land_y(i))
                END DO
             END IF
          END DO
          WHERE(mask/=1) out%SoilTemp(:,:,1,1) = 200.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SoilTemp,out%SoilTemp,start=(/1,1,1,ktau/), &
               count=(/xdimsize,ydimsize,ms,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          DO i = 1, mp
             out%SoilTemp(i,:,1,1)=ssoil%tgg(i,:)
             IF(check%ranges) THEN  ! Check ranges:
                DO j=1,ms
                   IF((out%SoilTemp(i,j,1,1)<ranges%SoilTemp(1)).OR. &
                        (out%SoilTemp(i,j,1,1)>ranges%SoilTemp(2))) &
                        CALL range_abort('SoilTemp out of specified ranges!', ktau, met, &
                        out%SoilTemp(i,j,1,1),ranges%SoilTemp,i)
                END DO
             END IF
          END DO
          status=NF90_PUT_VAR(ncid_out,ovid%SoilTemp,out%SoilTemp(:,:,:,1), &
               start=(/1,1,ktau/),count=(/mp,ms,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SoilTemp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! BaresoilT:
    IF(output%soil.OR.output%BaresoilT) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%BaresoilT(land_x(i),land_y(i),1) = ssoil%tgg(i,1)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%BaresoilT(land_x(i),land_y(i),1)<ranges%BaresoilT(1)).OR. &
                     (out%BaresoilT(land_x(i),land_y(i),1)>ranges%BaresoilT(2))) &
                     CALL range_abort('BaresoilT out of specified ranges!', ktau, met, &
                     out%BaresoilT(land_x(i),land_y(i),1),ranges%BaresoilT,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%BaresoilT(:,:,1) = 200.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%BaresoilT,out%BaresoilT,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%BaresoilT(:,1,1) = ssoil%tgg(:,1) 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%BaresoilT(i,1,1)<ranges%BaresoilT(1)).OR.(out%BaresoilT(i,1,1)>ranges%BaresoilT(2))) &
                     CALL range_abort('BaresoilT out of specified ranges!', ktau, met, &
                     out%BaresoilT(i,1,1),ranges%BaresoilT,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%BaresoilT,out%BaresoilT(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing BaresoilT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! SWE:
    IF(output%snow.OR.output%SWE) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SWE(land_x(i),land_y(i),1) = ssoil%snowd(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%SWE(land_x(i),land_y(i),1)<ranges%SWE(1)).OR. &
                     (out%SWE(land_x(i),land_y(i),1)>ranges%SWE(2))) &
                     CALL range_abort('SWE out of specified ranges!', ktau, met, &
                     out%SWE(land_x(i),land_y(i),1),ranges%SWE,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%SWE(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SWE,out%SWE,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%SWE(:,1,1) = ssoil%snowd 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%SWE(i,1,1)<ranges%SWE(1)).OR.(out%SWE(i,1,1)>ranges%SWE(2))) &
                     CALL range_abort('SWE out of specified ranges!', ktau, met, &
                     out%SWE(i,1,1),ranges%SWE,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%SWE,out%SWE(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWE variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! SnowT:
    IF(output%snow.OR.output%SnowT) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SnowT(land_x(i),land_y(i),1) = ssoil%tggsn(i,1)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%SnowT(land_x(i),land_y(i),1)<ranges%SnowT(1)).OR. &
                     (out%SnowT(land_x(i),land_y(i),1)>ranges%SnowT(2))) &
                     CALL range_abort('SnowT out of specified ranges!', ktau, met, &
                     out%SnowT(land_x(i),land_y(i),1),ranges%SnowT,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%SnowT(:,:,1) = 300.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SnowT,out%SnowT,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%SnowT(:,1,1) = ssoil%tggsn(:,1) 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%SnowT(i,1,1)<ranges%SnowT(1)).OR.(out%SnowT(i,1,1)>ranges%SnowT(2))) &
                     CALL range_abort('SnowT out of specified ranges!', ktau, met, &
                     out%SnowT(i,1,1),ranges%SnowT,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%SnowT,out%SnowT(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SnowT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! SnowDepth:
    IF(output%snow.OR.output%SnowDepth) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SnowDepth(land_x(i),land_y(i),1) = SUM(ssoil%sdepth(i,:))
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%SnowDepth(land_x(i),land_y(i),1)<ranges%SnowDepth(1)).OR. &
                     (out%SnowDepth(land_x(i),land_y(i),1)>ranges%SnowDepth(2))) &
                     CALL range_abort('SnowDepth out of specified ranges!', ktau, met, &
                     out%SnowDepth(land_x(i),land_y(i),1),ranges%SnowDepth,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%SnowDepth(:,:,1) = -0.1 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SnowDepth,out%SnowDepth,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%SnowDepth(:,1,1) = SUM(ssoil%sdepth(:,:),2) 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%SnowDepth(i,1,1)<ranges%SnowDepth(1)).OR.(out%SnowDepth(i,1,1)>ranges%SnowDepth(2))) &
                     CALL range_abort('SnowDepth out of specified ranges!', ktau, met, &
                     out%SnowDepth(i,1,1),ranges%SnowDepth,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%SnowDepth,out%SnowDepth(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SnowDepth variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    !----------------------------WRITE RADIATION DATA----------------------------------
    ! SWnet:
    IF(output%radiation.OR.output%SWnet) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%SWnet(land_x(i),land_y(i),1) = SUM(rad%qcan(i,:,1)) + &
                  SUM(rad%qcan(i,:,2)) + rad%qssabs(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%SWnet(land_x(i),land_y(i),1)<ranges%SWnet(1)).OR. &
                     (out%SWnet(land_x(i),land_y(i),1)>ranges%SWnet(2))) &
                     CALL range_abort('SWnet out of specified ranges!', ktau, met, &
                     out%SWnet(land_x(i),land_y(i),1),ranges%SWnet,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%SWnet(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%SWnet,out%SWnet,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%SWnet(:,1,1) = SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%SWnet(i,1,1)<ranges%SWnet(1)).OR.(out%SWnet(i,1,1)>ranges%SWnet(2))) &
                     CALL range_abort('SWnet out of specified ranges!', ktau, met, &
                     out%SWnet(i,1,1),ranges%SWnet,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%SWnet,out%SWnet(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWnet variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! LWnet:
    IF(output%radiation.OR.output%LWnet) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%LWnet(land_x(i),land_y(i),1) = met%fld(i) - sboltz*emleaf* &
               canopy%tv(i)**4 * (1-rad%transd(i)) - rad%flws(i)*rad%transd(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%LWnet(land_x(i),land_y(i),1)<ranges%LWnet(1)).OR. &
                     (out%LWnet(land_x(i),land_y(i),1)>ranges%LWnet(2))) &
                     CALL range_abort('LWnet out of specified ranges!', ktau, met, &
                     out%LWnet(land_x(i),land_y(i),1),ranges%LWnet,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%LWnet(:,:,1) = 100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%LWnet,out%LWnet,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%LWnet(:,1,1) = met%fld - sboltz*emleaf* &
               canopy%tv**4 * (1-rad%transd) - rad%flws*rad%transd
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%LWnet(i,1,1)<ranges%LWnet(1)).OR.(out%LWnet(i,1,1)>ranges%LWnet(2))) &
                     CALL range_abort('LWnet out of specified ranges!', ktau, met, &
                     out%LWnet(i,1,1),ranges%LWnet,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%LWnet,out%LWnet(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LWnet variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Rnet:
    IF(output%radiation.OR.output%Rnet) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Rnet(land_x(i),land_y(i),1) = & ! LWnet:
                  met%fld(i)-sboltz*emleaf*canopy%tv(i)**4 * &
             (1-rad%transd(i)) - rad%flws(i)*rad%transd(i) & 
                  + SUM(rad%qcan(i,:,1)) + & ! + SWnet
                  SUM(rad%qcan(i,:,2)) + rad%qssabs(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Rnet(land_x(i),land_y(i),1)<ranges%Rnet(1)).OR. &
                     (out%Rnet(land_x(i),land_y(i),1)>ranges%Rnet(2))) &
                     CALL range_abort('Rnet out of specified ranges!', ktau, met, &
                     out%Rnet(land_x(i),land_y(i),1),ranges%Rnet,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Rnet(:,:,1) = 100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Rnet,out%Rnet,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Rnet(:,1,1) = & ! LWnet + SWnet
               met%fld-sboltz*emleaf*canopy%tv**4 *(1-rad%transd)-rad%flws*rad%transd &
               + SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Rnet(i,1,1)<ranges%Rnet(1)).OR.(out%Rnet(i,1,1)>ranges%Rnet(2))) &
                     CALL range_abort('Rnet out of specified ranges!', ktau, met, &
                     out%Rnet(i,1,1),ranges%Rnet,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%Rnet,out%Rnet(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Rnet variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Albedo:
    IF(output%radiation.OR.output%Albedo) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Albedo(land_x(i),land_y(i),1) = (rad%albedo(i,1) + rad%albedo(i,2))*0.5
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%Albedo(land_x(i),land_y(i),1)<ranges%Albedo(1)).OR. &
                     (out%Albedo(land_x(i),land_y(i),1)>ranges%Albedo(2))) &
                     CALL range_abort('Albedo out of specified ranges!', ktau, met, &
                     out%Albedo(land_x(i),land_y(i),1),ranges%Albedo,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%Albedo(:,:,1) = -0.1 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Albedo,out%Albedo,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Albedo(:,1,1) = (rad%albedo(:,1) + rad%albedo(:,2))*0.5
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%Albedo(i,1,1)<ranges%Albedo(1)).OR.(out%Albedo(i,1,1)>ranges%Albedo(2))) &
                     CALL range_abort('Albedo out of specified ranges!', ktau, met, &
                     out%Albedo(i,1,1),ranges%Albedo,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%Albedo,out%Albedo(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Albedo variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! RadT:
    IF(output%radiation.OR.output%RadT) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%RadT(land_x(i),land_y(i),1) = (1.0-rad%transd(i))*emleaf*sboltz* &
               canopy%tv(i)**4 + rad%transd(i)*emsoil*sboltz* ( &
               (1-ssoil%isflag(i))*ssoil%tgg(i,1)+ssoil%isflag(i)*ssoil%tggsn(i,1)) ** 4
             out%RadT(land_x(i),land_y(i),1) = (out%RadT(land_x(i),land_y(i),1)/sboltz)**0.25
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%RadT(land_x(i),land_y(i),1)<ranges%RadT(1)).OR. &
                     (out%RadT(land_x(i),land_y(i),1)>ranges%RadT(2))) &
                     CALL range_abort('RadT out of specified ranges!', ktau, met, &
                     out%RadT(land_x(i),land_y(i),1),ranges%RadT,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%RadT(:,:,1) = 200.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%RadT,out%RadT,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%RadT(:,1,1) = (1.0-rad%transd)*emleaf*sboltz* &
               canopy%tv**4 + rad%transd*emsoil*sboltz* ( &
               (1-ssoil%isflag)*ssoil%tgg(:,1)+ssoil%isflag*ssoil%tggsn(:,1)) ** 4
          out%RadT(:,1,1) = (out%RadT(:,1,1)/sboltz)**0.25
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%RadT(i,1,1)<ranges%RadT(1)).OR.(out%RadT(i,1,1)>ranges%RadT(2))) &
                     CALL range_abort('RadT out of specified ranges!', ktau, met, &
                     out%RadT(i,1,1),ranges%RadT,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%RadT,out%RadT(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing RadT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! VegT:
    IF(output%veg.OR.output%VegT) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%VegT(land_x(i),land_y(i),1) = canopy%tv(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%VegT(land_x(i),land_y(i),1)<ranges%VegT(1)).OR. &
                     (out%VegT(land_x(i),land_y(i),1)>ranges%VegT(2))) &
                     CALL range_abort('VegT out of specified ranges!', ktau, met, &
                     out%VegT(land_x(i),land_y(i),1),ranges%VegT,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%VegT(:,:,1) = 200.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%VegT,out%VegT,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%VegT(:,1,1) = canopy%tv
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%VegT(i,1,1)<ranges%VegT(1)).OR.(out%VegT(i,1,1)>ranges%VegT(2))) &
                     CALL range_abort('VegT out of specified ranges!', ktau, met, &
                     out%VegT(i,1,1),ranges%VegT,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%VegT,out%VegT(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing VegT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! CanopInt:
    IF(output%veg.OR.output%CanopInt) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%CanopInt(land_x(i),land_y(i),1) = canopy%cansto(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%CanopInt(land_x(i),land_y(i),1)<ranges%CanopInt(1)).OR. &
                     (out%CanopInt(land_x(i),land_y(i),1)>ranges%CanopInt(2))) &
                     CALL range_abort('CanopInt out of specified ranges!', ktau, met, &
                     out%CanopInt(land_x(i),land_y(i),1),ranges%CanopInt,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%CanopInt(:,:,1) = -0.1 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%CanopInt,out%CanopInt,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%CanopInt(:,1,1) = canopy%cansto
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%CanopInt(i,1,1)<ranges%CanopInt(1)).OR.(out%CanopInt(i,1,1)>ranges%CanopInt(2))) &
                     CALL range_abort('CanopInt out of specified ranges!', ktau, met, &
                     out%CanopInt(i,1,1),ranges%CanopInt,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%CanopInt,out%CanopInt(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing CanopInt variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! LAI:
    IF(output%veg.OR.output%LAI) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%LAI(land_x(i),land_y(i),1) = veg%vlai(i)
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%LAI(land_x(i),land_y(i),1)<ranges%LAI(1)).OR. &
                     (out%LAI(land_x(i),land_y(i),1)>ranges%LAI(2))) &
                     CALL range_abort('LAI out of specified ranges!', ktau, met, &
                     out%LAI(land_x(i),land_y(i),1),ranges%LAI,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%LAI(:,:,1) = -1.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%LAI,out%LAI,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%LAI(:,1,1) = veg%vlai
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%LAI(i,1,1)<ranges%LAI(1)).OR.(out%LAI(i,1,1)>ranges%LAI(2))) &
                     CALL range_abort('LAI out of specified ranges!', ktau, met, &
                     out%LAI(i,1,1),ranges%LAI,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%LAI,out%LAI(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LAI variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Ebal:
    IF(output%balances.OR.output%Ebal) THEN ! If Ebal is requested output
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Ebal(land_x(i),land_y(i),1) = bal%ebal_tot(i)
          END DO
          WHERE(mask/=1) out%Ebal(:,:,1) = 0.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Ebal,out%Ebal,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Ebal(:,1,1) = bal%ebal_tot
          status=NF90_PUT_VAR(ncid_out,ovid%Ebal,out%Ebal(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       ! Check writing was successful:
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Ebal variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! Wbal:
    IF(output%balances.OR.output%Wbal) THEN ! If Wbal is requested output
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%Wbal(land_x(i),land_y(i),1) = bal%wbal_tot(i)
          END DO
          WHERE(mask/=1) out%Wbal(:,:,1) = 0.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%Wbal,out%Wbal,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%Wbal(:,1,1) = bal%wbal_tot
          status=NF90_PUT_VAR(ncid_out,ovid%Wbal,out%Wbal(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       ! Check writing was successful:
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Wbal variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! GPP:
    IF(output%carbon.OR.output%GPP) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%GPP(land_x(i),land_y(i),1) = -1.0*canopy%fpn(i)/1.201E-5
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%GPP(land_x(i),land_y(i),1)<ranges%GPP(1)).OR. &
                     (out%GPP(land_x(i),land_y(i),1)>ranges%GPP(2))) &
                     CALL range_abort('GPP out of specified ranges!', ktau, met, &
                     out%GPP(land_x(i),land_y(i),1),ranges%GPP,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%GPP(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%GPP,out%GPP,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%GPP(:,1,1) = -1.0*canopy%fpn/1.201E-5
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%GPP(i,1,1)<ranges%GPP(1)).OR.(out%GPP(i,1,1)>ranges%GPP(2))) &
                     CALL range_abort('GPP out of specified ranges!', ktau, met, &
                     out%GPP(i,1,1),ranges%GPP,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%GPP,out%GPP(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing GPP variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! AutoResp:
    IF(output%carbon.OR.output%AutoResp) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%AutoResp(land_x(i),land_y(i),1) = canopy%frp(i)/1.201E-5 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%AutoResp(land_x(i),land_y(i),1)<ranges%AutoResp(1)).OR. &
                     (out%AutoResp(land_x(i),land_y(i),1)>ranges%AutoResp(2))) &
                     CALL range_abort('AutoResp out of specified ranges!', ktau, met, &
                     out%AutoResp(land_x(i),land_y(i),1),ranges%AutoResp,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%AutoResp(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%AutoResp,out%AutoResp,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%AutoResp(:,1,1) = canopy%frp/1.201E-5 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%AutoResp(i,1,1)<ranges%AutoResp(1)).OR.(out%AutoResp(i,1,1)>ranges%AutoResp(2))) &
                     CALL range_abort('AutoResp out of specified ranges!', ktau, met, &
                     out%AutoResp(i,1,1),ranges%AutoResp,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%AutoResp,out%AutoResp(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing AutoResp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! HeteroResp:
    IF(output%carbon.OR.output%HeteroResp) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%HeteroResp(land_x(i),land_y(i),1) = canopy%frs(i)/1.201E-5 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%HeteroResp(land_x(i),land_y(i),1)<ranges%HeteroResp(1)).OR. &
                     (out%HeteroResp(land_x(i),land_y(i),1)>ranges%HeteroResp(2))) &
                     CALL range_abort('HeteroResp out of specified ranges!', ktau, met, &
                     out%HeteroResp(land_x(i),land_y(i),1),ranges%HeteroResp,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%HeteroResp(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%HeteroResp,out%HeteroResp,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%HeteroResp(:,1,1) = canopy%frs/1.201E-5 
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%HeteroResp(i,1,1)<ranges%HeteroResp(1)).OR. &
                     (out%HeteroResp(i,1,1)>ranges%HeteroResp(2))) &
                     CALL range_abort('HeteroResp out of specified ranges!', ktau, met, &
                     out%HeteroResp(i,1,1),ranges%HeteroResp,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%HeteroResp,out%HeteroResp(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HeteroResp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    ! NPP:
    IF(output%carbon.OR.output%NPP) THEN
       IF(gridType=='mask') THEN ! if land sea mask type grid
          DO i = 1, mp
             out%NPP(land_x(i),land_y(i),1) = out%GPP(land_x(i),land_y(i),1) - &
               out%AutoResp(land_x(i),land_y(i),1) 
             IF(check%ranges) THEN  ! Check ranges:
                IF((out%NPP(land_x(i),land_y(i),1)<ranges%NPP(1)).OR. &
                     (out%NPP(land_x(i),land_y(i),1)>ranges%NPP(2))) &
                     CALL range_abort('NPP out of specified ranges!', ktau, met, &
                     out%NPP(land_x(i),land_y(i),1),ranges%NPP,land_x(i),land_y(i))
             END IF
          END DO
          WHERE(mask/=1) out%NPP(:,:,1) = -100.0 ! not land
          status=NF90_PUT_VAR(ncid_out,ovid%NPP,out%NPP,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/)) ! write to file
       ELSE IF(gridType=='land') THEN ! else if land only grid
          out%NPP(:,1,1) = out%GPP(:,1,1) - out%AutoResp(:,1,1)
          IF(check%ranges) THEN  ! Check ranges:
             DO i = 1, mp
                IF((out%NPP(i,1,1)<ranges%NPP(1)).OR. &
                     (out%NPP(i,1,1)>ranges%NPP(2))) &
                     CALL range_abort('NPP out of specified ranges!', ktau, met, &
                     out%NPP(i,1,1),ranges%NPP,i)
             END DO
          END IF
          status=NF90_PUT_VAR(ncid_out,ovid%NPP,out%NPP(:,:,1), &
               start=(/1,ktau/),count=(/mp,1/))
       END IF
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing NPP variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
  END SUBROUTINE write_output
  !====================================================================
  SUBROUTINE close_output_file(filename_out,bal,air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)
    CHARACTER(LEN=*), INTENT(IN) :: filename_out
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (radiation_type),INTENT(INOUT)  :: rad
    TYPE (sum_flux_type), INTENT(INOUT)  :: sum_flux
    TYPE(balances_type),INTENT(INOUT) :: bal
    INTEGER :: i ! do loop counter

    ! Close file
    status = NF90_CLOSE(ncid_out)
    IF(status/=NF90_NOERR) CALL nc_abort('Error closing output file ' &
         //TRIM(filename_out)// '(SUBROUTINE close_output_file)')

    ! Deallocate output variables:
    IF(output%met) DEALLOCATE(out%SWdown,out%LWdown,out%Tair, &
         out%PSurf,out%Qair,out%CO2air,out%Rainf,out%Snowf,out%Wind) 
    IF(output%flux) DEALLOCATE(out%Qle,out%Qh,out%Qg,out%Qs,out%Qsb, &
         out%Evap,out%ECanop,out%TVeg,out%ESoil,out%HVeg,out%HSoil, &
         out%NEE) 
    IF(output%soil) DEALLOCATE(out%SoilMoist,out%SoilTemp,out%BaresoilT,&
         out%SWE,out%SnowT,out%SnowDepth) 
    IF(output%radiation) DEALLOCATE(out%SWnet,out%LWnet,out%Rnet,out%Albedo,out%RadT) 
    IF(output%veg) DEALLOCATE(out%VegT,out%LAI,out%CanopInt) 
    IF(output%carbon) THEN
       DEALLOCATE(out%NPP,out%GPP,out%AutoResp,out%HeteroResp)
       IF(.NOT.output%flux) DEALLOCATE(out%NEE)
    END IF
    IF(output%balances) THEN
       DEALLOCATE(out%Ebal,out%Wbal)
       IF(verbose) THEN
          WRITE(logn,*)
          DO i = 1, mp
             WRITE(logn,'(A33,I7,1X,A2,E12.4,A6)') ' Cumulative energy balance site #',&
                  i,'is',bal%ebal_tot(i),' W/m^2'
             WRITE(logn,'(A32,I7,1X,A2,E12.4,A3)') ' Cumulative water balance site #', &
                  i,'is',bal%wbal_tot(i),' mm'
          END DO
       END IF
    END IF
    IF(output%params) THEN
       DEALLOCATE(out%isoil,out%iveg)
    END IF
    DEALLOCATE(mask,timevar,land_x,land_y)
    DEALLOCATE(latitude,longitude,gdpt)

    ! Deallocate CABLE main variables:
    CALL dealloc_cbm_var(air, mp)
    CALL dealloc_cbm_var(bgc, mp)
    CALL dealloc_cbm_var(canopy, mp)
    CALL dealloc_cbm_var(met, mp)
    CALL dealloc_cbm_var(bal, mp)
    CALL dealloc_cbm_var(rad, mp)
    CALL dealloc_cbm_var(rough, mp)
    CALL dealloc_cbm_var(soil, mp)
    CALL dealloc_cbm_var(ssoil, mp)
    CALL dealloc_cbm_var(sum_flux, mp)
    CALL dealloc_cbm_var(veg, mp)

    ! Successful run!
    WRITE(logn,*)
    WRITE(logn,*) 'Run finished and output file closed.'
  END SUBROUTINE close_output_file
  !==========================================================================
  SUBROUTINE create_restart(filename_restart_out,filename_met,logn,kstart, &
       kend,soil,veg,ssoil,canopy,rough,bgc,bal)
    ! Creates a restart file for CABLE using a land only grid of size mp
    ! and CABLE's internal variable names.
    CHARACTER(LEN=*), INTENT(IN) :: filename_met 
    CHARACTER(LEN=*), INTENT(IN) :: filename_restart_out 
    INTEGER, INTENT(IN) :: logn ! log file number
    INTEGER, INTENT(IN) :: kstart ! starting time step number in run
    INTEGER, INTENT(IN) :: kend ! ending time step number in run
    TYPE (soil_parameter_type),INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type),INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (soil_snow_type),INTENT(IN)      :: ssoil  ! soil and snow variables
    TYPE (bgc_pool_type),INTENT(IN)       :: bgc    ! carbon pool variables
    TYPE (canopy_type),INTENT(IN)         :: canopy ! vegetation variables
    TYPE (roughness_type),INTENT(IN)      :: rough  ! roughness varibles
    TYPE (balances_type),INTENT(IN)  :: bal  ! energy and water balance variables
    TYPE(parID_type) :: rpid ! parameter IDs for restart nc file
    INTEGER :: ncid_restart ! netcdf restart file ID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
    INTEGER :: mpID,radID,soilID,soilcarbID,plantcarbID,tID,snowID ! dimension IDs
    INTEGER :: tvarID,latID,lonID ! time,lat,lon variable ID
    INTEGER :: tggID,wbID,wbiceID,tssID,ssdnnID,ssdnID,osnowdID, &
         smassID,sdepthID,snageID,snowdID,rtsoilID,isflagID,canstoID, &
         albsoilsnID,gammzzID,tggsnID,sghfluxID,ghfluxID,runoffID, &
         rnof1ID,rnof2ID,gaID,dgdtgID,fevID,fesID,fhsID,wbtot0ID,osnowd0ID, &
         cplantID,csoilID

    WRITE(logn,'(A24)') ' Writing restart file...'
    ! Create output file:
    status = NF90_CREATE(filename_restart_out,NF90_CLOBBER,ncid_restart)
    IF(status/=NF90_NOERR) CALL nc_abort('Error creating restart file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! Put the file in define mode:
    status = NF90_REDEF(ncid_restart)
    ! Define dimensions:
    status = NF90_DEF_DIM(ncid_restart,'mp',mp,mpID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining mp dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'soil',ms,soilID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining vertical soil dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'snow',3,snowID) ! number of snow layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining vertical snow dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'rad',nrb,radID) ! number of radiation bands
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining radiation dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'soil_carbon_pools',ncs,soilcarbID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining soil carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'plant_carbon_pools',ncp,plantcarbID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining plant carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_DEF_DIM(ncid_restart,'time',1,tID) 
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define "time" variable and its attributes:
    status=NF90_DEF_VAR(ncid_restart,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tvarID,'units',timeunits)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tvarID,'coordinate',time_coord)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define latitude and longitude variable:
    status=NF90_DEF_VAR(ncid_restart,'latitude',NF90_FLOAT,(/mpID/),latID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining latitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,latID,'units','degrees_north')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining latitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status=NF90_DEF_VAR(ncid_restart,'longitude',NF90_FLOAT,(/mpID/),lonID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining longitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,lonID,'units','degrees_east')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining longitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !=======begin defining state variables========================================
    !------------------soil states--------------------------------------
    ! Define tgg and units:
    status=NF90_DEF_VAR(ncid_restart,'tgg',NF90_FLOAT,(/mpID,soilID/),tggID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tgg variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tggID,'units','K')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tgg variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tggID,'long_name', &
         'Average layer soil temperature')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tgg variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define wb and units:
    status=NF90_DEF_VAR(ncid_restart,'wb',NF90_DOUBLE,(/mpID,soilID/),wbID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wb variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbID,'units','vol/vol')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wb variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbID,'long_name', &
         'Average layer volumetric soil moisture')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wb variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define wbice and units:
    status=NF90_DEF_VAR(ncid_restart,'wbice',NF90_DOUBLE,(/mpID,soilID/),wbiceID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbice variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbiceID,'units','vol/vol')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbice variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbiceID,'long_name', &
         'Average layer volumetric soil ice')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbice variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define tss and units:
    status=NF90_DEF_VAR(ncid_restart,'tss',NF90_FLOAT,(/mpID/),tssID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tss variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tssID,'units','K')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tss variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tssID,'long_name', &
         'Combined soil/snow temperature')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tss variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define albsoilsn and units:
    status=NF90_DEF_VAR(ncid_restart,'albsoilsn',NF90_FLOAT,(/mpID,radID/),albsoilsnID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining albsoilsn variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,albsoilsnID,'units','-')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining albsoilsn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,albsoilsnID,'long_name', &
         'Combined soil/snow albedo')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining albsoilsn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define rtsoil and units:
    status=NF90_DEF_VAR(ncid_restart,'rtsoil',NF90_FLOAT,(/mpID/),rtsoilID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rtsoil variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rtsoilID,'units','??')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rtsoil variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rtsoilID,'long_name', &
         'turbulent resistance for soil')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rtsoil variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define gammzz and units:
    status=NF90_DEF_VAR(ncid_restart,'gammzz',NF90_DOUBLE,(/mpID,soilID/),gammzzID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining gammzz variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,gammzzID,'units',' ')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining gammzz variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,gammzzID,'long_name', &
         'heat capacity for each soil layer')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining gammzz variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')  
    ! Define runoff and units:
    status=NF90_DEF_VAR(ncid_restart,'runoff',NF90_FLOAT,(/mpID/),runoffID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining runoff variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,runoffID,'units','mm/timestep')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining runoff variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,runoffID,'long_name', &
         'Total runoff')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining runoff variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define rnof1 and units:
    status=NF90_DEF_VAR(ncid_restart,'rnof1',NF90_FLOAT,(/mpID/),rnof1ID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof1 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rnof1ID,'units','mm/timestep')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof1 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rnof1ID,'long_name', &
         'Surface runoff')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof1 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define rnof2 and units:
    status=NF90_DEF_VAR(ncid_restart,'rnof2',NF90_FLOAT,(/mpID/),rnof2ID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof2 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rnof2ID,'units','mm/timestep')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof2 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rnof2ID,'long_name', &
         'Subsurface runoff')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rnof2 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !-----------------snow states--------------------------------------------
    ! Define tggsn and units:
    status=NF90_DEF_VAR(ncid_restart,'tggsn',NF90_FLOAT,(/mpID,snowID/),tggsnID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tggsn variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tggsnID,'units','K')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tggsn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,tggsnID,'long_name', &
         'Average layer snow temperature')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tggsn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define ssdnn and units:
    status=NF90_DEF_VAR(ncid_restart,'ssdnn',NF90_FLOAT,(/mpID/),ssdnnID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdnn variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ssdnnID,'units','kg/m^3')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdnn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ssdnnID,'long_name', &
         'Average snow density')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdnn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define ssdn and units:
    status=NF90_DEF_VAR(ncid_restart,'ssdn',NF90_FLOAT,(/mpID,snowID/),ssdnID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdn variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ssdnID,'units','kg/m^3')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ssdnID,'long_name', &
         'Average layer snow density')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssdn variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define snowd and units:
    status=NF90_DEF_VAR(ncid_restart,'snowd',NF90_FLOAT,(/mpID/),snowdID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snowd variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,snowdID,'units','mm?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snowd variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,snowdID,'long_name', &
         'Liquid water eqivalent snow depth')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snowd variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define snage and units:
    status=NF90_DEF_VAR(ncid_restart,'snage',NF90_FLOAT,(/mpID/),snageID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snage variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,snageID,'units','?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snage variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,snageID,'long_name', &
         'Snow age')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining snage variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define smass and units:
    status=NF90_DEF_VAR(ncid_restart,'smass',NF90_FLOAT,(/mpID,snowID/),smassID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining smass variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,smassID,'units','kg?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining smass variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,smassID,'long_name', &
         'Average layer snow mass')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining smass variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define sdepth and units:
    status=NF90_DEF_VAR(ncid_restart,'sdepth',NF90_FLOAT,(/mpID,snowID/),sdepthID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sdepth variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,sdepthID,'units','mm?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sdepth variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,sdepthID,'long_name', &
         'Snow layer depth')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sdepth variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define osnowd and units:
    status=NF90_DEF_VAR(ncid_restart,'osnowd',NF90_FLOAT,(/mpID/),osnowdID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,osnowdID,'units','mm?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,osnowdID,'long_name', &
         'Previous time step snow depth')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define isflag and units:
    status=NF90_DEF_VAR(ncid_restart,'isflag',NF90_INT,(/mpID/),isflagID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining isflag variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,isflagID,'units','mm?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining isflag variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,isflagID,'long_name', &
         'Snow flag')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining isflag variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !----------------define canopy states----------------------------------
    ! Define cansto and units:
    status=NF90_DEF_VAR(ncid_restart,'cansto',NF90_FLOAT,(/mpID/),canstoID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cansto variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,canstoID,'units','mm')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cansto variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,canstoID,'long_name', &
         'Canopy water storage')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cansto variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define ghflux and units:
    status=NF90_DEF_VAR(ncid_restart,'ghflux',NF90_FLOAT,(/mpID/),ghfluxID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ghflux variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ghfluxID,'units','W/m^2?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ghflux variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,ghfluxID,'long_name', &
         '????')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ghflux variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define sghflux and units:
    status=NF90_DEF_VAR(ncid_restart,'sghflux',NF90_FLOAT,(/mpID/),sghfluxID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sghflux variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,sghfluxID,'units','W/m^2?')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sghflux variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,sghfluxID,'long_name', &
         '????')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sghflux variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define ga and units:
    status=NF90_DEF_VAR(ncid_restart,'ga',NF90_FLOAT,(/mpID/),gaID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ga variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,gaID,'units','W/m^2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ga variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,gaID,'long_name', &
         'Ground heat flux')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ga variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define dgdtg and units:
    status=NF90_DEF_VAR(ncid_restart,'dgdtg',NF90_FLOAT,(/mpID/),dgdtgID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining dgdtg variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,dgdtgID,'units','W/m^2/K')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining dgdtg variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,dgdtgID,'long_name', &
         'Derivative of ground heat flux wrt soil temperature')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining dgdtg variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define fev and units:
    status=NF90_DEF_VAR(ncid_restart,'fev',NF90_FLOAT,(/mpID/),fevID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fev variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fevID,'units','W/m^2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fev variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fevID,'long_name', &
         'Latent heat flux from vegetation')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fev variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define fes and units:
    status=NF90_DEF_VAR(ncid_restart,'fes',NF90_FLOAT,(/mpID/),fesID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fes variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fesID,'units','W/m^2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fes variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fesID,'long_name', &
         'Latent heat flux from soil')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fes variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define fhs and units:
    status=NF90_DEF_VAR(ncid_restart,'fhs',NF90_FLOAT,(/mpID/),fhsID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fhs variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fhsID,'units','W/m^2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fhs variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,fhsID,'long_name', &
         'Sensible heat flux from soil')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining fhs variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !--------------biogeochemical variables------------------------
    ! Define cplant and units:
    status=NF90_DEF_VAR(ncid_restart,'cplant',NF90_FLOAT,(/mpID,plantcarbID/), &
         cplantID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cplant variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,cplantID,'units','gC/m2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cplant variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,cplantID,'long_name', &
         'Plant carbon stores')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining cplant variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define csoil and units:
    status=NF90_DEF_VAR(ncid_restart,'csoil',NF90_FLOAT,(/mpID,soilcarbID/), &
         csoilID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining csoil variable in restart file. (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,csoilID,'units','gC/m2')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining csoil variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,csoilID,'long_name', &
         'Plant carbon stores')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining csoil variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !-------------------others---------------------------------
    ! Define wbtot0 and units:
    status=NF90_DEF_VAR(ncid_restart,'wbtot0',NF90_FLOAT,(/mpID/),wbtot0ID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbtot0 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbtot0ID,'units','mm')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbtot0 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,wbtot0ID,'long_name', &
         'Initial time step soil water total')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining wbtot0 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define osnowd0 and units:
    status=NF90_DEF_VAR(ncid_restart,'osnowd0',NF90_FLOAT,(/mpID/),osnowd0ID)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd0 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,osnowd0ID,'units','mm')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd0 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,osnowd0ID,'long_name', &
         'Initial time step snow water total')
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining osnowd0 variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    !---------------------MODEL PARAMETERS---------------------------------   
    WRITE(logn,*) 'Writing model parameters to restart file'
    ! Allocate size of output parameter variables:
    ! iveg (Vegetation type):
    status=NF90_DEF_VAR(ncid_restart,'iveg',NF90_FLOAT,(/mpID/), &
         rpid%iveg)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining iveg variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%iveg,"long_name",&
         "Vegetation type")
    status = NF90_PUT_ATT(ncid_restart,rpid%iveg,"units","-")
    ! isoil (Soil type):
    status=NF90_DEF_VAR(ncid_restart,'isoil',NF90_FLOAT,(/mpID/), &
         rpid%isoil)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining isoil variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%isoil,"long_name",&
         "Soil type")
    status = NF90_PUT_ATT(ncid_restart,rpid%isoil,"units","-")
    ! clay (fraction of soil which is clay):
    status=NF90_DEF_VAR(ncid_restart,'clay',NF90_FLOAT,(/mpID/), &
         rpid%clay)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining clay variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%clay,"long_name",&
         "Fraction of soil which is clay")
    status = NF90_PUT_ATT(ncid_restart,rpid%clay,"units","-")
    ! sand:
    status=NF90_DEF_VAR(ncid_restart,'sand',NF90_FLOAT,(/mpID/), &
         rpid%sand)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sand variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%sand,"long_name",&
         "Fraction of soil which is sand")
    status = NF90_PUT_ATT(ncid_restart,rpid%sand,"units","-")
    ! silt :
    status=NF90_DEF_VAR(ncid_restart,'silt',NF90_FLOAT,(/mpID/), &
         rpid%silt)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining silt variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%silt,"long_name",&
         "Fraction of soil which is silt")
    status = NF90_PUT_ATT(ncid_restart,rpid%silt,"units","-")
    ! ssat (Volume of soil volume which is water @ saturation):
    status=NF90_DEF_VAR(ncid_restart,'ssat',NF90_FLOAT,(/mpID/), &
         rpid%ssat)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ssat variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%ssat,"long_name",&
         "Fraction of soil volume which is water @ saturation")
    status = NF90_PUT_ATT(ncid_restart,rpid%ssat,"units","-")
    ! sfc (Volume of soil volume which is water @ field capacity):
    status=NF90_DEF_VAR(ncid_restart,'sfc',NF90_FLOAT,(/mpID/), &
         rpid%sfc)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sfc variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%sfc,"long_name",&
         "Fraction of soil volume which is water @ field capacity")
    status = NF90_PUT_ATT(ncid_restart,rpid%sfc,"units","-")
    ! swilt (Volume of soil volume which is water @ wilting point):
    status=NF90_DEF_VAR(ncid_restart,'swilt',NF90_FLOAT,(/mpID/), &
         rpid%swilt)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining swilt variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%swilt,"long_name",&
         "Fraction of soil volume which is water @ wilting point")
    status = NF90_PUT_ATT(ncid_restart,rpid%swilt,"units","-")
    ! zse (depth of each soil layer):
    status=NF90_DEF_VAR(ncid_restart,'zse',NF90_FLOAT,(/soilID/), &
         rpid%zse)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining zse variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%zse,"long_name",&
         "Depth of each soil layer")
    status = NF90_PUT_ATT(ncid_restart,rpid%zse,"units","m")
    ! froot (fraction of roots in each soil layer):
    status=NF90_DEF_VAR(ncid_restart,'froot',NF90_FLOAT,(/mpID,soilID/), &
         rpid%froot)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining froot variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%froot,"long_name",&
         "Fraction of roots in each soil layer")
    status = NF90_PUT_ATT(ncid_restart,rpid%froot,"units","-")
    ! bch (Parameter b, Campbell eqn 1985):
    status=NF90_DEF_VAR(ncid_restart,'bch',NF90_FLOAT,(/mpID/), &
         rpid%bch)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining bch variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%bch,"long_name","Parameter b, Campbell eqn 1985")
    status = NF90_PUT_ATT(ncid_restart,rpid%bch,"units","-")
    ! hyds (Hydraulic conductivity @ saturation):
    status=NF90_DEF_VAR(ncid_restart,'hyds',NF90_FLOAT,(/mpID/), &
         rpid%hyds)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining hyds variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%hyds,"long_name",&
         "Hydraulic conductivity @ saturation")
    status = NF90_PUT_ATT(ncid_restart,rpid%hyds,"units","m/s")
    ! sucs (suction @ saturation):
    status=NF90_DEF_VAR(ncid_restart,'sucs',NF90_FLOAT,(/mpID/), &
         rpid%sucs)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining sucs variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%sucs,"long_name",&
         "Suction @ saturation")
    status = NF90_PUT_ATT(ncid_restart,rpid%sucs,"units","m")
    ! css (Heat capacity of soil minerals):
    status=NF90_DEF_VAR(ncid_restart,'css',NF90_FLOAT,(/mpID/), &
         rpid%css)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining css variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%css,"long_name",&
         "Heat capacity of soil minerals")
    status = NF90_PUT_ATT(ncid_restart,rpid%css,"units","J/kg/C")
    ! rhosoil (density of soil minerals):
    status=NF90_DEF_VAR(ncid_restart,'rhosoil',NF90_FLOAT,(/mpID/), &
         rpid%rhosoil)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rhosoil variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%rhosoil,"long_name",&
         "Density of soil minerals")
    status = NF90_PUT_ATT(ncid_restart,rpid%rhosoil,"units","kg/m^3")
    ! rs20 (Soil respiration coefficient at 20C):
    status=NF90_DEF_VAR(ncid_restart,'rs20',NF90_FLOAT,(/mpID/), &
         rpid%rs20)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rs20 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%rs20,"long_name",&
         "Soil respiration coefficient at 20C")
    status = NF90_PUT_ATT(ncid_restart,rpid%rs20,"units","-")
    ! albsoil (Soil reflectance):
    status=NF90_DEF_VAR(ncid_restart,'albsoil',NF90_FLOAT,(/mpID/), &
         rpid%albsoil)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining albsoil variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%albsoil,"long_name",&
         "Soil reflectance")
    status = NF90_PUT_ATT(ncid_restart,rpid%albsoil,"units","-")
    ! hc (height of canopy):
    status=NF90_DEF_VAR(ncid_restart,'hc',NF90_FLOAT,(/mpID/), &
         rpid%hc)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining hc variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%hc,"long_name",&
         "Height of canopy")
    status = NF90_PUT_ATT(ncid_restart,rpid%hc,"units","mm/LAI")
    ! canst1 (Max water intercepted by canopy):
    status=NF90_DEF_VAR(ncid_restart,'canst1',NF90_FLOAT,(/mpID/), &
         rpid%canst1)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining canst1 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%canst1,"long_name",&
         "Max water intercepted by canopy")
    status = NF90_PUT_ATT(ncid_restart,rpid%canst1,"units","mm/LAI")
    ! dleaf (Chararacteristic length of leaf):
    status=NF90_DEF_VAR(ncid_restart,'dleaf',NF90_FLOAT,(/mpID/), &
         rpid%dleaf)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining dleaf variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%dleaf,"long_name",&
         "Chararacteristic length of leaf")
    status = NF90_PUT_ATT(ncid_restart,rpid%dleaf,"units","m")
    ! frac4 (Fraction of plants which are C4):
    status=NF90_DEF_VAR(ncid_restart,'frac4',NF90_FLOAT,(/mpID/), &
         rpid%frac4)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining frac4 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%frac4,"long_name",&
         "Fraction of plants which are C4")
    status = NF90_PUT_ATT(ncid_restart,rpid%frac4,"units","-")
    ! ejmax (max pot. electron transport rate top leaf):
    status=NF90_DEF_VAR(ncid_restart,'ejmax',NF90_FLOAT,(/mpID/), &
         rpid%ejmax)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ejmax variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%ejmax,"long_name",&
         "Max potential electron transport rate top leaf")
    status = NF90_PUT_ATT(ncid_restart,rpid%ejmax,"units","mol/m^2/s")
    ! vcmax (Maximum RuBP carboxylation rate top leaf):
    status=NF90_DEF_VAR(ncid_restart,'vcmax',NF90_FLOAT,(/mpID/), &
         rpid%vcmax)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining vcmax variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%vcmax,"long_name",&
         "Maximum RuBP carboxylation rate top leaf")
    status = NF90_PUT_ATT(ncid_restart,rpid%vcmax,"units","mol/m^2/s")
    ! rp20 (Plant respiration coefficient at 20C):
    status=NF90_DEF_VAR(ncid_restart,'rp20',NF90_FLOAT,(/mpID/), &
         rpid%rp20)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rp20 variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%rp20,"long_name",&
         "Plant respiration coefficient at 20C")
    status = NF90_PUT_ATT(ncid_restart,rpid%rp20,"units","-")
    ! rpcoef (Temperature coef nonleaf plant respiration):
    status=NF90_DEF_VAR(ncid_restart,'rpcoef',NF90_FLOAT,(/mpID/), &
         rpid%rpcoef)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining rpcoef variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%rpcoef,"long_name",&
         "Temperature coef nonleaf plant respiration")
    status = NF90_PUT_ATT(ncid_restart,rpid%rpcoef,"units","1/C")
    ! shelrb (sheltering factor):
    status=NF90_DEF_VAR(ncid_restart,'shelrb',NF90_FLOAT,(/mpID/), &
         rpid%shelrb)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining shelrb variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%shelrb,"long_name",&
         "sheltering factor")
    status = NF90_PUT_ATT(ncid_restart,rpid%shelrb,"units","-")
    ! xfang (leaf angle parameter):
    status=NF90_DEF_VAR(ncid_restart,'xfang',NF90_FLOAT,(/mpID/), &
         rpid%xfang)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining xfang variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%xfang,"long_name",&
         "Leaf angle parameter")
    status = NF90_PUT_ATT(ncid_restart,rpid%xfang,"units","-")
    ! tminvj (Min temperature for the start of photosynthesis):
    status=NF90_DEF_VAR(ncid_restart,'tminvj',NF90_FLOAT,(/mpID/), &
         rpid%tminvj)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tminvj variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%tminvj,"long_name",&
         "Min temperature for the start of photosynthesis")
    status = NF90_PUT_ATT(ncid_restart,rpid%tminvj,"units","-")
    ! tmaxvj (Max temperature for the start of photosynthesis):
    status=NF90_DEF_VAR(ncid_restart,'tmaxvj',NF90_FLOAT,(/mpID/), &
         rpid%tmaxvj)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining tmaxvj variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%tmaxvj,"long_name",&
         "Max temperature for the start of photosynthesis")
    status = NF90_PUT_ATT(ncid_restart,rpid%tmaxvj,"units","-")
    ! vbeta (Stomatal sensitivity to soil water):
    status=NF90_DEF_VAR(ncid_restart,'vbeta',NF90_FLOAT,(/mpID/), &
         rpid%vbeta)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining vbeta variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%vbeta,"long_name",&
         "Stomatal sensitivity to soil water")
    status = NF90_PUT_ATT(ncid_restart,rpid%vbeta,"units","-")
    ! ratecp (Plant carbon rate constant):
    status=NF90_DEF_VAR(ncid_restart,'ratecp',NF90_FLOAT,(/plantcarbID/), &
         rpid%ratecp)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ratecp variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%ratecp,"long_name",&
         "Plant carbon rate constant")
    status = NF90_PUT_ATT(ncid_restart,rpid%ratecp,"units","1/year")
    ! ratecs (Soil carbon rate constant):
    status=NF90_DEF_VAR(ncid_restart,'ratecs',NF90_FLOAT,(/soilcarbID/), &
         rpid%ratecs)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining ratecs variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%ratecs,"long_name",&
         "Soil carbon rate constant")
    status = NF90_PUT_ATT(ncid_restart,rpid%ratecs,"units","1/year")
    ! meth:
    status=NF90_DEF_VAR(ncid_restart,'meth',NF90_FLOAT,(/mpID/), &
         rpid%meth)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining meth variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%meth,"long_name",&
         "Canopy turbulence parameterisation")
    status = NF90_PUT_ATT(ncid_restart,rpid%meth,"units","-")
    ! za:
    status=NF90_DEF_VAR(ncid_restart,'za',NF90_FLOAT,(/mpID/), &
         rpid%za)
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining za variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,rpid%za,"long_name",&
         "Reference height (lowest atm. model layer)")
    status = NF90_PUT_ATT(ncid_restart,rpid%za,"units","m")

    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    status = NF90_PUT_ATT(ncid_restart,NF90_GLOBAL,"Production", &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_restart_out)// ' (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,NF90_GLOBAL,"Source", &
         'CABLE LSM restart file')
     IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_restart_out)// ' (SUBROUTINE create_restart)')
    status = NF90_PUT_ATT(ncid_restart,NF90_GLOBAL,"CABLE_input_file", &
         TRIM(filename_met))
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing global details to file ' &
         //TRIM(filename_restart_out)// ' (SUBROUTINE create_restart)')
    
    ! End netcdf define mode:
    status = NF90_ENDDEF(ncid_restart)
    IF(status/=NF90_NOERR) CALL nc_abort('Error creating output file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
   
    ! Write time variable:
    status=NF90_PUT_VAR(ncid_restart,tvarID,timevar(kend-kstart+1))
    IF(status/=NF90_NOERR) CALL nc_abort('Error time variable to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')

    ! Write latitude and longitude variables:
    status=NF90_PUT_VAR(ncid_restart,latID,latitude)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing latitude variable to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,lonID,longitude)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing longitude variable to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')

    ! Write parameters to netcdf file:
    status=NF90_PUT_VAR(ncid_restart,rpid%iveg,veg%iveg)        
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing iveg parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%isoil,soil%isoilm)       
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing isoil parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%bch,soil%bch)       
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing bch parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%clay,soil%clay)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing clay parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%sand,soil%sand)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing sand parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%silt,soil%silt)           
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing silt parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%css,soil%css)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing css parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%rhosoil,soil%rhosoil)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rhosoil parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%hyds,soil%hyds)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing hyds parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%sucs,soil%sucs)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing sucs parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%rs20,soil%rs20)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rs20 parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%ssat,soil%ssat)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ssat parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%sfc,soil%sfc)         
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing sfc parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%swilt,soil%swilt)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing swilt parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%froot,veg%froot)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing froot parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%zse,soil%zse)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing zse parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%albsoil,soil%albsoil)         
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing albsoil parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%canst1,veg%canst1)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing canst1 parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%dleaf,veg%dleaf)         
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing dleaf parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%ejmax,veg%ejmax)        
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ejmax parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%vcmax,veg%vcmax)        
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing vcmax parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%frac4,veg%frac4)        
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing frac4 parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%hc,veg%hc)           
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing hc parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%rp20,veg%rp20)           
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rp20 parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%rpcoef,veg%rpcoef)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rpcoef parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%shelrb,veg%shelrb)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing shelrb parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%xfang,veg%xfang)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing xfang parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%tminvj,veg%tminvj)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing tminvj parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%tmaxvj,veg%tmaxvj)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing tmaxvj parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%vbeta,veg%vbeta)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing vbeta parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%ratecp,bgc%ratecp)
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecp parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%ratecs,bgc%ratecs)         
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ratecs parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%meth,veg%meth)          
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing meth parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    status=NF90_PUT_VAR(ncid_restart,rpid%za,rough%za)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing za parameter to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! Write states to restart file:
    ! tgg:
    status=NF90_PUT_VAR(ncid_restart,tggID,ssoil%tgg)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing tgg to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! wb:
    status=NF90_PUT_VAR(ncid_restart,wbID,ssoil%wb)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing wb to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! wbice:
    status=NF90_PUT_VAR(ncid_restart,wbiceID,ssoil%wbice)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing wbice to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! tss:
    status=NF90_PUT_VAR(ncid_restart,tssID,ssoil%tss)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing tss to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! ssdnn:
    status=NF90_PUT_VAR(ncid_restart,ssdnnID,ssoil%ssdnn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ssdnn to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! ssdn:
    status=NF90_PUT_VAR(ncid_restart,ssdnID,ssoil%ssdn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ssdn to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! osnowd:
    status=NF90_PUT_VAR(ncid_restart,osnowdID,ssoil%osnowd)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing osnowd to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! smass:
    status=NF90_PUT_VAR(ncid_restart,smassID,ssoil%smass)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing smass to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! sdepth:
    status=NF90_PUT_VAR(ncid_restart,sdepthID,ssoil%sdepth)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing sdepth to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! snage:
    status=NF90_PUT_VAR(ncid_restart,snageID,ssoil%snage)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing snage to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! snowd:
    status=NF90_PUT_VAR(ncid_restart,snowdID,ssoil%snowd)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing snowd to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! rtsoil:
    status=NF90_PUT_VAR(ncid_restart,rtsoilID,ssoil%rtsoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rtsoil to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! isflag:
    status=NF90_PUT_VAR(ncid_restart,isflagID,ssoil%isflag)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing isflag to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! cansto:
    status=NF90_PUT_VAR(ncid_restart,canstoID,canopy%cansto)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing cansto to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! albsoilsn:
    status=NF90_PUT_VAR(ncid_restart,albsoilsnID,ssoil%albsoilsn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing albsoilsn to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! gammzz:
    status=NF90_PUT_VAR(ncid_restart,gammzzID,ssoil%gammzz)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing gammzz to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! tggsn:
    status=NF90_PUT_VAR(ncid_restart,tggsnID,ssoil%tggsn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing tggsn to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! sghflux:
    status=NF90_PUT_VAR(ncid_restart,sghfluxID,canopy%sghflux)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing sghflux to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! ghflux:
    status=NF90_PUT_VAR(ncid_restart,ghfluxID,canopy%ghflux)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ghflux to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! runoff:
    status=NF90_PUT_VAR(ncid_restart,runoffID,ssoil%runoff)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing runoff to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! rnof1:
    status=NF90_PUT_VAR(ncid_restart,rnof1ID,ssoil%rnof1)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rnof1 to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! rnof2:
    status=NF90_PUT_VAR(ncid_restart,rnof2ID,ssoil%rnof2)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing rnof2 to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! ga:
    status=NF90_PUT_VAR(ncid_restart,gaID,canopy%ga)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing ga to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! dgdtg:
    status=NF90_PUT_VAR(ncid_restart,dgdtgID,canopy%dgdtg)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing dgdtg to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! fev:
    status=NF90_PUT_VAR(ncid_restart,fevID,canopy%fev)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing fev to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! fes:
    status=NF90_PUT_VAR(ncid_restart,fesID,canopy%fes)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing fes to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! fhs:
    status=NF90_PUT_VAR(ncid_restart,fhsID,canopy%fhs)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing fhs to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! cplant:
    status=NF90_PUT_VAR(ncid_restart,cplantID,bgc%cplant)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing cplant to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! csoil:
    status=NF90_PUT_VAR(ncid_restart,csoilID,bgc%csoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing csoil to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! wbtot0:
    status=NF90_PUT_VAR(ncid_restart,wbtot0ID,bal%wbtot0)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing wbtot0 to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')
    ! osnowd0:
    status=NF90_PUT_VAR(ncid_restart,osnowd0ID,bal%osnowd0)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error writing osnowd0 to file ' &
         //TRIM(filename_restart_out)// '(SUBROUTINE create_restart)')

    ! Close restart file
    status = NF90_CLOSE(ncid_restart)

    WRITE(logn,'(A34)') ' Restart file complete and closed.'

  END SUBROUTINE create_restart
  !========================================================================
  SUBROUTINE range_abort(message,ktau,met,value,var_range,xx,yy)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER(i_d), INTENT(IN) :: ktau ! time step
    TYPE(met_type),INTENT(IN) :: met  ! met data
    REAL(r_1),INTENT(IN) :: value ! value deemed to be out of range
    INTEGER(i_d),INTENT(IN) :: xx ! coordinates of erroneous grid square
    INTEGER(i_d),INTENT(IN),OPTIONAL ::yy ! coordinates of erroneous grid square
    REAL(r_1),DIMENSION(2),INTENT(IN) :: var_range ! appropriate var range 
    WRITE(*,*) message ! error from subroutine
    IF(PRESENT(yy)) THEN ! i.e. using rectangular land/sea grid
       WRITE(*,*) 'Site lat, lon:',lat_all(xx,yy),lon_all(xx,yy)
       WRITE(*,*) 'Timestep',ktau, &
             ', or ', met%hod,' hod, ',INT(met%doy),'doy, ',INT(met%year)
    ELSE ! i.e. using compressed land only grid
       WRITE(*,*) 'Site lat, lon:',latitude(xx),longitude(xx)
       WRITE(*,*) 'Timestep',ktau, &
         ', or ', met%hod(xx),' hod, ',INT(met%doy(xx)),'doy, ',INT(met%year(xx))
    END IF
    WRITE(*,*) 'Specified acceptable range (checks.f90):', var_range(1), &
         'to',var_range(2)
    WRITE(*,*) 'Value:',value 
    STOP
  END SUBROUTINE range_abort

END MODULE output_module
