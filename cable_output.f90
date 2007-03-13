!!$ cable_output.f90
!!$
!!$ Output module for CABLE land surface scheme offline netcdf driver; 
!!$
!!$ Gab Abramowitz 2006 CSIRO Marine and Atmospheric
!!$ Research/ Macquarie University; gabsun@gmail.com
!!$
!!$ Subroutines in this file open the output netcdf file for offline CABLE
!!$ and write to it, providing details about variable units and names. 

MODULE output_module
  USE input_module
  IMPLICIT NONE
  INTEGER :: ncid_out ! output data netcdf file ID
  TYPE output_switch_type
     LOGICAL :: met, flux, radiation, carbon, soil, veg, &
          params, balances
  END TYPE output_switch_type
  TYPE(output_switch_type) :: output
  INTEGER :: tvarID ! time variable ID
  TYPE out_varID_type ! output variable IDs in netcdf file
     INTEGER :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair,Qair,Rainf, &
          Snowf,CO2air,Qle,Qh,Qg,NEE,SWnet,LWnet,SoilMoist, &
          SoilTemp,Albedo,Qs,Qsb,Evap,BaresoilT,SWE,SnowT,RadT, &
          VegT,Ebal,Wbal,AutoResp,HeteroResp,GPP,NPP,LAI,ECanop, &
          TVeg,ESoil,CanopInt,SnowDepth,HVeg,HSoil
  END TYPE out_varID_type
  TYPE(out_varID_type) :: ovid ! netcdf variable IDs for output variables
  TYPE out_parID_type ! output model parameter IDs in netcdf file
     INTEGER :: bch,latitude,c3,clay,css,rhosoil,hyds,rs20,sand,sfc,silt, &
         ssat,sucs,swilt,froot,zse,canst1,dleaf,meth,za, &
         ejmax,frac4,hc,lai,rp20,rpcoef,shelrb, vbeta, &
         vcmax,xfang,ratecp,ratecs,refsbare,isoil,iveg,albsoil,&
         taul,refl,tauw,refw,tminvj,tmaxvj,veg_class,soil_class
  END TYPE out_parID_type
  TYPE(out_parID_type) :: opid ! netcdf variable IDs for output variables
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
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: HVeg         ! sensible heat from vegetation [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: HSoil        ! sensible heat from soil [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Ebal         ! cumulative energy balance [W/m2]
     REAL(r_1), POINTER, DIMENSION(:,:,:) :: Wbal         ! cumulative water balance [W/m2]
     ! Model parameters
     REAL(r_1),POINTER,DIMENSION(:,:) :: bch       ! parameter b in Campbell equation 1985
     REAL(r_1),POINTER,DIMENSION(:,:) :: latitude  ! site latitude
     ! Frac of lowest layer soil moisture above field capacity which drains
     REAL(r_1),POINTER,DIMENSION(:,:) :: c3        ! c3 drainage coeff (fraction)
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
     REAL(r_1),POINTER,DIMENSION(:,:) :: vbeta     ! ???
     REAL(r_1),POINTER,DIMENSION(:,:) :: iveg      ! vegetation type from global index
     REAL(r_1),POINTER,DIMENSION(:,:) :: isoil     ! soil type from global index
     REAL(r_1),POINTER,DIMENSION(:,:) :: meth      ! method for solving ??? in canopy scheme
     REAL(r_1),POINTER,DIMENSION(:,:) :: za        ! something to do with roughness ????
  END TYPE dataset_type
  TYPE (dataset_type) :: out ! CABLE output after units adjustment

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
    INTEGER :: xID,yID,zID,radID,soilID,soilcarbID,plantcarbID,tID ! dimension IDs
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
    status = NF90_DEF_DIM(ncid_out,'z',1,zID) ! number of soil layers
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error defining z dimension in output file. '// &
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
    ! Define meteorological variables in output file if requested:
    IF(output%met) THEN 
       WRITE(logn,*) 'Writing met input data to output file'
       ! Define SWdown and units:
       status=NF90_DEF_VAR(ncid_out,'SWdown',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%SWdown)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWdown variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWdown,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWdown,'long_name', &
            'Downward shortwave radiation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SWdown(xdimsize,ydimsize,1))
       ! Define LWdown and units:
       status=NF90_DEF_VAR(ncid_out,'LWdown',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%LWdown)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWdown variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWdown,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWdown,'long_name', &
            'Downward longwave radiation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWdown variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%LWdown(xdimsize,ydimsize,1))
       ! Define Tair and units:
       status=NF90_DEF_VAR(ncid_out,'Tair',NF90_FLOAT,(/xID,yID,zID,tID/), &
            ovid%Tair)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Tair variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Tair,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Tair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Tair,'long_name', &
            'Surface air temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Tair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Tair(xdimsize,ydimsize,1,1))
       ! Define Rainf and units:
       status=NF90_DEF_VAR(ncid_out,'Rainf',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Rainf)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rainf variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Rainf,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rainf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Rainf,'long_name', &
            'Rainfall AND snowfall')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Rainf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Rainf(xdimsize,ydimsize,1))
        ! Define Snowf and units:
       status=NF90_DEF_VAR(ncid_out,'Snowf',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Snowf)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Snowf variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Snowf,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Snowf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Snowf,'long_name', &
            'Snowfall')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Snowf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Snowf(xdimsize,ydimsize,1))
       ! Define Qair and units:
       status=NF90_DEF_VAR(ncid_out,'Qair',NF90_FLOAT,(/xID,yID,zID,tID/), &
            ovid%Qair)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qair variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qair,'units','kg/kg')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qair,'long_name', &
            'Surface specific humidity')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qair variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qair(xdimsize,ydimsize,1,1))
       ! Define Wind and units:
       status=NF90_DEF_VAR(ncid_out,'Wind',NF90_FLOAT,(/xID,yID,zID,tID/), &
            ovid%Wind)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wind variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Wind,'units','m/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wind variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Wind,'long_name', &
            'Scalar surface wind speed')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Wind variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Wind(xdimsize,ydimsize,1,1))
       ! Define PSurf and units:
       status=NF90_DEF_VAR(ncid_out,'PSurf',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%PSurf)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining PSurf variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%PSurf,'units','hPa')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining PSurf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%PSurf,'long_name', &
            'Surface air pressure')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining PSurf variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%PSurf(xdimsize,ydimsize,1))
       ! Define CO2air and units:
       status=NF90_DEF_VAR(ncid_out,'CO2air',NF90_FLOAT,(/xID,yID,zID,tID/), &
            ovid%CO2air)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CO2air variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%CO2air,'units','ppmv')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CO2air variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%CO2air,'long_name', &
            'Surface air CO2 concentration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CO2air variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%CO2air(xdimsize,ydimsize,1,1))
    END IF
    !----------------------FLUX VARIABLES------------------------------------
    IF(output%flux) THEN ! output fluxes
       WRITE(logn,*) 'Writing surface fluxes to output file'
       ! Define Qle and units (latent heat):
       status=NF90_DEF_VAR(ncid_out,'Qle',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Qle)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qle variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qle,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qle variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qle,'long_name', &
            'Surface latent heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qle variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qle(xdimsize,ydimsize,1))
       ! Define Qh and units (sensible heat):
       status=NF90_DEF_VAR(ncid_out,'Qh',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Qh)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qh variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qh,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qh variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qh,'long_name', &
            'Surface sensible heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qh variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qh(xdimsize,ydimsize,1))
       ! Define Qg and units (ground heat flux):
       status=NF90_DEF_VAR(ncid_out,'Qg',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Qg)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qg variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qg,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qg,'long_name', &
            'Ground heat flux')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qg(xdimsize,ydimsize,1))
       ! Define Qs and units (surface runoff):
       status=NF90_DEF_VAR(ncid_out,'Qs',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Qs)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qs variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qs,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qs,'long_name', &
            'Surface runoff')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qs variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qs(xdimsize,ydimsize,1))
       ! Define Qsb and units (subsurface runoff):
       status=NF90_DEF_VAR(ncid_out,'Qsb',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Qsb)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qsb variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Qsb,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qsb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%Qsb,'long_name', &
            'Subsurface runoff')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Qsb variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Qsb(xdimsize,ydimsize,1))
       ! Define Evap and units (total evapotranspiration):
       status=NF90_DEF_VAR(ncid_out,'Evap',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Evap)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Evap variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Evap,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Evap variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Evap,'long_name', &
            'Total evapotranspiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Evap variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Evap(xdimsize,ydimsize,1))
       ! Define ECanop and units (evaporation from wet canopy):
       status=NF90_DEF_VAR(ncid_out,'ECanop',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%ECanop)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ECanop variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%ECanop,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ECanop variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%ECanop,'long_name', &
            'Wet canopy evaporation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ECanop variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%ECanop(xdimsize,ydimsize,1))
       ! Define TVeg and units (vegetation transpiration):
       status=NF90_DEF_VAR(ncid_out,'TVeg',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%TVeg)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining TVeg variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%TVeg,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining TVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%TVeg,'long_name', &
            'Vegetation transpiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining TVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%TVeg(xdimsize,ydimsize,1))
       ! Define ESoil and units (soil evaporation):
       status=NF90_DEF_VAR(ncid_out,'ESoil',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%ESoil)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ESoil variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%ESoil,'units','kg/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ESoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%ESoil,'long_name', &
            'Evaporation from soil')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ESoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%ESoil(xdimsize,ydimsize,1))
       ! Define HVeg and units (sensible heat from vegetation):
       status=NF90_DEF_VAR(ncid_out,'HVeg',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%HVeg)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HVeg variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HVeg,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%HVeg,'long_name', &
            'Sensible heat from vegetation')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HVeg variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%HVeg(xdimsize,ydimsize,1))
       ! Define HSoil and units (sensible heat from soil):
       status=NF90_DEF_VAR(ncid_out,'HSoil',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%HSoil)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HSoil variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HSoil,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HSoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%HSoil,'long_name', &
            'Sensible heat from soil')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HSoil variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%HSoil(xdimsize,ydimsize,1))
       ! Define NEE and units (net ecosystem exchange):
       status=NF90_DEF_VAR(ncid_out,'NEE',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%NEE)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NEE variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%NEE,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NEE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%NEE,'long_name', &
            'Net ecosystem exchange of CO2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NEE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%NEE(xdimsize,ydimsize,1))
    END IF
    !------------------SOIL STATES--------------------------------------
    IF(output%soil) THEN ! output soil variables
        WRITE(logn,*) 'Writing soil data to output file'
       ! Define SoilMoist and units:
       status=NF90_DEF_VAR(ncid_out,'SoilMoist',NF90_FLOAT,(/xID,yID,soilID,tID/), &
            ovid%SoilMoist)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilMoist variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SoilMoist,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilMoist variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%SoilMoist,'long_name', &
            'Average layer soil moisture')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilMoist variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SoilMoist(xdimsize,ydimsize,ms,1))
       ! Define SoilTemp and units:
       status=NF90_DEF_VAR(ncid_out,'SoilTemp',NF90_FLOAT,(/xID,yID,soilID,tID/), &
            ovid%SoilTemp)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilTemp variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SoilTemp,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilTemp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
        status = NF90_PUT_ATT(ncid_out,ovid%SoilTemp,'long_name', &
            'Average layer soil temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SoilTemp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SoilTemp(xdimsize,ydimsize,ms,1))
       ! Define BaresoilT and units:
       status=NF90_DEF_VAR(ncid_out,'BaresoilT',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%BaresoilT)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining BaresoilT variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%BaresoilT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining BaresoilT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%BaresoilT,'long_name', &
            'Bare soil temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining BaresoilT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%BaresoilT(xdimsize,ydimsize,1))
       ! Define SWE and units (snow water equivalent):
       status=NF90_DEF_VAR(ncid_out,'SWE',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%SWE)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWE variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWE,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWE,'long_name', &
            'Snow water equivalent')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWE variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SWE(xdimsize,ydimsize,1))
       ! Define SnowT and units:
       status=NF90_DEF_VAR(ncid_out,'SnowT',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%SnowT)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowT variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowT,'long_name', &
            'Snow surface temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SnowT(xdimsize,ydimsize,1))
       ! Define SnowDepth and units:
       status=NF90_DEF_VAR(ncid_out,'SnowDepth',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%SnowDepth)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowDepth variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowDepth,'units','m')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowDepth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SnowDepth,'long_name', &
            'Snow depth')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SnowDepth variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SnowDepth(xdimsize,ydimsize,1))
    END IF
    !------------------RADIATION VARIABLES--------------------------------------
    IF(output%radiation) THEN ! output soil variables
       WRITE(logn,*) 'Writing radiation data to output file'
       ! Define SWnet and units:
       status=NF90_DEF_VAR(ncid_out,'SWnet',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%SWnet)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWnet variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWnet,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%SWnet,'long_name', &
            'Net shortwave radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining SWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%SWnet(xdimsize,ydimsize,1))
       ! Define LWnet and units:
       status=NF90_DEF_VAR(ncid_out,'LWnet',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%LWnet)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWnet variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWnet,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LWnet,'long_name', &
            'Net longwave radiation absorbed by surface')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LWnet variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%LWnet(xdimsize,ydimsize,1))
       ! Define Albedo and units:
       status=NF90_DEF_VAR(ncid_out,'Albedo',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Albedo)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Albedo variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Albedo,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Albedo variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Albedo,'long_name', &
            'Surface net albedo')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining Albedo variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Albedo(xdimsize,ydimsize,1))
       ! Define RadT and units:
       status=NF90_DEF_VAR(ncid_out,'RadT',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%RadT)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining RadT variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%RadT,'units','K')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining RadT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%RadT,'long_name', &
            'Radiative surface temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining RadT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%RadT(xdimsize,ydimsize,1))
    END IF
    !--------------------VEGETATION VARIABLES------------------------------
    IF(output%veg) THEN
       WRITE(logn,*) 'Writing vegetation data to output file'
       ! Define VegT and units:
       status=NF90_DEF_VAR(ncid_out,'VegT',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%VegT)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining VegT variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%VegT,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining VegT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%VegT,'long_name', &
            'Average vegetation temperature')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining VegT variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%VegT(xdimsize,ydimsize,1))
       ! Define CanopInt and units (canopy water storage):
       status=NF90_DEF_VAR(ncid_out,'CanopInt',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%CanopInt)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CanopInt variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%CanopInt,'units','kg/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CanopInt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%CanopInt,'long_name', &
            'Canopy water storage')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining CanopInt variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%CanopInt(xdimsize,ydimsize,1))
       ! Define LAI and units:
       status=NF90_DEF_VAR(ncid_out,'LAI',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%LAI)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LAI variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LAI,'units','-')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LAI variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%LAI,'long_name', &
            'Leaf area index')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining LAI variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%LAI(xdimsize,ydimsize,1))
    END IF
    !------------------------CONSERVATION VARIABLES---------------------------
    IF(output%balances) THEN
       WRITE(logn,*) 'Writing water and energy balance data to output file'
       ! Define Ebal and units:
       status=NF90_DEF_VAR(ncid_out,'Ebal',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Ebal)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining energy balance variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Ebal,'units','W/m^2')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining energy balance variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Ebal,'long_name', &
            'Cumulative energy balance')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining energy balance variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Ebal(xdimsize,ydimsize,1))
       ! Define Wbal and units:
       status=NF90_DEF_VAR(ncid_out,'Wbal',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%Wbal)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining water balance variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Wbal,'units','mm')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining water balance variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%Wbal,'long_name', &
            'Cumulative water balance')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining water balance variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%Wbal(xdimsize,ydimsize,1))
    END IF
    !--------------------CARBON VARIABLES------------------------------
    IF(output%carbon) THEN
       WRITE(logn,*) 'Writing carbon data to output file'
       ! Check if NEE is already being reported in "fluxes" output
       IF(.NOT.output%flux) THEN
          ! Define NEE and units (net ecosystem exchange):
          status=NF90_DEF_VAR(ncid_out,'NEE',NF90_FLOAT,(/xID,yID,tID/), &
               ovid%NEE)
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NEE variable in output file. '// &
               '(SUBROUTINE open_output_file)')
          status = NF90_PUT_ATT(ncid_out,ovid%NEE,'units','umol/m^2/s')
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NEE variable attributes in output file. '// &
               '(SUBROUTINE open_output_file)')
          status = NF90_PUT_ATT(ncid_out,ovid%NEE,'long_name', &
               'Net ecosystem exchange of CO2')
          IF (status /= NF90_NOERR) CALL nc_abort &
               ('Error defining NEE variable attributes in output file. '// &
               '(SUBROUTINE open_output_file)')
          ALLOCATE(out%NEE(xdimsize,ydimsize,1))
       END IF
       ! Define AutoResp and units:
       status=NF90_DEF_VAR(ncid_out,'AutoResp',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%AutoResp)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining AutoResp variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%AutoResp,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining AutoResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%AutoResp,'long_name', &
            'Autotrophic respiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining AutoResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%AutoResp(xdimsize,ydimsize,1))
       ! Define HeteroResp and units:
       status=NF90_DEF_VAR(ncid_out,'HeteroResp',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%HeteroResp)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HeteroResp variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HeteroResp,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HeteroResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%HeteroResp,'long_name', &
            'Heterotrophic respiration')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining HeteroResp variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%HeteroResp(xdimsize,ydimsize,1))
       ! Define GPP and units:
       status=NF90_DEF_VAR(ncid_out,'GPP',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%GPP)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining GPP variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%GPP,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining GPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%GPP,'long_name', &
            'Gross primary production')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining GPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%GPP(xdimsize,ydimsize,1))
       ! Define NPP and units:
       status=NF90_DEF_VAR(ncid_out,'NPP',NF90_FLOAT,(/xID,yID,tID/), &
            ovid%NPP)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NPP variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%NPP,'units','umol/m^2/s')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,ovid%NPP,'long_name', &
            'Net primary production')
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining NPP variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ALLOCATE(out%NPP(xdimsize,ydimsize,1))
    END IF
    !---------------------MODEL PARAMETERS---------------------------------
    IF(output%params) THEN
       WRITE(logn,*) 'Writing model parameters to output file'
       ! Allocate size of output parameter variables:
       ! iveg (Vegetation type):
       status=NF90_DEF_VAR(ncid_out,'iveg',NF90_FLOAT,(/xID,yID/), &
            opid%iveg)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining iveg variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%iveg,"long_name",&
            "Vegetation type")
       status = NF90_PUT_ATT(ncid_out,opid%iveg,"units","-")
       ALLOCATE(out%iveg(xdimsize,ydimsize))
       ! isoil (Soil type):
       status=NF90_DEF_VAR(ncid_out,'isoil',NF90_FLOAT,(/xID,yID/), &
            opid%isoil)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining isoil variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%isoil,"long_name",&
            "Soil type")
       status = NF90_PUT_ATT(ncid_out,opid%isoil,"units","-")
       ALLOCATE(out%isoil(xdimsize,ydimsize))
       ! clay (fraction of soil which is clay):
       status=NF90_DEF_VAR(ncid_out,'clay',NF90_FLOAT,(/xID,yID/), &
            opid%clay)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining clay variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%clay,"long_name",&
            "Fraction of soil which is clay")
       status = NF90_PUT_ATT(ncid_out,opid%clay,"units","-")
        ALLOCATE(out%clay(xdimsize,ydimsize))
       ! sand:
       status=NF90_DEF_VAR(ncid_out,'sand',NF90_FLOAT,(/xID,yID/), &
            opid%sand)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sand variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sand,"long_name",&
            "Fraction of soil which is sand")
       status = NF90_PUT_ATT(ncid_out,opid%sand,"units","-")
       ALLOCATE(out%sand(xdimsize,ydimsize))
       ! silt :
       status=NF90_DEF_VAR(ncid_out,'silt',NF90_FLOAT,(/xID,yID/), &
            opid%silt)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining silt variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%silt,"long_name",&
            "Fraction of soil which is silt")
       status = NF90_PUT_ATT(ncid_out,opid%silt,"units","-")
       ALLOCATE(out%silt(xdimsize,ydimsize))
       ! ssat (Volume of soil volume which is water @ saturation):
       status=NF90_DEF_VAR(ncid_out,'ssat',NF90_FLOAT,(/xID,yID/), &
            opid%ssat)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ssat variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ssat,"long_name",&
            "Fraction of soil volume which is water @ saturation")
       status = NF90_PUT_ATT(ncid_out,opid%ssat,"units","-")
       ALLOCATE(out%ssat(xdimsize,ydimsize))
       ! sfc (Volume of soil volume which is water @ field capacity):
       status=NF90_DEF_VAR(ncid_out,'sfc',NF90_FLOAT,(/xID,yID/), &
            opid%sfc)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sfc variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sfc,"long_name",&
            "Fraction of soil volume which is water @ field capacity")
       status = NF90_PUT_ATT(ncid_out,opid%sfc,"units","-")
       ALLOCATE(out%sfc(xdimsize,ydimsize))
       ! swilt (Volume of soil volume which is water @ wilting point):
       status=NF90_DEF_VAR(ncid_out,'swilt',NF90_FLOAT,(/xID,yID/), &
            opid%swilt)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining swilt variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%swilt,"long_name",&
            "Fraction of soil volume which is water @ wilting point")
       status = NF90_PUT_ATT(ncid_out,opid%swilt,"units","-")
       ALLOCATE(out%swilt(xdimsize,ydimsize))
       ! zse (depth of each soil layer):
       status=NF90_DEF_VAR(ncid_out,'zse',NF90_FLOAT,(/xID,yID,soilID/), &
            opid%zse)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining zse variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%zse,"long_name",&
            "Depth of each soil layer")
       status = NF90_PUT_ATT(ncid_out,opid%zse,"units","m")
       ALLOCATE(out%zse(xdimsize,ydimsize,ms))
       ! froot (fraction of roots in each soil layer):
       status=NF90_DEF_VAR(ncid_out,'froot',NF90_FLOAT,(/xID,yID,soilID/), &
            opid%froot)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining froot variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%froot,"long_name",&
            "Fraction of roots in each soil layer")
       status = NF90_PUT_ATT(ncid_out,opid%froot,"units","-")
       ALLOCATE(out%froot(xdimsize,ydimsize,ms))
       ! bch (Parameter b, Campbell eqn 1985):
       status=NF90_DEF_VAR(ncid_out,'bch',NF90_FLOAT,(/xID,yID/), &
            opid%bch)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining bch variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%bch,"long_name","Parameter b, Campbell eqn 1985")
       status = NF90_PUT_ATT(ncid_out,opid%bch,"units","-")
       ALLOCATE(out%bch(xdimsize,ydimsize))
       ! c3 (Fraction of bottom layer soil moisture > field capacity which drains):
       status=NF90_DEF_VAR(ncid_out,'c3',NF90_FLOAT,(/xID,yID/), &
            opid%c3)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining c3 variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%c3,"long_name",&
            "Fraction of bottom layer soil moisture > field capacity which drains")
       status = NF90_PUT_ATT(ncid_out,opid%c3,"units","-")
       ALLOCATE(out%c3(xdimsize,ydimsize))
       ! hyds (Hydraulic conductivity @ saturation):
       status=NF90_DEF_VAR(ncid_out,'hyds',NF90_FLOAT,(/xID,yID/), &
            opid%hyds)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hyds variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%hyds,"long_name",&
            "Hydraulic conductivity @ saturation")
       status = NF90_PUT_ATT(ncid_out,opid%hyds,"units","m/s")
       ALLOCATE(out%hyds(xdimsize,ydimsize))
       ! sucs (suction @ saturation):
       status=NF90_DEF_VAR(ncid_out,'sucs',NF90_FLOAT,(/xID,yID/), &
            opid%sucs)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining sucs variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%sucs,"long_name",&
            "Suction @ saturation")
       status = NF90_PUT_ATT(ncid_out,opid%sucs,"units","m")
       ALLOCATE(out%sucs(xdimsize,ydimsize))
       ! css (Heat capacity of soil minerals):
       status=NF90_DEF_VAR(ncid_out,'css',NF90_FLOAT,(/xID,yID/), &
            opid%css)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining css variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%css,"long_name",&
            "Heat capacity of soil minerals")
       status = NF90_PUT_ATT(ncid_out,opid%css,"units","J/kg/C")
       ALLOCATE(out%css(xdimsize,ydimsize))
       ! rhosoil (density of soil minerals):
       status=NF90_DEF_VAR(ncid_out,'rhosoil',NF90_FLOAT,(/xID,yID/), &
            opid%rhosoil)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rhosoil variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rhosoil,"long_name",&
            "Density of soil minerals")
       status = NF90_PUT_ATT(ncid_out,opid%rhosoil,"units","kg/m^3")
       ALLOCATE(out%rhosoil(xdimsize,ydimsize))
       ! rs20 (Soil respiration coefficient at 20C):
       status=NF90_DEF_VAR(ncid_out,'rs20',NF90_FLOAT,(/xID,yID/), &
            opid%rs20)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rs20 variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rs20,"long_name",&
            "Soil respiration coefficient at 20C")
       status = NF90_PUT_ATT(ncid_out,opid%rs20,"units","-")
       ALLOCATE(out%rs20(xdimsize,ydimsize))
       ! albsoil (Soil reflectance):
       status=NF90_DEF_VAR(ncid_out,'albsoil',NF90_FLOAT,(/xID,yID/), &
            opid%albsoil)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining albsoil variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%albsoil,"long_name",&
            "Soil reflectance")
       status = NF90_PUT_ATT(ncid_out,opid%albsoil,"units","-")
       ALLOCATE(out%albsoil(xdimsize,ydimsize))
       ! hc (height of canopy):
       status=NF90_DEF_VAR(ncid_out,'hc',NF90_FLOAT,(/xID,yID/), &
            opid%hc)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining hc variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%hc,"long_name",&
            "Height of canopy")
       status = NF90_PUT_ATT(ncid_out,opid%hc,"units","mm/LAI")
       ALLOCATE(out%hc(xdimsize,ydimsize))
       ! canst1 (Max water intercepted by canopy):
       status=NF90_DEF_VAR(ncid_out,'canst1',NF90_FLOAT,(/xID,yID/), &
            opid%canst1)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining canst1 variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%canst1,"long_name",&
            "Max water intercepted by canopy")
       status = NF90_PUT_ATT(ncid_out,opid%canst1,"units","mm/LAI")
       ALLOCATE(out%canst1(xdimsize,ydimsize))
       ! dleaf (Chararacteristic length of leaf):
       status=NF90_DEF_VAR(ncid_out,'dleaf',NF90_FLOAT,(/xID,yID/), &
            opid%dleaf)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining dleaf variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%dleaf,"long_name",&
            "Chararacteristic length of leaf")
       status = NF90_PUT_ATT(ncid_out,opid%dleaf,"units","m")
       ALLOCATE(out%dleaf(xdimsize,ydimsize))
       ! frac4 (Fraction of plants which are C4):
       status=NF90_DEF_VAR(ncid_out,'frac4',NF90_FLOAT,(/xID,yID/), &
            opid%frac4)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining frac4 variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%frac4,"long_name",&
            "Fraction of plants which are C4")
       status = NF90_PUT_ATT(ncid_out,opid%frac4,"units","-")
       ALLOCATE(out%frac4(xdimsize,ydimsize))
       ! ejmax (max pot. electron transport rate top leaf):
       status=NF90_DEF_VAR(ncid_out,'ejmax',NF90_FLOAT,(/xID,yID/), &
            opid%ejmax)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ejmax variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ejmax,"long_name",&
            "Max potential electron transport rate top leaf")
       status = NF90_PUT_ATT(ncid_out,opid%ejmax,"units","mol/m^2/s")
       ALLOCATE(out%ejmax(xdimsize,ydimsize))
       ! vcmax (Maximum RuBP carboxylation rate top leaf):
       status=NF90_DEF_VAR(ncid_out,'vcmax',NF90_FLOAT,(/xID,yID/), &
            opid%vcmax)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vcmax variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%vcmax,"long_name",&
            "Maximum RuBP carboxylation rate top leaf")
       status = NF90_PUT_ATT(ncid_out,opid%vcmax,"units","mol/m^2/s")
       ALLOCATE(out%vcmax(xdimsize,ydimsize))
       ! rp20 (Plant respiration coefficient at 20C):
       status=NF90_DEF_VAR(ncid_out,'rp20',NF90_FLOAT,(/xID,yID/), &
            opid%rp20)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rp20 variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rp20,"long_name",&
            "Plant respiration coefficient at 20C")
       status = NF90_PUT_ATT(ncid_out,opid%rp20,"units","-")
       ALLOCATE(out%rp20(xdimsize,ydimsize))
       ! rpcoef (Temperature coef nonleaf plant respiration):
       status=NF90_DEF_VAR(ncid_out,'rpcoef',NF90_FLOAT,(/xID,yID/), &
            opid%rpcoef)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining rpcoef variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%rpcoef,"long_name",&
            "Temperature coef nonleaf plant respiration")
       status = NF90_PUT_ATT(ncid_out,opid%rpcoef,"units","1/C")
       ALLOCATE(out%rpcoef(xdimsize,ydimsize))
       ! shelrb (sheltering factor):
       status=NF90_DEF_VAR(ncid_out,'shelrb',NF90_FLOAT,(/xID,yID/), &
            opid%shelrb)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining shelrb variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%shelrb,"long_name",&
            "sheltering factor")
       status = NF90_PUT_ATT(ncid_out,opid%shelrb,"units","-")
       ALLOCATE(out%shelrb(xdimsize,ydimsize))
       ! xfang (leaf angle parameter):
       status=NF90_DEF_VAR(ncid_out,'xfang',NF90_FLOAT,(/xID,yID/), &
            opid%xfang)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining xfang variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%xfang,"long_name",&
            "Leaf angle parameter")
       status = NF90_PUT_ATT(ncid_out,opid%xfang,"units","-")
       ALLOCATE(out%xfang(xdimsize,ydimsize))
       ! tminvj (Min temperature for the start of photosynthesis):
       status=NF90_DEF_VAR(ncid_out,'tminvj',NF90_FLOAT,(/xID,yID/), &
            opid%tminvj)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tminvj variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%tminvj,"long_name",&
            "Min temperature for the start of photosynthesis")
       status = NF90_PUT_ATT(ncid_out,opid%tminvj,"units","-")
       ALLOCATE(out%tminvj(xdimsize,ydimsize))
       ! tmaxvj (Max temperature for the start of photosynthesis):
       status=NF90_DEF_VAR(ncid_out,'tmaxvj',NF90_FLOAT,(/xID,yID/), &
            opid%tmaxvj)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining tmaxvj variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%tmaxvj,"long_name",&
            "Max temperature for the start of photosynthesis")
       status = NF90_PUT_ATT(ncid_out,opid%tmaxvj,"units","-")
       ALLOCATE(out%tmaxvj(xdimsize,ydimsize))
       ! vbeta (???):
       status=NF90_DEF_VAR(ncid_out,'vbeta',NF90_FLOAT,(/xID,yID/), &
            opid%vbeta)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining vbeta variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%vbeta,"long_name",&
            "???")
       status = NF90_PUT_ATT(ncid_out,opid%vbeta,"units","-")
       ALLOCATE(out%vbeta(xdimsize,ydimsize))
       ! ratecp (Plant carbon rate constant):
       status=NF90_DEF_VAR(ncid_out,'ratecp',NF90_FLOAT,(/xID,yID,plantcarbID/), &
            opid%ratecp)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecp variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ratecp,"long_name",&
            "Plant carbon rate constant")
       status = NF90_PUT_ATT(ncid_out,opid%ratecp,"units","1/year")
       ALLOCATE(out%ratecp(xdimsize,ydimsize,ncp))
       ! ratecs (Soil carbon rate constant):
       status=NF90_DEF_VAR(ncid_out,'ratecs',NF90_FLOAT,(/xID,yID,soilcarbID/), &
            opid%ratecs)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining ratecs variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%ratecs,"long_name",&
            "Soil carbon rate constant")
       status = NF90_PUT_ATT(ncid_out,opid%ratecs,"units","1/year")
       ALLOCATE(out%ratecs(xdimsize,ydimsize,ncs))
       ! meth:
       status=NF90_DEF_VAR(ncid_out,'meth',NF90_FLOAT,(/xID,yID/), &
            opid%meth)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining meth variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%meth,"long_name",&
            "Canopy turbulence parameterisation")
       status = NF90_PUT_ATT(ncid_out,opid%meth,"units","-")
       ALLOCATE(out%meth(xdimsize,ydimsize))
       ! za:
       status=NF90_DEF_VAR(ncid_out,'za',NF90_FLOAT,(/xID,yID/), &
            opid%za)
       IF (status /= NF90_NOERR) CALL nc_abort &
            ('Error defining za variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       status = NF90_PUT_ATT(ncid_out,opid%za,"long_name",&
            "Reference height (lowest atm. model layer)")
       status = NF90_PUT_ATT(ncid_out,opid%za,"units","m")
       ALLOCATE(out%za(xdimsize,ydimsize))
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
       ! Define global grid parameter variables
       DO i = 1, mp
          out%iveg(land_x(i),land_y(i)) = veg%iveg(i)
          out%isoil(land_x(i),land_y(i)) = soil%isoilm(i)
          out%bch(land_x(i),land_y(i)) = soil%bch(i)
          out%c3(land_x(i),land_y(i)) = soil%c3(i)
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
          out%froot(land_x(i),land_y(i),:) = soil%froot(i,:)
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
       status=NF90_PUT_VAR(ncid_out,opid%c3,out%c3,start=(/1,1/), &
            count=(/xdimsize,ydimsize/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing c3 parameter to file ' &
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

       ! Deallocate output parameter variables, since they've now been written to file:
       DEALLOCATE(out%bch,out%c3,out%clay,out%sand,out%silt,out%css,out%rhosoil, &
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
    IF(output%met) THEN
       ! Write to ALMA-style output met variables:
       DO i = 1, mp
          out%SWdown(land_x(i),land_y(i),1) = met%fsd(i)
          out%LWdown(land_x(i),land_y(i),1) = met%fld(i)
          out%Rainf(land_x(i),land_y(i),1) = met%precip(i)/dels
          out%Snowf(land_x(i),land_y(i),1) = met%precips(i)/dels
          out%PSurf(land_x(i),land_y(i),1) = met%pmb(i)
          out%Tair(land_x(i),land_y(i),1,1) = met%tk(i)
          out%Qair(land_x(i),land_y(i),1,1) = met%qv(i)
          out%Wind(land_x(i),land_y(i),1,1) = met%ua(i)
          out%CO2air(land_x(i),land_y(i),1,1) = met%ca(i)*1000000.0
       END DO
       WHERE(mask/=1) ! where not land:
          out%SWdown(:,:,1) = -9.0
          out%LWdown(:,:,1) = 100.0
          out%Rainf(:,:,1) = -0.01
          out%Snowf(:,:,1) = -0.01
          out%PSurf(:,:,1) = 1200.0
          out%Tair(:,:,1,1) = 180.0
          out%Qair(:,:,1,1) = 0.0
          out%Wind(:,:,1,1) = -1.0
          out%CO2air(:,:,1,1) = 100.0
       END WHERE
       ! Ranges of met variables have been checked when read in.
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%SWdown,out%SWdown,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWdown variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%LWdown,out%LWdown,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LWdown variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Rainf,out%Rainf,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Rainf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Snowf,out%Snowf,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Snowf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%PSurf,out%PSurf,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing PSurf variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Tair,out%Tair,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,1,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Tair variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Qair,out%Qair,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,1,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qair variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Wind,out%Wind,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,1,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Wind variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%CO2air,out%CO2air,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,1,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing CO2air variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    !--------------------------WRITE FLUX DATA-------------------------------------
    IF(output%flux) THEN ! If asked to write flux variables to file:
       ! Write to ALMA-style output flux variables:
       DO i = 1, mp
          out%Qle(land_x(i),land_y(i),1) = canopy%fe(i)
          out%Qh(land_x(i),land_y(i),1) = canopy%fh(i)
          out%Qg(land_x(i),land_y(i),1) = canopy%ga(i)
          out%Qs(land_x(i),land_y(i),1) = ssoil%rnof1(i)/dels
          out%Qsb(land_x(i),land_y(i),1) = ssoil%rnof2(i)/dels
          out%Evap(land_x(i),land_y(i),1) = canopy%fe(i)/air%rlam(i) ! W/m2 to kg/m2/s
          out%ECanop(land_x(i),land_y(i),1) = canopy%fevw(i)/air%rlam(i) ! W/m2 to kg/m2/s
          out%TVeg(land_x(i),land_y(i),1) = canopy%fevc(i)/air%rlam(i) ! W/m2 to kg/m2/s
          out%ESoil(land_x(i),land_y(i),1) = canopy%fes(i)/air%rlam(i) ! W/m2 to kg/m2/s
          out%HVeg(land_x(i),land_y(i),1) = canopy%fhv(i)
          out%HSoil(land_x(i),land_y(i),1) = canopy%fhs(i)
          out%NEE(land_x(i),land_y(i),1) = canopy%fnee(i)/1.201E-5 ! g/m2/s to umol/m2/s       
          IF(check%ranges) THEN  ! Check ranges:
             ! Latent heat:
             IF((out%Qle(land_x(i),land_y(i),1)<ranges%Qle(1)).OR. &
               (out%Qle(land_x(i),land_y(i),1)>ranges%Qle(2))) &
               CALL range_abort('Qle out of specified ranges!', ktau, met, &
               out%Qle(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Qle)
             ! sensible heat:
             IF((out%Qh(land_x(i),land_y(i),1)<ranges%Qh(1)).OR. &
               (out%Qh(land_x(i),land_y(i),1)>ranges%Qh(2))) &
               CALL range_abort('Qh out of specified ranges!', ktau, met, &
               out%Qh(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Qh)
             ! ground heat:
             IF((out%Qg(land_x(i),land_y(i),1)<ranges%Qg(1)).OR. &
               (out%Qg(land_x(i),land_y(i),1)>ranges%Qg(2))) &
               CALL range_abort('Qg out of specified ranges!', ktau, met, &
               out%Qg(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Qg)
             ! Surface runoff:
             IF((out%Qs(land_x(i),land_y(i),1)<ranges%Qs(1)).OR. &
               (out%Qs(land_x(i),land_y(i),1)>ranges%Qs(2))) &
               CALL range_abort('Qs out of specified ranges!', ktau, met, &
               out%Qs(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Qs)
             ! Subsurface runoff:
             IF((out%Qsb(land_x(i),land_y(i),1)<ranges%Qsb(1)).OR. &
               (out%Qsb(land_x(i),land_y(i),1)>ranges%Qsb(2))) &
               CALL range_abort('Qsb out of specified ranges!', ktau, met, &
               out%Qsb(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Qsb)
             ! Total evapotranspiration:
             IF((out%Evap(land_x(i),land_y(i),1)<ranges%Evap(1)).OR. &
               (out%Evap(land_x(i),land_y(i),1)>ranges%Evap(2))) &
               CALL range_abort('Evap out of specified ranges!', ktau, met, &
               out%Evap(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Evap)
             ! Wet canopy evaporation:
             IF((out%ECanop(land_x(i),land_y(i),1)<ranges%ECanop(1)).OR. &
               (out%ECanop(land_x(i),land_y(i),1)>ranges%ECanop(2))) &
               CALL range_abort('ECanop out of specified ranges!', ktau, met, &
               out%ECanop(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%ECanop)
             ! Vegetation transpiration:
             IF((out%TVeg(land_x(i),land_y(i),1)<ranges%TVeg(1)).OR. &
               (out%TVeg(land_x(i),land_y(i),1)>ranges%TVeg(2))) &
               CALL range_abort('TVeg out of specified ranges!', ktau, met, &
               out%TVeg(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%TVeg)
             ! Soil evaporation:
             IF((out%ESoil(land_x(i),land_y(i),1)<ranges%ESoil(1)).OR. &
               (out%ESoil(land_x(i),land_y(i),1)>ranges%ESoil(2))) &
               CALL range_abort('ESoil out of specified ranges!', ktau, met, &
               out%ESoil(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%ESoil)
             ! Sensible heat from canopy:
             IF((out%HVeg(land_x(i),land_y(i),1)<ranges%HVeg(1)).OR. &
               (out%HVeg(land_x(i),land_y(i),1)>ranges%HVeg(2))) &
               CALL range_abort('HVeg out of specified ranges!', ktau, met, &
               out%HVeg(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%HVeg)
             ! Sensible heat from soil:
             IF((out%HSoil(land_x(i),land_y(i),1)<ranges%HSoil(1)).OR. &
               (out%HSoil(land_x(i),land_y(i),1)>ranges%HSoil(2))) &
               CALL range_abort('HSoil out of specified ranges!', ktau, met, &
               out%HSoil(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%HSoil)
             ! Net ecosystem exchange:
             IF((out%NEE(land_x(i),land_y(i),1)<ranges%NEE(1)).OR. &
               (out%NEE(land_x(i),land_y(i),1)>ranges%NEE(2))) &
               CALL range_abort('NEE out of specified ranges!', ktau, met, &
               out%NEE(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%NEE)
          END IF
       END DO
       WHERE(mask/=1) ! where not land:
          out%Qle(:,:,1) = -1000.0
          out%Qh(:,:,1) = -1000.0
          out%Qg(:,:,1) = -1000.0
          out%Qs(:,:,1) = -0.01
          out%Qsb(:,:,1) = -0.01
          out%Evap(:,:,1) = -0.01
          out%ECanop(:,:,1) = -0.01
          out%TVeg(:,:,1) = -0.01
          out%ESoil(:,:,1) = -0.01  
          out%HVeg(:,:,1) = -1000.0
          out%HSoil(:,:,1) = -1000.0
          out%NEE(:,:,1) = -100.0
       END WHERE
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%Qle,out%Qle,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qle variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Qh,out%Qh,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qh variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Qg,out%Qg,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
        status=NF90_PUT_VAR(ncid_out,ovid%Qs,out%Qs,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qs variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Qsb,out%Qsb,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Qsb variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Evap,out%Evap,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Evap variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%ECanop,out%ECanop,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing ECanop variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%TVeg,out%TVeg,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing TVeg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%ESoil,out%ESoil,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing ESoil variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%HVeg,out%HVeg,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HVeg variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%HSoil,out%HSoil,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HSoil variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%NEE,out%NEE,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing NEE variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')    
    END IF
    !--------------------------WRITE SOIL DATA-------------------------------------
    IF(output%soil) THEN
       ! Write to ALMA-style output soil variables:
       DO i = 1, mp
          out%SoilMoist(land_x(i),land_y(i),:,1) = ssoil%wb(i,:)*soil%zse*1000.0
          out%SoilTemp(land_x(i),land_y(i),:,1) = ssoil%tgg(i,:)
          out%BaresoilT(land_x(i),land_y(i),1) = ssoil%tgg(i,1)
          out%SWE(land_x(i),land_y(i),1) = ssoil%snowd(i)
          out%SnowT(land_x(i),land_y(i),1) = ssoil%tggsn(i,1)
          out%SnowDepth(land_x(i),land_y(i),1) = SUM(ssoil%sdepth(i,:))
           IF(check%ranges) THEN  ! Check ranges:
              DO j=1,ms
                 ! SoilMoisture:
                 IF((out%SoilMoist(land_x(i),land_y(i),j,1)<ranges%SoilMoist(1)).OR. &
                      (out%SoilMoist(land_x(i),land_y(i),j,1)>ranges%SoilMoist(2))) &
                      CALL range_abort('SoilMoist out of specified ranges!', ktau, met, &
                      out%SoilMoist(land_x(i),land_y(i),j,1), land_x(i), land_y(i), &
                      ranges%SoilMoist)
                 ! Soil Temperature:
                 IF((out%SoilTemp(land_x(i),land_y(i),j,1)<ranges%SoilTemp(1)).OR. &
                      (out%SoilTemp(land_x(i),land_y(i),j,1)>ranges%SoilTemp(2))) &
                      CALL range_abort('SoilTemp out of specified ranges!', ktau, met, &
                      out%SoilTemp(land_x(i),land_y(i),j,1), land_x(i), land_y(i), &
                      ranges%SoilTemp)
              END DO
             ! Bare soil temperature:
             IF((out%BaresoilT(land_x(i),land_y(i),1)<ranges%BaresoilT(1)).OR. &
               (out%BaresoilT(land_x(i),land_y(i),1)>ranges%BaresoilT(2))) &
               CALL range_abort('BaresoilT out of specified ranges!', ktau, met, &
               out%BaresoilT(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%BaresoilT)
             ! Snow water equivalent:
             IF((out%SWE(land_x(i),land_y(i),1)<ranges%SWE(1)).OR. &
               (out%SWE(land_x(i),land_y(i),1)>ranges%SWE(2))) &
               CALL range_abort('SWE out of specified ranges!', ktau, met, &
               out%SWE(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%SWE)
             ! Snow surface temperature:
             IF((out%SnowT(land_x(i),land_y(i),1)<ranges%SnowT(1)).OR. &
               (out%SnowT(land_x(i),land_y(i),1)>ranges%SnowT(2))) &
               CALL range_abort('SnowT out of specified ranges!', ktau, met, &
               out%SnowT(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%SnowT)
              ! Snow depth:
             IF((out%SnowDepth(land_x(i),land_y(i),1)<ranges%SnowDepth(1)).OR. &
               (out%SnowDepth(land_x(i),land_y(i),1)>ranges%SnowDepth(2))) &
               CALL range_abort('SnowDepth out of specified ranges!', ktau, met, &
               out%SnowDepth(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%SnowDepth)
          END IF
       END DO
        WHERE(mask/=1) ! where not land:
          out%SoilMoist(:,:,1,1) = -100.0
          out%SoilMoist(:,:,2,1) = -100.0
          out%SoilMoist(:,:,3,1) = -100.0
          out%SoilMoist(:,:,4,1) = -100.0
          out%SoilMoist(:,:,5,1) = -100.0
          out%SoilMoist(:,:,6,1) = -100.0
          out%SoilTemp(:,:,1,1) = 200.0
          out%SoilTemp(:,:,2,1) = 200.0
          out%SoilTemp(:,:,3,1) = 200.0
          out%SoilTemp(:,:,4,1) = 200.0
          out%SoilTemp(:,:,5,1) = 200.0
          out%SoilTemp(:,:,6,1) = 200.0
          out%BaresoilT(:,:,1) = 200.0
          out%SWE(:,:,1) = -100.0
          out%SnowT(:,:,1) = 300.0
          out%SnowDepth(:,:,1) = -0.1
       END WHERE
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%SoilMoist,out%SoilMoist,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,ms,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SoilMoist variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%SoilTemp,out%SoilTemp,start=(/1,1,1,ktau/), &
            count=(/xdimsize,ydimsize,ms,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SoilTemp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%BaresoilT,out%BaresoilT,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing BaresoilT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)') 
       status=NF90_PUT_VAR(ncid_out,ovid%SWE,out%SWE,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWE variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)') 
       status=NF90_PUT_VAR(ncid_out,ovid%SnowT,out%SnowT,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SnowT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)') 
       status=NF90_PUT_VAR(ncid_out,ovid%SnowDepth,out%SnowDepth,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SnowDepth variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)') 
    END IF
    !----------------------------WRITE RADIATION DATA----------------------------------
    IF(output%radiation) THEN
       ! Write to ALMA-style output radiation variables:
       DO i = 1, mp
          out%SWnet(land_x(i),land_y(i),1) = SUM(rad%qcan(i,:,1)) + &
               SUM(rad%qcan(i,:,2)) + rad%qssabs(i)
          out%LWnet(land_x(i),land_y(i),1) = met%fld(1) - sboltz*emleaf* &
               canopy%tv(1)**4 * (1-rad%transd(1)) - rad%flws(1)*rad%transd(1)
          out%Albedo(land_x(i),land_y(i),1) = (rad%albedo(1,1) + rad%albedo(1,2))*0.5
          out%RadT(land_x(i),land_y(i),1) = (1.0-rad%transd(i))*emleaf*sboltz* &
               canopy%tv(i)**4 + rad%transd(i)*emsoil*sboltz* ( &
               (1-ssoil%isflag(i))*ssoil%tgg(i,1)+ssoil%isflag(i)*ssoil%tggsn(i,1) &
               ) ** 4
          out%RadT(land_x(i),land_y(i),1) = (out%RadT(land_x(i),land_y(i),1)/sboltz)**0.25
           IF(check%ranges) THEN  ! Check ranges:
             ! Net shortwave:
             IF((out%SWnet(land_x(i),land_y(i),1)<ranges%SWnet(1)).OR. &
               (out%SWnet(land_x(i),land_y(i),1)>ranges%SWnet(2))) &
               CALL range_abort('SWnet out of specified ranges!', ktau, met, &
               out%SWnet(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%SWnet)
             ! Net longwave:
             IF((out%LWnet(land_x(i),land_y(i),1)<ranges%LWnet(1)).OR. &
               (out%LWnet(land_x(i),land_y(i),1)>ranges%LWnet(2))) &
               CALL range_abort('LWnet out of specified ranges!', ktau, met, &
               out%LWnet(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%LWnet)
             ! Albedo:
             IF((out%Albedo(land_x(i),land_y(i),1)<ranges%Albedo(1)).OR. &
               (out%Albedo(land_x(i),land_y(i),1)>ranges%Albedo(2))) &
               CALL range_abort('Albedo out of specified ranges!', ktau, met, &
               out%Albedo(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%Albedo)
             ! Surface radiative temperature:
             IF((out%RadT(land_x(i),land_y(i),1)<ranges%RadT(1)).OR. &
               (out%RadT(land_x(i),land_y(i),1)>ranges%RadT(2))) &
               CALL range_abort('RadT out of specified ranges!', ktau, met, &
               out%RadT(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%RadT)
          END IF
       END DO
       WHERE(mask/=1) ! where not land:
          out%SWnet(:,:,1) = -100.0
          out%LWnet(:,:,1) = 100.0
          out%Albedo(:,:,1) = -0.1
          out%RadT(:,:,1) = 200.0
       END WHERE     
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%SWnet,out%SWnet,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing SWnet variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%LWnet,out%LWnet,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LWnet variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Albedo,out%Albedo,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing Albedo variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%RadT,out%RadT,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing RadT variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
    END IF
    !-----------------------WRITE VEGETATION DATA----------------------------------
    IF(output%veg) THEN
       ! Write to ALMA-style output vegetation variables:
       DO i = 1, mp
          out%VegT(land_x(i),land_y(i),1) = canopy%tv(i)
          out%CanopInt(land_x(i),land_y(i),1) = canopy%cansto(i)
          out%LAI(land_x(i),land_y(i),1) = veg%vlai(i)
          IF(check%ranges) THEN  ! Check ranges:
             ! Vegetation temperature:
             IF((out%VegT(land_x(i),land_y(i),1)<ranges%VegT(1)).OR. &
               (out%VegT(land_x(i),land_y(i),1)>ranges%VegT(2))) &
               CALL range_abort('VegT out of specified ranges!', ktau, met, &
               out%VegT(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%VegT)
             IF((out%CanopInt(land_x(i),land_y(i),1)<ranges%CanopInt(1)).OR. &
               (out%CanopInt(land_x(i),land_y(i),1)>ranges%CanopInt(2))) &
               CALL range_abort('CanopInt out of specified ranges!', ktau, met, &
               out%CanopInt(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%CanopInt)
             ! LAI is an input, so ranges are not checked here.
          END IF
       END DO
       WHERE(mask/=1) ! where not land:
          out%VegT(:,:,1) = 200.0
          out%CanopInt(:,:,1) = -0.1
          out%LAI(:,:,1) = -1.0
       END WHERE  
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%VegT,out%VegT,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing VegT variable to file ' &
            //TRIM(filename_out)// ' (SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%CanopInt,out%CanopInt,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing CanopInt variable to file ' &
            //TRIM(filename_out)// ' (SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%LAI,out%LAI,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing LAI variable to file ' &
            //TRIM(filename_out)// ' (SUBROUTINE write_output)')
    END IF
    !------------------------WRITE CONSERVATION DATA---------------------------------
    IF(output%balances) THEN
       ! Write energy and water balance variables:
       DO i = 1, mp
          out%Ebal(land_x(i),land_y(i),1) = bal%ebal_tot(i)
          out%Wbal(land_x(i),land_y(i),1) = bal%wbal_tot(i)
       END DO
       WHERE(mask/=1) ! where not land:
          out%Ebal(:,:,1) = 0.0
          out%Wbal(:,:,1) = 0.0
       END WHERE  
       ! Write variables to netcdf file:
       status=NF90_PUT_VAR(ncid_out,ovid%Ebal,out%Ebal,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing energy balance variable to file ' &
            //TRIM(filename_out)// ' (SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%Wbal,out%Wbal,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing water balance variable to file ' &
            //TRIM(filename_out)// ' (SUBROUTINE write_output)')
    END IF
    !------------------------WRITE CARBON DATA----------------------------------------------
    IF(output%carbon) THEN
       ! Write energy and water balance variables:
       DO i = 1, mp
          out%NEE(land_x(i),land_y(i),1) = canopy%fnee(i)/1.201E-5 ! g/m2/s to umol/m2/s 
          out%GPP(land_x(i),land_y(i),1) = -1.0*canopy%fpn(i)/1.201E-5 
          out%AutoResp(land_x(i),land_y(i),1) = canopy%frp(i)/1.201E-5 
          out%HeteroResp(land_x(i),land_y(i),1) = canopy%frs(i)/1.201E-5 
          out%NPP(land_x(i),land_y(i),1) = out%GPP(land_x(i),land_y(i),1) - &
               out%AutoResp(land_x(i),land_y(i),1)
          IF(check%ranges) THEN  ! Check ranges:
              IF((out%NEE(land_x(i),land_y(i),1)<ranges%NEE(1)).OR. &
               (out%NEE(land_x(i),land_y(i),1)>ranges%NEE(2))) &
               CALL range_abort('NEE out of specified ranges!', ktau, met, &
               out%NEE(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%NEE)
              IF((out%GPP(land_x(i),land_y(i),1)<ranges%GPP(1)).OR. &
               (out%GPP(land_x(i),land_y(i),1)>ranges%GPP(2))) &
               CALL range_abort('GPP out of specified ranges!', ktau, met, &
               out%GPP(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%GPP)
              IF((out%AutoResp(land_x(i),land_y(i),1)<ranges%AutoResp(1)).OR. &
               (out%AutoResp(land_x(i),land_y(i),1)>ranges%AutoResp(2))) &
               CALL range_abort('AutoResp out of specified ranges!', ktau, met, &
               out%AutoResp(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%AutoResp)
              IF((out%HeteroResp(land_x(i),land_y(i),1)<ranges%HeteroResp(1)).OR. &
               (out%HeteroResp(land_x(i),land_y(i),1)>ranges%HeteroResp(2))) &
               CALL range_abort('HeteroResp out of specified ranges!', ktau, met, &
               out%HeteroResp(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%HeteroResp)
              IF((out%NPP(land_x(i),land_y(i),1)<ranges%NPP(1)).OR. &
               (out%NPP(land_x(i),land_y(i),1)>ranges%NPP(2))) &
               CALL range_abort('NPP out of specified ranges!', ktau, met, &
               out%NPP(land_x(i),land_y(i),1), land_x(i), land_y(i), &
               ranges%NPP)
           END IF
       END DO
       WHERE(mask/=1) ! where not land:
          out%NEE(:,:,1) = -100.0
          out%GPP(:,:,1) = -100.0
          out%NPP(:,:,1) = -100.0
          out%AutoResp(:,:,1) = -100.0
          out%HeteroResp(:,:,1) = -100.0
       END WHERE  
       ! Write variables to netcdf file:
       IF(.NOT.output%flux) THEN ! make sure NEE hasn't already been written to file:
          status=NF90_PUT_VAR(ncid_out,ovid%NEE,out%NEE,start=(/1,1,ktau/), &
               count=(/xdimsize,ydimsize,1/))
          IF(status/=NF90_NOERR) CALL nc_abort('Error writing NEE variable to file ' &
               //TRIM(filename_out)// '(SUBROUTINE write_output)')
       END IF
       status=NF90_PUT_VAR(ncid_out,ovid%GPP,out%GPP,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing GPP variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%AutoResp,out%AutoResp,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing AutoResp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%HeteroResp,out%HeteroResp,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
       IF(status/=NF90_NOERR) CALL nc_abort('Error writing HeteroResp variable to file ' &
            //TRIM(filename_out)// '(SUBROUTINE write_output)')
       status=NF90_PUT_VAR(ncid_out,ovid%NPP,out%NPP,start=(/1,1,ktau/), &
            count=(/xdimsize,ydimsize,1/))
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
    IF(output%radiation) DEALLOCATE(out%SWnet,out%LWnet,out%Albedo,out%RadT) 
    IF(output%veg) DEALLOCATE(out%VegT,out%LAI,out%CanopInt) 
    IF(output%carbon) THEN
       DEALLOCATE(out%NPP,out%GPP,out%AutoResp,out%HeteroResp)
       IF(.NOT.output%flux) DEALLOCATE(out%NEE)
    END IF
    IF(output%balances) THEN
       DEALLOCATE(out%Ebal,out%Wbal)
       IF(mp<5) THEN
          WRITE(logn,*)
          WRITE(logn,*) 'Cumulative energy balance: ',bal%ebal_tot,'W/m^2'
          WRITE(logn,*) 'Cumulative water balance: ', bal%wbal_tot,'mm'
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
  !========================================================================
  SUBROUTINE range_abort(message,ktau,met,value,xx,yy,var_range)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER(i_d), INTENT(IN) :: ktau ! time step
    TYPE(met_type),INTENT(IN) :: met  ! met data
    REAL(r_1),INTENT(IN) :: value ! value deemed to be out of range
    INTEGER(i_d),INTENT(IN) :: xx
    INTEGER(i_d),INTENT(IN) ::yy ! coordinates of erroneous grid square
    REAL(r_1),DIMENSION(2),INTENT(IN) :: var_range ! appropriate var range 
    WRITE(*,*) message ! error from subroutine
    WRITE(*,44) 'Timestep',ktau,', or ', met%hod,'hod, ', &
         INT(met%doy),'doy, ',INT(met%year)
44  FORMAT(1X,A8,I8,A5,F5.2,1X,A5,I3,1X,A5,I4)
    WRITE(*,*) 'Site lat, lon:',lat_all(xx,yy),lon_all(xx,yy)
    WRITE(*,*) 'Specified acceptable range (checks.f90):', var_range(1), &
         'to',var_range(2)
    WRITE(*,*) 'Value:',value 
    STOP
  END SUBROUTINE range_abort

END MODULE output_module
