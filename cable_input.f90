!!$ cable_input.f90
!!$
!!$ Input module for CABLE land surface scheme offline driver; 
!!$
!!$ Gab Abramowitz 2007 University of New South Wales/
!!$ CSIRO Marine and Atmospheric Research; gabsun@gmail.com
!!$
!!$ The subroutines in this module are:
!!$ 
!!$ get_default_lai - reads LAI from default coarse grid netcdf file every time step
!!$ open_met_file - opens netcdf met forcing file and checks for variables
!!$ get_met_data - reads met (and LAI if present) from met file each time step
!!$ close_met_file - closes the netcdf met forcing file
!!$ get_restart_data - reads initialisations and parameters from restart file
!!$ load_parameters - decides where CABLE's parameters and init should load from
!!$ get_parameters_met - looks for and loads any of CABLE's parameters from met file
!!$ nc_abort - an abort subroutine that prints netcdf error message

MODULE input_module   
  USE checks_module
  USE abort_module
  USE canopy_module
  USE physical_constants
  USE parameter_module
  USE netcdf ! link must be made in cd to netcdf-x.x.x/src/f90/netcdf.mod
  IMPLICIT NONE
  INTEGER(i_d) :: ncid_met ! met data netcdf file ID
  INTEGER(i_d) :: ncid_par ! parameter data netcdf file ID
  INTEGER(i_d) :: ncid_lai ! lai data netcdf file ID
  INTEGER(i_d) :: ncid_rin ! input netcdf restart file ID
  INTEGER(i_d) :: status   ! netcdf error status
  INTEGER(i_d) :: logn     ! log file unit number
  CHARACTER(LEN=4) :: gridType ! Either 'land' or 'mask' - see ALMA compress by gathering
  INTEGER(i_d),POINTER,DIMENSION(:) :: landGrid ! for ALMA compressed variables
  INTEGER(i_d),POINTER,DIMENSION(:,:) :: mask ! land/sea mask from met file
  REAL(r_1),POINTER, DIMENSION(:,:) :: lat_all, lon_all ! lat and lon
  INTEGER(i_d),POINTER,DIMENSION(:) :: land_x,land_y ! indicies of land in mask
  REAL(r_1),POINTER,DIMENSION(:,:)  :: elevation ! site/grid cell elevation
  INTEGER(i_d) :: xdimsize,ydimsize ! sizes of x and y dimensions
  INTEGER(i_d) :: ngridcells ! number of gridcells in simulation
  LOGICAL :: leaps   ! use leap year timing?
  LOGICAL :: verbose ! print init and param details of all grid cells?
  CHARACTER(LEN=3) :: time_coord ! GMT or LOCal time variables
  CHARACTER(LEN=33) :: timeunits ! timing info read from nc file
  REAL(r_2),POINTER,DIMENSION(:) :: timevar ! time variable from file
  REAL(r_1)    :: shod ! start time hour-of-day
  INTEGER(i_d) :: sdoy,smoy,syear ! start time day-of-year month and year
  TYPE met_varID_type 
     INTEGER(i_d) :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair,Qair,Rainf, &
          Snowf,CO2air,Elev,LAI,iveg,isoil
  END TYPE met_varID_type
  TYPE(met_varID_type) :: id ! netcdf variable IDs for input met variables
  TYPE met_units_type
     CHARACTER(LEN=20) :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair, &
          Qair,Rainf,Snowf,CO2air,Elev
  END TYPE met_units_type
  TYPE(met_units_type) :: metunits ! units for meteorological variables
  TYPE convert_units_type
     REAL(r_1) :: PSurf,Tair,Qair,Rainf,CO2air,Elev
  END TYPE convert_units_type
  TYPE(convert_units_type) :: convert ! units change factors for met variables
  TYPE input_details_type
     LOGICAL :: Wind ! T => 'Wind' is present; F => use vector component wind
     LOGICAL :: LWdown ! T=> downward longwave is present in met file
     LOGICAL :: CO2air ! T=> air CO2 concentration is present in met file
     LOGICAL :: PSurf ! T=> surface air pressure is present in met file
     LOGICAL :: Snowf ! T=> snowfall variable is present in met file
     LOGICAL :: LAI   ! T=> LAI is present in the met file
     LOGICAL :: parameters ! TRUE if non-default parameters are found
     LOGICAL :: initial ! switched to TRUE when initialisation data are loaded
  END TYPE input_details_type
  TYPE(input_details_type) :: exists
  REAL(r_1) :: fixedCO2 ! CO2 level if CO2air not in met file
  TYPE checks_type
     LOGICAL :: ranges, energy_bal, mass_bal
  END TYPE checks_type
  TYPE(checks_type) :: check ! what types of checks to perform
  TYPE parID_type ! model parameter IDs in netcdf file
     INTEGER :: bch,latitude,clay,css,rhosoil,hyds,rs20,sand,sfc,silt, &
         ssat,sucs,swilt,froot,zse,canst1,dleaf,meth,za, &
         ejmax,frac4,hc,lai,rp20,rpcoef,shelrb, vbeta, &
         vcmax,xfang,ratecp,ratecs,refsbare,isoil,iveg,albsoil,&
         taul,refl,tauw,refw,tminvj,tmaxvj,veg_class,soil_class
  END TYPE parID_type
CONTAINS

! =================================== LAI ====================================
  SUBROUTINE get_default_lai(ktau,doy,veg,kend) 
    ! Fetches default LAI data from netcdf file, based on day-of-year (doy).
    INTEGER, INTENT(IN) :: ktau ! total run timestep
    ! gridpoint numbers for sites:
    REAL(r_1),DIMENSION(*),INTENT(IN) :: doy ! day of year
    TYPE(veg_parameter_type),INTENT(INOUT) :: veg ! LAI retrieved from file
    REAL(r_1),DIMENSION(1,1) :: templai2 ! LAI retrieved from file
    INTEGER(i_d), INTENT(IN) :: kend ! total number of timesteps in run
    INTEGER(i_d) :: month ! Current timestep's month
    INTEGER(i_d),SAVE :: ncid_lai ! LAI file ID
    INTEGER(i_d) :: laiID ! LAI variable ID
    INTEGER(i_d) :: a ! do loop counter
    
    IF(ktau==1) THEN
       ! Open netcdf file
       status = NF90_OPEN('surface_data/lai48.nc',0,ncid_lai) 
       IF (status /= NF90_NOERR) CALL nc_abort('Error opening LAI file.')
    END IF
   
    ! Allocate space for LAI variable:
    DO a=1,SIZE(gdpt)
       ! Determine month; don't bother about distinction between leap and 
       ! non-leap for LAI purposes:
       SELECT CASE(INT(doy(a)))
       CASE(0:31)
          month = 1
       CASE(32:59)
          month = 2
       CASE(60:90)
          month = 3
       CASE(91:120)
          month = 4
       CASE(121:151)
          month = 5
       CASE(152:181)
          month = 6
       CASE(182:212)
          month = 7
       CASE(213:243)
          month = 8
       CASE(244:273)
          month = 9
       CASE(274:304)
          month = 10
       CASE(305:334)
          month = 11
       CASE(335:366)
          month = 12
       CASE DEFAULT
          print*, 'Day of year at site #',a,', is ',doy(a)
          CALL abort('Unknown doy!')
       END SELECT
       ! Read LAI value:
       status = NF90_INQ_VARID(ncid_lai,'lai',laiID)
       IF (status /= NF90_NOERR)CALL nc_abort ('Error finding LAI variable.')
       status = NF90_GET_VAR(ncid_lai,laiID,templai2, &
            (/gdpt(a),month/),(/1,1/))
       veg%vlai(a) = templai2(1,1)
       veg%vlaimax = 1.1 * veg%vlai
    END DO
  
    ! Close netcdf file
    IF(ktau==kend) CLOSE(ncid_lai)

    ! Check for unrealistic LAI values:
    IF(ANY(veg%vlai<0.0).OR.ANY(veg%vlai>15.0)) THEN
       WRITE(*,*) 'Timestep: ',ktau
       CALL abort('Bad LAI value! (outside [0,15.0])')
    END IF
  END SUBROUTINE get_default_lai
!================================================================================
  SUBROUTINE open_met_file(filename_met,dels,kend)
    ! Opens netcdf file containing meteorological (LSM input) data
    ! and determines:
    ! 1. Spatial details - number of sites/grid cells, latitudes, longitudes
    ! 2. Timing details - time step size, number of timesteps, starting date,
    !    and whether time coordinate is local or GMT
    ! 3. Checks availability, including units issues, of all required
    !    meteorological input variables. Also checks whether or not LAI is 
    !    present, and fetches prescribed veg ans soil type if present.
    CHARACTER(LEN=*), INTENT(IN) :: filename_met ! name of file for met data
    REAL(r_1), INTENT(OUT) :: dels ! time step size
    INTEGER(i_d), INTENT(OUT) :: kend ! number of time steps in simulation
    INTEGER(i_d) :: timevarID ! time variable ID number
    INTEGER(i_d),DIMENSION(1) :: timedimID ! time dimension ID number
    INTEGER(i_d) :: xdimID,ydimID    ! x and y dimension ID numbers
    INTEGER(i_d) :: maskID    ! mask variable ID
    INTEGER(i_d) :: landID    ! land variable ID
    INTEGER(i_d) :: landdimID ! land dimension ID
    INTEGER(i_d) :: latitudeID, longitudeID ! lat and lon variable IDs
    REAL(r_1),POINTER, DIMENSION(:) :: lat_temp, lon_temp ! lat and lon
    INTEGER,POINTER,DIMENSION(:) ::land_xtmp,land_ytmp ! temp indicies
    REAL(r_1)    :: tshod        ! temporary variable
    INTEGER(i_d) :: tsdoy,tsyear ! temporary variables
    REAL(r_1)    :: ehod ! end time hour-of-day
    INTEGER(i_d) :: edoy,eyear ! end time day-of-year and year
    INTEGER(i_d) :: jump_days ! days made by first "time" entry
    INTEGER(i_d) :: sdoytmp ! used to determine start time hour-of-day
    INTEGER(i_d) :: mp_ctr ! counter for number of land points read from file
    INTEGER(i_d) :: mp_fromfile ! number of land points in file 
    LOGICAL :: all_met ! ALL required met in met file (no synthesis)?
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp log file
    REAL(r_1),DIMENSION(1,1) :: data2 ! temp variable for netcdf reading
    INTEGER(i_d) :: x,y,i,j ! do loop counters
    
    ! Initialise parameter loading switch - will be set to TRUE when 
    ! parameters are loaded:
    exists%parameters = .FALSE. ! initialise
    ! Initialise initialisation loading switch - will be set to TRUE when 
    ! initialisation data are loaded:
    exists%initial = .FALSE. ! initialise

    ! Write filename to log file:
    WRITE(logn,*) '============================================================'
    WRITE(logn,*) 'Log file for offline CABLE run:'
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    WRITE(logn,*) TRIM(nowtime),' ',TRIM(todaydate)
    WRITE(logn,*) '============================================================'
    WRITE(logn,*) 'Opening met data file: ', TRIM(filename_met)

    ! Open netcdf file:
    status = NF90_OPEN(filename_met,0,ncid_met) ! open met data file
    IF (status /= NF90_NOERR) CALL nc_abort &
         ('Error opening netcdf file '//TRIM(filename_met)// &
         ' (SUBROUTINE open_met_file)') 
    
    !!=====================VV Determine spatial details VV=================
    ! Determine number of sites/gridcells.
    ! Find size of 'x' or 'lat' dimension:
    status = NF90_INQ_DIMID(ncid_met,'x', xdimID)
    IF(status/=NF90_NOERR) THEN ! if failed
       ! Try 'lat' instead of x
       status = NF90_INQ_DIMID(ncid_met,'lat', xdimID)
       IF(status/=NF90_NOERR) CALL nc_abort &
            ('Error finding x dimension in '&
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    END IF
    status = NF90_INQUIRE_DIMENSION(ncid_met,xdimID,len=xdimsize)
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error determining size of x dimension in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Find size of 'y' dimension:
    status = NF90_INQ_DIMID(ncid_met,'y', ydimID)
    IF(status/=NF90_NOERR) THEN ! if failed
       ! Try 'lon' instead of y
       status = NF90_INQ_DIMID(ncid_met,'lon', ydimID)
       IF(status/=NF90_NOERR) CALL nc_abort &
            ('Error finding y dimension in ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    END IF
    status = NF90_INQUIRE_DIMENSION(ncid_met,ydimID,len=ydimsize)
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error determining size of y dimension in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Determine number of gridcells in netcdf file:
    ngridcells = xdimsize*ydimsize
    WRITE(logn,'(A28,I7)') 'Total number of gridcells: ', ngridcells
    
    ! Get all latitude and longitude values.
    ! Find latitude variable (try 'latitude' and 'nav_lat'(ALMA)):
    status = NF90_INQ_VARID(ncid_met, 'latitude', latitudeID)
    IF(status /= NF90_NOERR) THEN
       status = NF90_INQ_VARID(ncid_met, 'nav_lat', latitudeID)
       IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding latitude variable in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    END IF
    ! Allocate space for lat_all variable:
    ALLOCATE(lat_all(xdimsize,ydimsize))
    ! Get latitude values for entire region:
    status= NF90_GET_VAR(ncid_met,latitudeID,lat_all)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error reading latitude variable in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Find longitude variable (try 'longitude' and 'nav_lon'(ALMA)):
    status = NF90_INQ_VARID(ncid_met, 'longitude', longitudeID)
    IF(status /= NF90_NOERR) THEN
       status = NF90_INQ_VARID(ncid_met, 'nav_lon', longitudeID)
       IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding longitude variable in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    END IF
    ! Allocate space for lon_all variable:
    ALLOCATE(lon_all(xdimsize,ydimsize))
    ! Get longitude values for entire region:
    status= NF90_GET_VAR(ncid_met,longitudeID,lon_all)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error reading longitude variable in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')

    ! Check for "mask" variable or "land" variable to tell grid type
    ! (and allow neither if only one gridpoint). "mask" is a 2D variable
    ! with dims x,y and "land" is a 1D variable.
    status = NF90_INQ_VARID(ncid_met, 'mask', maskID) ! check for "mask"
    IF(status /= NF90_NOERR) THEN ! if error, i.e. no "mask" variable:
       ! Check for "land" variable:
       status = NF90_INQ_VARID(ncid_met, 'land', landID)
       IF(status /= NF90_NOERR) THEN ! ie no "land" or "mask"
          IF(ngridcells==1) THEN 
             ! Allow no explicit grid system if only one gridpoint
             ALLOCATE(mask(xdimsize,ydimsize)) ! Allocate "mask" variable
             gridType='mask' ! Use mask system, one gridpoint.
             mask = 1
             ALLOCATE(latitude(1),longitude(1))
             latitude = lat_all(1,1)
             longitude = lon_all(1,1)
             mp_fromfile=1
             ALLOCATE(land_x(mp_fromfile),land_y(mp_fromfile))
             land_x = 1
             land_y = 1
          ELSE
             ! Call abort if more than one gridcell and no
             ! recognised grid system:
             CALL nc_abort &
                  ('Error finding grid system ("mask" or "land") variable in ' &
                  //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
          END IF
       ELSE ! i.e. "land" variable exists
          gridType='land'
          ! Check size of "land" dimension:
          status = NF90_INQ_DIMID(ncid_met,'land', landdimID)
          IF(status/=NF90_NOERR) CALL nc_abort &
               ('Error finding land dimension in ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
          status = NF90_INQUIRE_DIMENSION(ncid_met,landdimID,len=mp_fromfile)
          IF(status/=NF90_NOERR) CALL nc_abort &
               ('Error determining size of land dimension in ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
          ! Allocate landGrid variable:
          ALLOCATE(landGrid(mp_fromfile))
          ! Get values of "land" variable from file:
          status= NF90_GET_VAR(ncid_met,landID,landGrid)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading "land" variable in ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
          ! Allocate latitude and longitude variables:
          ALLOCATE(latitude(mp_fromfile),longitude(mp_fromfile))
          ! Write to indicies of points in all-grid which are land
          ALLOCATE(land_x(mp_fromfile),land_y(mp_fromfile))
          ! Allocate "mask" variable:
          ALLOCATE(mask(xdimsize,ydimsize))
          ! Initialise all gridpoints as sea:
          mask = 0
          DO j=1, mp_fromfile ! over all land points
             ! Find x and y coords of current land point
             y = INT((landGrid(j)-1)/xdimsize)
             x = landGrid(j) - y * xdimsize
             y=y+1
             ! Write lat and lon to land-only lat/lon vars:
             latitude(j) = lat_all(x,y)
             longitude(j) = lon_all(x,y)
             ! Write to mask variable:
             mask(x,y)=1
             ! Save indicies:
             land_x(j) = x
             land_y(j) = y
          END DO
       END IF ! does "land" variable exist 
    ELSE ! i.e. "mask" variable exists
       ! Allocate "mask" variable:
       ALLOCATE(mask(xdimsize,ydimsize))
       gridType='mask' ! Use mask system
       ! Get mask values from file:
       status= NF90_GET_VAR(ncid_met,maskID,mask)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error reading "mask" variable in ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       ! Allocate space for extracting land lat/lon values:
       ALLOCATE(lat_temp(ngridcells),lon_temp(ngridcells))
       ! Allocate space for extracting index of mask which is land
       ALLOCATE(land_xtmp(ngridcells),land_ytmp(ngridcells))
       ! Cycle through all gridsquares:
       mp_ctr = 0 ! initialise
       DO y=1,ydimsize
          DO x=1,xdimsize
             IF(mask(x,y)==1) THEN ! If land
                mp_ctr = mp_ctr + 1
                ! Store lat and lon for land points
                lat_temp(mp_ctr) = lat_all(x,y)
                lon_temp(mp_ctr) = lon_all(x,y)
                ! Store indicies of points in mask which are land
                land_xtmp(mp_ctr) = x
                land_ytmp(mp_ctr) = y
             END IF
          END DO
       END DO
       ! Record number of land points
       mp_fromfile = mp_ctr
       ! Allocate latitude and longitude variables:
       ALLOCATE(latitude(mp_fromfile),longitude(mp_fromfile))
       ! Write to latitude and longitude variables:
       latitude = lat_temp(1:mp_fromfile)
       longitude = lon_temp(1:mp_fromfile)
       ! Write to indicies of points in mask which are land
       ALLOCATE(land_x(mp_fromfile),land_y(mp_fromfile))
       land_x = land_xtmp(1:mp_fromfile)
       land_y = land_ytmp(1:mp_fromfile)
       ! Clear lon_temp, lat_temp,land_xtmp,land_ytmp
       DEALLOCATE(lat_temp,lon_temp,land_xtmp,land_ytmp)
    END IF ! "mask" variable or no "mask" variable

    ! Write to number of land points to log file:
    WRITE(logn,'(24X,I7,A29)') mp_fromfile, ' of which are land grid cells'
 
    ! Set longitudes to be [-180,180]:
    WHERE(longitude>180.0) 
       longitude = longitude - 360.0
    END WHERE
    ! Check ranges for latitude and longitude:
    IF(ANY(longitude>180.0).OR.ANY(longitude<-180.0)) &
         CALL abort('Longitudes read from '//TRIM(filename_met)// &
         ' are not [-180,180] or [0,360]! Please set.')
    IF(ANY(latitude>90.0).OR.ANY(latitude<-90.0)) &
         CALL abort('Latitudes read from '//TRIM(filename_met)// &
         ' are not [-90,90]! Please set.')

    ! Set global mp value (number of land points), used to allocate
    ! all of CABLE's arrays:
    mp = mp_fromfile
    
    !!=================^^ End spatial details ^^========================

    !!=========VV Determine simulation timing details VV================
    ! Inquire 'time' variable's ID:
    status = NF90_INQ_VARID(ncid_met, 'time', timevarID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding time variable in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Get ID for dimension upon which time depends:
    status = NF90_INQUIRE_VARIABLE(ncid_met,timevarID,dimids=timedimID)
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error determining "time" dimension dimension in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Determine number of time steps:
    status = NF90_INQUIRE_DIMENSION(ncid_met,timedimID(1),len=kend)
    IF(status/=NF90_NOERR) CALL nc_abort &
         ('Error determining number of timesteps in ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Allocate time variable:
    ALLOCATE(timevar(kend))
    ! Fetch 'time' variable:
    status= NF90_GET_VAR(ncid_met,timevarID,timevar)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error reading time variable in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Set time step size:
    dels = REAL(timevar(2) - timevar(1),r_1)
    WRITE(logn,'(1X,A29,I8,A3,F10.3,A5)') 'Number of time steps in run: ',&
         kend,' = ', REAL(kend)/(3600/dels*24),' days'
    ! Write time step size to log file:
    WRITE(logn,'(1X,A17,F8.1,1X,A7)') 'Time step size:  ', dels, 'seconds'
    ! Get units for 'time' variable:
    status = NF90_GET_ATT(ncid_met,timevarID,'units',timeunits)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding time variable units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    WRITE(logn,*) 'Time variable units: ', timeunits
    ! Get coordinate field:
    status = NF90_GET_ATT(ncid_met,timevarID,'coordinate',time_coord)
    ! If error getting coordinate field (i.e. it doesn't exist):
    IF(status /= NF90_NOERR) THEN
       ! Assume default time coordinate:
       IF(mp_fromfile==1) THEN ! If single site, this is local time
          time_coord = 'LOC' ! 12am is 12am local time, at site/gridcell
       ELSE ! If multiple/global/regional, use GMT
          time_coord = 'GMT' ! 12am is GMT time, local time set by longitude
       END IF
    ELSE IF((status==NF90_NOERR.AND.time_coord=='LOC'.AND.mp_fromfile>1)) THEN
       ! Else if local time is selected for regional simulation, abort:
       CALL abort('"time" variable must be GMT for multiple site simulation!' &
            //' Check "coordinate" field in time variable.' &
            //' (SUBROUTINE open_met_file)')
    ELSE IF(time_coord/='LOC'.AND.time_coord/='GMT') THEN
       CALL abort('Meaningless time coordinate in met data file!' &
            // ' (SUBROUTINE open_met_file)')
    END IF
    
    ! Use internal files to convert "time" variable units (giving the run's 
    ! start time) from character to integer; calculate starting hour-of-day,
    ! day-of-year, year:
    READ(timeunits(15:18),*) syear
    READ(timeunits(20:21),*) smoy ! integer month
    READ(timeunits(23:24),*) sdoytmp ! integer day of that month
    READ(timeunits(26:27),*) shod  ! starting hour of day 
    ! Decide day-of-year for non-leap year:
    SELECT CASE(smoy)
    CASE(1) ! Jan
       sdoy=sdoytmp
    CASE(2) ! Feb
       sdoy=sdoytmp+31
    CASE(3) ! Mar
       sdoy=sdoytmp+59
    CASE(4)
       sdoy=sdoytmp+90
    CASE(5)
       sdoy=sdoytmp+120
    CASE(6)
       sdoy=sdoytmp+151
    CASE(7)
       sdoy=sdoytmp+181
    CASE(8)
       sdoy=sdoytmp+212
    CASE(9)
       sdoy=sdoytmp+243
    CASE(10)
       sdoy=sdoytmp+273
    CASE(11)
       sdoy=sdoytmp+304
    CASE(12) 
       sdoy=sdoytmp+334
    CASE DEFAULT
       CALL abort('Could not interpret month in "time" units from ' &
            //TRIM(filename_met)// '(SUBROUTINE open_met_file)')
    END SELECT
    ! If we're using leap year timing:
    IF(leaps) THEN
       ! If start year is a leap year and start month > Feb, add a day:
       IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & ! \
            (MOD(syear,4)==0.AND.MOD(syear,400)==0)) &   ! - leap year
            .AND.smoy>2) sdoy = sdoy + 1
       ! Number of days between start position and 1st timestep:
       jump_days = INT((timevar(1)/3600.0 + shod)/24.0)
       ! Cycle through days to find leap year inclusive starting date:
       DO i=1,jump_days
          sdoy = sdoy + 1
          ! IF leap year: 
          IF((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)) THEN
             ! Set moy:
             SELECT CASE(sdoy)
             CASE(1) ! Jan
                smoy = 1
             CASE(32) ! Feb
                smoy = 2
             CASE(61) ! Mar
                smoy = 3
             CASE(92)
                smoy = 4
             CASE(122)
                smoy = 5
             CASE(153)
                smoy = 6
             CASE(183)
                smoy = 7
             CASE(214)
                smoy = 8
             CASE(245)
                smoy = 9
             CASE(275)
                smoy = 10
             CASE(306)
                smoy = 11
             CASE(336) 
                smoy = 12
             CASE(367)! end of year; increment
                syear = syear + 1 
                smoy = 1
                sdoy = 1
             END SELECT
          ! ELSE IF not leap year:
          ELSE 
             ! Set moy:
             SELECT CASE(sdoy)
             CASE(1) ! Jan
                smoy = 1
             CASE(32) ! Feb
                smoy = 2
             CASE(60) ! Mar
                smoy = 3
             CASE(91)
                smoy = 4
             CASE(121)
                smoy = 5
             CASE(152)
                smoy = 6
             CASE(182)
                smoy = 7
             CASE(213)
                smoy = 8
             CASE(244)
                smoy = 9
             CASE(274)
                smoy = 10
             CASE(305)
                smoy = 11
             CASE(335) 
                smoy = 12
             CASE(366) ! end of year; increment
                syear = syear + 1 
                smoy = 1
                sdoy = 1
             END SELECT
          END IF
       END DO
       ! Update starting hour-of-day fot first time step's value
       shod = MOD(REAL(timevar(1)/3600.0 + shod),24.0)
    ELSE ! If not leap year timing,
       ! simply update starting times for first value of "time":
       tshod = MOD(REAL(timevar(1)/3600.0 + shod),24.0)
       tsdoy = MOD(INT((timevar(1)/3600.0 + shod)/24.0) + sdoy, 365)
       tsyear = INT(REAL(INT((timevar(1)/3600.0+shod)/24.0)+sdoy)/365.0)+syear
       shod=tshod  ! real valued
       sdoy=tsdoy  ! integer valued
       syear=tsyear ! integer valued
       ! Set moy:
       SELECT CASE(sdoy)
       CASE(1:31) ! Jan
          smoy = 1
       CASE(32:59) ! Feb
          smoy = 2
       CASE(60:90) ! Mar
          smoy = 3
       CASE(91:120)
          smoy = 4
       CASE(121:151)
          smoy = 5
       CASE(152:181)
          smoy = 6
       CASE(182:212)
          smoy = 7
       CASE(213:243)
          smoy = 8
       CASE(244:273)
          smoy = 9
       CASE(274:304)
          smoy = 10
       CASE(305:334)
          smoy = 11
       CASE(335:365) 
          smoy = 12
       END SELECT
    END IF
    ! Now all start time variables established, report to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') &
         'Run begins: ',shod,' hour-of-day, ',sdoy, ' day-of-year, ',&
         syear, time_coord, 'time'
    ! Determine ending time of run...
    IF(leaps) THEN ! If we're using leap year timing...
       ! Number of days between beginning and end of run:
       jump_days = INT(((timevar(kend)-timevar(1)+dels)/3600.0 + shod)/24.0)
       ! initialise:
       ehod = shod
       edoy = sdoy
       eyear = syear
       ! Cycle through days to find leap year inclusive ending date:
       DO i=1,jump_days
          edoy = edoy + 1
          ! IF leap year: 
          IF((MOD(eyear,4)==0.AND.MOD(eyear,100)/=0).OR. & 
               (MOD(eyear,4)==0.AND.MOD(eyear,400)==0)) THEN
             ! Set moy:
             SELECT CASE(edoy)
             CASE(367)! end of year; increment
                eyear = eyear + 1 
                edoy = 1
             END SELECT
          ! ELSE IF not leap year:
          ELSE 
             ! Set moy:
             SELECT CASE(edoy)
             CASE(366) ! end of year; increment
                eyear = eyear + 1 
                edoy = 1
             END SELECT
          END IF
       END DO
       ! Update starting hour-of-day fot first time step's value
       ehod = MOD(REAL((timevar(kend)-timevar(1)+dels)/3600.0 + shod),24.0)
    ELSE ! if not using leap year timing
       ! Update shod, sdoy, syear for first "time" value:
       ehod = MOD(REAL((timevar(kend)-timevar(1)+dels)/3600.0 + shod),24.0)
       edoy = MOD(INT(((timevar(kend)-timevar(1)+dels)/3600.0 + shod)/24.0) &
            + sdoy, 365)
       eyear = INT(REAL(INT(((timevar(kend)-timevar(1)+dels) &
            /3600.0+shod)/24.0)+sdoy)/365.0)+syear
    END IF
    ! Report finishing time to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') 'Run ends:   ',&
         ehod,' hour-of-day, ',edoy, &
         ' day-of-year, ', eyear, time_coord, 'time'
    !!===================^^ End timing details ^^==========================

    !!===================VV Look for met variables VV======================
    all_met = .TRUE. ! initialise
    ! Look for SWdown (essential):- - - - - - - - - - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'SWdown',id%SWdown)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding SWdown in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Get SWdown units and check okay:
    status = NF90_GET_ATT(ncid_met,id%SWdown,'units',metunits%SWdown)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding SWdown units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    IF(metunits%SWdown(1:4)/='W/m2'.AND.metunits%SWdown(1:5) &
         /='W/m^2'.AND.metunits%SWdown(1:5)/='Wm^-2' &
         .AND.metunits%SWdown(1:4)/='Wm-2') THEN
       WRITE(*,*) metunits%SWdown
       CALL abort('Unknown units for SWdown'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Tair (essential):- - - - - - - - - - - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'Tair',id%Tair)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Tair in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Get Tair units and check okay:
    status = NF90_GET_ATT(ncid_met,id%Tair,'units',metunits%Tair)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Tair units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Tair(1:1)=='C'.OR.metunits%Tair(1:1)=='c') THEN
       ! Change from celsius to kelvin:
       convert%Tair = tfrz
       WRITE(logn,*) 'Temperature will be converted from C to K'
    ELSE IF(metunits%Tair(1:1)=='K'.OR.metunits%Tair(1:1)=='k') THEN
       ! Units are correct
       convert%Tair = 0.0
    ELSE
       WRITE(*,*) metunits%Tair
       CALL abort('Unknown units for Tair'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Qair (essential):- - - - - - - - - - - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'Qair',id%Qair)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Qair in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Get Qair units:
    status = NF90_GET_ATT(ncid_met,id%Qair,'units',metunits%Qair)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Qair units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Qair(1:1)=='%'.OR.metunits%Qair(1:1)=='-') THEN
       ! Change from relative humidity to specific humidity:
       convert%Qair = -999.0
       WRITE(logn,*) 'Humidity will be converted from relative to specific'
    ELSE IF(metunits%Qair(1:3)=='g/g'.OR.metunits%Qair(1:5)=='kg/kg' &
         .OR.metunits%Qair(1:3)=='G/G'.OR.metunits%Qair(1:5)=='KG/KG') THEN
          ! Units are correct
       convert%Qair=1.0
    ELSE
       WRITE(*,*) metunits%Qair
       CALL abort('Unknown units for Qair'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Rainf (essential):- - - - - - - - - - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'Rainf',id%Rainf)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Rainf in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    ! Get Rainf units:
    status = NF90_GET_ATT(ncid_met,id%Rainf,'units',metunits%Rainf)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Rainf units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Rainf(1:8)=='kg/m^2/s'.OR.metunits%Rainf(1:10)== &
         'kgm^-2s^-1'.OR.metunits%Rainf(1:4)=='mm/s'.OR. &
         metunits%Rainf(1:6)=='mms^-1'.OR. &
         metunits%Rainf(1:7)=='kg/m^2s') THEN
       ! Change from mm/s to mm/time step:
       convert%Rainf = dels
    ELSE IF(metunits%Rainf(1:4)=='mm/h'.OR.metunits%Rainf(1:6)== &
         'mmh^-1') THEN
       ! Change from mm/h to mm/time step:
       convert%Rainf = dels/3600.0
    ELSE
       WRITE(*,*) metunits%Rainf
       CALL abort('Unknown units for Rainf'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Multiply acceptable Rainf ranges by time step size:
    ranges%Rainf = ranges%Rainf*dels ! range therefore depends on dels
    ! Look for Wind (essential):- - - - - - - - - - - - - - - - - - -
    status = NF90_INQ_VARID(ncid_met,'Wind',id%Wind)
    IF(status /= NF90_NOERR) THEN
       ! Look for vector wind:
       status = NF90_INQ_VARID(ncid_met,'Wind_N',id%Wind)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding Wind in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       status = NF90_INQ_VARID(ncid_met,'Wind_E',id%Wind_E)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding Wind_E in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       exists%Wind = .FALSE. ! Use vector wind when reading met
    ELSE
       exists%Wind = .TRUE. ! 'Wind' variable exists
    END IF
    ! Get Wind units:
    status = NF90_GET_ATT(ncid_met,id%Wind,'units',metunits%Wind)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding Wind units in met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Wind(1:3)/='m/s'.AND.metunits%Wind(1:2)/='ms') THEN
       WRITE(*,*) metunits%Wind
       CALL abort('Unknown units for Wind'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Now "optional" variables:
    ! Look for LWdown (can be synthesised):- - - - - - - - - - - - - - -
    status = NF90_INQ_VARID(ncid_met,'LWdown',id%LWdown)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       exists%LWdown = .TRUE. ! LWdown is present in met file
       ! Get LWdown units and check okay:
       status = NF90_GET_ATT(ncid_met,id%LWdown,'units',metunits%LWdown)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding LWdown units in met data file ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       IF(metunits%LWdown(1:4)/='W/m2'.AND.metunits%LWdown(1:5) &
            /='W/m^2'.AND.metunits%LWdown(1:5)/='Wm^-2' &
            .AND.metunits%LWdown(1:4)/='Wm-2') THEN
          WRITE(*,*) metunits%LWdown
          CALL abort('Unknown units for LWdown'// &
               ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE
       exists%LWdown = .FALSE. ! LWdown is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,*) 'LWdown not present in met file; ', &
            'values will be synthesised based on air temperature.'
    END IF
    ! Look for PSurf (can be synthesised):- - - - - - - - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'PSurf',id%PSurf)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       exists%PSurf = .TRUE. ! PSurf is present in met file
       ! Get PSurf units and check:
       status = NF90_GET_ATT(ncid_met,id%PSurf,'units',metunits%PSurf)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding PSurf units in met data file ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       IF(metunits%PSurf(1:2)=='Pa'.OR.metunits%PSurf(1:2)=='pa'.OR. &
            metunits%PSurf(1:2)=='PA' ) THEN
          ! Change from pa to mbar (cable uses mbar):
          convert%PSurf = 0.01
          WRITE(logn,*) 'Pressure will be converted from Pa to mb'
       ELSE IF(metunits%PSurf(1:2)=='KP'.OR.metunits%PSurf(1:2)=='kP' &
            .OR.metunits%PSurf(1:2)=='Kp'.OR.metunits%PSurf(1:2)=='kp') THEN
          ! convert from kPa to mb
          convert%PSurf = 10.0
          WRITE(logn,*) 'Pressure will be converted from kPa to mb'
       ELSE IF(metunits%PSurf(1:2)=='MB'.OR.metunits%PSurf(1:2)=='mB' &
            .OR.metunits%PSurf(1:2)=='Mb'.OR.metunits%PSurf(1:2)=='mb') THEN
          ! Units are correct
          convert%PSurf = 1.0
       ELSE
          WRITE(*,*) metunits%PSurf
          CALL abort('Unknown units for PSurf'// &
               ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! If PSurf not present
       exists%PSurf = .FALSE. ! PSurf is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Look for "elevation" variable to assume static pressure:
       status = NF90_INQ_VARID(ncid_met,'Elevation',id%Elev)
       IF(status == NF90_NOERR) THEN ! elevation present
          ! Get elevation units:
          status = NF90_GET_ATT(ncid_met,id%Elev,'units',metunits%Elev)
          CALL nc_abort &
               ('Error finding elevation units in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
          ! Allocate space for elevation variable:
          ALLOCATE(elevation(xdimsize,ydimsize))
          ! Get site elevation:
          status=NF90_GET_VAR(ncid_met,id%Elev,elevation)
          CALL nc_abort &
               ('Error reading elevation variable in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       ELSE
          ! If both PSurf and elevation aren't present, abort:
          CALL nc_abort &
               ('Error finding PSurf or Elevation in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       END IF
       ! Note static pressure based on elevation in log file:
       WRITE(logn,*) 'PSurf not present in met file; values will be ', &
            'synthesised based on elevation and temperature.'
    END IF
    ! Look for CO2air (can be assumed to be static):- - - - - - - - - - -
    status = NF90_INQ_VARID(ncid_met,'CO2air',id%CO2air)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       exists%CO2air = .TRUE. ! CO2air is present in met file
       ! Get CO2air units:
       status = NF90_GET_ATT(ncid_met,id%CO2air,'units',metunits%CO2air)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding CO2air units in met data file ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       IF(metunits%CO2air(1:3)/='ppm') THEN
          WRITE(*,*) metunits%CO2air
          CALL abort('Unknown units for CO2air'// &
               ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! CO2 not present
       exists%CO2air = .FALSE. ! CO2air is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,'(A32,A24,I4,A5)') ' CO2air not present in met file; ', &
            'values will be fixed at ',INT(fixedCO2),' ppmv'
    END IF
    ! Look for Snowf (could be part of Rainf variable):- - - - - - - - - - 
    status = NF90_INQ_VARID(ncid_met,'Snowf',id%Snowf)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       exists%Snowf = .TRUE. ! Snowf is present in met file
       ! Get Snowf units:
       status = NF90_GET_ATT(ncid_met,id%Snowf,'units',metunits%Snowf)
       IF(status /= NF90_NOERR) CALL nc_abort &
            ('Error finding Snowf units in met data file ' &
            //TRIM(filename_met)//' (SUBROUTINE open_met_file)')
       ! Make sure Snowf units are the same as Rainf units:
       IF(metunits%Rainf/=metunits%Snowf) CALL abort &
            ('Please ensure Rainf and Snowf units are the same'// &
            ' in '//TRIM(filename_met)//' (SUBROUTINE open_met_data)')
    ELSE
       exists%Snowf = .FALSE. ! Snowf is not present in met file
       !  all_met=.FALSE. not required; Snowf assumed to be in Rainf
       ! Note this in log file:
       WRITE(logn,*) 'Snowf not present in met file; ', &
            'Assumed to be contained in Rainf variable'
    END IF
    ! Look for LAI - - - - - - - - - - - - - - - - - - - - - - - - -
    status = NF90_INQ_VARID(ncid_met,'LAI',id%LAI)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       exists%LAI = .TRUE. ! LAI is present in met file
       ! LAI will be read in which ever land grid is used
    ELSE
       exists%LAI = .FALSE. ! LAI is not present in met file
       ! Report to log file
       WRITE(logn,*) 'LAI not present in met file; ', &
            'Will use MODIS coarse grid monthly LAI'
    END IF
    ! Look for veg type:
    status = NF90_INQ_VARID(ncid_met,'iveg',id%iveg)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       ! Allocate space for user-defined veg type variable:
       ALLOCATE(vegtype(mp))
       ! Get veg type from met file:
       IF(gridType=='mask') THEN
          DO i = 1, mp
             status= NF90_GET_VAR(ncid_met,id%iveg,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading iveg in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             vegtype(i)=data2(1,1)
          END DO
       ELSE IF(gridType=='land') THEN
          ! Collect data from land only grid in netcdf file:
          status= NF90_GET_VAR(ncid_met,id%iveg,vegtype)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading iveg in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')

       END IF
    END IF
    ! Look for soil type:
    status = NF90_INQ_VARID(ncid_met,'isoil',id%isoil)
    IF(status == NF90_NOERR) THEN ! If inquiry is okay
       ! Allocate space for user-defined soil type variable:
       ALLOCATE(soiltype(mp))
       ! Get soil type from met file:
       IF(gridType=='mask') THEN
          DO i = 1, mp
             status= NF90_GET_VAR(ncid_met,id%isoil,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading isoil in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             soiltype(i)=data2(1,1)
          END DO
       ELSE IF(gridType=='land') THEN
          ! Collect data from land only grid in netcdf file:
          status= NF90_GET_VAR(ncid_met,id%isoil,soiltype)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading isoil in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')

       END IF
    END IF

    ! Report finding met variables to log file:
    IF(all_met) THEN
       WRITE(logn,*) 'Found all met variables in met file.'
    ELSE
       WRITE(logn,*) 'Found all ESSENTIAL met variables in met file,', &
            ' some synthesised (as above).'
    END IF
    !!=================^^ End met variables search^^=======================
  END SUBROUTINE open_met_file
  !========================================================================
  SUBROUTINE get_met_data(ktau,filename_met,met,soil,rad,veg,kend,dels) 
    ! Fetches meteorological forcing data for a single time step,
    ! including LAI if it exists.
    INTEGER(i_d), INTENT(IN) :: ktau ! time step number
    CHARACTER(LEN=*),INTENT(IN) :: filename_met 
    TYPE(met_type),INTENT(OUT):: met ! meteorological data
    TYPE (soil_parameter_type),INTENT(IN)	:: soil	
    TYPE (radiation_type),INTENT(IN)	:: rad
    TYPE(veg_parameter_type),INTENT(INOUT) :: veg ! LAI retrieved from file
    INTEGER(i_d), INTENT(IN) :: kend ! total number of timesteps in run
    REAL(r_1),INTENT(IN) :: dels ! time step size
    REAL(r_1),DIMENSION(1,1,1) :: data3 ! temp variable for netcdf reading
    REAL(r_1),DIMENSION(1,1,1,1) :: data4 !  " " "
    REAL(r_1),DIMENSION(1,1)    :: data2 ! " " 
    INTEGER(i_d) :: i ! do loop counter

    DO i=1,mp ! over all land points/grid cells
       ! First set timing variables:
       IF(ktau==1) THEN ! initialise...
          SELECT CASE(time_coord)
          CASE('LOC')! i.e. use local time by default
             ! hour-of-day = starting hod 
             met%hod(i) = shod 
             met%doy(i) = sdoy
             met%moy(i) = smoy
             met%year(i) = syear
          CASE('GMT')! use GMT
             ! hour-of-day = starting hod + offset from GMT time:
             met%hod(i) = shod + (longitude(i)/180.0)*12.0
             met%doy(i) = sdoy
             met%moy(i) = smoy
             met%year(i) = syear
          CASE DEFAULT
             CALL abort('Unknown time coordinate! ' &
                  //' (SUBROUTINE get_met_data)')
          END SELECT
       ELSE
          ! increment hour-of-day by time step size:
          met%hod(i) = met%hod(i) + dels/3600.0
       END IF
       ! 
       IF(met%hod(i)<0.0) THEN ! may be -ve since longitude [-180,180]
          ! Reduce day-of-year by one and ammend hour-of-day:
          met%doy(i) = met%doy(i) - 1
          met%hod(i) = met%hod(i) + 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(i)))
             CASE(0) ! ie Dec previous year
                met%moy(i) = 12
                met%year(i) = met%year(i) - 1
                met%doy(i) = 365 ! prev year not a leap year as this is
             CASE(31) ! Jan
                met%moy(i) = 1
             CASE(60) ! Feb
                met%moy(i) = 2
             CASE(91) ! Mar
                met%moy(i) = 3
             CASE(121)
                met%moy(i) = 4
             CASE(152)
                met%moy(i) = 5
             CASE(182)
                met%moy(i) = 6
             CASE(213)
                met%moy(i) = 7
             CASE(244)
                met%moy(i) = 8
             CASE(274)
                met%moy(i) = 9
             CASE(305)
                met%moy(i) = 10
             CASE(335)
                met%moy(i) = 11
             END SELECT
          ELSE ! not a leap year or not using leap year timing
             SELECT CASE(INT(met%doy(i)))
             CASE(0) ! ie Dec previous year
                met%moy(i) = 12
                met%year(i) = met%year(i) - 1
                ! If previous year is a leap year
                IF((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
                     (MOD(syear,4)==0.AND.MOD(syear,400)==0)) THEN
                   met%doy(i) = 366
                ELSE
                   met%doy(i) = 365
                END IF
             CASE(31) ! Jan
                met%moy(i) = 1 
             CASE(59) ! Feb
                met%moy(i) = 2
             CASE(90)
                met%moy(i) = 3
             CASE(120)
                met%moy(i) = 4
             CASE(151)
                met%moy(i) = 5
             CASE(181)
                met%moy(i) = 6
             CASE(212)
                met%moy(i) = 7
             CASE(243)
                met%moy(i) = 8
             CASE(273)
                met%moy(i) = 9
             CASE(304)
                met%moy(i) = 10
             CASE(334) 
                met%moy(i) = 11
             CASE(36) ! end of year; increment
                syear = syear + 1 
             END SELECT
          END IF ! if leap year or not
       ELSE IF(met%hod(i)>=24.0) THEN ! increment or GMT adj has shifted day
          ! Adjust day-of-year and hour-of-day:
          met%doy(i) = met%doy(i) + 1
          met%hod(i) = met%hod(i) - 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(i)))
             CASE(32) ! Feb
                met%moy(i) = 2
             CASE(61) ! Mar
                met%moy(i) = 3
             CASE(92)
                met%moy(i) = 4
             CASE(122)
                met%moy(i) = 5
             CASE(153)
                met%moy(i) = 6
             CASE(183)
                met%moy(i) = 7
             CASE(214)
                met%moy(i) = 8
             CASE(245)
                met%moy(i) = 9
             CASE(275)
                met%moy(i) = 10
             CASE(306)
                met%moy(i) = 11
             CASE(336) 
                met%moy(i) = 12
             CASE(367)! end of year; increment
                met%year(i) = met%year(i) + 1 
                met%moy(i) = 1
                met%doy(i) = 1
             END SELECT
          ! ELSE IF not leap year and Dec 31st, increment year
          ELSE 
             SELECT CASE(INT(met%doy(i)))
             CASE(32) ! Feb
                met%moy(i) = 2
             CASE(60) ! Mar
                met%moy(i) = 3
             CASE(91)
                met%moy(i) = 4
             CASE(121)
                met%moy(i) = 5
             CASE(152)
                met%moy(i) = 6
             CASE(182)
                met%moy(i) = 7
             CASE(213)
                met%moy(i) = 8
             CASE(244)
                met%moy(i) = 9
             CASE(274)
                met%moy(i) = 10
             CASE(305)
                met%moy(i) = 11
             CASE(335) 
                met%moy(i) = 12
             CASE(366)! end of year; increment
                met%year(i) = met%year(i) + 1 
                met%moy(i) = 1
                met%doy(i) = 1
             END SELECT
          END IF ! if leap year or not
       END IF ! if increment has pushed hod to a different day
          
       IF(gridType=='mask') THEN
          ! Get SWdown data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%SWdown,data3, &
               start=(/land_x(i),land_y(i),ktau/),count=(/1,1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading SWdown in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable (no units change required):
          met%fsd(i) = data3(1,1,1)
          ! Get Tair data:- - - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Tair,data4, &
               start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Tair in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable with units change:
          met%tk(i) = data4(1,1,1,1) + convert%Tair
          met%tc(i) =  met%tk(i) - tfrz
          ! Get PSurf data:- - - - - - - - - - - - - - - - - - - - -
          IF(exists%PSurf) THEN ! IF PSurf is in met file:
             status= NF90_GET_VAR(ncid_met,id%PSurf,data4, &
                  start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading PSurf in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%pmb(i) = data4(1,1,1,1) * convert%PSurf
          ELSE ! PSurf must be fixed as a function of site elevation and T:
             met%pmb(i)=101.325*(met%tk(i)/(met%tk(i) + 0.0065* &
                  elevation(land_x(i),land_y(i))))**(9.80665/287.04/0.0065)
          END IF
          ! Get Qair data:- - - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Qair,data4, &
               start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Qair in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          IF(convert%Qair==-999.0) THEN
             CALL rh_sh(data4(1,1,1,1), met%tk(i), &
                  met%pmb(i),met%qv(i))
          ELSE
             met%qv(i) = data4(1,1,1,1)
          END IF
          ! Get Wind data:- - - - - - - - - - - - - - - - - - - - - -
          IF(exists%Wind) THEN ! Scalar Wind
             status= NF90_GET_VAR(ncid_met,id%Wind,data4, &
                  start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             ! Assign value to met data variable (no units change required):
             met%ua(i) = data4(1,1,1,1)
          ELSE ! Vector wind
             ! Get Wind_N:
             status= NF90_GET_VAR(ncid_met,id%Wind,data4, &
                  start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind_N in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%ua(i) = data4(1,1,1,1) ! only part of the wind variable
             status= NF90_GET_VAR(ncid_met,id%Wind_E,data4, &
                  start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind_E in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             ! Write final scalar Wind value:
             met%ua(i) = SQRT(met%ua(i)**2 + data4(1,1,1,1)**2)
          END IF
          ! Get Rainf and Snowf data:- - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Rainf,data3, &
               start=(/land_x(i),land_y(i),ktau/),count=(/1,1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Rainf in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          met%precip(i) = data3(1,1,1) ! store Rainf value
          IF(exists%Snowf) THEN
             status= NF90_GET_VAR(ncid_met,id%Snowf,data3, &
                  start=(/land_x(i),land_y(i),ktau/),count=(/1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Snowf in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             data3(1,1,1) = data3(1,1,1) + met%precip(i)
          END IF
          ! Convert units:
          met%precip(i) = data3(1,1,1) * convert%Rainf
          ! Get LWdown data:- - - - - - - - - - - - - - - - - - - - - - - - 
          IF(exists%LWdown) THEN ! If LWdown exists in met file
             status= NF90_GET_VAR(ncid_met,id%LWdown,data3, &
                  start=(/land_x(i),land_y(i),ktau/),count=(/1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading LWdown in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%fld(i)=data3(1,1,1)
          ELSE ! Synthesise LWdown based on temperature
             ! Use Swinbank formula:
             met%fld(i)=0.0000094*0.0000000567*(met%tk(i)**6.0)
          END IF
          ! Get CO2air data:- - - - - - - - - - - - - - - - - - - - - - - -
          IF(exists%CO2air) THEN ! If CO2air exists in met file
             status= NF90_GET_VAR(ncid_met,id%CO2air,data4, &
                  start=(/land_x(i),land_y(i),1,ktau/),count=(/1,1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading CO2air in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%ca(i) = data4(1,1,1,1)/1000000.0
          ELSE 
             ! Fix CO2 air concentration:
             met%ca(i) = fixedCO2 /1000000.0
          END IF
          ! Get LAI if it's present: - - - - - - - - - - - - - - - - -
          IF(exists%LAI) THEN ! If LAI exists in met file
             status= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                  start=(/land_x(i),land_y(i),ktau/),count=(/1,1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading LAI in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             veg%vlai(i) = data3(1,1,1)
          ELSE 
             ! If not in met file, load default LAI value:
             CALL get_default_lai(ktau,met%doy,veg,kend)
          END IF
          
          ! Set solid precip based on temp
          met%precips = 0.0
          IF( met%tc(i) <= 0.0 ) met%precips(i) = met%precip(i)

       ELSE IF(gridType=='land') THEN
          ! Collect data from alnd only grid in netcdf file:
          ! Get SWdown data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%SWdown,data2, &
               start=(/i,ktau/),count=(/1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading SWdown in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable (no units change required):
          met%fsd(i) = data2(1,1)
          ! Get Tair data:- - - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Tair,data2, &
               start=(/i,ktau/),count=(/1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Tair in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable with units change:
          met%tk(i) = data2(1,1) + convert%Tair
          met%tc(i) =  met%tk(i) - tfrz
          ! Get PSurf data:- - - - - - - - - - - - - - - - - - - - -
          IF(exists%PSurf) THEN ! IF PSurf is in met file:
             status= NF90_GET_VAR(ncid_met,id%PSurf,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading PSurf in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%pmb(i) = data2(1,1) * convert%PSurf
          ELSE ! PSurf must be fixed as a function of site elevation and T:
             met%pmb(i)=101.325*(met%tk(i)/(met%tk(i) + 0.0065* &
                  elevation(land_x(i),land_y(i))))**(9.80665/287.04/0.0065)
          END IF
          ! Get Qair data:- - - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Qair,data2, &
               start=(/i,ktau/),count=(/1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Qair in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          IF(convert%Qair==-999.0) THEN
             CALL rh_sh(data2(1,1), met%tk(i), &
                  met%pmb(i),met%qv(i))
          ELSE
             met%qv(i) = data2(1,1)
          END IF
          ! Get Wind data:- - - - - - - - - - - - - - - - - - - - - -
          IF(exists%Wind) THEN ! Scalar Wind
             status= NF90_GET_VAR(ncid_met,id%Wind,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             ! Assign value to met data variable (no units change required):
             met%ua(i) = data2(1,1)
          ELSE ! Vector wind
             ! Get Wind_N:
             status= NF90_GET_VAR(ncid_met,id%Wind,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind_N in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%ua(i) = data2(1,1) ! only part of the wind variable
             status= NF90_GET_VAR(ncid_met,id%Wind_E,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Wind_E in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             ! Write final scalar Wind value:
             met%ua(i) = SQRT(met%ua(i)**2 + data2(1,1)**2)
          END IF
          ! Get Rainf and Snowf data:- - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,id%Rainf,data2, &
               start=(/i,ktau/),count=(/1,1/))
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading Rainf in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
          met%precip(i) = data2(1,1) ! store Rainf value
          IF(exists%Snowf) THEN
             status= NF90_GET_VAR(ncid_met,id%Snowf,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading Snowf in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             data2(1,1) = data2(1,1) + met%precip(i)
          END IF
          ! Convert units:
          met%precip(i) = data2(1,1) * convert%Rainf
          ! Get LWdown data:- - - - - - - - - - - - - - - - - - - - - - - - 
          IF(exists%LWdown) THEN ! If LWdown exists in met file
             status= NF90_GET_VAR(ncid_met,id%LWdown,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading LWdown in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%fld(i)=data2(1,1)
          ELSE ! Synthesise LWdown based on temperature
             ! Use Swinbank formula:
             met%fld(i)=0.0000094*0.0000000567*(met%tk(i)**6.0)
          END IF
          ! Get CO2air data:- - - - - - - - - - - - - - - - - - - - - - - -
          IF(exists%CO2air) THEN ! If CO2air exists in met file
             status= NF90_GET_VAR(ncid_met,id%CO2air,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading CO2air in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             met%ca(i) = data2(1,1)/1000000.0
          ELSE 
             ! Fix CO2 air concentration:
             met%ca(i) = fixedCO2 /1000000.0
          END IF
          ! Get LAI data:- - - - - - - - - - - - - - - - - - - - - - - -
          IF(exists%LAI) THEN ! If LAI exists in met file
             status= NF90_GET_VAR(ncid_met,id%LAI,data2, &
                  start=(/i,ktau/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading LAI in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_met_data)')
             veg%vlai(i) = data2(1,1)
          ELSE 
             ! If not in met file, load default LAI value:
             CALL get_default_lai(ktau,met%doy,veg,kend)
          END IF

          ! Set solid precip based on temp
          met%precips = 0.0
          IF( met%tc(i) <= 0.0 ) met%precips(i) = met%precip(i)

       ELSE
          CALL abort('Unrecognised grid type')
       END IF ! grid type
    END DO
    ! Set cosine of zenith angle (provided by GCM when online):
    met%coszen = sinbet(met%doy, rad%latitude, met%hod)
    ! initialise within canopy air temp
    met%tvair = met%tk 
    met%tvrad = met%tk 
    IF(check%ranges) THEN
       ! Check ranges are okay:
       IF(ANY(met%fsd<ranges%SWdown(1)).OR.ANY(met%fsd>ranges%SWdown(2))) &
            CALL abort('SWdown out of specified ranges!')
       IF(ANY(met%fld<ranges%LWdown(1)).OR.ANY(met%fld>ranges%LWdown(2))) &
            CALL abort('LWdown out of specified ranges!')
       IF(ANY(met%qv<ranges%Qair(1)).OR.ANY(met%qv>ranges%Qair(2))) &
            CALL abort('Qair out of specified ranges!')
       IF(ANY(met%precip<ranges%Rainf(1)).OR.ANY(met%precip>ranges%Rainf(2))) &
            CALL abort('Rainf out of specified ranges!')
       IF(ANY(met%ua<ranges%Wind(1)).OR.ANY(met%ua>ranges%Wind(2))) &
            CALL abort('Wind out of specified ranges!')
       IF(ANY(met%tk<ranges%Tair(1)).OR.ANY(met%tk>ranges%Tair(2))) &
            CALL abort('Tair out of specified ranges!')
       IF(ANY(met%pmb<ranges%PSurf(1)).OR.ANY(met%pmb>ranges%PSurf(2))) &
            CALL abort('PSurf out of specified ranges!')
    END IF

  END SUBROUTINE get_met_data
!============================================================================
  SUBROUTINE close_met_file(filename_met)
    CHARACTER(LEN=*), INTENT(IN) :: filename_met
    status=NF90_CLOSE(ncid_met)
    IF(status /= NF90_NOERR) CALL nc_abort ('Error closing met data file ' &
         //TRIM(filename_met)//' (SUBROUTINE close_met_file)')
     ! Clear lat_all and lon_all variables
    DEALLOCATE(lat_all,lon_all)
  END SUBROUTINE close_met_file
  !=========================================================================
  SUBROUTINE get_restart_data(filename_restart_in,filename_met,logn, &
       ssoil,canopy,rough,bgc,bal,veg,soil,rad,sum_flux)
    ! Reads restart file, if available, and checks its compatibility.
    ! Initialisations and parameter values will be loaded.
    CHARACTER(LEN=*), INTENT(IN) :: filename_restart_in
    CHARACTER(LEN=*), INTENT(IN) :: filename_met
    INTEGER, INTENT(IN) :: logn ! log file number
    TYPE (soil_snow_type),INTENT(INOUT)      :: ssoil  ! soil and snow variables
    TYPE (bgc_pool_type),INTENT(INOUT)       :: bgc    ! carbon pool variables
    TYPE (canopy_type),INTENT(INOUT)         :: canopy ! vegetation variables
    TYPE (roughness_type),INTENT(INOUT)      :: rough  ! roughness varibles
    TYPE (balances_type),INTENT(INOUT) :: bal  ! energy and water balance variables
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(OUT) :: soil ! soil parameters
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (parID_type) :: ipid ! input parameter IDs in netcdf restart file
    REAL(r_1), POINTER,DIMENSION(:) :: lat_restart, lon_restart
    INTEGER :: mp_restart ! number of land points in restart file
    INTEGER :: mpID ! dimension IDs
    INTEGER :: latID,lonID ! time,lat,lon variable ID
    INTEGER :: tggID,wbID,wbiceID,tssID,ssdnnID,ssdnID,osnowdID, &
         smassID,sdepthID,snageID,snowdID,rtsoilID,isflagID,canstoID, &
         albsoilsnID,gammzzID,tggsnID,sghfluxID,ghfluxID,runoffID, &
         rnof1ID,rnof2ID,gaID,dgdtgID,fevID,fesID,fhsID,wbtot0ID, &
         osnowd0ID,cplantID,csoilID
    
    ! Write to screen the restart file is found:
    WRITE(*,*) 'Reading restart data from: ' ,TRIM(filename_restart_in)
    ! Check number of gridpoints in restart file is correct:
    status = NF90_INQ_DIMID(ncid_rin,'mp',mpID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding mp dimension in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status = NF90_INQUIRE_DIMENSION(ncid_rin,mpID,len=mp_restart)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding number of land points in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    IF(mp_restart /= mp) CALL abort('Number of land points in '// &
         'restart file '//TRIM(filename_restart_in)// &
         ' differs from number in met file '//TRIM(filename_met))
    ! check that lat/lon b/w run and restart are compatible:
    ALLOCATE(lat_restart(mp),lon_restart(mp))
    status = NF90_INQ_VARID(ncid_rin,'latitude',latID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding latitude in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status = NF90_INQ_VARID(ncid_rin,'longitude',lonID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding longitude in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,latID,lat_restart)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading latitude in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! Set rad%latitude parameter
    rad%latitude=lat_restart
    status=NF90_GET_VAR(ncid_rin,lonID,lon_restart)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading longitude in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    IF(ANY(lat_restart/=latitude)) CALL abort('Latitude of land points in '// &
         'restart file '//TRIM(filename_restart_in)// &
         ' differs from met file '//TRIM(filename_met))
    IF(ANY(lon_restart/=longitude)) CALL abort('Longitude of land points in '// &
         'restart file '//TRIM(filename_restart_in)// &
         ' differs from met file '//TRIM(filename_met))
    DEALLOCATE(lat_restart,lon_restart)
    ! Get variable initialisations =============================
    ! tgg:
    status = NF90_INQ_VARID(ncid_rin,'tgg',tggID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding tgg in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,tggID,ssoil%tgg)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading tgg in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! wb:
    status = NF90_INQ_VARID(ncid_rin,'wb',wbID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding wb in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,wbID,ssoil%wb)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading wb in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! wbice:
    status = NF90_INQ_VARID(ncid_rin,'wbice',wbiceID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding wbice in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,wbiceID,ssoil%wbice)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading wbice in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! tss:
    status = NF90_INQ_VARID(ncid_rin,'tss',tssID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding tss in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,tssID,ssoil%tss)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading tss in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ssdnn:
    status = NF90_INQ_VARID(ncid_rin,'ssdnn',ssdnnID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ssdnn in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ssdnnID,ssoil%ssdnn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ssdnn in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ssdn:
    status = NF90_INQ_VARID(ncid_rin,'ssdn',ssdnID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ssdn in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ssdnID,ssoil%ssdn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ssdn in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! osnowd:
    status = NF90_INQ_VARID(ncid_rin,'osnowd',osnowdID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding osnowd in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,osnowdID,ssoil%osnowd)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading osnowd in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! smass:
    status = NF90_INQ_VARID(ncid_rin,'smass',smassID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding smass in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,smassID,ssoil%smass)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading smass in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! sdepth:
    status = NF90_INQ_VARID(ncid_rin,'sdepth',sdepthID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding sdepth in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,sdepthID,ssoil%sdepth)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading sdepth in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! snage:
    status = NF90_INQ_VARID(ncid_rin,'snage',snageID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding snage in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,snageID,ssoil%snage)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading snage in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! snowd:
    status = NF90_INQ_VARID(ncid_rin,'snowd',snowdID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding snowd in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,snowdID,ssoil%snowd)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading snowd in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rtsoil:
    status = NF90_INQ_VARID(ncid_rin,'rtsoil',rtsoilID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rtsoil in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,rtsoilID,ssoil%rtsoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rtsoil in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! isflag:
    status = NF90_INQ_VARID(ncid_rin,'isflag',isflagID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding isflag in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,isflagID,ssoil%isflag)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading isflag in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! cansto:
    status = NF90_INQ_VARID(ncid_rin,'cansto',canstoID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding cansto in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,canstoID,canopy%cansto)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading cansto in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! albsoilsn:
    status = NF90_INQ_VARID(ncid_rin,'albsoilsn',albsoilsnID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding albsoilsn in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,albsoilsnID,ssoil%albsoilsn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading albsoilsn in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! gammzz:
    status = NF90_INQ_VARID(ncid_rin,'gammzz',gammzzID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding gammzz in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,gammzzID,ssoil%gammzz)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading gammzz in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! tggsn:
    status = NF90_INQ_VARID(ncid_rin,'tggsn',tggsnID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding tggsn in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,tggsnID,ssoil%tggsn)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading tggsn in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! sghflux:
    status = NF90_INQ_VARID(ncid_rin,'sghflux',sghfluxID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding sghflux in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,sghfluxID,canopy%sghflux)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading sghflux in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ghflux:
    status = NF90_INQ_VARID(ncid_rin,'ghflux',ghfluxID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ghflux in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ghfluxID,canopy%ghflux)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ghflux in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! runoff:
    status = NF90_INQ_VARID(ncid_rin,'runoff',runoffID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding runoff in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,runoffID,ssoil%runoff)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading runoff in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rnof1:
    status = NF90_INQ_VARID(ncid_rin,'rnof1',rnof1ID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rnof1 in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,rnof1ID,ssoil%rnof1)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rnof1 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rnof2:
    status = NF90_INQ_VARID(ncid_rin,'rnof2',rnof2ID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rnof2 in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,rnof2ID,ssoil%rnof2)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rnof2 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ga:
    status = NF90_INQ_VARID(ncid_rin,'ga',gaID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ga in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,gaID,canopy%ga)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ga in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! dgdtg:
    status = NF90_INQ_VARID(ncid_rin,'dgdtg',dgdtgID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding dgdtg in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,dgdtgID,canopy%dgdtg)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading dgdtg in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! fev:
    status = NF90_INQ_VARID(ncid_rin,'fev',fevID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding fev in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,fevID,canopy%fev)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading fev in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! fes:
    status = NF90_INQ_VARID(ncid_rin,'fes',fesID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding fes in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,fesID,canopy%fes)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading fes in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! fhs:
    status = NF90_INQ_VARID(ncid_rin,'fhs',fhsID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding fhs in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,fhsID,canopy%fhs)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading fhs in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! cplant:
    status = NF90_INQ_VARID(ncid_rin,'cplant',cplantID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding cplant in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,cplantID,bgc%cplant)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading cplant in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! csoil:
    status = NF90_INQ_VARID(ncid_rin,'csoil',csoilID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding csoil in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,csoilID,bgc%csoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading csoil in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! wbtot0:
    status = NF90_INQ_VARID(ncid_rin,'wbtot0',wbtot0ID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding wbtot0 in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,wbtot0ID,bal%wbtot0)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading wbtot0 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! osnowd0:
    status = NF90_INQ_VARID(ncid_rin,'osnowd0',osnowd0ID)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding osnowd0 in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,osnowd0ID,bal%osnowd0)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading osnowd0 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! Get model parameters =============================================
    ! rad%latitude set above in lat/lon checking section     
    ! iveg:
    status = NF90_INQ_VARID(ncid_rin,'iveg',ipid%iveg)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding iveg parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%iveg,veg%iveg)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading iveg in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! isoil:
    status = NF90_INQ_VARID(ncid_rin,'isoil',ipid%isoil)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding isoil parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%isoil,soil%isoilm)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading isoil in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! clay:
    status = NF90_INQ_VARID(ncid_rin,'clay',ipid%clay)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding clay parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%clay,soil%clay)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading clay in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! sand:
    status = NF90_INQ_VARID(ncid_rin,'sand',ipid%sand)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding sand parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%sand,soil%sand)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading sand in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! silt:
    status = NF90_INQ_VARID(ncid_rin,'silt',ipid%silt)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding silt parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%silt,soil%silt)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading silt in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ssat:
    status = NF90_INQ_VARID(ncid_rin,'ssat',ipid%ssat)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ssat parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%ssat,soil%ssat)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ssat in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! sfc:
    status = NF90_INQ_VARID(ncid_rin,'sfc',ipid%sfc)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding sfc parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%sfc,soil%sfc)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading sfc in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! swilt:
    status = NF90_INQ_VARID(ncid_rin,'swilt',ipid%swilt)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding swilt parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%swilt,soil%swilt)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading swilt in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! zse:
    status = NF90_INQ_VARID(ncid_rin,'zse',ipid%zse)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding zse parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%zse,soil%zse)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading zse in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! froot:
    status = NF90_INQ_VARID(ncid_rin,'froot',ipid%froot)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding froot parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%froot,veg%froot)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading froot in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! bch:
    status = NF90_INQ_VARID(ncid_rin,'bch',ipid%bch)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding bch parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%bch,soil%bch)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading bch in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! hyds:
    status = NF90_INQ_VARID(ncid_rin,'hyds',ipid%hyds)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding hyds parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%hyds,soil%hyds)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading hyds in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! sucs:
    status = NF90_INQ_VARID(ncid_rin,'sucs',ipid%sucs)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding sucs parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%sucs,soil%sucs)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading sucs in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! css:
    status = NF90_INQ_VARID(ncid_rin,'css',ipid%css)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding css parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%css,soil%css)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading css in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rhosoil:
    status = NF90_INQ_VARID(ncid_rin,'rhosoil',ipid%rhosoil)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rhosoil parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%rhosoil,soil%rhosoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rhosoil in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rs20:
    status = NF90_INQ_VARID(ncid_rin,'rs20',ipid%rs20)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rs20 parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%rs20,soil%rs20)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rs20 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! albsoil:
    status = NF90_INQ_VARID(ncid_rin,'albsoil',ipid%albsoil)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding albsoil parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%albsoil,soil%albsoil)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading albsoil in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! hc:
    status = NF90_INQ_VARID(ncid_rin,'hc',ipid%hc)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding hc parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%hc,veg%hc)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading hc in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! canst1:
    status = NF90_INQ_VARID(ncid_rin,'canst1',ipid%canst1)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding canst1 parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%canst1,veg%canst1)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading canst1 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! dleaf:
    status = NF90_INQ_VARID(ncid_rin,'dleaf',ipid%dleaf)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding dleaf parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%dleaf,veg%dleaf)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading dleaf in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! frac4:
    status = NF90_INQ_VARID(ncid_rin,'frac4',ipid%frac4)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding frac4 parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%frac4,veg%frac4)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading frac4 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ejmax:
    status = NF90_INQ_VARID(ncid_rin,'ejmax',ipid%ejmax)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ejmax parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%ejmax,veg%ejmax)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ejmax in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! vcmax:
    status = NF90_INQ_VARID(ncid_rin,'vcmax',ipid%vcmax)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding vcmax parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%vcmax,veg%vcmax)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading vcmax in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rp20:
    status = NF90_INQ_VARID(ncid_rin,'rp20',ipid%rp20)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rp20 parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%rp20,veg%rp20)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rp20 in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! rpcoef:
    status = NF90_INQ_VARID(ncid_rin,'rpcoef',ipid%rpcoef)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding rpcoef parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%rpcoef,veg%rpcoef)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading rpcoef in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! shelrb:
    status = NF90_INQ_VARID(ncid_rin,'shelrb',ipid%shelrb)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding shelrb parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%shelrb,veg%shelrb)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading shelrb in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! xfang:
    status = NF90_INQ_VARID(ncid_rin,'xfang',ipid%xfang)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding xfang parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%xfang,veg%xfang)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading xfang in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! tminvj:
    status = NF90_INQ_VARID(ncid_rin,'tminvj',ipid%tminvj)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding tminvj parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%tminvj,veg%tminvj)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading tminvj in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! tmaxvj:
    status = NF90_INQ_VARID(ncid_rin,'tmaxvj',ipid%tmaxvj)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding tmaxvj parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%tmaxvj,veg%tmaxvj)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading tmaxvj in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! vbeta:
    status = NF90_INQ_VARID(ncid_rin,'vbeta',ipid%vbeta)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding vbeta parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%vbeta,veg%vbeta)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading vbeta in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ratecp:
    status = NF90_INQ_VARID(ncid_rin,'ratecp',ipid%ratecp)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ratecp parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%ratecp,bgc%ratecp)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ratecp in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! ratecs:
    status = NF90_INQ_VARID(ncid_rin,'ratecs',ipid%ratecs)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding ratecs parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%ratecs,bgc%ratecs)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading ratecs in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! meth:
    status = NF90_INQ_VARID(ncid_rin,'meth',ipid%meth)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding meth parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%meth,veg%meth)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading meth in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')
    ! za:
    status = NF90_INQ_VARID(ncid_rin,'za',ipid%za)
    IF(status /= NF90_NOERR) CALL nc_abort &
         ('Error finding za parameter in restart file ' &
         //TRIM(filename_restart_in)//' (SUBROUTINE get_restart)')
    status=NF90_GET_VAR(ncid_rin,ipid%za,rough%za)            
    IF(status/=NF90_NOERR) CALL nc_abort('Error reading za in file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')

    ! Close restart file:
    status = NF90_CLOSE(ncid_rin)
    IF(status/=NF90_NOERR) CALL nc_abort('Error closing restart file ' &
         //TRIM(filename_restart_in)// '(SUBROUTINE get_restart)')

    ! Establish derivative parameters from those read in:
    soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints:
    soil%zshh(ms+1)=0.5*soil%zse(ms)
    soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
    soil%cnsd  = soil%sand*0.3 + soil%clay*0.25 &
         + soil%silt*0.265 ! set dry soil thermal conductivity [W/m/K]
    soil%hsbh  = soil%hyds*ABS(soil%sucs)*soil%bch        !difsat*etasat
    soil%ibp2  = NINT(soil%bch)+2
    soil%i2bp3 = 2*NINT(soil%bch)+3
    ! Initialise sum flux variables:
    sum_flux%sumpn = 0.0
    sum_flux%sumrp = 0.0
    sum_flux%sumrpw = 0.0
    sum_flux%sumrpr = 0.0
    sum_flux%sumrs = 0.0
    sum_flux%sumrd = 0.0
    sum_flux%dsumpn = 0.0
    sum_flux%dsumrp = 0.0
    sum_flux%dsumrd = 0.0
    ! Initialise conservation variables:
    bal%precip_tot = 0.0
    bal%rnoff_tot = 0.0
    bal%evap_tot = 0.0
    bal%wbal_tot = 0.0
    bal%ebal_tot = 0.0
    bal%drybal = 0.0
    bal%wetbal = 0.0

  END SUBROUTINE get_restart_data
  !============================================================================
  SUBROUTINE load_parameters(filename_restart_in,filename_met,met,air, &
       ssoil,veg,bgc,soil,canopy,rough,rad,sum_flux,bal,logn)
    ! Checks where parameters and initialisations should be loaded from.
    ! If they can be found in either the met file or restart file, they will 
    ! load from there, with the met file taking precedence. Otherwise, they'll
    ! be chosen from a coarse global grid of veg and soil types, based on 
    ! the lat/lon coordinates.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename_restart_in
    CHARACTER(LEN=*), INTENT(IN) :: filename_met
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    INTEGER,INTENT(IN) :: logn     ! log file unit number
    LOGICAL :: completeSet ! was a complete parameter set found?
    INTEGER :: i ! do loop variable

    WRITE(logn,*)
    WRITE(logn,*) 'Looking for parameters and initial states....'  
    ! Unless a restart file is found AND the met file has LAI,
    ! default_params will be called (although the initialisations 
    ! parameter values it writes may be overwritten by restart
    ! file values, if found).

    ! Look for restart file (which will have parameters):
    status = NF90_OPEN(filename_restart_in,0,ncid_rin) ! open restart file
    IF (status /= NF90_NOERR) THEN
       ! If no restart file, use default_parameters to load default
       ! parameters and initialisations (and get default grid coords 
       ! if default LAI is required). Parameter values will be over written 
       ! by any found in the met file.
       WRITE(logn,*) ' Could not find restart file: ',&
            TRIM(filename_restart_in)
       WRITE(logn,*) ' Loading initialisations ', &
            'from default grid. '
       ! Load default parameters and initialisations:
       CALL default_params(met,air,ssoil,veg,bgc,soil,canopy,rough, &
            rad,sum_flux,bal,logn)
       ! Check if there are any parameters in the met file; if so, 
       ! overwrite default values for those that are available:
       CALL get_parameters_met &
            (filename_met,soil,veg,bgc,rough,sum_flux,bal,completeSet)
       ! Results of looking for pars in the met file:
       IF(exists%parameters.AND.completeSet) THEN
          ! All pars were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from met input file: ', TRIM(filename_met)
       ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
          ! Only some pars were found in met file:
          WRITE(logn,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename_met), &
               ' the rest are default values - see below'
          WRITE(*,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename_met), &
               ' the rest are default values - check log file'
       ELSE
          ! No parameters were found in met file:
          WRITE(logn,*) ' Loaded parameters ', &
               'from default grid. '
       END IF
    ELSE ! RESTART FILE EXISTS
       ! If met file does not have LAI, then default_params will
       ! need to be called to load default grid to get LAI in the main
       ! time step loop. The parameters and initialisations it writes 
       ! will be overwritten by the restart (and possibly met file) values.
       IF(.NOT.exists%LAI) THEN
          ! Call default_params to get default grid details and 
          ! allocate CABLE's main variables:
          WRITE(logn,*) ' Loading default grid to use default LAI:'
          CALL default_params(met,air,ssoil,veg,bgc,soil,canopy,rough, &
               rad,sum_flux,bal,logn)
       ELSE ! i.e. LAI in met file
          ! Allocate CABLE's main variables:
          CALL alloc_cbm_var(air, mp)
          CALL alloc_cbm_var(bgc, mp)
          CALL alloc_cbm_var(canopy, mp)
          CALL alloc_cbm_var(met, mp)
          CALL alloc_cbm_var(bal, mp)
          CALL alloc_cbm_var(rad, mp)
          CALL alloc_cbm_var(rough, mp)
          CALL alloc_cbm_var(soil, mp)
          CALL alloc_cbm_var(ssoil, mp)
          CALL alloc_cbm_var(sum_flux, mp)
          CALL alloc_cbm_var(veg, mp)
       END IF ! if LAI in met file
       ! Restart file exists, parameters and init will be loaded from it.
       WRITE(logn,*) ' Loading initialisations ', &
            'from restart file: ', TRIM(filename_restart_in)
       ! Load initialisations and parameters from restart file:
       CALL get_restart_data(filename_restart_in,filename_met, &
            logn,ssoil,canopy,rough,bgc,bal,veg,&
            soil,rad,sum_flux)
       ! Overwrite any parameters found in met file:
       CALL get_parameters_met &
            (filename_met,soil,veg,bgc,rough,sum_flux,bal,completeSet)
       ! Results of looking for pars in the met file:
       IF(exists%parameters.AND.completeSet) THEN
          ! All pars were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from met input file: ', TRIM(filename_met)
       ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
          ! Only some pars were found in met file:
          WRITE(logn,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename_met), &
               ' the rest are from restart file: ', &
               TRIM(filename_restart_in),' - see below'
          WRITE(*,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename_met), &
               ' the rest are from restart file: ', &
               TRIM(filename_restart_in),' - check log file'
       ELSE
          ! No parameters were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from restart file: ',TRIM(filename_restart_in)
       END IF
    END IF ! if restart file exists
    WRITE(logn,*)

    ! Check for basic inconsistencies in parameter values:
    IF(ANY(soil%sand+soil%silt+soil%clay/=1.0)) THEN
       DO i=1,mp
          IF(soil%sand(i)+soil%silt(i)+soil%clay(i)>1.0001.OR. &
              soil%sand(i)+soil%silt(i)+soil%clay(i)<0.9999) THEN
             WRITE(*,*) 'SUBROUTINE load_parameters:'
             WRITE(*,*) 'At land point number', i
             CALL abort ('clay+sand+silt fraction does not sum to 1!')
          END IF
       END DO
    END IF
    IF(ANY(soil%swilt>soil%sfc).OR.ANY(soil%sfc>soil%ssat)) THEN
       DO i=1,mp
          IF(soil%swilt(i)>soil%sfc(i).OR.soil%sfc(i)>soil%ssat(i)) THEN
             WRITE(*,*) 'SUBROUTINE load_parameters:'
             WRITE(*,*) 'At land point number', i
             CALL abort ('Wilting pt < field capacity < saturation violated!')
          END IF
       END DO
    END IF

    ! Write parameter values to log file if requested:
    IF(verbose) CALL report_parameters(logn,soil,veg,bgc, &
         rough,ssoil,canopy)

  END SUBROUTINE load_parameters
  !============================================================================
  SUBROUTINE get_parameters_met(filename_met,soil,veg,bgc,rough, &
       sum_flux,bal,completeSet)
    ! This subroutine looks for parameters in the met file, and 
    ! loads those that are found.
    CHARACTER(LEN=*), INTENT(IN) :: filename_met
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    LOGICAL, INTENT(OUT) :: completeSet ! were all pars found?
    TYPE (parID_type) :: ipid ! parameter IDs in netcdf met file
    TYPE (parID_type) :: load ! switch to load a param - 1=>load it
    REAL(r_1),DIMENSION(1,1) :: data2 ! temporary for ncdf read in
    INTEGER(i_d),DIMENSION(1,1) :: data2i ! temporary for ncdf read in
    REAL(r_1),DIMENSION(1,1,1) :: data3 ! temporary for ncdf read in
    INTEGER :: i ! do loop variable

    completeSet=.TRUE. ! initialise

    ! Get parameter IDs:
    ! iveg:
    status = NF90_INQ_VARID(ncid_met,'iveg',ipid%iveg)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%iveg=0
    ELSE
       exists%parameters=.TRUE.
       load%iveg=1
    END IF
    ! isoil:
    status = NF90_INQ_VARID(ncid_met,'isoil',ipid%isoil)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%isoil=0
    ELSE
       exists%parameters=.TRUE.
       load%isoil=1
    END IF
    ! clay:
    status = NF90_INQ_VARID(ncid_met,'clay',ipid%clay)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%clay=0
    ELSE
       exists%parameters=.TRUE.
       load%clay=1
    END IF
    ! sand:
    status = NF90_INQ_VARID(ncid_met,'sand',ipid%sand)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%sand=0
    ELSE
       exists%parameters=.TRUE.
       load%sand=1
    END IF
    ! silt:
    status = NF90_INQ_VARID(ncid_met,'silt',ipid%silt)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%silt=0
    ELSE
       exists%parameters=.TRUE.
       load%silt=1
    END IF
    ! ssat:
    status = NF90_INQ_VARID(ncid_met,'ssat',ipid%ssat)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%ssat=0
    ELSE
       exists%parameters=.TRUE.
       load%ssat=1
    END IF
    ! sfc:
    status = NF90_INQ_VARID(ncid_met,'sfc',ipid%sfc)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%sfc=0
    ELSE
       exists%parameters=.TRUE.
       load%sfc=1
    END IF
    ! swilt:
    status = NF90_INQ_VARID(ncid_met,'swilt',ipid%swilt)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%swilt=0
    ELSE
       exists%parameters=.TRUE.
       load%swilt=1
    END IF
    ! zse:
    status = NF90_INQ_VARID(ncid_met,'zse',ipid%zse)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%zse=0
    ELSE
       exists%parameters=.TRUE.
       load%zse=1
    END IF
    ! froot:
    status = NF90_INQ_VARID(ncid_met,'froot',ipid%froot)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%froot=0
    ELSE
       exists%parameters=.TRUE.
       load%froot=1
    END IF
    ! bch:
    status = NF90_INQ_VARID(ncid_met,'bch',ipid%bch)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%bch=0
    ELSE
       exists%parameters=.TRUE.
       load%bch=1
    END IF
    ! hyds:
    status = NF90_INQ_VARID(ncid_met,'hyds',ipid%hyds)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%hyds=0
    ELSE
       exists%parameters=.TRUE.
       load%hyds=1
    END IF
    ! sucs:
    status = NF90_INQ_VARID(ncid_met,'sucs',ipid%sucs)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%sucs=0
    ELSE
       exists%parameters=.TRUE.
       load%sucs=1
    END IF
    ! css:
    status = NF90_INQ_VARID(ncid_met,'css',ipid%css)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%css=0
    ELSE
       exists%parameters=.TRUE.
       load%css=1
    END IF
    ! rhosoil:
    status = NF90_INQ_VARID(ncid_met,'rhosoil',ipid%rhosoil)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%rhosoil=0
    ELSE
       exists%parameters=.TRUE.
       load%rhosoil=1
    END IF
    ! rs20:
    status = NF90_INQ_VARID(ncid_met,'rs20',ipid%rs20)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%rs20=0
    ELSE
       exists%parameters=.TRUE.
       load%rs20=1
    END IF
    ! albsoil:
    status = NF90_INQ_VARID(ncid_met,'albsoil',ipid%albsoil)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%albsoil=0
    ELSE
       exists%parameters=.TRUE.
       load%albsoil=1
    END IF
    ! hc:
    status = NF90_INQ_VARID(ncid_met,'hc',ipid%hc)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%hc=0
    ELSE
       exists%parameters=.TRUE.
       load%hc=1
    END IF
    ! canst1:
    status = NF90_INQ_VARID(ncid_met,'canst1',ipid%canst1)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%canst1=0
    ELSE
       exists%parameters=.TRUE.
       load%canst1=1
    END IF
    ! dleaf:
    status = NF90_INQ_VARID(ncid_met,'dleaf',ipid%dleaf)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%dleaf=0
    ELSE
       exists%parameters=.TRUE.
       load%dleaf=1
    END IF
    ! frac4:
    status = NF90_INQ_VARID(ncid_met,'frac4',ipid%frac4)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%frac4=0
    ELSE
       exists%parameters=.TRUE.
       load%frac4=1
    END IF
    ! ejmax:
    status = NF90_INQ_VARID(ncid_met,'ejmax',ipid%ejmax)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%ejmax=0
    ELSE
       exists%parameters=.TRUE.
       load%ejmax=1
    END IF
    ! vcmax:
    status = NF90_INQ_VARID(ncid_met,'vcmax',ipid%vcmax)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%vcmax=0
    ELSE
       exists%parameters=.TRUE.
       load%vcmax=1
    END IF
    ! rp20:
    status = NF90_INQ_VARID(ncid_met,'rp20',ipid%rp20)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%rp20=0
    ELSE
       exists%parameters=.TRUE.
       load%rp20=1
    END IF
    ! rpcoef:
    status = NF90_INQ_VARID(ncid_met,'rpcoef',ipid%rpcoef)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%rpcoef=0
    ELSE
       exists%parameters=.TRUE.
       load%rpcoef=1
    END IF
    ! shelrb:
    status = NF90_INQ_VARID(ncid_met,'shelrb',ipid%shelrb)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%shelrb=0
    ELSE
       exists%parameters=.TRUE.
       load%shelrb=1
    END IF
    ! xfang:
    status = NF90_INQ_VARID(ncid_met,'xfang',ipid%xfang)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%xfang=0
    ELSE
       exists%parameters=.TRUE.
       load%xfang=1
    END IF
    ! tminvj:
    status = NF90_INQ_VARID(ncid_met,'tminvj',ipid%tminvj)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%tminvj=0
    ELSE
       exists%parameters=.TRUE.
       load%tminvj=1
    END IF
    ! tmaxvj:
    status = NF90_INQ_VARID(ncid_met,'tmaxvj',ipid%tmaxvj)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%tmaxvj=0
    ELSE
       exists%parameters=.TRUE.
       load%tmaxvj=1
    END IF
    ! vbeta:
    status = NF90_INQ_VARID(ncid_met,'vbeta',ipid%vbeta)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%vbeta=0
    ELSE
       exists%parameters=.TRUE.
       load%vbeta=1
    END IF
    ! ratecp:
    status = NF90_INQ_VARID(ncid_met,'ratecp',ipid%ratecp)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%ratecp=0
    ELSE
       exists%parameters=.TRUE.
       load%ratecp=1
    END IF
    ! ratecs:
    status = NF90_INQ_VARID(ncid_met,'ratecs',ipid%ratecs)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%ratecs=0
    ELSE
       exists%parameters=.TRUE.
       load%ratecs=1
    END IF
    ! meth:
    status = NF90_INQ_VARID(ncid_met,'meth',ipid%meth)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%meth=0
    ELSE
       exists%parameters=.TRUE.
       load%meth=1
    END IF
    ! za:
    status = NF90_INQ_VARID(ncid_met,'za',ipid%za)
    IF(status /= NF90_NOERR) THEN
       completeSet=.FALSE.
       load%za=0
    ELSE
       exists%parameters=.TRUE.
       load%za=1
    END IF

    ! Now fetch parameter values:
    IF(gridType=='mask') THEN
       DO i=1,mp ! over all land points/grid cells
          ! Get data from land/sea mask type grid:
          IF(load%iveg==1) THEN
             ! Get iveg data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%iveg,data2i, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading iveg in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%iveg(i)=data2i(1,1)
          END IF
          IF(load%isoil==1) THEN
             ! Get isoil data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%isoil,data2i, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading isoil in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%isoilm(i)=data2i(1,1)
          END IF
          IF(load%clay==1) THEN
             ! Get clay data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%clay,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading clay in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%clay(i)=data2(1,1)
          END IF
          IF(load%sand==1) THEN
             ! Get sand data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%sand,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading sand in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%sand(i)=data2(1,1)
          END IF
          IF(load%silt==1) THEN
             ! Get silt data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%silt,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading silt in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%silt(i)=data2(1,1)
          END IF
          IF(load%ssat==1) THEN
             ! Get ssat data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%ssat,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading ssat in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%ssat(i)=data2(1,1)
          END IF
          IF(load%sfc==1) THEN
             ! Get sfc data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%sfc,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading sfc in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%sfc(i)=data2(1,1)
          END IF
          IF(load%swilt==1) THEN
             ! Get swilt data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%swilt,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading swilt in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%swilt(i)=data2(1,1)
          END IF
          IF(load%zse==1) THEN
             ! Get zse data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%zse,soil%zse, &
                  start=(/1/),count=(/ms/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading zse in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
          END IF
          IF(load%froot==1) THEN
             ! Get froot data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%froot,data3, &
                  start=(/land_x(i),land_y(i),1/),count=(/1,1,ms/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading froot in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%froot(i,:)=data3(1,1,:)
          END IF
          IF(load%bch==1) THEN
             ! Get bch data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%bch,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading bch in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%bch(i)=data2(1,1)
          END IF
          IF(load%hyds==1) THEN
             ! Get hyds data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%hyds,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading hyds in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%hyds(i)=data2(1,1)
          END IF
          IF(load%sucs==1) THEN
             ! Get sucs data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%sucs,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading sucs in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%sucs(i)=data2(1,1)
          END IF
          IF(load%css==1) THEN
             ! Get css data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%css,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading css in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%css(i)=data2(1,1)
          END IF
          IF(load%rhosoil==1) THEN
             ! Get rhosoil data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%rhosoil,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading rhosoil in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%rhosoil(i)=data2(1,1)
          END IF
          IF(load%rs20==1) THEN
             ! Get rs20 data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%rs20,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading rs20 in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%rs20(i)=data2(1,1)
          END IF
          IF(load%albsoil==1) THEN
             ! Get albsoil data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%albsoil,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading albsoil in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             soil%albsoil(i)=data2(1,1)
          END IF
          IF(load%hc==1) THEN
             ! Get hc data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%hc,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading hc in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%hc(i)=data2(1,1)
          END IF
          IF(load%canst1==1) THEN
             ! Get canst1 data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%canst1,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading canst1 in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%canst1(i)=data2(1,1)
          END IF
          IF(load%dleaf==1) THEN
             ! Get dleaf data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%dleaf,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading dleaf in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%dleaf(i)=data2(1,1)
          END IF
          IF(load%frac4==1) THEN
             ! Get frac4 data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%frac4,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading frac4 in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%frac4(i)=data2(1,1)
          END IF
          IF(load%ejmax==1) THEN
             ! Get ejmax data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%ejmax,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading ejmax in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%ejmax(i)=data2(1,1)
          END IF
          IF(load%vcmax==1) THEN
             ! Get vcmax data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%vcmax,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading vcmax in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%vcmax(i)=data2(1,1)
          END IF
          IF(load%rp20==1) THEN
             ! Get rp20 data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%rp20,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading rp20 in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%rp20(i)=data2(1,1)
          END IF
          IF(load%rpcoef==1) THEN
             ! Get rpcoef data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%rpcoef,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading rpcoef in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%rpcoef(i)=data2(1,1)
          END IF
          IF(load%shelrb==1) THEN
             ! Get shelrb data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%shelrb,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading shelrb in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%shelrb(i)=data2(1,1)
          END IF
          IF(load%xfang==1) THEN
             ! Get xfang data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%xfang,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading xfang in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%xfang(i)=data2(1,1)
          END IF
          IF(load%tminvj==1) THEN
             ! Get tminvj data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%tminvj,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading tminvj in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%tminvj(i)=data2(1,1)
          END IF
          IF(load%tmaxvj==1) THEN
             ! Get tmaxvj data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%tmaxvj,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading tmaxvj in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%tmaxvj(i)=data2(1,1)
          END IF
          IF(load%vbeta==1) THEN
             ! Get vbeta data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%vbeta,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading vbeta in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%vbeta(i)=data2(1,1)
          END IF
          IF(load%ratecp==1) THEN
             ! Get ratecp data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%ratecp,bgc%ratecp, &
                  start=(/1/),count=(/ncp/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading ratecp in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
          END IF
          IF(load%ratecs==1) THEN
             ! Get ratecs data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%ratecs,bgc%ratecs, &
                  start=(/1/),count=(/ncs/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading ratecs in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
          END IF
          IF(load%meth==1) THEN
             ! Get meth data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%meth,data2i, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading meth in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
             veg%meth(i)=data2i(1,1)
          END IF
          IF(load%za==1) THEN
             ! Get za data:- - - - - - - - - - - - - - - - - - - -
             status= NF90_GET_VAR(ncid_met,ipid%za,data2, &
                  start=(/land_x(i),land_y(i)/),count=(/1,1/))
             IF(status /= NF90_NOERR) CALL nc_abort &
                  ('Error reading za in met data file ' &
                  //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)') 
             rough%za(i)=data2(1,1)
          END IF
       END DO
    ELSE IF(gridType=='land') THEN
       ! Collect data from land only grid in netcdf file:
       IF(load%iveg==1) THEN
          ! Get iveg data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%iveg,veg%iveg)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading iveg in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%isoil==1) THEN
          ! Get isoil data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%isoil,soil%isoilm)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading isoil in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%clay==1) THEN
          ! Get clay data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%clay,soil%clay)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading clay in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%sand==1) THEN
          ! Get sand data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%sand,soil%sand)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading sand in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%silt==1) THEN
          ! Get silt data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%silt,soil%silt)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading silt in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%ssat==1) THEN
          ! Get ssat data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%ssat,soil%ssat)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading ssat in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%sfc==1) THEN
          ! Get sfc data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%sfc,soil%sfc)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading sfc in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%swilt==1) THEN
          ! Get swilt data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%swilt,soil%swilt)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading swilt in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%zse==1) THEN
          ! Get zse data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%zse,soil%zse)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading zse in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%froot==1) THEN
          ! Get froot data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%froot,veg%froot)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading froot in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%bch==1) THEN
          ! Get bch data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%bch,soil%bch)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading bch in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%hyds==1) THEN
          ! Get hyds data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%hyds,soil%hyds)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading hyds in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%sucs==1) THEN
          ! Get sucs data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%sucs,soil%sucs)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading sucs in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%css==1) THEN
          ! Get css data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%css,soil%css)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading css in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%rhosoil==1) THEN
          ! Get rhosoil data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%rhosoil,soil%rhosoil)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading rhosoil in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%rs20==1) THEN
          ! Get rs20 data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%rs20,soil%rs20)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading rs20 in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%albsoil==1) THEN
          ! Get albsoil data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%albsoil,soil%albsoil)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading albsoil in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%hc==1) THEN
          ! Get hc data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%hc,veg%hc)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading hc in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%canst1==1) THEN
          ! Get canst1 data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%canst1,veg%canst1)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading canst1 in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%dleaf==1) THEN
          ! Get dleaf data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%dleaf,veg%dleaf)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading dleaf in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%frac4==1) THEN
          ! Get frac4 data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%frac4,veg%frac4)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading frac4 in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%ejmax==1) THEN
          ! Get ejmax data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%ejmax,veg%ejmax)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading ejmax in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%vcmax==1) THEN
          ! Get vcmax data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%vcmax,veg%vcmax)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading vcmax in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%rp20==1) THEN
          ! Get rp20 data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%rp20,veg%rp20)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading rp20 in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%rpcoef==1) THEN
          ! Get rpcoef data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%rpcoef,veg%rpcoef)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading rpcoef in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%shelrb==1) THEN
          ! Get shelrb data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%shelrb,veg%shelrb)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading shelrb in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%xfang==1) THEN
          ! Get xfang data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%xfang,veg%xfang)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading xfang in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%tminvj==1) THEN
          ! Get tminvj data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%tminvj,veg%tminvj)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading tminvj in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%tmaxvj==1) THEN
          ! Get tmaxvj data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%tmaxvj,veg%tmaxvj)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading tmaxvj in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%vbeta==1) THEN
          ! Get vbeta data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%vbeta,veg%vbeta)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading vbeta in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%ratecp==1) THEN
          status= NF90_GET_VAR(ncid_met,ipid%ratecp,bgc%ratecp)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading ratecp in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%ratecs==1) THEN
          ! Get ratecs data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%ratecs,bgc%ratecs)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading ratecs in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%meth==1) THEN
          ! Get meth data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%meth,veg%meth)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading meth in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
       IF(load%za==1) THEN
          ! Get za data:- - - - - - - - - - - - - - - - - - - -
          status= NF90_GET_VAR(ncid_met,ipid%za,rough%za)
          IF(status /= NF90_NOERR) CALL nc_abort &
               ('Error reading za in met data file ' &
               //TRIM(filename_met)//' (SUBROUTINE get_parameters_met)')
       END IF
    END IF ! grid type mask or land-only

    ! Establish derivative parameters from those read in:
    soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints:
    soil%zshh(ms+1)=0.5*soil%zse(ms)
    soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
    soil%cnsd  = soil%sand*0.3 + soil%clay*0.25 &
         + soil%silt*0.265 ! set dry soil thermal conductivity [W/m/K]
    soil%hsbh  = soil%hyds*ABS(soil%sucs)*soil%bch        !difsat*etasat
    soil%ibp2  = NINT(soil%bch)+2
    soil%i2bp3 = 2*NINT(soil%bch)+3
    ! Initialise sum flux variables:
    sum_flux%sumpn = 0.0
    sum_flux%sumrp = 0.0
    sum_flux%sumrpw = 0.0
    sum_flux%sumrpr = 0.0
    sum_flux%sumrs = 0.0
    sum_flux%sumrd = 0.0
    sum_flux%dsumpn = 0.0
    sum_flux%dsumrp = 0.0
    sum_flux%dsumrd = 0.0
    ! Initialise conservation variables:
    bal%precip_tot = 0.0
    bal%rnoff_tot = 0.0
    bal%evap_tot = 0.0
    bal%wbal_tot = 0.0
    bal%ebal_tot = 0.0
    bal%drybal = 0.0
    bal%wetbal = 0.0

  END SUBROUTINE get_parameters_met
  !============================================================================
  SUBROUTINE nc_abort(message)
    ! Error subroutine in case of fatal error:
    CHARACTER(LEN=*), INTENT(IN) :: message
    WRITE(*,*) message ! error from subroutine
    WRITE(*,*) NF90_STRERROR(status) ! netcdf error details
    STOP
  END SUBROUTINE nc_abort
  !==========================================================================
END MODULE input_module
