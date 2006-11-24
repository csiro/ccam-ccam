!!$ Netcdf offline driver for CABLE land surface scheme, 2006.
!!$ Gab Abramowitz, CSIRO Marine and Atmospheric Research/
!!$ Macquarie University; gabsun@gmail.com 

PROGRAM offline_driver
  USE cbm_module
  USE output_module
  USE parameter_module
  IMPLICIT NONE
  INTEGER(i_d)		:: kend ! no. of time steps in run
  TYPE (air_type)	:: air  ! air property variables
  TYPE (bgc_pool_type)	:: bgc	! carbon pool variables
  TYPE (canopy_type)	:: canopy ! vegetation variables
  TYPE (met_type) 	:: met  ! met input variables
  TYPE (balances_type)  :: bal  ! energy and water balance variables
  TYPE (radiation_type) :: rad  ! radiation variables
  TYPE (roughness_type) :: rough ! roughness varibles
  TYPE (soil_parameter_type) :: soil ! soil parameters	
  TYPE (soil_snow_type)	:: ssoil ! soil and snow variables
  TYPE (sum_flux_type)	:: sum_flux ! cumulative flux variables
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters	
  REAL(r_1)	        :: dels ! time step size in seconds
  REAL(r_1), ALLOCATABLE,DIMENSION(:) :: latitude, longitude
  INTEGER(i_d),ALLOCATABLE,DIMENSION(:) :: gdpt ! gridpoint number (default params)
  INTEGER(i_d) 	:: kstart ! start of simulation #
  INTEGER(i_d)  :: ktau	  ! index of time step = 1 ..  kend
  INTEGER(i_d)  :: ktauyear ! CCAM parameter; ignore
  CHARACTER(LEN=99) :: filename_params ! name of file for soil parameters
  CHARACTER(LEN=99) :: filename_met ! name of file for met. parameters
  CHARACTER(LEN=99) :: filename_out ! name of file for CABLE output
  CHARACTER(LEN=99) :: filename_log ! name of file for execution log

  !===================================================================!
  ! Filenames:
  filename_params = 'default'
  filename_met    = './sample_met/Tumbarumba_met.nc'
  !  filename_met    = './sample_met/Tharandt_met.nc'
  !  filename_met    = './sample_met/Bondville_met.nc'
  !  filename_met    = './sample_met/multi_site.nc' 
  filename_out    = 'out_cable.nc'
  filename_log    = 'log_cable.txt'
  
  ! Which variables should be in the output file?
  output%met = .TRUE.      ! input met data
  output%flux = .TRUE.     ! convective, runoff, NEE
  output%soil = .TRUE.     ! soil states
  output%radiation = .TRUE.! net rad, albedo
  output%carbon = .TRUE.   ! NEE, GPP, NPP, stores 
  output%veg = .TRUE.      ! vegetation states
  output%params = .TRUE.   ! input parameters used to produce run
  output%balances = .TRUE. ! energy and water balances

  ! Which checks to perform?
  check%ranges = .TRUE.     ! variable ranges, input and output
  check%energy_bal = .TRUE. ! energy balance
  check%mass_bal = .TRUE.   ! water/mass balance

  logn = 88 ! log file number - declared in input module
  leaps=.FALSE. ! calculate timing with leap years?
  fixedCO2 = 350.0 ! if not found in met file, in ppmv

  !=====================================================================!
  ! Open log file:
  OPEN(logn,FILE=filename_log)
  
  ! Open met data and get site information from netcdf file.
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites.
  CALL open_met_file(filename_met,dels,kend,latitude,longitude)

  ! Get land surface parameters and allocate main variables:
  IF(filename_params=='default') CALL default_params &
       (latitude,longitude,met,air,ssoil,veg,bgc,soil,canopy,rough, &
       rad,sum_flux,bal,gdpt,logn)    

  ! Open output file:
  CALL open_output_file(filename_out,dels,latitude,longitude, &
       soil,veg,bgc,rough)

  kstart = 1
  DO ktau = kstart, kend ! time step loop
    
     ! Get met data, set time variables:
     CALL get_met_data(ktau,filename_met,met,soil,rad,longitude,dels)
    
     ! Get this time step's LAI:
     CALL get_lai(ktau,met%doy,veg,kend,gdpt)

     ! CALL land surface scheme for this timestep, all grid points:
     CALL cbm(ktau, kstart, kend, ktauyear, dels, air, bgc, canopy, met, &
          bal, rad, rough, soil, ssoil, sum_flux, veg)

     ! Write time step's output to file:
     CALL write_output(ktau,dels,filename_out,met,canopy,ssoil, &
          rad,bal,air,soil,veg)

  END DO
  
  ! Close met data input file:
  CALL close_met_file(filename_met)
  ! Close output file and deallocate main variables:
  CALL close_output_file(filename_out,bal,latitude,longitude,gdpt,air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)
  
END PROGRAM offline_driver

