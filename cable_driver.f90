!!$ Netcdf offline driver for CABLE land surface scheme, 2007.
!!$ Gab Abramowitz, University of New South Wales/
!!$ CSIRO Marine and Atmospheric Research; gabsun@gmail.com 

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
  INTEGER(i_d) 	:: kstart ! start of simulation #
  INTEGER(i_d)  :: ktau	  ! index of time step = 1 ..  kend
  CHARACTER(LEN=99) :: filename_restart_in ! name of restart file to read
  CHARACTER(LEN=99) :: filename_restart_out ! name of restart file to read
  CHARACTER(LEN=99) :: filename_met ! name of file for met. parameters
  CHARACTER(LEN=99) :: filename_out ! name of file for CABLE output
  CHARACTER(LEN=99) :: filename_log ! name of file for execution log
  LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
  LOGICAL    :: spinConv ! has spinup converged?
  REAL(r_1)  :: delsoilM ! allowed variation in soil moisture for spin up
  REAL(r_1)  :: delsoilT ! allowed variation in soil temperature for spin up
  REAL(r_1),POINTER  :: soilMtemp(:,:) ! temporary storage for spin up
  REAL(r_1),POINTER  :: soilTtemp(:,:) ! temporary storage for spin up
  INTEGER(i_d) :: tstep  ! time step counter for spinup
  
  !===================================================================!
  ! Filenames:
  filename_met    = './sample_met/Tumbarumba.nc'
  ! filename_met    = './sample_met/Tharandt.nc'
  ! filename_met    = './sample_met/Bondville.nc'
  filename_restart_in = './restart_Tumbarumba.nc' ! will use defaults if not found
  filename_restart_out = './restart_cableOut.nc'
  filename_out    = 'out_cable.nc'
  filename_log    = 'log_cable.txt'
  
  ! Spin up details (not recommended for regional simulations);
  ! currently uses soil moisture and temperature only:
  spinup=.FALSE.     ! do we spin up the model?
  delsoilM=0.001    ! allowed variation in soil moisture for spin up
  delsoilT=0.01     ! allowed variation in soil temperature for spin up
  
  ! Which groups of variables should be written out?
  output%restart = .TRUE.   ! should a restart file be created?
  output%met = .TRUE.       ! input met data
  output%flux = .TRUE.      ! convective, runoff, NEE
  output%soil = .TRUE.      ! soil states
  output%snow = .TRUE.      ! snow states
  output%radiation = .TRUE. ! net rad, albedo
  output%carbon = .TRUE.    ! NEE, GPP, NPP, stores 
  output%veg = .TRUE.       ! vegetation states
  output%params = .TRUE.    ! input parameters used to produce run
  output%balances = .TRUE.  ! energy and water balances
  ! And any others to be included individually (see user guide for complete list):
  output%NEE = .TRUE.
  output%Qle = .TRUE.
  output%Qh = .TRUE.
  output%SWnet = .TRUE.
  output%LWnet = .TRUE.
  output%Qg = .TRUE.

  ! Which checks to perform?
  check%ranges = .TRUE.     ! variable ranges, input and output
  check%energy_bal = .TRUE. ! energy balance
  check%mass_bal = .TRUE.   ! water/mass balance
  
  verbose=.TRUE. ! write details of every grid cell init and params to log?
  leaps=.FALSE.  ! calculate timing with leap years?
  fixedCO2 = 350.0 ! if not found in met file, in ppmv
  logn = 88 ! log file output device number - declared in input module
  
  !=====================================================================!
  ! Open log file:
  OPEN(logn,FILE=filename_log)
  
  ! Open met data and get site information from netcdf file.
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites. 
  CALL open_met_file(filename_met,dels,kend)

  ! Checks where parameters and initialisations should be loaded from.
  ! If they can be found in either the met file or restart file, they will 
  ! load from there, with the met file taking precedence. Otherwise, they'll
  ! be chosen from a coarse global grid of veg and soil types, based on 
  ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
  CALL load_parameters(filename_restart_in,filename_met,met,air,ssoil, &
       veg,bgc,soil,canopy,rough,rad,sum_flux,bal,logn)

  ! Open output file:
  CALL open_output_file(filename_out,filename_met,dels,soil,veg,bgc,rough)

  kstart = 1
  tstep = 0          ! initialise
  spinConv = .FALSE. ! initialise
  ! spinup loop:
  DO
     ! time step loop:
     DO ktau = kstart, kend ! time step loop
        ! increment total timstep counter
        tstep = tstep + 1

        ! Get met data and LAI, set time variables:
        CALL get_met_data(ktau,filename_met,met,soil,rad,veg,kend,dels) 
        
        ! CALL land surface scheme for this timestep, all grid points:
        CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
             bal, rad, rough, soil, ssoil, sum_flux, veg)
        
        ! Write time step's output to file if either: we're not spinning up 
        ! or we're spinning up and the spinup has converged:
        IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
             (ktau,dels,filename_out,met,canopy,ssoil,rad,bal,air,soil,veg)
     END DO
     ! see if spinup (if conducting one) has converged:
     IF(spinup.AND..NOT.spinConv) THEN
        ! Write to screen and log file:
        WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
             ' of data set complete...'
        WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
             ' of data set complete...'
        ! IF not 1st run through whole dataset:
        IF(INT(tstep/kend)>1) THEN 
           ! evaluate spinup
           IF(ANY(ABS(ssoil%wb-soilMtemp)>delsoilM).OR. &
                ANY(ABS(ssoil%tgg-soilTtemp)>delsoilT)) THEN
              ! no convergence yet
           ELSE ! spinup has converged
              spinConv = .TRUE.
              ! Write to screen and log file:
              WRITE(*,'(A33)') ' Spinup has converged - final run'
              WRITE(logn,'(A52)') &
                   ' Spinup has converged - final run - writing all data'
              WRITE(logn,'(A37,F7.5,A28)') &
                   ' Criteria: Change in soil moisture < ', &
                   delsoilM, ' in any layer over whole run'
              WRITE(logn,'(A40,F7.5,A28)' ) & 
                   '           Change in soil temperature < ', &
                   delsoilT, ' in any layer over whole run'
           END IF
        ELSE ! allocate variables for storage
           ALLOCATE(soilMtemp(mp,ms),soilTtemp(mp,ms))
        END IF
        ! store soil moisture and temperature
        soilTtemp = ssoil%tgg
        soilMtemp = REAL(ssoil%wb,r_1)
     ELSE
        ! if not spinning up, or spin up has converged, exit:
        EXIT
     END IF
  END DO

  ! Write restart file if requested:
  IF(output%restart) CALL create_restart(filename_restart_out, &
  filename_met,logn,kstart,kend,soil,veg,ssoil,canopy,rough,bgc,bal)

  ! Close met data input file:
  CALL close_met_file(filename_met)
  ! Close output file and deallocate main variables:
  CALL close_output_file(filename_out, bal, air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)
  
END PROGRAM offline_driver

