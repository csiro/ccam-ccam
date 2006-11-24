!!$ Sample text offline driver for CABLE land surface scheme, 2006.
!!$ Gab Abramowitz, CSIRO Marine and Atmospheric Research/
!!$ Macquarie University; gabsun@gmail.com 

PROGRAM offline_driver
  USE cbm_module
  USE canopy_module
  USE physical_constants
  USE parameter_module
  IMPLICIT NONE
  INTEGER(i_d)		:: kend ! no. of time steps in run
  TYPE (air_type)	:: air
  TYPE (bgc_pool_type)	:: bgc	
  TYPE (canopy_type)	:: canopy
  TYPE (met_type) 	:: met
  TYPE (balances_type)  :: bal
  TYPE (radiation_type) :: rad
  TYPE (roughness_type) :: rough
  TYPE (soil_parameter_type) :: soil	
  TYPE (soil_snow_type)	:: ssoil
  TYPE (sum_flux_type)	:: sum_flux
  TYPE (veg_parameter_type)	:: veg	
  REAL(r_1)	        :: dels ! time step size in seconds
  REAL(r_1),ALLOCATABLE,DIMENSION(:) :: latitude, longitude
  INTEGER(i_d),ALLOCATABLE,DIMENSION(:) :: gdpt ! gridpoint number (default params)
  INTEGER(i_d) :: logn     ! log file unit number
  INTEGER(i_d) 	:: kstart   ! start of simulation #
  INTEGER(i_d)  :: ktau	  ! index of time step = 1 ..  kend
  INTEGER(i_d)  :: ktauyear ! CCAM parameter; ignore
  CHARACTER(LEN=99) :: filename_met ! name of file for met. parameters
  CHARACTER(LEN=99) :: filename_out ! name of file for CABLE output
  CHARACTER(LEN=99) :: filename_log ! name of file for execution log
  INTEGER(i_d) :: u=30, k_in, ios,j
  LOGICAL :: is_open
  REAL(r_1) :: xdum1, xdum2

  !======================!
  ! Number of grid points:
  mp = 1 
  !======================!
  ! Allocate latitude, longitude 
  ALLOCATE(latitude(mp),longitude(mp))
  filename_log    = 'log_cable.txt'
  filename_out    = 'out_cable.txt'
  logn = 88 ! log file unit number
 
  !=============== SINGLE SITE EXAMPLES ===============
  ! Tumbarumba (Australian eucalypt site):
  filename_met    = 'sample_met/Tumbarumba_met.txt'
  dels = 3600.0       ! set time step size in seconds
  kend = 17520         ! set number of time steps
  latitude = -35.6557 ! set latitude
  longitude = 148.152 ! set longitude
!!$  ! Tharandt (German coniferous site):
!!$  filename_met    = 'sample_met/Tharandt_met.txt'
!!$  dels = 1800.0       ! set time step size in seconds
!!$  kend = 87600        ! set number of time steps 
!!$  latitude = 50.9670  ! set latitude
!!$  longitude = 13.633  ! set longitude
!!$  ! Bondville (USA crop site):
!!$  filename_met    = 'sample_met/Bondville_met.txt'
!!$  dels = 1800.0       ! set time step size in seconds
!!$  kend = 52560        ! set number of time steps 
!!$  latitude = 40.006   ! set latitude
!!$  longitude = -88.292 ! set longitude
  !========================================================

  ! Open log file:
  OPEN(logn,FILE=filename_log)

  ! Get default land surface parameters and allocate main variables:
  CALL default_params(latitude,longitude,met,air,ssoil,veg,bgc,soil,canopy, &
       rough,rad,sum_flux,bal,gdpt,logn)    

  !=====================================================================!
  ! user-forced initialisations/parameter values (see user guide):
  
  
  !=====================================================================!

  kstart = 1
  DO ktau = kstart, kend ! time step loop
    
     ! Get met data:
     INQUIRE (u, opened=is_open)
     IF (.not. is_open) THEN
        OPEN(u, file=filename_met, iostat=ios)
     END IF
     ! read met data for one site:
     READ(u, *) k_in, met%doy(1), met%hod(1), met%fsd(1), met%fld(1), &
          met%precip(1), met%tc(1), met%ua(1), met%qv(1), met%pmb(1), &
          met%ca(1), rad%latitude(1),  &
          met%coszen(1), xdum1,xdum2,veg%vlai(1)
     
     ! Initialise other met variables:
     met%precips = 0.0
     WHERE(met%tc <= 0.0) met%precips = met%precip
     met%tk = met%tc + tfrz
     met%tvair = met%tk 
     met%tvrad = met%tk 
     met%coszen = sinbet(met%doy, rad%latitude, met%hod)
   
     ! CALL land surface scheme for this timestep
     CALL cbm(ktau,kstart,kend,ktauyear,dels,air,bgc,canopy,met, &
          bal,rad,rough,soil,ssoil,sum_flux,veg)

  !   CALL text_output(ktau,kstart,kend,ktauyear,dels,air,bgc, &
  !        canopy,met,bal,rad,rough,soil,ssoil,sum_flux,veg, &
  !        filename_out)
  END DO

  print*, ssoil%wb
  print*, ssoil%tgg

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

END PROGRAM offline_driver
