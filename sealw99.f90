                           module sealw99_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  This code provides the core functionality of FMS longwave
!  radiation. It is based on exchange method with prescribed
!  coefficients embedded in the code.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!
!  shared modules:

!  radiation package shared modules:

use rad_utilities_mod,   only: rad_utilities_init, Lw_control, &
                               cldrad_properties_type, &
                               cld_specification_type, &
                               lw_output_type, longwave_tables1_type,  &
                               longwave_tables2_type,  &
                               longwave_tables3_type, atmos_input_type,&
                               Cldrad_control, radiative_gases_type, &
                               aerosol_type, aerosol_properties_type, &
                               aerosol_diagnostics_type, &
                               optical_path_type, gas_tf_type, &
                               lw_table_type, Lw_parameters, &
                               lw_diagnostics_type, lw_clouds_type, &
                               locate_in_table, looktab, mass_1,  &
                               temp_1, Rad_control, time_type
use longwave_params_mod, only: longwave_params_init, NBCO215, &
                               NBLY_CKD, NBLY_RSB, longwave_params_end

! radiation package modules:

use longwave_clouds_mod, only: longwave_clouds_init, cldtau, cloud,   &
                               lw_clouds_dealloc, longwave_clouds_end, &
                               thickcld 
use longwave_fluxes_mod, only: longwave_fluxes_ks,    &
                               longwave_fluxes_init, &
                               longwave_fluxes_end, &
                               longwave_fluxes_k_down,  &
                               longwave_fluxes_KE_KEp1,  &
                               longwave_fluxes_diag,   &
                               longwave_fluxes_sum
use longwave_tables_mod, only: longwave_tables_init,  &
                               longwave_tables_end
use optical_path_mod,    only: optical_path_setup, &
                               optical_path_init,  &
                               optical_trans_funct_from_KS, &
                               optical_trans_funct_k_down, &
                               optical_trans_funct_KE, &
                               optical_trans_funct_diag, &
                               get_totvo2, get_totch2o, get_totch2obd,&
                               optical_dealloc, optical_path_end
use gas_tf_mod,          only: co2coef, transcolrow, transcol, &
                               get_control_gas_tf, gas_tf_init, &
                               gas_tf_dealloc, gas_tf_end,   &
                               trans_sfc, trans_nearby
use lw_gases_stdtf_mod,  only: lw_gases_stdtf_init, cfc_indx8,  &
                               cfc_indx8_part, cfc_exact, &
                               lw_gases_stdtf_time_vary, &
                               lw_gases_stdtf_dealloc, co2_lblinterp, &
                               ch4_lblinterp, n2o_lblinterp, &
                               lw_gases_stdtf_end

!------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    sealw99_mod is the internal driver for the 1999 sea longwave 
!    radiation code.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

    character(len=128)  :: version =  '$Id: sealw99.F90,v 18.0.2.1.2.1.2.1 2010/08/30 20:33:33 wfc Exp $'
    character(len=128)  :: tagname =  '$Name: testing $'
    logical, save       :: module_is_initialized = .false.

!---------------------------------------------------------------------
!-------  interfaces --------

public       &
         sealw99_init,  sealw99_time_vary,  sealw99,  &
         sealw99_endts, sealw99_end

private   &
          check_tf_interval, obtain_gas_tfs, &
          sealw99_alloc, &
          cool_to_space_approx,   cool_to_space_exact, &
          e1e290, e290, esfc, enear, &
          co2_source_calc, nlte, co2curt, &
          co2_time_vary, ch4_time_vary, n2o_time_vary


!---------------------------------------------------------------------
!-------- namelist  ---------

logical, save      ::    &
                 do_thick = .false. ! perform "pseudo-convective  
                                    ! adjustment" for maximally 
                                    ! overlapped clouds ?
logical, save      ::    &
             do_lwcldemiss = .true. ! use multiple bands to calculate
                                    ! lw cloud emissivites ? 
logical, save      ::    &
         do_ch4lbltmpint  = .false. ! perform and save intermediate
                                    ! flux calculations for ch4?
logical, save      ::    &
         do_n2olbltmpint  = .false. ! perform and save intermediate
                                    ! flux calculations for n2o?
logical, save      ::   &
                  do_nlte = .false. ! there is a non-local thermodynamic
                                    ! equilibrium region at top 
                                    ! of model ?
character(len=16), save  ::  &
          continuum_form = 'ckd2.1' ! continuum specification; either
                                    ! 'ckd2.1', 'ckd2.4', 'mt_ckd1.0', 
                                    ! 'rsb' or 'none'
character(len=16), save  ::  &
   linecatalog_form = 'hitran_2000' ! line catalog specification; either
                                    ! 'hitran_1992' or 'hitran_2000'
real, save         ::  &
        co2_tf_calc_intrvl = 1.0E6  ! interval between recalculating co2
                                    ! transmission functions, relevant
                                    ! for time-varying co2 cases 
                                    ! [ hours ]
logical, save      ::  &
        calc_co2_tfs_on_first_step = .true. 
                                    ! always calculate co2 tfs on 
                                    ! first time step of job ?
logical, save      ::  &
        use_current_co2_for_tf = .false.  
                                    ! use current co2 mixing ratio for  
                                    ! calculating tfs ?
real, save         ::  &
        co2_tf_time_displacement = 0.0 
                                    ! time displacement from job start 
                                    ! to the point in time where co2 
                                    ! tfs are to be valid -- may be (-),
                                    ! 0.0 or (+); used only when
                                    ! calc_co2_tfs_on_first_step is true
                                    ! [ hours ]
real, save         ::  &
        ch4_tf_calc_intrvl = 1.0E6  ! interval between recalculating ch4
                                    ! transmission functions, relevant
                                    ! for time-varying ch4 cases 
                                    ! [ hours ]
logical, save      ::  &
        calc_ch4_tfs_on_first_step = .true. 
                                    ! always calculate ch4 tfs on 
                                    ! first time step of job ?
logical, save      ::  &
        use_current_ch4_for_tf = .false. 
                                    ! use current ch4 mixing ratio for  
                                    ! calculating tfs ?
real, save         ::  &
        ch4_tf_time_displacement = 0.0 
                                    ! time displacement from job start 
                                    ! to the point in time where ch4 
                                    ! tfs are to be valid -- may be (-),
                                    ! 0.0 or (+); used only when
                                    ! calc_ch4_tfs_on_first_step is true
                                    ! [ hours ]
real, save         ::  &
        n2o_tf_calc_intrvl = 1.0E6  ! interval between recalculating n2o
                                    ! transmission functions, relevant
                                    ! for time-varying n2o cases 
                                    ! [ hours ]
logical, save      ::  &
        calc_n2o_tfs_on_first_step = .true. 
                                    ! always calculate n2o tfs on 
                                    ! first time step of job ?
logical, save      ::  &
        use_current_n2o_for_tf = .false. 
                                    ! use current n2o mixing ratio for  
                                    ! calculating tfs ?
real, save         ::  &
        n2o_tf_time_displacement = 0.0 
                                    ! time displacement from job start 
                                    ! to the point in time where n2o 
                                    ! tfs are to be valid -- may be (-),
                                    ! 0.0 or (+); used only when
                                    ! calc_n2o_tfs_on_first_step is true
                                    ! [ hours ]
integer, save      ::  &
        verbose = 0                 ! verbosity level, ranges from 0
                                    ! (min output) to 5 (max output)
logical, save      ::  &
       calc_co2_tfs_monthly = .false.
logical, save      ::  &
       calc_ch4_tfs_monthly = .false.
logical, save      ::  &
       calc_n2o_tfs_monthly = .false.
integer, save      ::      &
       no_h2o_bands_1200_1400 = 1 ! number of bands in the lw par-
                                  ! ameterization between 1200 and 1400
                                  ! cm (-1); 0 and 1 have been
                                  ! tested. other potentially available
                                  ! values are 2, 4, 10  and 20. 
logical, save      ::      &
    use_bnd1_cldtf_for_h2o_bands = .false. ! the 1200-1400 cm(-1) band
                                             ! uses the same radiative
                                             ! properties as the 0-160,
                                             ! 1400-2200 band. needed 
                                             ! for backward compatibil-
                                             ! ity


!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

real, parameter :: GRAV             = 9.80
real, parameter :: DIFFAC           = 1.660000E+00
real, parameter :: RDGAS            = 287.04 
real, parameter :: KAPPA            = 2./7.  
real, parameter :: CP_AIR           = RDGAS/KAPPA
real, parameter :: SECONDS_PER_DAY  = 8.640000E+04
real, parameter :: RADCON_MKS       = (GRAV/CP_AIR)*SECONDS_PER_DAY
real, parameter :: RADCON           = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY

!---------------------------------------------------------------------
!     apcm, bpcm    capphi coefficients for NBLY bands.
!     atpcm, btpcm  cappsi coefficients for NBLY bands.
!     acomb         random "a" parameter for NBLY bands.
!     bcomb         random "b" parameter for NBLY bands.
!---------------------------------------------------------------------
real, dimension (:), allocatable, save    ::  apcm, bpcm, atpcm, btpcm, &
                                              acomb, bcomb

!-------------------------------------------------------------------
!    the following longwave tables are retained for the life of the
!    run.
!-------------------------------------------------------------------
type (longwave_tables3_type), save       :: tabsr
type (longwave_tables1_type), save       :: tab1, tab2, tab3, tab1w
type (longwave_tables2_type), save       :: tab1a, tab2a, tab3a


!-------------------------------------------------------------------
!
!--------------------------------------------------------------------
integer, parameter                        ::  no_combined_bands = 8
real, dimension(no_combined_bands), save  ::  band_no_start, band_no_end
integer, dimension(NBLY_CKD-1), save      ::  cld_indx_table

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
real, dimension (:),    allocatable, save    ::  c1b7, c2b7
integer, dimension (:), allocatable, save    ::  cld_indx

!---------------------------------------------------------------------
!    miscellaneous variables.
!---------------------------------------------------------------------
integer, save    ::  nbly      ! number of frequency bands for exact
                               ! cool-to-space computations.
integer, save    ::  nbtrge, nbtrg
integer, save    ::  ixprnlte
integer, save    ::  ks, ke
logical, save    ::  do_co2_tf_calc = .true.
logical, save    ::  do_ch4_tf_calc = .true.
logical, save    ::  do_n2o_tf_calc = .true.
logical, save    ::  do_co2_tf_calc_init = .true.
logical, save    ::  do_ch4_tf_calc_init = .true.
logical, save    ::  do_n2o_tf_calc_init = .true.

integer, save    ::  month_of_co2_tf_calc = 0
integer, save    ::  month_of_ch4_tf_calc = 0
integer, save    ::  month_of_n2o_tf_calc = 0

!----------------------------------------------------------------------
!----------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="sealw99_init">
!  <OVERVIEW>
!   Subroutine to initialize longwave radiation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine initializes longwave radiation. It includes the
!   prescribed gas band coefficients, initializes gas optical depth, 
!   longwave tables, and allocate cloud related variables in the
!   longwave spectrum.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sealw99_init (latb, lonb, pref, Lw_tables)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell boundaries [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell boundaries [radians]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   array containing two reference pressure profiles [pascals]
!  </IN>
!  <INOUT NAME="Lw_tables" TYPE="lw_table_type">
!   lw_tables_type variable containing various longwave
!                 table specifiers needed by radiation_diag_mod.
!  </INOUT>
! </SUBROUTINE>
!
subroutine sealw99_init (pref, Lw_tables)
 
!---------------------------------------------------------------------
!    sealw99_init is the constructor for sealw99_mod.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in) :: pref
type(lw_table_type), intent(inout) :: Lw_tables

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb           array of model longitudes on cell boundaries
!                      [ radians ]
!       latb           array of model latitudes at cell boundaries
!                      [ radians ]
!       pref           array containing two reference pressure profiles 
!                      for use in defining transmission functions
!                      [ Pa ]
!
!   intent(inout)    
!
!       Lw_tables      lw_table_type variable which holds much of the
!                      relevant data used in the longwave radiation
!                      parameterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

       real, dimension (no_combined_bands)  :: band_no_start_rsb,   &
                                           band_no_start_ckd,       &
                                             band_no_end_rsb,       &
                                           band_no_end_ckd

       data band_no_start_ckd / 1, 25, 41, 42, 43, 44, 46, 47 /
       data band_no_end_ckd   / 24,40, 41, 42, 43, 45, 46, 47 /

       data band_no_start_rsb / 1, 5, 9, 10, 11, 12, 14, 15 /
       data band_no_end_rsb   / 4, 8, 9, 10, 11, 13, 14, 15 /

       integer, dimension (NBLY_CKD-1) :: cld_indx_table_lwclde, &
                                          cld_indx_table_rsb

       data cld_indx_table_lwclde /40*1, 2, 2, 2, 3, 4, 5, 6 /

       data cld_indx_table_rsb   / 47*1 /

       real, dimension(NBLY_RSB)   :: apcm_n, bpcm_n,      &
                                      atpcm_n, btpcm_n,   &
                                      acomb_n, bcomb_n
       real, dimension(NBLY_CKD) :: apcm_c, bpcm_c, atpcm_c,  &
                                    btpcm_c, acomb_c, bcomb_c 
       real, dimension(size(pref,1) ) :: plm
       real, dimension (NBCO215) :: cent, del

       integer         :: k, n,  nn
       integer         :: ioffset
       real            :: prnlte
       integer         :: kmax, kmin
       character(len=4)  :: gas_name

!---------------------------------------------------------------------
!  local variables:
!
!     band_no_start_rsb
!     band_no_start_ckd
!     band_no_end_rsb
!     band_no_end_ckd
!     cld_indx_table_lwclde
!     cld_indx_table_rsb
!     apcm_n 
!     bpcm_n
!     atpcm_n 
!     btpcm_n
!     acomb_n 
!     bcomb_n
!     apcm_c 
!     bpcm_c 
!     atpcm_c 
!     btpcm_c 
!     acomb_c 
!     bcomb_c 
!     plm
!     cent 
!     del
!     unit 
!     ierr 
!     io 
!     k,n,m,nn
!     ioffset
!     prnlte
!     kmax 
!     kmin
!     inrad
!     dum
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call rad_utilities_init

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kmax = size(pref,1) - 1   ! radiation grid size
      ks = 1
      ke = kmax

!---------------------------------------------------------------------
!    be sure that the radiation time step has been defined before using
!    it.
!---------------------------------------------------------------------
      if (.not.Rad_control%rad_time_step_iz) then
        write(6,*) "ERROR: must define rad_time_step before using it"
        stop
      endif

!---------------------------------------------------------------------
!    call check_tf_interval to verify that the namelist input variables
!    related to co2 are consistent.
!---------------------------------------------------------------------
      gas_name = 'co2 '
      call check_tf_interval (gas_name, co2_tf_calc_intrvl, &
                              calc_co2_tfs_on_first_step,   &
                              calc_co2_tfs_monthly, &
                              use_current_co2_for_tf)

!--------------------------------------------------------------------
!    define the radiation_control_type variable components related to
!    co2 tf calculation.
!--------------------------------------------------------------------
      Rad_control%co2_tf_calc_intrvl = co2_tf_calc_intrvl
      Rad_control%use_current_co2_for_tf = use_current_co2_for_tf
      Rad_control%calc_co2_tfs_monthly       =   &
                                             calc_co2_tfs_monthly      
      Rad_control%calc_co2_tfs_on_first_step =   &
                                             calc_co2_tfs_on_first_step
      Rad_control%co2_tf_time_displacement = co2_tf_time_displacement
      Rad_control%co2_tf_calc_intrvl_iz = .true.             
      Rad_control%use_current_co2_for_tf_iz = .true.
      Rad_control%calc_co2_tfs_on_first_step_iz = .true.
      Rad_control%calc_co2_tfs_monthly_iz    =  .true.
      Rad_control%co2_tf_time_displacement_iz = .true.

!---------------------------------------------------------------------
!    call check_tf_interval to verify that the namelist input variables
!    related to ch4 are consistent.
!---------------------------------------------------------------------
      gas_name = 'ch4 '
      call check_tf_interval (gas_name, ch4_tf_calc_intrvl, &
                              calc_ch4_tfs_on_first_step,   &
                              calc_ch4_tfs_monthly, &
                              use_current_ch4_for_tf)

!--------------------------------------------------------------------
!    define the radiation_control_type variable components related to
!    ch4 tf calculation.
!--------------------------------------------------------------------
      Rad_control%ch4_tf_calc_intrvl = ch4_tf_calc_intrvl
      Rad_control%use_current_ch4_for_tf = use_current_ch4_for_tf
      Rad_control%calc_ch4_tfs_on_first_step =   &
                                             calc_ch4_tfs_on_first_step
      Rad_control%calc_ch4_tfs_monthly       =   &
                                             calc_ch4_tfs_monthly      
      Rad_control%ch4_tf_time_displacement = ch4_tf_time_displacement
      Rad_control%ch4_tf_calc_intrvl_iz = .true.             
      Rad_control%use_current_ch4_for_tf_iz = .true.
      Rad_control%calc_ch4_tfs_on_first_step_iz = .true.
      Rad_control%calc_ch4_tfs_monthly_iz    =  .true.
      Rad_control%ch4_tf_time_displacement_iz = .true.

!---------------------------------------------------------------------
!    call check_tf_interval to verify that the namelist input variables
!    related to n2o are consistent.
!---------------------------------------------------------------------
      gas_name = 'n2o '
      call check_tf_interval (gas_name, n2o_tf_calc_intrvl, &
                              calc_n2o_tfs_on_first_step,   &
                              calc_n2o_tfs_monthly, &
                              use_current_n2o_for_tf)

!--------------------------------------------------------------------
!    define the radiation_control_type variable components related to
!    n2o tf calculation.
!--------------------------------------------------------------------
      Rad_control%n2o_tf_calc_intrvl = n2o_tf_calc_intrvl
      Rad_control%use_current_n2o_for_tf = use_current_n2o_for_tf
      Rad_control%calc_n2o_tfs_on_first_step =   &
                                           calc_n2o_tfs_on_first_step
      Rad_control%calc_n2o_tfs_monthly       =   &
                                             calc_n2o_tfs_monthly      
      Rad_control%n2o_tf_time_displacement = n2o_tf_time_displacement
      Rad_control%n2o_tf_calc_intrvl_iz = .true.             
      Rad_control%use_current_n2o_for_tf_iz = .true.
      Rad_control%calc_n2o_tfs_on_first_step_iz = .true.
      Rad_control%calc_n2o_tfs_monthly_iz    =  .true.
      Rad_control%n2o_tf_time_displacement_iz = .true.

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
    call longwave_params_init

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(continuum_form) == 'ckd2.1'      .or.   &
          trim(continuum_form) == 'ckd2.4'      .or.   &
          trim(continuum_form) == 'mt_ckd1.0'   .or.   &
          trim(continuum_form) == 'rsb'         .or.   &
          trim(continuum_form) == 'none'        )      then
      else
        write(6,*) "continuum_form is not specified correctly"
        stop
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(linecatalog_form) == 'hitran_1992'  .or.  &
          trim(linecatalog_form) == 'hitran_2000' )  then
      else
        write(6,*) "linecatalog_form is not specified correctly"
        stop
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
    Lw_control%continuum_form     = continuum_form
    Lw_control%linecatalog_form   = linecatalog_form
    Lw_control%do_ch4lbltmpint = do_ch4lbltmpint
    Lw_control%do_n2olbltmpint = do_n2olbltmpint
    Lw_control%do_ch4lbltmpint_iz  = .true.
    Lw_control%do_n2olbltmpint_iz  = .true.

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then
      Lw_parameters%offset = 32
      NBLY = NBLY_CKD
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then
     Lw_parameters%offset = 0
      NBLY = NBLY_RSB  
   endif
   Lw_parameters%offset_iz = .true.

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_1992' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          !  ckd rndm coeff for 40 bands (160-560) and 8 wide bands (560-1400)
          acomb_c=(/ 0.849130E+03,  0.135587E+05,  0.286836E+04,  0.169580E+04,  0.208642E+05, &
                     0.126034E+04,  0.109494E+05,  0.335111E+03,  0.488806E+04,  0.860045E+04, &
                     0.537333E+03,  0.437769E+04,  0.345836E+04,  0.129538E+03,  0.463562E+04, &
                     0.251489E+03,  0.256264E+04,  0.485091E+03,  0.889881E+03,  0.116598E+04, &
                     0.125244E+03,  0.457264E+03,  0.142197E+03,  0.444551E+03,  0.301446E+02, &
                     0.392750E+03,  0.436426E+02,  0.347449E+02,  0.612509E+02,  0.142506E+03, &
                     0.103643E+02,  0.721701E+02,  0.315040E+02,  0.941759E+01,  0.540473E+02, &
                     0.350084E+02,  0.300816E+02,  0.379020E+01,  0.125727E+02,  0.545869E+01, &
                     0.788160E+01,  0.207769E+01,  0.433980E+00,  0.499718E-01,  0.166671E-01, &
                     0.213073E-01,  0.210228E+00,  0.458804E-02 /)
          bcomb_c=(/ 0.175174E+00,  0.176667E+00,  0.109512E+00,  0.111893E+00,  0.145289E+00, &
                     0.203190E+00,  0.151547E+00,  0.911103E-01,  0.151444E+00,  0.850764E-01, &
                     0.756520E-01,  0.100377E+00,  0.171557E+00,  0.125429E+00,  0.105915E+00, &
                     0.816331E-01,  0.149472E+00,  0.857054E-01,  0.107092E+00,  0.185458E+00, &
                     0.753818E-01,  0.108639E+00,  0.123934E+00,  0.178712E+00,  0.833855E-01, &
                     0.119886E+00,  0.133082E+00,  0.935851E-01,  0.156848E+00,  0.166457E+00, &
                     0.162215E+00,  0.114845E+00,  0.724304E-01,  0.740525E-01,  0.734090E-01, &
                     0.141319E+00,  0.359408E-01,  0.833614E-01,  0.128919E+00,  0.996329E-01, &
                     0.676888E-01,  0.690856E-01,  0.660847E-01,  0.730849E-01,  0.654438E-01, &
                     0.848294E-01,  0.852442E-01,  0.242787E+00 /)
          apcm_c=(/ 0.549325E-02, -0.150653E-02,  0.268788E-02,  0.138495E-01, -0.714528E-03, &
                    0.112319E-01,  0.113418E-02,  0.215116E-01,  0.388898E-02,  0.398385E-02, &
                    0.931768E-02,  0.655185E-02,  0.735642E-02,  0.190346E-01,  0.104483E-01, &
                    0.917671E-02,  0.108668E-01,  0.305797E-02,  0.163975E-01,  0.147718E-01, &
                    0.485502E-02,  0.223258E-01,  0.567357E-02,  0.197808E-01,  0.245634E-01, &
                    0.116045E-01,  0.269989E-01,  0.176298E-01,  0.128961E-01,  0.134788E-01, &
                    0.391238E-01,  0.117165E-01,  0.691808E-02,  0.202443E-01,  0.137798E-01, &
                    0.215153E-01,  0.154358E-01,  0.850256E-02,  0.111306E-01,  0.185757E-01, &
                    0.168938E-01,  0.211805E-01,  0.224487E-01,  0.257983E-01,  0.282393E-01, &
                    0.293678E-01,  0.204784E-01,  0.383510E-01 /)
          bpcm_c=(/ -0.305151E-04, -0.244741E-05, -0.203093E-04, -0.736015E-04, -0.158662E-04, &
                    -0.381826E-04, -0.197166E-04, -0.984160E-04, -0.222455E-04, -0.346880E-04, &
                    -0.395593E-04, -0.426165E-04, -0.410312E-04, -0.848479E-04, -0.597304E-04, &
                    -0.318474E-04, -0.450295E-04, -0.284497E-04, -0.772035E-04, -0.545821E-04, &
                    -0.242272E-04, -0.105653E-03, -0.854473E-05, -0.672510E-04, -0.109627E-03, &
                    -0.330446E-04, -0.682675E-04,  0.479154E-04,  0.411211E-04, -0.554504E-04, &
                    -0.145967E-03, -0.425913E-04,  0.413272E-05, -0.531586E-04, -0.429814E-04, &
                    -0.847248E-04, -0.733456E-04,  0.403362E-05, -0.389712E-04, -0.531450E-04, &
                    -0.560515E-04, -0.645812E-04, -0.559626E-04, -0.456728E-04, -0.789113E-04, &
                    -0.105722E-03, -0.811760E-04, -0.103490E-03 /)
          atpcm_c=(/ 0.541966E-02, -0.153876E-02,  0.158818E-02,  0.133698E-01, -0.270746E-02, &
                     0.691660E-02,  0.485749E-05,  0.199036E-01,  0.319826E-02,  0.220802E-02, &
                     0.985921E-02,  0.462151E-02,  0.554947E-02,  0.149315E-01,  0.911982E-02, &
                     0.696417E-02,  0.892579E-02,  0.222579E-02,  0.123105E-01,  0.129295E-01, &
                     0.745423E-02,  0.203878E-01,  0.597427E-02,  0.170838E-01,  0.231443E-01, &
                     0.127692E-01,  0.222678E-01,  0.165331E-01,  0.141333E-01,  0.141307E-01, &
                     0.312569E-01,  0.114137E-01,  0.126050E-01,  0.242966E-01,  0.145196E-01, &
                     0.224802E-01,  0.150399E-01,  0.173815E-01,  0.103438E-01,  0.188690E-01, &
                     0.167743E-01,  0.195884E-01,  0.228968E-01,  0.274106E-01,  0.305775E-01, &
                     0.292968E-01,  0.201658E-01,  0.373711E-01 /)
          btpcm_c=(/ -0.202513E-04,  0.948663E-06, -0.130243E-04, -0.678688E-04, -0.986181E-05, &
                     -0.258818E-04, -0.139292E-04, -0.916890E-04, -0.138148E-04, -0.275165E-04, &
                     -0.298451E-04, -0.310005E-04, -0.314745E-04, -0.695971E-04, -0.486158E-04, &
                     -0.198383E-04, -0.421494E-04, -0.102396E-04, -0.591516E-04, -0.575729E-04, &
                     -0.238464E-04, -0.938688E-04, -0.885666E-05, -0.728295E-04, -0.938897E-04, &
                     -0.448622E-04, -0.642904E-04, -0.102699E-04, -0.348743E-05, -0.568427E-04, &
                     -0.104122E-03, -0.313943E-04,  0.939109E-05, -0.631332E-04, -0.325944E-04, &
                     -0.757699E-04, -0.398066E-04,  0.285026E-04, -0.355222E-04, -0.266443E-04, &
                     -0.510390E-04, -0.643921E-04, -0.722986E-04, -0.742483E-04, -0.903545E-04, &
                     -0.102420E-03, -0.702446E-04, -0.120756E-03 /)
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          !  rsb rndm coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          acomb_n=(/ 0.151896E+05,  0.331743E+04,  0.526549E+03,  0.162813E+03,  0.268531E+03, &
                     0.534040E+02,  0.267777E+02,  0.123003E+02,  0.790655E+01,  0.208292E+01, &
                     0.435231E+00,  0.501124E-01,  0.167058E-01,  0.213731E-01,  0.211154E+00, &
                     0.458804E-02 /)
          bcomb_n=(/ 0.153898E+00,  0.116352E+00,  0.102814E+00,  0.966280E-01,  0.128586E+00, &
                     0.119393E+00,  0.898891E-01,  0.645340E-01,  0.676982E-01,  0.691090E-01, &
                     0.660834E-01,  0.730815E-01,  0.654709E-01,  0.848528E-01,  0.852302E-01, &
                     0.242787E+00 /)
          apcm_n=(/ -0.531605E-03,  0.679588E-02,  0.145133E-01,  0.928230E-02,  0.120727E-01, &
                     0.159210E-01,  0.175268E-01,  0.150281E-01,  0.168815E-01,  0.211690E-01, &
                     0.224373E-01,  0.257869E-01,  0.282310E-01,  0.293583E-01,  0.204657E-01, &
                     0.383510E-01 /)
          bpcm_n=(/ -0.122002E-04, -0.335081E-04, -0.449999E-04, -0.231290E-04, -0.385298E-04, &
                    -0.157170E-04,  0.183434E-04, -0.501718E-04, -0.560024E-04, -0.645285E-04, &
                    -0.559339E-04, -0.456444E-04, -0.788833E-04, -0.105692E-03, -0.811350E-04, &
                    -0.103490E-03 /)
          atpcm_n=(/ -0.156433E-02,  0.611348E-02,  0.134469E-01,  0.871421E-02,  0.133191E-01, &
                      0.164697E-01,  0.199640E-01,  0.163219E-01,  0.167628E-01,  0.195779E-01, &
                      0.228872E-01,  0.274006E-01,  0.305678E-01,  0.292862E-01,  0.201527E-01, &
                      0.373711E-01 /)
          btpcm_n=(/ -0.698856E-05, -0.295269E-04, -0.503674E-04, -0.268392E-04, -0.496663E-04, &
                     -0.333122E-04, -0.337346E-04, -0.258706E-04, -0.510036E-04, -0.643516E-04, &
                     -0.722669E-04, -0.742132E-04, -0.903238E-04, -0.102388E-03, -0.702017E-04, &
                     -0.120756E-03 /)
        endif
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          !  ckd rndm coeff for 40 bands (160-560) and 8 wide bands (560-1400)
          acomb_c=(/ 0.852443E+03,  0.135680E+05,  0.287214E+04,  0.169786E+04,  0.208774E+05, &
                     0.126347E+04,  0.109575E+05,  0.335659E+03,  0.489190E+04,  0.860683E+04, &
                     0.537720E+03,  0.438100E+04,  0.346085E+04,  0.129930E+03,  0.463856E+04, &
                     0.251866E+03,  0.256429E+04,  0.485501E+03,  0.890461E+03,  0.116674E+04, &
                     0.125349E+03,  0.457574E+03,  0.142291E+03,  0.444840E+03,  0.301673E+02, &
                     0.392997E+03,  0.436753E+02,  0.347728E+02,  0.612912E+02,  0.142595E+03, &
                     0.103720E+02,  0.722146E+02,  0.315233E+02,  0.942401E+01,  0.513072E+02, &
                     0.333356E+02,  0.282452E+02,  0.364795E+01,  0.120566E+02,  0.524435E+01, &
                     0.766557E+01,  0.206654E+01,  0.434517E+00,  0.484834E-01,  0.149319E-01, &
                     0.208170E-01,  0.213094E+00,  0.411377E-02 /)
          bcomb_c=(/ 0.198784E+00,  0.180681E+00,  0.118658E+00,  0.119064E+00,  0.146603E+00, &
                     0.220406E+00,  0.155103E+00,  0.976744E-01,  0.156269E+00,  0.877632E-01, &
                     0.779564E-01,  0.103170E+00,  0.175775E+00,  0.138026E+00,  0.107161E+00, &
                     0.886369E-01,  0.151599E+00,  0.898578E-01,  0.108788E+00,  0.188239E+00, &
                     0.785380E-01,  0.111297E+00,  0.126167E+00,  0.181797E+00,  0.864333E-01, &
                     0.121608E+00,  0.136989E+00,  0.976423E-01,  0.159293E+00,  0.167830E+00, &
                     0.165480E+00,  0.115992E+00,  0.731373E-01,  0.759051E-01,  0.850304E-01, &
                     0.157206E+00,  0.468525E-01,  0.928124E-01,  0.154442E+00,  0.113983E+00, &
                     0.801016E-01,  0.925560E-01,  0.871279E-01,  0.971572E-01,  0.948153E-01, &
                     0.978814E-01,  0.863140E-01,  0.246144E+00 /)
          apcm_c=(/ 0.546539E-02, -0.150952E-02,  0.268314E-02,  0.138377E-01, -0.717680E-03, &
                    0.112103E-01,  0.113131E-02,  0.214898E-01,  0.388605E-02,  0.398117E-02, &
                    0.931482E-02,  0.654904E-02,  0.735371E-02,  0.190129E-01,  0.104452E-01, &
                    0.917983E-02,  0.108637E-01,  0.305806E-02,  0.163946E-01,  0.147688E-01, &
                    0.485485E-02,  0.223223E-01,  0.567185E-02,  0.197775E-01,  0.245593E-01, &
                    0.116015E-01,  0.269930E-01,  0.176268E-01,  0.128932E-01,  0.134757E-01, &
                    0.391168E-01,  0.117135E-01,  0.691524E-02,  0.202408E-01,  0.138023E-01, &
                    0.215678E-01,  0.154495E-01,  0.854992E-02,  0.111647E-01,  0.185593E-01, &
                    0.169321E-01,  0.212939E-01,  0.229009E-01,  0.269002E-01,  0.290811E-01, &
                    0.295474E-01,  0.204218E-01,  0.397540E-01 /)
          bpcm_c=(/ -0.305711E-04, -0.267166E-05, -0.205244E-04, -0.737412E-04, -0.160900E-04, &
                    -0.383000E-04, -0.199417E-04, -0.984509E-04, -0.224718E-04, -0.349128E-04, &
                    -0.397845E-04, -0.428414E-04, -0.412565E-04, -0.849355E-04, -0.599544E-04, &
                    -0.320908E-04, -0.452531E-04, -0.286591E-04, -0.774278E-04, -0.548061E-04, &
                    -0.244370E-04, -0.105872E-03, -0.875835E-05, -0.674738E-04, -0.109824E-03, &
                    -0.332675E-04, -0.684616E-04,  0.476755E-04,  0.408936E-04, -0.556738E-04, &
                    -0.146146E-03, -0.428146E-04,  0.391166E-05, -0.533799E-04, -0.429577E-04, &
                    -0.848461E-04, -0.735139E-04,  0.475512E-05, -0.391524E-04, -0.528768E-04, &
                    -0.566290E-04, -0.659191E-04, -0.583359E-04, -0.489007E-04, -0.799401E-04, &
                    -0.108513E-03, -0.817415E-04, -0.124641E-03 /)
          atpcm_c=(/ 0.511647E-02, -0.149589E-02,  0.168155E-02,  0.130626E-01, -0.267983E-02, &
                     0.672956E-02,  0.320920E-04,  0.194706E-01,  0.319870E-02,  0.227324E-02, &
                     0.987638E-02,  0.464699E-02,  0.557515E-02,  0.145732E-01,  0.912502E-02, &
                     0.727025E-02,  0.892842E-02,  0.246321E-02,  0.123671E-01,  0.129566E-01, &
                     0.758092E-02,  0.203546E-01,  0.613865E-02,  0.170817E-01,  0.230765E-01, &
                     0.128080E-01,  0.221615E-01,  0.166057E-01,  0.141910E-01,  0.141719E-01, &
                     0.311608E-01,  0.114747E-01,  0.126978E-01,  0.242290E-01,  0.143966E-01, &
                     0.222507E-01,  0.145647E-01,  0.163848E-01,  0.105427E-01,  0.185355E-01, &
                     0.169812E-01,  0.200214E-01,  0.233600E-01,  0.280627E-01,  0.313602E-01, &
                     0.293442E-01,  0.201411E-01,  0.386191E-01 /)
          btpcm_c=(/ -0.195751E-04,  0.614318E-06, -0.131839E-04, -0.659056E-04, -0.100020E-04, &
                     -0.254160E-04, -0.142042E-04, -0.884293E-04, -0.141881E-04, -0.279176E-04, &
                     -0.300722E-04, -0.313475E-04, -0.318919E-04, -0.679308E-04, -0.488824E-04, &
                     -0.218240E-04, -0.423678E-04, -0.108001E-04, -0.595518E-04, -0.578809E-04, &
                     -0.244881E-04, -0.937841E-04, -0.908428E-05, -0.729140E-04, -0.927663E-04, &
                     -0.450776E-04, -0.636769E-04, -0.114250E-04, -0.409253E-05, -0.570420E-04, &
                     -0.103645E-03, -0.316653E-04,  0.907677E-05, -0.627932E-04, -0.343476E-04, &
                     -0.757366E-04, -0.433972E-04,  0.284288E-04, -0.359012E-04, -0.326903E-04, &
                     -0.517480E-04, -0.645067E-04, -0.709235E-04, -0.729186E-04, -0.904949E-04, &
                     -0.971933E-04, -0.682280E-04, -0.131039E-03 /)
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          !  rsb rndm coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          acomb_n=(/ 0.151998E+05,  0.332049E+04,  0.526971E+03,  0.163056E+03,  0.268698E+03, &
                     0.525512E+02,  0.267972E+02,  0.117653E+02,  0.769374E+01,  0.207347E+01, &
                     0.436111E+00,  0.486557E-01,  0.149839E-01,  0.208859E-01,  0.213994E+00, &
                     0.411377E-02 /)
          bcomb_n=(/ 0.156626E+00,  0.121214E+00,  0.106041E+00,  0.102696E+00,  0.130106E+00, &
                     0.125977E+00,  0.924532E-01,  0.758848E-01,  0.801078E-01,  0.925801E-01, &
                     0.871212E-01,  0.971485E-01,  0.948122E-01,  0.979028E-01,  0.863083E-01, &
                     0.246144E+00 /)
          apcm_n=(/ -0.534632E-03,  0.679130E-02,  0.145095E-01,  0.928332E-02,  0.120697E-01, &
                     0.159191E-01,  0.175239E-01,  0.150708E-01,  0.169191E-01,  0.212791E-01, &
                     0.228860E-01,  0.268843E-01,  0.290691E-01,  0.295375E-01,  0.204095E-01, &
                     0.397540E-01 /)
          bpcm_n=(/ -0.124246E-04, -0.337232E-04, -0.452213E-04, -0.233748E-04, -0.387528E-04, &
                    -0.152874E-04,  0.181160E-04, -0.499820E-04, -0.565675E-04, -0.658376E-04, &
                    -0.582841E-04, -0.488377E-04, -0.799114E-04, -0.108471E-03, -0.816965E-04, &
                    -0.124641E-03 /)
          atpcm_n=(/ -0.152854E-02,  0.607607E-02,  0.134403E-01,  0.882810E-02,  0.133579E-01, &
                      0.164403E-01,  0.199834E-01,  0.159606E-01,  0.169686E-01,  0.200086E-01, &
                      0.233479E-01,  0.280500E-01,  0.313482E-01,  0.293334E-01,  0.201281E-01, &
                      0.386191E-01 /)
          btpcm_n=(/ -0.722716E-05, -0.295919E-04, -0.505120E-04, -0.278532E-04, -0.498659E-04, &
                     -0.340396E-04, -0.341553E-04, -0.281647E-04, -0.517000E-04, -0.644504E-04, &
                     -0.708746E-04, -0.728627E-04, -0.904533E-04, -0.971527E-04, -0.681793E-04, &
                     -0.131039E-03 /)
        endif
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the number of separate bands in the 1200 - 1400 cm-1
!    region. this region may be broken into 1, 2 , 4, 10 or 20 individ-
!    ual bands, or it may remain part of the 0-160, 1200-2200 band. 
!---------------------------------------------------------------------
      nbtrge = no_h2o_bands_1200_1400 
      nbtrg  = no_h2o_bands_1200_1400 
      Lw_parameters%NBTRG  = nbtrg
      Lw_parameters%NBTRGE = nbtrge

!---------------------------------------------------------------------
!    set flag indicating these values have been initialized.
!---------------------------------------------------------------------
      Lw_parameters%NBTRG_iz = .true.
      Lw_parameters%NBTRGE_iz = .true.

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      call lw_gases_stdtf_init(pref)
      call optical_path_init(pref)
      call longwave_tables_init (Lw_tables, tabsr, &
                                 tab1, tab2, tab3, tab1w, &
                                 tab1a, tab2a, tab3a)
      if (Lw_control%do_co2_iz) then
        if (do_nlte .and. .not. Lw_control%do_co2) then
         write(6,*) "cannot activate nlte when co2 not active as radiative gas"
         stop
        endif 
      else  ! (do_co2_iz)
        write(6,*) "do_co2 not yet defined"
        stop
      endif ! (do_co2_iz)

!---------------------------------------------------------------------
!    define pressure-dependent index values used in the infrared
!    radiation code. by this manner, the coefficients are defined
!    at execution time (not dependent on the choice of vertical
!    layers)
!      prnlte : pressure (mb) below which non-LTE code (Nlte.F) affects
!               CO2 transmissivities
!---------------------------------------------------------------------
      if (do_nlte) then
        prnlte = 0.1

!--------------------------------------------------------------------
!    abort execution if trying to run with modified radiation grid
!    code must be added to properly map plm (on the model grid)
!    to the radiation grid. if simply dropping upper levels, then is
!    fairly easy.  abort here so that problem may be addressed and to
!    prevent uncertified results.
!    solution is likely to be to pass in plm on radiation grid, whatever
!    that is.
!--------------------------------------------------------------------
        kmin = 1
        plm (kmin) = 0.
        do k=kmin+1,kmax
          plm (k) = 0.5*(pref (k-1,1) + pref (k,1))
        end do
        plm (kmax+1) = pref (kmax+1,1)

!--------------------------------------------------------------------
!    convert pressure specification for bottom (flux) pressure level
!    for nlte calculation into an index (ixprnlte)
!!! CAN THE MODEL TOP BE AT A PRESSURE THAT IS NOT ZERO ??
!!! if not, then must define plm(ks) always to be 0.0
!!! implications for lw_gas tf calcs ??
!      kmax = size(plm,1) - 1
!-------------------------------------------------------------------
        ixprnlte = 1 
        do k=ks+1, kmax
          if (plm(k)*1.0E-02  .LT. prnlte) then
            ixprnlte = k-ks+1 
          else
            exit
          endif
        end do

!---------------------------------------------------------------------
!   allocate and obtain elements of the source function for bands in 
!   the 15 um range (used in nlte)
!---------------------------------------------------------------------
        ioffset = Lw_parameters%offset
        allocate ( c1b7    (NBCO215) )
        allocate ( c2b7    (NBCO215) )
        do n=1,NBCO215 
          cent(n) = 0.5E+00*(Lw_tables%bdlocm(n+8+ioffset) +   &
                             Lw_tables%bdhicm(n+8+ioffset))
          del (n) = Lw_tables%bdhicm(n+8+ioffset) -    &
                    Lw_tables%bdlocm(n+8+ioffset)
          c1b7(n) = (3.7412E-05)*cent(n)*cent(n)*cent(n)*del(n) 
          c2b7(n) = (1.4387E+00)*cent(n)
        end do
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      call gas_tf_init (pref)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      call longwave_clouds_init

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        band_no_start = band_no_start_ckd
        band_no_end   = band_no_end_ckd
        NBLY = NBLY_CKD
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        band_no_start = band_no_start_rsb  
        band_no_end   = band_no_end_rsb  
        NBLY = NBLY_RSB 
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      Lw_parameters%NBLY = NBLY
      Lw_parameters%NBLY_iz = .true.
      Lw_control%do_lwcldemiss = do_lwcldemiss
      Lw_control%do_lwcldemiss_iz = .true.

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_lwcldemiss) then
        cld_indx_table = cld_indx_table_lwclde
        Cldrad_control%nlwcldb = 7
      else
        cld_indx_table = cld_indx_table_rsb
        Cldrad_control%nlwcldb = 1
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      allocate (apcm(NBLY))
      allocate (bpcm(NBLY))
      allocate (atpcm(NBLY))
      allocate (btpcm(NBLY))
      allocate (acomb(NBLY))
      allocate (bcomb(NBLY))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        apcm = apcm_c
        atpcm = atpcm_c
        bpcm = bpcm_c
        btpcm = btpcm_c
        acomb = acomb_c
        bcomb = bcomb_c
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        apcm = apcm_n
        atpcm = atpcm_n
        bpcm = bpcm_n
        btpcm = btpcm_n
        acomb = acomb_n
        bcomb = bcomb_n
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      allocate (cld_indx(6+NBTRGE) )
      do nn=1,6+NBTRGE
        if (nn > Cldrad_control%nlwcldb) then
          cld_indx(nn) = Cldrad_control%nlwcldb
        else
          cld_indx(nn) = nn
        endif
      end do
      if (NBTRGE == 0 .and. .not. use_bnd1_cldtf_for_h2o_bands) then
        write(6,*) "must use band1 cld tfs for the 1200-1400 cm(-1) bands when"
        write(6,*) "they are included in the 1200-2200 cm(-1) band"
        stop
      endif 

      if (NBTRGE  > 0 .and. use_bnd1_cldtf_for_h2o_bands) then
        cld_indx(7:) = 1
      endif

!--------------------------------------------------------------------
!
!---------------------------------------------------------------------
      call longwave_fluxes_init

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
     module_is_initialized = .true.
!---------------------------------------------------------------------


end subroutine sealw99_init

!####################################################################

subroutine sealw99_time_vary (Rad_time, Rad_gases)

type(time_type), intent(in) :: Rad_time
type(radiative_gases_type),    intent(inout)    ::  Rad_gases   

         logical                    :: calc_co2, calc_n2o, calc_ch4

         character(len=4)           :: gas_name
!---------------------------------------------------------------------------

      call get_control_gas_tf (calc_co2, calc_ch4, calc_n2o)

!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
      call lw_gases_stdtf_time_vary

!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
      if (Rad_gases%time_varying_ch4 .or.    &
          Rad_gases%time_varying_n2o) then
        if (Rad_gases%time_varying_ch4 .and. .not. calc_ch4) then
          !call error_mesg ('sealw99_mod', &
          !' if ch4 amount is to vary in time, ch4 tfs must be '//&
          !                                'recalculated', FATAL)
          stop
        endif
        if (Rad_gases%time_varying_n2o        .and. .not. calc_n2o) then
          !call error_mesg ('sealw99_mod', &
          !' if n2o amount is to vary in time, n2o tfs must be '//&
          !                                   'recalculated', FATAL)
          stop
        endif

      endif

!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
      if (Rad_gases%time_varying_co2) then 
        if (.not. calc_co2) then
          !call error_mesg ('sealw99_mod', &
          !' if co2 amount is to vary in time, co2 tfs must be '//&
          !                                  'recalculated', FATAL)
          stop
        endif
      endif

!----------------------------------------------------------------------
!    if ch4 is activated in this job, varying in time, and 
!    calculation of ch4 tfs are requested, call obtain_gas_tfs to
!    define the tfs.
!--------------------------------------------------------------------
      if (Rad_gases%time_varying_ch4) then
        if (Lw_control%do_ch4) then
          if (do_ch4_tf_calc) then
            gas_name = 'ch4 '
            call obtain_gas_tfs (gas_name, Rad_time,   &
                                 Rad_gases%Ch4_time,  &
                                 ch4_tf_calc_intrvl,&
                                 Rad_gases%ch4_tf_offset,  &
                                 calc_ch4_tfs_on_first_step, &
                                 calc_ch4_tfs_monthly, &
                                 month_of_ch4_tf_calc, &
                                 Rad_gases%ch4_for_next_tf_calc,  &
                                 Rad_gases%ch4_for_last_tf_calc, &
                                 do_ch4_tf_calc, do_ch4_tf_calc_init)
          endif  ! (do_ch4_tf_calc)

        endif ! (do_ch4)

!---------------------------------------------------------------------
!    if ch4 is not time-varying and it is the initial call to sealw99,
!    call ch4_time_vary to calculate the tfs. set flags to indicate
!    the calculation has been done.
!---------------------------------------------------------------------
      else ! (time_varying_ch4)
        if (Lw_control%do_ch4 .and. do_ch4_tf_calc ) then
          call ch4_time_vary (Rad_gases%rrvch4,0)
       !   do_ch4_tf_calc = .false.
       !   do_ch4_tf_calc_init = .false.
       ! else if (.not. Lw_control%do_ch4) then
       !   do_ch4_tf_calc = .false.
       !   do_ch4_tf_calc_init = .false.
        endif
      endif  ! (time_varying_ch4)

!----------------------------------------------------------------------
!    if n2o is activated in this job, varying in time, and 
!    calculation of n2o tfs are requested, call obtain_gas_tfs to
!    define the tfs.
!--------------------------------------------------------------------
      if (Rad_gases%time_varying_n2o) then
        if (Lw_control%do_n2o) then
          if (do_n2o_tf_calc) then
            gas_name = 'n2o '
            call obtain_gas_tfs (gas_name, Rad_time,   &
                                 Rad_gases%N2o_time,  &
                                 n2o_tf_calc_intrvl,&
                                 Rad_gases%n2o_tf_offset,  &
                                 calc_n2o_tfs_on_first_step, &
                                 calc_n2o_tfs_monthly, &
                                 month_of_n2o_tf_calc, &
                                 Rad_gases%n2o_for_next_tf_calc,  &
                                 Rad_gases%n2o_for_last_tf_calc, &
                                 do_n2o_tf_calc, do_n2o_tf_calc_init)
          endif  ! (do_n2o_tf_calc)
        endif ! (do_n2o)

!---------------------------------------------------------------------
!    if n2o is not time-varying and it is the initial call to sealw99,
!    call n2o_time_vary to calculate the tfs. set flags to indicate
!    the calculation has been done.
!---------------------------------------------------------------------
      else
        if (Lw_control%do_n2o .and. do_n2o_tf_calc) then
          call n2o_time_vary (Rad_gases%rrvn2o,0)
       !   do_n2o_tf_calc = .false.
       !   do_n2o_tf_calc_init = .false.
       ! else if (.not. Lw_control%do_n2o) then
       !   do_n2o_tf_calc = .false.
       !   do_n2o_tf_calc_init = .false.
        endif
      endif  ! (time_varying_n2o)


!----------------------------------------------------------------------
!    if co2 is activated in this job, varying in time, and 
!    calculation of co2 tfs are requested, call obtain_gas_tfs to
!    define the tfs.
!--------------------------------------------------------------------
      if (Rad_gases%time_varying_co2) then
        if (Lw_control%do_co2) then
          if (do_co2_tf_calc) then
            gas_name = 'co2 '
            call obtain_gas_tfs (gas_name, Rad_time,  &
                                 Rad_gases%Co2_time,  &
                                 co2_tf_calc_intrvl,&
                                 Rad_gases%co2_tf_offset,  &
                                 calc_co2_tfs_on_first_step, &
                                 calc_co2_tfs_monthly, &
                                 month_of_co2_tf_calc, &
                                 Rad_gases%co2_for_next_tf_calc,  &
                                 Rad_gases%co2_for_last_tf_calc, &
                                 do_co2_tf_calc, do_co2_tf_calc_init)
          endif  ! (do_co2_tf_calc)
        endif ! (do_co2)

!---------------------------------------------------------------------
!    if co2 is not time-varying and it is the initial call to sealw99,
!    call co2_time_vary to calculate the tfs. set flags to indicate
!    the calculation has been done.
!---------------------------------------------------------------------
      else
! interactive co2 mod for radiation calculation
! here it's hardcoded to recompute co2 TF on the 1st of each month
         if (Rad_gases%use_model_supplied_co2) then
            write(6,*) "ERROR: use_model_supplied_co2 not supported"
            stop
            !call get_date (Rad_time, year, month, day, hour, minute,&
            !     second)
            !if (day == 1 .and. hour == 0 .and. minute == 0 .and. &
            !     second == 0) then
            !   call co2_time_vary (Rad_gases%rrvco2)
            !   Rad_gases%co2_for_last_tf_calc = Rad_gases%rrvco2
            !   do_co2_tf_calc_init = .false.
            !else
            !   if (do_co2_tf_calc_init) then
            !      call co2_time_vary (Rad_gases%co2_for_last_tf_calc)
            !      do_co2_tf_calc_init = .false.
            !   endif
            !endif
         else  !(Rad_gases%use_model_supplied_co2)
            if (Lw_control%do_co2 .and. do_co2_tf_calc) then
               call co2_time_vary (Rad_gases%rrvco2,0)
           !    do_co2_tf_calc = .false.
           !    do_co2_tf_calc_init = .false.
           ! else if (.not. Lw_control%do_co2) then
           !    do_co2_tf_calc = .false.
           !    do_co2_tf_calc_init = .false.
            endif
         endif  !(Rad_gases%use_model_supplied_co2)
      endif  ! (time_varying_co2)

      if (.not.Rad_gases%time_varying_ch4) then
        if (Lw_control%do_ch4 .and. do_ch4_tf_calc ) then
          call ch4_time_vary (Rad_gases%rrvch4,1)
          do_ch4_tf_calc = .false.
          do_ch4_tf_calc_init = .false.
        else if (.not. Lw_control%do_ch4) then
          do_ch4_tf_calc = .false.
          do_ch4_tf_calc_init = .false.
        endif
      endif  ! (time_varying_ch4)

      if (.not.Rad_gases%time_varying_n2o) then
        if (Lw_control%do_n2o .and. do_n2o_tf_calc) then
          call n2o_time_vary (Rad_gases%rrvn2o,1)
          do_n2o_tf_calc = .false.
          do_n2o_tf_calc_init = .false.
        else if (.not. Lw_control%do_n2o) then
          do_n2o_tf_calc = .false.
          do_n2o_tf_calc_init = .false.
        endif
      endif  ! (time_varying_n2o)
      
      if (.not.Rad_gases%time_varying_co2) then
         if (.not.Rad_gases%use_model_supplied_co2) then
            if (Lw_control%do_co2 .and. do_co2_tf_calc) then
               call co2_time_vary (Rad_gases%rrvco2,1)
               do_co2_tf_calc = .false.
               do_co2_tf_calc_init = .false.
            else if (.not. Lw_control%do_co2) then
               do_co2_tf_calc = .false.
               do_co2_tf_calc_init = .false.
            endif
         endif  !(Rad_gases%use_model_supplied_co2)
      endif  ! (time_varying_co2)

!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
      if ((Lw_control%do_co2 .and. calc_co2) .or. &
          (Lw_control%do_ch4 .and. calc_ch4) .or. &
          (Lw_control%do_n2o .and. calc_n2o)) then
        call lw_gases_stdtf_dealloc
      endif
 
!------------------------------------------------------------------------


end subroutine sealw99_time_vary


!#####################################################################
 
subroutine sealw99_endts (Rad_gases_tv)

type(radiative_gases_type), intent(in) :: Rad_gases_tv

       if (Rad_gases_tv%time_varying_ch4) then
         if (Lw_control%do_ch4) then
           if (.not. calc_ch4_tfs_on_first_step) then
             do_ch4_tf_calc = .true.
           endif
          endif
        endif

       if (Rad_gases_tv%time_varying_n2o) then
         if (Lw_control%do_n2o) then
          if (.not. calc_n2o_tfs_on_first_step) then
            do_n2o_tf_calc = .true.
         endif
        endif
      endif

     if (Rad_gases_tv%time_varying_co2) then
       if (Lw_control%do_co2) then
         if (.not. calc_co2_tfs_on_first_step) then
           do_co2_tf_calc = .true.
         endif
       endif
      endif

end subroutine sealw99_endts


!#####################################################################
! <SUBROUTINE NAME="sealw99">
!  <OVERVIEW>
!   Subroutine to calculate longwave radiation flux and heating rate.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine calculates longwave radiation flux and heating rate
!   based on the simplified exchange method. It also provides diagnostics
!   for the longwave radiation and cloud.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sealw99 (is, ie, js, je, Atmos_input, Rad_gases, &
!                 Aerosol, Aerosol_props, Cldrad_props, Cld_spec, &
!                Lw_output, Lw_diagnostics)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting subdomain i indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending subdomain i indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting subdomain j indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending subdomain j indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to longwave radiation
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties
!  </INOUT>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_TYPE">
!   Cloud specification type contains cloud microphysical, geometrical,
!   and distribution properties in a model column.
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </INOUT>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  </INOUT>
! </SUBROUTINE>
!
subroutine sealw99 (is, ie, js, je, Rad_time, Atmos_input, Rad_gases, &
                    Aerosol, Aerosol_props, Cldrad_props, Cld_spec,   &
                    Aerosol_diags, Lw_output, Lw_diagnostics,         &
                    including_aerosols, Lw_clouds, Optical, Gas_tf)

!---------------------------------------------------------------------
!    sealw99 is the longwave driver subroutine.
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1991), 9075-9096.
!
!     (2) schwarzkopf, m. d., and s. b. fels, "improvements to the
!         algorithm for computing co2 transmissivities and cooling
!         rates," journal geophysical research, 90 (1985) 10541-10550.
!
!     (3) fels, s.b., "simple strategies for inclusion of voigt
!         effects in infrared cooling calculations," application
!         optics, 18 (1979), 2634-2637.
!
!     (4) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 10/7/93
!
!     certified:  radiation version 1.0
!---------------------------------------------------------------------

integer,                        intent(in)    ::  is, ie, js, je
type(time_type),                intent(in)    ::  Rad_time
type(atmos_input_type),         intent(in)    ::  Atmos_input  
type(radiative_gases_type),     intent(inout)    ::  Rad_gases   
type(aerosol_type),             intent(in)    ::  Aerosol      
type(aerosol_properties_type),  intent(inout) ::  Aerosol_props      
type(aerosol_diagnostics_type), intent(inout) ::  Aerosol_diags      
type(cldrad_properties_type),   intent(in)    ::  Cldrad_props
type(cld_specification_type),   intent(in)    ::  Cld_spec
type(lw_output_type),           intent(inout) ::  Lw_output   
type(lw_diagnostics_type),      intent(inout) ::  Lw_diagnostics
logical,                        intent(in)    :: including_aerosols  

!-----------------------------------------------------------------------
!  intent(in) variables:
!
!    is,js,ie,je 
!    Atmos_input
!    Rad_gases
!    Aerosol
!    Cldrad_props
!    Cld_spec
!
!  intent(inout) variables:
!
!    Aerosol_props
!    Lw_output
!    Lw_diagnostics
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3),    &
                       Cldrad_control%nlwcldb) ::    &
           cldtf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3)) ::    &
           cnttaub1, cnttaub2, cnttaub3, co21c,         &
           co21r,                                       &
           heatem, overod, tmp1,                        &
           to3cnt

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), 2)  ::    &
            emspec

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2))   ::  &
           co21c_KEp1, co21r_KEp1   

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3), 3) ::    &
           contdg

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3) -1) ::    &
           pdflux

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2))   ::  &
           flx1e1cf, gxctscf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3)) ::    &
           e1ctw1, e1ctw2,                              &
           emisdg, flxcf,                               &
           heatemcf, flx, to3dg,                        &
           tcfc8
                       
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3) - 1) ::    &
           cts_sum, cts_sumcf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3), NBTRGE) ::    &
           emisdgf, tch4n2oe

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3),   &
                                               8)::  &
           sorc

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3),   &
                                       6+NBTRGE     ) ::    &
           source_band, dsrcdp_band

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3) ,  &
                                        6+NBTRGE     ) ::  &
           trans_band1, trans_band2

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                                        6+NBTRGE     ) ::  &
           trans_b2d1, trans_b2d2

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), 2, NBTRGE) ::    &
           emspecf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),  &     
                       size(Atmos_input%press,3)-1 )  :: &
            pdfinv, heatra_save
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),  &
                       size(Atmos_input%press,3) )  :: flxnet_save

      type(lw_clouds_type), intent(inout)    :: Lw_clouds ! MJT moved to outside sealw99
      type(optical_path_type), intent(inout) :: Optical   ! MJT moved to outside sealw99 
      type(gas_tf_type), intent(inout)       :: Gas_tf    ! MJT moved to outside sealw99

      integer                       :: ix, jx, kx
      integer                       :: k, kp, m, j
      integer                       :: kk, i, l
      integer                       :: nprofiles, nnn

!---------------------------------------------------------------------
!  local variables:
!
!     cldtf
!     cnttaub1
!     cnttaub2      
!     cnttaub3      
!     co21c       transmission function for the 560-800 cm-1 band 
!                 from levels k through KE+1 to level k. used for flux
!                 at level k arising from exchange with levels k 
!                 through KE+1. includes co2 (from tables), h2o (after
!                 multiplication with over).
!     co21diag    transmission function for co2 only in the 560-800
!                 cm-1 band, from levels k to k, where k ranges from
!                 KS to KE+1. 
!     co21r       transmission function for the 560-800 cm-1 band 
!                 from level k to levels k+1 through KE+1. used for 
!                 flux at levels k+1 through KE+1 arising from
!                 exchange with level k. includes co2 (from tables),
!                 h2o (after multiplication with over).
!     dsorc15
!     dsorc93
!     dsorcb1
!     dsorcb2
!     dsorcb3
!     dt4       
!     emiss
!     heatem
!     overod
!     sorc15      planck function for 560-800 cm-1 bands (sum over
!                 bands 9 and 10).
!     t4
!     tmp1        temporary array, used for computational purposes
!                 in various places. should have no consequences
!                 beyond the immediate module wherein it is defined.
!     tmp2        temporary arrays, similar to tmp1
!     to3cnt      transmission function for the 990-1070 cm-1 band
!                 including o3(9.6 um) + h2o continuum (no lines) 
!                 and possibly cfcs.
!     ch41c
!     n2o1c
!     n2o17c
!     emspec
!     s1a
!     flxge1
!     co21c_KEp1
!     co21r_KEp1
!     contdg
!     cfc_tf
!     pdflux
!     flx1e1cf
!     flxge1cf
!     gxctscf
!     emissb 
!     e1cts1
!     e1cts2
!     e1ctw1
!     e1ctw2
!     soe2
!     soe3
!     soe4
!     soe5
!     emisdg
!     flxcf
!     heatemcf
!     flx
!     to3dg
!     taero8
!     taero8kp
!     totaer_tmp
!     tcfc8
!     cts_sum
!     cts_sumcf
!     emissbf
!     e1cts1f
!     e1cts2f
!     emisdgf
!     emissf
!     tch4n2oe
!     flx1e1fcf
!     flxge1f
!     flxge1fcf
!     sorc        planck function, at model temperatures, for all
!                 bands;  used in cool-to-space calculations.
!     source_band
!     dsrcdp_band
!     trans_band1
!     trans_band2
!     trans_b2d1
!     trans_b2d2
!     emspecf
!     pdfinv
!     Lw_clouds
!     Optical
!     Gas_tf
!     calc_co2
!     calc_n2o
!     calc_ch4
!     ch4_vmr
!     n2o_vmr
!     co2_vmr
!     ix,jx,kx
!     n,k,kp,m,j
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        write(6,*) "ERROR: sealw99 module has not been initialized"
        stop
      endif

!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
      kx = size(Atmos_input%press,3) - 1
      ix = ie-is+1
      jx = je-js+1

!--------------------------------------------------------------------
!    call sealw99_alloc to allocate component arrays of 
!    lw_diagnostics_type variables.
!----------------------------------------------------------------------
      if ( .not. allocated(Lw_diagnostics%flx1e1)) then
        call sealw99_alloc (ix, jx, kx, Lw_diagnostics)
      else
        Lw_diagnostics%flx1e1   = 0.
        Lw_diagnostics%cts_out    = 0.
        Lw_diagnostics%cts_outcf = 0.
        Lw_diagnostics%gxcts    = 0.
        Lw_diagnostics%excts  = 0.
        Lw_diagnostics%exctsn   = 0.
        Lw_diagnostics%fctsg   = 0.
        Lw_diagnostics%fluxn  = 0.
        if (Rad_control%do_totcld_forcing) then
          Lw_diagnostics%fluxncf = 0.
        endif
        Lw_diagnostics%flx1e1f  = 0.
      endif
      
!----------------------------------------------------------------------
!
!--------------------------------------------------------------------
       call optical_path_setup (is, ie, js, je, Atmos_input, Rad_gases, &
                               Aerosol, Aerosol_props, Aerosol_diags, &
                               Optical, including_aerosols) 
        
!--------------------------------------------------------------------
!    call co2coef to compute some co2 temperature and pressure   
!    interpolation quantities and to compute temperature-corrected co2
!    transmission functions (co2spnb and co2nbl). 
!-------------------------------------------------------------------
      call co2coef (Atmos_input, Gas_tf)

!----------------------------------------------------------------------
!    call co2_source_calc to calculate the source function.
!-----------------------------------------------------------------------
      call co2_source_calc (Atmos_input, Rad_gases, sorc, Gas_tf, &
                            source_band, dsrcdp_band)
      
      if (Cldrad_control%do_ica_calcs) then
        nprofiles = Cldrad_control%nlwcldb
        heatra_save = 0.
        flxnet_save = 0.0
      else
        nprofiles = 1
      endif
 
      do nnn = 1, nprofiles

!  reinitialize these values (done in sealw99_alloc)
       Lw_diagnostics%flx1e1    = 0.
       Lw_diagnostics%cts_out   = 0.
       Lw_diagnostics%cts_outcf = 0.
       Lw_diagnostics%gxcts     = 0.
       Lw_diagnostics%excts     = 0.
       Lw_diagnostics%exctsn    = 0.
       Lw_diagnostics%fctsg     = 0.

       Lw_diagnostics%fluxn(:,:,:,:) = 0.0

       if (Rad_control%do_totcld_forcing) then
         Lw_diagnostics%fluxncf(:,:,:,:) = 0.0
       endif

       if (NBTRGE > 0) then
         Lw_diagnostics%flx1e1f  = 0.
       end if

!----------------------------------------------------------------------
!    call cldtau to compute cloud layer transmission functions for 
!    all layers.
!-------------------------------------------------------------------
      call cldtau (nnn, Cldrad_props, Cld_spec, Lw_clouds)


!---------------------------------------------------------------------
!    BEGIN LEVEL KS CALCULATIONS
!---------------------------------------------------------------------

!-----------------------------------------------------------------
!    compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!    mission functions and n2o 560-670 cm-1 transmission functions
!    appropriate for level KS. 
!---------------------------------------------------------------------
      call transcolrow (Gas_tf, KS, KS, KS, KE+1, KS+1, KE+1,   &
                        co21c, co21r, tch4n2oe)

!---------------------------------------------------------------------
!    go into optical_path_mod to obtain the optical path functions 
!    needed for use from level KS.
!    to3cnt contains values in the 990 - 1070 cm-1 range (ozone, water
!    vapor continuum, aerosol, cfc).
!    overod contains values in the 560 - 800 cm-1 range (water vapor 
!    lines, continuum, aerosol, 17 um n2o band, cfc).
!    cnttaub1 is continuum band 4 (water vapor continuum, aerosol, cfc)
!    cnttaub2 is continuum band 5 (water vapor continuum, aerosol, cfc)
!    cnttaub3 is continuum band 7 (water vapor continuum, aerosol, cfc)
!    the 15um band transmission functions between levels KS and KS+1
!    are stored in overod and co2nbl; they will not be overwritten,
!    as long as calculations are made for pressure levels increasing
!    from KS.
!---------------------------------------------------------------------
       call optical_trans_funct_from_KS (Gas_tf, to3cnt, overod,  &
                                        Optical, cnttaub1, cnttaub2, &
                                        cnttaub3, including_aerosols)   

!-----------------------------------------------------------------------
!    compute cloud transmission functions between level KS and all
!    other levels.
!-----------------------------------------------------------------------
      call cloud (KS, Cldrad_props, Cld_spec, Lw_clouds, cldtf)

!-----------------------------------------------------------------------
!    obtain exact cool-to-space for water and co2, and approximate 
!    cool-to-space for co2 and o3.
!----------------------------------------------------------------------
        call cool_to_space_exact (cldtf, Atmos_input, Optical, Gas_tf, &
                                 sorc, to3cnt, Lw_diagnostics, cts_sum, &
                                 cts_sumcf, gxctscf, including_aerosols)

!----------------------------------------------------------------------
!    compute the emissivity fluxes for k=KS.
!! trans_band1:
!    index 1 = e1flx
!    index 2 = co21r*overod
!    index 3 = cnttaub1
!    index 4 = cnttaub2
!    index 5 = to3cnt
!    index 6 = cnttaub3
!    index 7 = e1flxf
!! trans_band2:
!    index 1 = emiss
!    index 2 = co21c*overod
!    index 3 = cnttaub1
!    index 4 = cnttaub2
!    index 5 = to3cnt
!    index 6 = cnttaub3
!    index 7 = emissf
!----------------------------------------------------------------------
       call e1e290 (Atmos_input,  e1ctw1, e1ctw2, trans_band1,   &
                    trans_band2,  Optical, tch4n2oe, &
                    source_band(:,:,:,1), Lw_diagnostics, cldtf, &
                    cld_indx, flx1e1cf,  tcfc8, including_aerosols)  

     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,2) = co21r   (i,j,kk)*  &
                                               overod(i,j,kk)
           end do
        end do
     end do
     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,3) = cnttaub1(i,j,kk)
           end do
        end do
     end do
     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,4) = cnttaub2(i,j,kk)
           end do
        end do
     end do
     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,5) = to3cnt  (i,j,kk)
           end do
        end do
     end do
     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,6) = cnttaub3(i,j,kk)
           end do
        end do
     end do

!! trans_band2:
!   index 1 = emiss
!   index 2 = co21c*overod
!   index 3 = cnttaub1
!   index 4 = cnttaub2
!   index 5 = to3cnt
!   index 6 = cnttaub3
!   index 7 = emissf


     do kk = KS+1,KE+1
        do j = 1,size(trans_band2(:,:,:,:),2)
           do i = 1,size(trans_band2(:,:,:,:),1)
              trans_band2(i,j,kk,2) = co21c(i,j,kk)*   &
                                       overod(i,j,kk)
           end do
        end do
     end do
     do l = 3,6
        do kk = KS+1,KE+1
           do j = 1,size(trans_band2(:,:,:,:),2)
              do i = 1,size(trans_band2(:,:,:,:),1)
                 trans_band2(i,j,kk,l) = trans_band1(i,j,kk,l)
              end do
           end do
        end do
     end do

!----------------------------------------------------------------------
!     the following is a rewrite of the original code largely to
!     eliminate three-dimensional arrays.  the code works on the
!     following principles.  let k be a fixed flux level and kp be
!     a varying flux level, then
! 
!     flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=KS,KE+1.
!
!     if we assume a symmetrical array tau(k,kp)=tau(kp,k), we can
!     break down the calculations for k=KS,KE+1 as follows:
! 
!     flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=k+1,KE+1            (1)
!
!     flux(kp) =   (deltab(k )*tau(kp,k)) for kp=k+1,KE+1.           (2)
!
!     plus deltab(k)*tau(k,k) for all k.
!
!     if we compute a one-dimensional array tauod(kp) for 
!     kp=k+1,KE+1, equations (1) and (2) become:
!
!     tauod(kp) = tau(kp,k)                                          (3)
!
!     flux (k ) = sum(deltab(kp)*tauod(kp)) for kp=k+1,KE+1          (4)
!
!     flux (kp) =    (deltab(k )*tauod(kp)) for kp=k+1,KE+1          (5)
!
!     where tau(k,k) and nearby layer terms are handled separately.
!
!     compute fluxes at level k = KS
!     compute the terms for flux at levels KS+1 to KE+1 from level KS.
!     compute terms for flux at level KS from level KS.
!     compute the terms for flux at level KS due to levels KP from KS+1 
!     to KE+1.
!
!-----------------------------------------------------------------------
    call longwave_fluxes_ks (source_band, trans_band1, dsrcdp_band,  &
                             trans_band2, cldtf, cld_indx , &
                                          Lw_diagnostics)  

     do kk = KS,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,1) = e1ctw2(i,j,kk)
           end do
        end do
     end do

!-----------------------------------------------------------------------
!     compute approximate cool-to-space heating rates for 1 wide band
!     in the 15um  range (560-800 cm-1) (ctsco2) and for 1 band in 
!     the 9.6 um band (ctso3).
!----------------------------------------------------------------------
    call cool_to_space_approx ( Atmos_input%pflux,  source_band, &
                                trans_band1, cldtf, cld_indx ,   &
                                Lw_diagnostics, e1ctw1)

!-----------------------------------------------------------------------
!     perform flux calculations for the flux levels KS+1 to KE-1. calcu-
!     lations for flux levels KE and KE+1 are done separately, as all
!     calculations are special cases or nearby layers.
!----------------------------------------------------------------------

    do k=KS+1,KE-1

!--------------------------------------------------------------------
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level k. 
!--------------------------------------------------------------------
        call transcolrow (Gas_tf, k, k, k, KE+1, k+1, KE+1,  &
                          co21c, co21r, tch4n2oe)

!-------------------------------------------------------------------
!     the 15 um band transmission functions between levels k and k+1
!     are stored in overod and co2nbl; they will not be overwritten,
!     as long as calculations are made for pressure levels increasing
!     from k.
!---------------------------------------------------------------------
       call optical_trans_funct_k_down (Gas_tf, k, to3cnt, overod,  &
                                         Optical, including_aerosols)  

!-----------------------------------------------------------------------
!     compute cloud transmission functions between level k and all
!     other levels greater or equal to k.
!---------------------------------------------------------------------
      call cloud (k,Cldrad_props,Cld_spec,  Lw_clouds, cldtf)

!-----------------------------------------------------------------------
!     compute the exchange terms in the flux equation (except the 
!     nearby layer (k,k) terms, done later).
!! trans_band1:
!   index 1 = emissb
!   index 2 = co21r*overod
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emissbf
!! trans_band2:
!   index 1 = emiss
!   index 2 = co21c*overod
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emissf
!---------------------------------------------------------------------
       call e290 (Atmos_input, k, trans_band2, trans_band1,   &
                  Optical,  tch4n2oe,  tcfc8, including_aerosols)     
       do kp=k,KE
      do j = 1,size(trans_band1(:,:,:,:),2)
         do i = 1,size(trans_band1(:,:,:,:),1)
            trans_band1(i,j,kp+1,3) = cnttaub1(i,j,kp+1  )/cnttaub1(i,j,k  )
         end do
      end do
      do j = 1,size(trans_band1(:,:,:,:),2)
         do i = 1,size(trans_band1(:,:,:,:),1)
            trans_band1(i,j,kp+1,4) = cnttaub2(i,j,kp+1  )/cnttaub2(i,j,k  )
         end do
      end do
      do j = 1,size(trans_band1(:,:,:,:),2)
         do i = 1,size(trans_band1(:,:,:,:),1)
            trans_band1(i,j,kp+1,6) = cnttaub3(i,j,kp+1  )/cnttaub3(i,j,k  )
         end do
      end do
       end do


     do kk = k+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,2) = co21r(i,j,kk)*   &
                                       overod(i,j,kk)
           end do
        end do
     end do
     do kk = k+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,5) = to3cnt(i,j,kk)
           end do
        end do
     end do


     do kk = k+1,KE+1
        do j = 1,size(trans_band2(:,:,:,:),2)
           do i = 1,size(trans_band2(:,:,:,:),1)
              trans_band2(i,j,kk,2) = co21c(i,j,kk)*   &
                                            overod(i,j,kk)
           end do
        end do
     end do
     do l = 3,6
        do kk = k+1,KE+1
           do j = 1,size(trans_band2(:,:,:,:),2)
              do i = 1,size(trans_band2(:,:,:,:),1)
                 trans_band2(i,j,kk,l) = trans_band1(i,j,kk,l)
              end do
           end do
        end do
     end do

!-----------------------------------------------------------------------
!     compute the terms for flux at levels k+1 to KE+1 from level k.
!     compute the terms for flux at level k due to levels 
!     kp from k+1 to KE+1.
!----------------------------------------------------------------------
      call longwave_fluxes_k_down  (k, dsrcdp_band, trans_band1, &
                                    trans_band2, cldtf, cld_indx,  &
                                    Lw_diagnostics          )  
    end do   ! (end of k=KS+1,KE-1 loop)

!-----------------------------------------------------------------------
!     compute remaining flux terms. these include:
!       1) the (k,k) terms, for pressure levels k from KS+1 to KE-1 
!          (the KS,KS term was handled earlier);
!       2) terms for pressure level KE. these include the (KE,KE) term,
!          computed as in (1), and the (KE,KE+1) and (KE+1,KE) terms,
!          computed somewhat differently from the similar terms at
!          higher levels, owing to the proximity to the surface layer
!          KE+1;
!       3) the term for pressure level KE+1 (the (KE+1,KE+1 term).
!
!     compute k=KE case.  since the kp loop is length one, many 
!     simplifications occur.  the co2 quantities and the emissivity
!     quantities) are computed in the nbl section. therefore, we want
!     to compute over, to3cnt, and contod; according to our notation    
!     over(:,:,KE), to3cnt(:,:,KE), and contod(:,:,KE).  the boundary
!     layer and nearby layer corrections to the transmission functions 
!     are obtained above.  the following ratios are used in various nbl
!     nbl calculations.  the remaining calculations are for:
!
!       1) the (k,k) terms, k=KS+1,KE-1;
!       2) the (KE,KE    ) term;
!       3) the (KE,KE+1  ) term;
!       4) the (KE+1,KE  ) term;
!       5) the (KE+1,KE+1) term.
!
!     each is uniquely handled.  different flux terms are computed
!     differently the fourth section obtains water transmission 
!     functions used in q(approximate) calculations and also makes nbl 
!     corrections:
!  
!       1) emiss (:,:) is the transmission function matrix obtained 
!          using E2spec;
! 
!       2) "nearby layer" corrections (emiss(i,i)) are obtained
!          using E3v88;
! 
!       3) special values at the surface (emiss(KE,KE+1),
!          emiss(KE+1,KE), emiss(KE+1,KE+1)) are calculated.
!
!
!     compute temperature and/or scaled amount indices and residuals 
!     for nearby layer and special layer lookup tables.
!
!          calculation for special cases (KE,KE+1) and (KE+1,KE)
!
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level KE. if activated, save ch4n2o tf term.
!-------------------------------------------------------------------
      call transcolrow (Gas_tf, KE, KE, KE, KE+1, KE+1, KE+1,  &
                        co21c, co21r, tch4n2oe)

!----------------------------------------------------------------------
!     get optical path terms for KE
!----------------------------------------------------------------------
      call optical_trans_funct_KE (Gas_tf, to3cnt, Optical, overod, &
                                   including_aerosols)        

   do j = 1,size(trans_b2d1(:,:,:),2)
      do i = 1,size(trans_b2d1(:,:,:),1)
         trans_b2d1(i,j,3  ) = cnttaub1(i,j,KE+1)/cnttaub1(i,j,KE  )
      end do
   end do
   do j = 1,size(trans_b2d1(:,:,:),2)
      do i = 1,size(trans_b2d1(:,:,:),1)
         trans_b2d1(i,j,4  ) = cnttaub2(i,j,KE+1)/cnttaub2(i,j,KE  )
      end do
   end do
   do j = 1,size(trans_b2d1(:,:,:),2)
      do i = 1,size(trans_b2d1(:,:,:),1)
         trans_b2d1(i,j,6  ) = cnttaub3(i,j,KE+1)/cnttaub3(i,j,KE  )
      end do
   end do

!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE and KE and
!     KE+1
!----------------------------------------------------------------------
      call cloud (KE, Cldrad_props, Cld_spec, Lw_clouds, cldtf)
 
!-------------------------------------------------------------------- 
!     compute mean temperature in the "nearby layer" between a flux
!     level and the first data level below the flux level (tpl1) or the
!     first data level above the flux level (tpl2)
!---------------------------------------------------------------------
      call esfc  (Atmos_input, emspec, Optical, emspecf, &
                  tch4n2oe, tcfc8)

!----------------------------------------------------------------------
!     compute nearby layer transmission functions for 15 um band, cont-
!     inuum bands, and 9.3 um band in subroutine Nearbylyrtf. trans-
!     mission functions for the special cases (KE,KE+1) and (KE+1,KE)
!     are also computed for the 15 um band.
!! trans_band1:
!    index 1 = emspec(KS+1)
!    index 2 = co21c_KEp1  
!    index 3 = contodb1
!    index 4 = contodb2
!    index 5 = to3cnt
!    index 6 = contodb3
!    index 7 = emspecf(KS+1)

!! trans_band2:
!    index 1 = emspec(KS)
!    index 2 = co21r_KEp1   
!    index 3 = contodb1
!    index 4 = contodb2
!    index 5 = to3cnt
!    index 6 = contodb3
!    index 7 = emspecf(KS)
!----------------------------------------------------------------------
    call trans_sfc    (Gas_tf, Atmos_input, overod, co21c_KEp1, &
                       co21r_KEp1)

     do j = 1,size(trans_b2d1(:,:,:),2)
        do i = 1,size(trans_b2d1(:,:,:),1)
           trans_b2d1(i,j,1) = emspec(i,j,KS+1)
        end do
     end do
     do j = 1,size(trans_b2d1(:,:,:),2)
        do i = 1,size(trans_b2d1(:,:,:),1)
           trans_b2d1(i,j,2) = co21c_KEp1(i,j)
        end do
     end do
     do j = 1,size(trans_b2d1(:,:,:),2)
        do i = 1,size(trans_b2d1(:,:,:),1)
           trans_b2d1(i,j,5) = to3cnt(i,j,KE+1)
        end do
     end do

     do m=1,NBTRGE
     do j = 1,size(trans_b2d1(:,:,:),2)
        do i = 1,size(trans_b2d1(:,:,:),1)
           trans_b2d1(i,j,6+m) = emspecf(i,j,KS+1,m)
        end do
     end do
     end do

     do j = 1,size(trans_b2d2(:,:,:),2)
        do i = 1,size(trans_b2d2(:,:,:),1)
           trans_b2d2(i,j,1) = emspec(i,j,KS)
        end do
     end do
     do j = 1,size(trans_b2d2(:,:,:),2)
        do i = 1,size(trans_b2d2(:,:,:),1)
           trans_b2d2(i,j,2) = co21r_KEP1(i,j)
        end do
     end do
     do kk = 3,6
        do j = 1,size(trans_b2d2(:,:,:),2)
           do i = 1,size(trans_b2d2(:,:,:),1)
              trans_b2d2(i,j,kk) = trans_b2d1(i,j,kk)
           end do
        end do
     end do
     do m=1,NBTRGE
     do j = 1,size(trans_b2d2(:,:,:),2)
        do i = 1,size(trans_b2d2(:,:,:),1)
           trans_b2d2(i,j,6+m) = emspecf(i,j,KS,m)
        end do
     end do
     end do

!-----------------------------------------------------------------------
!     obtain fluxes for the two terms (KE,KE+1) and (KE+1,KE), both 
!     using the same cloud transmission functions (from layer KE)
!----------------------------------------------------------------------
    call longwave_fluxes_KE_KEp1 (dsrcdp_band, trans_b2d1, &
                                  trans_b2d2, cldtf, cld_indx,  &
                                  Lw_diagnostics )

!---------------------------------------------------------------------
!     call enear to calculate emissivity arrays
!----------------------------------------------------------------------
      call enear (Atmos_input, emisdg, Optical, emisdgf , tch4n2oe, &
                  tcfc8)

!-------------------------------------------------------------------
!     obtain optical path transmission functions for diagonal terms
!------------------------------------------------------------------
    call optical_trans_funct_diag (Atmos_input, contdg, to3dg, &
                                   Optical)
 
!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE+1 and KE+1
!------------------------------------------------------------------
     call cloud (KE+1, Cldrad_props, Cld_spec, Lw_clouds, cldtf)
 
!----------------------------------------------------------------------
!     compute nearby layer transmission functions for 15 um band, cont-
!     inuum bands, and 9.3 um band in subroutine Nearbylyrtf. trans-
!     mission functions for the special cases (KE,KE+1) and (KE+1,KE)
!     are also computed for the 15 um band.
!----------------------------------------------------------------------
     call trans_nearby (Gas_tf, Atmos_input, overod,  co21c)


     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,1) = emisdg(i,j,kk)
           end do
        end do
     end do
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,2) = co21c(i,j,kk)
           end do
        end do
     end do
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,3) = contdg(i,j,kk,1)
           end do
        end do
     end do
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,4) = contdg(i,j,kk,2)
           end do
        end do
     end do
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,5) = to3dg(i,j,kk)
           end do
        end do
     end do
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,6) = contdg(i,j,kk,3)
           end do
        end do
     end do

     do m=1,NBTRGE
     do kk = ks+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,6+m) = emisdgf(i,j,kk,m)
           end do
        end do
     end do
     end do

!-----------------------------------------------------------------------
!     obtain fluxes for the diagonal terms at all levels.
!-----------------------------------------------------------------------
    call longwave_fluxes_diag (dsrcdp_band, trans_band1, cldtf , &
                               cld_indx, Lw_diagnostics )
    
!--------------------------------------------------------------------
!      sum up fluxes over bands
!-----------------------------------------------------------------------
    if (Rad_control%do_totcld_forcing) then
      call longwave_fluxes_sum (is, ie, js, je, flx, NBTRGE, &
                                Lw_diagnostics, flxcf)
    else
      call longwave_fluxes_sum (is, ie, js, je, flx, NBTRGE, &
                                Lw_diagnostics)
    endif
    Lw_output%bdy_flx(:,:,1) = Lw_output%bdy_flx(:,:,1) + &
                               Lw_diagnostics%fluxn(:,:,1,3) + &
                               Lw_diagnostics%fluxn(:,:,1,4) + &
                               Lw_diagnostics%fluxn(:,:,1,5) + &
                               Lw_diagnostics%fluxn(:,:,1,6) 
    Lw_output%bdy_flx(:,:,2) = Lw_output%bdy_flx(:,:,2) + &
                               Lw_diagnostics%fluxn(:,:,1,4) 
    Lw_output%bdy_flx(:,:,3) = Lw_output%bdy_flx(:,:,3) + &
                               Lw_diagnostics%fluxn(:,:,ke+1,3) + &
                               Lw_diagnostics%fluxn(:,:,ke+1,4) + &
                               Lw_diagnostics%fluxn(:,:,ke+1,5) + &
                               Lw_diagnostics%fluxn(:,:,ke+1,6) 
    Lw_output%bdy_flx(:,:,4) = Lw_output%bdy_flx(:,:,4) + &
                               Lw_diagnostics%fluxn(:,:,ke+1,4) 
!   Lw_output%bdy_flx = 1.0E-03*Lw_output%bdy_flx

    if (nnn == 1) then ! need do only once
    if (Rad_control%do_totcld_forcing) then
      Lw_output%bdy_flx_clr(:,:,1) = Lw_diagnostics%fluxncf(:,:,1,3) + &
                                 Lw_diagnostics%fluxncf(:,:,1,4) + &
                                 Lw_diagnostics%fluxncf(:,:,1,5) + &
                                 Lw_diagnostics%fluxncf(:,:,1,6) 
      Lw_output%bdy_flx_clr(:,:,2) = Lw_diagnostics%fluxncf(:,:,1,4) 
      Lw_output%bdy_flx_clr(:,:,3) = Lw_diagnostics%fluxncf(:,:,ke+1,3) + &
                                 Lw_diagnostics%fluxncf(:,:,ke+1,4) + &
                                 Lw_diagnostics%fluxncf(:,:,ke+1,5) + &
                                 Lw_diagnostics%fluxncf(:,:,ke+1,6) 
      Lw_output%bdy_flx_clr(:,:,4) = Lw_diagnostics%fluxncf(:,:,ke+1,4) 
      Lw_output%bdy_flx_clr = 1.0E-03*Lw_output%bdy_flx_clr
    endif
    endif

!-----------------------------------------------------------------------
!     compute emissivity heating rates.
!-----------------------------------------------------------------------

     do kk = ks,ke
        do j = 1,size(pdfinv(:,:,:),2)
           do i = 1,size(pdfinv(:,:,:),1)
              pdfinv(i,j,kk) = 1.0/(Atmos_input%pflux(i,j,kk+ks+1-(ks)) -  &
                            Atmos_input%pflux(i,j,kk))
           end do
        end do
     end do



    do kk = KS,KE
       do j = 1,size(heatem(:,:,:),2)
          do i = 1,size(heatem(:,:,:),1)
             heatem(i,j,kk) = (radcon_mks*(flx(i,j,kk+KS+1-(KS)) -    &
                        flx(i,j,kk))*pdfinv(i,j,kk))*1.0e-03
          end do
       end do
    end do
    if (Rad_control%do_totcld_forcing) then                    
      do kk = KS,KE
         do j = 1,size(heatemcf(:,:,:),2)
            do i = 1,size(heatemcf(:,:,:),1)
               heatemcf(i,j,kk) = (radcon_mks*(flxcf(i,j,kk+KS+1-(KS)) -    &
                            flxcf(i,j,kk))*pdfinv(i,j,kk)*1.0e-03)
            end do
         end do
      end do
    endif

!-----------------------------------------------------------------------
!     compute total heating rates.
!-----------------------------------------------------------------------
!--------------------------------------------------------------------
!     cts_sum is the sum of the values from cool_to_space_exact and
!     the values defined here in cool_to_space_approx. it will be used
!     by longwave_driver_mod.
!--------------------------------------------------------------------
    do kk = 1,size(cts_sum(:,:,:),3)
       do j = 1,size(cts_sum(:,:,:),2)
          do i = 1,size(cts_sum(:,:,:),1)
             cts_sum(i,j,kk) = ((((((cts_sum(i,j,kk) -   &
               Lw_diagnostics%cts_out(i,j,kk,2)) -  & 
               Lw_diagnostics%cts_out(i,j,kk,5) )-  &
               Lw_diagnostics%cts_out(i,j,kk,1))  - &
               Lw_diagnostics%cts_out(i,j,kk,3)) - &
               Lw_diagnostics%cts_out(i,j,kk,4))- &
               Lw_diagnostics%cts_out(i,j,kk,6))
          end do
       end do
    end do

    if (Rad_control%do_totcld_forcing) then
    do kk = 1,size(cts_sumcf(:,:,:),3)
       do j = 1,size(cts_sumcf(:,:,:),2)
          do i = 1,size(cts_sumcf(:,:,:),1)
             cts_sumcf(i,j,kk) = ((((((cts_sumcf(i,j,kk) -  &
                            Lw_diagnostics%cts_outcf(i,j,kk,2)) - &
                            Lw_diagnostics%cts_outcf(i,j,kk,5)) - &
                            Lw_diagnostics%cts_outcf(i,j,kk,1)) - &
                            Lw_diagnostics%cts_outcf(i,j,kk,3)) - &
                            Lw_diagnostics%cts_outcf(i,j,kk,4) ) - &
                            Lw_diagnostics%cts_outcf(i,j,kk,6))
          end do
       end do
    end do
    endif
      do kk = KS,KE
         do j = 1,size(Lw_output%heatra(:,:,:),2)
            do i = 1,size(Lw_output%heatra(:,:,:),1)
               Lw_output%heatra(i,j,kk) = heatem(i,j,kk) +   &
               cts_sum  (i,j,kk)  
            end do
         end do
      end do

    if (nnn == 1) then ! only need to do once
    if (Rad_control%do_totcld_forcing) then                    
      do kk = KS,KE
         do j = 1,size(Lw_output%heatracf(:,:,:),2)
            do i = 1,size(Lw_output%heatracf(:,:,:),1)
               Lw_output%heatracf(i,j,kk) = heatemcf(i,j,kk) +   &
                           cts_sumcf(i,j,kk) 
            end do
         end do
      end do
    endif
    endif !  (nnn == 1)

!-----------------------------------------------------------------------
!     compute the flux at each flux level using the flux at the
!     top (flx1e1 + gxcts) and the integral of the heating rates.
!---------------------------------------------------------------------
      do kk = KS,KE
         do j = 1,size(pdflux(:,:,:),2)
            do i = 1,size(pdflux(:,:,:),1)
               pdflux(i,j,kk) = Atmos_input%pflux(i,j,kk+KS+1-(KS)) -   &
                              Atmos_input%pflux(i,j,kk)
            end do
         end do
      end do

      do j = 1,size(Lw_output%flxnet(:,:,:),2)
         do i = 1,size(Lw_output%flxnet(:,:,:),1)
            Lw_output%flxnet(i,j,KS   ) = Lw_diagnostics%flx1e1(i,j) + Lw_diagnostics%gxcts(i,j)
         end do
      end do


!---------------------------------------------------------------------
! convert values to mks (1.0e-03 factor) 
!---------------------------------------------------------------------
      do j = 1,size(Lw_diagnostics%gxcts(:,:),2)
         do i = 1,size(Lw_diagnostics%gxcts(:,:),1)
            Lw_diagnostics%gxcts(i,j) = 1.0e-03*Lw_diagnostics%gxcts(i,j)
         end do
      end do
      do j = 1,size(Lw_diagnostics%flx1e1(:,:),2)
         do i = 1,size(Lw_diagnostics%flx1e1(:,:),1)
            Lw_diagnostics%flx1e1(i,j) = 1.0e-03*Lw_diagnostics%flx1e1(i,j)
         end do
      end do
      if (NBTRGE> 0) then
      do kk = 1,size(Lw_diagnostics%flx1e1f(:,:,:),3)
         do j = 1,size(Lw_diagnostics%flx1e1f(:,:,:),2)
            do i = 1,size(Lw_diagnostics%flx1e1f(:,:,:),1)
               Lw_diagnostics%flx1e1f(i,j,kk) = &
                    1.0e-03*Lw_diagnostics%flx1e1f(i,j,kk)
            end do
         end do
      end do
      endif


!---------------------------------------------------------------------
!    convert mks values to cgs (1.0e03 factor) so can be summed with
!    cgs value.
!---------------------------------------------------------------------
      do kk = KS,KE
         do j = 1,size(tmp1(:,:,:),2)
            do i = 1,size(tmp1(:,:,:),1)
               tmp1(i,j,kk) = 1.0e03*Lw_output%heatra(i,j,kk)*pdflux(i,j,kk)/radcon_mks
            end do
         end do
      end do
    do k=KS+1,KE+1
      do j = 1,size(Lw_output%flxnet(:,:,:),2)
         do i = 1,size(Lw_output%flxnet(:,:,:),1)
            Lw_output%flxnet(i,j,k) = Lw_output%flxnet(i,j,k-1) + tmp1(i,j,k-1)
         end do
      end do
    enddo

   if (nnn == 1) then   ! only need to do once
   if (Rad_control%do_totcld_forcing) then                    
     do j = 1,size(Lw_output%flxnetcf(:,:,:),2)
        do i = 1,size(Lw_output%flxnetcf(:,:,:),1)
           Lw_output%flxnetcf(i,j,KS   ) = flx1e1cf(i,j) + gxctscf(i,j)
        end do
     end do
!---------------------------------------------------------------------
!    convert mks values to cgs (1.0e03 factor) so can be summed 
!    with cgs value.
!---------------------------------------------------------------------
     do kk = KS,KE
        do j = 1,size(tmp1(:,:,:),2)
           do i = 1,size(tmp1(:,:,:),1)
              tmp1(i,j,kk) = 1.0e03*Lw_output%heatracf(i,j,kk)*pdflux(i,j,kk)/radcon_mks 
           end do
        end do
     end do
     do k=KS+1,KE+1
       do j = 1,size(Lw_output%flxnetcf(:,:,:),2)
          do i = 1,size(Lw_output%flxnetcf(:,:,:),1)
             Lw_output%flxnetcf(i,j,k) = Lw_output%flxnetcf(i,j,k-1) + tmp1(i,j,k-1)
          end do
       end do
     enddo
   endif
   endif  ! (nnn == 1) 
 
   if (Cldrad_control%do_ica_calcs) then
     heatra_save = heatra_save + Lw_output%heatra
     flxnet_save = flxnet_save + Lw_output%flxnet
   endif
 
   end do  ! (profiles loop)
 
    if (Cldrad_control%do_ica_calcs) then
      Lw_output%heatra = heatra_save / Float(nprofiles)
      Lw_output%flxnet = flxnet_save / Float(nprofiles)
      Lw_output%bdy_flx = 1.0e-03*Lw_output%bdy_flx/Float(nprofiles)
    endif

!-----------------------------------------------------------------------
!    call thickcld to perform "pseudo-convective adjustment" for
!    maximally overlapped clouds, if desired.
!-----------------------------------------------------------------------
   if (do_thick) then
       call thickcld (Atmos_input%pflux, Cldrad_props, Cld_spec, &
                      Lw_output)
   endif  ! (do_thick)

!--------------------------------------------------------------------
!   convert lw fluxes to mks units.
!---------------------------------------------------------------------
     do kk = 1,size(Lw_output%flxnet(:,:,:),3)
        do j = 1,size(Lw_output%flxnet(:,:,:),2)
           do i = 1,size(Lw_output%flxnet(:,:,:),1)
              Lw_output%flxnet(i,j,kk) = 1.0E-03*Lw_output%flxnet(i,j,kk)
           end do
        end do
     end do
     if (Rad_control%do_totcld_forcing) then
        do kk = 1,size(Lw_output%flxnetcf(:,:,:),3)
           do j = 1,size(Lw_output%flxnetcf(:,:,:),2)
              do i = 1,size(Lw_output%flxnetcf(:,:,:),1)
                 Lw_output%flxnetcf(i,j,kk) = 1.0E-03*Lw_output%flxnetcf(i,j,kk)
              end do
           end do
        end do
      endif

!--------------------------------------------------------------------
!    call lw_clouds_dealloc to deallocate component arrays of Lw_clouds.
!--------------------------------------------------------------------
      call lw_clouds_dealloc (Lw_clouds)

!--------------------------------------------------------------------
!    call gas_tf_dealloc to deallocate component arrays of Gas_tf.
!--------------------------------------------------------------------
      call gas_tf_dealloc (Gas_tf)

!--------------------------------------------------------------------
!    call optical_dealloc to deallocate component arrays of Optical.
!--------------------------------------------------------------------
      call optical_dealloc (Optical, including_aerosols)      
      
!--------------------------------------------------------------------
      
end subroutine sealw99 

!#####################################################################
! <SUBROUTINE NAME="sealw99_end">
!  <OVERVIEW>
!   sealw99_end is the destructor for sealw99_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   sealw99_end is the destructor for sealw99_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sealw99_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine sealw99_end                  

!---------------------------------------------------------------------
!    sealw99_end is the destructor for sealw99_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg( 'sealw99_mod',  &
        !     'module has not been initialized', FATAL )
        stop
      endif

!---------------------------------------------------------------------
!    call the destructor routines for the modules initialized by 
!    sealw99.
!---------------------------------------------------------------------
      call gas_tf_end
      call optical_path_end
      call lw_gases_stdtf_end
      call longwave_clouds_end
      call longwave_fluxes_end
      call longwave_tables_end
      call longwave_params_end

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------


end subroutine sealw99_end                  




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                 
!#####################################################################
! <SUBROUTINE NAME="check_tf_interval">
!  <OVERVIEW>
!   check_tf_interval verifies that requested tf calculation intervals
!   are compatible with radiation time step
!  </OVERVIEW>
!
!  <DESCRIPTION>
!   check_tf_interval verifies the following relationships:
!     1) that the tf calculation interval is no smaller than the
!        radiation time step;
!     2) that the tf calculation interval is an integral multiple of
!        the radiation time step;
!     3) that the specification for calculating tfs on the first step
!        of the job is done properly.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_tf_interval (gas, gas_tf_calc_intrvl, &
!                           calc_gas_tfs_on_first_step,  &
!                           use_current_gas_for_tf)
!   call sealw99_alloc (ix, jx, kx, Lw_diagnostics)  
!  </TEMPLATE>
!  <IN NAME="gas">    
!   name associated with the gas
!  </IN>
!  <IN NAME="gas_tf_calc_intrvl">    
!   time interval between recalculating transmission fumctions [ hours ]
!  </IN>
!  <IN NAME="calc_gas_tfs_on_first_step">  
!   flag indicating if tfs are to be calculated only on first step 
!   of job
!  </IN>
!  <IN NAME="use_current_gas_for_tf">  
!   flag indicating if gas mixing ratio at current time is to be used
!   for calculation of gas tfs 
!  </IN>
! </SUBROUTINE>
! <PUBLICROUTINE>
!
!NOTE: THIS IS A PRIVATE SUBROUTINE.
!

subroutine check_tf_interval (gas, gas_tf_calc_intrvl, &
                              calc_gas_tfs_on_first_step,  &
                              calc_gas_tfs_monthly,        &
                              use_current_gas_for_tf)

!--------------------------------------------------------------------
character(len=4), intent(in) :: gas
real,             intent(in) :: gas_tf_calc_intrvl
logical,          intent(in) :: calc_gas_tfs_on_first_step,  &
                                calc_gas_tfs_monthly,        &
                                use_current_gas_for_tf

! </PUBLICROUTINE>

!---------------------------------------------------------------------
!    if tfs are not being calculated on the first step, the requested 
!    gas transmission function recalculation interval must be greater
!    than the radiation time step. 
!---------------------------------------------------------------------
      if (.not. calc_gas_tfs_on_first_step .and. &
          .not. calc_gas_tfs_monthly) then
        if (INT(3600.0*gas_tf_calc_intrvl) <   &
            Rad_control%rad_time_step) then
          write(6,*) trim(gas)//' tf calculation interval must be greater'
          write(6,*) 'than or equal to the radiation time step'
          stop
        endif

!---------------------------------------------------------------------
!    be sure that the tf calculation interval is an integral multiple 
!    of the radiation timestep.
!---------------------------------------------------------------------
        if (mod(INT(3600.0*gas_tf_calc_intrvl),   &
                Rad_control%rad_time_step) /= 0) then
          write(6,*) trim(gas)//' transmission function calculation interval'
          write(6,*) 'must be integral multiple of radiation time step'
          stop
        endif
      endif ! (.not. calc_gas_tfs_on_first_step)

!---------------------------------------------------------------------
!    to calculate the tfs using the gas value at the start of the run,
!    one must set use_current_gas_for_tf to .false, and set the 
!    gas_tf_time_displacement to 0.0, rather than setting 
!    use_current_gas_for_tf to .true.
!---------------------------------------------------------------------
      if (calc_gas_tfs_on_first_step) then
        if (use_current_gas_for_tf) then
          write(6,*) 'cannot specify use of current '//trim(gas)//' value'
          write(6,*) 'for tfs when calculating tfs on first step; instead'
          write(6,*) 'set use_current_'//trim(gas)//'_for_tf to false and set'
          write(6,*) trim(gas)//'_tf_time_displacement =    0.0'
          stop
        endif
      endif

      if (calc_gas_tfs_on_first_step .and. &
          calc_gas_tfs_monthly) then
        write(6,*) "cannot request calc of tfs both on first step and monthly"
        stop
      endif

!---------------------------------------------------------------------


end subroutine check_tf_interval



!####################################################################
! <SUBROUTINE NAME="obtain_gas_tfs">
!  <OVERVIEW>
!   obtain_gas_tfs obtains the transmission functions for the requested
!   gas
!  </OVERVIEW>
!
!  <DESCRIPTION>
!   obtain_gas_tfs performs the following functions:
!     a) if time variation of the gas has begun at the current time:
!        1) defines how long the gas has been varying and whether
!           the tfs are due to be recalculated at the current time;
!        2) if the tfs are not to be always recalculated on the first
!           step:
!           a) if this is a recalculation step;
!             1) call the routine to calculate the tfs for the input 
!                gas;
!             2) redefine the value of the gas mixing ratio used fro the
!                last tf calculation to be the one just used;
!             3) set the flag indicating the need to initially calculate
!                the tfs to .false.
!           b) if this is not a recalculation step:
!             1) if this is the initial step, call the routine to calc-
!                ulate the tfs for the input gas;
!             2) set the flag indicating the need to initially calculate
!                the tfs to .false.
!        3) if the tfs are to be always calculated on the first step:
!           a) call the routine to calculate the tfs for the input gas;
!           b) redefine the value of the gas mixing ratio used from the
!              last tf calculation to be the one just used;
!           c) set the flag indicating the need to initially calculate
!              the tfs to .false.
!     b) if time variation of the gas has not begun at the current time:
!         1) if this is the initial call of the job, call the routine
!            to calculate the tfs for the input gas;
!         2) set a flag to indicate that the initial call has been made.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_gas_tfs (gas, Rad_time, Gas_time, gas_tf_calc_intrvl,&
!                       gas_tf_offset, calc_gas_tfs_on_first_step, &
!                       gas_for_next_tf_calc, gas_for_last_tf_calc, &
!                       do_gas_tf_calc, do_gas_tf_calc_init)
!  </TEMPLATE>
!  <IN NAME="gas">    
!   name associated with the gas
!  </IN>
!  <IN NAME="Rad_time">
!   current model time [ time_type ]
!  </IN>
!  <IN NAME="Gas_time">
!   time since time variation of gas began [ time_type ]
!  </IN>
!  <IN NAME="gas_tf_calc_intrvl">    
!   time interval between recalculating transmission fumctions [ hours ]
!  </IN>
!  <IN NAME="gas_tf_offset">
!   time difference between current time and the time for which the 
!   tfs are calculated
!  </IN>
!  <IN NAME="calc_gas_tfs_on_first_step">  
!   flag indicating if tfs are to be calculated only on first step 
!   of job
!  </IN>
!  <INOUT NAME="gas_for_next_tf_calc">  
!   value of gas mixing ratio to be used when tfs are next calculated
!   [ no. / no. ]
!  </INOUT>
!  <INOUT NAME="gas_for_last_tf_calc">  
!   value of gas mixing ratio to be used when tfs were last calculated
!   [ no. / no. ]
!  </INOUT>
!  <INOUT NAME="do_gas_tf_calc">
!   if true, calculate gas tfs when alarm again goes off
!  </INOUT>
!  <INOUT NAME="do_gas_tf_calc_init">
!   this variable is true initially to force calculation of the tfs on 
!   the first call of the job; it is then set to false
!  </INOUT>
!</SUBROUTINE>
!
!NOTE: THIS IS A PRIVATE SUBROUTINE.
!
subroutine obtain_gas_tfs (gas, Rad_time, Gas_time, gas_tf_calc_intrvl,&
                           gas_tf_offset, calc_gas_tfs_on_first_step, &
                           calc_gas_tfs_monthly,month_of_gas_tf_calc, &
                           gas_for_next_tf_calc, gas_for_last_tf_calc, &
                           do_gas_tf_calc, do_gas_tf_calc_init)

!----------------------------------------------------------------------
character(len=4),       intent(in)    :: gas
type(time_type),        intent(in)    :: Rad_time, Gas_time
real,                   intent(in)    :: gas_tf_calc_intrvl,  &
                                         gas_tf_offset
logical,                intent(in)    :: calc_gas_tfs_on_first_step
integer,                intent(inout) :: month_of_gas_tf_calc
logical,                intent(in)    :: calc_gas_tfs_monthly       
real,                   intent(inout) :: gas_for_next_tf_calc, &
                                         gas_for_last_tf_calc
logical,                intent(inout) :: do_gas_tf_calc,&
                                         do_gas_tf_calc_init

!---------------------------------------------------------------------
!  local variables:
!
!     Time_since_gas_start   length of time that gas has been varying
!                            [ time_type ]
!     seconds                seconds component of Time_since_gas_start
!                            [ seconds ]
!     days                   days component of Time_since_gas_start
!                            [ days ]
!     minutes_from_start     minutes since time variation started
!                            [ minutes ]
!     alarm                  if alarm = 0, it is time to calculate the
!                            tfs [ minutes ]
!     year                   year component of current time
!     month                  month component of current time
!     day                    day component of current time
!     hour                   hour component of current time
!     minute                 minute component of current time
!     second                 second component of current time
!     chvers, chversx        characters used to output model variables
!                            through the error_mesg interface
!      
!---------------------------------------------------------------------


!--------------------------------------------------------------------
!    if gas variation is underway, define how long it has been since
!    that started, and whether the current time is an integral multiple
!    of the calculation interval from that gas starting time. put
!    the current time into character variables for output via the
!    error_mesg interface.
!--------------------------------------------------------------------
      write(6,*) "ERROR: unsupported option"
      stop
!      if (Rad_time >= Gas_time) then 
!        Time_since_gas_start = Rad_time - Gas_time
!        call get_time (Time_since_gas_start, seconds, days)
!        call get_date (Rad_time, year, month, day, hour, minute, second)
!        write (chvers, '(i4)') year
!        write (chvers2, '(i4)') month
!        write (chvers3, '(i4)') day  
!        write (chvers4, '(i4)') hour 
!        write (chvers5, '(i4)') minute
!        write (chvers6, '(i4)') second
!        write (chvers10, '( f9.3)') gas_tf_offset
!
!!---------------------------------------------------------------------
!!    if tfs are not automatically calculated on the first step of the
!!    job and if the current time is a desired recalculation time, call 
!!    xxx_time_vary to do the calculation.
!!---------------------------------------------------------------------
!!       if (.not. calc_gas_tfs_on_first_step) then
!        if (.not. calc_gas_tfs_on_first_step .and. &
!            .not. calc_gas_tfs_monthly) then
!          minutes_from_start = INT(days*1440.0 + real(seconds)/60.)
!          if (gas_tf_calc_intrvl /= 0.0) then
!            alarm = MOD (minutes_from_start,   &
!                         INT(gas_tf_calc_intrvl*60.0))
!          endif
!          if (alarm == 0) then
!            if (trim(gas) == 'ch4') then
!              call ch4_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(4pe15.7)')  gas_for_next_tf_calc
!            else if (trim(gas) == 'n2o') then
!              call n2o_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!            else if (trim(gas) == 'co2') then
!              call co2_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(3pe15.7)')  gas_for_next_tf_calc
!            endif
!
!!--------------------------------------------------------------------
!!    redefine the value for the gas mixing ratio used fro the last tf 
!!    calculation.
!!--------------------------------------------------------------------
!            gas_for_last_tf_calc = gas_for_next_tf_calc
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!            if (verbose >= 1) then
!              if (gas_tf_offset /= 0.0) then
!                !call error_mesg ('sealw99_mod',  &
!                ! 'calculating '//trim(gas)//' transmission functions&
!                ! & at time '//chvers//chvers2//chvers3//chvers4//  &
!                ! chvers5// chvers6// ', using '//trim(gas)//' &
!                ! &mixing ratio of:' // chvers7 //', which&
!                ! & is the value '//chvers10// 'hours from current &
!                !                                         & time.', NOTE)
!                stop
!              else
!                !call error_mesg ('sealw99_mod',  &
!                ! 'calculating '//trim(gas)//' transmission functions&
!                ! & at time ' //chvers//chvers2//chvers3//chvers4//&
!                ! chvers5// chvers6// ', using '//trim(gas)//'  &
!                ! &mixing ratio of:' // chvers7 //', which&
!                !             & is the value at the current time', NOTE)
!                stop
!              endif
!            endif ! (verbose)
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the initial tf calculation has been 
!!    completed.
!!---------------------------------------------------------------------
!            do_gas_tf_calc_init = .false.
!
!!----------------------------------------------------------------------
!!    if alarm is not 0 and the tfs have not yet been calculated, call 
!!    xxx_time_vary to do the calculation. set the flag appropriately.
!!---------------------------------------------------------------------
!          else
!            if (do_gas_tf_calc_init) then
!              if (trim(gas) == 'ch4') then
!                call ch4_time_vary (gas_for_last_tf_calc)
!                write (chvers7, '(4pe15.7)') gas_for_last_tf_calc
!              else if (trim(gas) == 'n2o') then
!                call n2o_time_vary (gas_for_last_tf_calc)
!                write (chvers7, '(3pe15.7)') gas_for_last_tf_calc
!              else if (trim(gas) == 'co2') then
!                call co2_time_vary (gas_for_last_tf_calc)
!                write (chvers7, '(3pe15.7)') gas_for_last_tf_calc
!              endif
!              do_gas_tf_calc_init = .false.
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!              if (verbose >= 1) then
!                if (gas_tf_offset /= 0.0) then
!                  !call error_mesg ('sealw99_mod',  &
!                  !     'initial '//gas//' transmission function  &
!                  !      &calculation uses '//trim(gas)//' mixing ratio &
!                  !      &of:'//chvers7//'.', NOTE) 
!                  stop
!                else
!                  !call error_mesg ('sealw99_mod',  &
!                  !      'initial '//trim(gas)//' transmission function&
!                  !       & calculation uses '//trim(gas)//' mixing &
!                  !      &ratio of:' // chvers7 //', which is the value &
!                  !      &at the current time.', NOTE)
!                  stop
!                endif
!              endif
!            endif
!          endif !(alarm == 0)
!
!!---------------------------------------------------------------------
!!    if it is desired that the tfs be calculated only on the first 
!!    step, call the appropriate subroutines to do so. redefine the
!!    gas values used for the last tf calculation  redefine the
!!    gas values used for the last tf calculation. 
!!---------------------------------------------------------------------
!!       else ! (.not. calc_gas_tfs_on_first_step)
!        else  if (calc_gas_tfs_on_first_step) then
!          if (trim(gas) == 'ch4') then
!            call ch4_time_vary (gas_for_next_tf_calc)
!            write (chvers7, '(4pe15.7)') gas_for_next_tf_calc
!          else if (trim(gas) == 'n2o') then
!            call n2o_time_vary (gas_for_next_tf_calc)
!            write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!          else if (trim(gas) == 'co2') then
!            call co2_time_vary (gas_for_next_tf_calc)
!            write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!          endif
!          gas_for_last_tf_calc = gas_for_next_tf_calc
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!          if (verbose >= 1) then
!            if (gas_tf_offset /= 0.0) then
!              !call error_mesg ('sealw99_mod',  &
!              !    'calculating '//trim(gas)//' transmission functions&
!              !    & at time '  //chvers//chvers2//chvers3//chvers4// &
!              !    chvers5//chvers6// ', using '//trim(gas)//' mixing &
!              !    &ratio of:' // chvers7 //', which is the value ' &
!              !    //chvers10// 'hours from current time.', NOTE)
!              stop
!            else
!              !call error_mesg ('sealw99_mod',  &
!              !     'calculating '//trim(gas)//' transmission functions&
!              !     & at time ' //chvers//chvers2//chvers3//chvers4//&
!              !     chvers5//chvers6// ', using '//trim(gas)//' mixing &
!              !     &ratio of:' // chvers7 //', which is the value at &
!              !     &the current time.', NOTE)
!              stop
!            endif
!          endif
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the initial tf calculation has been 
!!    completed.
!!---------------------------------------------------------------------
!          do_gas_tf_calc_init = .false.
!        else if (calc_gas_tfs_monthly) then
!          if (do_gas_tf_calc_init)  then          
!            if (trim(gas) == 'ch4') then
!              call ch4_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(4pe15.7)') gas_for_next_tf_calc
!            else if (trim(gas) == 'n2o') then
!              call n2o_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!            else if (trim(gas) == 'co2') then
!              call co2_time_vary (gas_for_next_tf_calc)
!              write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!            endif
!            gas_for_last_tf_calc = gas_for_next_tf_calc
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!            if (verbose >= 1) then
!              if (gas_tf_offset /= 0.0) then
!                !call error_mesg ('sealw99_mod',  &
!                !  'calculating '//trim(gas)//' transmission functions&
!                !  & at time '  //chvers//chvers2//chvers3//chvers4// &
!                !  chvers5//chvers6// ', using '//trim(gas)//' mixing &
!                !  &ratio of:' // chvers7 //', which is the value ' &
!                !  //chvers10// 'hours from current time  .', NOTE)
!                stop
!              else
!                !call error_mesg ('sealw99_mod',  &
!                !   'calculating '//trim(gas)//' transmission functions&
!                !   & at time ' //chvers//chvers2//chvers3//chvers4//&
!                !   chvers5//chvers6// ', using '//trim(gas)//' mixing &
!                !   &ratio of:' // chvers7 //', which is the value at &
!                !   &the current time.', NOTE)
!                stop
!              endif
!            endif
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the initial tf calculation has been 
!!    completed.
!!---------------------------------------------------------------------
!            do_gas_tf_calc_init = .false.
!            month_of_gas_tf_calc = month
!          else ! (do_gas_tf_calc_init)
!            if (month /= month_of_gas_tf_calc) then
!              if (trim(gas) == 'ch4') then
!                call ch4_time_vary (gas_for_next_tf_calc)
!                write (chvers7, '(4pe15.7)') gas_for_next_tf_calc
!              else if (trim(gas) == 'n2o') then
!                call n2o_time_vary (gas_for_next_tf_calc)
!                write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!              else if (trim(gas) == 'co2') then
!                call co2_time_vary (gas_for_next_tf_calc)
!                write (chvers7, '(3pe15.7)') gas_for_next_tf_calc
!              endif
!              gas_for_last_tf_calc = gas_for_next_tf_calc
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!              if (verbose >= 1) then
!                if (gas_tf_offset /= 0.0) then
!                  !call error_mesg ('sealw99_mod',  &
!                  !'calculating '//trim(gas)//' transmission functions&
!                  !& at time '  //chvers//chvers2//chvers3//chvers4// &
!                  !chvers5//chvers6// ', using '//trim(gas)//' mixing &
!                  !&ratio of:' // chvers7 //', which is the value ' &
!                  !//chvers10// 'hours from current time.', NOTE)
!                  stop
!                else
!                  !call error_mesg ('sealw99_mod',  &
!                  ! 'calculating '//trim(gas)//' transmission functions&
!                  ! & at time ' //chvers//chvers2//chvers3//chvers4//&
!                  ! chvers5//chvers6// ', using '//trim(gas)//' mixing &
!                  ! &ratio of:' // chvers7 //', which is the value at &
!                  ! &the current time.', NOTE)
!                  stop
!                endif
!              endif
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the initial tf calculation has been 
!!    completed.
!!---------------------------------------------------------------------
!              month_of_gas_tf_calc = month
!            endif
!          endif ! (do_gas_tf_calc_init)
!        endif ! (.not. calc_gas_tfs_on_first_step)
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the tf calculation has been 
!!    completed on the current timestep.
!!---------------------------------------------------------------------
!        do_gas_tf_calc = .false.
!
!!---------------------------------------------------------------------
!!    if the time variation of the gas has not yet begun, and it is the 
!!    initial call of the job, call the appropriate subroutines to
!!    define the transmission functions.
!!---------------------------------------------------------------------
!      else ! (Rad_time >= Gas_time)
!        if (do_gas_tf_calc_init) then
!          if (trim(gas) == 'ch4') then
!            call ch4_time_vary (gas_for_last_tf_calc)
!            write (chvers7, '(4pe15.7)') gas_for_last_tf_calc
!          else if (trim(gas) == 'n2o') then
!            call n2o_time_vary (gas_for_last_tf_calc)
!            write (chvers7, '(3pe15.7)') gas_for_last_tf_calc
!          else if (trim(gas) == 'co2') then
!            call co2_time_vary (gas_for_last_tf_calc)
!            write (chvers7, '(3pe15.7)') gas_for_last_tf_calc
!          endif
!
!!---------------------------------------------------------------------
!!    set the flag to indicate that the initial tf calculation has been 
!!    completed.
!!---------------------------------------------------------------------
!          do_gas_tf_calc_init = .false.
!
!!---------------------------------------------------------------------
!!    if a record of the tf calculation path is desired, print out the
!!    relevant data.
!!---------------------------------------------------------------------
!          if (verbose >= 1) then
!            write (chvers10, '( f9.3)') gas_tf_offset    
!            if (gas_tf_offset /= 0.0) then
!              !call error_mesg ('sealw99_mod',  &
!              !   'initial '//trim(gas)//' transmission function  &
!              !   &calculation uses '//trim(gas)//' mixing ratio of:'  &
!              !   // chvers7 //', which is the value '//chvers10//  &
!              !   'hours from current time.', NOTE)
!              stop
!            else
!              !call error_mesg ('sealw99_mod',  &
!              !   'initial '//trim(gas)//' transmission function   &
!              !   &calculation uses '//trim(gas)//' mixing ratio of:'  &
!              !   // chvers7 //', which is the value at the current  &
!              !                                          &time.', NOTE)
!              stop
!            endif
!          endif
!        endif ! (do_gas_tf_calc_init)
!      endif   ! (Rad_time >= Gas_time) 
!
!!--------------------------------------------------------------------


end subroutine obtain_gas_tfs 




                                  
!#####################################################################
! <SUBROUTINE NAME="sealw99_alloc">
!  <OVERVIEW>
!   Subroutine to allocate variables needed for longwave diagnostics
!  </OVERVIEW>
!  <TEMPLATE>
!   call sealw99_alloc (ix, jx, kx, Lw_diagnostics)  
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   Dimension 1 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   Dimension 2 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   Dimension 3 length of radiation arrays to be allocated
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_output_type">
!   lw_diagnostics_type variable containing longwave 
!                   radiation output data
!  </INOUT>
! </SUBROUTINE>
!
subroutine sealw99_alloc (ix, jx, kx, Lw_diagnostics)

!--------------------------------------------------------------------
!    sealw99_alloc allocates and initializes the components of the 
!    lw_diagnostics_type variable Lw_diagnostics which holds diagnostic
!    output generated by sealw99_mod.  
!--------------------------------------------------------------------

integer,                   intent(in)    :: ix, jx, kx
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      ix,jx,kx     (i,j,k) lengths of radiation arrays to be allocated
!
!
!   intent(inout) variables:
!
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer ::  NBTRGE, NBLY

      NBTRGE = Lw_parameters%NBTRGE
      NBLY   = Lw_parameters%NBLY
      
      ! gol124: allocate only if not done previously!
      IF (.NOT. allocated(Lw_diagnostics%flx1e1)) THEN
!---------------------------------------------------------------------
!    allocate (and initialize where necessary) lw_diagnostics_type 
!    component arrays.
!---------------------------------------------------------------------

      allocate ( Lw_diagnostics%flx1e1   (ix, jx                ) )
      allocate ( Lw_diagnostics%fluxn    (ix, jx, kx+1, 6+NBTRGE) )
      allocate (Lw_diagnostics%cts_out   (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%cts_outcf (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%gxcts     (ix, jx                ) )
      allocate (Lw_diagnostics%excts     (ix, jx, kx            ) )
      allocate (Lw_diagnostics%exctsn    (ix, jx, kx,   NBLY    ) )
      allocate (Lw_diagnostics%fctsg     (ix, jx,       NBLY    ) )

      allocate( Lw_diagnostics%flx1e1f  (ix, jx,       NBTRGE  ) )

      ! gol124: end if not allocated
      END IF

      Lw_diagnostics%flx1e1    = 0.
      Lw_diagnostics%fluxn     = 0.
      Lw_diagnostics%cts_out   = 0.
      Lw_diagnostics%cts_outcf = 0.
      Lw_diagnostics%gxcts     = 0.
      Lw_diagnostics%excts     = 0.
      Lw_diagnostics%exctsn    = 0.
      Lw_diagnostics%fctsg     = 0.

      Lw_diagnostics%flx1e1f   = 0.

      if (Rad_control%do_totcld_forcing) then
        if (.NOT. allocated(Lw_diagnostics%fluxncf )) then
          allocate ( Lw_diagnostics%fluxncf (ix, jx, kx+1, 6+NBTRGE) )
        endif
        Lw_diagnostics%fluxncf = 0.
      endif

      


!--------------------------------------------------------------------

end subroutine sealw99_alloc



!####################################################################
! <SUBROUTINE NAME="cool_to_space_approx">
!  <OVERVIEW>
!   Subroutine the calculate the cool to space approximation longwave
!   radiation.
!  </OVERVIEW>
!  <TEMPLATE>
!   call cool_to_space_approx (     pflux_in,        source,  &
!                                 trans,      cld_trans, cld_ind, &
!                                 Lw_diagnostics, &
!                                 trans2      )
!  </TEMPLATE>
!  <IN NAME="pflux_in" TYPE="real">
!   pressure values at flux levels
!  </IN>
!  <IN NAME="source" TYPE="real">
!   band integrated longwave source function of each model layer
!  </IN>
!  <IN NAME="trans" TYPE="real">
!   clear sky longwave transmission
!  </IN>
!  <IN NAME="cld_trans" TYPE="real">
!   cloud transmission
!  </IN>
!  <IN NAME="cld_ind" TYPE="real">
!   cloud type index
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_dignostics_type">
!   longwave diagnostics output
!  </INOUT>
!  <IN NAME="trans2" TYPE="real">
!   optional input alternative transmission profile
!  </IN>
! </SUBROUTINE>
! 
subroutine cool_to_space_approx ( pflux_in, source, trans, cld_trans, &
                                  cld_ind, Lw_diagnostics, trans2      )

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)           :: pflux_in
real, dimension (:,:,:,:), intent(in)           :: source, trans, &
                                                   cld_trans
integer, dimension (:),    intent(in)           :: cld_ind
type(lw_diagnostics_type), intent(inout)        :: Lw_diagnostics
real, dimension (:,:,:),   intent(in), optional :: trans2

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     pflux_in
!     source
!     trans
!     cld_trans
!     cld_ind
!     
!  intent(inout) variables:
!
!     Lw_diagnostics
!
!  intent(in),optional:
!     trans2
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

    real, dimension(size(pflux_in,1), size(pflux_in,2), &
                    size(pflux_in,3)-1) :: pdfinv
    integer  ::  i,j,kk
    integer  :: index, nbands

!---------------------------------------------------------------------
!   local variables:
!
!      pdfinv      
!      index
!      nbands
!      j
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
     nbands = size(source,4) - NBTRGE
     do kk = KS,KE
        do j = 1,size(pdfinv(:,:,:),2)
           do i = 1,size(pdfinv(:,:,:),1)
              pdfinv(i,j,kk) = 1.0/(pflux_in(i,j,kk+KS+1-(KS)) - pflux_in(i,j,kk))
           end do
        end do
     end do
!---------------------------------------------------------------------
!
!--------------------------------------------------------------------
     do index=1,nbands
    if (index == 1     ) then
      do kk = KS,KE
         do j = 1,size(Lw_diagnostics%cts_out(:,:,:,:),2)
            do i = 1,size(Lw_diagnostics%cts_out(:,:,:,:),1)
               Lw_diagnostics%cts_out(i,j,kk,index) = (radcon_mks*pdfinv(i,j,kk)*     &
                             source(i,j,kk,index)*  &
                             (trans(i,j,kk,index)*   &
                                cld_trans(i,j,kk+KS+1-(KS), cld_ind(index)) -   &
                             trans2(i,j,kk)*     &
                             cld_trans(i,j,kk, cld_ind(index)))*1.0e-03)
            end do
         end do
      end do
    else
      do kk = KS,KE
         do j = 1,size(Lw_diagnostics%cts_out(:,:,:,:),2)
            do i = 1,size(Lw_diagnostics%cts_out(:,:,:,:),1)
               Lw_diagnostics%cts_out(i,j,kk,index) = (radcon_mks*pdfinv(i,j,kk)*     &
                             source(i,j,kk, index)*  &
                             (trans(i,j,kk+KS+1-(KS),index    )*   &
                             cld_trans(i,j,kk+KS+1-(KS), cld_ind(index)) -   &
                             trans(i,j,kk,index    )*     &
                             cld_trans(i,j,kk, cld_ind(index)))*1.0e-03)
            end do
         end do
      end do
    endif
    if (Rad_control%do_totcld_forcing) then
    if (index == 1     ) then
        do kk = KS,KE
           do j = 1,size(Lw_diagnostics%cts_outcf(:,:,:,:),2)
              do i = 1,size(Lw_diagnostics%cts_outcf(:,:,:,:),1)
                 Lw_diagnostics%cts_outcf(i,j,kk,index) = (radcon_mks*pdfinv(i,j,kk)*     &
                                 source(i,j,kk,index)* &
                                 (trans(i,j,kk,index      ) -     &
                                  trans2(i,j,kk))*1.0e-03)
              end do
           end do
        end do
      else
        do kk = KS,KE
           do j = 1,size(Lw_diagnostics%cts_outcf(:,:,:,:),2)
              do i = 1,size(Lw_diagnostics%cts_outcf(:,:,:,:),1)
                 Lw_diagnostics%cts_outcf(i,j,kk,index) = (radcon_mks*pdfinv(i,j,kk)*     &
                                 source(i,j,kk,index)* &
                                 (trans(i,j,kk+KS+1-(KS), index    ) -     &
                                  trans(i,j,kk, index    ))*1.0e-03)
              end do
           end do
        end do
      endif
    endif
   end do  ! (index loop)

!--------------------------------------------------------------------


end subroutine cool_to_space_approx



!####################################################################
! <SUBROUTINE NAME="cool_space_exact">
!  <OVERVIEW>
!   cool_to_space calculates the cool-to-space cooling rate for 
!   a band n.
!  </OVERVIEW>
!  <TEMPLATE>
!   call  cool_to_space_exact (                cldtf,          &
!                             Atmos_input, Optical, Gas_tf,  &
!                             sorc,        to3cnt, Lw_diagnostics, &
!                             cts_sum, cts_sumcf, &
!                             gxctscf) 
!  </TEMPLATE>
!  <IN NAME="cldtf" TYPE="real">
!   cloud transmission function between levels k level KS.
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input to the cool to space approximation method
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth of atmospheric layers and clouds
!  </INOUT>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission function
!  </INOUT>
!  <IN NAME="sorc" TYPE="real">
!   band-integrated Planck function, for each combined
!   band in the 160-1200 cm-1 region.
!  </IN>
!  <IN NAME="to3cnt" TYPE="real">
!   transmission functions between levels k and
!   level KS for the 990-1070 cm-1 range.
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   Longwave diagnostics 
!  </INOUT>
!  <INOUT NAME="cts_sum" TYPE="real">
!   Cool to space heating rates
!  </INOUT>
!  <INOUT NAME="cts_sumcf" TYPE="real">
!   Cool to space heating rates due to cloud forcing
!  </INOUT>
!  <INOUT NAME="gxctscf" TYPE="real">
!   gxcts is the "exact" surface flux accumulated over
!   the frequency bands in the 160-1200 cm-1 range.
!  </INOUT>
! </SUBROUTINE>
subroutine cool_to_space_exact (cldtf, Atmos_input, Optical, Gas_tf,  &
                                sorc, to3cnt, Lw_diagnostics, &
                                cts_sum, cts_sumcf, gxctscf,  &
                                including_aerosols)     

!-----------------------------------------------------------------------
!    cool_to_space calculates the cool-to-space cooling rate for 
!    a band n.
!-----------------------------------------------------------------------

real, dimension (:,:,:,:), intent(in)     :: cldtf
type(atmos_input_type),    intent(in)     :: Atmos_input
type(optical_path_type),   intent(inout)  :: Optical
type(gas_tf_type),         intent(inout)  :: Gas_tf 
real, dimension (:,:,:,:), intent(in)     :: sorc
real, dimension (:,:,:),   intent(in)     :: to3cnt
type(lw_diagnostics_type), intent(inout)  :: Lw_diagnostics
real, dimension(:,:,:),    intent(inout)  :: cts_sum, cts_sumcf
real, dimension(:,:),      intent(inout)  :: gxctscf
logical,                   intent(in)  :: including_aerosols   

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     cldtf        cloud transmission function between levels k and 
!                  level KS.
!     Atmos_input
!     sorc          band-integrated Planck function, for each combined
!                   band in the 160-1200 cm-1 region.
!     to3cnt        transmission functions between levels k and
!                   level KS for the 990-1070 cm-1 range.
!
!  intent(inout) variables:
!
!     Optical
!     Gas_tf
!     Lw_diagnostics
!     cts_sum
!     cts_sumcf
!     gxctscf
!
!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
      integer        :: i,kk
    real, dimension(size(Atmos_input%pflux,1), &
                    size(Atmos_input%pflux,2), &
                    size(Atmos_input%pflux,3)-1) :: pdfinv, pdfinv2

    real, dimension(size(Atmos_input%pflux,1), &
                    size(Atmos_input%pflux,2), &
                    size(Atmos_input%pflux,3)  ) :: &
                                              dte1, press, temp, pflux

    integer, dimension(size(Atmos_input%pflux,1), &
                       size(Atmos_input%pflux,2), &
                       size(Atmos_input%pflux,3)  ) :: ixoe1 

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2))   ::   &
                                                      pfac1, pfac2

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)) :: &
                              sorc_tmp, ctmp, totch2o_tmp, totaer_tmp

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    2:size(Atmos_input%pflux,3)) :: &
                                                    totvo2_tmp

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)-1) :: &
                                  exctscf, tt, x, y, topm, &
                                  topphi, phitmp, psitmp, ag, &
                                  agg, f, ff, tmp1, tmp2, fac1, &
                                  fac2, cfc_tf

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)-1, NBLY) :: &
                            exctsncf

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                                                 NBLY) :: &
                                fctsgcf
      integer        :: n, k, j, ioffset

!-----------------------------------------------------------------------
!  local variables
!
!     pdfinv      inverse of pressure difference between flux levels.
!     pdfinv2
!     dte1
!     press       pressure at data levels of model.
!     temp        temperature at data levels of model.
!     pflux       pressure at flux levels of model.
!     ixoe1
!     pfac1
!     pfac2
!     sorc_tmp
!     ctmp
!     totch2o_tmp
!     totaer_tmp
!     totvo2_tmp
!     exctscf
!     tt
!     x
!     y
!     topm
!     topphi
!     phitmp
!     psitmp
!     ag
!     agg
!     f 
!     ff
!     tmp1
!     tmp2
!     fac1
!     fac2
!     cfc_tf
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  convert press and pflux to cgs.
        do kk = 1,size(press(:,:,:),3)
           do j = 1,size(press(:,:,:),2)
              do i = 1,size(press(:,:,:),1)
                 press(i,j,kk) = 10.0*Atmos_input%press(i,j,kk)
              end do
           end do
        end do
       do kk = 1,size(pflux(:,:,:),3)
          do j = 1,size(pflux(:,:,:),2)
             do i = 1,size(pflux(:,:,:),1)
                pflux(i,j,kk) = 10.0*Atmos_input%pflux(i,j,kk)
             end do
          end do
       end do
      do kk = 1,size(temp(:,:,:),3)
         do j = 1,size(temp(:,:,:),2)
            do i = 1,size(temp(:,:,:),1)
               temp(i,j,kk) = Atmos_input%temp(i,j,kk)
            end do
         end do
      end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------


      do kk = KS,KE
         do j = 1,size(pdfinv2(:,:,:),2)
            do i = 1,size(pdfinv2(:,:,:),1)
               pdfinv2(i,j,kk) = 1.0/(pflux(i,j,kk+KS+1-(KS)) - pflux(i,j,kk))
            end do
         end do
      end do
     do kk = KS,KE
        do j = 1,size(pdfinv(:,:,:),2)
           do i = 1,size(pdfinv(:,:,:),1)
              pdfinv(i,j,kk) = 1.0/(Atmos_input%pflux(i,j,kk+KS+1-(KS)) -   &
                               Atmos_input%pflux(i,j,kk))
           end do
        end do
     end do
!----------------------------------------------------------------------
      ioffset = Lw_parameters%offset

!-----------------------------------------------------------------------
!     initialize quantities.
!-----------------------------------------------------------------------
      do kk = KS,KE
         do j = 1,size(Lw_diagnostics%excts(:,:,:),2)
            do i = 1,size(Lw_diagnostics%excts(:,:,:),1)
               Lw_diagnostics%excts(i,j,kk) = 0.0E+00
            end do
         end do
      end do
      do j = 1,size(Lw_diagnostics%gxcts(:,:),2)
         do i = 1,size(Lw_diagnostics%gxcts(:,:),1)
            Lw_diagnostics%gxcts(i,j) = 0.0E+00
         end do
      end do

      if (Rad_control%do_totcld_forcing) then
        do kk = KS,KE
           do j = 1,size(exctscf(:,:,:),2)
              do i = 1,size(exctscf(:,:,:),1)
                 exctscf(i,j,kk) = 0.0E+00
              end do
           end do
        end do
        do j = 1,size(gxctscf(:,:),2)
           do i = 1,size(gxctscf(:,:),1)
              gxctscf(i,j) = 0.0E+00
           end do
        end do
      endif

!-----------------------------------------------------------------------
!     compute temperature quantities.
!-----------------------------------------------------------------------
      do kk = KS,KE
         do j = 1,size(x(:,:,:),2)
            do i = 1,size(x(:,:,:),1)
               x(i,j,kk) = temp(i,j,kk) - 2.5E+02
            end do
         end do
      end do
      do kk = KS,KE
         do j = 1,size(y(:,:,:),2)
            do i = 1,size(y(:,:,:),1)
               y(i,j,kk) = x(i,j,kk)*x(i,j,kk)
            end do
         end do
      end do
      do j = 1,size(ctmp(:,:,:),2)
         do i = 1,size(ctmp(:,:,:),1)
            ctmp(i,j,KS) = 1.0E+00
         end do
      end do


!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
       call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      do j = 1,size(Lw_diagnostics%fctsg(:,:,:),2)
         do i = 1,size(Lw_diagnostics%fctsg(:,:,:),1)
            Lw_diagnostics%fctsg(i,j,NBLY) = 0.0
         end do
      end do
      do n=1,NBLY-1
!-----------------------------------------------------------------------
!     obtain temperature correction capphi, cappsi, then multiply
!     by optical path var1, var2 to compute temperature-corrected
!     optical path and mean pressure for a layer: phitmp, psitmp.
!-----------------------------------------------------------------------
        do kk = KS,KE
           do j = 1,size(f(:,:,:),2)
              do i = 1,size(f(:,:,:),1)
                 f(i,j,kk) = 0.44194E-01*(apcm (n)*x(i,j,kk) +  &
                            bpcm (n)*y(i,j,kk)) 
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(ff(:,:,:),2)
              do i = 1,size(ff(:,:,:),1)
                 ff(i,j,kk) = 0.44194E-01*(atpcm(n)*x(i,j,kk) +   &
                            btpcm(n)*y(i,j,kk))
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(ag(:,:,:),2)
              do i = 1,size(ag(:,:,:),1)
                 ag(i,j,kk) = (1.418191E+00 + f (i,j,kk))*   &
                            f (i,j,kk) + 1.0E+00
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(agg(:,:,:),2)
              do i = 1,size(agg(:,:,:),1)
                 agg(i,j,kk) = (1.418191E+00 + ff(i,j,kk))*   &
                            ff(i,j,kk) + 1.0E+00 
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(phitmp(:,:,:),2)
              do i = 1,size(phitmp(:,:,:),1)
                 phitmp(i,j,kk) = Optical%var1(i,j,kk)*     &
                            ((((ag (i,j,kk)*        &
                            ag (i,j,kk))**2)**2)**2)
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(psitmp(:,:,:),2)
              do i = 1,size(psitmp(:,:,:),1)
                 psitmp(i,j,kk) = Optical%var2(i,j,kk)*     &
                            ((((agg(i,j,kk)*       &
                            agg(i,j,kk))**2)**2)**2)
              end do
           end do
        end do

!-----------------------------------------------------------------------
!     obtain optical path and mean pressure from the top of the 
!     atmosphere to the level k.
!-----------------------------------------------------------------------
        do j = 1,size(topm(:,:,:),2)
           do i = 1,size(topm(:,:,:),1)
              topm(i,j,KS) = phitmp(i,j,KS) 
           end do
        end do
        do j = 1,size(topphi(:,:,:),2)
           do i = 1,size(topphi(:,:,:),1)
              topphi(i,j,KS) = psitmp(i,j,KS) 
           end do
        end do
        do k=KS+1,KE
          do j = 1,size(topm(:,:,:),2)
             do i = 1,size(topm(:,:,:),1)
                topm(i,j,k) = topm  (i,j,k-1) + phitmp(i,j,k) 
             end do
          end do
          do j = 1,size(topphi(:,:,:),2)
             do i = 1,size(topphi(:,:,:),1)
                topphi(i,j,k) = topphi(i,j,k-1) + psitmp(i,j,k) 
             end do
          end do
        enddo

!-----------------------------------------------------------------------
!     tt is the cloud-free h2o cool-to-space transmission function.
!-----------------------------------------------------------------------

      if (Lw_control%do_h2o) then
        do kk = KS,KE
           do j = 1,size(fac1(:,:,:),2)
              do i = 1,size(fac1(:,:,:),1)
                 fac1(i,j,kk) = acomb(n)*topm(i,j,kk)
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(fac2(:,:,:),2)
              do i = 1,size(fac2(:,:,:),1)
                 fac2(i,j,kk) = fac1(i,j,kk)*topm(i,j,kk)/   &
                              (bcomb(n)*topphi(i,j,kk))
              end do
           end do
        end do
        do kk = KS,KE
           do j = 1,size(tmp1(:,:,:),2)
              do i = 1,size(tmp1(:,:,:),1)
                 tmp1(i,j,kk) = fac1(i,j,kk)/SQRT(1.0E+00 +     &
                          fac2(i,j,kk))
              end do
           end do
        end do

      else
        do kk = KS,KE
          do j = 1,size(tmp1(:,:,:),2)
            do i = 1,size(tmp1(:,:,:),1)
              tmp1(i,j,kk) = 0.0
            end do
          end do
        end do
      endif
        


!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        if (n >= band_no_start(1) .and. n <= band_no_end(1)) then
!                       160-400 cm-1 region (h2o)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! bands 1-24.
            call get_totch2o (n, Optical, totch2o_tmp, dte1, ixoe1)
            do kk = KS,KE
               do j = 1,size(tt(:,:,:),2)
                  do i = 1,size(tt(:,:,:),1)
                     tt(i,j,kk) = EXP(-1.0*(tmp1(i,j,kk) + diffac*   &
                            totch2o_tmp(i,j,kk+KS+1-(KS))))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! bands 1-4.
            do kk = KS,KE
               do j = 1,size(tt(:,:,:),2)
                  do i = 1,size(tt(:,:,:),1)
                     tt(i,j,kk) = EXP(-1.0*tmp1(i,j,kk)) 
                  end do
               end do
            end do
          endif

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(2) .and. n <= band_no_end(2)) then
!                       400-560 cm-1 region (h2o)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! bands 25-40.
            call get_totch2o (n, Optical, totch2o_tmp, dte1, ixoe1)
            do kk = KS,KE
               do j = 1,size(tt(:,:,:),2)
                  do i = 1,size(tt(:,:,:),1)
                     tt(i,j,kk) = EXP(-1.0*(tmp1(i,j,kk) + diffac*   &
                            totch2o_tmp(i,j,kk+KS+1-(KS))))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! bands 5-8.
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tt(:,:,:),2)
                  do i = 1,size(tt(:,:,:),1)
                     tt(i,j,kk) = EXP(-1.0*(tmp1(i,j,kk) +  &
                                totvo2_tmp(i,j,kk+KS+1-(KS)))) 
                  end do
               end do
            end do
          endif

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(3) .and. n <= band_no_end(3)) then
!                       560-630 cm-1 region (h2o, co2, n2o, aerosol)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! band 41.
            call get_totch2obd (n-40, Optical, totch2o_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) + diffac*     &
                            totch2o_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! band 9
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) +               &
                                totvo2_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          if (including_aerosols) then                    
            do kk = 1,size(totaer_tmp(:,:,:),3)
               do j = 1,size(totaer_tmp(:,:,:),2)
                  do i = 1,size(totaer_tmp(:,:,:),1)
                     totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,1)
                  end do
               end do
            end do
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp2(i,j,kk) +        &
                            totaer_tmp  (i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp2(i,j,kk))*    &
                                (Gas_tf%co2spnb(i,j,kk+KS+1-(KS),1)* &
                                Gas_tf%tn2o17(i,j,kk+KS+1-(KS)))
                end do
             end do
          end do

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(4) .and. n <= band_no_end(4)) then
!                       630-700 cm-1 region (h2o, co2, aerosol)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! band 42.
            call get_totch2obd (n-40, Optical, totch2o_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) + diffac*     &
                            totch2o_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! band 10
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) +               &
                                totvo2_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          if (including_aerosols) then      
            do kk = 1,size(totaer_tmp(:,:,:),3)
               do j = 1,size(totaer_tmp(:,:,:),2)
                  do i = 1,size(totaer_tmp(:,:,:),1)
                     totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,2)
                  end do
               end do
            end do
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp2(i,j,kk) +          &
                            totaer_tmp   (i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp2(i,j,kk))*     &
                                Gas_tf%co2spnb(i,j,kk+KS+1-(KS),2)
                end do
             end do
          end do

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(5) .and. n <= band_no_end(5)) then
!                       700-800 cm-1 region (h2o, co2, aerosol)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! band 43.
            call get_totch2obd (n-40, Optical, totch2o_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) + diffac*      &
                            totch2o_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! band 11
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) +                &
                                totvo2_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          if (including_aerosols) then           
            do kk = 1,size(totaer_tmp(:,:,:),3)
               do j = 1,size(totaer_tmp(:,:,:),2)
                  do i = 1,size(totaer_tmp(:,:,:),1)
                     totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,3)
                  end do
               end do
            end do
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp2(i,j,kk) +       &
                            totaer_tmp   (i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp2(i,j,kk))*     &
                                Gas_tf%co2spnb(i,j,kk+KS+1-(KS),3)
                end do
             end do
          end do

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(6) .and. n <= band_no_end(6)) then
!                       800-990 cm-1 region (h2o, aerosol)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! bands 44-45.
            call get_totch2obd (n-40, Optical, totch2o_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) + diffac*    &
                            totch2o_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! bands 12-13
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) +               &
                                totvo2_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          if (including_aerosols) then    
            do kk = 1,size(totaer_tmp(:,:,:),3)
               do j = 1,size(totaer_tmp(:,:,:),2)
                  do i = 1,size(totaer_tmp(:,:,:),1)
                     totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,n-8-ioffset)
                  end do
               end do
            end do
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp2(i,j,kk) +         &
                            totaer_tmp   (i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp2(i,j,kk))
                end do
             end do
          end do

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(7) .and. n <= band_no_end(7)) then
!                       990-1070 cm-1 region (h2o(lines), o3)
!                       band 46 (ckd2.1) or 14 (rsb)

          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp1(i,j,kk))*      &
                               to3cnt (i,j,kk+KS+1-(KS))
                end do
             end do
          end do

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
else if (n >= band_no_start(8) .and. n <= band_no_end(8)) then
!                       1070-1200 cm-1 region (h2o, n2o)
  if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
      trim(Lw_control%continuum_form) == 'ckd2.4' ) then 
! band 47.
            call get_totch2obd (n-40, Optical, totch2o_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) + diffac*    &
                            totch2o_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
  else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
! band 15
            call get_totvo2 (n, Optical, totvo2_tmp)
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp1(i,j,kk) +              &
                                totvo2_tmp(i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          if (including_aerosols) then      
            do kk = 1,size(totaer_tmp(:,:,:),3)
               do j = 1,size(totaer_tmp(:,:,:),2)
                  do i = 1,size(totaer_tmp(:,:,:),1)
                     totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,7)
                  end do
               end do
            end do
            do kk = KS,KE
               do j = 1,size(tmp2(:,:,:),2)
                  do i = 1,size(tmp2(:,:,:),1)
                     tmp2(i,j,kk) = tmp2(i,j,kk)    +                &
                            totaer_tmp   (i,j,kk+KS+1-(KS))
                  end do
               end do
            end do
          endif
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = EXP(-1.0E+00*tmp2(i,j,kk))*   &
                                Gas_tf%n2o9c  (i,j,kk+KS+1-(KS))
                end do
             end do
          end do

        endif
!--------------------------------------------------------------------
!     calculate or retrieve the source function for the current band.
!--------------------------------------------------------------------
        if (n <= 8 + ioffset) then
          call looktab (tabsr, ixoe1, dte1, sorc_tmp, KS, KE+1, n)
        else
         do kk = 1,size(sorc_tmp(:,:,:),3)
            do j = 1,size(sorc_tmp(:,:,:),2)
               do i = 1,size(sorc_tmp(:,:,:),1)
                  sorc_tmp(i,j,kk) = sorc(i,j,kk,n-8-ioffset)
               end do
            end do
         end do
        endif

!---------------------------------------------------------------------
!     retrieve the cfc effect if cfcs are activated.
!---------------------------------------------------------------------
        if (Lw_control%do_cfc .and. n >= 9+ioffset) then
          call cfc_exact(n-8-ioffset, Optical, cfc_tf)
          do kk = KS,KE
             do j = 1,size(tt(:,:,:),2)
                do i = 1,size(tt(:,:,:),1)
                   tt(i,j,kk) = tt(i,j,kk)*cfc_tf(i,j,kk)
                end do
             end do
          end do
        endif

!------------------------------------------------------------------
!     define some near-surface pressure functions that are needed
!------------------------------------------------------------------
        do j = 1,size(pfac1(:,:),2)
           do i = 1,size(pfac1(:,:),1)
              pfac1(i,j) = 0.5*pdfinv2(i,j,KE)*(pflux(i,j,KE+1) - &
                     press(i,j,KE))*tt(i,j,KE-1)
           end do
        end do
        do j = 1,size(pfac2(:,:),2)
           do i = 1,size(pfac2(:,:),1)
              pfac2(i,j) = 0.5*pdfinv2(i,j,KE)*(pflux(i,j,KE+1) +    &
                     press(i,j,KE) - 2.0*pflux(i,j,KE))*tt(i,j,KE)
           end do
        end do

!--------------------------------------------------------------------
!     calculate the ground fluxes (?)
!--------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          do j = 1,size(fctsgcf(:,:,:),2)
             do i = 1,size(fctsgcf(:,:,:),1)
                fctsgcf(i,j,n) = tt(i,j,KE)*sorc_tmp(i,j,KE) +   &
                           (pfac1(i,j) + pfac2(i,j))*   &
                           (sorc_tmp(i,j,KE+1) - sorc_tmp(i,j,KE))
             end do
          end do
    do j = 1,size(Lw_diagnostics%fctsg(:,:,:),2)
       do i = 1,size(Lw_diagnostics%fctsg(:,:,:),1)
          Lw_diagnostics%fctsg(i,j,n) = cldtf(i,j,KE+1,cld_indx_table(n+32-ioffset))* &
                         fctsgcf(i,j,n)
       end do
    end do
          do j = 1,size(gxctscf(:,:),2)
             do i = 1,size(gxctscf(:,:),1)
                gxctscf(i,j) = gxctscf(i,j) + fctsgcf(i,j,n)
             end do
          end do
          do j = 1,size(Lw_diagnostics%gxcts(:,:),2)
             do i = 1,size(Lw_diagnostics%gxcts(:,:),1)
                Lw_diagnostics%gxcts(i,j) = Lw_diagnostics%gxcts(i,j) +   &
                      Lw_diagnostics%fctsg(i,j,n)
             end do
          end do
        else
do j = 1,size(Lw_diagnostics%fctsg(:,:,:),2)
   do i = 1,size(Lw_diagnostics%fctsg(:,:,:),1)
      Lw_diagnostics%fctsg(i,j,n) = cldtf(i,j,KE+1,cld_indx_table(n+32-ioffset))* &
                         (tt(i,j,KE)*sorc_tmp(i,j,KE) +  &
                         (pfac1(i,j) + pfac2(i,j))*   &
                         (sorc_tmp(i,j,KE+1) - sorc_tmp(i,j,KE)))
   end do
end do
          do j = 1,size(Lw_diagnostics%gxcts(:,:),2)
             do i = 1,size(Lw_diagnostics%gxcts(:,:),1)
                Lw_diagnostics%gxcts(i,j) = Lw_diagnostics%gxcts(i,j) +   &
          Lw_diagnostics%fctsg(i,j,n)
             end do
          end do
        endif
        do j = 1,size(Lw_diagnostics%fctsg(:,:,:),2)
           do i = 1,size(Lw_diagnostics%fctsg(:,:,:),1)
              Lw_diagnostics%fctsg(i,j,n) = 1.0e-03*   &
                                      Lw_diagnostics%fctsg(i,j,n)
           end do
        end do

!--------------------------------------------------------------------
!    include the effect of the cloud transmission function.
!--------------------------------------------------------------------
        do kk = KS+1,KE+1
           do j = 1,size(ctmp(:,:,:),2)
              do i = 1,size(ctmp(:,:,:),1)
                 ctmp(i,j,kk) = tt(i,j,kk+KS-(KS+1))*        &
                      cldtf(i,j,kk,cld_indx_table(n+32-ioffset)) 
              end do
           end do
        end do

!---------------------------------------------------------------------
!    if diagnostics is on, save each band's contribution separately.
!    exctsn is the cool-to-space heating rate for each frequency
!    band. fctsg is the "exact" surface flux for each frequency
!    band in the 160-1200 cm-1 range.
!---------------------------------------------------------------------
!! the following array only needed when diagnostics on
          do kk = KS,KE
             do j = 1,size(Lw_diagnostics%exctsn(:,:,:,:),2)
                do i = 1,size(Lw_diagnostics%exctsn(:,:,:,:),1)
                   Lw_diagnostics%exctsn(i,j,kk,n) = sorc_tmp(i,j,kk)*     &
                                (ctmp(i,j,kk+KS+1-(KS)) - ctmp(i,j,kk))
                end do
             end do
          end do
          if (Rad_control%do_totcld_forcing) then
            do j = 1,size(exctsncf(:,:,:,:),2)
               do i = 1,size(exctsncf(:,:,:,:),1)
                  exctsncf(i,j,KS,n) = sorc_tmp(i,j,KS)*         &
                                  (tt(i,j,KS) - 1.0E+00)
               end do
            end do
            do kk = KS+1,KE
               do j = 1,size(exctsncf(:,:,:,:),2)
                  do i = 1,size(exctsncf(:,:,:,:),1)
                     exctsncf(i,j,kk,n) = sorc_tmp(i,j,kk)*     &
                                     (tt(i,j,kk) - tt(i,j,kk+KS-(KS+1)))
                  end do
               end do
            end do
          endif

!-----------------------------------------------------------------------
!     excts is the cool-to-space cooling rate accumulated over
!     frequency bands.
!-----------------------------------------------------------------------
        do kk = KS,KE
           do j = 1,size(Lw_diagnostics%excts(:,:,:),2)
              do i = 1,size(Lw_diagnostics%excts(:,:,:),1)
                 Lw_diagnostics%excts(i,j,kk) = &
               Lw_diagnostics%excts(i,j,kk) +   &
               Lw_diagnostics%exctsn(i,j,kk,n)  
              end do
           end do
        end do
        if (Rad_control%do_totcld_forcing) then
          do j = 1,size(exctscf(:,:,:),2)
             do i = 1,size(exctscf(:,:,:),1)
                exctscf(i,j,KS) = exctscf(i,j,KS) +     &
                             exctsncf(i,j,KS,n)
             end do
          end do
          do kk = KS+1,KE
             do j = 1,size(exctscf(:,:,:),2)
                do i = 1,size(exctscf(:,:,:),1)
                   exctscf(i,j,kk) = exctscf(i,j,kk) +     &
                             exctsncf(i,j,kk,n)
                end do
             end do
          end do
        endif
      end do

!-----------------------------------------------------------------------
!     gxcts is the "exact" surface flux accumulated over
!     the frequency bands in the 160-1200 cm-1 range.
!     obtain cool-to-space flux at the top by integration of heating
!     rates and using cool-to-space flux at the bottom (current value 
!     of gxcts).  note that the pressure quantities and conversion
!     factors have not been included either in excts or in gxcts.
!     these cancel out, thus reducing computations.
!-----------------------------------------------------------------------
      do k=KS,KE
        do j = 1,size(Lw_diagnostics%gxcts(:,:),2)
           do i = 1,size(Lw_diagnostics%gxcts(:,:),1)
              Lw_diagnostics%gxcts(i,j) = Lw_diagnostics%gxcts(i,j) - Lw_diagnostics%excts(i,j,k)
           end do
        end do
      enddo
      if (Rad_control%do_totcld_forcing) then
        do k=KS,KE
          do j = 1,size(gxctscf(:,:),2)
             do i = 1,size(gxctscf(:,:),1)
                gxctscf(i,j) = gxctscf(i,j) - exctscf(i,j,k)
             end do
          end do
        enddo  
      endif

!-----------------------------------------------------------------------
!     now scale the cooling rate excts by including the pressure 
!     factor pdfinv and the conversion factor radcon.
!-----------------------------------------------------------------------
!! the following array only needed when diagnostics on
        do n=1,NBLY-1
          do kk = KS,KE
             do j = 1,size(Lw_diagnostics%exctsn(:,:,:,:),2)
                do i = 1,size(Lw_diagnostics%exctsn(:,:,:,:),1)
                   Lw_diagnostics%exctsn(i,j,kk,n) = 1.0e-03*(Lw_diagnostics%exctsn(i,j,kk,n)*radcon_mks*     &
                                pdfinv(i,j,kk))
                end do
             end do
          end do
        enddo
        if (Rad_control%do_totcld_forcing) then
          do n=1,NBLY-1
            do kk = KS,KE
               do j = 1,size(exctsncf(:,:,:,:),2)
                  do i = 1,size(exctsncf(:,:,:,:),1)
                     exctsncf(i,j,kk,n) = 1.0e-03*(exctsncf(i,j,kk,n)*radcon_mks*    &
                                    pdfinv(i,j,kk) )
                  end do
               end do
            end do
          enddo
        endif

      do kk = KS,KE
         do j = 1,size(Lw_diagnostics%excts(:,:,:),2)
            do i = 1,size(Lw_diagnostics%excts(:,:,:),1)
               Lw_diagnostics%excts(i,j,kk) = (Lw_diagnostics%excts(i,j,kk)*radcon_mks*pdfinv(i,j,kk))*1.0e-03
            end do
         end do
      end do

      if (Rad_control%do_totcld_forcing) then
        do kk = KS,KE
           do j = 1,size(exctscf(:,:,:),2)
              do i = 1,size(exctscf(:,:,:),1)
                 exctscf(i,j,kk) = (exctscf(i,j,kk)*radcon_mks*pdfinv(i,j,kk))*1.0e-03
              end do
           end do
        end do
      endif

!--------------------------------------------------------------------
!    save the heating rates to be later sent to longwave_driver_mod.
!--------------------------------------------------------------------
      do kk = 1,size(cts_sum(:,:,:),3)
         do j = 1,size(cts_sum(:,:,:),2)
            do i = 1,size(cts_sum(:,:,:),1)
               cts_sum(i,j,kk) = Lw_diagnostics%excts(i,j,kk)
            end do
         end do
      end do
      do kk = 1,size(cts_sumcf(:,:,:),3)
         do j = 1,size(cts_sumcf(:,:,:),2)
            do i = 1,size(cts_sumcf(:,:,:),1)
               cts_sumcf(i,j,kk) = exctscf(i,j,kk)
            end do
         end do
      end do


!--------------------------------------------------------------------


end  subroutine cool_to_space_exact




!####################################################################
! <SUBROUTINE NAME="e1e290">
!  <OVERVIEW>
!   Subroutine to compute thermal exchange terms and emissivities used
!   to obtain the cool-to-space heating rates for all pressure layers.
!  </OVERVIEW>
!  <DESCRIPTION>
!   !     E1e290 computes two different quantities.
!     
!     1) emissivities used to compute the exchange terms for flux at the
!     top of the atmosphere (level KS). (the top layer, isothermal by
!     assumption, does not contribute to photon exchanges with other
!     layers). these terms are obtained using precomputed e2 functions
!     (see ref. (2)).
!
!     2) emissivities used to obtain the cool-to-space heating rates
!     for all pressure layers. these are obtained using precomputed
!     e1 functions (see ref. (2)).
!
!     the frequency ranges for the e2 calculations are 0-560 and 1200-
!     2200 cm-1. the CTS calculations also require calculations in the
!     160-560 cm-1 range. (see refs. (1) and (2)).
!ifdef ch4n2o
!
!     if ch4 and n2o are included, the frequency range for emissivities
!     is 1400-2200 cm-1, with separate emissivity quantities for the
!     1200-1400 cm-1 range.
!endif ch4n2o
!
!     the reason for combining these calculations is that both use
!     the same scaled h2o amount (avephi) as input, thus reducing
!     some calculation time for obtaining index quantities.
!   
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!   </DESCRIPTION>
!   <TEMPLATE>
!    call e1e290 (Atmos_input,                 e1ctw1, e1ctw2,   &
!                 trans_band1, trans_band2, Optical, tch4n2oe, &
!                 t4, Lw_diagnostics, cldtf, cld_indx, flx1e1cf, &
!                 tcfc8)
!   </TEMPLATE>
!   <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmospheric input data to the thermal exchange method
!   </IN>
!   <OUT NAME="e1ctw1" TYPE="real">
!    Cool to space thermal exchange terms in 0-560 band
!   </OUT>
!   <OUT NAME="e1ctw2" TYPE="real">
!    Cool to space thermal exchange terms in 1200-2200  band
!   </OUT>
!   <OUT NAME="trans_band1" TYPE="real">
!    transmission functions in band 1
!   </OUT>
!   <OUT NAME="trans_band2" TYPE="real">
!    transmission function in band 2
!   </OUT>
!   <INOUT NAME="Optical" TYPE="optical_path_type">
!    thermal layer optical path 
!   </INOUT>
!   <IN NAME="tch4n2oe" TYPE="real">
!    CH4 and N2O transmission functions
!   </IN>
!   <IN NAME="t4" TYPE="real">
!    source function of the top most layer
!   </IN>
!   <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!    longwave diagnostics variable
!   </IN>
!   <IN NAME="cld_tf" TYPE="real">
!    Cloud transmission function
!   </IN>
!   <IN NAME="cld_indx" TYPE="real">
!    Cloud type index
!   </IN>
!   <OUT NAME="flx1e1cf" TYPE="real">
!    TOA flux due to cloud forcing
!   </OUT>
!   <INOUT NAME="tcfc8" TYPE="real">
!    CFC transmission function (chloroflurocarbons)
!   </INOUT>
! </SUBROUTINE>
!
subroutine e1e290 (Atmos_input, e1ctw1, e1ctw2, trans_band1, &
                   trans_band2, Optical, tch4n2oe, t4, Lw_diagnostics, &
                   cldtf, cld_indx, flx1e1cf, tcfc8, including_aerosols) 

!-----------------------------------------------------------------------
!
!     E1e290 computes two different quantities.
!     
!     1) emissivities used to compute the exchange terms for flux at the
!     top of the atmosphere (level KS). (the top layer, isothermal by
!     assumption, does not contribute to photon exchanges with other
!     layers). these terms are obtained using precomputed e2 functions
!     (see ref. (2)).
!
!     2) emissivities used to obtain the cool-to-space heating rates
!     for all pressure layers. these are obtained using precomputed
!     e1 functions (see ref. (2)).
!
!     the frequency ranges for the e2 calculations are 0-560 and 1200-
!     2200 cm-1. the CTS calculations also require calculations in the
!     160-560 cm-1 range. (see refs. (1) and (2)).
!ifdef ch4n2o
!
!     if ch4 and n2o are included, the frequency range for emissivities
!     is 1400-2200 cm-1, with separate emissivity quantities for the
!     1200-1400 cm-1 range.
!endif ch4n2o
!
!     the reason for combining these calculations is that both use
!     the same scaled h2o amount (avephi) as input, thus reducing
!     some calculation time for obtaining index quantities.
!   
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------

type(atmos_input_type),    intent(in)    :: Atmos_input
real, dimension (:,:,:),   intent(out)   :: e1ctw1, e1ctw2
real, dimension (:,:,:,:), intent(out)   :: trans_band1, trans_band2  
type(optical_path_type),   intent(inout) :: Optical
real, dimension (:,:,:,:), intent(in)    ::  tch4n2oe                  
real, dimension(:,:,:),    intent(in)    :: t4
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
real, dimension(:,:,:,:),  intent(in)    :: cldtf
integer, dimension(:),     intent(in)    :: cld_indx
real, dimension(:,:),      intent(out)   :: flx1e1cf
real, dimension (:,:,:),   intent(inout) ::  tcfc8           
logical,                   intent(in)            :: including_aerosols  
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input
!      tch4n2oe
!      t4
!      cldtf
!      cld_indx
!
!   intent(inout) variables:
!
!      Optical
!      Lw_diagnostics
!      tcfc8
!
!   intent(out) variables:
!
!      e1ctw1
!      e1ctw2
!      trans_band1
!      trans_band2
!      flx1e1cf
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      real, dimension (size(trans_band2,1), &
                       size(trans_band2,2), &
                       size(trans_band2,3)) :: dte1, dte2,ttmp,ttmp0

      integer, dimension (size(trans_band2,1), &
                          size(trans_band2,2), &
                          size(trans_band2,3)) :: ixoe1, ixoe2

      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::   &
                                     temp, tflux, totaer_tmp, taero8, &
                                     tmp1, tmp2, e1cts1, e1cts2, &
                                     avphilog, dt1, du, dup, du1, dup1

      integer, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::   &
                                         ixo1, iyo, iyop, iyo1, iyop1 

      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2)) :: &
                        s1a, flxge1, flxge1cf

      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), NBTRGE) :: &
                        flx1e1fcf, flxge1f, flxge1fcf
      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3), NBTRGE) ::   &
                         e1cts1f, e1cts2f

      integer  :: i,j,k,kk,m
!---------------------------------------------------------------------
!  local variables
!
!     dte1
!     dte2
!     ixoe1
!     ixoe2  
!     temp
!     tflux
!     totaer_tmp
!     taero8
!     tmp1
!     tmp2
!     e1cts1
!     e1cts2
!     avphilog
!     dt1
!     du
!     dup
!     ixo1
!     iyo
!     iyop
!     s1a
!     flxge1
!     flxge1cf
!     flx1e1fcf
!     flxge1f
!     flxge1fcf
!     e1cts1f
!     e1cts2f
!     k,m
!
!----------------------------------------------------------------------

     do kk = 1,size(tflux(:,:,:),3)
        do j = 1,size(tflux(:,:,:),2)
           do i = 1,size(tflux(:,:,:),1)
              tflux(i,j,kk) = Atmos_input%tflux(i,j,kk)
           end do
        end do
     end do
     do kk = 1,size(temp(:,:,:),3)
        do j = 1,size(temp(:,:,:),2)
           do i = 1,size(temp(:,:,:),1)
              temp(i,j,kk) = Atmos_input%temp(i,j,kk)
           end do
        end do
     end do

!---------------------------------------------------------------------
!     obtain the "exchange" emissivities as a function of temperature 
!     (fxo) and water amount (fyo). the temperature indices have
!     been obtained in longwave_setup_mod.
!-------------------------------------------------------------------
  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  do kk = KS,KE
     do j = 1,size(ixoe2(:,:,:),2)
        do i = 1,size(ixoe2(:,:,:),1)
           ixoe2(i,j,kk) = ixoe2(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  do kk = KS,KE
     do j = 1,size(dte2(:,:,:),2)
        do i = 1,size(dte2(:,:,:),1)
           dte2(i,j,kk) = dte2 (i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  do j = 1,size(ixoe2(:,:,:),2)
     do i = 1,size(ixoe2(:,:,:),1)
        ixoe2(i,j,KE+1) = ixoe1(i,j,KE)
     end do
  end do
  do j = 1,size(dte2(:,:,:),2)
     do i = 1,size(dte2(:,:,:),1)
        dte2(i,j,KE+1) = dte1 (i,j,KE)
     end do
  end do



   if (Lw_control%do_h2o) then
      do kk = KS,KE+1
         do j = 1,size(avphilog(:,:,:),2)
            do i = 1,size(avphilog(:,:,:),1)
               avphilog(i,j,kk) = LOG10(Optical%avephi(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, avphilog, du, iyo, KS, KE+1)
      iyo(:,:,KS:KE+1) = iyo(:,:,KS:KE+1) + 1
      call looktab (tab2, ixoe2, iyo, dte2, du, &
                    trans_band2(:,:,:,1), KS, KE+1)
   else
     iyo1(:,:,KS:KE+1) = 1
     du1(:,:,KS:KE+1) = 0.0
     call looktab (tab2, ixoe2, iyo1, dte2, du1, &
                   trans_band2(:,:,:,1), KS, KE+1)
   endif

!-----------------------------------------------------------------------
!     the special case emiss(:,:,KE+1) for layer KE+1 is obtained by 
!     averaging the values for KE and KE+1.
!---------------------------------------------------------------------
      do j = 1,size(trans_band2(:,:,:,:),2)
         do i = 1,size(trans_band2(:,:,:,:),1)
            trans_band2(i,j,KE+1,1) = 0.5E+00*(trans_band2(i,j,KE,1) +  &
                          trans_band2(i,j,KE+1, 1))
         end do
      end do
      ttmp(:,:,:) = trans_band2(:,:,:,1)
      do kk = KS+1,KE
         do j = 1,size(trans_band2(:,:,:,:),2)
            do i = 1,size(trans_band2(:,:,:,:),1)
               trans_band2(i,j,kk,1) = ttmp(i,j,kk-1)
            end do
         end do
      end do
 
!---------------------------------------------------------------------
!     perform calculations for the e1 function. the terms involving top
!     layer du are not known.  we use index two to represent index one
!     in previous calculations. (now 3)
!--------------------------------------------------------------------
      do j = 1,size(iyop(:,:,:),2)
         do i = 1,size(iyop(:,:,:),1)
            iyop(i,j,KS) = 2
         end do
      end do
      do kk = KS+1,KE+1
         do j = 1,size(iyop(:,:,:),2)
            do i = 1,size(iyop(:,:,:),1)
               iyop(i,j,kk) = iyo(i,j,kk+KS-(KS+1))
            end do
         end do
      end do
      do j = 1,size(dup(:,:,:),2)
         do i = 1,size(dup(:,:,:),1)
            dup(i,j,KS) = 0.0E+00
         end do
      end do
      do kk = KS+1,KE+1
         do j = 1,size(dup(:,:,:),2)
            do i = 1,size(dup(:,:,:),1)
               dup(i,j,kk) = du (i,j,kk+KS-(KS+1))
            end do
         end do
      end do
      do k=KS,KE+1
        do j = 1,size(ixo1(:,:,:),2)
           do i = 1,size(ixo1(:,:,:),1)
              ixo1(i,j,k) = ixoe1(i,j,KS)
           end do
        end do
        do j = 1,size(dt1(:,:,:),2)
           do i = 1,size(dt1(:,:,:),1)
              dt1(i,j,k) = dte1 (i,j,KS)
           end do
        end do
      enddo

    if (Lw_control%do_h2o) then
!-----------------------------------------------------------------------
!     e1flx(:,:,KS) equals e1cts1(:,:,KS).
!-----------------------------------------------------------------------
      call looktab (tab1, ixoe1, iyop, dte1, dup, e1cts1, KS, KE+1)
      call looktab (tab1, ixoe1, iyo, dte1, du, e1cts2, KS, KE)
      call looktab (tab1, ixo1, iyop, dt1, dup, &
                    trans_band1(:,:,:,1), KS, KE+1)
      call looktab (tab1w, ixoe1, iyop, dte1, dup, e1ctw1, KS, KE+1)
      call looktab (tab1w, ixoe1, iyo, dte1, du, e1ctw2, KS, KE)
    else
      iyop1(:,:,:) = 1
      dup1 (:,:,:) = 0.0
      call looktab (tab1, ixoe1, iyop1, dte1, dup1, e1cts1, KS, KE+1)
      call looktab (tab1, ixoe1, iyo1, dte1, du1, e1cts2, KS, KE)
      call looktab (tab1, ixo1, iyop1, dt1, dup1, &
                    trans_band1(:,:,:,1), KS, KE+1)
      call looktab (tab1w, ixoe1, iyop1, dte1, dup1, e1ctw1, KS, KE+1)
      call looktab (tab1w, ixoe1, iyo1, dte1, du1, e1ctw2, KS, KE)
    endif

!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o.
!--------------------------------------------------------------------
      if (NBTRGE > 0) then
        do m=1,NBTRGE
          if (Lw_control%do_h2o) then
          do kk = KS,KE+1
             do j = 1,size(avphilog(:,:,:),2)
                do i = 1,size(avphilog(:,:,:),1)
                   avphilog(i,j,kk) = LOG10(Optical%avephif(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, avphilog, du, iyo, KS , KE+1)
          iyo(:,:,KS:KE+1) = iyo(:,:,KS:KE+1) + 1
          do j = 1,size(iyop(:,:,:),2)
             do i = 1,size(iyop(:,:,:),1)
                iyop(i,j,KS) = 2
             end do
          end do
          do kk = KS+1,KE+1
             do j = 1,size(iyop(:,:,:),2)
                do i = 1,size(iyop(:,:,:),1)
                   iyop(i,j,kk) = iyo(i,j,kk+KS-(KS+1))
                end do
             end do
          end do
          do j = 1,size(dup(:,:,:),2)
             do i = 1,size(dup(:,:,:),1)
                dup(i,j,KS) = 0.0E+00
             end do
          end do
          do kk = KS+1,KE+1
             do j = 1,size(dup(:,:,:),2)
                do i = 1,size(dup(:,:,:),1)
                   dup(i,j,kk) = du (i,j,kk+KS-(KS+1))
                end do
             end do
          end do
          call looktab (tab2a, ixoe2, iyo, dte2, du, &
                        trans_band2(:,:,:,6+m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyop, dte1, dup,    &
                        e1cts1f(:,:,:,m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyo, dte1, du, e1cts2f(:,:,:,m),&
                        KS, KE, m)
          call looktab (tab1a, ixo1, iyop, dt1, dup,   &
                        trans_band1(:,:,:,6+m), KS, KE+1, m)
         else
          iyo1(:,:,KS:KE+1) = 1
          iyop1(:,:,KS:KE+1) = 1
          du1(:,:,KS:KE+1) = 0.0
          dup1(:,:,KS:KE+1) = 0.0
          call looktab (tab2a, ixoe2, iyo1, dte2, du1, &
                        trans_band2(:,:,:,6+m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyop1, dte1, dup1,    &
                        e1cts1f(:,:,:,m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyo1, dte1, du1, &
                        e1cts2f(:,:,:,m), KS, KE, m)
          call looktab (tab1a, ixo1, iyop1, dt1, dup1,   &
                        trans_band1(:,:,:,6+m), KS, KE+1, m)
        endif
        enddo

!--------------------------------------------------------------------
!     the special case emissf(:,:,KE+1,m) for layer KE+1 is obtained by 
!     averaging the values for KE and KE+1.
!--------------------------------------------------------------------
        do m=1,NBTRGE
         ttmp(:,:,:) = trans_band2(:,:,:,6+m)
         ttmp0(:,:,:) = trans_band2(:,:,:,6+m)
          do j = 1,size(trans_band2(:,:,:,:),2)
             do i = 1,size(trans_band2(:,:,:,:),1)
                trans_band2(i,j,KE+1,6+m) = 0.5E+00*    &
                             (ttmp(i,j,KE) +  &
                                   ttmp0(i,j,KE+1))
             end do
          end do
      ttmp(:,:,:)=trans_band2(:,:,:,6+m)
      do kk = KS+1,KE
         do j = 1,size(trans_band2(:,:,:,:),2)
            do i = 1,size(trans_band2(:,:,:,:),1)
               trans_band2(i,j,kk,6+m) = ttmp(i,j,kk+KS-(KS+1))
            end do
         end do
      end do
     enddo
    endif

!---------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
   if (NBTRGE > 0) then
   if (Lw_control%do_ch4 .or. Lw_control%do_n2o) then
     do kk = KS+1,KE+1
        do j = 1,size(trans_band1(:,:,:,:),2)
           do i = 1,size(trans_band1(:,:,:,:),1)
              trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*  &
                                tch4n2oe(i,j,kk,1)
           end do
        end do
     end do
     do kk = KS,KE+1
        do j = 1,size(e1cts1f(:,:,:,:),2)
           do i = 1,size(e1cts1f(:,:,:,:),1)
              e1cts1f(i,j,kk,1) = e1cts1f(i,j,kk,1)*  &
                                tch4n2oe(i,j,kk,1)
           end do
        end do
     end do
     do kk = KS,KE
        do j = 1,size(e1cts2f(:,:,:,:),2)
           do i = 1,size(e1cts2f(:,:,:,:),1)
              e1cts2f(i,j,kk,  1) = e1cts2f(i,j,kk,1)*   &
                                tch4n2oe(i,j,kk+KS+1-(KS),1)
           end do
        end do
     end do
     do kk = KS+1,KE+1
        do j = 1,size(trans_band2(:,:,:,:),2)
           do i = 1,size(trans_band2(:,:,:,:),1)
              trans_band2(i,j,kk,  6+1) = trans_band2(i,j,kk,6+1)*   &
                                tch4n2oe(i,j,kk,1)
           end do
        end do
     end do
   endif
 
!----------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range ( presently index 8) are present.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_indx8 (8, Optical, tcfc8)
       do kk = KS+1,KE+1
          do j = 1,size(trans_band1(:,:,:,:),2)
             do i = 1,size(trans_band1(:,:,:,:),1)
                trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*  &
                                  tcfc8(i,j,kk)
             end do
          end do
       end do
       do kk = KS,KE+1
          do j = 1,size(e1cts1f(:,:,:,:),2)
             do i = 1,size(e1cts1f(:,:,:,:),1)
                e1cts1f(i,j,kk,1) = e1cts1f(i,j,kk,1)*  &
                                  tcfc8(i,j,kk)
             end do
          end do
       end do
       do kk = KS,KE
          do j = 1,size(e1cts2f(:,:,:,:),2)
             do i = 1,size(e1cts2f(:,:,:,:),1)
                e1cts2f(i,j,kk,  1) = e1cts2f(i,j,kk,1)*  &
                                  tcfc8(i,j,kk+KS+1-(KS))
             end do
          end do
       end do
       do kk = KS+1,KE+1
          do j = 1,size(trans_band2(:,:,:,:),2)
             do i = 1,size(trans_band2(:,:,:,:),1)
                trans_band2(i,j,kk,6+1) = trans_band2(i,j,kk,6+1)*   &
                                  tcfc8(i,j,kk)
             end do
          end do
       end do
     endif 

!----------------------------------------------------------------------
!    compute aerosol transmission function for 1200-1400 cm-1 region
!----------------------------------------------------------------------
      if (including_aerosols) then
        do kk = 1,size(totaer_tmp(:,:,:),3)
           do j = 1,size(totaer_tmp(:,:,:),2)
              do i = 1,size(totaer_tmp(:,:,:),1)
                 totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,8)
              end do
           end do
        end do
       do kk = KS,KE+1
          do j = 1,size(taero8(:,:,:),2)
             do i = 1,size(taero8(:,:,:),1)
                taero8(i,j,kk) = EXP(-1.0E+00*totaer_tmp(i,j,kk))
             end do
          end do
       end do
       do kk = KS+1,KE+1
          do j = 1,size(trans_band1(:,:,:,:),2)
             do i = 1,size(trans_band1(:,:,:,:),1)
                trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*   &
                                  taero8(i,j,kk)
             end do
          end do
       end do
       do kk = KS,KE+1
          do j = 1,size(e1cts1f(:,:,:,:),2)
             do i = 1,size(e1cts1f(:,:,:,:),1)
                e1cts1f(i,j,kk,1) = e1cts1f(i,j,kk,1)*  &
                                  taero8(i,j,kk)
             end do
          end do
       end do
       do kk = KS,KE
          do j = 1,size(e1cts2f(:,:,:,:),2)
             do i = 1,size(e1cts2f(:,:,:,:),1)
                e1cts2f(i,j,kk,  1) = e1cts2f(i,j,kk,1)*   &
                                  taero8(i,j,kk+KS+1-(KS))
             end do
          end do
       end do
       do kk = KS+1,KE+1
          do j = 1,size(trans_band2(:,:,:,:),2)
             do i = 1,size(trans_band2(:,:,:,:),1)
                trans_band2(i,j,kk,6+1) = trans_band2(i,j,kk,6+1)*   &
                                  taero8(i,j,kk)
             end do
          end do
       end do
     endif
   endif  ! (NBTRGE > 0) 

!-----------------------------------------------------------------------
!     obtain the flux at the top of the atmosphere in the 0-160, 
!     1200-2200 cm-1 frequency ranges, where heating rates and fluxes
!     are derived from h2o emissivity calculations (flx1e1) by:
!     1) obtaining the surface flux (flxge1); 2) summing the
!     emissivity flux divergence for these ranges (tmp1) over all 
!     pressure layers.
!#ifdef ch4n2o
!     if the 1200-1400 cm-1 range is computed separately, flux calcu-
!     lations are done separately in this range, then combined with
!     those from the other frequency range.
!#endif ch4n2o
!----------------------------------------------------------------------
    do j = 1,size(s1a(:,:),2)
       do i = 1,size(s1a(:,:),1)
          s1a(i,j) = t4(i,j,KE+1)*(e1cts1(i,j,KE+1) - e1ctw1(i,j,KE+1))
       end do
    end do
    do j = 1,size(flxge1(:,:),2)
       do i = 1,size(flxge1(:,:),1)
          flxge1(i,j) = s1a(i,j)*cldtf(i,j,KE+1,1)
       end do
    end do
    do kk = KS,KE
       do j = 1,size(tmp1(:,:,:),2)
          do i = 1,size(tmp1(:,:,:),1)
             tmp1(i,j,kk) = t4(i,j,kk)*    &
                      (e1cts1(i,j,kk) - e1ctw1(i,j,kk)) 
          end do
       end do
    end do
    do kk = KS,KE
       do j = 1,size(tmp2(:,:,:),2)
          do i = 1,size(tmp2(:,:,:),1)
             tmp2(i,j,kk) = t4(i,j,kk)*   &
                       (e1cts2(i,j,kk) - e1ctw2(i,j,kk))
          end do
       end do
    end do
    do j = 1,size(Lw_diagnostics%flx1e1(:,:),2)
       do i = 1,size(Lw_diagnostics%flx1e1(:,:),1)
          Lw_diagnostics%flx1e1(i,j) = flxge1(i,j)
       end do
    end do
    do k=KS,KE
      do j = 1,size(Lw_diagnostics%flx1e1(:,:),2)
         do i = 1,size(Lw_diagnostics%flx1e1(:,:),1)
            Lw_diagnostics%flx1e1(i,j) = Lw_diagnostics%flx1e1(i,j) + tmp1(i,j,k)*cldtf(i,j,k,1) -   &
                    tmp2(i,j,k)*cldtf(i,j,k+1,1)
         end do
      end do
    enddo
    if (Rad_control%do_totcld_forcing) then
      do j = 1,size(flxge1cf(:,:),2)
         do i = 1,size(flxge1cf(:,:),1)
            flxge1cf(i,j) = s1a(i,j)
         end do
      end do
      do j = 1,size(flx1e1cf(:,:),2)
         do i = 1,size(flx1e1cf(:,:),1)
            flx1e1cf(i,j) = flxge1cf(i,j)
         end do
      end do
      do k=KS,KE
        do j = 1,size(flx1e1cf(:,:),2)
           do i = 1,size(flx1e1cf(:,:),1)
              flx1e1cf(i,j) = flx1e1cf(i,j) + tmp1(i,j,k) - tmp2(i,j,k)
           end do
        end do
      enddo
    endif
    if (NBTRGE > 0) then
      do m=1,NBTRGE
        do j = 1,size(s1a(:,:),2)
           do i = 1,size(s1a(:,:),1)
              s1a(i,j) = t4(i,j,KE+1)*e1cts1f(i,j,KE+1,m)
           end do
        end do
        do j = 1,size(flxge1f(:,:,:),2)
           do i = 1,size(flxge1f(:,:,:),1)
              flxge1f(i,j,m) = s1a(i,j)*cldtf(i,j,KE+1,cld_indx(7))
           end do
        end do
        do j = 1,size(Lw_diagnostics%flx1e1f(:,:,:),2)
           do i = 1,size(Lw_diagnostics%flx1e1f(:,:,:),1)
              Lw_diagnostics%flx1e1f(i,j,m) = flxge1f(i,j,m)
           end do
        end do
        do k=KS,KE
          do j = 1,size(tmp1(:,:,:),2)
             do i = 1,size(tmp1(:,:,:),1)
                tmp1(i,j,k) = t4(i,j,k)*e1cts1f(i,j,k,m)
             end do
          end do
          do j = 1,size(tmp2(:,:,:),2)
             do i = 1,size(tmp2(:,:,:),1)
                tmp2(i,j,k) = t4(i,j,k)*e1cts2f(i,j,k,m)
             end do
          end do
          do j = 1,size(Lw_diagnostics%flx1e1f(:,:,:),2)
             do i = 1,size(Lw_diagnostics%flx1e1f(:,:,:),1)
                Lw_diagnostics%flx1e1f(i,j,m) =  &
                       Lw_diagnostics%flx1e1f(i,j,m) + tmp1(i,j,k)*   &
                           cldtf(i,j,k,cld_indx(7)) - tmp2(i,j,k)*  &
                           cldtf(i,j,k+1,cld_indx(7))
             end do
          end do
        end do
      end do
      do m=1,NBTRGE
        do j = 1,size(Lw_diagnostics%flx1e1(:,:),2)
           do i = 1,size(Lw_diagnostics%flx1e1(:,:),1)
              Lw_diagnostics%flx1e1(i,j) =    &
                                 Lw_diagnostics%flx1e1(i,j) +  &
                                 Lw_diagnostics%flx1e1f(i,j,m)
           end do
        end do
      enddo
      if (Rad_control%do_totcld_forcing) then
        do m=1,NBTRGE
          do j = 1,size(flxge1fcf(:,:,:),2)
             do i = 1,size(flxge1fcf(:,:,:),1)
                flxge1fcf(i,j,m) = s1a(i,j)
             end do
          end do
          do j = 1,size(flx1e1fcf(:,:,:),2)
             do i = 1,size(flx1e1fcf(:,:,:),1)
                flx1e1fcf(i,j,m) = s1a(i,j)
             end do
          end do
          do k=KS,KE
            do j = 1,size(flx1e1fcf(:,:,:),2)
               do i = 1,size(flx1e1fcf(:,:,:),1)
                  flx1e1fcf(i,j,m) = flx1e1fcf(i,j,m) + tmp1(i,j,k) -   &
                                                    tmp2(i,j,k)
               end do
            end do
          end do
        end do
        do m=1,NBTRGE
          do j = 1,size(flx1e1cf(:,:),2)
             do i = 1,size(flx1e1cf(:,:),1)
                flx1e1cf(i,j) = flx1e1cf(i,j) + flx1e1fcf(i,j,m)
             end do
          end do
        enddo
      endif
    endif  ! (ntrge > 0)

!----------------------------------------------------------------------



end  subroutine e1e290



!###################################################################
! <SUBROUTINE NAME="e290">
!  <OVERVIEW>
!   e290 computes the exchange terms in the flux equation for longwave
!     radiation for all terms except the exchange with the top of the
!     atmosphere.
!  </OVERVIEW>
!  <DESCRIPTION>
!     e290 computes the exchange terms in the flux equation for longwave
!     radiation for all terms except the exchange with the top of the
!     atmosphere.  the method is a table lookup on a pre-computed e2
!     function (defined in reference (2)).  calculation are done in the
!     frequency range: 0-560, 1200-2200 cm-1 for q(approximate).
!     motivation for these calculations is in references (1) and (2).
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call e290 (Atmos_input, k, trans_band2, trans_band1, Optical,  tch4n2oe, &
!              tcfc8)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data
!  </IN>
!  <IN NAME="k" TYPE="integer">
!   Starting vertical level k to compute exchange terms
!  </IN>
!  <INOUT NAME="trans_band2" TYPE="real">
!   Transmission funciton in band 1200-2200
!  </INOUT>
!  <INOUT NAME="trans_band" TYPE="real">
!   Transmission function in band 0-560
!  </INOUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth of the thermal layers
!  </INOUT>
!  <IN NAME="tch4n2oe" TYPE="real">
!   CH4 and N2O transmission function
!  </IN>
!  <INOUT NAME="tcfc8" TYPE="real">
!   CFC transmission function
!  </INOUT>
! </SUBROUTINE>
!
subroutine e290 (Atmos_input, k, trans_band2, trans_band1, Optical,  &
                 tch4n2oe, tcfc8, including_aerosols)   

!-----------------------------------------------------------------------
!
!     e290 computes the exchange terms in the flux equation for longwave
!     radiation for all terms except the exchange with the top of the
!     atmosphere.  the method is a table lookup on a pre-computed e2
!     function (defined in reference (2)).  calculation are done in the
!     frequency range: 0-560, 1200-2200 cm-1 for q(approximate).
!     motivation for these calculations is in references (1) and (2).
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!---------------------------------------------------------------------

type(atmos_input_type),   intent(in)    ::  Atmos_input
integer,                  intent(in)    ::  k
real, dimension(:,:,:,:), intent(inout) :: trans_band1, trans_band2
type(optical_path_type),  intent(inout) ::  Optical
real, dimension(:,:,:,:), intent(in)    :: tch4n2oe       
real, dimension(:,:,:),   intent(inout) :: tcfc8          
logical,                   intent(in)            :: including_aerosols  
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!   intent(in) variables:
!
!       Atmos_input
!       k
!       tch4n2oe
!
!   intent(inout) variables:
!
!       trans_band1
!       trans_band2
!       Optical
!       tcfc8
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables
!-------------------------------------------------------------------
      integer      :: i,j,kk,kp, m

      real,    dimension (size(Atmos_input%temp,1),    &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)) ::    &
                                  temp, tflux, totaer_tmp, taero8, &
                                  taero8kp, avphilog, dtk, du, du1

      integer, dimension (size(Atmos_input%temp,1),    &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)) ::              &
                                                  ixok, iyo, iyo1    

      real,    dimension (size(trans_band2,1), &
                          size(trans_band2,2), &
                          size(trans_band2,3)) :: dte1, dte2

      integer, dimension (size(trans_band2,1), &
                          size(trans_band2,2), &
                          size(trans_band2,3)) :: ixoe1, ixoe2
      real,    dimension (size(trans_band2,1), &
                          size(trans_band2,2), &
                          size(trans_band2,3)) ::ttmp
      integer,    dimension (size(trans_band2,1), &
                          size(trans_band2,2), &
                          size(trans_band2,3)) ::ittmp      

      real,    dimension (size(trans_band1,1), &
                          size(trans_band1,2), &
                          size(trans_band1,3)) :: ttmp0

!-------------------------------------------------------------------
!  local variables:
!
!      temp
!      tflux
!      totaer_tmp
!      taero8
!      taero8kp
!      avphilog
!      dtk
!      du
!      du1
!      ixok
!      iyo
!      iyo1
!      dte1
!      dte2
!      ixoe1
!      ixoe2
!      kp,m
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      temp = Atmos_input%temp
      tflux = Atmos_input%tflux

!-----------------------------------------------------------------------
!     obtain the "exchange" emissivities as a function of temperature 
!     (fxo) and water amount (avephi). the temperature indices have
!     been obtained in Lwrad. calculations are for flux level k, with
!     kp = k+1 to KE+1. the case k=KS is excluded (done in E1e290).
!     calculations are also made for flux levels k to KE, for
!     contributions from flux level k. in this case, the temperature
!     index (ixok) represents tflux(:,:,k-1); the water index (iyo)
!     has the same values as in the case with varying kp.
!---------------------------------------------------------------------

!----------------------------------------------------------------------

  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ittmp(:,:,:) = ixoe2(:,:,:)
  do kk = KS,KE
     do j = 1,size(ixoe2(:,:,:),2)
        do i = 1,size(ixoe2(:,:,:),1)
           ixoe2(i,j,kk) = ittmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  ttmp(:,:,:)=dte2 (:,:,:)
  do kk = KS,KE
     do j = 1,size(dte2(:,:,:),2)
        do i = 1,size(dte2(:,:,:),1)
           dte2(i,j,kk) = ttmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  do j = 1,size(ixoe2(:,:,:),2)
     do i = 1,size(ixoe2(:,:,:),1)
        ixoe2(i,j,KE+1) = ixoe1(i,j,KE)
     end do
  end do
  do j = 1,size(dte2(:,:,:),2)
     do i = 1,size(dte2(:,:,:),1)
        dte2(i,j,KE+1) = dte1 (i,j,KE)
     end do
  end do


      do kp=k,KE
        do j = 1,size(ixok(:,:,:),2)
           do i = 1,size(ixok(:,:,:),1)
              ixok(i,j,kp) = ixoe2(i,j,k-1)
           end do
        end do
        do j = 1,size(dtk(:,:,:),2)
           do i = 1,size(dtk(:,:,:),1)
              dtk(i,j,kp) = dte2 (i,j,k-1)
           end do
        end do
      end do

   if (Lw_control%do_h2o) then
      do kk = k,KE+1
         do j = 1,size(avphilog(:,:,:),2)
            do i = 1,size(avphilog(:,:,:),1)
               avphilog(i,j,kk) = LOG10(Optical%avephi(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, avphilog, du, iyo,k, KE+1)
      iyo(:,:,k:KE+1) = iyo(:,:,k:KE+1) + 1
      call looktab (tab2, ixoe2, iyo, dte2, du, &
                    trans_band2(:,:,:,1), k, KE+1)
      call looktab (tab2, ixok, iyo, dtk, du, &
                    trans_band1(:,:,:,1), k, KE)
   else
     iyo1(:,:,k:KE+1) = 1
     du1(:,:,k:KE+1) = 0.0
     call looktab (tab2, ixoe2, iyo1, dte2, du1, &
                   trans_band2(:,:,:,1), k, KE+1)
     call looktab (tab2, ixok, iyo1, dtk, du1, &
                   trans_band1(:,:,:,1), k, KE)
   endif
       ttmp0(:,:,:)=trans_band1(:,:,:,1)
       do kk = k+1,KE+1
          do j = 1,size(trans_band1(:,:,:,:),2)
             do i = 1,size(trans_band1(:,:,:,:),1)
                trans_band1(i,j,kk,1) = ttmp0(i,j,kk+k-(k+1))
             end do
          end do
       end do

!--------------------------------------------------------------------
!     the special case emiss(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emiss(:,:,KE+1) 
!     is not useful after this point.
!-------------------------------------------------------------------
      do j = 1,size(trans_band2(:,:,:,:),2)
         do i = 1,size(trans_band2(:,:,:,:),1)
            trans_band2(i,j,KE+1,1) = 0.5E+00*(trans_band2(i,j,KE,1) +  &
                                               trans_band2(i,j,KE+1,1))
         end do
      end do
      ttmp(:,:,:)=trans_band2(:,:,:,1)
      do kk = k+1,KE
         do j = 1,size(trans_band2(:,:,:,:),2)
            do i = 1,size(trans_band2(:,:,:,:),1)
               trans_band2(i,j,kk,1) = ttmp(i,j,kk+k-(k+1))
            end do
         end do
      end do
 
!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o. reqults are in emissf (flux level k) and
!     emissbf (other levels).
!-------------------------------------------------------------------
      if (nbtrge > 0) then
        do m=1,NBTRGE
   if (Lw_control%do_h2o) then
          do kk = k,KE+1
             do j = 1,size(avphilog(:,:,:),2)
                do i = 1,size(avphilog(:,:,:),1)
                   avphilog(i,j,kk) = LOG10(Optical%avephif(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, avphilog, du, iyo, k, KE+1)
          iyo(:,:,k:KE+1) = iyo(:,:,k:KE+1) + 1
          call looktab (tab2a, ixoe2, iyo, dte2, du, &
                        trans_band2(:,:,:,6+m), k, KE+1, m)
          call looktab (tab2a, ixok, iyo, dtk, du, &
                        trans_band1(:,:,:,6+m), k, KE, m)
   else
      iyo1(:,:,k:KE+1) = 1
      du1(:,:,k:KE+1) = 0.0
      call looktab (tab2a, ixoe2, iyo1, dte2, du1, &
                    trans_band2(:,:,:,6+m), k, KE+1, m)
      call looktab (tab2a, ixok, iyo1, dtk, du1, &
                    trans_band1(:,:,:,6+m), k, KE, m)
   endif
       ttmp0(:,:,:)=trans_band1(:,:,:,6+m)
       do kk = k+1,KE+1
          do j = 1,size(trans_band1(:,:,:,:),2)
             do i = 1,size(trans_band1(:,:,:,:),1)
                trans_band1(i,j,kk,6+m) = ttmp0(i,j,kk+k-(k+1))
             end do
          end do
       end do
      enddo

!----------------------------------------------------------------------
!     the special case emissf(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emissf(:,:,KE+1,m)
!     is not useful after this point.
!----------------------------------------------------------------------
        do m=1,NBTRGE
          do j = 1,size(trans_band2(:,:,:,:),2)
             do i = 1,size(trans_band2(:,:,:,:),1)
                trans_band2(i,j,KE+1,6+m) = 0.5E+00*  &
                             (trans_band2(i,j,KE,6+m) +   &
                               trans_band2(i,j,KE+1,6+m))
             end do
          end do
          ttmp(:,:,:)=trans_band2(:,:,:,6+m)
          do kk = k+1,KE
             do j = 1,size(trans_band2(:,:,:,:),2)
                do i = 1,size(trans_band2(:,:,:,:),1)
                   trans_band2(i,j,kk,6+m) = ttmp(i,j,kk+k-(k+1))
                end do
             end do
          end do
        enddo
      endif

!----------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!---------------------------------------------------------------------
    if (nbtrge > 0) then
      if (Lw_control%do_ch4 .or. Lw_control%do_n2o) then
        do kk = k+1,KE+1
           do j = 1,size(trans_band1(:,:,:,:),2)
              do i = 1,size(trans_band1(:,:,:,:),1)
                 trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*  &
                         tch4n2oe(i,j,kk,1)
              end do
           end do
        end do
        do kk = k+1,KE+1
           do j = 1,size(trans_band2(:,:,:,:),2)
              do i = 1,size(trans_band2(:,:,:,:),1)
                 trans_band2(i,j,kk,6+1) = trans_band2(i,j,kk,6+1)*tch4n2oe(i,j,kk,1)
              end do
           end do
        end do
      endif 

!--------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!--------------------------------------------------------------------
        if (Lw_control%do_cfc) then
          call cfc_indx8_part (8, Optical, tcfc8, k)
          do kk = k+1,KE+1
             do j = 1,size(trans_band1(:,:,:,:),2)
                do i = 1,size(trans_band1(:,:,:,:),1)
                   trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*tcfc8(i,j,kk)
                end do
             end do
          end do
          do kk = k+1,KE+1
             do j = 1,size(trans_band2(:,:,:,:),2)
                do i = 1,size(trans_band2(:,:,:,:),1)
                   trans_band2(i,j,kk,6+1) = trans_band2(i,j,kk,6+1)*tcfc8(i,j,kk)
                end do
             end do
          end do
        endif

!--------------------------------------------------------------------
!     compute aerosol transmission function for 1200-1400 cm-1 region
!    (as quotient of 2 exponentials)
!     taero8kp(k) contains the (k+1,k) transmissivities for all k
!     in the 1200-1400 cm-1 frequency range.
!---------------------------------------------------------------------
      if (including_aerosols) then
        do kk = 1,size(totaer_tmp(:,:,:),3)
           do j = 1,size(totaer_tmp(:,:,:),2)
              do i = 1,size(totaer_tmp(:,:,:),1)
                 totaer_tmp(i,j,kk) = Optical%totaerooptdep(i,j,kk,8)
              end do
           end do
        end do
       do kk = KS,KE+1
          do j = 1,size(taero8(:,:,:),2)
             do i = 1,size(taero8(:,:,:),1)
                taero8(i,j,kk) = EXP(-1.0E+00*totaer_tmp(i,j,kk))
             end do
          end do
       end do
          do kp = k+1,KE+1
            do j = 1,size(taero8kp(:,:,:),2)
               do i = 1,size(taero8kp(:,:,:),1)
                  taero8kp(i,j,kp) = taero8(i,j,kp)/taero8(i,j,k)
               end do
            end do
          enddo
          do kk = k+1,KE+1
             do j = 1,size(trans_band1(:,:,:,:),2)
                do i = 1,size(trans_band1(:,:,:,:),1)
                   trans_band1(i,j,kk,6+1) = trans_band1(i,j,kk,6+1)*  &
                  taero8kp(i,j,kk)
                end do
             end do
          end do
          do kk = k+1,KE+1
             do j = 1,size(trans_band2(:,:,:,:),2)
                do i = 1,size(trans_band2(:,:,:,:),1)
                   trans_band2(i,j,kk,6+1) = trans_band2(i,j,kk,6+1)*  &
                  taero8kp(i,j,kk)
                end do
             end do
          end do
        endif
      endif  ! (nbtrge > 0)

!--------------------------------------------------------------------

end subroutine e290



!####################################################################
! <SUBROUTINE NAME="esfc">
!  <OVERVIEW>
!   Subroutine to compute thermal layer emissivity using pre computed
!   look up tables
!  </OVERVIEW>
!  <TEMPLATE>
!   call esfc  (Atmos_input,         emspec,             Optical, &
!               emspecf, tch4n2oe, tcfc8 )
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data such as temperature and flux level temp
!  </IN>
!  <OUT NAME="emspec" TYPE="real">
!   Emissivity of thermal layers
!  </OUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth of thermal layers
!  </INOUT>
!  <OUT NAME="emspecf" TYPE="real">
!   Emissivity of thermal layers including effects of minor gas species
!  </OUT>
!  <IN NAME="tch4n2oe" TYPE="real">
!   CH4 and N2O transmission function
!  </IN>
!  <INOUT NAME="tcfc8" TYPE="real">
!   CFC transmission function
!  </INOUT>
! </SUBROUTINE>
!
subroutine esfc (Atmos_input, emspec, Optical, emspecf, tch4n2oe, tcfc8 ) 
   
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

type(atmos_input_type),    intent(in)    :: Atmos_input
real, dimension (:,:,:),   intent(out)   :: emspec
type(optical_path_type),   intent(inout) :: Optical
real, dimension (:,:,:,:), intent(out)   :: emspecf
real, dimension (:,:,:,:), intent(in)    :: tch4n2oe
real, dimension (:,:,:),   intent(inout) :: tcfc8   
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  intent(in)variables:
!    
!     Atmos_input
!     tch4n2oe
!
!  intent(inout) variables:
!
!     Optical
!     tcfc8
!
!  intent(out) variables:
!
!     emspec
!     emspecf
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      integer :: i,j,k,kk,m



      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) :: &
                                     temp, tflux, tpl1, tpl2, dte1, &
                                     dte2, dxsp, ylog, dysp, emiss, &
                                     emd1, emd2, dysp1

      integer, dimension (size(Atmos_input%temp,1),   &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)) :: &
                                        ixsp, iysp, iysp1, ixoe1, ixoe2
      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2),   &
                       size(Atmos_input%temp,3)) :: ttmp
      integer, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2),   &
                       size(Atmos_input%temp,3)) :: ittmp
      
      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3), NBTRGE) ::    &
                                                  emissf, emd2f,  emd1f

!--------------------------------------------------------------------
!   local variables:
!
!      temp
!      tflux
!      tpl1
!      tpl2
!      dte1
!      dte2 
!      dxsp
!      ylog
!      dysp
!      dysp1
!      emiss
!      emd1
!      emd2
!      ixsp
!      iysp
!      iysp1
!      ixoe1
!      ixoe2
!      emissf
!      emd2f
!      emd1f
!      m,k
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

     do kk = 1,size(tflux(:,:,:),3)
        do j = 1,size(tflux(:,:,:),2)
           do i = 1,size(tflux(:,:,:),1)
              tflux(i,j,kk) = Atmos_input%tflux(i,j,kk)
           end do
        end do
     end do
      do kk = 1,size(temp(:,:,:),3)
         do j = 1,size(temp(:,:,:),2)
            do i = 1,size(temp(:,:,:),1)
               temp(i,j,kk) = Atmos_input%temp(i,j,kk)
            end do
         end do
      end do
 


      do j = 1,size(tpl1(:,:,:),2)
         do i = 1,size(tpl1(:,:,:),1)
            tpl1(i,j,KS) = temp(i,j,KE)
         end do
      end do
     do kk = KS+1,KE
        do j = 1,size(tpl1(:,:,:),2)
           do i = 1,size(tpl1(:,:,:),1)
              tpl1(i,j,kk) = tflux(i,j,kk)
           end do
        end do
     end do
     do j = 1,size(tpl1(:,:,:),2)
        do i = 1,size(tpl1(:,:,:),1)
           tpl1(i,j,KE+1) = 0.5E+00*(tflux(i,j,KE+1) +   &
                                     temp(i,j,KE))
        end do
     end do
  do kk = KS+1,KE
     do j = 1,size(tpl2(:,:,:),2)
        do i = 1,size(tpl2(:,:,:),1)
           tpl2(i,j,kk) = tflux(i,j,kk)
        end do
     end do
  end do
    do j = 1,size(tpl2(:,:,:),2)
       do i = 1,size(tpl2(:,:,:),1)
          tpl2(i,j,KE+1) = 0.5E+00*(tflux(i,j,KE) +    &
                                temp(i,j,KE))
       end do
    end do


!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ittmp(:,:,:)=ixoe2(:,:,:)
  do kk = KS,KE
     do j = 1,size(ixoe2(:,:,:),2)
        do i = 1,size(ixoe2(:,:,:),1)
           ixoe2(i,j,kk) = ittmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  ttmp(:,:,:)=dte2 (:,:,:)
  do kk = KS,KE
     do j = 1,size(dte2(:,:,:),2)
        do i = 1,size(dte2(:,:,:),1)
           dte2(i,j,kk) = ttmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  do j = 1,size(ixoe2(:,:,:),2)
     do i = 1,size(ixoe2(:,:,:),1)
        ixoe2(i,j,KE+1) = ixoe1(i,j,KE)
     end do
  end do
  do j = 1,size(dte2(:,:,:),2)
     do i = 1,size(dte2(:,:,:),1)
        dte2(i,j,KE+1) = dte1 (i,j,KE)
     end do
  end do




      do j = 1,size(ixsp(:,:,:),2)
         do i = 1,size(ixsp(:,:,:),1)
            ixsp(i,j,KE) = ixoe2(i,j,KE-1)
         end do
      end do
      do j = 1,size(ixsp(:,:,:),2)
         do i = 1,size(ixsp(:,:,:),1)
            ixsp(i,j,KE+1) = ixoe1(i,j,KE-1)
         end do
      end do
      do j = 1,size(dxsp(:,:,:),2)
         do i = 1,size(dxsp(:,:,:),1)
            dxsp(i,j,KE) = dte2(i,j,KE-1)
         end do
      end do
      do j = 1,size(dxsp(:,:,:),2)
         do i = 1,size(dxsp(:,:,:),1)
            dxsp(i,j,KE+1) = dte1(i,j,KE-1)
         end do
      end do

    if (Lw_control%do_h2o) then
      do j = 1,size(ylog(:,:,:),2)
         do i = 1,size(ylog(:,:,:),1)
            ylog(i,j,KE  ) = ALOG10(Optical%var2(i,j,KE))
         end do
      end do
      do j = 1,size(ylog(:,:,:),2)
         do i = 1,size(ylog(:,:,:),1)
            ylog(i,j,KE+1) = ALOG10(Optical%var2(i,j,KE) + Optical%empl1(i,j,KE))
         end do
      end do

      call locate_in_table (mass_1, ylog, dysp, iysp, KE, KE+1)
      iysp(:,:,KE:KE+1) = iysp(:,:,KE:KE+1) + 1

!--------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!--------------------------------------------------------------------
      call looktab (tab2, ixsp, iysp, dxsp, dysp, emiss, KE, KE+1)

    else
      iysp1(:,:,KE:KE+1) = 1
      dysp1(:,:,KE:KE+1) = 0.0
      call looktab (tab2, ixsp, iysp1, dxsp, dysp1, emiss, KE, KE+1)
    endif

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (nbtrge > 0) then
        do m=1,NBTRGE
    if (Lw_control%do_h2o) then
          do j = 1,size(ylog(:,:,:),2)
             do i = 1,size(ylog(:,:,:),1)
                ylog(i,j,KE  ) = ALOG10(Optical%vrpfh2o(i,j,KE,m))
             end do
          end do
          do j = 1,size(ylog(:,:,:),2)
             do i = 1,size(ylog(:,:,:),1)
                ylog(i,j,KE+1) = ALOG10(Optical%vrpfh2o(i,j,KE,m) + Optical%empl1f(i,j,KE,m))
             end do
          end do

          call locate_in_table (mass_1, ylog, dysp, iysp, KE, KE+1)
          iysp(:,:,KE:KE+1) = iysp(:,:,KE:KE+1) + 1

!-----------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!---------------------------------------------------------------------
          call looktab (tab2a, ixsp, iysp, dxsp, dysp, &
                        emissf(:,:,:,m), KE, KE+1, m)
   else
        iysp1(:,:,KE:KE+1) = 1
        dysp1(:,:,KE:KE+1) = 0.0
        call looktab (tab2a, ixsp, iysp1, dxsp, dysp1, &
                      emissf(:,:,:,m), KE, KE+1, m)
   endif
        enddo
      endif
!-----------------------------------------------------------------------
!     compute nearby layer transmissivities for h2o.
!--------------------------------------------------------------------
    if (Lw_control%do_h2o) then
      call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
      do kk = KS,KE+1
         do j = 1,size(ylog(:,:,:),2)
            do i = 1,size(ylog(:,:,:),1)
               ylog(i,j,kk) = ALOG10(Optical%empl1(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
      iysp(:,:,KS:KE+1) = iysp(:,:,KS:KE+1) + 1
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd1, KS, KE+1)
   else
     call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
     iysp1(:,:,KS:KE+1) = 1
     dysp1(:,:,KS:KE+1) = 0.0
     call looktab (tab3, ixsp, iysp1, dxsp, dysp1, emd1, KS, KE+1)
   endif

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!------------------------------------------------------------------
      if (nbtrge > 0) then
        do m=1,NBTRGE
   if (Lw_control%do_h2o) then
          do kk = KS,KE+1
             do j = 1,size(ylog(:,:,:),2)
                do i = 1,size(ylog(:,:,:),1)
                   ylog(i,j,kk) = ALOG10(Optical%empl1f(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
          iysp(:,:,KS:KE+1) = iysp(:,:,KS:KE+1) + 1
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, &
                        emd1f(:,:,:,m), KS, KE+1, m)
    else
      iysp1(:,:,KS:KE+1) = 1
      dysp1(:,:,KS:KE+1) = 0.0
      call looktab (tab3a, ixsp, iysp1, dxsp, dysp1, &
                    emd1f(:,:,:,m), KS, KE+1, m)
    endif
        enddo
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  if (Lw_control%do_h2o) then
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      do kk = KS+1,KE+1
         do j = 1,size(ylog(:,:,:),2)
            do i = 1,size(ylog(:,:,:),1)
               ylog(i,j,kk) = ALOG10(Optical%empl2(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
      iysp(:,:,KS+1:KE+1) = iysp(:,:,KS+1:KE+1) + 1
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd2, KS+1, KE+1)
  else
     call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
     iysp1(:,:,KS+1:KE+1) = 1
     dysp1(:,:,KS+1:KE+1) = 0.0
     call looktab (tab3, ixsp, iysp1, dxsp, dysp1, emd2, KS+1, KE+1)
   endif

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (nbtrge > 0 ) then
        do m=1,NBTRGE
    if (Lw_control%do_h2o) then
          do kk = KS+1,KE+1
             do j = 1,size(ylog(:,:,:),2)
                do i = 1,size(ylog(:,:,:),1)
                   ylog(i,j,kk) = ALOG10(Optical%empl2f(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
          iysp(:,:,KS+1:KE+1) = iysp(:,:,KS+1:KE+1) + 1
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, &
                        emd2f(:,:,:,m), KS+1, KE+1, m)
    else
      iysp1(:,:,KS+1:KE+1) = 1
      dysp1(:,:,KS+1:KE+1) = 0.0
      call looktab (tab3a, ixsp, iysp1, dxsp, dysp1, &
                    emd2f(:,:,:,m), KS+1, KE+1, m)
    endif
        enddo
      endif

!---------------------------------------------------------------------- 
!     compute nearby layer and special-case transmissivities for
!     emissivity using methods for h2o given in reference (4).
!-------------------------------------------------------------------- 
     if (Lw_control%do_h2o) then   
      do j = 1,size(emspec(:,:,:),2)
         do i = 1,size(emspec(:,:,:),1)
            emspec(i,j,KS     ) = (emd1(i,j,KS)*Optical%empl1(i,j,KS) -    &
                             emd1(i,j,KE+1)*Optical%empl1(i,j,KE+1))/  &
                             Optical%emx1(i,j) + 0.25E+00*(emiss(i,j,KE) +   &
                             emiss(i,j,KE+1))
         end do
      end do
      do j = 1,size(emspec(:,:,:),2)
         do i = 1,size(emspec(:,:,:),1)
            emspec(i,j,KS+1) = 2.0E+00*(emd1(i,j,KS)*Optical%empl1(i,j,KS) -    &
                         emd2(i,j,KE+1)*Optical%empl2(i,j,KE+1))/  &
                           Optical%emx2(i,j)
         end do
      end do
   else
     do j = 1,size(emspec(:,:,:),2)
       do i = 1,size(emspec(:,:,:),1)
         emspec(i,j,KS     ) = (emd1(i,j,KS) - emd1(i,j,KE+1)) /  &
                              (Atmos_input%press(i,j,KE) - Atmos_input%pflux(i,j,KE)) + &
                               0.25E+00*(emiss(i,j,KE) + emiss(i,j,KE+1))
        end do
      end do
      do j = 1,size(emspec(:,:,:),2)
        do i = 1,size(emspec(:,:,:),1)
            emspec(i,j,KS+1) = 2.0E+00*(emd1(i,j,KS) -  emd2(i,j,KE+1)) /  &
                      (Atmos_input%pflux(i,j,KE+1) - Atmos_input%press(i,j,KE))
        end do
      end do
    endif

     if (nbtrge > 0) then
       do m=1,NBTRGE
     if (Lw_control%do_h2o) then
         do j = 1,size(emspecf(:,:,:,:),2)
            do i = 1,size(emspecf(:,:,:,:),1)
               emspecf(i,j,KS,m   ) = (emd1f(i,j,KS,m)*Optical%empl1f(i,j,KS,m) -   &
                              emd1f(i,j,KE+1,m)*Optical%empl1f(i,j,KE+1,m))/   &
                   Optical%emx1f(i,j,m) + 0.25E+00*(emissf(i,j,KE,m) +  &
                             emissf(i,j,KE+1,m))
            end do
         end do
         do j = 1,size(emspecf(:,:,:,:),2)
            do i = 1,size(emspecf(:,:,:,:),1)
               emspecf(i,j,KS+1,m) = 2.0E+00*    &
                                 (emd1f(i,j,KS,m)*Optical%empl1f(i,j,KS,m) -  &
                              emd2f(i,j,KE+1,m)*Optical%empl2f(i,j,KE+1,m)) / &
                              Optical%emx2f(i,j,m)
            end do
         end do
     else
       do j = 1,size(emspecf(:,:,:,:),2)
           do i = 1,size(emspecf(:,:,:,:),1)
               emspecf(i,j,KS,m   ) = (emd1f(i,j,KS,m) - emd1f(i,j,KE+1,m)) / &
                            (Atmos_input%press(i,j,KE) - Atmos_input%pflux(i,j,KE)) + &
                               0.25E+00*(emissf(i,j,KE,m) +  emissf(i,j,KE+1,m))
            end do
          end do
          do j = 1,size(emspecf(:,:,:,:),2)
            do i = 1,size(emspecf(:,:,:,:),1)
              emspecf(i,j,KS+1,m) = 2.0E+00*    &
                                  (emd1f(i,j,KS,m) - emd2f(i,j,KE+1,m)) / &
                          (Atmos_input%pflux(i,j,KE+1) - Atmos_input%press(i,j,KE))
             end do
          end do
        endif

       enddo
     endif

!--------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
  if (nbtrge > 0) then
    if (Lw_control%do_ch4 .or. Lw_control%do_n2o) then
      do k=KS,KS+1
        do j = 1,size(emspecf(:,:,:,:),2)
           do i = 1,size(emspecf(:,:,:,:),1)
              emspecf(i,j,K,1) = emspecf(i,j,K,1)*tch4n2oe(i,j,KE+1,1)
           end do
        end do
      end do
     endif 

!--------------------------------------------------------------------
!     add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_indx8_part (8, Optical, tcfc8, KE)
        do k=KS,KS+1
          do j = 1,size(emspecf(:,:,:,:),2)
             do i = 1,size(emspecf(:,:,:,:),1)
                emspecf(i,j,K,1) = emspecf(i,j,K,1)*tcfc8(i,j,KE+1)
             end do
          end do
        end do
      endif
    endif ! (nbtrge > 0)

!------------------------------------------------------------------



end subroutine esfc 


!######################################################################
! <SUBROUTINE NAME="enear">
!  <OVERVIEW>
!   Subroutine to compute thermal layer emissivity using pre computed
!   look up tables
!  </OVERVIEW>
!  <TEMPLATE>
!   call enear  (Atmos_input,         emisdg,             Optical, &
!               emisdgf, tch4n2oe, tcfc8 )
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data such as temperature and flux level temp
!  </IN>
!  <OUT NAME="emisdg" TYPE="real">
!   Emissivity of thermal layers
!  </OUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth of thermal layers
!  </INOUT>
!  <OUT NAME="emisdgf" TYPE="real">
!   Emissivity of thermal layers including effects of minor gas species
!  </OUT>
!  <IN NAME="tch4n2oe" TYPE="real">
!   CH4 and N2O transmission function
!  </IN>
!  <INOUT NAME="tcfc8" TYPE="real">
!   CFC transmission function
!  </INOUT>
! </SUBROUTINE>
!
subroutine enear (Atmos_input, emisdg, Optical, emisdgf, tch4n2oe,  &
                  tcfc8) 
   
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

type(atmos_input_type),    intent(in)    ::  Atmos_input
real, dimension (:,:,:),   intent(out)   ::  emisdg 
type(optical_path_type),   intent(inout) ::  Optical
real, dimension (:,:,:,:), intent(out)   ::  emisdgf
real, dimension (:,:,:,:), intent(in)    ::  tch4n2oe
real, dimension (:,:,:),   intent(inout) ::  tcfc8       
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     Atmos_input
!     tch4n2oe
!
!  intent(inout) variables:
!
!     Optical
!     tcfc8
!
!  intent(out) variables:
!
!     emisdg
!     emisdgf
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer   :: i,j,kk,m

      real, dimension (size(emisdg,1), &
                       size(emisdg,2), &
                       size(emisdg,3)) :: dte1, dte2

      integer, dimension (size(emisdg,1), &
                          size(emisdg,2), &
                          size(emisdg,3)) :: ixoe1, ixoe2

      real, dimension (size(emisdg,1), &
                       size(emisdg,2), &
                       size(emisdg,3)) :: ttmp 

      integer, dimension (size(emisdg,1), &
                       size(emisdg,2), &
                       size(emisdg,3)) :: ittmp 
      
      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) :: &
                             temp, tflux, tpl1, tpl2, &
                             dxsp, ylog, dysp, dysp1, emd1, emd2

      integer, dimension (size(Atmos_input%temp,1),   &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)) :: &
                                                       ixsp, iysp, iysp1

      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3), NBTRGE) ::    &
                            emd2f,  emd1f

!--------------------------------------------------------------------
!   local variables:
!
!      dte1
!      dte2
!      ixoe1
!      ixoe2 
!      temp
!      tflux
!      tpl1
!      tpl2
!      dxsp
!      ylog
!      dysp
!      dysp1
!      emiss
!      emd1
!      emd2
!      ixsp
!      iysp
!      iysp1
!      emissf
!      emd2f
!      emd1f
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
     do kk = 1,size(tflux(:,:,:),3)
        do j = 1,size(tflux(:,:,:),2)
           do i = 1,size(tflux(:,:,:),1)
              tflux(i,j,kk) = Atmos_input%tflux(i,j,kk)
           end do
        end do
     end do
      do kk = 1,size(temp(:,:,:),3)
         do j = 1,size(temp(:,:,:),2)
            do i = 1,size(temp(:,:,:),1)
               temp(i,j,kk) = Atmos_input%temp(i,j,kk)
            end do
         end do
      end do
 


      do j = 1,size(tpl1(:,:,:),2)
         do i = 1,size(tpl1(:,:,:),1)
            tpl1(i,j,KS) = temp(i,j,KE)
         end do
      end do
     do kk = KS+1,KE
        do j = 1,size(tpl1(:,:,:),2)
           do i = 1,size(tpl1(:,:,:),1)
              tpl1(i,j,kk) = tflux(i,j,kk)
           end do
        end do
     end do
     do j = 1,size(tpl1(:,:,:),2)
        do i = 1,size(tpl1(:,:,:),1)
           tpl1(i,j,KE+1) = 0.5E+00*(tflux(i,j,KE+1) +   &
                                     temp(i,j,KE))
        end do
     end do
  do kk = KS+1,KE
     do j = 1,size(tpl2(:,:,:),2)
        do i = 1,size(tpl2(:,:,:),1)
           tpl2(i,j,kk) = tflux(i,j,kk)
        end do
     end do
  end do
    do j = 1,size(tpl2(:,:,:),2)
       do i = 1,size(tpl2(:,:,:),1)
          tpl2(i,j,KE+1) = 0.5E+00*(tflux(i,j,KE) +    &
                                temp(i,j,KE))
       end do
    end do


!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ittmp(:,:,:)=ixoe2(:,:,:)
  do kk = KS,KE
     do j = 1,size(ixoe2(:,:,:),2)
        do i = 1,size(ixoe2(:,:,:),1)
           ixoe2(i,j,kk) = ittmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  ttmp(:,:,:)=dte2 (:,:,:)
  do kk = KS,KE
     do j = 1,size(dte2(:,:,:),2)
        do i = 1,size(dte2(:,:,:),1)
           dte2(i,j,kk) = ttmp(i,j,kk+KS+1-(KS))
        end do
     end do
  end do
  do j = 1,size(ixoe2(:,:,:),2)
     do i = 1,size(ixoe2(:,:,:),1)
        ixoe2(i,j,KE+1) = ixoe1(i,j,KE)
     end do
  end do
  do j = 1,size(dte2(:,:,:),2)
     do i = 1,size(dte2(:,:,:),1)
        dte2(i,j,KE+1) = dte1 (i,j,KE)
     end do
  end do




      do j = 1,size(ixsp(:,:,:),2)
         do i = 1,size(ixsp(:,:,:),1)
            ixsp(i,j,KE) = ixoe2(i,j,KE-1)
         end do
      end do
      do j = 1,size(ixsp(:,:,:),2)
         do i = 1,size(ixsp(:,:,:),1)
            ixsp(i,j,KE+1) = ixoe1(i,j,KE-1)
         end do
      end do
      do j = 1,size(dxsp(:,:,:),2)
         do i = 1,size(dxsp(:,:,:),1)
            dxsp(i,j,KE) = dte2(i,j,KE-1)
         end do
      end do
      do j = 1,size(dxsp(:,:,:),2)
         do i = 1,size(dxsp(:,:,:),1)
            dxsp(i,j,KE+1) = dte1(i,j,KE-1)
         end do
      end do

!--------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute nearby layer transmissivities for h2o.
!--------------------------------------------------------------------
   if (Lw_control%do_h2o) then
      call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
      do kk = KS,KE+1
         do j = 1,size(ylog(:,:,:),2)
            do i = 1,size(ylog(:,:,:),1)
               ylog(i,j,kk) = ALOG10(Optical%empl1(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
      iysp(:,:,KS:KE+1) = iysp(:,:,KS:KE+1) + 1
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd1, KS, KE+1)
   else
      call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
      iysp1(:,:,KS:KE+1) = 1
      dysp1(:,:,KS:KE+1) = 0.0
      call looktab (tab3, ixsp, iysp1, dxsp, dysp1, emd1, KS, KE+1)
   endif

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!------------------------------------------------------------------
      if (nbtrge > 0) then
        do m=1,NBTRGE
    if (Lw_control%do_h2o) then
          do kk = KS,KE+1
             do j = 1,size(ylog(:,:,:),2)
                do i = 1,size(ylog(:,:,:),1)
                   ylog(i,j,kk) = ALOG10(Optical%empl1f(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
          iysp(:,:,KS:KE+1) = iysp(:,:,KS:KE+1) + 1
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, &
                        emd1f(:,:,:,m), KS, KE+1, m)
    else
          iysp1(:,:,KS:KE+1) = 1
          dysp1(:,:,KS:KE+1) = 0.0
          call looktab (tab3a, ixsp, iysp1, dxsp, dysp1, &
                        emd1f(:,:,:,m), KS, KE+1, m)
    endif
        enddo
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
   if (Lw_control%do_h2o) then
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      do kk = KS+1,KE+1
         do j = 1,size(ylog(:,:,:),2)
            do i = 1,size(ylog(:,:,:),1)
               ylog(i,j,kk) = ALOG10(Optical%empl2(i,j,kk))
            end do
         end do
      end do
      call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
      iysp(:,:,KS+1:KE+1) = iysp(:,:,KS+1:KE+1) + 1
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd2, KS+1, KE+1)
   else
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      iysp1(:,:,KS+1:KE+1) = 1
      dysp1(:,:,KS+1:KE+1) = 0.0
      call looktab (tab3, ixsp, iysp1, dxsp, dysp1, emd2, KS+1, KE+1)
   endif

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (nbtrge > 0) then
        do m=1,NBTRGE
   if (Lw_control%do_h2o) then
          do kk = KS+1,KE+1
             do j = 1,size(ylog(:,:,:),2)
                do i = 1,size(ylog(:,:,:),1)
                   ylog(i,j,kk) = ALOG10(Optical%empl2f(i,j,kk,m))
                end do
             end do
          end do
          call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
          iysp(:,:,KS+1:KE+1) = iysp(:,:,KS+1:KE+1) + 1
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, &
                        emd2f(:,:,:,m), KS+1, KE+1, m)
    else
          iysp1(:,:,KS+1:KE+1) = 1
          dysp1(:,:,KS+1:KE+1) = 0.0
          call looktab (tab3a, ixsp, iysp1, dxsp, dysp1, &
                        emd2f(:,:,:,m), KS+1, KE+1, m)
    endif
        enddo
      endif

!---------------------------------------------------------------------- 
!     compute nearby layer and special-case transmissivities for
!     emissivity using methods for h2o given in reference (4).
!-------------------------------------------------------------------- 
      do kk = KS+1,KE
         do j = 1,size(emisdg(:,:,:),2)
            do i = 1,size(emisdg(:,:,:),1)
               emisdg(i,j,kk) = emd2(i,j,kk) + emd1(i,j,kk)
            end do
         end do
      end do
      do j = 1,size(emisdg(:,:,:),2)
         do i = 1,size(emisdg(:,:,:),1)
            emisdg(i,j,KE+1) = 2.0E+00*emd1(i,j,KE+1)
         end do
      end do

     if (nbtrge > 0) then
       do m=1,NBTRGE
         do kk = KS+1,KE
            do j = 1,size(emisdgf(:,:,:,:),2)
               do i = 1,size(emisdgf(:,:,:,:),1)
                  emisdgf(i,j,kk,m) = &
                            emd2f(i,j,kk,m) + emd1f(i,j,kk,m)
               end do
            end do
         end do
         do j = 1,size(emisdgf(:,:,:,:),2)
            do i = 1,size(emisdgf(:,:,:,:),1)
               emisdgf(i,j,KE+1,m) = 2.0E+00*emd1f(i,j,KE+1,m)
            end do
         end do
       enddo
     endif

!--------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
  if (nbtrge > 0) then
    if (Lw_control%do_ch4 .or. Lw_control%do_n2o) then
      do kk = KS+1,KE+1
         do j = 1,size(emisdgf(:,:,:,:),2)
            do i = 1,size(emisdgf(:,:,:,:),1)
               emisdgf(i,j,kk,1) = emisdgf(i,j,kk,1) *   &
                                 tch4n2oe(i,j,kk,1)
            end do
         end do
      end do
     endif 

!--------------------------------------------------------------------
!     add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_indx8_part (8, Optical, tcfc8, KE)
        do kk = KS+1,KE+1
           do j = 1,size(emisdgf(:,:,:,:),2)
              do i = 1,size(emisdgf(:,:,:,:),1)
                 emisdgf(i,j,kk,1) = emisdgf(i,j,kk,1) *   &
                                   tcfc8(i,j,kk)
              end do
           end do
        end do
      endif
    endif ! (nbtrge > 0)

!------------------------------------------------------------------



end subroutine enear



!####################################################################
! <SUBROUTINE NAME="co2_source_calc">
!  <OVERVIEW>
!   Subroutine to calculate CO2 source function
!  </OVERVIEW>
!  <TEMPLATE>
!   call co2_source_calc (Atmos_input, Rad_gases, sorc,  Gas_tf, &
!                         source_band, dsrcdp_band)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties
!  </IN>
!  <OUT NAME="sorc" TYPE="real">
!   CO2 source function results
!  </OUT>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission functions
!  </IN>
!  <OUT NAME="source_band" TYPE="real">
!   CO2 source function bands
!  </OUT>
!  <OUT NAME="dsrcdp_band" TYPE="real">
!   Difference of source function between nearby thermal layers
!  </OUT>
! </SUBROUTINE>
!
subroutine co2_source_calc (Atmos_input, Rad_gases, sorc,  Gas_tf, &
                            source_band, dsrcdp_band)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

type(atmos_input_type),     intent(in)   ::  Atmos_input
type(radiative_gases_type), intent(in)   ::  Rad_gases
real, dimension(:,:,:,:),   intent(out)  ::  sorc               
type(gas_tf_type),          intent(in)   ::  Gas_tf
real, dimension(:,:,:,:),   intent(out)  ::  source_band, dsrcdp_band  

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     Atmos_input
!     Rad_gases
!     Gas_tf
!
!  intent(out) variables:
!
!     sorc
!     source_band
!     dsrcdp_band
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::   &
                                              dte1, press, pflux, temp

      integer, dimension (size(Atmos_input%temp,1),   &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)) ::          ixoe1

      integer            ::   i,j,kk
      integer            ::   n, ioffset, m
      integer            :: nbly                    
      real               :: rrvco2

!---------------------------------------------------------------------
!  local variables:
!
!     dte1
!     press
!     pflux
!     temp
!     ixoe1
!     n
!     ioffset
!     m
!     nbly
!     rrvco2
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      ioffset = Lw_parameters%offset
      nbly = 16+ioffset

!--------------------------------------------------------------------
 !  convert press and pflux to cgs.
        do kk = 1,size(press(:,:,:),3)
           do j = 1,size(press(:,:,:),2)
              do i = 1,size(press(:,:,:),1)
                 press(i,j,kk) = 10.0*Atmos_input%press(i,j,kk)
              end do
           end do
        end do
        do kk = 1,size(pflux(:,:,:),3)
           do j = 1,size(pflux(:,:,:),2)
              do i = 1,size(pflux(:,:,:),1)
                 pflux(i,j,kk) = 10.0*Atmos_input%pflux(i,j,kk)
              end do
           end do
        end do
     do kk = 1,size(temp(:,:,:),3)
        do j = 1,size(temp(:,:,:),2)
           do i = 1,size(temp(:,:,:),1)
              temp(i,j,kk) = Atmos_input%temp(i,j,kk)
           end do
        end do
     end do
      rrvco2 = Rad_gases%rrvco2

!----------------------------------------------------------------------
!     compute source function for frequency bands (9+ioffset to NBLY-1) 
!     at layer temperatures using table lookup.
!----------------------------------------------------------------------
      call locate_in_table(temp_1, temp, dte1, ixoe1, KS, Ke+1)
      do n=9+ioffset,NBLY-1
        call looktab (tabsr, ixoe1, dte1,   & 
                      sorc(:,:,:,n-ioffset-8), KS, ke+1, n)
      enddo

!-----------------------------------------------------------------------
!     compute the nlte source function for co2.
!-----------------------------------------------------------------------
      if (do_nlte) then
        call nlte (pflux, press, rrvco2, sorc, Gas_tf)
      endif

!----------------------------------------------------------------------
!    define "source function" appropriate for emissivity calculations
!    (temp**4), source functions for selected ranges including more 
!    than 1 frequency band (sorc15 for 15 um co2 band)
!    and differences in source functions (deltab) over
!    pressure layers.
!  
!    note: the values of sorc, sorc15, sorcwin, and derivatives 
!    depend on the no. of freq. bands!
!-----------------------------------------------------------------------
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,1) = Atmos_input%temp (i,j,kk+KS-(1))**4
            end do
         end do
      end do
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,2) = sorc(i,j,kk, 1) + &
                             sorc(i,j,kk, 2) + &
                             sorc(i,j,kk, 3 )
            end do
         end do
      end do
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,3) = sorc(i,j,kk,4 )
            end do
         end do
      end do
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,4) = sorc(i,j,kk,5 )
            end do
         end do
      end do
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,5) = sorc(i,j,kk, 6 )
            end do
         end do
      end do
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,6) = sorc(i,j,kk,7 )
            end do
         end do
      end do
      do m=1,NBTRGE
      do kk = 1,size(source_band(:,:,:,:),3)
         do j = 1,size(source_band(:,:,:,:),2)
            do i = 1,size(source_band(:,:,:,:),1)
               source_band(i,j,kk,6+m) = source_band(i,j,kk,1)
            end do
         end do
      end do
      end do

      do n=1, 6+NBTRGE       
       do kk = KS+1,KE+1
          do j = 1,size(dsrcdp_band(:,:,:,:),2)
             do i = 1,size(dsrcdp_band(:,:,:,:),1)
                dsrcdp_band(i,j,kk,n) = source_band(i,j,kk,n) - &
                               source_band(i,j,kk+KS-(KS+1),n)
             end do
          end do
       end do
      end do

!-------------------------------------------------------------------


end subroutine co2_source_calc




!#####################################################################
! <SUBROUTINE NAME="nlte">
!  <OVERVIEW>
!   nlte is the present formulation of an nlte calculation of the 
!     source function in the 15 um region (two bands).
!  </OVERVIEW>
!  <DESCRIPTION>
!   nlte is the present formulation of an nlte calculation of the 
!     source function in the 15 um region (two bands).
!
!     the essential theory is:
!
!           phi = C*j
!             j = b + E*phi
!
!     where
!             C = Curtis matrix
!              E = NLTE contribution (diagonal matrix)
!           phi = heating rate vector
!             b = LTE source function vector
!             j = NLTE source function vector
!
!             j = b (by assumption) for pressure layers > ixnltr
!             j = b (by assumption) for pressure layers > ixprnlte
!      E is obtained using a formulation devised by Fels (denoted
!      Ri in his notes).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call nlte (pflux, press, rrvco2, sorc, Gas_tf)
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   pressure values at flux levels.
!  </IN>
!  <IN NAME="press" TYPE="real">
!   pressure cordinates
!  </IN>
!  <IN NAME="rrvco2" TYPE="real">
!   CO2 volumn mixing ratio
!  </IN>
!  <INOUT NAME="sorc" TYPE="real">
!   CO2 source function to be calculated
!  </INOUT>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission function 
!  </IN>
! </SUBROUTINE>
!
subroutine nlte (pflux, press, rrvco2, sorc, Gas_tf)

!-----------------------------------------------------------------------
!     nlte is the present formulation of an nlte calculation of the 
!     source function in the 15 um region (two bands).
!
!     the essential theory is:
!
!           phi = C*j
!             j = b + E*phi
!
!     where
!             C = Curtis matrix
!              E = NLTE contribution (diagonal matrix)
!           phi = heating rate vector
!             b = LTE source function vector
!             j = NLTE source function vector
!
!             j = b (by assumption) for pressure layers > ixnltr
!             j = b (by assumption) for pressure layers > ixprnlte
!      E is obtained using a formulation devised by Fels (denoted
!      Ri in his notes).
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)    ::  pflux, press
real,                      intent(in)    ::  rrvco2
real, dimension (:,:,:,:), intent(inout) ::  sorc               
type(gas_tf_type),         intent(in)    :: Gas_tf

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     pflux
!     press
!     rrvco2
!     Gas_tf
!
!  intent(inout) variables:
!
!     sorc
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(press,1), size(press,2), &
                       ixprnlte) ::  &
                                ag, az, bdenom, cdiag, &
                                                tcoll, phifx, phivar

      real, dimension (size(press,1), size(press,2), &
                       ixprnlte, NBCO215) ::  &
                                    fnlte
      real, dimension (size(press,1), size(press,2), &
                       size(press,3)-1, ixprnlte ) ::  &
                                     cmtrx

      real                                   :: degen = 0.5
      integer                                :: i,j,kk
      integer                                :: n, k, inb, kp, ioffset

!---------------------------------------------------------------------
!  local variables:
!
!     ag
!     az
!     bdenom
!     cdiag
!     tcoll
!     phifx     fixed portion of PHI (contributions from
!               layers > ixnltr, where j(k) = b(k))
!               layers > ixprnlte, where j(k) = b(k))
!     phivar    varying portion of PHI (contributions
!               from layers <= ixprnlte).
!               from layers <= ixnltr).
!     fnlte     NLTE contribution: (E in above notes)
!     cmtrx
!     degen     degeneracy factor (= 0.5)
!     n
!     k
!     inb
!     kp
!     ioffset
!
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---------------------------------------------------------------------
      ioffset =  Lw_parameters%offset

!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute curtis matrix for both frequency bands.
!-----------------------------------------------------------------------
      call co2curt (pflux, cmtrx, Gas_tf)

      do k=KS,ixprnlte
        do j = 1,size(cdiag(:,:,:),2)
           do i = 1,size(cdiag(:,:,:),1)
              cdiag(i,j,k) = cmtrx(i,j,k,k)
           end do
        end do
      end do

!-----------------------------------------------------------------------
!   collisional relaxation time (see fels notes for "tcoll")
!-----------------------------------------------------------------------
      do k=KS,ixprnlte
        do j = 1,size(tcoll(:,:,:),2)
           do i = 1,size(tcoll(:,:,:),1)
              tcoll(i,j,k) = degen*1.5E-05*press(i,j,KE+1)/   &
                       (seconds_per_day*press(i,j,k)) 
           end do
        end do
      end do

!-----------------------------------------------------------------------
!   compute NLTE contribution for each band at each pressure level
!   <= ixprnlte. fnlte = zero by assumption at other levels.
!-----------------------------------------------------------------------
      do n=1,NBCO215
        do kk = KS,ixprnlte
           do j = 1,size(fnlte(:,:,:,:),2)
              do i = 1,size(fnlte(:,:,:,:),1)
                 fnlte(i,j,kk,n) = 3.5E+00*tcoll(i,j,kk)*  &
                                    c1b7(n)/(rrvco2*c2b7(n)) 
              end do
           end do
        end do
      enddo

!-----------------------------------------------------------------------
!     begin computations for (NBCO215) bands in 15um range.
!-----------------------------------------------------------------------
      do inb = 1,NBCO215
        do kk = KS,ixprnlte
           do j = 1,size(bdenom(:,:,:),2)
              do i = 1,size(bdenom(:,:,:),1)
                 bdenom(i,j,kk) = 1.0E+00/   &
              (1.0E+00 - fnlte(i,j,kk,inb)*   &
                        cdiag(i,j,kk))
              end do
           end do
        end do
        do kk = KS,ixprnlte
           do j = 1,size(phifx(:,:,:),2)
              do i = 1,size(phifx(:,:,:),1)
                 phifx(i,j,kk) = 0.0E+00
              end do
           end do
        end do
        do k=KS,ixprnlte
          do kp=ixprnlte+1,KE
            do j = 1,size(phifx(:,:,:),2)
               do i = 1,size(phifx(:,:,:),1)
                  phifx(i,j,k) = phifx(i,j,k) +   &
                           cmtrx(i,j,kp,k)*sorc(i,j,kp,inb           )
               end do
            end do
          end do
        end do
        do kk = KS,ixprnlte
           do j = 1,size(az(:,:,:),2)
              do i = 1,size(az(:,:,:),1)
                 az(i,j,kk) = sorc (i,j,kk,inb           ) +  &
                     fnlte(i,j,kk,inb)*phifx(i,j,kk)
              end do
           end do
        end do

!----------------------------------------------------------------------
!     first iteration. (J(k) = B(k)) as initial guess)
!-----------------------------------------------------------------------
        do kk = KS,ixprnlte
           do j = 1,size(phivar(:,:,:),2)
              do i = 1,size(phivar(:,:,:),1)
                 phivar(i,j,kk) = 0.0E+00
              end do
           end do
        end do
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            do j = 1,size(phivar(:,:,:),2)
               do i = 1,size(phivar(:,:,:),1)
                  phivar(i,j,k) = phivar(i,j,k) +   &
                            cmtrx(i,j,kp,k)*sorc(i,j,kp,inb           )
               end do
            end do
          end do
        end do
        do kk = KS,ixprnlte
           do j = 1,size(ag(:,:,:),2)
              do i = 1,size(ag(:,:,:),1)
                 ag(i,j,kk) = fnlte(i,j,kk,inb)*   &
                                (phivar(i,j,kk) -   &
                                 cdiag(i,j,kk)*  &
                                 sorc(i,j,kk,inb           ))
              end do
           end do
        end do

        do kk = KS,ixprnlte
           do j = 1,size(sorc(:,:,:,:),2)
              do i = 1,size(sorc(:,:,:,:),1)
                 sorc(i,j,kk,inb           ) = bdenom(i,j,kk)*&
                                               (az(i,j,kk) + &
                                                ag(i,j,kk)) 
              end do
           end do
        end do

!-----------------------------------------------------------------------
!     second iteration.  (J(k) = result of first iteration as guess)
!-----------------------------------------------------------------------
        do kk = KS,ixprnlte
           do j = 1,size(phivar(:,:,:),2)
              do i = 1,size(phivar(:,:,:),1)
                 phivar(i,j,kk) = 0.0E+00
              end do
           end do
        end do
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            do j = 1,size(phivar(:,:,:),2)
               do i = 1,size(phivar(:,:,:),1)
                  phivar(i,j,k) = phivar(i,j,k) +    &
                            cmtrx(i,j,kp,k)*sorc(i,j,kp,inb           )
               end do
            end do
          end do
        end do
        do kk = KS,ixprnlte
           do j = 1,size(ag(:,:,:),2)
              do i = 1,size(ag(:,:,:),1)
                 ag(i,j,kk) = fnlte(i,j,kk,inb)*   &
                        (phivar(i,j,kk) -   &
            cdiag(i,j,kk)*sorc(i,j,kk,inb           ))
              end do
           end do
        end do

        do kk = KS,ixprnlte
           do j = 1,size(sorc(:,:,:,:),2)
              do i = 1,size(sorc(:,:,:,:),1)
                 sorc(i,j,kk,inb           ) = bdenom(i,j,kk)*&
                                               (az(i,j,kk) +  &
                                                ag(i,j,kk)) 
              end do
           end do
        end do
      enddo

!-----------------------------------------------------------------------


end subroutine nlte



!#####################################################################
! <SUBROUTINE NAME="co2curt">
!  <OVERVIEW>
!   co2curt computes Curtis matrix elements derived from co2
!     transmission functions.
!  </OVERVIEW>
!  <TEMPLATE>
!   call co2curt (pflux, cmtrx, Gas_tf)
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   pressure values at flux levels
!  </IN>
!  <OUT NAME="cmtrx" TYPE="real">
!   Curtis matrix elements
!  </OUT>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   gas transmission function
!  </IN>
! </SUBROUTINE>
!
subroutine co2curt (pflux, cmtrx, Gas_tf)

!----------------------------------------------------------------------
!     co2curt computes Curtis matrix elements derived from co2
!     transmission functions.
!     functions.
!
!     author: m. d. schwarzkopf
!
!     revised: 8/18/94
!
!     certified:  radiation version 1.0
!
!---------------------------------------------------------------------
real, dimension(:,:, :),   intent(in)  :: pflux                  
real, dimension(:,:, :,:), intent(out) :: cmtrx                  
type(gas_tf_type),         intent(in)  :: Gas_tf

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     pflux
!     Gas_tf
!
!  intent(out) variables:
!
!     cmtrx    cutris matrix.
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! local variables:

      real, dimension (size(pflux,1),  &
                       size(pflux,2), &
                       size(pflux,3)-1) ::            pdfinv

      real, dimension (size(pflux,1),   &
                       size(pflux,2), &
                       size(pflux,3)) ::              co2row, co2rowp

     integer   :: k, krow, kp
     integer   :: i, j, kk

!---------------------------------------------------------------------
! local variables:
!
!     pdfinv
!     co2row
!     co2rowp
!     k
!     krow
!     kp
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute co2 transmission functions.
!-----------------------------------------------------------------------
      do kk = KS,KE+1
         do j = 1,size(co2row(:,:,:),2)
            do i = 1,size(co2row(:,:,:),1)
               co2row(i,j,kk) = 1.0E+00
            end do
         end do
      end do
      do kk = KS,KE+1
         do j = 1,size(co2rowp(:,:,:),2)
            do i = 1,size(co2rowp(:,:,:),1)
               co2rowp(i,j,kk) = 1.0E+00
            end do
         end do
      end do

!-----------------------------------------------------------------------
!    compute curtis matrix for rows from KS to ixprnlte
!-----------------------------------------------------------------------
      do k = KS,ixprnlte
        krow = k
        do j = 1,size(pdfinv(:,:,:),2)
           do i = 1,size(pdfinv(:,:,:),1)
              pdfinv(i,j,k) = 1.0/(pflux(i,j,k+1) - pflux(i,j,k))
           end do
        end do

        call transcol ( KS, krow, KS, KE+1, co2row, Gas_tf)
        call transcol ( KS, krow+1, KS, KE+1, co2rowp, Gas_tf)
        do kp=KS,KE-1 
          do j = 1,size(cmtrx(:,:,:,:),2)
             do i = 1,size(cmtrx(:,:,:,:),1)
                cmtrx(i,j,kp,k) = radcon*pdfinv(i,j,k)*   &
                            (co2rowp(i,j,kp) - co2rowp(i,j,kp+1) -  &
                             co2row(i,j,kp) + co2row(i,j,kp+1)) 
             end do
          end do
        end do

        do j = 1,size(cmtrx(:,:,:,:),2)
           do i = 1,size(cmtrx(:,:,:,:),1)
              cmtrx(i,j,KE,k) = radcon*pdfinv(i,j,k)*   &
                          (co2rowp(i,j,KE) - co2row(i,j,KE)) 
           end do
        end do
      enddo

!--------------------------------------------------------------------


end subroutine co2curt




!####################################################################

! <SUBROUTINE NAME="co2_time_vary">
!  <OVERVIEW>
!   Calculate CO2 absorption coefficient based on its volume
!   mixing ratio using precomputed lbl tables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculate CO2 absorption coefficient based on its volume
!   mixing ratio using precomputed lbl tables
!  </DESCRIPTION>
!  <TEMPLATE>
!   call co2_time_vary ( rrvco2 )
!  </TEMPLATE>
!  <IN NAME="rrvco2" TYPE="real">
!   CO2 volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine co2_time_vary ( rrvco2, phase )

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, intent(in   )    ::  rrvco2
integer, intent(in), optional  ::  phase

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     rrvco2
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real    ::    co2_vmr   !
      integer ::    lphase

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
        if ( present(phase) ) then
          lphase = phase
        else
          lphase = -1
        end if
        co2_vmr = rrvco2*1.0E+06

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
        call co2_lblinterp  (co2_vmr, lphase)

!--------------------------------------------------------------------


end subroutine co2_time_vary



!####################################################################
! <SUBROUTINE NAME="ch4_time_vary">
!  <OVERVIEW>
!   Calculate CH4 and N2O absorption coefficients from their
!   mixing ratios using precomputed lbl tables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculate CH4 and N2O absorption coefficients from their
!   mixing ratios using precomputed lbl tables
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ch4_n2o_time_vary (rrvch4, rrvn2o)
!  </TEMPLATE>
!  <IN NAME="rrvch4" TYPE="real">
!   ch4 volume mixing ratio
!  </IN>
!  <IN NAME="rrvn2o" TYPE="real">
!   n2o volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine ch4_time_vary (rrvch4,phase)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, intent(in) :: rrvch4         
integer, intent(in), optional :: phase

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      rrvch4
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
 
     real   ::  ch4_vmr !
     integer :: lphase

!---------------------------------------------------------------------
!  the ch4 volume mixing ratio is set to the initial value (rch4) and 
!  the mass mixing ratio is defined on the first access of this
!  routine. then the lbl transmission function is calculated. after 
!  first access, this routine does nothing. 
!--------------------------------------------------------------------
        if ( present(phase) ) then
          lphase = phase
        else
          lphase = -1
        end if
        ch4_vmr = rrvch4*1.0E+09
        call Ch4_lblinterp  (ch4_vmr, lphase)

!----------------------------------------------------------------------

end subroutine ch4_time_vary


!####################################################################
! <SUBROUTINE NAME="n2o_time_vary">
!  <OVERVIEW>
!   Calculate CH4 and N2O absorption coefficients from their
!   mixing ratios using precomputed lbl tables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculate CH4 and N2O absorption coefficients from their
!   mixing ratios using precomputed lbl tables
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ch4_n2o_time_vary (rrvch4, rrvn2o)
!  </TEMPLATE>
!  <IN NAME="rrvch4" TYPE="real">
!   ch4 volume mixing ratio
!  </IN>
!  <IN NAME="rrvn2o" TYPE="real">
!   n2o volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine n2o_time_vary (rrvn2o,phase)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, intent(in) :: rrvn2o            
integer, intent(in), optional :: phase

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      rrvch4
!      rrvn2o
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
 
     real   ::  n2o_vmr !
     integer :: lphase

!---------------------------------------------------------------------
!  the n2o volume mixing ratio is set to initial value (rn2o) and the 
!  mass mixing ratio is defined on the first access of this routine. 
!  routines are called to calculate the lbl transmission functions for 
!  n2o. after first access, this routine does nothing. 
!--------------------------------------------------------------------
        if ( present(phase) ) then
          lphase = phase
        else
          lphase = -1
        end if
        n2o_vmr = rrvn2o*1.0E+09
        call N2o_lblinterp (n2o_vmr, lphase)

!----------------------------------------------------------------------

end subroutine n2o_time_vary



!####################################################################




                  end module sealw99_mod
