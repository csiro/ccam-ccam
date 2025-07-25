                  module lw_gases_stdtf_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Module that computes longwave gas transmission functions
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!

!  shared modules:
     
!use fms_mod,              only: open_namelist_file, fms_init, &
!                                mpp_pe, mpp_root_pe, stdlog, &
!                                file_exist, write_version_number, &
!                                check_nml_error, error_mesg, &
!                                FATAL, NOTE, WARNING, close_file, &
!                                open_direct_file, mpp_error
!use fms_io_mod,           only: read_data
!  shared radiation package modules:

use rad_utilities_mod,    only: rad_utilities_init, Lw_parameters, &
                                Lw_control, optical_path_type

!  radiation package modules

use gas_tf_mod,           only: gas_tf_init,  &
                                put_co2_stdtf_for_gas_tf, &
                                put_co2_nbltf_for_gas_tf, &
                                put_ch4_stdtf_for_gas_tf, &
                                put_n2o_stdtf_for_gas_tf, &
                                get_control_gas_tf, &
                                process_co2_input_file, &
                                process_ch4_input_file, &
                                process_n2o_input_file

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    lw_gases_stdf_mod computes line-by-line transmission 
!    functions for co2, ch4 and n2o for a usstd temperature 
!    profile and (if needed) that profile +/- 25 degrees, for 
!    the vertical layer structure of the atmospheric model 
!    with surface pressures of 1013.25 hPa and 0.8*1013.25 hPa.
!    options are taken from namelist(s).
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: lw_gases_stdtf.F90,v 17.0.4.1 2010/08/30 20:33:32 wfc Exp $'
character(len=128)  :: tagname =  '$Name: testing $'


!---------------------------------------------------------------------
!-------  interfaces --------

public         &
        lw_gases_stdtf_init, lw_gases_stdtf_time_vary,      &
        ch4_lblinterp, co2_lblinterp, n2o_lblinterp, &
        lw_gases_stdtf_dealloc,  lw_gases_stdtf_end, &
        cfc_exact, cfc_exact_part, cfc_indx8, cfc_indx8_part, &
        cfc_overod, cfc_overod_part

private        &
        std_lblpressures, approx_fn, approx_fn_std, &
        gasins, gasint, coeint, intcoef_1d, intcoef_2d, &
        intcoef_2d_std, interp_error, interp_error_r, &
        pathv1, rctrns, read_lbltfs, allocate_interp_arrays, &
        deallocate_interp_arrays


!---------------------------------------------------------------------
!-------- namelist  ---------

logical, save    :: do_co2_bug = .false.
logical, save    :: do_coeintdiag = .false.
integer, save    :: NSTDCO2LVLS = 496 ! # of levels at which lbl tfs exist
integer, save    :: llco2 = 1  ! #no iterations in coeint for co2
integer, save    :: llch4 = 1  ! #no iterations in coeint for ch4
integer, save    :: lln2o = 1  ! #no iterations in coeint for n2o

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


real,    dimension (:,:),   allocatable, save  ::                      &   
                                            pressint_hiv_std_pt1,   & 
                                            pressint_lov_std_pt1,   & 
                                            pressint_hiv_std_pt2,   &
                                            pressint_lov_std_pt2
integer, dimension (:,:),   allocatable, save  ::                      &   
                                            indx_pressint_hiv_std_pt1, &
                                            indx_pressint_lov_std_pt1, &
                                            indx_pressint_hiv_std_pt2, &
                                            indx_pressint_lov_std_pt2

!--------------------------------------------------------------------
!    xa, ca, dop_core, uexp, sexp are coefficients for
!    the approximation function (Eq. (4) in Ref. (2)) used trns_std_hi_nf
!    the co2 interpolation algorithm. the nomenclature is:
!
!      this code           Ref. (2)
!      ---------           --------
!       xa                  X (see Eq. A1) 
!       ca                  C (see Eq. A1)
!       uexp                delta (see Eq. A6b)
!       sexp                gamma (see Eq. A6c)
!       dop_core            core (see Eq. A6a)
!----------------------------------------------------------------------
real,    dimension (:),     allocatable, save  :: xa, ca, uexp, sexp, &
                                            press_lo, press_hi
real,    dimension(:,:),   allocatable, save  :: trns_std_hi, &
                                                  trns_std_lo

!---------------------------------------------------------------------
!   pa          = pressure levels where line-by-line co2 transmission
!                 functions have been calculated
!---------------------------------------------------------------------
real, dimension(:), allocatable, save   :: pa

!----------------------------------------------------------------------
!    ch4 data
!----------------------------------------------------------------------
integer, parameter                              ::  number_std_ch4_vmrs = 9
real,    dimension(number_std_ch4_vmrs), save   ::   ch4_std_vmr
data ch4_std_vmr / 0., 300., 700., 1250., 1750., 2250., 2800., 4000., 6000. /

integer, parameter                              ::  nfreq_bands_sea_ch4 = 1

logical, dimension(nfreq_bands_sea_ch4), save   ::  do_lyrcalc_ch4_nf, &
                                                    do_lvlcalc_ch4_nf, &
                                                    do_lvlctscalc_ch4_nf
data   do_lyrcalc_ch4_nf    /  .true.  /
data   do_lvlcalc_ch4_nf    /  .true.  /
data   do_lvlctscalc_ch4_nf /  .false. /

integer, dimension(nfreq_bands_sea_ch4), save   ::  ntbnd_ch4
data   ntbnd_ch4      /  3  /

!----------------------------------------------------------------------
!    n2o data
!----------------------------------------------------------------------
integer, parameter                              ::  number_std_n2o_vmrs = 9
real,    dimension(number_std_n2o_vmrs), save   ::  n2o_std_vmr
data n2o_std_vmr / 0., 180., 275., 310., 340., 375., 500., 600., 800. /

integer, parameter                              ::  nfreq_bands_sea_n2o = 3
logical, dimension(nfreq_bands_sea_n2o), save   ::  do_lyrcalc_n2o_nf, &
                                                    do_lvlcalc_n2o_nf, &
                                                    do_lvlctscalc_n2o_nf
data do_lyrcalc_n2o_nf    / .true., .true., .true./
data do_lvlcalc_n2o_nf    / .true., .true., .true./
data do_lvlctscalc_n2o_nf / .false., .false., .false./

integer, dimension(nfreq_bands_sea_n2o), save   ::  ntbnd_n2o
data ntbnd_n2o /  3, 3, 3/

!----------------------------------------------------------------------
!    co2 data
!----------------------------------------------------------------------
integer, parameter                        ::  number_std_co2_vmrs = 14
real,    dimension(number_std_co2_vmrs)   ::  co2_std_vmr
data co2_std_vmr / 0., 165.0, 300.0, 330.0, 348.0, 356.0, 360.0,  &
                   600.0, 660.0, 1320.0, 1600.0, 2000.0, 5000.0, 10000.0 /

integer, parameter                              ::  nfreq_bands_sea_co2 = 8
logical, dimension(nfreq_bands_sea_co2), save   ::  do_lyrcalc_co2_nf, &
                                                    do_lvlcalc_co2_nf, &
                                                    do_lvlctscalc_co2_nf
data do_lyrcalc_co2_nf    / .true., .false., .false., .false., .true., .true., .true., .true./
data do_lvlcalc_co2_nf    / .true., .true., .true., .true., .true., .true., .true., .true./
data do_lvlctscalc_co2_nf / .false., .true., .true., .true., .false., .true., .true., .true./

integer, dimension(nfreq_bands_sea_co2), save   ::  ntbnd_co2
data ntbnd_co2 / 3, 3, 3, 3, 1, 3, 3, 3/

integer, parameter :: nfreq_bands_sea_all = nfreq_bands_sea_ch4 &
                                           +nfreq_bands_sea_n2o &
                                           +nfreq_bands_sea_co2

real,  dimension (:,:,:), allocatable, save   :: dgasdt8_lvl, dgasdt10_lvl, &
                                                 d2gast8_lvl, d2gast10_lvl, &
                                                 gasp10_lvl, gasp8_lvl,   &  
                                                 dgasdt8_lyr, dgasdt10_lyr, &
                                                 d2gast8_lyr, d2gast10_lyr, &
                                                 gasp10_lyr, gasp8_lyr
real,  dimension (:,:), allocatable, save :: dgasdt8_lvlcts, dgasdt10_lvlcts, &
                                             d2gast8_lvlcts, d2gast10_lvlcts, &
                                             gasp10_lvlcts, gasp8_lvlcts

real,  dimension (:,:), allocatable, save   :: trns_interp_lyr_ps, &
                                               trns_interp_lyr_ps8, &
                                               trns_interp_lvl_ps, &
                                               trns_interp_lvl_ps8
real,  dimension (:,:,:), allocatable, save :: trns_interp_lyr_ps_nf, &
                                               trns_interp_lyr_ps8_nf, &
                                               trns_interp_lvl_ps_nf, &
                                               trns_interp_lvl_ps8_nf
 
 
real, dimension(:), allocatable, save  :: plm, plm8, pd, pd8

!!$integer             :: k, kp, nf, nt
integer, save       :: ndimkp, ndimk, nlev
real, parameter     :: dop_core0 = 25.0
real, save          :: dop_core 
logical, save       :: do_calcstdco2tfs
logical, save       :: do_calcstdch4tfs
logical, save       :: do_calcstdn2otfs

!--------------------------------------------------------------------
!       NBLWCFC =  number of frequency bands with cfc band strengths
!                  included. The bands have the same frequency ranges
!                  as those used for h2o calculations
!--------------------------------------------------------------------
integer, parameter :: NBLWCFC = 8


!--------------------------------------------------------------------
!   data for averaged f11 band strength
!--------------------------------------------------------------------
real, dimension(nblwcfc), save :: strf11

data  strf11 /       &
         0.000000E+00,  0.000000E+00,  0.527655E+02,  0.297523E+04,  &
         0.134488E+03,  0.247279E+03,  0.710717E+03,  0.000000E+00/

!--------------------------------------------------------------------
!   data for averaged f12 band strength
!--------------------------------------------------------------------
real, dimension(nblwcfc), save :: strf12

data strf12 /       &
         0.552499E+01,  0.136436E+03,  0.243867E+02,  0.612532E+03, &
         0.252378E+04,  0.438226E+02,  0.274950E+04,  0.000000E+00/

!--------------------------------------------------------------------
!   data for averaged f113 band strength
!--------------------------------------------------------------------
real, dimension(nblwcfc), save :: strf113

data strf113 /     &
         0.627223E+01,  0.690936E+02,  0.506764E+02,  0.122039E+04,  &
         0.808762E+03,  0.742843E+03,  0.109485E+04,  0.194768E+03/

!--------------------------------------------------------------------
!   data for averaged f22 band strength
!--------------------------------------------------------------------
real, dimension(nblwcfc), save :: strf22

data strf22 /    &
         0.301881E+02,  0.550826E+01,  0.397496E+03,  0.124802E+04,  &
         0.190285E+02,  0.460065E+02,  0.367359E+04,  0.508838E+03/

!--------------------------------------------------------------------
!   data for averaged f11 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf1115=0.219856E+02

!--------------------------------------------------------------------
!   data for averaged f12 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf1215=0.515665E+02

!--------------------------------------------------------------------
!   data for averaged f113 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf11315=0.430969E+02

!--------------------------------------------------------------------
!   data for averaged f22 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf2215=0.176035E+03

!--------------------------------------------------------------------
!   data for averaged f11 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf11ct=0.125631E+04

!--------------------------------------------------------------------
!   data for averaged f12 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf12ct=0.201821E+04

!--------------------------------------------------------------------
!   data for averaged f113 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf113ct=0.105362E+04

!--------------------------------------------------------------------
!   data for averaged f22 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf22ct=0.188775E+04

integer, save :: ksrad, kerad
logical, save   :: module_is_initialized = .false.

logical, save :: rctrns_first = .true.
logical, save :: intcoef_2d_std_first = .true.

!---------------------------------------------------------------------
!---------------------------------------------------------------------



                        contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="lw_gases_stdtf_init">
!  <OVERVIEW>
!   Subroutine to initialize longwave gas transmission function 
!   calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize longwave gas transmission function 
!   calculation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_gases_stdtf_init ( pref)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference level pressure array
!  </IN>
! </SUBROUTINE>
!
subroutine lw_gases_stdtf_init ( pref)

!-------------------------------------------------------------------
!
!-------------------------------------------------------------------

real,  dimension(:,:), intent(in) :: pref

!--------------------------------------------------------------------
!  intent(in)  variables:
!
!    pref
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer :: kmin, kmax, k

!---------------------------------------------------------------------
!  local variables:
!
!    unit
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
      call gas_tf_init (pref)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      allocate (pd   (size(pref,1)   ))
      allocate (plm  (size(pref,1)  ))
      allocate (pd8  (size(pref,1)  ))
      allocate (plm8 (size(pref,1 ) ))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kmin = 1
      kmax = size(pref,1) - 1
      pd (:) = pref(:,1)
      pd8(:) = pref(:,2)
      plm (kmin) = 0.
      plm8(kmin) = 0.
      do k=kmin+1,kmax
        plm (k) = 0.5*(pd (k-1) + pd (k))
        plm8(k) = 0.5*(pd8(k-1) + pd8(k))
      enddo
      plm (kmax+1) = pd (kmax+1)
      plm8(kmax+1) = pd8(kmax+1)
      pd = pd*1.0E-02
      pd8 = pd8*1.0E-02
      plm =plm*1.0E-02           
      plm8 = plm8*1.0E-02           

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ksrad = 1
      kerad = kmax

!--------------------------------------------------------------------
!    define the standard pressure levels for use in calculating the 
!    transmission functions.
!--------------------------------------------------------------------- 
      call std_lblpressures

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine lw_gases_stdtf_init



!#####################################################################
! <SUBROUTINE NAME="lw_gases_stdtf_time_vary">
!  <OVERVIEW>
!   Allocate transmission function memory tables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate transmission function memory tables
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_gases_stdtf_time_vary
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_gases_stdtf_time_vary

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!-------------------------------------------------------------------
!    determine if the tfs are to be calculated.
!-------------------------------------------------------------------
      call get_control_gas_tf (calc_co2=do_calcstdco2tfs, &
                               calc_n2o=do_calcstdn2otfs, &
                               calc_ch4=do_calcstdch4tfs)

!-------------------------------------------------------------------
!    define the number of levels being used in the calculation and
!    the dimension extents of the interpolation arrays:
!-------------------------------------------------------------------
      if (do_calcstdco2tfs .or. do_calcstdn2otfs .or.   &
          do_calcstdch4tfs) then
        nlev = KERAD - KSRAD + 1
        ndimkp = nlev + 1
        ndimk  = nlev + 1

!---------------------------------------------------------------------
!    allocate module variables.
!---------------------------------------------------------------------
        if ( .not.allocated(xa) ) then
          allocate (xa          (NSTDCO2LVLS) )
          allocate (ca          (NSTDCO2LVLS) )
          allocate (uexp        (NSTDCO2LVLS) )
          allocate (sexp        (NSTDCO2LVLS) )
          allocate (press_lo    (NSTDCO2LVLS) )
          allocate (press_hi    (NSTDCO2LVLS) )

          allocate (  dgasdt8_lvl (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  dgasdt10_lvl(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast8_lvl (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast10_lvl(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp10_lvl(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp8_lvl (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  dgasdt8_lvlcts (KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  dgasdt10_lvlcts(KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast8_lvlcts (KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast10_lvlcts(KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp10_lvlcts(KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp8_lvlcts (KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  dgasdt8_lyr (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  dgasdt10_lyr(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast8_lyr (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  d2gast10_lyr(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp10_lyr(KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
          allocate (  gasp8_lyr (KSRAD:KERAD+1,KSRAD:KERAD+1,nfreq_bands_sea_all) )
        end if
      endif

!------------------------------------------------------------------


end subroutine lw_gases_stdtf_time_vary



!###################################################################
! <SUBROUTINE NAME="ch4_lblinterp">
!  <OVERVIEW>
!   Subroutine to interpolate ch4 transmission function to user
!   specified pressure levels and ch4 concentration
!  </OVERVIEW>
!  <DESCRIPTION>
!   this routine is 1) a standalone program for a ch4 interpolation
!     to user-specified pressure levels and ch4 concentration;
!     2) an interface between a GCM and ch4 interpolation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ch4_lblinterp (ch4_vmr)
!  </TEMPLATE>
!  <IN NAME="ch4_vmr" TYPE="real">
!   ch4 volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine ch4_lblinterp (ch4_vmr, phase)

use cc_mpi
use newmpar_m, only : nproc

!-----------------------------------------------------------------
!    ch4_lblinterp is 
!    1) a standalone program for a ch4 interpolation
!       to user-specified pressure levels and ch4 concentration;
!    2) an interface between a GCM and ch4 interpolation
!
!    input files:
!
!         1)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (ch4_vmr).
!         2)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (ch4_vmr). may not be
!                  used, depending on value of (ch4_vmr).
!
!    output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these 
!                  files are created if ifdef (writeinterpch4) is
!                  on. otherwise, it is assumed the data is fed
!                 directly back to the parent model.
!-----------------------------------------------------------------

real,              intent(in)  :: ch4_vmr
integer,           intent(in)  :: phase

!-----------------------------------------------------------------
!  intent(in) variables:
!
!     ch4_vmr
!
!-----------------------------------------------------------------

        
!---------------------------------------------------------------------
!  local variables:

      logical                   ::  callrctrns_ch4
      logical                   ::  do_lyrcalc_ch4
      logical                   ::  do_lvlcalc_ch4
      logical                   ::  do_lvlctscalc_ch4
      real                      ::  ch4_std_lo, ch4_std_hi
      integer                   ::  n, nf, nt
      integer                   ::  nstd_ch4_lo, nstd_ch4_hi
      integer                   ::  origin, nf_offset
      integer, parameter        ::  offset = 0
      character(len=8)          ::  gas_type = 'ch4'
      real, dimension(:,:,:), allocatable, save :: trns_std_hi_nf, trns_std_lo_nf


!---------------------------------------------------------------------
!  local variables:
!
!     callrctrns_ch4
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
      
      allocate( trns_std_hi_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
      allocate( trns_std_lo_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
      
!---------------------------------------------------------------------
!    this routine does nothing unless the calculation of tfs is desired.
!---------------------------------------------------------------------
      if (do_calcstdch4tfs) then

        if ( phase==0 .or. phase==-1 ) then

          rctrns_first = .true.
          intcoef_2d_std_first = .true.
          
          allocate (trns_std_hi(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_std_lo(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_interp_lyr_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lyr_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          allocate (trns_interp_lvl_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lvl_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          
!--------------------------------------------------------------------
!    using the value of the ch4 volume mixing ratio (ch4_vmr) and
!    the available standard ch4 mixing ratios (with lbl 
!    transmission functions) obtain the two standard mixing ratios
!    which bracket (ch4_vmr). if (as in a fixed ch4 experiment) the
!    difference between (ch4_vmr) and the higher of the standard
!    mixing ratios is less than a tolerance (taken as 0.1 ppbv)
!    we will assume that that standard ch4 transmissivity applies
!    to (ch4_vmr) without interpolation. otherwise, interpolation
!    to (ch4_vmr) will be performed, in rctrns.F
!--------------------------------------------------------------------
          if (ch4_vmr .LT. ch4_std_vmr(1) .OR.               &
              ch4_vmr .GT. ch4_std_vmr(number_std_ch4_vmrs)) then
            write(6,*) "ch4 volume mixing ratio is out of range"
            stop
          endif

          if (ch4_vmr .EQ. ch4_std_vmr(1)) then
            ch4_std_lo = ch4_std_vmr(1)
            ch4_std_hi = ch4_std_vmr(1)
            nstd_ch4_lo = 1
            nstd_ch4_hi = 1
          else 
            do n=1,number_std_ch4_vmrs-1
              if (ch4_vmr .GT. ch4_std_vmr(n) .AND.            &
                  ch4_vmr .LE. ch4_std_vmr(n+1)) then
                ch4_std_lo = ch4_std_vmr(n)
                ch4_std_hi = ch4_std_vmr(n+1)
                nstd_ch4_lo = n
                nstd_ch4_hi = n+1
                exit
              endif
            enddo
          endif

!--------------------------------------------------------------------
!    ch4_std_lo, nstd_ch4_lo have arbitrary definitions, since they
!    will not be used, as callrctrns will be false in this case.
!--------------------------------------------------------------------
          if (ABS(ch4_vmr - ch4_std_hi) .LE. 1.0e-1) then
            callrctrns_ch4 = .false.
          else
            callrctrns_ch4 = .true.
          endif

!-------------------------------------------------------------------
!    allocate pressure, index arrays used in rctrns (if needed)
!-------------------------------------------------------------------
          if (callrctrns_ch4) then
            call allocate_interp_arrays
          endif
          
        end if ! phase=0 .or. phase=1
 
!---------------------------------------------------------------------
!    loop on frequency bands. in the 1996 SEA formulation, there are
!    1 frequency ranges for lbl ch4 transmissions:
!    nf = 1:  lbl transmissions over 1200-1400 cm-1    
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!    read in ch4 transmissivities at this point-----
!    data is read for all temperature profiles required for the
!    frequency band (at 1 or 2 appropriate concentrations). the 
!    number of temperature profiles for band (nf) is ntbnd(nf).
!    in the 1996 SEA formulation, the profiles required are 3 
!    (USSTD,1976; USSTD,1976 +- 25).
!----------------------------------------------------------------------
        do nf = 1,nfreq_bands_sea_ch4
            
          origin=mod(nf-1+offset,nproc)
          nf_offset=nf+offset
          
          if ( phase==0 .or. phase==-1 ) then
              
            if ( myid==origin ) then  
              if (trim(Lw_control%linecatalog_form) == 'hitran_2012' ) then  
                call read_lbltfs ('ch4', callrctrns_ch4, nstd_ch4_lo, &
                                  nstd_ch4_hi, nf, ntbnd_ch4, & 
                                  trns_std_hi_nf, trns_std_lo_nf )
              else
                call read_lbltfs_old ('ch4', callrctrns_ch4, nstd_ch4_lo, &
                                  nstd_ch4_hi, nf, ntbnd_ch4, & 
                                  trns_std_hi_nf, trns_std_lo_nf )
              end if  
              do_lyrcalc_ch4    = do_lyrcalc_ch4_nf(nf)
              do_lvlcalc_ch4    = do_lvlcalc_ch4_nf(nf)
              do_lvlctscalc_ch4 = do_lvlctscalc_ch4_nf(nf)
 
!---------------------------------------------------------------------
!    load in appropriate ch4 transmission functions
!---------------------------------------------------------------------
              if (ch4_vmr /= 0.0) then
                do nt = 1,ntbnd_ch4(nf)   ! temperature structure loop.
                  trns_std_hi = trns_std_hi_nf(:,:,nt)
                  if (callrctrns_ch4) then
                    trns_std_lo = trns_std_lo_nf(:,:,nt)
                  endif
                  call gasint(gas_type,           &
                              ch4_vmr, ch4_std_lo, ch4_std_hi,   &
                              callrctrns_ch4,   &
                              do_lvlcalc_ch4, do_lvlctscalc_ch4,    &
                              do_lyrcalc_ch4, nf, nt)
                  trns_interp_lyr_ps_nf(:,:,nt) = trns_interp_lyr_ps(:,:)
                  trns_interp_lyr_ps8_nf(:,:,nt) = trns_interp_lyr_ps8(:,:)
                  trns_interp_lvl_ps_nf(:,:,nt) = trns_interp_lvl_ps(:,:)
                  trns_interp_lvl_ps8_nf(:,:,nt) = trns_interp_lvl_ps8(:,:)
                enddo   ! temperature structure loop
              endif
 
!--------------------------------------------------------------------
!    perform final processing for each frequency band.
!--------------------------------------------------------------------
              if (ch4_vmr /= 0.0) then
                call gasins(gas_type, do_lvlcalc_ch4, do_lvlctscalc_ch4,      &
                            do_lyrcalc_ch4, nf, ntbnd_ch4(nf), ndimkp, ndimk, &
                            dgasdt10_lvl(:,:,nf_offset), dgasdt10_lvlcts(:,nf_offset), dgasdt10_lyr(:,:,nf_offset),      &
                            gasp10_lvl(:,:,nf_offset), gasp10_lvlcts(:,nf_offset), gasp10_lyr(:,:,nf_offset),            &
                            d2gast10_lvl(:,:,nf_offset), d2gast10_lvlcts(:,nf_offset), d2gast10_lyr(:,:,nf_offset),      &
                            dgasdt8_lvl(:,:,nf_offset),  dgasdt8_lvlcts(:,nf_offset),  dgasdt8_lyr(:,:,nf_offset) ,      &
                            gasp8_lvl(:,:,nf_offset),  gasp8_lvlcts(:,nf_offset),  gasp8_lyr(:,:,nf_offset) ,            &
                            d2gast8_lvl(:,:,nf_offset),  d2gast8_lvlcts(:,nf_offset),  d2gast8_lyr(:,:,nf_offset) )
 
              else
!---------------------------------------------------------------------
!    define arrays for the SEA module. the SEA model nomenclature
!    has been used here and the values of do_lvlcalc, do_lvlctscalc,
!    and do_lyrcalc are assumed to be from the data statement.
!---------------------------------------------------------------------
!15
                gasp10_lyr(:,:,nf_offset) = 1.0 
                gasp8_lyr(:,:,nf_offset) = 1.0 
                dgasdt10_lyr(:,:,nf_offset) = 0.0
                dgasdt8_lyr(:,:,nf_offset) = 0.0
                d2gast10_lyr(:,:,nf_offset) = 0.0
                d2gast8_lyr(:,:,nf_offset) = 0.0
              endif
             
            end if ! myid==origin
         
          end if ! phase==0
          
          if ( phase==1 .or. phase==-1 ) then  
             
            call ccmpi_bcastr8(gasp10_lyr(:,:,nf_offset),origin,comm_world) 
            call ccmpi_bcastr8(gasp8_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(dgasdt10_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(dgasdt8_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(d2gast10_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(d2gast8_lyr(:,:,nf_offset),origin,comm_world)
           
            call put_ch4_stdtf_for_gas_tf (gasp10_lyr(:,:,nf_offset), gasp8_lyr(:,:,nf_offset),      &
                                           dgasdt10_lyr(:,:,nf_offset), dgasdt8_lyr(:,:,nf_offset),  &
                                           d2gast10_lyr(:,:,nf_offset),  d2gast8_lyr(:,:,nf_offset))
          end if ! phase==1 .or. phase==-1
         
        enddo  !  frequency band loop
        
!--------------------------------------------------------------------
!    deallocate pressure, index arrays used in rctrns (if needed)
!--------------------------------------------------------------------
        
        if ( phase==0 .or. phase==-1 ) then
        
          if (callrctrns_ch4) then
            call deallocate_interp_arrays
          endif

          deallocate(trns_interp_lyr_ps)
          deallocate(trns_interp_lyr_ps8)
          deallocate(trns_interp_lvl_ps)
          deallocate(trns_interp_lvl_ps8)
          deallocate(trns_interp_lyr_ps_nf)
          deallocate(trns_interp_lyr_ps8_nf)
          deallocate(trns_interp_lvl_ps_nf)
          deallocate(trns_interp_lvl_ps8_nf)
          deallocate(trns_std_hi)
          deallocate(trns_std_lo)
          
        end if ! phase==0 .or. phase==-1

!-----------------------------------------------------------------
!    pass necessary data to gas_tf in case stdtf file is to be written
!-----------------------------------------------------------------
        if ( phase==1 .or. phase==-1 ) then
          call process_ch4_input_file (gas_type, ch4_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase==1 .or. phase==-1
          
!---------------------------------------------------------------------
!    if not calculating tfs, read them in
!---------------------------------------------------------------------
      else                                      
          
        if ( phase==0 .or. phase==-1 ) then          
          call process_ch4_input_file (gas_type, ch4_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase==0 .or. phase==-1

      endif   ! (do_calcstdch4tfs)

      
      deallocate( trns_std_hi_nf )
      deallocate( trns_std_lo_nf )
!-------------------------------------------------------------------
 
 
end subroutine ch4_lblinterp



!####################################################################
! <SUBROUTINE NAME="co2_lblinterp">
!  <OVERVIEW>
!   Subroutine to interpolate co2 transmission function to user
!   specified pressure levels and co2 concentration
!  </OVERVIEW>
!  <DESCRIPTION>
!   this routine is 1) a standalone program for a co2 interpolation
!     to user-specified pressure levels and co2 concentration;
!     2) an interface between a GCM and co2 interpolation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call co2_lblinterp (co2_vmr)
!  </TEMPLATE>
!  <IN NAME="co2_vmr" TYPE="real">
!   co2 volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine co2_lblinterp (co2_vmr, phase)

use cc_mpi
use newmpar_m, only : nproc

!--------------------------------------------------------------------
!    this routine is 
!    1) a standalone program for a co2 interpolation
!       to user-specified pressure levels and co2 concentration;
!    2) an interface between a GCM and co2 interpolation
!
!    input files:
!
!         1)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (co2_vmr).
!         2)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (co2_vmr). may not be
!                  used, depending on value of (co2_vmr).
!
!    output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these 
!                  files are created if ifdef (writeinterpco2) is
!                  on. otherwise, it is assumed the data is fed
!                  directly back to the parent model.
!--------------------------------------------------------------------

real,             intent(in)     ::  co2_vmr
integer,          intent(in)     :: phase

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    co2_vmr
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
      logical              ::  callrctrns_co2
      logical              ::  do_lyrcalc_co2
      logical              ::  do_lvlcalc_co2
      logical              ::  do_lvlctscalc_co2
      real                 ::  co2_std_lo, co2_std_hi
      integer              ::  n, nf, nt
      integer              ::  nstd_co2_lo, nstd_co2_hi
      integer              ::  origin, nf_offset
      integer              ::  nfco2
      integer, parameter   ::  offset=nfreq_bands_sea_ch4
      character(len=8)     ::  gas_type = 'co2'
      real, dimension(:,:,:), allocatable, save :: trns_std_hi_nf, trns_std_lo_nf
      real, dimension(:,:), allocatable, save :: dum_lvlcts
      
!---------------------------------------------------------------------
!  local variables
!
!     callrctrns_co2
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod',   &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
      allocate( trns_std_hi_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
      allocate( trns_std_lo_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
      
!---------------------------------------------------------------------
!    this routine does nothing unless the calculation of tfs is desired.
!---------------------------------------------------------------------
      if (do_calcstdco2tfs) then

        if ( phase==0 .or. phase==-1 ) then

          rctrns_first = .true.
          intcoef_2d_std_first = .true.
          
          allocate (trns_std_hi(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_std_lo(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_interp_lyr_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lyr_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          allocate (trns_interp_lvl_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lvl_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          
!--------------------------------------------------------------------
!    using the value of the co2 volume mixing ratio (co2_vmr) and
!    the available standard co2 mixing ratios (with lbl 
!    transmission functions) obtain the two standard mixing ratios
!    which bracket (co2_vmr). if (as in a fixed co2 experiment) the
!    difference between (co2_vmr) and the higher of the standard
!    mixing ratios is less than a tolerance (taken as .0001 ppmv)
!    we will assume that that standard co2 transmissivity applies
!    to (co2_vmr) without interpolation. otherwise, interpolation
!    to (co2_vmr) will be performed, in rctrns.F
!--------------------------------------------------------------------- 
          if (co2_vmr .LT. co2_std_vmr(1) .OR.                     &
            co2_vmr .GT. co2_std_vmr(number_std_co2_vmrs)) then
            write(6,*) "lw_gases_stdtf_mod: co2 volume mixing ratio is out of range"
            stop
          endif

          if (co2_vmr .EQ. co2_std_vmr(1)) then
            co2_std_lo = co2_std_vmr(1)
            co2_std_hi = co2_std_vmr(1)
            nstd_co2_lo = 1
            nstd_co2_hi = 1
          else 
            do n=1,number_std_co2_vmrs-1
              if ( co2_vmr .GT. co2_std_vmr(n) .AND.                 &  
                   co2_vmr .LE. co2_std_vmr(n+1)) then
                co2_std_lo = co2_std_vmr(n)
                co2_std_hi = co2_std_vmr(n+1)
                nstd_co2_lo = n
                nstd_co2_hi = n+1
                exit
              endif
            enddo
          endif 

!-------------------------------------------------------------------
!    co2_std_lo, nstd_co2_lo have arbitrary definitions, since they
!    will not be used, as callrctrns will be false in this case.
!-------------------------------------------------------------------
          if (ABS(co2_vmr - co2_std_hi) .LE. 1.0e-4) then
            callrctrns_co2 = .false.
          else
            callrctrns_co2 = .true.
          endif

!-------------------------------------------------------------------
!    allocate pressure, index arrays used in rctrns (if needed)
!-------------------------------------------------------------------
          if (callrctrns_co2) then
            call allocate_interp_arrays
          endif
        
        end if ! phase==0 .or. phase==-1

!--------------------------------------------------------------------
!    loop on frequency bands. in the 1996 SEA formulation, there are
!    8 frequency ranges for lbl co2 transmissions:
!    nf = 1:  lbl transmissions over 490-850 cm-1    
!    nf = 2:  lbl transmissions over 490-630 cm-1    
!    nf = 3:  lbl transmissions over 630-700 cm-1    
!    nf = 4:  lbl transmissions over 700-800 cm-1    
!    nf = 5:  lbl transmissions over 2270-2380 cm-1    
!    nf = 6:  lbl transmissions over 990-1070 cm-1    
!    nf = 7:  lbl transmissions over 900-990 cm-1    
!    nf = 8:  lbl transmissions over 1070-1200 cm-1
!---------------------------------------------------------------------
        allocate( dum_lvlcts(KSRAD:KERAD+1,1:6) )
        
        if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
           nfco2 = 5  ! only the first 5 bands are used, not the 10 um band
        else 
           nfco2 = nfreq_bands_sea_co2
        endif
        
        do nf = 1,nfco2
 
!---------------------------------------------------------------------
!    read in co2 transmissivities at this point-----
!    data is read for all temperature profiles required for the
!    frequency band (at 1 or 2 appropriate concentrations). the 
!    number of temperature profiles for band (nf) is ntbnd_co2(nf).
!    in the 1996 SEA formulation, the profiles required are 3 
!    (USSTD,1976; USSTD,1976 +- 25) except for the 4.3 um band (nf=2)
!    where the number is one.
!---------------------------------------------------------------------

          origin=mod(nf-1+offset,nproc)
          nf_offset=nf+offset

          if ( phase==0 .or. phase==-1 ) then

            if ( myid==origin ) then
              if (trim(Lw_control%linecatalog_form) == 'hitran_2012' ) then
                call read_lbltfs('co2',                                     &
                                 callrctrns_co2, nstd_co2_lo, nstd_co2_hi,  &
                                 nf, ntbnd_co2,                             &
                                 trns_std_hi_nf, trns_std_lo_nf )
              else
                call read_lbltfs_old('co2',                                     &
                                 callrctrns_co2, nstd_co2_lo, nstd_co2_hi,  &
                                 nf, ntbnd_co2,                             &
                                 trns_std_hi_nf, trns_std_lo_nf )
              end if    
              do_lyrcalc_co2 = do_lyrcalc_co2_nf(nf)
              do_lvlcalc_co2 = do_lvlcalc_co2_nf(nf)
              do_lvlctscalc_co2 = do_lvlctscalc_co2_nf(nf)
 
!--------------------------------------------------------------------
!    load in appropriate co2 transmission functions
!--------------------------------------------------------------------
              if (co2_vmr /= 0.0) then
                do nt = 1,ntbnd_co2(nf)    !  temperature structure loop.
                  trns_std_hi = trns_std_hi_nf(:,:,nt)
                  if (callrctrns_co2) then
                    trns_std_lo = trns_std_lo_nf(:,:,nt)
                  endif
                  call gasint(gas_type,                             &
                              co2_vmr, co2_std_lo, co2_std_hi,      &
                              callrctrns_co2,                       & 
                              do_lvlcalc_co2, do_lvlctscalc_co2,    &
                              do_lyrcalc_co2,   &
                              nf, nt)
                  
                  trns_interp_lyr_ps_nf(:,:,nt) = trns_interp_lyr_ps(:,:)
                  trns_interp_lyr_ps8_nf(:,:,nt) = trns_interp_lyr_ps8(:,:)
                  trns_interp_lvl_ps_nf(:,:,nt) = trns_interp_lvl_ps(:,:)
                  trns_interp_lvl_ps8_nf(:,:,nt) = trns_interp_lvl_ps8(:,:)
                enddo        !  temperature structure loop
 
              endif
!--------------------------------------------------------------------
!    perform final processing for each frequency band.
!--------------------------------------------------------------------
              if (co2_vmr /= 0.0) then
                call gasins('co2',                                 &
                            do_lvlcalc_co2, do_lvlctscalc_co2,    &
                            do_lyrcalc_co2, &
                            nf, ntbnd_co2(nf),                            & 
                            ndimkp,ndimk,                                 & 
                            dgasdt10_lvl(:,:,nf_offset), dgasdt10_lvlcts(:,nf_offset), dgasdt10_lyr(:,:,nf_offset),   & 
                            gasp10_lvl(:,:,nf_offset), gasp10_lvlcts(:,nf_offset), gasp10_lyr(:,:,nf_offset),        & 
                            d2gast10_lvl(:,:,nf_offset), d2gast10_lvlcts(:,nf_offset), d2gast10_lyr(:,:,nf_offset),   & 
                            dgasdt8_lvl(:,:,nf_offset),  dgasdt8_lvlcts(:,nf_offset),  dgasdt8_lyr(:,:,nf_offset) ,    & 
                            gasp8_lvl(:,:,nf_offset),  gasp8_lvlcts(:,nf_offset),  gasp8_lyr(:,:,nf_offset) ,        & 
                            d2gast8_lvl(:,:,nf_offset),  d2gast8_lvlcts(:,nf_offset),  d2gast8_lyr(:,:,nf_offset) )
              else
                dgasdt10_lvl(:,:,nf_offset) = 0.
                dgasdt10_lvlcts(:,nf_offset)  = 0.
                dgasdt10_lyr(:,:,nf_offset) = 0.
                gasp10_lvl(:,:,nf_offset)  = 1.
                gasp10_lvlcts(:,nf_offset) = 1.
                gasp10_lyr(:,:,nf_offset) = 1.
                d2gast10_lvl(:,:,nf_offset) = 0.
                d2gast10_lvlcts(:,nf_offset)   = 0.
                d2gast10_lyr(:,:,nf_offset) = 0.
                dgasdt8_lvl(:,:,nf_offset) = 0.
                dgasdt8_lvlcts(:,nf_offset) = 0.
                dgasdt8_lyr(:,:,nf_offset)      = 0.
                gasp8_lvl(:,:,nf_offset)  = 1.
                gasp8_lvlcts(:,nf_offset)  = 1.
                gasp8_lyr(:,:,nf_offset)  = 1.
                d2gast8_lvl(:,:,nf_offset)  = 0.
                d2gast8_lvlcts(:,nf_offset)  = 0.
                d2gast8_lyr(:,:,nf_offset) = 0.
              endif
              
            end if ! myid==origin
            
          end if ! phase==0 .or. phase==-1
 
!--------------------------------------------------------------------
!    define arrays for the SEA module. the SEA model nomenclature
!    has been used here and the values of do_lvlcalc, do_lvlctscalc,
!    and do_lyrcalc are assumed to be from the data statement.
!----------------------------------------------------------------------
          
          if ( phase==1 .or. phase==-1 ) then
          
            if (nf == 1 .or. nf == 5 .or. nf == 6 .or. nf == 7 .or. nf == 8 ) then
              call ccmpi_bcastr8(gasp10_lyr(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(gasp8_lyr(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(dgasdt10_lyr(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(dgasdt8_lyr(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(d2gast10_lyr(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(d2gast8_lyr(:,:,nf_offset),origin,comm_world)
              call put_co2_stdtf_for_gas_tf (nf, gasp10_lyr(:,:,nf_offset), gasp8_lyr(:,:,nf_offset), &
                                             dgasdt10_lyr(:,:,nf_offset), dgasdt8_lyr(:,:,nf_offset), &
                                             d2gast10_lyr(:,:,nf_offset),  d2gast8_lyr(:,:,nf_offset))
            endif
            if (nf <= 4 .or. nf == 6 .or. nf == 7 .or. nf == 8 ) then
              call ccmpi_bcastr8(gasp10_lvl(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(dgasdt10_lvl(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(d2gast10_lvl(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(gasp8_lvl(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(dgasdt8_lvl(:,:,nf_offset),origin,comm_world)
              call ccmpi_bcastr8(d2gast8_lvl(:,:,nf_offset),origin,comm_world)
              if ( myid==origin ) then
                dum_lvlcts(:,1) = gasp10_lvlcts(:,nf_offset)
                dum_lvlcts(:,2) = gasp8_lvlcts(:,nf_offset)
                dum_lvlcts(:,3) = dgasdt10_lvlcts(:,nf_offset)
                dum_lvlcts(:,4) = dgasdt8_lvlcts(:,nf_offset)
                dum_lvlcts(:,5) = d2gast10_lvlcts(:,nf_offset)
                dum_lvlcts(:,6) = d2gast8_lvlcts(:,nf_offset)
              end if
              call ccmpi_bcastr8(dum_lvlcts(:,1:6),origin,comm_world)
              gasp10_lvlcts(:,nf_offset)   = dum_lvlcts(:,1)
              gasp8_lvlcts(:,nf_offset)    = dum_lvlcts(:,2)
              dgasdt10_lvlcts(:,nf_offset) = dum_lvlcts(:,3)
              dgasdt8_lvlcts(:,nf_offset)  = dum_lvlcts(:,4)
              d2gast10_lvlcts(:,nf_offset) = dum_lvlcts(:,5)
              d2gast8_lvlcts(:,nf_offset)  = dum_lvlcts(:,6)
              call put_co2_nbltf_for_gas_tf (nf, gasp10_lvl(:,:,nf_offset),   &
                                           dgasdt10_lvl(:,:,nf_offset), d2gast10_lvl(:,:,nf_offset), &
                                           gasp8_lvl(:,:,nf_offset), dgasdt8_lvl(:,:,nf_offset),  &
                                           d2gast8_lvl(:,:,nf_offset), gasp10_lvlcts(:,nf_offset), &
                                           gasp8_lvlcts(:,nf_offset),&
                                           dgasdt10_lvlcts(:,nf_offset),  &
                                           dgasdt8_lvlcts(:,nf_offset),&
                                           d2gast10_lvlcts(:,nf_offset), &
                                           d2gast8_lvlcts(:,nf_offset))
            endif 
            
          end if ! phase==1 .or. phase==-1  
            
        enddo  ! frequency band loop
        deallocate( dum_lvlcts )

!-----------------------------------------------------------------
!    deallocate pressure, index arrays used in rctrns (if needed)
!-----------------------------------------------------------------

        if ( phase==0 .or. phase==-1 ) then
        
          if (callrctrns_co2) then
            call deallocate_interp_arrays
          endif

          deallocate(trns_interp_lyr_ps)
          deallocate(trns_interp_lyr_ps8)
          deallocate(trns_interp_lvl_ps)
          deallocate(trns_interp_lvl_ps8)
          deallocate(trns_interp_lyr_ps_nf)
          deallocate(trns_interp_lyr_ps8_nf)
          deallocate(trns_interp_lvl_ps_nf)
          deallocate(trns_interp_lvl_ps8_nf)
          deallocate(trns_std_hi)
          deallocate(trns_std_lo)
          
        end if ! phase==0 .or. phase==-1
        
!-----------------------------------------------------------------
!    pass necessary data to gas_tf in case stdtf file is to be written
!-----------------------------------------------------------------
        if ( phase==1 .or. phase==-1 ) then
          call process_co2_input_file (gas_type, co2_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase==1 .or. phase==-1

!---------------------------------------------------------------------
!    if not calculating tfs, read them in
!---------------------------------------------------------------------
      else                              
        if ( phase==0 .or. phase==-1 ) then 
          call process_co2_input_file (gas_type, co2_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase ==0 .or. phase==-1
      endif    !   (do_calcstdco2tfs)

!-------------------------------------------------------------------
      deallocate( trns_std_hi_nf )
      deallocate( trns_std_lo_nf )


end subroutine co2_lblinterp



!#####################################################################
! <SUBROUTINE NAME="n2o_lblinterp">
!  <OVERVIEW>
!   Subroutine to interpolate n2o transmission function to user
!   specified pressure levels and n2o concentration
!  </OVERVIEW>
!  <DESCRIPTION>
!   this routine is 1) a standalone program for a n2o interpolation
!     to user-specified pressure levels and n2o concentration;
!     2) an interface between a GCM and n2o interpolation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call n2o_lblinterp (n2o_vmr)
!  </TEMPLATE>
!  <IN NAME="n2o_vmr" TYPE="real">
!   n2o volume mixing ratio
!  </IN>
! </SUBROUTINE>
!
subroutine n2o_lblinterp (n2o_vmr, phase)

use cc_mpi
use newmpar_m, only : nproc

!---------------------------------------------------------------------
!    n2o_lblinterp is 
!    1) a standalone program for a n2o interpolation
!       to user-specified pressure levels and n2o concentration;
!    2) an interface between a GCM and n2o interpolation
!
!    input files:
!
!         1)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (n2o_vmr).
!         2)     : gas transmission function at higher of 2 
!                  standard mixing ratios, for a specified frequency
!                  range, determined by choice of (n2o_vmr). may not be
!                  used, depending on value of (n2o_vmr).
!
!    output files:
!
!         id2,   : interpolated gas transmission fctns and derivatives
!       id2nb      saved in format suitable as input to operational
!                  radiation program, for the desired gas mixing ratio
!                  and frequency range. The number of records will
!                  vary, depending on the frequency range. these 
!                  files are created if ifdef (writeinterpn2o) is
!                  on. otherwise, it is assumed the data is fed
!                 directly back to the parent model.
!---------------------------------------------------------------------

real,             intent(in)   :: n2o_vmr
integer,          intent(in)   :: phase

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     n2o_vmr
!
!--------------------------------------------------------------------
        
!--------------------------------------------------------------------
!    local variables

      logical                   ::  callrctrns_n2o
      logical                   ::  do_lyrcalc_n2o
      logical                   ::  do_lvlcalc_n2o
      logical                   ::  do_lvlctscalc_n2o
      real                      ::  n2o_std_lo, n2o_std_hi
      integer                   ::  n, nf, nt
      integer                   ::  nstd_n2o_lo, nstd_n2o_hi
      integer                   ::  origin, nf_offset
      integer, parameter        ::  offset=nfreq_bands_sea_ch4+nfreq_bands_sea_co2
      character(len=8)          ::  gas_type = 'n2o'
      real, dimension(:,:,:), allocatable, save :: trns_std_hi_nf, trns_std_lo_nf

!--------------------------------------------------------------------
!    local variables
!
!      callrctrns_n2o
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
      
      allocate( trns_std_hi_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
      allocate( trns_std_lo_nf(NSTDCO2LVLS,NSTDCO2LVLS,3) )
 
!---------------------------------------------------------------------
!    this routine does nothing unless the calculation of tfs is desired.
!---------------------------------------------------------------------
      if (do_calcstdn2otfs) then
          
        if ( phase==0 .or. phase==-1 ) then

          rctrns_first = .true.
          intcoef_2d_std_first = .true.

          allocate (trns_std_hi(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_std_lo(NSTDCO2LVLS, NSTDCO2LVLS) )
          allocate (trns_interp_lyr_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lvl_ps8(KSRAD:KERAD+1, KSRAD:KERAD+1) )
          allocate (trns_interp_lyr_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lyr_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          allocate (trns_interp_lvl_ps_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3) )
          allocate (trns_interp_lvl_ps8_nf(KSRAD:KERAD+1,KSRAD:KERAD+1,3))
          
!--------------------------------------------------------------------
!    using the value of the n2o volume mixing ratio (n2o_vmr) and
!    the available standard n2o mixing ratios (with lbl 
!    transmission functions) obtain the two standard mixing ratios
!    which bracket (n2o_vmr). if (as in a fixed n2o experiment) the
!    difference between (n2o_vmr) and the higher of the standard
!    mixing ratios is less than a tolerance (taken as 0.1 ppbv)
!    we will assume that that standard n2o transmissivity applies
!    to (n2o_vmr) without interpolation. otherwise, interpolation
!    to (n2o_vmr) will be performed, in rctrns.F
!---------------------------------------------------------------------
          if (n2o_vmr .LT. n2o_std_vmr(1) .OR.     &
               n2o_vmr .GT. n2o_std_vmr(number_std_n2o_vmrs)) then
            !call error_mesg ('lw_gases_stdtf_mod', &
            !     'n2o volume mixing ratio is out of range', FATAL)
            write(6,*) "lw_gases_stdtf_mod n2o volume mixing ratio is out of range"
            write(6,*) n2o_vmr,n2o_std_vmr(1),n2o_std_vmr(number_std_n2o_vmrs)
            stop
          endif

          if (n2o_vmr .EQ. n2o_std_vmr(1)) then
            n2o_std_lo = n2o_std_vmr(1)
            n2o_std_hi = n2o_std_vmr(1)
            nstd_n2o_lo = 1
            nstd_n2o_hi = 1
          else 
            do n=1,number_std_n2o_vmrs-1
              if ( n2o_vmr .GT. n2o_std_vmr(n) .AND.               &
                   n2o_vmr .LE. n2o_std_vmr(n+1)) then
                n2o_std_lo = n2o_std_vmr(n)
                n2o_std_hi = n2o_std_vmr(n+1)
                nstd_n2o_lo = n
                nstd_n2o_hi = n+1
                exit
              endif
            enddo
          endif

!-------------------------------------------------------------------
!    n2o_std_lo, nstd_n2o_lo have arbitrary definitions, since they
!    will not be used, as callrctrns will be false in this case.
!-------------------------------------------------------------------
          if (ABS(n2o_vmr - n2o_std_hi) .LE. 1.0e-1) then
            callrctrns_n2o = .false.
          else
           callrctrns_n2o = .true.
          endif
 
!--------------------------------------------------------------------
!    allocate pressure, index arrays used in rctrns (if needed)
!--------------------------------------------------------------------
          if (callrctrns_n2o) then
            call allocate_interp_arrays
          endif
          
        end if ! phase==0 .or. phase==-1
          
!---------------------------------------------------------------------
!    loop on frequency bands. in the 1996 SEA formulation, there are
!    3 frequency ranges for lbl n2o transmissions:
!    nf = 1:  lbl transmissions over 1200-1400 cm-1    
!    nf = 2:  lbl transmissions over 1070-1200 cm-1    
!    nf = 3:  lbl transmissions over 560-630 cm-1    
!---------------------------------------------------------------------
        do nf = 1,nfreq_bands_sea_n2o
 
!---------------------------------------------------------------------
!    read in n2o transmissivities at this point-----
!    data is read for all temperature profiles required for the
!    frequency band (at 1 or 2 appropriate concentrations). the 
!    number of temperature profiles for band (nf) is ntbnd(nf).
!    in the 1996 SEA formulation, the profiles required are 3 
!    (USSTD,1976; USSTD,1976 +- 25).
!----------------------------------------------------------------------
            
          origin=mod(nf-1+offset,nproc)
          nf_offset=nf+offset

          if ( phase==0 .or. phase==-1 ) then

            if ( myid==origin ) then 
              if (trim(Lw_control%linecatalog_form) == 'hitran_2012' ) then
                call read_lbltfs('n2o',                                       &   
                                 callrctrns_n2o, nstd_n2o_lo, nstd_n2o_hi, nf,&
                                 ntbnd_n2o,                                   &
                                 trns_std_hi_nf, trns_std_lo_nf )
              else
                call read_lbltfs_old('n2o',                                       &   
                                 callrctrns_n2o, nstd_n2o_lo, nstd_n2o_hi, nf,&
                                 ntbnd_n2o,                                   &
                                 trns_std_hi_nf, trns_std_lo_nf )
              end if     
              do_lyrcalc_n2o = do_lyrcalc_n2o_nf(nf)
              do_lvlcalc_n2o = do_lvlcalc_n2o_nf(nf)
              do_lvlctscalc_n2o = do_lvlctscalc_n2o_nf(nf)
 
!----------------------------------------------------------------------
!    load in appropriate n2o transmission functions
!----------------------------------------------------------------------
              if (n2o_vmr /= 0.0) then
                do nt = 1,ntbnd_n2o(nf) ! temperature structure loop
                  trns_std_hi = trns_std_hi_nf(:,:,nt)
                  if (callrctrns_n2o) then
                    trns_std_lo = trns_std_lo_nf(:,:,nt)
                  endif
                  call gasint(gas_type, n2o_vmr, n2o_std_lo, n2o_std_hi,    &
                              callrctrns_n2o,     &
                              do_lvlcalc_n2o, do_lvlctscalc_n2o,    &
                              do_lyrcalc_n2o, nf, nt)
                  trns_interp_lyr_ps_nf(:,:,nt) = trns_interp_lyr_ps(:,:)
                  trns_interp_lyr_ps8_nf(:,:,nt) = trns_interp_lyr_ps8(:,:)
                  trns_interp_lvl_ps_nf(:,:,nt) = trns_interp_lvl_ps(:,:)
                  trns_interp_lvl_ps8_nf(:,:,nt) = trns_interp_lvl_ps8(:,:)
                enddo    ! temperature structure loop
              endif 

!--------------------------------------------------------------------
!    perform final processing for each frequency band.
!--------------------------------------------------------------------
              if (n2o_vmr /= 0.0) then
                call gasins('n2o',                                        &
                            do_lvlcalc_n2o, do_lvlctscalc_n2o,   &
                            do_lyrcalc_n2o, nf, ntbnd_n2o(nf),   &
                            ndimkp,ndimk,                               &
                            dgasdt10_lvl(:,:,nf_offset), dgasdt10_lvlcts(:,nf_offset), dgasdt10_lyr(:,:,nf_offset),   &
                            gasp10_lvl(:,:,nf_offset), gasp10_lvlcts(:,nf_offset), gasp10_lyr(:,:,nf_offset),      &
                            d2gast10_lvl(:,:,nf_offset), d2gast10_lvlcts(:,nf_offset), d2gast10_lyr(:,:,nf_offset),  &  
                            dgasdt8_lvl(:,:,nf_offset),  dgasdt8_lvlcts(:,nf_offset),  dgasdt8_lyr(:,:,nf_offset) ,    &
                            gasp8_lvl(:,:,nf_offset),  gasp8_lvlcts(:,nf_offset),  gasp8_lyr(:,:,nf_offset) ,      &
                            d2gast8_lvl(:,:,nf_offset),  d2gast8_lvlcts(:,nf_offset),  d2gast8_lyr(:,:,nf_offset) )
              else
!15
               gasp10_lyr(:,:,nf_offset) = 1.0 
               gasp8_lyr(:,:,nf_offset) = 1.0 
               dgasdt10_lyr(:,:,nf_offset) = 0.0
               dgasdt8_lyr(:,:,nf_offset) = 0.0
               d2gast10_lyr(:,:,nf_offset) = 0.0
               d2gast8_lyr(:,:,nf_offset) = 0.0
             endif
              
           end if ! myid==origin   
              
         end if  ! phase==0 .or phase==-1   
 
!--------------------------------------------------------------------
!    define arrays for the SEA module. the SEA model nomenclature
!    has been used here and the values of do_lvlcalc, do_lvlctscalc,
!    and do_lyrcalc are assumed to be from the data statement.
!--------------------------------------------------------------------
         
         if ( phase==1 .or. phase==-1 ) then 

            call ccmpi_bcastr8(gasp10_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(gasp8_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(dgasdt10_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(dgasdt8_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(d2gast10_lyr(:,:,nf_offset),origin,comm_world)
            call ccmpi_bcastr8(d2gast8_lyr(:,:,nf_offset),origin,comm_world)
             
           call put_n2o_stdtf_for_gas_tf (nf, gasp10_lyr(:,:,nf_offset), gasp8_lyr(:,:,nf_offset),   &
                                          dgasdt10_lyr(:,:,nf_offset), dgasdt8_lyr(:,:,nf_offset), &
                                          d2gast10_lyr(:,:,nf_offset), d2gast8_lyr(:,:,nf_offset))
           
          end if ! phase==1 .or. phase==-1
           
        enddo    ! frequency band loop
        
!---------------------------------------------------------------------
!    deallocate pressure, index arrays used in rctrns (if needed)
!--------------------------------------------------------------------

        if ( phase==0 .or. phase==-1 ) then
        
          if (callrctrns_n2o) then
            call deallocate_interp_arrays
          endif
          
          deallocate(trns_interp_lyr_ps)
          deallocate(trns_interp_lyr_ps8)
          deallocate(trns_interp_lvl_ps)
          deallocate(trns_interp_lvl_ps8)
          deallocate(trns_interp_lyr_ps_nf)
          deallocate(trns_interp_lyr_ps8_nf)
          deallocate(trns_interp_lvl_ps_nf)
          deallocate(trns_interp_lvl_ps8_nf)
          deallocate(trns_std_hi)
          deallocate(trns_std_lo)
          
        end if ! phase==0 .or. phase==-1
        
!-----------------------------------------------------------------
!    pass necessary data to gas_tf in case stdtf file is to be written
!-----------------------------------------------------------------
        if ( phase==1 .or. phase==-1 ) then
          call process_n2o_input_file (gas_type, n2o_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase==1 .or. phase==-1

!---------------------------------------------------------------------
!   if not calculating tfs, read them in
!---------------------------------------------------------------------
      else   
          
        if ( phase==0 .or. phase==-1 ) then  
          call process_n2o_input_file (gas_type, n2o_vmr, NSTDCO2LVLS, &
                                       KSRAD, KERAD, pd, plm, pa)
        end if ! phase==0 .or. phase==-1
        
      endif
 
!------------------------------------------------------------------
      deallocate( trns_std_hi_nf )
      deallocate( trns_std_lo_nf )


end subroutine n2o_lblinterp





!#####################################################################
! <SUBROUTINE NAME="lw_gases_stdtf_dealloc">
!  <OVERVIEW>
!   Subroutine to deallocate long wave gas transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to deallocate long wave gas transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_gases_stdtf_dealloc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_gases_stdtf_dealloc

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif

end subroutine lw_gases_stdtf_dealloc



!####################################################################
! <SUBROUTINE NAME="cfc_exact">
!  <OVERVIEW>
!   cfc_exact computes exact cool-to-space transmission function 
!   for cfc for the desired band (given by index). 
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_exact computes exact cool-to-space transmission function 
!   for cfc for the desired band (given by index).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_exact (index, Optical, cfc_tf)
!  </TEMPLATE>
!  <IN NAME="index" TYPE="integer">
!   the spectral index where exact CTS transmision function is computed
!  </IN>
!  <IN NAME="Optical" TYPE="optical_depth_type">
!   The CFC gas optical depth
!  </IN>
!  <OUT NAME="cfc_tf" TYPE="real">
!   exact CTS transmission function output
!  </OUT>
! </SUBROUTINE>
!
subroutine cfc_exact (index, Optical, cfc_tf)

!----------------------------------------------------------------------
!    cfc_exact computes exact cool-to-space transmission function 
!    for cfc for the desired band (given by index). 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     index
!     Optical
!
!  intent(out) variables:
!
!     cfc_tf
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer :: kx   ! do-loop  index

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod',   &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kx = size (Optical%totf11,3) 
      cfc_tf(:,:,:          ) = 1.0E+00 -    &
                      strf113(index)*Optical%totf113(:,:,2:kx) -   &
                      strf22 (index)*Optical%totf22 (:,:,2:kx) -   &
                      strf11 (index)*Optical%totf11 (:,:,2:kx) -   &
                      strf12 (index)*Optical%totf12 (:,:,2:kx)    

!---------------------------------------------------------------------

end subroutine cfc_exact




!####################################################################
! <SUBROUTINE NAME="cfc_exact_part">
!  <OVERVIEW>
!   cfc_exact computes exact cool-to-space transmission function 
!   at levels below klevel for cfc for the band given by index. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_exact computes exact cool-to-space transmission function 
!   at levels below klevel for cfc for the band given by index.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_exact_part(index, Optical, cfc_tf, klevel)
!  </TEMPLATE>
!  <IN NAME="index" TYPE="integer">
!   The spectral index where exact CTS transmision function is computed
!  </IN>
!  <IN NAME="Optical" TYPE="optical_depth_type">
!   The CFC gas optical depth
!  </IN>
!  <OUT NAME="cfc_tf" TYPE="real">
!   exact CTS transmission function output
!  </OUT>
!  <IN NAME="klevel" TYPE="integer">
!   The level below which exact CTS transmision function is computed
!  </IN>
! </SUBROUTINE>
!
subroutine cfc_exact_part (index, Optical, cfc_tf, klevel)

!----------------------------------------------------------------------
!    cfc_exact computes exact cool-to-space transmission function 
!    at levels below klevel for cfc for the band given by index. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     index
!     klevel
!     Optical
!
!  intent(out) variables:
!
!     cfc_tf
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer     ::  k     ! do-loop index
      integer     ::  kx    !

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod',   &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        cfc_tf(:,:,k) = 1.0E+00 - strf113(index)*    &
                        (Optical%totf113(:,:,k+1) -  &
                         Optical%totf113(:,:,klevel)) -   &
                        strf22 (index)*   &
                        (Optical%totf22(:,:,k+1) -   &
                         Optical%totf22(:,:,klevel)) -   &
                        strf11 (index)*                          &
                        (Optical%totf11(:,:,k+1) -   &
                         Optical%totf11(:,:,klevel)) -   &
                        strf12 (index)*      &   
                        (Optical%totf12(:,:,k+1) -   &
                         Optical%totf12(:,:,klevel)) 
      end do

!--------------------------------------------------------------------

end subroutine cfc_exact_part



!####################################################################
! <SUBROUTINE NAME="cfc_indx8">
!  <OVERVIEW>
!   cfc_indx8 computes transmission function for cfc for the band 8. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_indx8 computes transmission function for cfc for the band 8. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_indx8(index, Optical, tcfc8)
!  </TEMPLATE>
!  <IN NAME="index" TYPE="integer">
!   The spectral index where exact CTS transmision function is computed
!  </IN>
!  <IN NAME="Optical" TYPE="optical_depth_type">
!   The CFC gas optical depth
!  </IN>
!  <OUT NAME="tcfc8" TYPE="real">
!   exact CTS transmission function output for the band 8
!  </OUT>
! </SUBROUTINE>
!
subroutine cfc_indx8 (index, Optical, tcfc8)

!----------------------------------------------------------------------
!     cfc_indx8 computes transmission function for cfc for the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     index
!     Optical
!
!  intent(out) variables:
!
!     tcfc8 
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      tcfc8 (:,:,:) = 1.0E+00 -    &
                      strf113(index)*Optical%totf113 -   &
                      strf22 (index)*Optical%totf22  

!---------------------------------------------------------------------

end subroutine cfc_indx8



!####################################################################
! <SUBROUTINE NAME="cfc_indx8_part">
!  <OVERVIEW>
!   cfc_indx8 computes transmission function for cfc for the band 8. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_indx8 computes transmission function for cfc for the band 8. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_indx8_part(index, Optical, tcfc8, klevel)
!  </TEMPLATE>
!  <IN NAME="index" TYPE="integer">
!   The spectral index where exact CTS transmision function is computed
!  </IN>
!  <IN NAME="Optical" TYPE="optical_depth_type">
!   The CFC gas optical depth
!  </IN>
!  <OUT NAME="tcfc8" TYPE="real">
!   exact CTS transmission function output for the band 8
!  </OUT>
!  <IN NAME="klevel" TYPE="integer">
!   The level below which exact CTS transmision function is computed
!  </IN>
! </SUBROUTINE>
!
subroutine cfc_indx8_part (index, Optical, tcfc8, klevel)

!----------------------------------------------------------------------
!     cfc_indx8_part computes transmission function for cfc for 
!     the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     index
!     klevel
!     Optical
!
!  intent(out) variables:
!
!     tcfc8 
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer     :: kx
      integer     :: k

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        tcfc8 (:,:,k+1) = 1.0E+00 -  strf113(index)*    &
                          (Optical%totf113(:,:,k+1) -  &
                           Optical%totf113(:,:,klevel)) -   &
                          strf22 (index)*  &
                          (Optical%totf22(:,:,k+1) -   &
                           Optical%totf22(:,:,klevel)) 
      end do

!--------------------------------------------------------------------

end subroutine cfc_indx8_part




!####################################################################
! <SUBROUTINE NAME="cfc_overod">
!  <OVERVIEW>
!   cfc_overod computes transmission function for cfc that is used   
!   with overod variable in the 15 um (560-800 cm-1) band.
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_overod computes transmission function for cfc that is used   
!   with overod variable in the 15 um (560-800 cm-1) band.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_overod (Optical, cfc_tf)
!  </TEMPLATE>
!  <IN NAME="Optical" TYPE="optical_path_type">
!   CFC optical depth values
!  </IN>
!  <OUT NAME="cfc_tf" TYPE="real">
!   CFC transmission function
!  </OUT>
! </SUBROUTINE>
! 
subroutine cfc_overod (Optical, cfc_tf)

!----------------------------------------------------------------------
!     cfc_overod computes transmission function for cfc that is used   
!     with overod variable.
!----------------------------------------------------------------------

type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     Optical
!
!  intent(out) variables:
!
!     cfc_tf
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer :: kx

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kx = size (Optical%totf11,3) 
      cfc_tf(:,:,:) = 1.0E+00 -    &
                      sf11315*Optical%totf113(:,:,2:kx) -   &
                      sf2215 *Optical%totf22 (:,:,2:kx) -  &
                      sf1115*Optical%totf11  (:,:,2:kx) -  &
                      sf1215*Optical%totf12  (:,:,2:kx)

!---------------------------------------------------------------------


end subroutine cfc_overod




!####################################################################
! <SUBROUTINE NAME="cfc_overod_part">
!  <OVERVIEW>
!   cfc_overod computes transmission function for cfc that is used   
!   with overod variable in the 15 um (560-800 cm-1) band from klevel down.
!  </OVERVIEW>
!  <DESCRIPTION>
!   cfc_overod computes transmission function for cfc that is used   
!   with overod variable in the 15 um (560-800 cm-1) band from klevel down.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_overod_part (Optical, cfc_tf, klevel)
!  </TEMPLATE>
!  <IN NAME="Optical" TYPE="optical_path_type">
!   CFC optical depth values
!  </IN>
!  <OUT NAME="cfc_tf" TYPE="real">
!   CFC transmission function
!  </OUT>
!  <IN NAME="klevel" TYPE="integer">
!   The level below which exact CTS transmision function is computed
!  </IN>
! </SUBROUTINE>
!  
subroutine cfc_overod_part (Optical, cfc_tf, klevel)

!----------------------------------------------------------------------
!    cfc_overod_part computes transmission function for cfc that is 
!    used with overod variable from klevel down.
!----------------------------------------------------------------------

type(optical_path_type), intent(in)    :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf
integer,                 intent(in)    :: klevel

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     Optical
!     klevel
!
!  intent(out) variables:
!
!     cfc_tf
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer     ::      kx
      integer     ::      k

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        cfc_tf(:,:,k) = 1.0E+00 - sf11315*      &
                        (Optical%totf113(:,:,k+1) -   &
                         Optical%totf113(:,:,klevel)) - &
                        sf2215 *  &
                        (Optical%totf22 (:,:,k+1) -   &
                         Optical%totf22 (:,:,klevel)) -   & 
                        sf1115*   &
                        (Optical%totf11 (:,:,k+1) -   &
                         Optical%totf11 (:,:,klevel)) -   & 
                        sf1215*  &
                        (Optical%totf12 (:,:,k+1) -   &
                         Optical%totf12 (:,:,klevel))   
      end do

!---------------------------------------------------------------------


end subroutine cfc_overod_part


!#####################################################################
! <SUBROUTINE NAME="lw_gases_stdtf_end">
!  <OVERVIEW>
!   lw_gases_stdtf_end is the destructor for lw_gases_stdtf_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   lw_gases_stdtf_end is the destructor for lw_gases_stdtf_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_gases_stdtf_end           
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_gases_stdtf_end

!--------------------------------------------------------------------
!    lw_gases_stdtf_end is the destructor for lw_gases_stdtf_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('lw_gases_stdtf_mod', &
        !      'module has not been initialized', FATAL )
        write(6,*) "lw_gases_stdtf_mod has not been initialized"
        stop
      endif
 
!-------------------------------------------------------------------
!    deallocate module variables.
!--------------------------------------------------------------------
      deallocate (pd, plm, pd8, plm8)
      deallocate (pa)

      deallocate (xa                )
      deallocate (ca                )
      deallocate (uexp              )
      deallocate (sexp              )
      deallocate (press_lo          )
      deallocate (press_hi          )

      deallocate (  trns_interp_lyr_ps )
      deallocate (  trns_interp_lyr_ps8 )
      deallocate (  trns_interp_lvl_ps )
      deallocate (  trns_interp_lvl_ps8 )
 
      deallocate (  trns_interp_lyr_ps_nf )
      deallocate (  trns_interp_lyr_ps8_nf )
      deallocate (  trns_interp_lvl_ps_nf )
      deallocate (  trns_interp_lvl_ps8_nf )

      deallocate (  dgasdt8_lvl )
      deallocate (  dgasdt10_lvl )
      deallocate (  d2gast8_lvl  )
      deallocate (  d2gast10_lvl )
      deallocate (  gasp10_lvl )
      deallocate (  gasp8_lvl  )
      deallocate (  dgasdt8_lvlcts  )
      deallocate (  dgasdt10_lvlcts )
      deallocate (  d2gast8_lvlcts  )
      deallocate (  d2gast10_lvlcts )
      deallocate (  gasp10_lvlcts )
      deallocate (  gasp8_lvlcts  )
      deallocate (  dgasdt8_lyr  )
      deallocate (  dgasdt10_lyr )
      deallocate (  d2gast8_lyr  )
      deallocate (  d2gast10_lyr )
      deallocate (  gasp10_lyr )
      deallocate (  gasp8_lyr  )

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

end subroutine lw_gases_stdtf_end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!######################################################################
! <SUBROUTINE NAME="std_lblpressures">
!  <OVERVIEW>
!   calculation of pa -- the "table" of (NSTDCO2LVLS) grid pressures
!  </OVERVIEW>
!  <DESCRIPTION>
!   calculation of pa -- the "table" of (NSTDCO2LVLS) grid pressures
!  </DESCRIPTION>
!  <TEMPLATE>
!   call std_lblpressures
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine std_lblpressures
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
      real      ::  fact15, fact30

!---------------------------------------------------------------------
!  local variables
!
!     fact15
!
!-------------------------------------------------------------------
     integer :: k
!---------------------------------------------------------------------
!    calculation of pa -- the "table" of (NSTDCO2LVLS) grid pressures
!    note-this code must not be changed by the user!!!!!!!!!
!---------------------------------------------------------------------
      allocate ( pa (NSTDCO2LVLS) )

      if (NSTDCO2LVLS .EQ. 109) then
        pa(1)=0.
        fact15=10.**(1./15.)
        fact30=10.**(1./30.)
        pa(2)=1.0e-3
        do k=2,76
          pa(k+1)=pa(k)*fact15
        enddo
        do k=77,108
          pa(k+1)=pa(k)*fact30
        enddo
      else if (NSTDCO2LVLS .EQ. 496) then
        pa=(/ 0.0000000000E+00, 0.1000000000E-05, 0.1079775162E-05, 0.1165914401E-05, 0.1258925412E-05, &
              0.1359356391E-05, 0.1467799268E-05, 0.1584893192E-05, 0.1711328304E-05, 0.1847849797E-05, &
              0.1995262315E-05, 0.2154434690E-05, 0.2326305067E-05, 0.2511886432E-05, 0.2712272579E-05, &
              0.2928644565E-05, 0.3162277660E-05, 0.3414548874E-05, 0.3686945065E-05, 0.3981071706E-05, &
              0.4298662347E-05, 0.4641588834E-05, 0.5011872336E-05, 0.5411695265E-05, 0.5843414134E-05, &
              0.6309573445E-05, 0.6812920691E-05, 0.7356422545E-05, 0.7943282347E-05, 0.8576958986E-05, &
              0.9261187281E-05, 0.1000000000E-04, 0.1079775162E-04, 0.1165914401E-04, 0.1258925412E-04, &
              0.1359356391E-04, 0.1467799268E-04, 0.1584893192E-04, 0.1711328304E-04, 0.1847849797E-04, &
              0.1995262315E-04, 0.2154434690E-04, 0.2326305067E-04, 0.2511886432E-04, 0.2712272579E-04, &
              0.2928644565E-04, 0.3162277660E-04, 0.3414548874E-04, 0.3686945065E-04, 0.3981071706E-04, &
              0.4298662347E-04, 0.4641588834E-04, 0.5011872336E-04, 0.5411695265E-04, 0.5843414134E-04, &
              0.6309573445E-04, 0.6812920691E-04, 0.7356422545E-04, 0.7943282347E-04, 0.8576958986E-04, &
              0.9261187281E-04, 0.1000000000E-03, 0.1079775162E-03, 0.1165914401E-03, 0.1258925412E-03, &
              0.1359356391E-03, 0.1467799268E-03, 0.1584893192E-03, 0.1711328304E-03, 0.1847849797E-03, &
              0.1995262315E-03, 0.2154434690E-03, 0.2326305067E-03, 0.2511886432E-03, 0.2712272579E-03, &
              0.2928644565E-03, 0.3162277660E-03, 0.3414548874E-03, 0.3686945065E-03, 0.3981071706E-03, &
              0.4298662347E-03, 0.4641588834E-03, 0.5011872336E-03, 0.5411695265E-03, 0.5843414134E-03, &
              0.6309573445E-03, 0.6812920691E-03, 0.7356422545E-03, 0.7943282347E-03, 0.8576958986E-03, &
              0.9261187281E-03, 0.1000000000E-02, 0.1079775162E-02, 0.1165914401E-02, 0.1258925412E-02, &
              0.1359356391E-02, 0.1467799268E-02, 0.1584893192E-02, 0.1711328304E-02, 0.1847849797E-02, &
              0.1995262315E-02, 0.2154434690E-02, 0.2326305067E-02, 0.2511886432E-02, 0.2712272579E-02, &
              0.2928644565E-02, 0.3162277660E-02, 0.3414548874E-02, 0.3686945065E-02, 0.3981071706E-02, &
              0.4298662347E-02, 0.4641588834E-02, 0.5011872336E-02, 0.5411695265E-02, 0.5843414134E-02, &
              0.6309573445E-02, 0.6812920691E-02, 0.7356422545E-02, 0.7943282347E-02, 0.8576958986E-02, &
              0.9261187281E-02, 0.1000000000E-01, 0.1079775162E-01, 0.1165914401E-01, 0.1258925412E-01, &
              0.1359356391E-01, 0.1467799268E-01, 0.1584893192E-01, 0.1711328304E-01, 0.1847849797E-01, &
              0.1995262315E-01, 0.2154434690E-01, 0.2326305067E-01, 0.2511886432E-01, 0.2712272579E-01, &
              0.2928644565E-01, 0.3162277660E-01, 0.3414548874E-01, 0.3686945065E-01, 0.3981071706E-01, &
              0.4298662347E-01, 0.4641588834E-01, 0.5011872336E-01, 0.5411695265E-01, 0.5843414134E-01, &
              0.6309573445E-01, 0.6812920691E-01, 0.7356422545E-01, 0.7943282347E-01, 0.8576958986E-01, &
              0.9261187281E-01, 0.1000000000E+00, 0.1079775162E+00, 0.1165914401E+00, 0.1258925412E+00, &
              0.1359356391E+00, 0.1467799268E+00, 0.1584893192E+00, 0.1711328304E+00, 0.1847849797E+00, &
              0.1995262315E+00, 0.2154434690E+00, 0.2326305067E+00, 0.2511886432E+00, 0.2712272579E+00, &
              0.2928644565E+00, 0.3162277660E+00, 0.3414548874E+00, 0.3686945065E+00, 0.3981071706E+00, &
              0.4298662347E+00, 0.4641588834E+00, 0.5011872336E+00, 0.5411695265E+00, 0.5843414134E+00, &
              0.6309573445E+00, 0.6812920691E+00, 0.7356422545E+00, 0.7943282347E+00, 0.8576958986E+00, &
              0.9261187281E+00, 0.1000000000E+01, 0.1079775162E+01, 0.1165914401E+01, 0.1258925412E+01, &
              0.1359356391E+01, 0.1467799268E+01, 0.1584893192E+01, 0.1711328304E+01, 0.1847849797E+01, &
              0.1995262315E+01, 0.2154434690E+01, 0.2326305067E+01, 0.2511886432E+01, 0.2712272579E+01, &
              0.2928644565E+01, 0.3162277660E+01, 0.3414548874E+01, 0.3686945065E+01, 0.3981071706E+01, &
              0.4298662347E+01, 0.4641588834E+01, 0.5011872336E+01, 0.5411695265E+01, 0.5843414134E+01, &
              0.6309573445E+01, 0.6812920691E+01, 0.7356422545E+01, 0.7943282347E+01, 0.8576958986E+01, &
              0.9261187281E+01, 0.1000000000E+02, 0.1079775162E+02, 0.1165914401E+02, 0.1258925412E+02, &
              0.1359356391E+02, 0.1467799268E+02, 0.1584893192E+02, 0.1711328304E+02, 0.1847849797E+02, &
              0.1995262315E+02, 0.2154434690E+02, 0.2326305067E+02, 0.2511886432E+02, 0.2712272579E+02, &
              0.2928644565E+02, 0.3162277660E+02, 0.3414548874E+02, 0.3686945065E+02, 0.3981071706E+02, &
              0.4298662347E+02, 0.4641588834E+02, 0.5011872336E+02, 0.5411695265E+02, 0.5843414134E+02, &
              0.6309573445E+02, 0.6812920691E+02, 0.7356422545E+02, 0.7943282347E+02, 0.8576958986E+02, &
              0.9261187281E+02, 0.1000000000E+03, 0.1075000000E+03, 0.1150000000E+03, 0.1225000000E+03, &
              0.1300000000E+03, 0.1375000000E+03, 0.1450000000E+03, 0.1525000000E+03, 0.1600000000E+03, &
              0.1675000000E+03, 0.1750000000E+03, 0.1825000000E+03, 0.1900000000E+03, 0.1975000000E+03, &
              0.2050000000E+03, 0.2125000000E+03, 0.2200000000E+03, 0.2275000000E+03, 0.2350000000E+03, &
              0.2425000000E+03, 0.2500000000E+03, 0.2575000000E+03, 0.2650000000E+03, 0.2725000000E+03, &
              0.2800000000E+03, 0.2875000000E+03, 0.2950000000E+03, 0.3025000000E+03, 0.3100000000E+03, &
              0.3175000000E+03, 0.3250000000E+03, 0.3325000000E+03, 0.3398345186E+03, 0.3470140487E+03, &
              0.3540480193E+03, 0.3609449404E+03, 0.3677125236E+03, 0.3743577834E+03, 0.3808871224E+03, &
              0.3873064033E+03, 0.3936210106E+03, 0.3998359038E+03, 0.4059556626E+03, 0.4119845264E+03, &
              0.4179264289E+03, 0.4237850281E+03, 0.4295637322E+03, 0.4352657234E+03, 0.4408939782E+03, &
              0.4464512851E+03, 0.4519402615E+03, 0.4573633676E+03, 0.4627229193E+03, 0.4680211000E+03, &
              0.4732599708E+03, 0.4784414802E+03, 0.4835674720E+03, 0.4886396934E+03, 0.4936598019E+03, &
              0.4986293714E+03, 0.5035498982E+03, 0.5084228063E+03, 0.5132494520E+03, 0.5180311284E+03, &
              0.5227690695E+03, 0.5274644538E+03, 0.5321184079E+03, 0.5367320095E+03, 0.5413062904E+03, &
              0.5458422391E+03, 0.5503408035E+03, 0.5548028929E+03, 0.5592293805E+03, 0.5636211050E+03, &
              0.5679788728E+03, 0.5723034597E+03, 0.5765956122E+03, 0.5808560493E+03, 0.5850854638E+03, &
              0.5892845238E+03, 0.5934538735E+03, 0.5975941348E+03, 0.6017059082E+03, 0.6057897738E+03, &
              0.6098462921E+03, 0.6138760054E+03, 0.6178794381E+03, 0.6218570977E+03, 0.6258094758E+03, &
              0.6297370483E+03, 0.6336402765E+03, 0.6375196075E+03, 0.6413754751E+03, 0.6452082997E+03, &
              0.6490184897E+03, 0.6528064415E+03, 0.6565725398E+03, 0.6603171586E+03, 0.6640406614E+03, &
              0.6677434013E+03, 0.6714257219E+03, 0.6750879572E+03, 0.6787304325E+03, 0.6823534641E+03, &
              0.6859573602E+03, 0.6895424207E+03, 0.6931089380E+03, 0.6966571969E+03, 0.7001874749E+03, &
              0.7037000426E+03, 0.7071951640E+03, 0.7106730964E+03, 0.7141340911E+03, 0.7175783929E+03, &
              0.7210062413E+03, 0.7244178697E+03, 0.7278135063E+03, 0.7311933739E+03, 0.7345576900E+03, &
              0.7379066675E+03, 0.7412405143E+03, 0.7445594335E+03, 0.7478636239E+03, 0.7511532800E+03, &
              0.7544285917E+03, 0.7576897452E+03, 0.7609369225E+03, 0.7641703017E+03, 0.7673900573E+03, &
              0.7705963600E+03, 0.7737893770E+03, 0.7769692722E+03, 0.7801362061E+03, 0.7832903357E+03, &
              0.7864318152E+03, 0.7895607956E+03, 0.7926774249E+03, 0.7957818482E+03, 0.7988742079E+03, &
              0.8019546434E+03, 0.8050232916E+03, 0.8080802869E+03, 0.8111257609E+03, 0.8141598430E+03, &
              0.8171826601E+03, 0.8201943367E+03, 0.8231949951E+03, 0.8261847554E+03, 0.8291637353E+03, &
              0.8321320508E+03, 0.8350898155E+03, 0.8380371412E+03, 0.8409741375E+03, 0.8439009124E+03, &
              0.8468175719E+03, 0.8497242200E+03, 0.8526209592E+03, 0.8555078901E+03, 0.8583851117E+03, &
              0.8612527213E+03, 0.8641108147E+03, 0.8669594858E+03, 0.8697988273E+03, 0.8726289303E+03, &
              0.8754498843E+03, 0.8782617776E+03, 0.8810646968E+03, 0.8838587274E+03, 0.8866439533E+03, &
              0.8894204574E+03, 0.8921883209E+03, 0.8949476242E+03, 0.8976984460E+03, 0.9004408642E+03, &
              0.9031749554E+03, 0.9059007948E+03, 0.9086184568E+03, 0.9113280145E+03, 0.9140295400E+03, &
              0.9167231043E+03, 0.9194087774E+03, 0.9220866283E+03, 0.9247567248E+03, 0.9274191339E+03, &
              0.9300739218E+03, 0.9327211534E+03, 0.9353608929E+03, 0.9379932036E+03, 0.9406181478E+03, &
              0.9432357871E+03, 0.9458461820E+03, 0.9484493924E+03, 0.9510454774E+03, 0.9536344950E+03, &
              0.9562165027E+03, 0.9587915571E+03, 0.9613597142E+03, 0.9639210289E+03, 0.9664755558E+03, &
              0.9690233485E+03, 0.9715644600E+03, 0.9740989426E+03, 0.9766268479E+03, 0.9791482268E+03, &
              0.9816631296E+03, 0.9841716060E+03, 0.9866737049E+03, 0.9891694749E+03, 0.9916589636E+03, &
              0.9941422182E+03, 0.9966192854E+03, 0.9990902111E+03, 0.1001555041E+04, 0.1004013820E+04, &
              0.1006466592E+04, 0.1008913401E+04, 0.1011354290E+04, 0.1013789303E+04, 0.1016218480E+04, &
              0.1018641865E+04, 0.1021059499E+04, 0.1023471421E+04, 0.1025877673E+04, 0.1028278294E+04, &
              0.1030673324E+04, 0.1033062801E+04, 0.1035446763E+04, 0.1037825250E+04, 0.1040198298E+04, &
              0.1042565945E+04, 0.1044928227E+04, 0.1047285181E+04, 0.1049636842E+04, 0.1051983246E+04, &
              0.1054324428E+04, 0.1056660423E+04, 0.1058991265E+04, 0.1061316988E+04, 0.1063637626E+04, &
              0.1065953212E+04, 0.1068263778E+04, 0.1070569358E+04, 0.1072869983E+04, 0.1075165685E+04, &
              0.1077456496E+04, 0.1079742446E+04, 0.1082023567E+04, 0.1084299889E+04, 0.1086571443E+04, &
              0.1088838257E+04, 0.1091100362E+04, 0.1093357787E+04, 0.1095610560E+04, 0.1097858711E+04, &
              0.1100102268E+04 /)
      endif

!---------------------------------------------------------------------
 
end subroutine std_lblpressures




!#####################################################################
! <SUBROUTINE NAME="approx_fn">
!  <OVERVIEW>
!   Subroutine to compute co2 approximation function
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute co2 approximation function
!  </DESCRIPTION>
!  <TEMPLATE>
!   call approx_fn (press_hi_app, press_lo_app, do_triangle,     &  
!                   nklo, nkhi, nkplo, nkphi,                      &
!                   ca_app, sexp_app, xa_app, uexp_app, approx)
!  </TEMPLATE>
!  <IN NAME="press_hi_app" TYPE="real">
!   high standard pressure array
!  </IN>
!  <IN NAME="press_lo_app" TYPE="real">
!   low standard pressure array
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   state variable of interpolation scheme
!  </IN>
!  <IN NAME="nklo, nkhi" TYPE="integer">
!   vertical level pairs: the upper and lower level index
!  </IN>
!  <IN NAME="nkplo, nkphi" TYPE="integer">
!   pressure level pairs: the upper and lower level index
!  </IN>
!  <IN NAME="ca_app, sexp_app, xa_app, uexp_app" TYPE="real">
!   The interpolation coefficients
!  </IN>
!  <OUT NAME="approx" TYPE="real">
!   co2 approximation function
!  </OUT>
! </SUBROUTINE>
!
subroutine approx_fn (press_hi_app, press_lo_app, do_triangle,     &  
                      nklo, nkhi, nkplo, nkphi,                      &
                      ca_app, sexp_app, xa_app, uexp_app, approx)

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    approx_fn computes the co2 approximation function
!          A(press_hi_app(i), press_lo_app(j))  (Eq.(4), Ref. (2))
!    for a particular co2 amount and a standard pressure grid (pa).
!    the calculation is performed for all (press_hi(k),press_lo(k')
!    pairs which are possible according to the state of (do_triangle).
!    the path function (upathv) is evaluated using the expression
!    given in Eqs. (5) and (A5) in Ref. (2) for the co2 interpolation
!    program between a lower model pressure (press_lo_app) and
!    a higher model pressure (press_hi_app) using the interpolation
!    coefficients (ca_app, sexp_app, xa_app, uexp_app) computed in
!    subroutine coeint.
!         the output is in (approx).
!----------------------------------------------------------------------

real, dimension(:,:), intent(in)  :: press_hi_app, press_lo_app, &
                                     ca_app, sexp_app, xa_app,   &
                                     uexp_app
logical,              intent(in)  :: do_triangle 
integer,              intent(in)  :: nklo, nkhi, nkplo, nkphi
real, dimension(:,:), intent(out) :: approx

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     press_hi_app
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables

       real, dimension(size(press_hi_app,1),size(press_hi_app,2)) :: upathv

       integer   :: k, kp, kp0
       integer   :: k1, k2

!----------------------------------------------------------------------
!    obtain array extents for internal arrays  and allocate these arrays
!----------------------------------------------------------------------
      k1 = size(press_hi_app,1)       ! this corresponds to ndimkp
      k2 = size(press_hi_app,2)       ! this corresponds to ndimk
 
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi

!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(alog(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!          upathv(kp,k) = (press_hi_app(kp,k) -    &
!   press_lo_app(kp,k))**(1./sexp_app(kp,k))* &
!                         (press_hi_app(kp,k) + press_lo_app(kp,k) +  &
!          upathv(kp,k) = upathv(kp,k)**uexp_app(kp,k)
!          approx(kp,k) = (ca_app(kp,k)*           &
!                          LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))**   &
!                          (sexp_app(kp,k)/uexp_app(kp,k))
!------------------------------------------------------------------

          upathv(kp,k) = (press_hi_app(kp,k) -    &
                          press_lo_app(kp,k))**(1./sexp_app(kp,k))* &
                          (press_hi_app(kp,k) + press_lo_app(kp,k) +  &
                          dop_core)
           upathv(kp,k) = upathv(kp,k)**uexp_app(kp,k)
           approx(kp,k) = (ca_app(kp,k)*           &
                           LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))**   &
                           (sexp_app(kp,k)/uexp_app(kp,k))

!          upathv(kp,k) = EXP((1.0/sexp_app(kp,k))* ALOG(  &
!  (press_hi_app(kp,k) - press_lo_app(kp,k))))* &
!                         (press_hi_app(kp,k) + press_lo_app(kp,k) +  &
!                                     dop_core)
!          upathv(kp,k) = EXP(uexp_app(kp,k)*ALOG(upathv(kp,k)))
!          approx(kp,k) = EXP( ((sexp_app(kp,k)/uexp_app(kp,k)) *   &
!   ALOG( ca_app(kp,k)*  &
!                          LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))))
          enddo
        enddo

end subroutine approx_fn



!#####################################################################
! <SUBROUTINE NAME="approx_fn_std">
!  <OVERVIEW>
!   Subroutine to compute co2 approximation function
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute co2 approximation function
!  </DESCRIPTION>
!  <TEMPLATE>
!   call approx_fn_std (press_hi_app, press_lo_app, do_triangle,     &  
!                   ca_app, sexp_app, xa_app, uexp_app, approx)
!  </TEMPLATE>
!  <IN NAME="press_hi_app" TYPE="real">
!   high standard pressure array
!  </IN>
!  <IN NAME="press_lo_app" TYPE="real">
!   low standard pressure array
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   state variable of interpolation scheme
!  </IN>
!  <IN NAME="ca_app, sexp_app, xa_app, uexp_app" TYPE="real">
!   The interpolation coefficients
!  </IN>
!  <OUT NAME="approx" TYPE="real">
!   co2 approximation function
!  </OUT>
! </SUBROUTINE>
!
subroutine approx_fn_std (press_hi_pa, press_lo_pa, do_triangle, &
                          ca_app, sexp_app, xa_app, uexp_app,  &
                          approx)
 
!---------------------------------------------------------------------
!    approx_fn_std computes the co2 approximation function
!               A(press_hi(i), press_lo(j))  (Eq.(4), Ref. (2))
!    for a particular co2 amount and a standard pressure grid (pa).
!    the calculation is performed for all (press_hi(k),press_lo(k')
!    pairs which are possible according to the state of (do_triangle).
!    the path function (upathv) is evaluated using the expression
!    given in Eqs. (5) and (A5) in Ref. (2) for the co2 interpolation
!    program between a lower standard pressure (press_lo_pa) and
!    a higher standard pressure (press_hi_pa) using the interpolation
!    coefficients (ca_app, sexp_app, xa_app, uexp_app) computed in
!    subroutine coeint.
!         the output is in (approx).
!--------------------------------------------------------------------

real, dimension(:,:), intent(in)    :: ca_app, sexp_app, xa_app, &
                                       uexp_app
real, dimension(:,:), intent(in)    :: press_hi_pa, press_lo_pa
logical,              intent(in)    :: do_triangle
real, dimension(:,:), intent(out)   :: approx

!--------------------------------------------------------------------
!  local variables
      
      real, dimension(NSTDCO2LVLS,NSTDCO2LVLS) :: upathv
      integer         :: k, kp, kp0, NSTDCO2LVLS_loop

!--------------------------------------------------------------------
!  local variables
!   
!     upathv
!
!--------------------------------------------------------------------
 
      if (do_co2_bug) then
        NSTDCO2LVLS_loop = NSTDCO2LVLS
      else 
        NSTDCO2LVLS_loop = NSTDCO2LVLS -1 
      end if 

      do k=1,NSTDCO2LVLS_loop  
        if (do_triangle) then
          kp0 = k + 1
        else
          kp0 = 1
        endif
        do kp=kp0,NSTDCO2LVLS

!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(alog(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in 
!    radiag file
!         upathv(kp,k) = (press_hi_pa(kp,k) -    &
!  press_lo_pa(kp,k))**(1./sexp_app(kp,k)) *  &
!                         (press_hi_pa(kp,k) + press_lo_pa(kp,k) +    &
!    dop_core)
!         upathv(kp,k) = upathv(kp,k)**uexp_app(kp,k)
!         approx(kp,k) = (ca_app(kp,k)*                              &
!                         LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))**    &
!                             (sexp_app(kp,k)/uexp_app(kp,k))
!---------------------------------------------------------------------
            
          upathv(kp,k) = EXP((1.0/sexp_app(kp,k))* log(  &
                        (press_hi_pa(kp,k) - press_lo_pa(kp,k))))* &
                        (press_hi_pa(kp,k) + press_lo_pa(kp,k) +  &
                         dop_core)
          upathv(kp,k) = EXP(uexp_app(kp,k)*log(upathv(kp,k)))
          approx(kp,k) = EXP( ((sexp_app(kp,k)/uexp_app(kp,k)) *  &
                         log( ca_app(kp,k)*   &
                         LOG(1.0 + xa_app(kp,k)*upathv(kp,k)))))
        enddo
      enddo

!--------------------------------------------------------------------


end subroutine approx_fn_std


!#####################################################################
! <SUBROUTINE NAME="gasins">
!  <OVERVIEW>
!   gasins processes transmission functions to produce 
!     "consolidated" functions over the specific frequency band
!     ranges needed by the SEA code, and the derivatives needed
!     by the SEA algorithm.
!  </OVERVIEW>
!  <DESCRIPTION>
!    gasins processes transmission functions to produce 
!     "consolidated" functions over the specific frequency band
!     ranges needed by the SEA code, and the derivatives needed
!     by the SEA algorithm.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call gasins (gas_type, do_lvlcalc, do_lvlctscalc, do_lyrcalc,  &
!                   nf, ntbnd, ndimkp, ndimk,  &
!                   dgasdt10_lvl, dgasdt10_lvlcts, dgasdt10_lyr,  &
!                   gasp10_lvl, gasp10_lvlcts, gasp10_lyr,  &
!                   d2gast10_lvl, d2gast10_lvlcts, d2gast10_lyr,  &
!                   dgasdt8_lvl,  dgasdt8_lvlcts,  dgasdt8_lyr ,  &
!                   gasp8_lvl,  gasp8_lvlcts,  gasp8_lyr ,  &
!                   d2gast8_lvl,  d2gast8_lvlcts,  d2gast8_lyr )
!  </TEMPLATE>
!  <IN NAME="gas_type" TYPE="character">
!   Gas type information
!  </IN>
!  <IN NAME="do_lvlcalc, do_lvlctscalc, do_lyrcalc" TYPE="logical">
!   State variables that determines calculation paths
!        do_lvlcalc     : compute level co2 transmissivities if true.
!        do_lyrcalc     : compute layer co2 transmissivities if true.
!        do_lvlctscalc  : compute cts level co2 transmissivities if true
!  </IN>
!  <IN NAME="nf" TYPE="integer">
!   frequency band number
!  </IN>
!  <IN NAME="ntbnd" TYPE="integer">
!   temperature index of the frequency band
!  </IN>
!  <IN NAME="ndimkp, ndimk" TYPE="integer">
!   extents of dimensions for output interpolation transmissivity arrays.
!  </IN>
!  <OUT NAME="dgasdt10_lvl gasp10_lvl d2gast10_lvl dgasdt8_lvl gasp8_lvl d2gast8_lvl" TYPE="real">
!   variables used in do_lvlcalc calculation path
!  </OUT>
!  <OUT NAME="dgasdt10_lvlcts gasp10_lvlcts d2gast10_lvlcts dgasdt8_lvlcts gasp8_lvlcts d2gast8_lvlcts" TYPE="real">
!   variables used in do_lvlctscalc calculation path
!  </OUT>
!  <OUT NAME="dgasdt10_lyr gasp10_lyr d2gast10_lyr dgasdt8_lyr gasp8_lyr d2gast8_lyr" TYPE="real">
!   variables used in do_lyrcalc calculation path
!  </OUT>
! </SUBROUTINE>
!

subroutine gasins (gas_type, do_lvlcalc, do_lvlctscalc, do_lyrcalc,  &
                   nf, ntbnd, ndimkp, ndimk,  &
                   dgasdt10_lvl, dgasdt10_lvlcts, dgasdt10_lyr,  &
                   gasp10_lvl, gasp10_lvlcts, gasp10_lyr,  &
                   d2gast10_lvl, d2gast10_lvlcts, d2gast10_lyr,  &
                   dgasdt8_lvl,  dgasdt8_lvlcts,  dgasdt8_lyr ,  &
                   gasp8_lvl,  gasp8_lvlcts,  gasp8_lyr ,  &
                   d2gast8_lvl,  d2gast8_lvlcts,  d2gast8_lyr )
 
!-------------------------------------------------------------------
!    gasins processes transmission functions to produce 
!    "consolidated" functions over the specific frequency band
!    ranges needed by the SEA code, and the derivatives needed
!    by the SEA algorithm. writing to a file, formerly done in
!    this module, is now done (if needed) in write_seaco2fcns.F
!-------------------------------------------------------------------

character(len=*),     intent(in)   :: gas_type
logical,              intent(in)   :: do_lvlcalc, do_lyrcalc, &
                                      do_lvlctscalc
integer,              intent(in)   :: nf, ntbnd
integer,              intent(in)   :: ndimkp, ndimk
real, dimension(:,:), intent(out)  :: &
                                      dgasdt10_lvl, dgasdt10_lyr,  &
                                      gasp10_lvl, gasp10_lyr,  &
                                      d2gast10_lvl, d2gast10_lyr,  &
                                      dgasdt8_lvl,  dgasdt8_lyr ,  &
                                      gasp8_lvl,  gasp8_lyr ,  &
                                      d2gast8_lvl,  d2gast8_lyr
real, dimension(:)  , intent(out)  :: &
                                      dgasdt10_lvlcts, &
                                      gasp10_lvlcts, &
                                      d2gast10_lvlcts, &
                                      dgasdt8_lvlcts,  &
                                      gasp8_lvlcts,  &
                                      d2gast8_lvlcts   

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    gas_type 
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables

      integer      :: k1, k2, k, kp
      real         :: c1,c2

!-------------------------------------------------------------------
!  local variables
!
!     k1
!
!-------------------------------------------------------------------

!--------------------------------------------------------------------
!    obtain array extents for internal arrays and allocate these arrays
!--------------------------------------------------------------------
      k1 = size(trns_interp_lvl_ps_nf,1) ! this corresponds to ndimkp
      k2 = size(trns_interp_lvl_ps_nf,2) ! this corresponds to ndimk
 
!---------------------------------------------------------------------
!    the following code is rewritten so that the radiative bands
!    are: 
!        nf=1    560-800     (consol.=490-850)
!        nf=2    560-630      consol=490-630
!        nf=3    630-700      consol=630-700
!        nf=4    700-800      consol=700-850
!        nf=5   2270-2380     consol=2270-2380
!        nf=6    990-1070     consol=990-1070
!        nf=7    900-990     consol=900-990
!        nf=8    1070-1200     consol=1070-1200
!    the following loop obtains transmission functions for bands
!    used in radiative model calculations,with the equivalent
!    widths kept from the original consolidated co2 tf's.
!---------------------------------------------------------------------
      if (gas_type .EQ. 'co2') then
        if (nf.eq.1) then
          c1=1.5
          c2=0.5
        else if (nf.eq.2) then
          c1=2.0
          c2=1.0
        else if (nf.eq.3) then
          c1=1.0
          c2=0.0
        else if (nf.eq.4) then
          c1=1.5
          c2=0.5
        else if (nf.eq.5) then
          c1=1.0
          c2=0.0
        else if (nf.eq.6) then
          c1=1.0
          c2=0.0
        else if (nf.eq.7) then
          c1=1.0
          c2=0.0
        else if (nf.eq.8) then
          c1=1.0
          c2=0.0
        else
          write(6,*) 'illegal value of nf for co2'
          stop
        endif

!--------------------------------------------------------------------
!    the following code is rewritten so that the radiative bands
!    are: 
!        nf=1    1200-1400    consol=1200-1400
!    the following loop obtains transmission functions for bands
!    used in radiative model calculations,with the equivalent
!    widths kept from the original consolidated ch4 tf's.
!--------------------------------------------------------------------
      else if (gas_type .EQ. 'ch4') then
        if (nf.eq.1) then
          c1=1.0
          c2=0.0
        else 
          write(6,*) 'illegal value of nf for ch4'
          stop
        endif

!--------------------------------------------------------------------
!    the following code is rewritten so that the radiative bands
!    are: 
!        nf=1    1200-1400    consol=1200-1400
!        nf=2    1070-1200    consol=1070-1200
!        nf=3    560-630    consol=560-630
!    the following loop obtains transmission functions for bands
!    used in radiative model calculations,with the equivalent
!    widths kept from the original consolidated n2o tf's.
!--------------------------------------------------------------------
      else if (gas_type .EQ. 'n2o') then
        if (nf.eq.1) then
          c1=1.0
          c2=0.0
        else if (nf.eq.2) then
          c1=1.0
          c2=0.0
        else if (nf.eq.3) then
          c1=1.0
          c2=0.0
        else
          write(6,*) 'illegal value of nf for n2o'
          stop
        endif
      else 
        write(6,*) 'radiative gas type unrecognized in lw_gases_stdtf'
        stop
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_lvlcalc) then
        do k=1,k2
          do kp=1,k1
            gasp10_lvl(kp,k) = c1*trns_interp_lvl_ps_nf(kp,k,1) - c2
            gasp8_lvl(kp,k) = c1*trns_interp_lvl_ps8_nf(kp,k,1) - c2
          enddo
        enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (ntbnd .EQ. 3) then
          do k=1,k2
            do kp=1,k1
              dgasdt10_lvl(kp,k) = .02*(trns_interp_lvl_ps_nf(kp,k,2) -&
                                        trns_interp_lvl_ps_nf(kp,k,3))*&
                                   100.
              dgasdt8_lvl(kp,k) = .02*(trns_interp_lvl_ps8_nf(kp,k,2) -&
                                       trns_interp_lvl_ps8_nf(kp,k,3))*&
                                   100.
              d2gast10_lvl(kp,k) = .0016*    &
                                   (trns_interp_lvl_ps_nf(kp,k,2) +&
                                    trns_interp_lvl_ps_nf(kp,k,3) -&
                                    2.0*   &
                                    trns_interp_lvl_ps_nf(kp,k,1))*&
                                    1000.
              d2gast8_lvl(kp,k)  =   &
                                    .0016*  &
                                    (trns_interp_lvl_ps8_nf(kp,k,2) +&
                                     trns_interp_lvl_ps8_nf(kp,k,3) - &
                                     2.0* &
                                     trns_interp_lvl_ps_nf(kp,k,1))*  &
                                     1000.
            enddo
          enddo
        endif
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_lvlctscalc) then
        do kp=1,k1
          gasp10_lvlcts(kp) = c1*trns_interp_lvl_ps_nf(kp,1,1) - c2
          gasp8_lvlcts(kp) = c1*trns_interp_lvl_ps8_nf(kp,1,1) - c2
        enddo
        if (ntbnd .EQ. 3) then
          do kp=1,k1
            dgasdt10_lvlcts(kp) = .02*(trns_interp_lvl_ps_nf(kp,1,2) - &
                                       trns_interp_lvl_ps_nf(kp,1,3))* &
                                  100.
            dgasdt8_lvlcts(kp)  = .02*(trns_interp_lvl_ps8_nf(kp,1,2) -&
                                       trns_interp_lvl_ps8_nf(kp,1,3))*&
                                  100.
            d2gast10_lvlcts(kp) = .0016*  &
                                  (trns_interp_lvl_ps_nf(kp,1,2) +   &
                                   trns_interp_lvl_ps_nf(kp,1,3) -  &
                                   2.0*trns_interp_lvl_ps_nf(kp,1,1))*&
                                   1000.
            d2gast8_lvlcts(kp)  = .0016*  &
                                  (trns_interp_lvl_ps8_nf(kp,1,2) +  &
                                   trns_interp_lvl_ps8_nf(kp,1,3) - &
                                  2.0*trns_interp_lvl_ps_nf(kp,1,1))*  &
                                  1000.
          enddo
        endif
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_lyrcalc) then
        do k=1,k2
          do kp=1,k1
            gasp10_lyr(kp,k) = c1*trns_interp_lyr_ps_nf(kp,k,1) - c2
            gasp8_lyr(kp,k) = c1*trns_interp_lyr_ps8_nf(kp,k,1) - c2
          enddo
        enddo
        if (ntbnd .EQ. 3) then
          do k=1,k2
            do kp=1,k1
              dgasdt10_lyr(kp,k) = .02*(trns_interp_lyr_ps_nf(kp,k,2) -&
                                        trns_interp_lyr_ps_nf(kp,k,3))*&
                                    100.
              dgasdt8_lyr(kp,k) = .02*(trns_interp_lyr_ps8_nf(kp,k,2) -&
                                       trns_interp_lyr_ps8_nf(kp,k,3))*&
                                                   100.
              d2gast10_lyr(kp,k) = .0016*  &
                                   (trns_interp_lyr_ps_nf(kp,k,2) +   &
                                    trns_interp_lyr_ps_nf(kp,k,3) -  &
                                    2.0*trns_interp_lyr_ps_nf(kp,k,1))*&
                                    1000.
              d2gast8_lyr(kp,k)  = .0016*  &
                                   (trns_interp_lyr_ps8_nf(kp,k,2) +   &
                                    trns_interp_lyr_ps8_nf(kp,k,3) - &
                                    2.0* &
                                    trns_interp_lyr_ps8_nf(kp,k,1))*  &
                                    1000.
            enddo
          enddo
        endif
      endif

!--------------------------------------------------------------------

 
end subroutine gasins




!#####################################################################
! <SUBROUTINE NAME="gasint">
!  <OVERVIEW>
!   gasint interpolates carbon dioxide transmission functions
!   from the standard level grid,for which the transmission functions
!   have been pre-calculated, to the grid structure specified by the
!   user.
!  </OVERVIEW>
!  <DESCRIPTION>
!   gasint interpolates carbon dioxide transmission functions
!   from the standard level grid,for which the transmission functions
!   have been pre-calculated, to the grid structure specified by the
!   user.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call gasint (gas_type, co2_vmr, co2_std_lo, co2_std_hi,  &
!                callrctrns, do_lvlcalc, do_lvlctscalc, do_lyrcalc,  &
!                nf, nt)
!  </TEMPLATE>
!  <IN NAME="gas_type" TYPE="character">
!   Gas type information
!  </IN>
!  <IN NAME="do_lvlcalc, do_lvlctscalc, do_lyrcalc" TYPE="logical">
!   State variables that determine calculation paths
!        do_lvlcalc     : compute level co2 transmissivities if true.
!        do_lyrcalc     : compute layer co2 transmissivities if true.
!        do_lvlctscalc  : compute cts level co2 transmissivities if true
!  </IN>
!  <IN NAME="nf" TYPE="integer">
!   frequency band number
!  </IN>
!  <IN NAME="nt" TYPE="integer">
!   temperature index of the frequency band
!  </IN>
!  <IN NAME="co2_vmr" TYPE="real">
!   co2 volume mixing ratio
!  </IN>
!  <IN NAME="co2_std_lo, co2_std_hi" TYPE="real">
!   standard co2 high and low volume mixing ratio (ppmv) pair
!  </IN>
!  <IN NAME="callrctrns" TYPE="logical">
!   state variable that determines calculation path
!  </IN>
! </SUBROUTINE>
!
subroutine gasint (gas_type, co2_vmr, co2_std_lo, co2_std_hi,  &
                   callrctrns, do_lvlcalc, do_lvlctscalc, do_lyrcalc,  &
                   nf, nt)
 
!---------------------------------------------------------------------
!    gasint interpolates carbon dioxide transmission functions
!    from the standard level grid,for which the transmission functions
!    have been pre-calculated, to the grid structure specified by the
!    user.
!
!        method: 
!
!    gasint is employable for two purposes: 1) to obtain transmis-
!    sivities between any 2 of an array of user-defined pressures; and
!    2) to obtain layer-mean transmissivities between any 2 of an array
!    of user-defined pressure layers.to clarify these two purposes,see
!    the diagram and discussion below.
!
!    let p be an array of user-defined pressures
!    and plm the array of user-defined level pressures 
!    and pd the extent of the interpolation layer.
!    for many purposes,plm will be chosen to be the average
!    pressure in the interpolation layer -i.e.,
!    plm(i) = 0.5*(pd(i-1)+pd(i)).
!
!       - - - - - - - - -   pd(i-1)  -----------!
!                                               !
!       -----------------   plm(i), p(k)-----!  !  interpolation layer
!                                            !  !
!       - - - - - - - - -   pd(i)    -----------!       model layer (i)
!                                            !
!       -----------------   plm(i+1)---------!
!            ...
!            ...                          (the notation used is
!            ...                          consistent with the code)
!            ...
!       - - - - - - - - -   pd(j-1)  -----------!
!                                               !
!       -----------------   plm(j), p(k')----!  !  interpolation layer
!                                            !  !
!       - - - - - - - - -   pd(j)    -----------!       model layer (j)
!                                            !
!       -----------------   plm(j+1)---------!
!
!    purpose 1:   the transmissivity between specific pressures
!    p(k) and p(k') ,tau(p(k),p(k'))  is computed by this program.
!    in this mode,there is no reference to layer pressures pd
!
!    purpose 2:   the transmissivity between a pressure p(k) and
!    the interpolation layer centered at p(k') (taulm(p(k),p(k'))
!    is obtained. it is computed by the integral
!
!                           pd(j)
!                           ----
!             1             !
!        -------------  *   !   tau ( p',plm(i) )  dp'
!        pd(j)-pd(j-1)      !
!                        ----
!                        pd(j-1)
!
!    the level pressures (plm) and layer-mean pressures (pd) are
!    both inputted in for this purpose.
!
!    in general, this integral is done using simpson's rule with
!    7 points. however , when plm(i) = pjm(j) (the nearby layer
!    case) a 51-point quadrature is used for greater accuracy.
!    note that taulm(k,k') is not a symmetrical matrix. also, no
!    quadrature is performed over the layer between the smallest 
!    nonzero pressure and zero pressure;
!    taulm is taulm(0,plm(j)) in this case,and taulm(0,0)=1.

!
!    the following paragraphs depict the utilization of this
!    code when used to compute transmissivities between specific
!    pressures. later paragraphs describe additional features needed
!    for layer-mean transmissivities.
!
!    for a given co2 mixing ratio and standard temperature
!    profile,a table of transmission functions for a fixed grid
!    of atmospheric pressures has been pre-calculated.
!    the standard temperature profile is from the us
!    standard atmosphere (1977) table.additionally, the
!    same transmission functions have been pre-calculated for a
!    temperature profile increased and decreased (at all levels)
!    by 25 degrees.
!    this program reads in the prespecified transmission functions
!    and a user-supplied pressure grids (p(k)) and calculates trans-
!    mission functions ,tau(p(k),p(k')), for all (k,k') using one
!    of the above methods, for one of the three temperature profiles.
!    outputs are tables of transmission functions. 
!
!
!    this code may be considered to be version 2 of the
!    interpolator. differences between this code and version 1, 
!    written in ~1983, are as follows:
!
!    1) the code is written using arrays (previous code was entirely
!       scalar)
!    2) double precision quantities are removed. it is assumed that
!       this code is being run on a 64-bit machine. if not, the
!       user will have to reinstate double precisioning, or
!       the appropriate KIND statement in Fortran 90.
!    3) many redundant calculations were eliminated
!    4) the error function is defined exactly as in Ref. (2).
!
!    as a result, this version of the code runs 100-200 times as fast
!    as version 1, and is suitable for on-line calculation of the
!    co2 transmission function.
!     
!                differences in the answers:
!
!    1) as part of the elimination of redundant calculation, the
!       following algorithmic change was performed:
!       calculation of the error function (error_guess1) at standard
!       pressures is done WITHOUT reference to model (user) pressures.
!       answers (in fractional absorptivity change) differ by 10-3 (max)
!       to 10-5. the new method should be "better", as there is no 
!       reason why the error function at standard pressures should care
!       about the pressures where it will be interpolated to.
!
!    2) in the "closely spaced" case (model pressures p,p' both
!       between standard pressures (pa(k),pa(k+1))) the coefficients 
!       (c,x,eta,sexp) are interpolated, not the approx function. this
!       is consistent with all other cases. fractional absorptivity 
!       changes are ~3x10-5 (max) and only for a(p,p') with p~ or = p'.
!
!    3) double precision to single precision change causes fractional
!       absorptivity changes of < 1x10-6.
!
!             references: 
!
!    (1): s.b.fels and m.d.schwarzkopf,"an efficient,accurate
!    algorithm for calculating co2 15 um band cooling rates",journal
!    of geophysical research,vol.86,no. c2, pp.1205-1232,1981.
!    (2): m.d. schwarzkopf and s.b. fels, "Improvements to the
!    algorithm for computing co2 transmissivities and cooling rates",
!    JGR, vol. 90, no. C10, pp10541-10550, 1985.
!
!            author:    m.daniel schwarzkopf
!
!            date:      14 july 1996
!
!            address: 
!
!                      GFDL
!                      p.o.box 308
!                      princeton,n.j.08542
!                      u.s.a.
!            telephone:  (609) 452-6521
!
!            e-mail:   Dan.Schwarzkopf@noaa.gov
!
!
!
!    NOTE: the comment line below is part of the original version
!    of this code, written in the late '70s by Stephen B. Fels. 
!    although the code has been extensively rewritten, and might
!    be unrecognizable to Steve, this line is kept as a tribute
!    to him.
!
!      ************   function interpolator routine  *****
!
!--------------------------------------------------------------------

logical,              intent(in)  ::  do_lvlcalc, do_lyrcalc,   &
                                      do_lvlctscalc, callrctrns
real,                 intent(in)  ::  co2_vmr, co2_std_lo, co2_std_hi
integer,              intent(in)  ::  nf, nt
character(len=*),     intent(in)  ::  gas_type
 
!--------------------------------------------------------------------
! miscellaneous variables: 
!
!    trns_std_lo : array of co2 transmission functions at the
!                  lower of two standard co2 concentrations
!                  at a given temperature profile.
!                  used if interpolation to the actual co2
!                  mixing ratio is required (callrctrns = true).
!                  dimensions: (NSTDCO2LVLS,NSTDCO2LVLS)
!    trns_std_hi : array of co2 transmission functions at the
!                  higher of two standard co2 concentrations
!                  at a given temperature profile.
!                  dimensions: (NSTDCO2LVLS,NSTDCO2LVLS)
!    co2_vmr  : actual co2 concentration (in ppmv)
!    co2_std_lo  : co2 concentration (ppmv) of lower of
!                  two standard concentrations.
!    co2_std_hi  : co2 concentration (ppmv) of higher of
!                  two standard concentrations.
!    callrctrns  : call rctrns.F if true.
! pressint_hiv_std_pt1  : allocated array used for rctrns hi pressure
! pressint_lov_std_pt1  : allocated array used for rctrns low pressure
! pressint_hiv_std_pt2  : allocated array used for rctrns hi pressure
! pressint_lov_std_pt2  : allocated array used for rctrns low pressure
! indx_pressint_hiv_std_pt1  : allocated index array used in rctrns
! indx_pressint_lov_std_pt1  : allocated index array used in rctrns
! indx_pressint_hiv_std_pt2  : allocated index array used in rctrns
! indx_pressint_lov_std_pt2  : allocated index array used in rctrns
!           do_lvlcalc  : compute level co2 transmissivities if true.
!           do_lyrcalc  : compute layer co2 transmissivities if true.
!        do_lvlctscalc  : compute cts level co2 transmissivities if true
!                   nf  : frequency band index
!                   nt  : temperature index (for the freq band)
!        ndimkp, ndimk  : extents of dimensions for output interp
!                         transmissivity arrays.
!              pd, plm  : see description below. note that the
!                         present limits on pressures (from the lbl
!                         computations require that the top level
!                         be 0 mb, and the bottom level pressure
!                         not exceed 1165 mb.
!             pd8, plm8 : same as pd, plm; highest pressure is 0.8*
!                         (highest pressure in plm).
!                         
!
!
!     outputs: 
!        trns_interp_lyr_ps : array of interpolated layer transmission
!           do_lvlcalc  : compute level co2 transmissivities if true.
!           do_lyrcalc  : compute layer co2 transmissivities if true.
!        do_lvlctscalc  : compute cts level co2 transmissivities if true
!                   nf  : frequency band index
!                   nt  : temperature index (for the freq band)
!        ndimkp, ndimk  : extents of dimensions for output interp
!                         transmissivity arrays.
!              pd, plm  : see description below. note that the
!                         present limits on pressures (from the lbl
!                         computations require that the top level
!                         be 0 mb, and the bottom level pressure
!                         not exceed 1165 mb.
!             pd8, plm8 : same as pd, plm; highest pressure is 0.8*
!                         (highest pressure in plm).
!                         
!
!
!     outputs: 
!        trns_interp_lyr_ps : array of interpolated layer transmission
!                             functions for the pressure array (pd).
!        trns_interp_lyr_ps8: array of interpolated layer transmission
!                             functions for the pressure array (pd8).
!        trns_interp_lyr_ps : array of interpolated level transmission
!                             functions for the pressure array (plm).
!        trns_interp_lyr_ps8: array of interpolated level transmission
!                             functions for the pressure array (plm8).
!

!--------------------------------------------------------------------
!  local variables

      real, dimension(NSTDCO2LVLS,NSTDCO2LVLS) :: trns_vmr,      &
                                                  approx_guess1, &
                                                  error_guess1,  &
                                                  press_hiv,     &
                                                  press_lov

      real, dimension(KSRAD:KERAD+1,KSRAD:KERAD+1) :: &
                                                  approxint_guess1, &
                                                  errorint_guess1,  &
                                                  pressint_lov,     &
                                                  pressint_hiv

      integer, dimension(KSRAD:KERAD+1,KSRAD:KERAD+1) :: &
                                                  indx_pressint_hiv, &
                                                  indx_pressint_lov

      real,    dimension(:,:), allocatable :: caintv, uexpintv, &
                                              sexpintv, xaintv

      real,    dimension(51,KERAD+1) :: sexpnblv, uexpnblv,  &
                                        canblv, xanblv,      &
                                        pressnbl_lov,        &
                                        pressnbl_hiv,        &
                                        approxnbl_guess1,    &
                                        errornbl_guess1

      integer, dimension(51,KERAD+1) :: indx_pressnbl_hiv, &
                                        indx_pressnbl_lov
 
      real, dimension(7)    ::  wgt_lyr
      real, dimension(51)   ::  wgt_nearby_lyr
      logical               ::  do_triangle
      integer               ::  n, k, kp, nklo, nkhi, nkplo, nkphi,   &
                                nq, nprofile

!--------------------------------------------------------------------
!    compute the layer weights for transmissivities. used only if
!    layer transmissivities are needed (do_lyrcalc = true)
!--------------------------------------------------------------------
      if (do_lyrcalc) then
        wgt_lyr(1) = 1./18.
        wgt_lyr(7) = 1./18.
        do n=1,3
          wgt_lyr(2*n) = 4./18.
        enddo
        do n=1,2
          wgt_lyr(2*n+1) = 2./18.
        enddo
        wgt_nearby_lyr(1) = 1./150.
        wgt_nearby_lyr(51) = 1./150.
        do n=1,25
          wgt_nearby_lyr(2*n) = 4./150.
        enddo
        do n=1,24
          wgt_nearby_lyr(2*n+1) = 2./150.
        enddo
      endif

!-------------------------------------------------------------------
!    define transmission function array for (co2_vmr) over
!    standard pressures (pa), using a call to rctrns if necessary.
!-------------------------------------------------------------------
      if (callrctrns) then
        call rctrns (gas_type, co2_std_lo, co2_std_hi, co2_vmr,  &
                     nf, nt, trns_vmr)
      else
        trns_vmr = trns_std_hi
      endif
                         
      do k=1,NSTDCO2LVLS
        trns_vmr(k,k)=1.0
      enddo
 
!-------------------------------------------------------------------
!    compute co2 transmission functions for actual co2 concentration
!    using method of section 5, Ref. (2).
!-------------------------------------------------------------------
      call coeint (gas_type, nf, trns_vmr, ca, sexp, xa, uexp)
 
!-------------------------------------------------------------------
!    compute the interpolation. 
!-------------------------------------------------------------------
      do_triangle = .true.
 
!-------------------------------------------------------------------
!    1) compute approx function at standard (pa) pressures
!-------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        press_hi(k) = pa(k)
        press_lo(k) = pa(k)
      enddo
 
!-------------------------------------------------------------------
!    allocate the 2-d input and output arrays needed to obtain the
!    approx function
!-------------------------------------------------------------------
      allocate ( caintv  (NSTDCO2LVLS,NSTDCO2LVLS))
      allocate ( sexpintv(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate ( xaintv  (NSTDCO2LVLS,NSTDCO2LVLS))
      allocate ( uexpintv(NSTDCO2LVLS,NSTDCO2LVLS))
   
!-------------------------------------------------------------------
!    compute the 2-d input arrays
!-------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        do kp=k,NSTDCO2LVLS
          press_hiv(kp,k) = pa(kp)
          press_lov(kp,k) = pa(k)
          caintv(kp,k) = ca(kp)
          sexpintv(kp,k) = sexp(kp)
          xaintv(kp,k) = xa(kp)
          uexpintv(kp,k) = uexp(kp)
        enddo
      enddo

!-------------------------------------------------------------------
!    the call (and calculations) to pathv2_std has been subsumed into
!    the subroutine approx_fn_std
!-------------------------------------------------------------------
      call approx_fn_std (press_hiv, press_lov, do_triangle, &
                          caintv, sexpintv, xaintv, uexpintv,  &
                          approx_guess1)

      deallocate (uexpintv)
      deallocate (xaintv)
      deallocate (sexpintv)
      deallocate (caintv)

!-------------------------------------------------------------------
!    2) compute error function at standard (pa) pressures
!-------------------------------------------------------------------
 
      do k=1,NSTDCO2LVLS
        do kp=k+1,NSTDCO2LVLS
          error_guess1(kp,k) = 1.0 - trns_vmr(kp,k) -  &
                               approx_guess1(kp,k)
        enddo
        error_guess1(k,k) = 0.0
      enddo
        
!-------------------------------------------------------------------
!    define the actual extents of the level interpolation calculation.
!    this depends on the type of calculation desired (whether
!    do_lvlcalc, do_lvlctscalc is true).
!-------------------------------------------------------------------
      if (do_lvlcalc) then
        nklo = KSRAD
        nkhi = KERAD + 1
        nkplo = KSRAD
        nkphi = KERAD + 1
      elseif (do_lvlctscalc) then
        nklo = KSRAD
        nkhi = KSRAD
        nkplo = KSRAD
        nkphi = KERAD + 1
      endif
 
!-------------------------------------------------------------------
!    allocate arrays with user-defined k-extents, which are used
!    in the remainder of the subroutine
!-------------------------------------------------------------------
      allocate ( caintv  (KSRAD:KERAD+1, KSRAD:KERAD+1) )
      allocate ( sexpintv(KSRAD:KERAD+1, KSRAD:KERAD+1) )
      allocate ( xaintv  (KSRAD:KERAD+1, KSRAD:KERAD+1) )
      allocate ( uexpintv(KSRAD:KERAD+1, KSRAD:KERAD+1) )
 
      if (do_lvlctscalc .OR. do_lvlcalc) then
        do k=KSRAD,KERAD+1
          trns_interp_lvl_ps(k,k) = 1.0
          trns_interp_lvl_ps8(k,k) = 1.0
        enddo
 
!-------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!-------------------------------------------------------------------
        do_triangle = .true.
        do nprofile = 1,2
          if (nprofile .EQ. 1) then
            do k=nklo,nkhi
              do kp=k+nkplo,nkphi
                pressint_hiv(kp,k) = plm(kp)
                pressint_lov(kp,k) = plm(k)
              enddo
            enddo
          else
            do k=nklo,nkhi
              do kp=k+nkplo,nkphi
                pressint_hiv(kp,k) = plm8(kp)
                pressint_lov(kp,k) = plm8(k)
              enddo
            enddo
          endif
          if (gas_type == 'co2' .and. (nf .ge. 6 .and. nf .le. 8)) then  !  co2 10um interpolation
            call intcoef_2d_10um (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv, uexpintv)
          else
            call intcoef_2d (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv, uexpintv)
          endif

!-------------------------------------------------------------------
!    4) interpolate error function to (pressint_hiv, pressint_lov)
!       for relevant (k',k)
!-------------------------------------------------------------------
          call interp_error (error_guess1, pressint_hiv, pressint_lov, &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             do_triangle, nklo, nkhi, nkplo, nkphi,  &
                             errorint_guess1)

!-------------------------------------------------------------------
!    5) compute approx function for (pressint_hiv, pressint_lov)
!       the call (and calculations) to pathv2 has been subsumed 
!       into subroutine approx_fn
!-------------------------------------------------------------------
          call approx_fn (pressint_hiv, pressint_lov, do_triangle,  &
                          nklo, nkhi, nkplo, nkphi,  &
                          caintv, sexpintv, xaintv, uexpintv,  &
                          approxint_guess1)
 
!-------------------------------------------------------------------
!    6) compute interp transmission function using Eq.(3),
!       Ref.(2).
!-------------------------------------------------------------------
          if (nprofile .EQ. 1) then
            do k=nklo,nkhi
              do kp=k+nkplo,nkphi
                trns_interp_lvl_ps(kp,k) = 1.0 -  &
                        (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                trns_interp_lvl_ps(k,kp) = trns_interp_lvl_ps(kp,k)
              enddo
            enddo
          else
            do k=nklo,nkhi
              do kp=k+nkplo,nkphi
                trns_interp_lvl_ps8(kp,k) = 1.0 -  &
                        (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                trns_interp_lvl_ps8(k,kp) = trns_interp_lvl_ps8(kp,k)
              enddo
            enddo
          endif  
        enddo  ! (nprofile loop)
      endif
 
      if (do_lyrcalc) then
!-------------------------------------------------------------------
!    A): calculate, for (kp,k) pairs with kp > k, a set of 7 transmis-
!    sivities, with the values of p'(kp) encompassing the layer bounded
!    by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!    averaged transmissivity (trns_interp_lyr_ps(8)(kp,k)).
!    B): calculate, for (kp,k) pairs with kp < k, a set of 7 transmis-
!    sivities, with the values of p'(kp) encompassing the layer bounded
!    by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!    averaged transmissivity (trns_interp_lyr_ps(8)(kp,k)).
!    C): calculate, for pairs (kp,kp) with kp > 1, a set of 51 transmis-
!    sivities, with the values of p'(kp) encompassing the layer bounded
!    by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!    averaged transmissivity (trns_interp_lyr_ps(8)(kp,k)).
!
!    note: one of the 7 (or 51) transmissivities equals the level
!    transmissivity (trns_interp_lvl_ps(8))
!
!    initialize the layer transmissivities to zero (summing later)
!    except the (1,1), which are set to 1
!-------------------------------------------------------------------
        trns_interp_lyr_ps = 0.0
        trns_interp_lyr_ps8 = 0.0
        trns_interp_lyr_ps(1,1) = 1.0
        trns_interp_lyr_ps8(1,1) = 1.0
 
!-------------------------------------------------------------------
!   case A): (kp) levels are at higher pressure, hence are used for
!            pressint_hiv. the (fixed) (k) levels are used for
!            pressint_lov
!-------------------------------------------------------------------
        do_triangle = .true.
        nklo = KSRAD
        nkhi = KERAD
        nkplo = KSRAD
        nkphi = KERAD + 1
 
!-------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!-------------------------------------------------------------------
        do nprofile = 1,2
          do nq = 1,7
            if (nprofile .EQ. 1) then
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  pressint_hiv(kp,k) = pd(kp-1) + (nq-1)*  &
                                       (pd(kp) - pd(kp-1))/6
                  pressint_lov(kp,k) = plm(k)
                enddo
              enddo
            else
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  pressint_hiv(kp,k) = pd8(kp-1) + (nq-1)*  &
                                       (pd8(kp) - pd8(kp-1))/6
                  pressint_lov(kp,k) = plm8(k)
                enddo
              enddo
            endif
            if (gas_type == 'co2' .and. (nf .ge. 6 .and. nf .le. 8)) then  !  co2 10um interpolation
              call intcoef_2d_10um (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv, uexpintv)
            else
              call intcoef_2d (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv, uexpintv)
            endif

!-------------------------------------------------------------------
!    4) interpolate error function to (pressint_hiv, pressint_lov)
!       for relevant (k',k)
!-------------------------------------------------------------------
            call interp_error (error_guess1, pressint_hiv, &
                               pressint_lov, indx_pressint_hiv,  &
                               indx_pressint_lov, do_triangle,  &
                               nklo, nkhi, nkplo, nkphi,  &
                               errorint_guess1)
 
!-------------------------------------------------------------------
!    5) compute approx function for (pressint_hiv, pressint_lov)
!       the call (and calculations) to pathv2 has been subsumed 
!       into subroutine approx_fn
!-------------------------------------------------------------------
            call approx_fn (pressint_hiv, pressint_lov, do_triangle,  &
                            nklo, nkhi, nkplo, nkphi,  &
                            caintv, sexpintv, xaintv, uexpintv,  &
                            approxint_guess1)
 
!-------------------------------------------------------------------
!    6) compute interp transmission function using Eq.(3),
!       Ref.(2).
!-------------------------------------------------------------------
            if (nprofile .EQ. 1) then
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  trns_interp_lyr_ps(kp,k) = trns_interp_lyr_ps(kp,k) +&
                                             wgt_lyr(nq)*(1.0 -   &
                                             (errorint_guess1(kp,k) +  &
                                             approxint_guess1(kp,k)))
 
!-------------------------------------------------------------------
!    for the case (nq=4), where  (pressint_hiv(kp,k) = plm(kp)) use
!    the (kp,1) unweighted values (errorint + approxint) for
!    the (1,kp) transmissivity, otherwise uncalculated. (exception:
!    when kp = nkphi, the (nq=7) value must be used)
!-------------------------------------------------------------------
                  if (nq .EQ. 4 .AND. k .EQ. nklo) then
                    trns_interp_lyr_ps(nklo,kp) = 1.0 -  &
                       (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                  endif
                enddo
              enddo
              if (nq .EQ. 7) then
                trns_interp_lyr_ps(nklo,nkphi) = 1.0 -  &
                                     (errorint_guess1(nkphi,nklo) +   &
                                      approxint_guess1(nkphi,nklo))
              endif
            else
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  trns_interp_lyr_ps8(kp,k) = &
                                         trns_interp_lyr_ps8(kp,k) +  &
                                         wgt_lyr(nq)*(1.0 -  &
                       (errorint_guess1(kp,k) + approxint_guess1(kp,k)))
 
!-------------------------------------------------------------------
!    for the case (nq=4), where  (pressint_hiv(kp,k) = plm(kp)) use
!    the (kp,1) unweighted values (errorint + approxint) for
!    the (1,kp) transmissivity, otherwise uncalculated. (exception:
!    when kp = nkphi, the (nq=7) value must be used)
!
!-------------------------------------------------------------------
                  if (nq .EQ. 4 .AND. k .EQ. nklo) then
                    trns_interp_lyr_ps8(nklo,kp) = 1.0 -  &
                        (errorint_guess1(kp,k) + approxint_guess1(kp,k))
                  endif
                enddo
              enddo
              if (nq .EQ. 7) then
                trns_interp_lyr_ps8(nklo,nkphi) = 1.0 -  &
                                     (errorint_guess1(nkphi,nklo) +   &
                                      approxint_guess1(nkphi,nklo))
              endif
            endif
          enddo
        enddo

!-------------------------------------------------------------------
!    case B): (k) levels are at higher pressure, hence are used for
!             pressint_hiv. the (variable) (kp) levels are used for
!             pressint_lov. (kp,k) calculations are loaded into
!             (k,kp) array locations to keep calculations into the
!             "upper sandwich". results are then put into their proper
!             array locations (before weighting function is applied).
!             also, no calculations are made for (1,k). these values
!             are obtained from level calcs for (k,1), nq=4.
!-------------------------------------------------------------------
        do_triangle = .true.
        nklo = KSRAD+1
        nkhi = KERAD
        nkplo = KSRAD
        nkphi = KERAD + 1

!-------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!-------------------------------------------------------------------
        do nprofile = 1,2
          do nq = 1,7
            if (nprofile .EQ. 1) then
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  pressint_hiv(kp,k) = plm(kp)
                  pressint_lov(kp,k) = pd(k-1) + (nq-1)*  &
                                       (pd(k) - pd(k-1))/6
                enddo
              enddo
            else
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  pressint_hiv(kp,k) = plm8(kp)
                  pressint_lov(kp,k) = pd8(k-1) + (nq-1)*  &
                                       (pd8(k) - pd8(k-1))/6
                enddo
              enddo
            endif
            if (gas_type == 'co2' .and. (nf .ge. 6 .and. nf .le. 8)) then  !  co2 10um interpolation
              call intcoef_2d_10um (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv, uexpintv)
            else
              call intcoef_2d (pressint_hiv, pressint_lov, do_triangle,  &
                             nklo, nkhi, nkplo, nkphi,  &
                             indx_pressint_hiv, indx_pressint_lov,  &
                             caintv, sexpintv, xaintv,uexpintv)
            endif

!-------------------------------------------------------------------
!    4) interpolate error function to (pressint_hiv, pressint_lov)
!       for relevant (k',k)
!-------------------------------------------------------------------
            call interp_error (error_guess1, pressint_hiv,    &
                               pressint_lov, indx_pressint_hiv,   &
                               indx_pressint_lov, do_triangle,  &
                               nklo, nkhi, nkplo, nkphi,  &
                               errorint_guess1)
 
!-------------------------------------------------------------------
!    5) compute approx function for (pressint_hiv, pressint_lov)
!       the call (and calculations) to pathv2 has been subsumed 
!       into subroutine approx_fn
!-------------------------------------------------------------------
            call approx_fn (pressint_hiv, pressint_lov, do_triangle,  &
                            nklo, nkhi, nkplo, nkphi,  &
                            caintv, sexpintv, xaintv, uexpintv,  &
                            approxint_guess1)
 
!-------------------------------------------------------------------
!    6) compute interp transmission function using Eq.(3),
!       Ref.(2).
!-------------------------------------------------------------------
            if (nprofile .EQ. 1) then
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  trns_interp_lyr_ps(k,kp) = trns_interp_lyr_ps(k,kp) +&
                                             wgt_lyr(nq)*(1.0 -  &
                       (errorint_guess1(kp,k) + approxint_guess1(kp,k)))
                enddo
              enddo
            else
              do k=nklo,nkhi
                do kp=k+nkplo,nkphi
                  trns_interp_lyr_ps8(k,kp) =    &
                                         trns_interp_lyr_ps8(k,kp) +  &
                                         wgt_lyr(nq)*(1.0 -  &
                       (errorint_guess1(kp,k) + approxint_guess1(kp,k)))
                enddo
              enddo
            endif
          enddo ! (nq loop)
        enddo   ! (nprofile loop)

!-------------------------------------------------------------------
!    C): calculate, for pairs (kp,kp) with kp > 1, a set of 51 transmis-
!    sivities, with the values of p'(kp) encompassing the layer bounded
!    by (pd(kp-1),pd(kp)). the weighted average of these is the layer-
!    averaged transmissivity (trns_interp_lyr_ps(8)(kp,k)).
!    case C): (kp) levels are at higher pressure, hence are used for
!           pressint_hiv. the (fixed) (k) levels are used for
!           pressint_lov
!-------------------------------------------------------------------
        do_triangle = .false.
        nklo = KSRAD + 1
        nkhi = KERAD + 1 
        nkplo = 1
        nkphi = 51

!-------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!-------------------------------------------------------------------
        do nprofile = 1,2
          if (nprofile .EQ. 1) then
            do k=nklo,nkhi-1
              do kp=1,25
                pressnbl_lov(kp,k) = pd(k-1) + (kp-1)*  &
                                     (pd(k) - pd(k-1))/50
                pressnbl_hiv(kp,k) = plm(k)
              enddo
              pressnbl_lov(26,k) = plm(k)
              pressnbl_hiv(26,k) = plm(k) + 1.0E-13*plm(k)
              do kp=27,51
                pressnbl_hiv(kp,k) = pd(k-1) + (kp-1)*  &
                                     (pd(k) - pd(k-1))/50
                pressnbl_lov(kp,k) = plm(k)
              enddo
            enddo
            do kp=1,50
              pressnbl_lov(kp,nkhi) = pd(nkhi-1) + (kp-1)*  &
                                      (pd(nkhi) - pd(nkhi-1))/50
              pressnbl_hiv(kp,nkhi) = plm(nkhi)
            enddo
            pressnbl_lov(51,nkhi) = plm(nkhi)
            pressnbl_hiv(51,nkhi) = plm(nkhi) + 1.0E-13*plm(nkhi)
          else
            do k=nklo,nkhi-1
              do kp=1,25
                pressnbl_lov(kp,k) = pd8(k-1) + (kp-1)*  &
                                     (pd8(k) - pd8(k-1))/50
                pressnbl_hiv(kp,k) = plm8(k)
              enddo
              pressnbl_lov(26,k) = plm8(k)
              pressnbl_hiv(26,k) = plm8(k) + 1.0E-13*plm8(k)
              do kp=27,51
                pressnbl_hiv(kp,k) = pd8(k-1) + (kp-1)*  &
                                     (pd8(k) - pd8(k-1))/50
                pressnbl_lov(kp,k) = plm8(k)
              enddo
            enddo
            do kp=1,50
              pressnbl_lov(kp,nkhi) = pd8(nkhi-1) + (kp-1)*  &
                                      (pd8(nkhi) - pd8(nkhi-1))/50
              pressnbl_hiv(kp,nkhi) = plm8(nkhi)
            enddo
            pressnbl_lov(51,nkhi) = plm8(nkhi)
            pressnbl_hiv(51,nkhi) = plm8(nkhi) + 1.0E-13*plm8(nkhi)
          endif
          if (gas_type == 'co2' .and. (nf .ge. 6 .and. nf .le. 8)) then  !  co2 10um interpolationn
            call intcoef_2d_10um (pressnbl_hiv, pressnbl_lov, do_triangle,  &
                           nklo, nkhi, nkplo, nkphi,  &
                           indx_pressnbl_hiv, indx_pressnbl_lov,  &
                           canblv, sexpnblv, xanblv, uexpnblv)
          else
            call intcoef_2d (pressnbl_hiv, pressnbl_lov, do_triangle,  &
                           nklo, nkhi, nkplo, nkphi,  &
                           indx_pressnbl_hiv, indx_pressnbl_lov,  &
                           canblv, sexpnblv, xanblv, uexpnblv)
          endif

!-------------------------------------------------------------------
!    4) interpolate error function to (pressnbl_hiv, pressnbl_lov)
!       for relevant (k',k)
!-------------------------------------------------------------------
          call interp_error (error_guess1, pressnbl_hiv, pressnbl_lov,&
                             indx_pressnbl_hiv, indx_pressnbl_lov,  &
                             do_triangle, nklo, nkhi, nkplo, nkphi,  &
                             errornbl_guess1)
 
!-------------------------------------------------------------------
!    5) compute approx function for (pressnbl_hiv, pressnbl_lov)
!       the call (and calculations) to pathv2 has been subsumed 
!       into subroutine approx_fn
!-------------------------------------------------------------------
          call approx_fn (pressnbl_hiv, pressnbl_lov, do_triangle,  &
                          nklo, nkhi, nkplo, nkphi,  &
                          canblv, sexpnblv, xanblv, uexpnblv,  &
                          approxnbl_guess1)
 
!-------------------------------------------------------------------
!    6) compute interp transmission function using Eq.(3),
!       Ref.(2).
!-------------------------------------------------------------------
          if (nprofile .EQ. 1) then
            do k=nklo,nkhi
              do kp=1,51
                trns_interp_lyr_ps(k,k) = trns_interp_lyr_ps(k,k) +  &
                                           wgt_nearby_lyr(kp)*   &
                (1.0 - (errornbl_guess1(kp,k) + approxnbl_guess1(kp,k)))
              enddo
            enddo
          else
            do k=nklo,nkhi
              do kp=1,51
                trns_interp_lyr_ps8(k,k) = trns_interp_lyr_ps8(k,k) +  &
                                            wgt_nearby_lyr(kp)*   &
                (1.0 - (errornbl_guess1(kp,k) + approxnbl_guess1(kp,k)))
              enddo
            enddo
          endif
        enddo    ! (nprofile loop)

      endif  ! do_lyrcalc

!-------------------------------------------------------------------
!      deallocate arrays with user-defined k-extents
!-------------------------------------------------------------------
      deallocate ( caintv )
      deallocate ( sexpintv )
      deallocate ( xaintv )
      deallocate ( uexpintv )

!---------------------------------------------------------------------


end subroutine gasint


!#####################################################################
! <SUBROUTINE NAME="coeint">
!  <OVERVIEW>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method
!  </DESCRIPTION>
!  <TEMPLATE>
!   call coeint (gas_type, nf, trns_val, ca, sexp, xa, uexp)
!  </TEMPLATE>
!  <IN NAME="gas_type" TYPE="character">
!   Gas type information
!  </IN>
!  <IN NAME="nf" TYPE="integer">
!   number of frequency band
!  </IN>
!  <IN NAME="trns_val" TYPE="real">
!   transmission function array
!  </IN>
!  <OUT NAME="ca, xa, sexp, uexp" TYPE="real">
!   coefficients in the transmission function between two pressure
!   levels
!  </OUT>
! </SUBROUTINE>
!
subroutine coeint (gas_type, nf, trns_val, ca, sexp, xa, uexp)

!-------------------------------------------------------------------
!    the transmission function between p1 and p2 is assumed to
!    have the  functional form
!         tau(p1,p2)= 1.0-(c*log(1.0+x*path**delta))**(gamma/delta),
!    where
!         path(p1,p2)=(p1-p2)**2)*(p1+p2+dop_core)
!         and p2 is the larger of the two pressures (p1,p2).

!    the coefficients c and x are functions of p2, while dop_core,
!    gamma and delta are predetermined coefficients.
!    (delta,gamma are uexp,sexp in this code).
!    subroutine coeint determines c(i) and x(i) by using actual
!    values of tau(p(i-2),p(i)) and tau(p(i-1),p(i)), obtained
!    from line-by-line calculations.
!    define: 
!             patha=(path(p(i),p(i-2),dop_core)**delta
!             pathb=(path(p(i),p(i-1),dop_core)**delta;
!    then
!         r=(1-tau(p(i),p(i-2)))/(1-tau(p(i),p(i-1)))
!          = (log(1+x(p(i))*patha)/log(1+x(p(i))*pathb))**(gamma/delta),
!    since   c(p(i)) cancels out
!    so that
!           r**(delta/gamma)= log(1+x(p(i))*patha)/log(1+x(p(i))*pathb).
!    this equation is solved by newton's method for x and then the
!    result used to find c. this is repeated for each value of i 
!    greater than 2 to give the arrays x(i), c(i).
!    there are several possible pitfalls: 
!       1) in the course of iteration, x may reach a value which makes
!          1+x*patha negative; in this case the iteration is stopped,
!          and an error message is printed out.
!       2) even if (1) does not occur, it is still possible that x may
!          be negative and large enough to make
!          1+x*path(p(i),0,dop_core) negative. this is checked in
!          a final loop, and if true,a warning is printed out.
!-----------------------------------------------------------------

character(len=*),       intent(in)  :: gas_type
real, dimension(:,:),   intent(in)  :: trns_val
integer,                intent(in)  :: nf
real, dimension(:),     intent(out) :: ca, xa, sexp, uexp

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      gas_type
!
!-------------------------------------------------------------------

!-----------------------------------------------------------------
!   local variables
      real, dimension(NSTDCO2LVLS) :: upath0, upatha, upathb,    &
                                      pam1, pam2, pa0, pr_hi, r, &
                                      rexp, f, f1, f2, fprime,   &
                                      ftest1, ftest2, xx, xxlog, &
                                      pa2, xx0, xxtest
      integer     :: k, ll
      integer     :: llmax
      real        :: check
      real, dimension(size(trns_val,1), size(trns_val,2)) :: abs_val
!-----------------------------------------------------------------
!   local variables
! 
!     upath0
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!    the following specifications for dop_core, sexp and uexp follow
!    "try9", which has (as of 5/27/97) been found to produce the
!    most accurate co2 40 level 490-850 cm-1 transmissivities, when
!    compared to LBL calculations over the same frequencies and 
!    vertical structure.
!--------------------------------------------------------------------
      if (gas_type .EQ. 'co2') then
        if (nf .eq. 1) dop_core = dop_core0
        if (nf .eq. 2) dop_core = dop_core0*560./670.
        if (nf .eq. 3) dop_core = dop_core0*665./670.
        if (nf .eq. 4) dop_core = dop_core0*775./670.
        if (nf .eq. 5) dop_core = dop_core0*2325./670.
        if (nf .eq. 6) dop_core = dop_core0*1030./670.
        if (nf .eq. 7) dop_core = dop_core0*945./670.
        if (nf .eq. 8) dop_core = dop_core0*1135./670.
      endif
      if (gas_type .EQ. 'ch4') then
        if (nf .eq. 1) dop_core = dop_core0*1300./670.
      endif
      if (gas_type .EQ. 'n2o') then
        if (nf .eq. 1) dop_core = dop_core0*1300./670.
        if (nf .eq. 2) dop_core = dop_core0*1135./670.
        if (nf .eq. 3) dop_core = dop_core0*595./670.
      endif
 
      do k=1,NSTDCO2LVLS
        pa2(k)=pa(k)*pa(k)
        sexp(k)=.505+2.0e-5*pa(k)+.035*(pa2(k)-.25)/(pa2(k)+.25)
        uexp(k) = sexp(k)*(1.0 + 0.33*pa2(k)/(pa2(k) + 40000.))
      enddo
 
      do k=1,NSTDCO2LVLS
        pr_hi(k) = pa(k)
      enddo
      do k=3,NSTDCO2LVLS
        pam1(k) = pa(k-1)
        pam2(k) = pa(k-2)
        pa0(k) = 0.0
      enddo
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      call pathv1 (pr_hi, pam1, 3, NSTDCO2LVLS, upathb)
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      call pathv1 (pr_hi, pam2, 3, NSTDCO2LVLS, upatha)
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      abs_val = MAX(1.0E-20,1.0 - trns_val)
      
      do k=3,NSTDCO2LVLS

        if (do_co2_bug) then 
          r(k) = (1.0 -trns_val(k,k-2))/(1.0 -trns_val(k,k-1))
        else
          r(k) = abs_val(k,k-2)/abs_val(k,k-1)
          if (r(k) .le. 1.0) then
            r(k) = 1.0
          endif
        endif 
!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(log(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!       rexp(k) = r(k)**(uexp(k)/sexp(k))
!       upatha(k) = upatha(k)**uexp(k)
!       upathb(k) = upathb(k)**uexp(k)
!------------------------------------------------------------------
        rexp(k) = EXP((uexp(k)/sexp(k))*log(r(k)))
        upatha(k) = EXP(uexp(k)*log(upatha(k)))
        upathb(k) = EXP(uexp(k)*log(upathb(k)))
        xx(k) = 2.0*(upathb(k)*rexp(k) - upatha(k))/   &
                (upathb(k)*upathb(k)*rexp(k) - upatha(k)*upatha(k))
        xx0(k) = xx(k)
        if (.not. do_co2_bug) then
          if (xx0(k) .le. 0.0) then
            xx0(k) = 0.0
          endif
        end if 
      enddo
!!    do ll=1,20
      if (gas_type .eq. 'co2') then
        llmax = llco2
      else if (gas_type .eq. 'ch4') then
        llmax = llch4
      else if (gas_type .eq. 'n2o') then
        llmax = lln2o
      endif
      do ll=1,llmax
        do k=3,NSTDCO2LVLS
          ftest1(k) =xx(k)*upatha(k)
          ftest2(k) =xx(k)*upathb(k)
!--------------------------------------------------------------------
!    end iteration and solve if ftest1 is small or ftest2 is large
!--------------------------------------------------------------------
          if (ftest1(k) .LE. 1.0E-10) then
            xa(k)=1.0
    
!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(log(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!           ca(k)=(1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/upatha(k)
!------------------------------------------------------------------
            if (do_co2_bug) then 
              ca(k)=EXP((uexp(k)/max(sexp(k),1.e-99))*   &
                  log((1.0 - trns_val(k,k-2))))/upatha(k)
            else     
              ca(k)=EXP((uexp(k)/max(sexp(k),1.e-99))*   &
                  log((abs_val(k,k-2))))/upatha(k)
              f(k) = 0.0
              fprime(k) = 1.0
            endif  
  
          else if (ftest2(k) .GE. 1.0E+8) then
            xxlog(k) = (LOG(upatha(k)) - rexp(k)*LOG(upathb(k)))/  &
                        (rexp(k)-1.0 )
            xa(k) = EXP(xxlog(k))

!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(log(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!           ca(k) = (1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/  &
!                   (xxlog(k) + LOG(upatha(k)))
!------------------------------------------------------------------
            if (do_co2_bug) then 
              ca(k)=EXP((uexp(k)/max(sexp(k),1.e-99))*   &
                 log((1.0 - trns_val(k,k-2))))/   &
                      (xxlog(k) + LOG(upatha(k)))
            else 
              ca(k)=EXP((uexp(k)/max(sexp(k),1.e-99))*   &
                  log((abs_val(k,k-2))))/   &
                       (xxlog(k) + LOG(upatha(k)))
              f(k) = 0.0
              fprime(k) = 1.0  
            end if 
          else
            f1(k) = LOG(1.0 + xx(k)*upatha(k))
            f2(k) = LOG(1.0 + xx(k)*upathb(k))
            f(k) = f1(k)/f2(k) - rexp(k)
            fprime(k) = (f2(k)*upatha(k)/(1.0 + xx(k)*upatha(k)) -  &
                         f1(k)*upathb(k)/(1.0 + xx(k)*upathb(k)))/  &
                         (f2(k)*f2(k))
 
            if (do_co2_bug) then
              xx(k) = xx(k) - f(k)/fprime(k)
            else
              xxtest(k) = xx(k) - f(k)/fprime(k)
              if (xx(k) .lt. 0.0) then
                xx(k) = xx0(k)
              else if (xxtest(k) .lt. 0.0) then
                xx(k) = xx0(k)
              else
                xx(k) = xxtest(k)
              endif
            endif
          endif
        enddo
      enddo

!--------------------------------------------------------------------
!    the following if loop is diagnostic only
!--------------------------------------------------------------------
      if (do_coeintdiag) then
        do k=3,NSTDCO2LVLS
          check=1.0 +xx(k)*upatha(k)
          if (check .le. 0.0) then
            write (     *, 360)  k, check
360         format ('check le zero, k=',i3, ' check =',f20.10)
            write(6,*) ' error, check le zero'
            stop
          endif
        enddo
      endif

      do k=3,NSTDCO2LVLS
        if (ftest1(k) .GT. 1.0E-10 .AND. ftest2(k) .LT. 1.0E+8) then
          if (do_co2_bug) then  
            ca(k) = (1.0 - trns_val(k,k-2))**(uexp(k)/sexp(k))/  &
                    (LOG(1.0 + xx(k)*upatha(k)) + 1.0e-20)
          else
            ca(k) = (abs_val(k,k-2))**(uexp(k)/sexp(k))/  &
                    (LOG(1.0 + xx(k)*upatha(k)) + 1.0e-20)
          endif
          xa(k) = xx(k)
       endif
      enddo

!----------------------------------------------------------------------
!    by assumption, ca, xa for the first two levels  are
!    equal to the values for the third level.
!----------------------------------------------------------------------
      xa(2)=xa(3)
      xa(1)=xa(3)
      ca(2)=ca(3)
      ca(1)=ca(3)
 
!--------------------------------------------------------------------
!    the following if loop is diagnostic only
!--------------------------------------------------------------------
      if (do_coeintdiag) then
        call pathv1 (pr_hi, pa0, 3, NSTDCO2LVLS, upath0)
        do k=3,NSTDCO2LVLS

!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(log(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!         upath0(k)=upath0(k)**uexp(k)
!------------------------------------------------------------------
          upath0(k)=EXP(uexp(k)*log(upath0(k)))
          upath0(k)=1.0 +xa(k)*upath0(k)
          if (upath0(k).lt.0.)   then
            write (     *, 361) k, upath0(k), xa(k) 
361         format (' 1+xa*path(pa(i),0) is negative,i= ',i3,/  &
                    20x,'upath0(i)=',f16.6,' xa(i)=',f16.6)
            write(6,*) '1+xa*path(pa(i),0) is negative'
            stop
          endif
        enddo
      endif

!--------------------------------------------------------------------


end subroutine coeint


!#####################################################################
! <SUBROUTINE NAME="intcoef_1d">
!  <OVERVIEW>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (1 dimensional)
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (1 dimensional)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call intcoef_1d (press_hi, press_lo, cav, sexpv, xav, uexpv)
!  </TEMPLATE>
!  <IN NAME="press_hi, press_lo" TYPE="real">
!   high and low pressure pair
!  </IN>
!  <OUT NAME="cav, xav, sexpv, uexpv" TYPE="real">
!   coefficients in the transmission function between two pressure
!   levels
!  </OUT>
! </SUBROUTINE>
!
subroutine intcoef_1d (press_hi, press_lo, cav, sexpv, xav, uexpv)
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

real, dimension (:), intent(in)  :: press_hi, press_lo
real, dimension (:), intent(out) :: cav, sexpv, xav, uexpv

!--------------------------------------------------------------------
!  intent(in) vaiables:
!
!       press_hi
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      real,    dimension(NSTDCO2LVLS) :: caxa, ca_hi, prod_hi, &
                                         sexp_hi, uexp_hi, xa_hi
      integer, dimension(NSTDCO2LVLS) :: indx_press_hi, indx_press_lo

      integer         :: k, kp, kpp

!-------------------------------------------------------------------
!  local variables:
!
!      caxa
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    compute the index of press_hi and press_lo  corresponding to pa
!---------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        if (press_hi(k) .LT. pa(1)) then
          indx_press_hi(k) = 1
        endif
        if (press_hi(k) .GE. pa(NSTDCO2LVLS)) then
          indx_press_hi(k) = NSTDCO2LVLS - 1
        endif
        if (press_lo(k) .LT. pa(1)) then
          indx_press_lo(k) = 1
        endif
        if (press_lo(k) .GE. pa(NSTDCO2LVLS)) then
          indx_press_lo(k) = NSTDCO2LVLS - 1
        endif
      enddo
      do k=1,NSTDCO2LVLS
        do kpp=1,NSTDCO2LVLS - 1
          if (press_hi(k) .GE. pa(kpp) .AND.  &
              press_hi(k) .LT. pa(kpp+1)) then
            indx_press_hi(k) = kpp
          endif
          if (press_lo(k) .GE. pa(kpp) .AND.  &
              press_lo(k) .LT. pa(kpp+1)) then
            indx_press_lo(k) = kpp
          endif
        enddo
      enddo

!--------------------------------------------------------------------
!    interpolate values of cint, xint, sexp, for the pressures
!    (press_hi)
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        caxa(k) = ca(k)*xa(k)
      enddo
      do k=1,NSTDCO2LVLS
        sexp_hi(k) = sexp(indx_press_hi(k)) +   &
                     (sexp(indx_press_hi(k)+1) -    &
                      sexp(indx_press_hi(k))) /  &
                     (pa  (indx_press_hi(k)+1) -    &
                      pa  (indx_press_hi(k))) *  &
                     (press_hi(k) - pa(indx_press_hi(k)))
        uexp_hi(k) = uexp(indx_press_hi(k)) +   &
                     (uexp(indx_press_hi(k)+1) -    &
                     uexp(indx_press_hi(k))) /  &
                     (pa  (indx_press_hi(k)+1) -    &
                      pa  (indx_press_hi(k))) *  &
                     (press_hi(k) - pa(indx_press_hi(k)))
        prod_hi(k) = caxa(indx_press_hi(k)) +   &
                     (caxa(indx_press_hi(k)+1) -   &
                     caxa(indx_press_hi(k))) /  &
                     (pa  (indx_press_hi(k)+1) -    &
                      pa  (indx_press_hi(k))) *  &
                     (press_hi(k) - pa(indx_press_hi(k)))
        xa_hi(k) = xa(indx_press_hi(k)) +   &
                   (xa(indx_press_hi(k)+1) - xa(indx_press_hi(k))) /  &
                   (pa  (indx_press_hi(k)+1) -     &
                    pa  (indx_press_hi(k))) *  &
                   (press_hi(k) - pa(indx_press_hi(k)))
                   ca_hi(k) = prod_hi(k)/xa_hi(k)
      enddo
 
!-------------------------------------------------------------------
!
!-------------------------------------------------------------------
      do kp=k, NSTDCO2LVLS
        sexpv(kp)     = sexp_hi(kp)
        uexpv(kp)     = uexp_hi(kp)
        cav(kp)     = ca_hi(kp)
        xav(kp)     = xa_hi(kp)
      enddo

!-------------------------------------------------------------------


end subroutine intcoef_1d


!#####################################################################
! <SUBROUTINE NAME="intcoef_2d">
!  <OVERVIEW>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (2 dimensional)
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (2 dimensional)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call intcoef_2d (press_hiv, press_lov, do_triangle,  &
!                    nklo, nkhi, nkplo, nkphi,  &
!                    indx_hiv, indx_lov,  &
!                    caintv,  sexpintv, xaintv, uexpintv)
!  </TEMPLATE>
!  <IN NAME="press_hiv, press_lov" TYPE="real">
!   high and low pressure pair
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   State variable of triangle interpolation scheme
!  </IN>
!  <IN NAME="nklo, nkhi, nkplo, nkphi" TYPE="integer">
!   the high and low level and pressure pairs
!  </IN>
!  <IN NAME="indx_hiv, indx_lov" TYPE="integer">
!   the high and low index pair
!  </IN>
!  <OUT NAME="caintv, xaintv, sexpintv, uexpintv" TYPE="real">
!   coefficients in the transmission function between two pressure
!   levels
!  </OUT>
! </SUBROUTINE>
!
subroutine intcoef_2d (press_hiv, press_lov, do_triangle,  &
                       nklo, nkhi, nkplo, nkphi,  &
                       indx_hiv, indx_lov,  &
                       caintv,  sexpintv, xaintv, uexpintv)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real,    dimension(:,:), intent(out) :: sexpintv, &
                                        uexpintv, caintv, xaintv
integer, dimension(:,:), intent(out) :: indx_hiv, indx_lov
integer,                 intent(in)  :: nklo, nkhi, nkplo, nkphi
real,    dimension(:,:), intent(in)  :: press_hiv, press_lov
logical,                 intent(in)  :: do_triangle

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    sexpintv
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      real, dimension(NSTDCO2LVLS) :: caxa
      real, dimension(size(press_hiv,1),size(press_hiv,2)) :: &
                                  sexp_hiv, uexp_hiv, ca_hiv, &
                                  prod_hiv, xa_hiv, d1kp,     &
                                  d2kp, bkp, akp, delp_hi

      integer    :: k, kp, kp0, kpp

!-------------------------------------------------------------------
!  local variables:
!
!      caxa
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    compute the index of the inputted pressures (press_hiv,
!    press_lov) corresponding to the standard (pa) pressures.
!---------------------------------------------------------------------
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          if (press_hiv(kp,k) .LT. pa(1)) then
            indx_hiv(kp,k) = 1
          endif
          if (press_hiv(kp,k) .GE. pa(NSTDCO2LVLS)) then
            indx_hiv(kp,k) = NSTDCO2LVLS - 1
          endif
          if (press_lov(kp,k) .LT. pa(1)) then
            indx_lov(kp,k) = 1
          endif
          if (press_lov(kp,k) .GE. pa(NSTDCO2LVLS)) then
            indx_lov(kp,k) = NSTDCO2LVLS - 1
          endif
        enddo
      enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          do kpp=1,NSTDCO2LVLS - 1
            if (press_hiv(kp,k) .GE. pa(kpp) .AND.  &
                press_hiv(kp,k) .LT. pa(kpp+1)) then
              indx_hiv(kp,k) = kpp
              exit
            endif
          enddo
          do kpp=1,NSTDCO2LVLS - 1
            if (press_lov(kp,k) .GE. pa(kpp) .AND.  &
                press_lov(kp,k) .LT. pa(kpp+1)) then
              indx_lov(kp,k) = kpp
              exit
            endif
          enddo
        enddo
      enddo
 
!---------------------------------------------------------------------
!    interpolate values of cint, xint, sexp, uexp for the pressures
!    (press_hiv)
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        caxa(k) = ca(k)*xa(k)
      enddo
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          sexp_hiv(kp,k) = sexp(indx_hiv(kp,k)) +   &
                   (sexp(indx_hiv(kp,k)+1) - sexp(indx_hiv(kp,k))) /  &
                   (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                   (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
           uexp_hiv(kp,k) = uexp(indx_hiv(kp,k)) +   &
                   (uexp(indx_hiv(kp,k)+1) - uexp(indx_hiv(kp,k))) /  &
                   (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                   (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
 
!--------------------------------------------------------------------
!    use 3-point interpolation: (indx_hiv of 1 or 2 are excluded
!    since ca and xa were arbitrarily set to ca(3),xa(3))
!--------------------------------------------------------------------
          if (indx_hiv(kp,k) .GT. 2 .AND.                        &
              indx_hiv(kp,k) .LT. NSTDCO2LVLS - 1) then     
            delp_hi(kp,k) =                           &
                 press_hiv(kp,k) - pa(indx_hiv(kp,k)+1)

!------------------------------------------------------------------
!    interpolate xa
!------------------------------------------------------------------
            d1kp(kp,k) =   &
              (xa(indx_hiv(kp,k)+2) - xa(indx_hiv(kp,k)+1)) /  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp,k) =   &
              (xa(indx_hiv(kp,k)+1) -  xa(indx_hiv(kp,k) )) /  &
              (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            xa_hiv(kp,k) =   &
              xa(indx_hiv(kp,k)+1) +                     &
              delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))

!----------------------------------------------------------------------
!    if xa_hiv is negative or zero, the interpolation fails and
!    the model may bomb. to avoid this, use 2-point interpolation
!    in this case. the 3-point interpolation for prod_hiv is
!    stable, so there is no need to change this calculation.
!----------------------------------------------------------------------
            if (xa_hiv(kp,k) .LE. 0.0E+00) then                 
              xa_hiv(kp,k) = xa(indx_hiv(kp,k)) +                  &
                 (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /     &
                 (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *   &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            endif
 
!----------------------------------------------------------------------
!    interpolate caxa
!----------------------------------------------------------------------
            d1kp(kp,k) =   &
              (caxa(indx_hiv(kp,k)+2) - caxa(indx_hiv(kp,k)+1)) /  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp,k) =   &
              (caxa(indx_hiv(kp,k)+1) -  caxa(indx_hiv(kp,k) )) /  &
              (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            prod_hiv(kp,k) =   &
              caxa(indx_hiv(kp,k)+1) +  &
                delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))
          else
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)) +   &
               (caxa(indx_hiv(kp,k)+1) - caxa(indx_hiv(kp,k))) /  &
               (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            xa_hiv(kp,k) = xa(indx_hiv(kp,k)) +   &
               (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /  &
               (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          endif

!---------------------------------------------------------------------
!    compute ca
!---------------------------------------------------------------------
          if (do_co2_bug) then
            ca_hiv(kp,k) = prod_hiv(kp,k)/xa_hiv(kp,k)
          else
            ca_hiv(kp,k) = prod_hiv(kp,k)/xa_hiv(kp,k)
            if (ca_hiv(kp,k) .lt. 0.0) then
              ca_hiv(kp,k) = 0!xa_hiv(kp,k) ** (-1./uexp_hiv(kp,k)) 
            endif
          endif
        enddo
      enddo
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          sexpintv(kp,k)     = sexp_hiv(kp,k)
          uexpintv(kp,k)     = uexp_hiv(kp,k)
          caintv(kp,k)     = ca_hiv(kp,k)
          xaintv(kp,k)     = xa_hiv(kp,k)
        enddo
      enddo

!--------------------------------------------------------------------- 

end subroutine intcoef_2d

subroutine intcoef_2d_10um (press_hiv, press_lov, do_triangle,  &
                       nklo, nkhi, nkplo, nkphi,  &
                       indx_hiv, indx_lov,  &
                       caintv,  sexpintv, xaintv, uexpintv)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real,    dimension(:,:), intent(out) :: sexpintv, &
                                        uexpintv, caintv, xaintv
integer, dimension(:,:), intent(out) :: indx_hiv, indx_lov
integer,                 intent(in)  :: nklo, nkhi, nkplo, nkphi
real,    dimension(:,:), intent(in)  :: press_hiv, press_lov
logical,                 intent(in)  :: do_triangle

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    sexpintv
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      real, dimension(NSTDCO2LVLS) :: caxa
      real, dimension(size(press_hiv,1),size(press_hiv,2)) :: &
                                  sexp_hiv, uexp_hiv, ca_hiv, &
                                  prod_hiv, xa_hiv, d1kp,     &
                                  d2kp, bkp, akp, delp_hi

      integer    :: k, kp, kp0, kpp
      integer    :: k1, k2

!-------------------------------------------------------------------
!  local variables:
!
!      caxa
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    compute the index of the inputted pressures (press_hiv,
!    press_lov) corresponding to the standard (pa) pressures.
!---------------------------------------------------------------------
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          if (press_hiv(kp,k) .LT. pa(1)) then
            indx_hiv(kp,k) = 1
          endif
          if (press_hiv(kp,k) .GE. pa(NSTDCO2LVLS)) then
            indx_hiv(kp,k) = NSTDCO2LVLS - 1
          endif
          if (press_lov(kp,k) .LT. pa(1)) then
            indx_lov(kp,k) = 1
          endif
          if (press_lov(kp,k) .GE. pa(NSTDCO2LVLS)) then
            indx_lov(kp,k) = NSTDCO2LVLS - 1
          endif
        enddo
      enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          do kpp=1,NSTDCO2LVLS - 1
            if (press_hiv(kp,k) .GE. pa(kpp) .AND.  &
                press_hiv(kp,k) .LT. pa(kpp+1)) then
              indx_hiv(kp,k) = kpp
              exit
            endif
          enddo
          do kpp=1,NSTDCO2LVLS - 1
            if (press_lov(kp,k) .GE. pa(kpp) .AND.  &
                press_lov(kp,k) .LT. pa(kpp+1)) then
              indx_lov(kp,k) = kpp
              exit
            endif
          enddo
        enddo
      enddo
 
!---------------------------------------------------------------------
!    interpolate values of cint, xint, sexp, uexp for the pressures
!    (press_hiv)
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        caxa(k) = ca(k)*xa(k)
      enddo
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          sexp_hiv(kp,k) = sexp(indx_hiv(kp,k)) +   &
                   (sexp(indx_hiv(kp,k)+1) - sexp(indx_hiv(kp,k))) /  &
                   (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                   (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          uexp_hiv(kp,k) = uexp(indx_hiv(kp,k)) +   &
                   (uexp(indx_hiv(kp,k)+1) - uexp(indx_hiv(kp,k))) /  &
                   (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                   (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
 
!--------------------------------------------------------------------
!    use 2-point interpolation: (indx_hiv of 1 or 2 are excluded
!    since ca and xa were arbitrarily set to ca(3),xa(3))
!--------------------------------------------------------------------
          prod_hiv(kp,k) = caxa(indx_hiv(kp,k)) +   &
               (caxa(indx_hiv(kp,k)+1) - caxa(indx_hiv(kp,k))) /  &
               (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          xa_hiv(kp,k) = xa(indx_hiv(kp,k)) +   &
               (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /  &
               (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))

!---------------------------------------------------------------------
!    compute ca
!---------------------------------------------------------------------
          if (do_co2_bug) then
            ca_hiv(kp,k) = prod_hiv(kp,k)/xa_hiv(kp,k)
          else
            ca_hiv(kp,k) = prod_hiv(kp,k)/xa_hiv(kp,k)
            if (ca_hiv(kp,k) .lt. 0.0) then
              ca_hiv(kp,k) = 0!xa_hiv(kp,k) ** (-1./uexp_hiv(kp,k)) 
            endif
          endif
        enddo
      enddo
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          sexpintv(kp,k)     = sexp_hiv(kp,k)
          uexpintv(kp,k)     = uexp_hiv(kp,k)
          caintv(kp,k)     = ca_hiv(kp,k)
          xaintv(kp,k)     = xa_hiv(kp,k)
        enddo
      enddo

!---------------------------------------------------------------------

end subroutine intcoef_2d_10um                       

!#####################################################################
! <SUBROUTINE NAME="intcoef_2d_std">
!  <OVERVIEW>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (2 dimensional)
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to inverse coefficients from transmission functions
!   using newton method (2 dimensional)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call intcoef_2d (press_hiv, press_lov, nf, nt, do_triangle,  &
!                    indx_hiv, indx_lov,  &
!                    caintv,  sexpintv, xaintv, uexpintv)
!  </TEMPLATE>
!  <IN NAME="press_hiv, press_lov" TYPE="real">
!   high and low pressure pair
!  </IN>
!  <IN NAME="nf" TYPE="integer">
!   number of frequency bands
!  </IN>
!  <IN NAME="nt" TYPE="integer">
!   number of temperature values
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   State variable of triangle interpolation scheme
!  </IN>
!  <IN NAME="nklo, nkhi, nkplo, nkphi" TYPE="integer">
!   the high and low level and pressure pairs
!  </IN>
!  <IN NAME="indx_hiv, indx_lov" TYPE="integer">
!   the high and low index pair
!  </IN>
!  <OUT NAME="caintv, xaintv, sexpintv, uexpintv" TYPE="real">
!   coefficients in the transmission function between two pressure
!   levels
!  </OUT>
! </SUBROUTINE>
!
subroutine intcoef_2d_std (press_hiv, press_lov, nf, nt, do_triangle,  &
                           indx_hiv, indx_lov,  &
                           caintv,  sexpintv, xaintv, uexpintv)

!---------------------------------------------------------------------

integer,                 intent(in)   :: nf, nt
real,    dimension(:,:), intent(in)   :: press_hiv, press_lov
logical,                 intent(in)   :: do_triangle
real,    dimension(:,:), intent(out)  :: sexpintv, uexpintv, caintv,  &
                                         xaintv
integer, dimension(:,:), intent(out)  :: indx_hiv, indx_lov

!----------------------------------------------------------------------
!  intent(in) variables:
!
!      nf
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:

      real, dimension(NSTDCO2LVLS,NSTDCO2LVLS) :: prod_hiv
      real, dimension(NSTDCO2LVLS) :: d1kp, d2kp, bkp, akp, &
                                      delp_hi, caxa
      integer         :: k, kp, kp0, kpp
 
!--------------------------------------------------------------------
!  local variables:
!
!     prod_hiv
!
!--------------------------------------------------------------------
     indx_hiv = 0   
      
!--------------------------------------------------------------------
!    compute the index of the inputted pressures (press_hiv,
!    press_lov) corresponding to the standard (pa) pressures.
!    (only calculate if nf = 1, nt = 1)
!--------------------------------------------------------------------
      !if (nf .EQ. 1 .AND. nt .EQ. 1) then
        do k=1,NSTDCO2LVLS
          if (do_triangle) then
            kp0 = k + 1
          else
            kp0 = 1
          endif
          do kp=kp0,NSTDCO2LVLS
            if (press_hiv(kp,k) .LT. pa(1)) then
              indx_hiv(kp,k) = 1
            endif
            if (press_hiv(kp,k) .GE. pa(NSTDCO2LVLS)) then
              indx_hiv(kp,k) = NSTDCO2LVLS - 1
            endif
            if (press_lov(kp,k) .LT. pa(1)) then
              indx_lov(kp,k) = 1
            endif
            if (press_lov(kp,k) .GE. pa(NSTDCO2LVLS)) then
              indx_lov(kp,k) = NSTDCO2LVLS - 1
            endif
          enddo
        enddo

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
        do k=1,NSTDCO2LVLS
          if (do_triangle) then
            kp0 = k + 1
          else
            kp0 = 1
          endif
          do kp=kp0,NSTDCO2LVLS
            do kpp=1,NSTDCO2LVLS - 1
              if (press_hiv(kp,k) .GE. pa(kpp) .AND.  &
                  press_hiv(kp,k) .LT. pa(kpp+1)) then
                indx_hiv(kp,k) = kpp
                exit
              endif
            enddo
            do kpp=1,NSTDCO2LVLS - 1
              if (press_lov(kp,k) .GE. pa(kpp) .AND.  &
                  press_lov(kp,k) .LT. pa(kpp+1)) then
                indx_lov(kp,k) = kpp
                exit
              endif
            enddo
          enddo
        enddo
      !endif
 
!--------------------------------------------------------------------
!    interpolate values of cint, xint, sexp, uexp for the pressures
!    (press_hiv) (for all values of nf, nt)
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        caxa(k) = ca(k)*xa(k)
      enddo

      do k=1,NSTDCO2LVLS
        if (do_triangle) then
          kp0 = k + 1
        else
          kp0 = 1
        endif
        do kp=kp0,NSTDCO2LVLS
          sexpintv(kp,k) = sexp(indx_hiv(kp,k)) +   &
                  (sexp(indx_hiv(kp,k)+1) - sexp(indx_hiv(kp,k))) /  &
                  (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                        (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          uexpintv(kp,k) = uexp(indx_hiv(kp,k)) +   &
                  (uexp(indx_hiv(kp,k)+1) - uexp(indx_hiv(kp,k))) /  &
                  (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                        (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
 
!--------------------------------------------------------------------
!    use 3-point interpolation: (indx_hiv of 1 or 2 are excluded
!    since ca and xa were arbitrarily set to ca(3),xa(3))
!--------------------------------------------------------------------
           if (indx_hiv(kp,k) .GT. 2 .AND.  &
              indx_hiv(kp,k) .LT. NSTDCO2LVLS - 1) then     
             delp_hi(kp) =                    &
                 press_hiv(kp,k) - pa(indx_hiv(kp,k)+1)

!---------------------------------------------------------------------
!    interpolate xa
!---------------------------------------------------------------------
            d1kp(kp) =   &
              (xa(indx_hiv(kp,k)+2) - xa(indx_hiv(kp,k)+1)) /  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp) =   &
              (xa(indx_hiv(kp,k)+1) -  xa(indx_hiv(kp,k) )) /  &
              (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp) = d1kp(kp) - bkp(kp)*  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            xaintv(kp,k) =   &
              xa(indx_hiv(kp,k)+1) +  &
                delp_hi(kp)*(akp(kp) + delp_hi(kp)*bkp(kp))

!--------------------------------------------------------------------
!    if xaintv is negative or zero, the interpolation fails and
!    the model may bomb. to avoid this, use 2-point interpolation
!    in this case. the 3-point interpolation for prod_hiv is
!    stable, so there is no need to change this calculation.
!--------------------------------------------------------------------
            if (xaintv(kp,k) .LE. 0.0E+00) then                 
              xaintv(kp,k) = xa(indx_hiv(kp,k)) +                  &
                 (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /     &
                  (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *   &
                            (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            endif
 
!-------------------------------------------------------------------
!    interpolate caxa
!-------------------------------------------------------------------
            d1kp(kp) =   &
              (caxa(indx_hiv(kp,k)+2) - caxa(indx_hiv(kp,k)+1)) /  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            d2kp(kp) =   &
              (caxa(indx_hiv(kp,k)+1) -  caxa(indx_hiv(kp,k) )) /  &
              (pa(indx_hiv(kp,k)+1) - pa(indx_hiv(kp,k)  ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)  ))
            akp(kp) = d1kp(kp) - bkp(kp)*  &
              (pa(indx_hiv(kp,k)+2) - pa(indx_hiv(kp,k)+1))
            prod_hiv(kp,k) =   &
              caxa(indx_hiv(kp,k)+1) +  &
                delp_hi(kp)*(akp(kp) + delp_hi(kp)*bkp(kp))
          else
            prod_hiv(kp,k) = caxa(indx_hiv(kp,k)) +   &
                 (caxa(indx_hiv(kp,k)+1) - caxa(indx_hiv(kp,k))) /  &
                  (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
            xaintv(kp,k) = xa(indx_hiv(kp,k)) +   &
                 (xa(indx_hiv(kp,k)+1) - xa(indx_hiv(kp,k))) /  &
                  (pa  (indx_hiv(kp,k)+1) - pa  (indx_hiv(kp,k))) *  &
                          (press_hiv(kp,k) - pa(indx_hiv(kp,k)))
          endif

!------------------------------------------------------------------
!    compute ca
!------------------------------------------------------------------
          if (do_co2_bug) then
            caintv(kp,k) = prod_hiv(kp,k)/xaintv(kp,k)
          else
            caintv(kp,k) = prod_hiv(kp,k)/xaintv(kp,k)
            if (caintv(kp,k) .lt. 0.0) then
              caintv(kp,k) = xaintv(kp,k) ** (-1./uexpintv(kp,k))
            endif 
          endif 
        enddo  ! (kp loop)
      enddo   ! (k loop)
 
!---------------------------------------------------------------------


end subroutine intcoef_2d_std


!#####################################################################
! <SUBROUTINE NAME="interp_error">
!  <OVERVIEW>
!   Subroutine to examine error associated with interpolation onto
!   pressure grids.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to examine error associated with interpolation onto
!   pressure grids.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call interp_error (error, pressint_hiv, pressint_lov,  &
!                         indx_press_hiv, indx_press_lov,  &
!                         do_triangle, nklo, nkhi, nkplo, nkphi,  &
!                         errorint)
!  </TEMPLATE>
!  <IN NAME="error" TYPE="real">
!   interpolation error at standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!  </IN>
!  <IN NAME="pressint_hiv, pressint_lov" TYPE="real">
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     pressint_lov = pressure of low (kp) interpolated pressure
!  </IN>
!  <IN NAME="indx_press_hiv, indx_press_lov" TYPE="real">
!   indx_press_hiv = pressure on std pa grid of high (kp) pressure
!   indx_press_lov = pressure on std pa grid of low  (kp) pressure
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   state variable that determines the interpolation scheme
!  </IN>
!  <IN NAME="nkl, nkhi, nkplo, nkphi" TYPE="integer">
!   The index of level and pressure high/low pair
!  </IN>
!  <OUT NAME="errorint" TYPE="real">
!   error at interpolated grid
!  </OUT>
! </SUBROUTINE>
!
subroutine interp_error (error, pressint_hiv, pressint_lov,  &
                         indx_press_hiv, indx_press_lov,  &
                         do_triangle, nklo, nkhi, nkplo, nkphi,  &
                         errorint)

!--------------------------------------------------------------------
!
!---------------------------------------------------------------------

real,    dimension(:,:), intent(in)   :: error,   &
                                         pressint_hiv, pressint_lov
integer, dimension(:,:), intent(in)   :: indx_press_hiv, indx_press_lov
logical,                 intent(in)   :: do_triangle
integer,                 intent(in)   :: nklo, nkhi, nkplo, nkphi
real,    dimension(:,:), intent(out)  :: errorint

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     error
!     press_hiv = pressure on std pa grid of high (kp) pressure
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     error = error ot standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!     errorint = error at interpolated grid
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension(size(pressint_hiv,1),size(pressint_hiv,2)) :: &
                                         delp_lo, delp_hi,          &
                                         d1kp, d2kp, bkp, akp, fkp, &
                                         fkp1, fkp2
      integer        :: k, kp, kp0

!---------------------------------------------------------------------
!  local variables:
!
!    delp_lo
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
 
      do k=nklo,nkhi
        if (do_triangle) then
          kp0 = k + nkplo
        else
          kp0 = nkplo
        endif
        do kp=kp0,nkphi
          if (indx_press_hiv(kp,k) - indx_press_lov(kp,k) .GE. 3 .AND. &
              indx_press_hiv(kp,k) .LT. NSTDCO2LVLS - 1         ) then

!---------------------------------------------------------------------
!    use quadratic interpolation:
!---------------------------------------------------------------------
            delp_lo(kp,k) =             &
                 pressint_lov(kp,k) - pa(indx_press_lov(kp,k)+1)

!--------------------------------------------------------------------
!    1) for fixed (kp), varying (k)
!--------------------------------------------------------------------
            d1kp(kp,k) =   &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+2) -  &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1)  )/  &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =   &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -  &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  )  )/  &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/  &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*  &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp(kp,k) =   &
              error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) +  &
                delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!--------------------------------------------------------------------
!    2) for fixed (kp+1), varying (k)
!--------------------------------------------------------------------
            d1kp(kp,k) =   &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+2) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1)  )/&
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =   &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  )  )/&
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  )) 
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp1(kp,k) =     &
              error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) +   &
                delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!----------------------------------------------------------------------
!    3) for fixed (kp+2), varying (k)
!----------------------------------------------------------------------
            d1kp(kp,k) =     &
              (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+2) -  &
               error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1)  )/&
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp(kp,k) =     &
              (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) - &
               error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)  )  )/&
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp2(kp,k) =     &
              error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) +  &
                delp_lo(kp,k)*(akp(kp,k) + delp_lo(kp,k)*bkp(kp,k))

!---------------------------------------------------------------------
!    4) finally, varying (kp) using (fkp,fkp1,fkp2)
!---------------------------------------------------------------------
            delp_hi(kp,k) =     &
                 pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k)+1)
            d1kp(kp,k) =     &
              (fkp2(kp,k) - fkp1(kp,k)) /    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            d2kp(kp,k) =     &
              (fkp1(kp,k) - fkp (kp,k)) /    &
              (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)+0))
            bkp(kp,k) = (d1kp(kp,k) - d2kp(kp,k))/    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)  ))
            akp(kp,k) = d1kp(kp,k) - bkp(kp,k)*    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            errorint(kp,k) =     &
                             fkp1(kp,k) +     &
              delp_hi(kp,k)*(akp(kp,k) + delp_hi(kp,k)*bkp(kp,k))
 
          elseif (indx_press_hiv(kp,k) .GT. indx_press_lov(kp,k)) then

!--------------------------------------------------------------------
!    use linear interpolation:
!--------------------------------------------------------------------
            delp_lo(kp,k) =     &
                 pressint_lov(kp,k) - pa(indx_press_lov(kp,k))

!--------------------------------------------------------------------
!    1) for fixed (kp), varying (k)
!--------------------------------------------------------------------
            d2kp(kp,k) =     &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -    &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  )  ) / &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            fkp(kp,k) =     &
                  error(indx_press_hiv(kp,k),indx_press_lov(kp,k)) +  &
                              delp_lo(kp,k)*d2kp(kp,k)

!--------------------------------------------------------------------
!    2) for fixed (kp+1), varying (k)
!--------------------------------------------------------------------
            d2kp(kp,k) =     &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k) ))/  &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k) ))
            fkp1(kp,k) =     &
                  error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)) +&
                              delp_lo(kp,k)*d2kp(kp,k)

!--------------------------------------------------------------------
!    3) linear interpolate (fkp,fkp1):
!--------------------------------------------------------------------
            errorint(kp,k) =     &
              (fkp(kp,k)*    &
              (pa(indx_press_hiv(kp,k)+1) - pressint_hiv(kp,k)) +    &
               fkp1(kp,k)*    &
              (pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k))) ) /    &
              (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)))

          else

!---------------------------------------------------------------------
!    the error function for closely-spaced pressures equals zero
!    (section 3.2, Ref. (2))
!---------------------------------------------------------------------
            errorint(kp,k) = 0.0
          endif
        enddo
      enddo

!---------------------------------------------------------------------


end subroutine interp_error



!#####################################################################
! <SUBROUTINE NAME="interp_error_r">
!  <OVERVIEW>
!   Subroutine to examine error associated with interpolation onto
!   pressure grids.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to examine error associated with interpolation onto
!   pressure grids.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call interp_error_r (error, pressint_hiv, pressint_lov,  &
!                         indx_press_hiv, indx_press_lov,  &
!                         do_triangle, errorint)
!  </TEMPLATE>
!  <IN NAME="error" TYPE="real">
!   interpolation error at standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!  </IN>
!  <IN NAME="pressint_hiv, pressint_lov" TYPE="real">
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     pressint_lov = pressure of low (kp) interpolated pressure
!  </IN>
!  <IN NAME="indx_press_hiv, indx_press_lov" TYPE="real">
!   indx_press_hiv = pressure on std pa grid of high (kp) pressure
!   indx_press_lov = pressure on std pa grid of low  (kp) pressure
!  </IN>
!  <IN NAME="do_triangle" TYPE="logical">
!   state variable that determines the interpolation scheme
!  </IN>
!  <OUT NAME="errorint" TYPE="real">
!   error at interpolated grid
!  </OUT>
! </SUBROUTINE>
!
subroutine interp_error_r (error, pressint_hiv, pressint_lov,    &
                           indx_press_hiv, indx_press_lov,    &
                           do_triangle, errorint)

!-------------------------------------------------------------------
!
!-------------------------------------------------------------------

logical,                 intent(in)   :: do_triangle
real,    dimension(:,:), intent(in)   :: error, pressint_hiv,   &
                                         pressint_lov
integer, dimension(:,:), intent(in)   :: indx_press_hiv, indx_press_lov
real,    dimension(:,:), intent(out)  :: errorint

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    do_triangle
!     press_hiv = pressure on std pa grid of high (kp) pressure
!     pressint_hiv = pressure of high(kp) interpolated pressure
!     error = error at standard pa grid. evaluated on
!             a (NSTDCO2LVLS,NSTDCO2LVLS) grid when kp ge k).
!     errorint = error at interpolated grid
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:

      real, dimension(NSTDCO2LVLS) :: delp_lo, delp_hi, d1kp, d2kp, &
                                      bkp, akp, fkp, d1kp1, d2kp1,  &
                                      bkp1, akp1, fkp1, d1kp2,   &
                                      d2kp2, bkp2, akp2, fkp2,   &
                                      d1kpf, d2kpf, bkpf, akpf
      integer     :: k, kp, kp0
 
!-------------------------------------------------------------------
!   local variables:
!
!     delp_lo
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------

      do k=1,NSTDCO2LVLS
        if (do_triangle) then
          kp0 = k + 1
        else
          kp0 = 1
        endif
        do kp=kp0,NSTDCO2LVLS
          if (indx_press_hiv(kp,k) - indx_press_lov(kp,k) .GE. 3 .AND. &
              indx_press_hiv(kp,k) .LT. NSTDCO2LVLS - 1         ) then

!---------------------------------------------------------------------
!    use quadratic interpolation:
!---------------------------------------------------------------------
            delp_lo(kp) =     &
                 pressint_lov(kp,k) - pa(indx_press_lov(kp,k)+1)
 
!---------------------------------------------------------------------
!    1) for fixed (kp), varying (k)
!---------------------------------------------------------------------
            d1kp(kp) =     &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+2) -    &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1)  ) / &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1)) 
            d2kp(kp) =     &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -    &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  )  ) / &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp(kp) = (d1kp(kp) - d2kp(kp))/    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp(kp) = d1kp(kp) - bkp(kp)*    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp(kp) =     &
              error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) +    &
                delp_lo(kp)*(akp(kp) + delp_lo(kp)*bkp(kp))

!---------------------------------------------------------------------
!    2) for fixed (kp+1), varying (k)
!---------------------------------------------------------------------
            d1kp1(kp) =     &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+2) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1))/  &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp1(kp) =     &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k) ) )/  &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp1(kp) = (d1kp1(kp) - d2kp1(kp))/    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp1(kp) = d1kp1(kp) - bkp1(kp)*    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp1(kp) =     &
              error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) +  &
                delp_lo(kp)*(akp1(kp) + delp_lo(kp)*bkp1(kp))
 
!---------------------------------------------------------------------
!    3) for fixed (kp+2), varying (k)
!---------------------------------------------------------------------
            d1kp2(kp) =     &
              (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+2) -  &
               error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1)  )/&
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            d2kp2(kp) =     &
              (error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) - & 
               error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)  )  )/&
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            bkp2(kp) = (d1kp2(kp) - d2kp2(kp))/    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)  ))
            akp2(kp) = d1kp2(kp) - bkp2(kp)*    &
              (pa(indx_press_lov(kp,k)+2) - pa(indx_press_lov(kp,k)+1))
            fkp2(kp) =     &
              error(indx_press_hiv(kp,k)+2,indx_press_lov(kp,k)+1) +  &
                delp_lo(kp)*(akp2(kp) + delp_lo(kp)*bkp2(kp))
 
!---------------------------------------------------------------------
!    4) finally, varying (kp) using (fkp,fkp1,fkp2)
!---------------------------------------------------------------------
            delp_hi(kp) =     &
                 pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k)+1)
            d1kpf(kp) =     &
              (fkp2(kp) - fkp1(kp)) /    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            d2kpf(kp) =     &
              (fkp1(kp) - fkp (kp)) /    &
              (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)+0))
            bkpf(kp) = (d1kpf(kp) - d2kpf(kp))/    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)  ))
            akpf(kp) = d1kpf(kp) - bkpf(kp)*    &
              (pa(indx_press_hiv(kp,k)+2) - pa(indx_press_hiv(kp,k)+1))
            errorint(kp,k) =     &
                             fkp1(kp) +     &
              delp_hi(kp)*(akpf(kp) + delp_hi(kp)*bkpf(kp))

          elseif (indx_press_hiv(kp,k) .GT. indx_press_lov(kp,k)) then

!---------------------------------------------------------------------
!    use linear interpolation:
!---------------------------------------------------------------------
            delp_lo(kp) =     &
                 pressint_lov(kp,k) - pa(indx_press_lov(kp,k))
 
!---------------------------------------------------------------------
!    1) for fixed (kp), varying (k)
!---------------------------------------------------------------------
            d2kp(kp) =     &
              (error(indx_press_hiv(kp,k),indx_press_lov(kp,k)+1) -    &
               error(indx_press_hiv(kp,k),indx_press_lov(kp,k)  )  ) / &
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            fkp(kp) =   &
                  error(indx_press_hiv(kp,k),indx_press_lov(kp,k)) +  &
                              delp_lo(kp)*d2kp(kp)

!---------------------------------------------------------------------
!    2) for fixed (kp+1), varying (k)
!---------------------------------------------------------------------
            d2kp1(kp) =   &
              (error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)+1) -  &
               error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)  )  )/&
              (pa(indx_press_lov(kp,k)+1) - pa(indx_press_lov(kp,k)  ))
            fkp1(kp) =   &
                  error(indx_press_hiv(kp,k)+1,indx_press_lov(kp,k)) + &
                              delp_lo(kp)*d2kp1(kp)

!---------------------------------------------------------------------
!    3) linear interpolate (fkp,fkp1):
!---------------------------------------------------------------------
            errorint(kp,k) =   &
              (fkp(kp)*  &
              (pa(indx_press_hiv(kp,k)+1) - pressint_hiv(kp,k)) +  &
               fkp1(kp)*  &
              (pressint_hiv(kp,k) - pa(indx_press_hiv(kp,k))) ) /  &
              (pa(indx_press_hiv(kp,k)+1) - pa(indx_press_hiv(kp,k)))
          else
 
!---------------------------------------------------------------------
!    the error function for closely-spaced pressures equals zero
!    (section 3.2, Ref. (2))
!---------------------------------------------------------------------
            errorint(kp,k) = 0.0
          endif
        enddo
      enddo
 
!---------------------------------------------------------------------
 

end subroutine interp_error_r



!#####################################################################
! <SUBROUTINE NAME="pathv1">
!  <OVERVIEW>
!   Subroutine to compute the path function for the co2 interpolation pgm. 
!   between a pressure (press_lo) and a variable pressure (press_hi)
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute the path function for the co2 interpolation pgm. 
!   between a pressure (press_lo) and a variable pressure (press_hi)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call pathv1 (press_hi, press_lo, ndimlo, ndimhi, upath)
!  </TEMPLATE>
!  <IN NAME="press_hi, press_lo" TYPE="real">
!   The reference pressure levels
!  </IN>
!  <IN NAME="ndimlo, ndimhi" TYPE="integer">
!   the index of pressure level bound
!  </IN>
!  <OUT NAME="upath" TYPE="real">
!   The path function for the co2 interpolation pgm.
!  </OUT>
! </SUBROUTINE>
!
subroutine pathv1 (press_hi, press_lo, ndimlo, ndimhi, upath)
 
!--------------------------------------------------------------------
!    pathv1 computes the path function given in Eqs. (5) and (A5) in
!    Ref. (2) for the co2 interpolation pgm. between a 
!    pressure (press_lo) and a variable pressure (press_hi). This
!    has been modified on 5/27/97.
!--------------------------------------------------------------------
 
real,     dimension (:), intent(in)       :: press_hi, press_lo
real,     dimension (:), intent(out)      :: upath
integer,                 intent(in)       :: ndimlo, ndimhi

!-------------------------------------------------------------------
!   intent(in) variables:
!
!     press_hi
!
!--------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables

      integer      :: k   ! do-loop index

!------------------------------------------------------------------
!
!------------------------------------------------------------------
      do k = ndimlo,ndimhi

!------------------------------------------------------------------
!    all  a(**)b code replaced with exp(b*(log(a)) code below for 
!    overall ~ 10% speedup in standalone code -- no change in radiag 
!    file
!       upath(k) = (press_hi(k) - press_lo(k))**(1./sexp(k))*   &
!                  (press_hi(k) + press_lo(k) + dop_core)
!------------------------------------------------------------------

        upath(k) = EXP((1./sexp(k))*log((press_hi(k) - press_lo(k))))*&
                   (press_hi(k) + press_lo(k) + dop_core)
      enddo

!---------------------------------------------------------------------

 
end subroutine pathv1




!#####################################################################
! <SUBROUTINE NAME="rctrns">
!  <OVERVIEW>
!   Subroutine to compute co2 transmission functions for actual co2 
!   concentration
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute co2 transmission functions for actual co2 
!   concentration
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rctrns (gas_type, co2_std_lo, co2_std_hi, co2_vmr,  &
!                nf, nt, trns_vmr)
!  </TEMPLATE>
! </SUBROUTINE>
subroutine rctrns (gas_type, co2_std_lo, co2_std_hi, co2_vmr,  &
                   nf, nt, trns_vmr)

!-------------------------------------------------------------------
!    rctrns computes co2 transmission functions for actual co2 
!    concentration using method of section 5, Ref. (2).
!-------------------------------------------------------------------

character(len=*),        intent(in)    :: gas_type
integer,                 intent(in)    :: nf,nt
real,                    intent(in)    :: co2_vmr, co2_std_lo,   &
                                          co2_std_hi
real,    dimension(:,:), intent(inout) :: trns_vmr

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     gas_type
!      co2_std_hi = value of higher std co2 concentration in ppmv
!      co2_std_lo = value of lower std co2 concentration in ppmv
!      co2_vmr   = value of actual co2 concentration in ppmv
!
!-------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension(NSTDCO2LVLS,NSTDCO2LVLS) :: &
                                           approx_guess1,    &
                                           approxint_guess1, &
                                           approxint_guess2, &
                                           error_guess1,     &
                                           errorint_guess1,  &
                                           errorint_guess2,  &
                                           trans_guess1,     &
                                           trans_guess2
      real, dimension(NSTDCO2LVLS,NSTDCO2LVLS) :: &
                                           caintv, uexpintv, &
                                           sexpintv, xaintv, &
                                           press_hiv, press_lov
      logical do_triangle

!--------------------------------------------------------------------
!  local variables:
!
!     approx_guess1
!
!---------------------------------------------------------------------
     integer :: k, kp
!----------------------------------------------------------------------
!    the first part of the method is to obtain a first guess co2
!    transmission function for the desired concentration using only the
!    co2 tf's for the higher standard concentration.
!----------------------------------------------------------------------
      call coeint (gas_type, nf, trns_std_hi, ca, sexp, xa, uexp)
 
!-------------------------------------------------------------------
!    compute the interpolation. 
!-------------------------------------------------------------------
      do_triangle = .true.

!--------------------------------------------------------------------
!    1) compute approx function at standard (pa) pressures
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        press_hi(k) = pa(k)
        press_lo(k) = pa(k)
      enddo
 
!-------------------------------------------------------------------
!    compute the 2-d input arrays
!-------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        do kp=k,NSTDCO2LVLS
          press_hiv(kp,k) = pa(kp)
          press_lov(kp,k) = pa(k)
          caintv(kp,k) = ca(kp)
          sexpintv(kp,k) = sexp(kp)
          xaintv(kp,k) = xa(kp)
          uexpintv(kp,k) = uexp(kp)
        enddo
      enddo

!-------------------------------------------------------------------
!    the call (and calculations) to pathv2_std has been subsumed into
!    the subroutine approx_fn_std
!-------------------------------------------------------------------
      call approx_fn_std (press_hiv, press_lov, do_triangle, &
                          caintv, sexpintv, xaintv, uexpintv,  &
                          approx_guess1)

!--------------------------------------------------------------------
!    2) compute error function at standard (pa) pressures
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        do kp=k+1,NSTDCO2LVLS
          error_guess1(kp,k) = 1.0 - trns_std_hi(kp,k) -  &
                               approx_guess1(kp,k)
        enddo
        error_guess1(k,k) = 0.0
      enddo
        
!---------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!---------------------------------------------------------------------
      !if (nf .EQ. 1 .AND. nt .EQ. 1) then
        do k=1,NSTDCO2LVLS
          do kp=k+1,NSTDCO2LVLS
            pressint_hiv_std_pt1(kp,k) = ((co2_vmr+co2_std_hi)*pa(kp) +&
                                         (co2_std_hi-co2_vmr)*pa(k))  /&
                                         (2.*co2_std_hi)
            pressint_lov_std_pt1(kp,k) = ((co2_std_hi-co2_vmr)*pa(kp) +&
                                         (co2_vmr+co2_std_hi)*pa(k))  /&
                                         (2.*co2_std_hi)
          enddo
        enddo
      !endif
      call intcoef_2d_std (pressint_hiv_std_pt1, pressint_lov_std_pt1, &
                           nf, nt, do_triangle,  &
                           indx_pressint_hiv_std_pt1,   &
                           indx_pressint_lov_std_pt1,  &
                           caintv, sexpintv, xaintv,uexpintv)

!----------------------------------------------------------------------
!    4) interpolate error function to (pressint_hiv, pressint_lov)
!    for all (k,k')
!----------------------------------------------------------------------
      call interp_error_r (error_guess1, pressint_hiv_std_pt1,  &
                           pressint_lov_std_pt1,  &
                           indx_pressint_hiv_std_pt1,   &
                           indx_pressint_lov_std_pt1, do_triangle,  &
                           errorint_guess1)

!---------------------------------------------------------------------
!    5) compute approx function for (pressint_hiv, pressint_lov)
!---------------------------------------------------------------------
      
!--------------------------------------------------------------------
!    the call (and calculations) to pathv2_std has been subsumed into
!    the subroutine approx_fn_std
!--------------------------------------------------------------------
      call approx_fn_std (pressint_hiv_std_pt1, pressint_lov_std_pt1,  &
                          do_triangle, caintv, sexpintv, xaintv,   &
                          uexpintv, approxint_guess1)

!---------------------------------------------------------------------
!    6) compute first guess transmission function using Eq.(3),
!    Ref.(2).
!---------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        do kp=k+1,NSTDCO2LVLS
          trans_guess1(kp,k) = 1.0 -  &
                       (errorint_guess1(kp,k) + approxint_guess1(kp,k))
        enddo
      enddo
 
!---------------------------------------------------------------------
!    the second part of the method is to obtain a second guess co2
!    transmission function for the lower standard  concentration using
!    only the co2 tf's for the higher standard concentration.
!    the coeint call and steps (1-2) of part (1) need not be repeated.
!---------------------------------------------------------------------
  if (co2_std_lo .GT. 0.0 .or. do_co2_bug) then
 !---------------------------------------------------------------------
!    3) derive the pressures for interpolation using Eqs. (8a-b)
!       in Ref.(2).
!---------------------------------------------------------------------
      !if (nf .EQ. 1 .AND. nt .EQ. 1) then
        do k=1,NSTDCO2LVLS
          do kp=k+1,NSTDCO2LVLS
            pressint_hiv_std_pt2(kp,k) = ((co2_std_lo+co2_std_hi)*  &
                                           pa(kp) +  &
                                         (co2_std_hi-co2_std_lo)*  &
                                           pa(k))/  &
                                          (2.*co2_std_hi)
            pressint_lov_std_pt2(kp,k) = ((co2_std_hi-co2_std_lo)* &
                                           pa(kp) +  &
                                         (co2_std_lo+co2_std_hi)* &
                                           pa(k))/  &
                                         (2.*co2_std_hi)
          enddo
        enddo
      !endif
      call intcoef_2d_std (pressint_hiv_std_pt2, pressint_lov_std_pt2, &
                           nf, nt, do_triangle,  &
                           indx_pressint_hiv_std_pt2,    &
                           indx_pressint_lov_std_pt2,  &
                           caintv, sexpintv, xaintv,uexpintv)

!---------------------------------------------------------------------
!    4) interpolate error function to (pressint_hiv, pressint_lov)
!       for all (k,k')
!---------------------------------------------------------------------
      call interp_error_r (error_guess1, pressint_hiv_std_pt2,   &
                           pressint_lov_std_pt2,  &
                           indx_pressint_hiv_std_pt2,   &
                           indx_pressint_lov_std_pt2,  &
                           do_triangle,  errorint_guess2)

!---------------------------------------------------------------------
!    5) compute approx function for (pressint_hiv, pressint_lov)
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    the call (and calculations) to pathv2_std has been subsumed into
!    the subroutine approx_fn_std
!---------------------------------------------------------------------
      call approx_fn_std (pressint_hiv_std_pt2, pressint_lov_std_pt2,  &
                          do_triangle, caintv, sexpintv, xaintv,   &
                          uexpintv, approxint_guess2)
 
!--------------------------------------------------------------------
!    6) compute second guess transmission function using Eq.(3),
!       Ref.(2).
!--------------------------------------------------------------------
      do k=1,NSTDCO2LVLS 
        do kp=k+1,NSTDCO2LVLS
          trans_guess2(kp,k) = 1.0 -  &
            (errorint_guess2(kp,k) + approxint_guess2(kp,k))
        enddo
      enddo

!---------------------------------------------------------------------
!    finally, obtain transmission function for (co2_vmr) using
!    Eq.(9), Ref. (2).
!---------------------------------------------------------------------
      do k=1,NSTDCO2LVLS
        do kp=k+1,NSTDCO2LVLS
          trns_vmr(kp,k) = trans_guess1(kp,k) +  &
                           (co2_std_hi - co2_vmr)/  &
                           (co2_std_hi - co2_std_lo)*  &
                          (trns_std_lo(kp,k) - trans_guess2(kp,k))
          trns_vmr(k,kp) = trns_vmr(kp,k)
        enddo
        trns_vmr(k,k) = 1.0
      enddo

    else
      do k=1,NSTDCO2LVLS
        do kp=k+1,NSTDCO2LVLS
          trns_vmr(kp,k) = trans_guess1(kp,k)
          trns_vmr(k,kp) = trns_vmr(kp,k)
        enddo
          trns_vmr(k,k) = 1.0
      enddo
    
    endif

!---------------------------------------------------------------------
       

end subroutine rctrns



!#####################################################################
! <SUBROUTINE NAME="read_lbltfs">
!  <OVERVIEW>
!   Subroutine to read gas transmission functions from input file
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to read gas transmission functions from input file
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_lbltfs (gas_type, callrctrns, nstd_lo, nstd_hi, nf,   &
!                     ntbnd, trns_std_hi_nf, trns_std_lo_nf )
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_lbltfs (gas_type, callrctrns, nstd_lo, nstd_hi, nf,   &
                        ntbnd, trns_std_hi_nf, trns_std_lo_nf )
 
use cc_mpi
use filnames_m
use infile
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),           intent(in)   :: gas_type
logical,                    intent(in)   :: callrctrns
integer,                    intent(in)   :: nstd_lo, nstd_hi, nf
integer, dimension(:),      intent(in)   :: ntbnd
real,    dimension(NSTDCO2LVLS,NSTDCO2LVLS,3),  intent(out)  :: &
                             trns_std_hi_nf, trns_std_lo_nf

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     gas_type
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      character(len=80) input_lblco2name(nfreq_bands_sea_co2,number_std_co2_vmrs)
      character(len=80) input_lblch4name(nfreq_bands_sea_ch4,number_std_ch4_vmrs)
      character(len=80) input_lbln2oname(nfreq_bands_sea_n2o,number_std_n2o_vmrs)
      character(len=80) name_lo
      character(len=80) name_hi
      character(len=1024) filename, ncname

      integer        :: n, nt, ierr
      
      integer, dimension(3) :: startpos, npos
      integer ncid, ncstatus, varid
      logical tst
 
      data (input_lblco2name(n,1),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_0_HITRAN2012_490850_495lyr      ', 'cns_co2_0_HITRAN2012_490630_495lyr      ', &
        'cns_co2_0_HITRAN2012_630700_495lyr      ', 'cns_co2_0_HITRAN2012_700850_495lyr      ', &
        'cns_co2_0_43um                          ', 'cns_co2_0_HITRAN2012_9901070_495lyr     ', &
        'cns_co2_0_HITRAN2012_900990_495lyr      ', 'cns_co2_0_HITRAN2012_10701200_495lyr    '  /
      data (input_lblco2name(n,2),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_165_HITRAN2012_490850_495lyr    ', 'cns_co2_165_HITRAN2012_490630_495lyr    ', &
        'cns_co2_165_HITRAN2012_630700_495lyr    ', 'cns_co2_165_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_165_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_165_HITRAN2012_900990_495lyr    ', 'cns_co2_165_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,3),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_300_HITRAN2012_490850_495lyr    ', 'cns_co2_300_HITRAN2012_490630_495lyr    ', &
        'cns_co2_300_HITRAN2012_630700_495lyr    ', 'cns_co2_300_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_300_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_300_HITRAN2012_900990_495lyr    ', 'cns_co2_300_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,4),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_330_HITRAN2012_490850_495lyr    ', 'cns_co2_330_HITRAN2012_490630_495lyr    ', &
        'cns_co2_330_HITRAN2012_630700_495lyr    ', 'cns_co2_330_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_330_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_330_HITRAN2012_900990_495lyr    ', 'cns_co2_330_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,5),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_348_HITRAN2012_490850_495lyr    ', 'cns_co2_348_HITRAN2012_490630_495lyr    ', &
        'cns_co2_348_HITRAN2012_630700_495lyr    ', 'cns_co2_348_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_348_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_348_HITRAN2012_900990_495lyr    ', 'cns_co2_348_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,6),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_356_HITRAN2012_490850_495lyr    ', 'cns_co2_356_HITRAN2012_490630_495lyr    ', &
        'cns_co2_356_HITRAN2012_630700_495lyr    ', 'cns_co2_356_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_356_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_356_HITRAN2012_900990_495lyr    ', 'cns_co2_356_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,7),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_360_HITRAN2012_490850_495lyr    ', 'cns_co2_360_HITRAN2012_490630_495lyr    ', &
        'cns_co2_360_HITRAN2012_630700_495lyr    ', 'cns_co2_360_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_360_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_360_HITRAN2012_900990_495lyr    ', 'cns_co2_360_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,8),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_600_HITRAN2012_490850_495lyr    ', 'cns_co2_600_HITRAN2012_490630_495lyr    ', &
        'cns_co2_600_HITRAN2012_630700_495lyr    ', 'cns_co2_600_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_600_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_600_HITRAN2012_900990_495lyr    ', 'cns_co2_600_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,9),n=1,nfreq_bands_sea_co2)/            &
        'cns_co2_660_HITRAN2012_490850_495lyr    ', 'cns_co2_660_HITRAN2012_490630_495lyr    ', &
        'cns_co2_660_HITRAN2012_630700_495lyr    ', 'cns_co2_660_HITRAN2012_700850_495lyr    ', &
        'cnsco2_0_43um                           ', 'cns_co2_660_HITRAN2012_9901070_495lyr   ', &
        'cns_co2_660_HITRAN2012_900990_495lyr    ', 'cns_co2_660_HITRAN2012_10701200_495lyr  '  /
      data (input_lblco2name(n,10),n=1,nfreq_bands_sea_co2)/           &
        'cns_co2_1320_HITRAN2012_490850_495lyr   ', 'cns_co2_1320_HITRAN2012_490630_495lyr   ', &
        'cns_co2_1320_HITRAN2012_630700_495lyr   ', 'cns_co2_1320_HITRAN2012_700850_495lyr   ', &
        'cnsco2_0_43um                           ', 'cns_co2_1320_HITRAN2012_9901070_495lyr  ', &
        'cns_co2_1320_HITRAN2012_900990_495lyr   ', 'cns_co2_1320_HITRAN2012_10701200_495lyr '  /
      data (input_lblco2name(n,11),n=1,nfreq_bands_sea_co2)/           &
        'cns_co2_1600_HITRAN2012_490850_495lyr   ', 'cns_co2_1600_HITRAN2012_490630_495lyr   ', &
        'cns_co2_1600_HITRAN2012_630700_495lyr   ', 'cns_co2_1600_HITRAN2012_700850_495lyr   ', &
        'cnsco2_0_43um                           ', 'cns_co2_1600_HITRAN2012_9901070_495lyr  ', &
        'cns_co2_1600_HITRAN2012_900990_495lyr   ', 'cns_co2_1600_HITRAN2012_10701200_495lyr '  /
      data (input_lblco2name(n,12),n=1,nfreq_bands_sea_co2)/           &
        'cns_co2_2000_HITRAN2012_490850_495lyr   ', 'cns_co2_2000_HITRAN2012_490630_495lyr   ', &
        'cns_co2_2000_HITRAN2012_630700_495lyr   ', 'cns_co2_2000_HITRAN2012_700850_495lyr   ', &
        'cnsco2_0_43um                           ', 'cns_co2_2000_HITRAN2012_9901070_495lyr  ', &
        'cns_co2_2000_HITRAN2012_900990_495lyr   ', 'cns_co2_2000_HITRAN2012_10701200_495lyr '  /
      data (input_lblco2name(n,13),n=1,nfreq_bands_sea_co2)/           &
        'cns_co2_5000_HITRAN2012_490850_495lyr   ', 'cns_co2_5000_HITRAN2012_490630_495lyr   ', &
        'cns_co2_5000_HITRAN2012_630700_495lyr   ', 'cns_co2_5000_HITRAN2012_700850_495lyr   ', &
        'cnsco2_0_43um                           ', 'cns_co2_5000_HITRAN2012_9901070_495lyr  ', &
        'cns_co2_5000_HITRAN2012_900990_495lyr   ', 'cns_co2_5000_HITRAN2012_10701200_495lyr '  /
      data (input_lblco2name(n,14),n=1,nfreq_bands_sea_co2)/           &
        'cns_co2_10000_HITRAN2012_490850_495lyr  ', 'cns_co2_10000_HITRAN2012_490630_495lyr  ', &
        'cns_co2_10000_HITRAN2012_630700_495lyr  ', 'cns_co2_10000_HITRAN2012_700850_495lyr  ', &
        'cnsco2_0_43um                           ', 'cns_co2_10000_HITRAN2012_9901070_495lyr ', &
        'cns_co2_10000_HITRAN2012_900990_495lyr  ', 'cns_co2_10000_HITRAN2012_10701200_495lyr'  /

 
      data (input_lblch4name(n,1),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_0_HITRAN2012_12001400_495lyr   '/
      data (input_lblch4name(n,2),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_300_HITRAN2012_12001400_495lyr '/
      data (input_lblch4name(n,3),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_700_HITRAN2012_12001400_495lyr '/
      data (input_lblch4name(n,4),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_1250_HITRAN2012_12001400_495lyr'/
      data (input_lblch4name(n,5),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_1750_HITRAN2012_12001400_495lyr'/
      data (input_lblch4name(n,6),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_2250_HITRAN2012_12001400_495lyr'/
      data (input_lblch4name(n,7),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_2800_HITRAN2012_12001400_495lyr'/
      data (input_lblch4name(n,8),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_4000_HITRAN2012_12001400_495lyr'/
      data (input_lblch4name(n,9),n=1,nfreq_bands_sea_ch4)/          &
        'cns_ch4_6000_HITRAN2012_12001400_495lyr'/
 
      data (input_lbln2oname(n,1),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_0_HITRAN2012_12001400_495lyr  ', 'cns_n2o_0_HITRAN2012_10701200_495lyr  ', &
        'cns_n2o_0_HITRAN2012_560630_495lyr    '/
      data (input_lbln2oname(n,2),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_180_HITRAN2012_12001400_495lyr', 'cns_n2o_180_HITRAN2012_10701200_495lyr', &
        'cns_n2o_180_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,3),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_275_HITRAN2012_12001400_495lyr', 'cns_n2o_275_HITRAN2012_10701200_495lyr', &
        'cns_n2o_275_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,4),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_310_HITRAN2012_12001400_495lyr', 'cns_n2o_310_HITRAN2012_10701200_495lyr', &
        'cns_n2o_310_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,5),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_340_HITRAN2012_12001400_495lyr', 'cns_n2o_340_HITRAN2012_10701200_495lyr', &
        'cns_n2o_340_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,6),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_375_HITRAN2012_12001400_495lyr', 'cns_n2o_375_HITRAN2012_10701200_495lyr', &
        'cns_n2o_375_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,7),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_500_HITRAN2012_12001400_495lyr', 'cns_n2o_500_HITRAN2012_10701200_495lyr', &
        'cns_n2o_500_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,8),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_600_HITRAN2012_12001400_495lyr', 'cns_n2o_600_HITRAN2012_10701200_495lyr', &
        'cns_n2o_600_HITRAN2012_560630_495lyr  '/
      data (input_lbln2oname(n,9),n=1,nfreq_bands_sea_n2o)/           &
        'cns_n2o_800_HITRAN2012_12001400_495lyr', 'cns_n2o_800_HITRAN2012_10701200_495lyr', &
        'cns_n2o_800_HITRAN2012_560630_495lyr  '/

!--------------------------------------------------------------------
!  local variables:
!
!     input_lblco2name    
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (gas_type .EQ. 'co2') then
        name_lo = input_lblco2name(nf,nstd_lo)
        name_hi = input_lblco2name(nf,nstd_hi)
      endif
      if (gas_type .EQ. 'ch4') then
        name_lo = input_lblch4name(nf,nstd_lo)
        name_hi = input_lblch4name(nf,nstd_hi)
      endif
      if (gas_type .EQ. 'n2o') then
        name_lo = input_lbln2oname(nf,nstd_lo)
        name_hi = input_lbln2oname(nf,nstd_hi)
      endif

!-------------------------------------------------------------------
!    read in tfs of higher std gas concentration
!-------------------------------------------------------------------

      filename = trim(cnsdir) // '/' // trim(name_hi)
      ncname = trim(filename) // '.nc'

      startpos(:) = 1
      npos(1) = NSTDCO2LVLS
      npos(2) = NSTDCO2LVLS
      npos(3) = ntbnd(nf)
      call ccnf_open(ncname,ncid,ncstatus)
      if ( ncstatus==0 ) then
        !write(6,*) "Reading ",trim(ncname)
        call ccnf_inq_varid(ncid,"trns_std_nf",varid,tst)
        if ( tst ) then
          write(6,*) "trns_std_nf not found in ",trim(ncname)
          call ccmpi_abort(-1)
        end if
        call ccnf_get_vara(ncid,varid,startpos,npos,trns_std_hi_nf(:,:,1:ntbnd(nf)))
        call ccnf_close(ncid)
      else
        !write(6,*) "Reading ",trim(filename)  
        open(11,file=filename,access='DIRECT',recl=NSTDCO2LVLS*NSTDCO2LVLS*8,form='UNFORMATTED',action='READ',iostat=ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot locate ",trim(ncname)," or ",trim(filename)
          call ccmpi_abort(-1)
        end if    
        do nt = 1,ntbnd(nf)
          read(11,rec=nt) trns_std_hi_nf(:,:,nt)  
        end do
        close(11)
      end if    
      
      if ( callrctrns ) then
          
!--------------------------------------------------------------------
!    if necessary, read in tfs of lower standard gas concentration
!-------------------------------------------------------------------
        filename = trim(cnsdir) // '/' // trim(name_lo )
        ncname = trim(filename) // '.nc'

        startpos(:) = 1
        npos(1) = NSTDCO2LVLS
        npos(2) = NSTDCO2LVLS
        npos(3) = ntbnd(nf)
        call ccnf_open(ncname,ncid,ncstatus)
        if ( ncstatus == 0 ) then ! Netcdf file
          !write(6,*) "Reading ",trim(ncname)
          call ccnf_inq_varid(ncid,"trns_std_nf",varid,tst)
          if ( tst ) then
            write(6,*) "trns_std_nf not found in ",trim(ncname)
            call ccmpi_abort(-1)
          end if
          call ccnf_get_vara(ncid,varid,startpos,npos,trns_std_lo_nf(:,:,1:ntbnd(nf)))
          call ccnf_close(ncid)
        else
          !write(6,*) "Reading ",trim(filename)    
          open(11,file=filename,access='DIRECT',recl=NSTDCO2LVLS*NSTDCO2LVLS*8,form='UNFORMATTED',action='READ',iostat=ierr)
          if ( ierr/=0 ) then
            write(6,*) "ERROR: Cannot locate ",trim(ncname)," or ",trim(filename)
            call ccmpi_abort(-1)
          end if           
          do nt = 1,ntbnd(nf)
            read(11,rec=nt) trns_std_lo_nf(:,:,nt)  
          end do
          close(11)
        end if    
          
      end if
      
!--------------------------------------------------------------------


end subroutine read_lbltfs

subroutine read_lbltfs_old (gas_type, callrctrns, nstd_lo, nstd_hi, nf,   &
                        ntbnd, trns_std_hi_nf, trns_std_lo_nf )
 
use cc_mpi
use filnames_m
use infile
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),           intent(in)   :: gas_type
logical,                    intent(in)   :: callrctrns
integer,                    intent(in)   :: nstd_lo, nstd_hi, nf
integer, dimension(:),      intent(in)   :: ntbnd
real,    dimension(NSTDCO2LVLS,NSTDCO2LVLS,3),  intent(out)  :: &
                             trns_std_hi_nf, trns_std_lo_nf

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     gas_type
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      character(len=24) input_lblco2name(5,11)
      character(len=24) input_lblch4name(nfreq_bands_sea_ch4,8)
      character(len=24) input_lbln2oname(nfreq_bands_sea_n2o,7)
      character(len=24) name_lo
      character(len=24) name_hi
      character(len=110) filename, ncname

      integer        :: n
      
      integer, dimension(3) :: startpos, npos
      integer ncid, ncstatus, varid, ierr, nt
      logical tst
 
      data (input_lblco2name(n,1),n=1,5)/            &
        'cns_0_490850   ', 'cns_0_490630   ', 'cns_0_630700   ',       &
        'cns_0_700850   ', 'cns_0_43um     '/
      data (input_lblco2name(n,2),n=1,5)/            &
        'cns_165_490850   ', 'cns_165_490630   ', 'cns_165_630700   ', &
        'cns_165_700850   ', 'cns_165_43um     '/
      data (input_lblco2name(n,3),n=1,5)/            &
        'cns_300_490850   ', 'cns_300_490630   ', 'cns_300_630700   ', &
        'cns_300_700850   ', 'cns_300_43um     '/
      data (input_lblco2name(n,4),n=1,5)/            &
        'cns_330_490850   ', 'cns_330_490630   ', 'cns_330_630700   ', &
        'cns_330_700850   ', 'cns_330_43um     '/
      data (input_lblco2name(n,5),n=1,5)/            &
        'cns_348_490850   ', 'cns_348_490630   ', 'cns_348_630700   ', &
        'cns_348_700850   ', 'cns_348_43um     '/
      data (input_lblco2name(n,6),n=1,5)/            &
        'cns_356_490850   ', 'cns_356_490630   ', 'cns_356_630700   ', &
        'cns_356_700850   ', 'cns_356_43um     '/
      data (input_lblco2name(n,7),n=1,5)/            &
        'cns_360_490850   ', 'cns_360_490630   ', 'cns_360_630700   ', &
        'cns_360_700850   ', 'cns_360_43um     '/
      data (input_lblco2name(n,8),n=1,5)/            &
        'cns_600_490850   ', 'cns_600_490630   ', 'cns_600_630700   ', &
        'cns_600_700850   ', 'cns_600_43um     '/
      data (input_lblco2name(n,9),n=1,5)/            &
        'cns_660_490850   ', 'cns_660_490630   ', 'cns_660_630700   ', &
        'cns_660_700850   ', 'cns_660_43um     '/
      data (input_lblco2name(n,10),n=1,5)/           &
        'cns_1320_490850  ', 'cns_1320_490630  ', 'cns_1320_630700  ', &
        'cns_1320_700850  ', 'cns_1320_43um    '/
      data (input_lblco2name(n,11),n=1,5)/           &
        'cns_1600_490850  ', 'cns_1600_490630  ', 'cns_1600_630700  ', &
        'cns_1600_700850  ', 'cns_1600_43um    '/
 
      data (input_lblch4name(n,1),n=1,nfreq_bands_sea_ch4)/          &
        'cns_0_12001400'/
      data (input_lblch4name(n,2),n=1,nfreq_bands_sea_ch4)/          &
        'cns_300_12001400'/
      data (input_lblch4name(n,3),n=1,nfreq_bands_sea_ch4)/          &
        'cns_700_12001400'/
      data (input_lblch4name(n,4),n=1,nfreq_bands_sea_ch4)/          &
        'cns_1250_12001400'/
      data (input_lblch4name(n,5),n=1,nfreq_bands_sea_ch4)/          &
        'cns_1750_12001400'/
      data (input_lblch4name(n,6),n=1,nfreq_bands_sea_ch4)/          &
        'cns_2250_12001400'/
      data (input_lblch4name(n,7),n=1,nfreq_bands_sea_ch4)/          &
        'cns_2800_12001400'/
      data (input_lblch4name(n,8),n=1,nfreq_bands_sea_ch4)/          &
        'cns_4000_12001400'/
 
      data (input_lbln2oname(n,1),n=1,nfreq_bands_sea_n2o)/           &
        'cns_0_12001400 ', 'cns_0_10701200 ', 'cns_0_560630   '/
      data (input_lbln2oname(n,2),n=1,nfreq_bands_sea_n2o)/           &
        'cns_180_12001400 ', 'cns_180_10701200 ', 'cns_180_560630   '/
      data (input_lbln2oname(n,3),n=1,nfreq_bands_sea_n2o)/           &
        'cns_275_12001400 ', 'cns_275_10701200 ', 'cns_275_560630   '/
      data (input_lbln2oname(n,4),n=1,nfreq_bands_sea_n2o)/           &
        'cns_310_12001400 ', 'cns_310_10701200 ', 'cns_310_560630   '/
      data (input_lbln2oname(n,5),n=1,nfreq_bands_sea_n2o)/           &
        'cns_340_12001400 ', 'cns_340_10701200 ', 'cns_340_560630   '/
      data (input_lbln2oname(n,6),n=1,nfreq_bands_sea_n2o)/           &
        'cns_375_12001400 ', 'cns_375_10701200 ', 'cns_375_560630   '/
      data (input_lbln2oname(n,7),n=1,nfreq_bands_sea_n2o)/           &
        'cns_500_12001400 ', 'cns_500_10701200 ', 'cns_500_560630   '/

!--------------------------------------------------------------------
!  local variables:
!
!     input_lblco2name    
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (gas_type .EQ. 'co2') then
        name_lo = input_lblco2name(nf,nstd_lo)
        name_hi = input_lblco2name(nf,nstd_hi)
      endif
      if (gas_type .EQ. 'ch4') then
        name_lo = input_lblch4name(nf,nstd_lo)
        name_hi = input_lblch4name(nf,nstd_hi)
      endif
      if (gas_type .EQ. 'n2o') then
        name_lo = input_lbln2oname(nf,nstd_lo)
        name_hi = input_lbln2oname(nf,nstd_hi)
      endif

!-------------------------------------------------------------------
!    read in tfs of higher std gas concentration
!-------------------------------------------------------------------

      filename = trim(cnsdir) // '/' // trim(name_hi)
      ncname = trim(filename) // '.nc'

      startpos(:) = 1
      npos(1) = NSTDCO2LVLS
      npos(2) = NSTDCO2LVLS
      npos(3) = ntbnd(nf)
      call ccnf_open(ncname,ncid,ncstatus)
      if ( ncstatus == 0 ) then
        !write(6,*) "Reading ",trim(ncname)
        call ccnf_inq_varid(ncid,"trns_std_nf",varid,tst)
        if ( tst ) then
          write(6,*) "trns_std_nf not found in ",trim(ncname)
          call ccmpi_abort(-1)
        end if
        call ccnf_get_vara(ncid,varid,startpos,npos,trns_std_hi_nf(:,:,1:ntbnd(nf)))
        call ccnf_close(ncid)
      else
        !write(6,*) "Reading ",trim(filename)  
        open(11,file=filename,access='DIRECT',recl=NSTDCO2LVLS*NSTDCO2LVLS*8,form='UNFORMATTED',action='READ',iostat=ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot locate ",trim(ncname)," or ",trim(filename)
          call ccmpi_abort(-1)
        end if    
        do nt = 1,ntbnd(nf)
          read(11,rec=nt) trns_std_hi_nf(:,:,nt)  
        end do
        close(11)
      end if            
      
      if ( callrctrns ) then
          
!--------------------------------------------------------------------
!    if necessary, read in tfs of lower standard gas concentration
!-------------------------------------------------------------------
        filename = trim(cnsdir) // '/' // trim(name_lo )
        ncname = trim(filename) // '.nc'

        startpos(:) = 1
        npos(1) = NSTDCO2LVLS
        npos(2) = NSTDCO2LVLS
        npos(3) = ntbnd(nf)
        call ccnf_open(ncname,ncid,ncstatus)
        if ( ncstatus == 0 ) then
          !write(6,*) "Reading ",trim(ncname)
          call ccnf_inq_varid(ncid,"trns_std_nf",varid,tst)
          if ( tst ) then
            write(6,*) "trns_std_nf not found in ",trim(ncname)
            call ccmpi_abort(-1)
          end if
          call ccnf_get_vara(ncid,varid,startpos,npos,trns_std_lo_nf(:,:,1:ntbnd(nf)))
          call ccnf_close(ncid)
      else
        !write(6,*) "Reading ",trim(filename)  
        open(11,file=filename,access='DIRECT',recl=NSTDCO2LVLS*NSTDCO2LVLS*8,form='UNFORMATTED',action='READ',iostat=ierr)
        if ( ierr/=0 ) then
          write(6,*) "ERROR: Cannot locate ",trim(ncname)," or ",trim(filename)
          call ccmpi_abort(-1)
        end if    
        do nt = 1,ntbnd(nf)
          read(11,rec=nt) trns_std_lo_nf(:,:,nt)  
        end do
        close(11)
      end if  
          
      end if
      
!--------------------------------------------------------------------


end subroutine read_lbltfs_old


!#####################################################################
! <SUBROUTINE NAME="allocate_interp_arrays">
!  <OVERVIEW>
!   Subroutine to allocate interpolation arrays
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to allocate interpolation arrays
!  </DESCRIPTION>
!  <TEMPLATE>
!   call allocate_interp_arrays
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine allocate_interp_arrays

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

      allocate (pressint_hiv_std_pt1(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (pressint_lov_std_pt1(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (pressint_hiv_std_pt2(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (pressint_lov_std_pt2(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (indx_pressint_hiv_std_pt1(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (indx_pressint_lov_std_pt1(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (indx_pressint_hiv_std_pt2(NSTDCO2LVLS,NSTDCO2LVLS))
      allocate (indx_pressint_lov_std_pt2(NSTDCO2LVLS,NSTDCO2LVLS))

!-------------------------------------------------------------------


end subroutine allocate_interp_arrays


!####################################################################

subroutine deallocate_interp_arrays

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

      deallocate (pressint_hiv_std_pt1)
      deallocate (pressint_lov_std_pt1)
      deallocate (pressint_hiv_std_pt2)
      deallocate (pressint_lov_std_pt2)
      deallocate (indx_pressint_hiv_std_pt1)
      deallocate (indx_pressint_lov_std_pt1)
      deallocate (indx_pressint_hiv_std_pt2)
      deallocate (indx_pressint_lov_std_pt2)

!-------------------------------------------------------------------


end subroutine deallocate_interp_arrays


!####################################################################



             end module lw_gases_stdtf_mod


