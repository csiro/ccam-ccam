                 module esfsw_parameters_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to initialize shortwave parameters and access flux data.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes shortwave radiation calculation parameters such as
!  solar flux input at top of the atmosphere, number of shortwave bands
!  depending on the spectral resolution used, number of frequency points
!  in the gaussian quadrature algorithm, the number of streams used in
!  multiple stream flux algorithm, and the number of water vapor bands.
!
!  The code also provides two access methods: get and put solar flux data
!  
! </DESCRIPTION>

!    shared modules:

!use fms_mod,           only: open_namelist_file, fms_init, &
!                             mpp_pe, mpp_root_pe, stdlog, &
!                             file_exist, write_version_number, &
!                             check_nml_error, error_mesg, &
!                             FATAL, NOTE, WARNING, close_file

!  shared radiation package modules:

use rad_utilities_mod, only:  solar_spectrum_type, rad_utilities_init


!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!     esfsw_parameters_mod defines parameters for esf shortwave code,
!     including a description of the band structure  used to define the
!     solar spectrum.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: esfsw_parameters.F90,v 18.0.2.1 2010/08/30 20:33:32 wfc Exp $'
character(len=128)  :: tagname =  '$Name: testing $'

!--------------------------------------------------------------------
!----- interfaces ------

public       &
         esfsw_parameters_init, &
         esfsw_parameters_end,  &
         sw_resolution,         &
         sw_diff_streams

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16), save :: sw_resolution = 'low' ! either 'high' or 'low'
integer, save           :: sw_diff_streams = 1   ! number of streams of
                                                 ! diffuse radiation that
                                                 ! are considered (1 or 4)

!-------------------------------------------------------------------
!----- public data --------

!---------------------------------------------------------------------
!    TOT_WVNUMS     number of wavenumbers included in the parameter-
!                   ization of the solar spectrum 
!    Solar_spect    solar_spectrum_type variable defining the nature
!                   of the solar spectral paramaterization
!---------------------------------------------------------------------
integer, parameter                      :: TOT_WVNUMS  = 57600
type(solar_spectrum_type), public, save :: Solar_spect


!-------------------------------------------------------------------
!----- private data --------

logical, save :: module_is_initialized = .false.  ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                      contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PUBLIC SUBROUTINES
!            
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
!
! <SUBROUTINE NAME="esfsw_parameters_init">
!  <OVERVIEW>
!   Subroutine that initializes and set up shortwave radiation.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that initializes shortwave radiation calculation parameters such as
!   solar flux input at top of the atmosphere, number of shortwave bands
!   depending on the spectral resolution used, number of frequency points
!   in the gaussian quadrature algorithm, the number of streams used in
!   multiple stream flux algorithm, and the number of water vapor bands.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_parameters_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine esfsw_parameters_init

!------------------------------------------------------------------
!    esfsw_parameters_init is the constructor for esfsw_parameters_mod.
!------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!--------------------------------------------------------------------
!    process the namelist entries to obtain the parameters specifying
!    the solar spectral parameterization.
!--------------------------------------------------------------------
      if (trim(sw_resolution) == 'high') then
        Solar_spect%nbands = 25
        Solar_spect%nfrqpts = 72
        Solar_spect%nh2obands = 14
      else if (trim(sw_resolution) == 'low') then
        Solar_spect%nbands = 18
        Solar_spect%nfrqpts = 38
        Solar_spect%nh2obands = 9
      else
       ! call error_mesg ( 'esfsw_parameters_mod',   &
       !' sw_resolution must be specified as "high" or "low".', FATAL)
       write(6,*) "sw_resolution must be specified as high or low."
       stop
      endif
      if (sw_diff_streams == 4) then
        Solar_spect%nstreams = 4
      else if (sw_diff_streams == 1) then
        Solar_spect%nstreams = 1
      else
        !call error_mesg ( 'esfsw_parameters_mod',   &
        !  ' sw_diff_streams must be specified as either 1 or 4.', FATAL)
        write(6,*) "sw_diff_streams must be specified as either 1 or 4."
        stop
      endif

!---------------------------------------------------------------------
!    include the total number of wavenumbers in the solar parameter-
!    ization in the solar_spectrum_type variable.
!---------------------------------------------------------------------
      Solar_spect%tot_wvnums = TOT_WVNUMS

!-------------------------------------------------------------------
!    indicate that visible_band_indx has not yet been defined.
!-------------------------------------------------------------------
      Solar_spect%visible_band_indx = -10000000
      Solar_spect%visible_band_indx_iz = .false.

!-------------------------------------------------------------------
!    allocate space for the array components of the solar_spect_type
!    variable.
!-------------------------------------------------------------------
      allocate (Solar_spect%solflxband (Solar_spect%nbands) )
      allocate (Solar_spect%solflxbandref (Solar_spect%nbands) )
      allocate (Solar_spect%endwvnbands (0:Solar_spect%nbands) )
      allocate (Solar_spect%solarfluxtoa (Solar_spect%tot_wvnums))

!------------------------------------------------------------------
!    mark the module as initialized.
!------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------


end subroutine esfsw_parameters_init



!####################################################################
!
! <SUBROUTINE NAME="esfsw_parameters_end">
!  <OVERVIEW>
!   Subroutine that is the destructor for esfsw_parameters_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that deallocates module variables and marks the module  
!   as uninitialized.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_parameters_end
!  </TEMPLATE>
! </SUBROUTINE>
!

subroutine esfsw_parameters_end

!--------------------------------------------------------------------
!    esfsw_parameters_end is the destructor for esfsw_parameters_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('esfsw_parameters_mod',   &
        !     'module has not been initialized', FATAL )
        stop
      endif

!--------------------------------------------------------------------
!    deallocate the components of the solar_spect_type variable.
!---------------------------------------------------------------------
      deallocate (Solar_spect%solflxband, &
                  Solar_spect%solflxbandref, &
                  Solar_spect%endwvnbands, &
                  Solar_spect%solarfluxtoa)

!-------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

end subroutine esfsw_parameters_end



!###################################################################

      end module esfsw_parameters_mod
