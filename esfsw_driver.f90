
module esfsw_driver_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
!
! <OVERVIEW>
!  Code that initializes and calculates shortwave radiative quantities
!  such as flux and heating rate.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes the necessary shortwave radiative parameters
!  in the initialization subroutine. It then uses delta-eddington approximation
!  and doubling and adding technique to calculate solar flux and
!  heating rate.
! </DESCRIPTION>
!

!  shared radiation package modules:

use esfsw_parameters_mod, only:  Solar_spect, esfsw_parameters_init, &
                                 esfsw_parameters_end
use rad_utilities_mod,    only:  Rad_control, rad_utilities_init, &
                                 cldrad_properties_type, &
                                 cld_specification_type, &
                                 astronomy_type, &
                                 aerosol_diagnostics_type, &
                                 radiative_gases_type, &
                                 solar_spectrum_type,  &
                                 aerosol_type, aerosol_properties_type,&
                                 Cldrad_control, &
                                 atmos_input_type, surface_type, &
                                 sw_output_type, Sw_control
!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    esfsw_driver_mod is the internal driver for the esf shortwave
!    package.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: esfsw_driver.F90,v 18.0.2.1.2.1.2.1.2.1 2011/01/25 10:28:27 Richard.Hemler Exp $'
character(len=128)  :: tagname =  '$Name: testing $'

real, parameter :: PI               = 3.14159265358979323846
real, parameter :: GRAV             = 9.80
real, parameter :: RDGAS            = 287.04 
real, parameter :: KAPPA            = 2./7.  
real, parameter :: CP_AIR           = RDGAS/KAPPA
real, parameter :: SECONDS_PER_DAY  = 8.640000E+04
real, parameter :: RADCON_MKS       = (GRAV/CP_AIR)*SECONDS_PER_DAY
real, parameter :: O2MIXRAT         = 2.0953E-01
real, parameter :: RHOAIR           = 1.292269
real, parameter :: PSTD_MKS         = 101325.0
real, parameter :: WTMAIR           = 2.896440E+01

!---------------------------------------------------------------------
!-------  interfaces --------

public    &
         esfsw_driver_init, swresf,   &
         esfsw_driver_end

private     &
!   called from swresf:
         adding, deledd


!---------------------------------------------------------------------
!-------- namelist  ---------

 logical, save ::  do_ica_calcs = .false.         ! do independent column
                                                  ! calculations when sto-
                                                  ! chastic clouds are
                                                  ! active ?
logical, save ::  do_rayleigh_all_bands = .true.  ! rayleigh scattering 
                                                  ! calculated in all sw 
                                                  ! bands ?
logical, save ::  do_herzberg = .false.           ! include the herzberg 
                                                  ! effect on the o2 
                                                  ! optical depth ?
logical, save ::  do_quench = .false.             ! include the quenching
                                                  ! effect of non-LTE 
                                                  ! processes on the co2 
                                                  ! optical depth ?
logical, save ::  do_ch4_sw_effects = .true.      ! the shortwave effects
                                                  ! of ch4 are included ?
logical, save ::  do_n2o_sw_effects = .true.      ! the shortwave effects
                                                  ! of n2o are included ?
!logical, save ::  do_coupled_stratozone = .false. ! include the coupled
!                                                  ! stratospheric ozone effects?
!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!    variables associated with absorptivity and sw transmission for
!    various gaseous atmospheric components
!
! powph2o      = the scaling factor used in the fit of the h2o          
!                transmission function                                  
!                                                                       
! p0h2o        = the reference pressure (mb) used in the fit of the     
!                h2o transmission function                              
!                                                                       
! c(n)co2(str) = coefficients for the absorptivity expression for co2   
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! c(n)o2(str)  = coefficients for the absorptivity expression for o2    
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! ""(schrun)   = coefficients for the absorptivity expression for the   
!                Schuman-Runge o2 band (non-scaled only)                
!                                                                       
! kh2o         =  the psuedo-absorption coefficients in cm2/gm for h2o  
!                                                                       
! ko3          = the absorption coefficients in cm2/gm for o3           
!                                                                       
! wtfreq       = the weight associated with each exponential term       
! strterm      = logical flag to indicate whether or not a h2o pseudo-  
!                absorption coefficient is assigned a non-scaled        
!                (true) or pressure-scaled (false) gas amount           
!---------------------------------------------------------------------
real, dimension (:), allocatable, save :: c1co2, c1co2str, c1o2, c1o2str, &
                                          c2co2, c2co2str, c2o2, c2o2str, &
                                          c3co2, c3co2str, c3o2, c3o2str, &
                                          c4co2, c4co2str, c4o2, c4o2str, &
                                          c1ch4, c1ch4str, c2ch4,         &
                                          c2ch4str, c3ch4, c3ch4str,      &
                                          c4ch4, c4ch4str,                &
                                          c1n2o, c1n2ostr, c2n2o,         &
                                          c2n2ostr, c3n2o, c3n2ostr,      &
                                          c4n2o, c4n2ostr,                &
                                          powph2o, p0h2o
real, save                             :: c1o2strschrun, c2o2strschrun,   &
                                          c3o2strschrun, c4o2strschrun
real, dimension (:), allocatable, save    :: kh2o, ko3, wtfreq
logical, dimension(:), allocatable, save  :: strterm

!---------------------------------------------------------------------
!    quantities associated with solar spectral parameterization
!                                                                       
! firstrayband = the first band number where the contribution by        
!                rayleigh scattering is included in the solar           
!                calculations                                           
!                                                                       
! nirbands     = the number of bands in the near-infrared (used in      
!                assigning the value of the surface albedo for the      
!                near-infrared, and the visible and ultraviolet         
!                regions, separately)                                   
! nfreqpts     = the number of pseudo-monochromatic frequencies         
! solflxband   = the solar flux in each parameterization band           
! solflxbandref = the solar flux in each parameterization band, used for
!                 defining band-averaged optical parameters. If the
!                 solar constant is time-invariant, it is also the solar
!                 flux in each parameterization band (solflxband).
! vis_wvnum    = the wavenumber of visible light (corresponds to
!                wavelength of 0.55 microns) [ cm **(-1) ]
!---------------------------------------------------------------------
real, save                                :: refquanray, solflxtotal
integer, save                             :: firstrayband, nirbands
integer, dimension (:), allocatable, save :: nfreqpts
real,    dimension(:), allocatable, save  :: solflxband
real,    dimension(:), allocatable, save  :: solflxbandref
real, dimension(:), allocatable, save     :: wtstr, cosangstr
real, dimension(4), parameter             :: wtstr_4 =             &
                                       (/0.347854845, 0.652145155, &
                                         0.347854845, 0.652145155/)

integer, save   :: nbands, tot_wvnums, nfrqpts, nh2obands, nstreams
logical, save   :: nstr4 = .false.
real, parameter :: vis_wvnum = 1.0E+04/0.55
real, parameter :: wvnum_870 = 1.0E+04/0.87
real, parameter :: wvnum_340 = 1.0E+04/0.34
real, parameter :: wvnum_380 = 1.0E+04/0.38
real, parameter :: wvnum_440 = 1.0E+04/0.44
real, parameter :: wvnum_670 = 1.0E+04/0.67
real, parameter :: one_micron_wvnum = 1.0E+04/1.00
real, parameter :: onepsix_micron_wvnum = 1.0E+04/1.61
integer, save   :: onepsix_band_indx

!---------------------------------------------------------------------
!    variables associated with rayleigh scattering
!---------------------------------------------------------------------
real, dimension (:), allocatable, save    :: betaddensitymol

!----------------------------------------------------------------------
!    variables associated with total optical path of species ? - smf
!----------------------------------------------------------------------
real, save                            :: toto2strmaxschrun
real, dimension(:), allocatable, save :: totco2max, totco2strmax, &
                                         toto2max, toto2strmax,   &
                                         totch4max, totch4strmax, &
                                         totn2omax, totn2ostrmax

!----------------------------------------------------------------------
!    variables associated with the herzberg effect. wtmo2 is the mol-
!    ecular weight of o2. herzberg_fac is a factor used in the last 
!    shortwave band, so that herzberg_fac*wo2 yields a correction for 
!    the o2 optical depth to account for the herzberg o2 heating. this 
!    is done only when do_herzberg is true.
!----------------------------------------------------------------------
real, parameter   :: wtmo2        = 3.19988E+01  
real, parameter   :: herzberg_fac = 9.9488377E-3

!----------------------------------------------------------------------
!    variables associated with the quenching effect. co2_quenchfac is a
!    multiplication factor that reduces the co2 gas optical depth, and 
!    hence, the solar heating in the upper atmosphere, to account for 
!    "quenching" due to non-LTE processes. co2_quenchfac_height is the 
!    reference height for co2_quenchfac [ meters ].
!----------------------------------------------------------------------
real, dimension(30), save :: co2_quenchfac
data co2_quenchfac /1.0,.954,.909,.853,.800,.747,.693,.637,.583, .526,&
                    .467,.416,.368,.325,.285,.253,.229,.206,.186,.170,&
                    .163,.156,.151,.144,.138,.132,.127,.124,.068,.037/

real, dimension(30), save :: co2_quenchfac_height
data co2_quenchfac_height /67304.,68310.,69303.,70288.,71267.,72245.,&
                           73221.,74195.,75169.,76141.,77112.,78082.,&
                           79051.,80018.,80985.,81950.,82914.,83876.,&
                           84838.,85798.,86757.,87715.,88672.,89627.,&
                           90582.,91535.,92487.,93438.,94387.,106747./



!---------------------------------------------------------------------
!    miscellaneous variables
!---------------------------------------------------------------------
integer, parameter :: NSOLWG = 1
real, dimension(NSOLWG), save :: gausswt
logical, save      :: module_is_initialized = .false.
logical, save      :: do_esfsw_band_diagnostics = .false.

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! <SUBROUTINE NAME="esfsw_driver_init">
!  <OVERVIEW>
!   Subroutine that defines the time-independent quantities associated
!   with the incoming shortwave radiation in the multiple-band solar
!   radiation parameterization.
!  </OVERVIEW>
!  <DESCRIPTION>
!   It first reads in the input namelist and then allocates gas absorption
!   coefficient variables. It then reads in the shortwave input namelist
!   file and assigns the gas absorption coefficients. Rayleigh scattering
!   coefficient is also calculated based on the temperature and pressure
!   structure of the atmosphere.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_driver_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine esfsw_driver_init
 
!---------------------------------------------------------------------- 
!    esfsw_driver_init is the constructor for esfsw_driver_mod. it
!    defines the time-independent quantities associated with the 
!    incoming shortwave radiation in the multiple-band solar radiation 
!    parameterization.                                    
!---------------------------------------------------------------------- 

!---------------------------------------------------------------------
!  local variables:

      real,    dimension(Solar_spect%nbands)    :: freqnu
      real,    dimension(Solar_spect%nstreams)  :: ptstr 
      integer, dimension(0:Solar_spect%nbands) :: endwvnbands

      integer, dimension(:), allocatable, save  :: nwvnsolar
      real   , dimension(:), allocatable, save  :: solint   

      real, dimension(4), parameter :: ptstr_4 = (/-0.861136312,  &
                                                   -0.339981044,  &
                                                    0.861136312,  &
                                                    0.339981044 /)
      real, parameter :: ptstr_1     = 0.2
      real, parameter :: temprefray  = 288.15
      real, parameter :: pressrefray = 101325.       ! MKS units
      real, parameter :: densmolref  = 2.54743E+19
      real, parameter :: convfac     = 1.0E+18
      real            :: corrfac, gamma, f1, f2, f3, pinteg, &
                         twopiesq, densmolrefsqt3, wavelength,  &
                         freqsq, ristdm1, ri
      integer         :: nband, ni, nw, nw1, nw2, nintsolar
      integer         :: i
      integer         :: n
      real, parameter :: input_flag = 1.0e-99

!---------------------------------------------------------------------
!  local variables:
!                                                                       
!      freqnu   
!      ptstr          gaussian points and weights for evaluation of the
!                     diffuse beam.
!      nwvnsolar      the number of wavenumbers in each region where the
!                     solar flux is constant                         
!      solint         the solar flux in watts per meter**2 in each      
!                     wavenumber region where it is constant    
!      endwvnbands    the wavenumber value for the band limits   
!      file_name
!      ptstr_4    
!      ptstr_1     
!      temprefray     reference temperature used in defining rayleigh
!                     optical depth
!      pressrefray    reference pressure used in defining rayleigh
!                     optical depth [ Pa ]
!      densmolref     reference density used in defining rayleigh
!                     optical depth
!      convfac     
!      corrfac
!      gamma
!      f1
!      f2
!      f3
!      pinteg
!      twopiesq
!      densmolrefsqt3
!      wavelength
!      freqsq
!      ristdm1
!      ri
!      iounit
!      nband
!      nf
!      ni
!      nw
!      nw1
!      nw2
!      nintsolar      the number of wavenumber regions where the  
!                     solar flux is constant.   
!      unit
!      io
!      ierr
!      i
!                                                                       
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call rad_utilities_init
      call esfsw_parameters_init

!---------------------------------------------------------------------
!    define flag indicating if ICA calculations being done.
!---------------------------------------------------------------------
      Cldrad_control%do_ica_calcs = do_ica_calcs
      Cldrad_control%do_ica_calcs_iz = .true.

!---------------------------------------------------------------------
!    allocate module variables
!---------------------------------------------------------------------
      nbands = Solar_spect%nbands       
      tot_wvnums = Solar_spect%tot_wvnums
      nfrqpts = Solar_spect%nfrqpts
      nstreams = Solar_spect%nstreams
      nh2obands = Solar_spect%nh2obands
      allocate ( betaddensitymol (nbands) )
      allocate ( c1co2   (nh2obands),  &
                 c1co2str(nh2obands),  &
                 c1o2    (nh2obands),  &
                 c1o2str (nh2obands),  &
                 c2co2   (nh2obands),  &
                 c2co2str(nh2obands),  &
                 c2o2    (nh2obands),  &
                 c2o2str (nh2obands),  &
                 c3co2   (nh2obands),  &
                 c3co2str(nh2obands),  &
                 c3o2    (nh2obands),  &
                 c3o2str (nh2obands),  &
                 c4co2   (nh2obands),  &
                 c4co2str(nh2obands),  &
                 c4o2    (nh2obands),  &
                 c4o2str (nh2obands),  &
                 c1ch4   (nh2obands),  &
                 c1ch4str(nh2obands),  &
                 c2ch4   (nh2obands),  &
                 c2ch4str(nh2obands),  &
                 c3ch4   (nh2obands),  &
                 c3ch4str(nh2obands),  &
                 c4ch4   (nh2obands),  &
                 c4ch4str(nh2obands),  &
                 c1n2o   (nh2obands),  &
                 c1n2ostr(nh2obands),  &
                 c2n2o   (nh2obands),  &
                 c2n2ostr(nh2obands),  &
                 c3n2o   (nh2obands),  &
                 c3n2ostr(nh2obands),  &
                 c4n2o   (nh2obands),  &
                 c4n2ostr(nh2obands),  &
                 powph2o (nh2obands),  &
                 p0h2o   (nh2obands)    )
      allocate ( nfreqpts        (nbands) )
      allocate ( solflxband      (nbands) )
      allocate ( solflxbandref   (nbands) )
      allocate ( kh2o            (nfrqpts),   & 
                 ko3             (nfrqpts),   &
                 wtfreq          (nfrqpts),   &
                 strterm         (nfrqpts)   )
      allocate ( wtstr           (nstreams),   & 
                 cosangstr       (nstreams)  )
      allocate ( totco2max    (nh2obands),     &
                 totco2strmax (nh2obands),     &
                 toto2max     (nh2obands),     &
                 toto2strmax  (nh2obands),   &
                 totch4max    (nh2obands),   &
                 totch4strmax (nh2obands),   &
                 totn2omax    (nh2obands),   &
                 totn2ostrmax (nh2obands)    )

      betaddensitymol = 0.0 ; c1co2    = 0.0 ; c1co2str = 0.0
      c1o2     = 0.0 ; c1o2str  = 0.0 ; c2co2    = 0.0
      c2co2str = 0.0 ; c2o2     = 0.0 ; c2o2str  = 0.0
      c3co2    = 0.0 ; c3co2str = 0.0 ; c3o2     = 0.0
      c3o2str  = 0.0 ; c4co2    = 0.0 ; c4co2str = 0.0
      c4o2     = 0.0 ; c4o2str  = 0.0 ; powph2o  = 0.0
      p0h2o    = 0.0 ; nfreqpts        = 0 ; solflxband      = 0.0
      solflxbandref   = 0.0 ; kh2o            = 0.0 ; ko3             = 0.0
      wtfreq          = 0.0 ; strterm     = .FALSE. ; wtstr           = 0.0
      cosangstr       = 0.0 ; totco2max     = 0.0 ; totco2strmax  = 0.0
      toto2max      = 0.0 ; toto2strmax   = 0.0
      c1ch4    = 0.0 ; c1ch4str = 0.0; c2ch4    = 0.0
      c2ch4str = 0.0 ; c3ch4    = 0.0 ; c3ch4str = 0.0
      c4ch4    = 0.0 ; c4ch4str = 0.0
      totch4max     = 0.0 ; totch4strmax  = 0.0
      c1n2o    = 0.0 ; c1n2ostr = 0.0; c2n2o    = 0.0
      c2n2ostr = 0.0 ; c3n2o    = 0.0 ; c3n2ostr = 0.0
      c4n2o    = 0.0 ; c4n2ostr = 0.0
      totn2omax     = 0.0 ; totn2ostrmax  = 0.0

!---------------------------------------------------------------------
!    allocate local variables.
!---------------------------------------------------------------------
      if (nstreams == 4) then
        ptstr(:) = ptstr_4(:)
        wtstr(:) = wtstr_4(:)
        nstr4 = .true.
      else if (nstreams == 1) then
        ptstr(1) = ptstr_1
        wtstr(1) = 1.0
        nstr4 = .false.
      endif

!---------------------------------------------------------------------
!    read input file for band positions, solar fluxes and band
!    strengths.
!---------------------------------------------------------------------
      if (nbands == 25 .and. nfrqpts == 72) then 
        !file_name = 'INPUT/esf_sw_input_data_n72b25'
        solflxbandref(1:nbands)=(/ 12.1587,  6.5070,   10.7300,  23.8226,  19.2689, 43.7116, 35.7886, 135.0955, 239.2806, &
                                  222.9263, 138.7890, 182.3105, 101.2186, 72.2298, 48.5104, 28.2587, 15.4827,  6.0424,   &
                                  3.7148 ,  3.0384,   1.7734,   1.9695,   3.1789,  1.0869,  1.0672/)
        nfreqpts(1:nbands)=(/ 6, 1, 5, 6, 1, 9, 5, 9, 7, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
        endwvnbands(1:nbands)=(/ 2500,  2900,  3400,  4200,  4700,  5600,  6200,  8200,  11500, 14600, 16700, 20000, 22300, 24600, &
                                27500, 30000, 31900, 33000, 33800, 34500, 35300, 36500, 40000, 43300, 57600 /)
        FIRSTRAYBAND=9
        NIRBANDS=10
        powph2o(1:NH2OBANDS)=(/ 0.84, 0.00, 0.96, 0.60, 0.00, 0.88, 0.91, 0.80, 0.49, 0.38, 0.00, 0.00, 0.00, 0.00 /)
        p0h2o(1:NH2OBANDS)=(/  3.00, 1013.25,  500.00,   50.00, 1013.25,   40.00,  700.00, 1013.25, 1013.25, &
                             700.00, 1013.25, 1013.25, 1013.25, 1013.25 /)
        c1co2(1:NH2OBANDS)=(/ 5.4E+02, 1.0E-99, 5.0E-02, 1.6E-01, 1.0E-99, 1.3E-02, 2.4E-06, 3.3E-04, 2.8E+00, 1.0E-99, 1.0E-99, &
                              1.0E-99, 1.0E-99, 1.0E-99 /)
        c1co2str(1:NH2OBANDS)=(/ 7.2E-03, 1.0E-99, 1.3E-04, 7.1E-03, 1.0E-99, 3.2E-02, 6.6E-01, 2.8E-04, 7.0E-02, 1.0E-99, &
                                 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2co2(1:NH2OBANDS)=(/ 7.2E-04, 1.0E-99, 1.0E+03, 1.0E-02, 1.0E-99, 4.2E-01, 9.9E+03, 6.6E+00, 5.8E+05, 1.0E-99, 1.0E-99, &
                              1.0E-99, 1.0E-99, 1.0E-99/)
        c2co2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E+02, 4.4E-03, 1.0E-99, 6.3E-01, 7.2E+02, 9.5E+00, 1.0E+03, 1.0E-99, &
                                 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3co2(1:NH2OBANDS)=(/ 1.8E-05, 1.0E-99, 6.8E-02, 1.0E-01, 1.0E-99, 4.1E-01, 9.9E-01, 5.6E-01, 2.8E-02, 1.0E-99, 1.0E-99, &
                              1.0E-99, 1.0E-99, 1.0E-99 /)
        c3co2str(1:NH2OBANDS)=(/ 4.0E-01, 1.0E-99, 5.0E-01, 9.7E-02, 1.0E-99, 3.3E-02, 4.9E-03, 4.6E-01, 9.4E-03, 1.0E-99, &
                                 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 3.4E-02, 1.0E-99, 3.9E-03, 1.4E-07, &
                             1.0E-99, 1.0E-99, 1.0E-99 /)
        c1o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 7.9E-06, 1.0E-99, 1.9E-03,  &
                                4.9E-06, 1.0E-99, 1.0E-99, 1.0E-99/)
        c1o2strschrun=2.7E-02
        c2o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 4.5E+05, 1.0E-99, 3.3E+03, 5.2E+08, &
                             1.0E-99, 1.0E-99, 1.0E-99 /)
        c2o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 9.8E+03, 1.0E-99, 1.5E+02, 7.0E+04, &
                                1.0E-99, 1.0E-99, 1.0E-99 /)
        c2o2strschrun=1.2E-01
        c3o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 8.4E-02, 1.0E-99, 2.1E-01, 7.8E-01, &
                             1.0E-99, 1.0E-99, 1.0E-99 /)
        c3o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 4.7E-01, 1.0E-99, 9.7E-02, 4.6E-01, &
                                1.0E-99, 1.0E-99, 1.0E-99 /)
        c3o2strschrun=1.5E-01
        !c1ch4(1:NH2OBANDS)=(/ /)
        !c1ch4str(1:NH2OBANDS)=(/ /)
        !c2ch4(1:NH2OBANDS)=(/ /)
        !c2ch4str(1:NH2OBANDS)=(/ /)
        !c3ch4(1:NH2OBANDS)=(/ /)
        !c3ch4str(1:NH2OBANDS)=(/ /)
        !c1n2o(1:NH2OBANDS)=(/ /)
        !c1n2ostr(1:NH2OBANDS)=(/ /)
        !c2n2o(1:NH2OBANDS)=(/ /)
        !c2n2ostr(1:NH2OBANDS)=(/ /)
        !c3n2o(1:NH2OBANDS)=(/ /)
        !c3n2ostr(1:NH2OBANDS)=(/ /)
        ! MJT fix
        c1ch4(1:NH2OBANDS)=1.e-99
        c1ch4str(1:NH2OBANDS)=1.e-99
        c2ch4(1:NH2OBANDS)=1.e-99
        c2ch4str(1:NH2OBANDS)=1.e-99
        c3ch4(1:NH2OBANDS)=1.e-99
        c3ch4str(1:NH2OBANDS)=1.e-99
        c1n2o(1:NH2OBANDS)=1.e-99
        c1n2ostr(1:NH2OBANDS)=1.e-99
        c2n2o(1:NH2OBANDS)=1.e-99
        c2n2ostr(1:NH2OBANDS)=1.e-99
        c3n2o(1:NH2OBANDS)=1.e-99
        c3n2ostr(1:NH2OBANDS)=1.e-99
        wtfreq(1:nfrqpts)=(/ 9.312504E-02, 1.626168E-01, 1.433068E-01, 3.071999E-01, 2.737514E-01, 2.000000E-02, 1.000000E+00, &
                             7.874545E-02, 2.476534E-01, 3.794250E-01, 2.755277E-01, 1.864845E-02, 1.921193E-01, 2.751345E-01, &
                             2.480785E-01, 2.262776E-01, 1.839013E-02, 4.000000E-02, 1.000000E+00, 3.837498E-02, 8.791451E-02, &
                             1.637795E-01, 1.189705E-01, 2.020078E-01, 1.509591E-01, 1.693291E-01, 6.466451E-02, 4.000000E-03, &
                             1.605849E-02, 8.165883E-02, 6.004498E-02, 3.943264E-01, 4.479113E-01, 4.476276E-02, 1.119772E-01, &
                             1.180453E-01, 6.352831E-02, 8.756711E-02, 5.027652E-02, 1.652990E-01, 3.485438E-01, 1.000000E-02, &
                             5.315826E-03, 2.793918E-02, 4.138399E-02, 1.803244E-01, 1.955244E-01, 1.832192E-01, 3.662930E-01, &
                             3.601121E-03, 1.291813E-02, 4.997123E-02, 8.941965E-02, 1.016408E-01, 1.407492E-01, 2.304239E-01, &
                             3.712796E-01, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, &
                             1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, &
                             1.000000E+00, 1.000000E+00 /)
        kh2o(1:nfrqpts)=(/ 4.000000E+00, 2.546100E-01, 1.489350E-02, 5.228619E-04, 0.000000E+00, 1.000000E+03, 3.524000E-02, &
                           4.000000E+01, 2.923541E+00, 3.820742E-01, 4.629296E-02, 0.000000E+00, 5.000000E+01, 4.756542E+00, &
                           3.740251E-01, 2.046388E-02, 0.000000E+00, 1.000000E+03, 1.443000E-02, 1.000000E+02, 6.131070E+00, &
                           8.325077E-01, 2.336348E-01, 5.893526E-02, 4.474983E-03, 3.560266E-04, 0.000000E+00, 4.000000E+03, &
                           4.000000E+00, 2.855524E-01, 5.188692E-02, 6.285804E-03, 0.000000E+00, 1.000000E+02, 1.124617E+01, &
                           2.061909E+00, 6.405800E-01, 2.300286E-01, 8.792879E-02, 1.722225E-02, 0.000000E+00, 5.000000E+02, &
                           5.000000E+01, 7.355553E+00, 2.057103E+00, 4.088761E-01, 6.036039E-02, 1.437820E-02, 0.000000E+00, &
                           3.000000E+00, 7.057908E-01, 1.566582E-01, 3.665539E-02, 1.463247E-02, 5.343198E-03, 1.353690E-03, &
                           0.000000E+00, 1.925000E-03, 2.025000E-03, 3.664000E-05, 6.239000E-05, 0.000000E+00, 0.000000E+00, &
                           0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                           0.000000E+00, 0.000000E+00 /)
        ko3(1:nfrqpts)=(/ 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          5.041992E+00, 5.041992E+00, 5.041992E+00, 5.041992E+00, 5.041992E+00, 5.041992E+00, 5.041992E+00, &
                          5.041992E+00, 3.745368E+01, 4.021207E+01, 7.062558E+00, 9.620166E-01, 0.000000E+00, 9.738058E+00, &
                          2.378498E+02, 1.567498E+03, 5.394976E+03, 1.207743E+04, 2.706951E+04, 5.277161E+04, 1.177364E+05, &
                          1.035882E+05, 2.475921E+04 /)
        strterm(1:nfrqpts)=(/ .false., .false., .false., .false., .false., .true.,  .false., .false., .false., .false., .false., &
                              .false., .false., .false., .false., .false., .false., .true.,  .false., .false., .false., .false., &
                              .false., .false., .false., .false., .false., .true.,  .false., .false., .false., .false., .false., &
                              .false., .false., .false., .false., .false., .false., .false., .false., .true.,  .false., .false., &
                              .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
                              .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
                              .false., .false., .false., .false., .false., .false. /)
        nintsolar=151
        allocate ( nwvnsolar (nintsolar) )
        allocate ( solint    (nintsolar) )
        nwvnsolar(1:nintsolar)=(/  100,   11,   14,   18,   24,   33,   50,   83,   12,   12,   13,   15,   15,   17,   18,   20, &
                                    21,   24,   26,   30,   32,   37,   42,   47,   55,   64,   76,   91,  111,  139,  179,  238, &
                                   333,   41,   42,   45,   46,   48,   51,   53,   55,   58,   61,   64,   68,   71,   75,   79, &
                                    84,   89,   95,  101,  107,  115,  123,  133,  142,  154,  167,  181,  197,  217,  238,  263, &
                                   293,  326,  368,  417,  476,  549,  641,  758,  909,  101,  103,  105,  108,  109,  112,  115, &
                                   117,  119,  122,  125,  128,  130,  134,  137,  140,  143,  147,  151,  154,  158,  163,  166, &
                                   171,  175,  181,  185,  190,  196,  201,  207,  213,  219,  227,  233,  240,  248,  256,  264, &
                                   274,  282,  292,  303,  313,  325,  337,  349,  363,  377,  392,  408,  425,  444,  462,  483, &
                                   505,  529,  554,  580,  610,  641,  675,  711,  751,  793,  841,  891,  947, 1008, 1075, 1150, &
                                  1231, 1323, 1425, 1538, 1667, 1633,14300 /)
         solint(1:nintsolar)=(/ 1.60000E-06, 2.88000E-05, 3.60000E-05, 4.59200E-05, 6.13200E-05, 8.55000E-05, 1.28600E-04, &
                                2.16000E-04, 2.90580E-04, 3.10184E-04, 3.34152E-04, 3.58722E-04, 3.88050E-04, 4.20000E-04, &
                                4.57056E-04, 4.96892E-04, 5.45160E-04, 6.00600E-04, 6.53600E-04, 7.25040E-04, 7.98660E-04, &
                                9.11200E-04, 1.03680E-03, 1.18440E-03, 1.36682E-03, 1.57560E-03, 1.87440E-03, 2.25500E-03, &
                                2.74500E-03, 3.39840E-03, 4.34000E-03, 5.75400E-03, 7.74000E-03, 9.53050E-03, 9.90192E-03, &
                                1.02874E-02, 1.06803E-02, 1.11366E-02, 1.15830E-02, 1.21088E-02, 1.26420E-02, 1.32250E-02, &
                                1.38088E-02, 1.44612E-02, 1.51164E-02, 1.58878E-02, 1.66500E-02, 1.75140E-02, 1.84450E-02, &
                                1.94106E-02, 2.04864E-02, 2.17248E-02, 2.30640E-02, 2.44470E-02, 2.59840E-02, 2.75940E-02, &
                                2.94138E-02, 3.13950E-02, 3.34800E-02, 3.57696E-02, 3.84054E-02, 4.13490E-02, 4.46880E-02, &
                                4.82220E-02, 5.22918E-02, 5.70078E-02, 6.19888E-02, 6.54720E-02, 6.69060E-02, 6.81226E-02, &
                                6.97788E-02, 7.12668E-02, 7.27100E-02, 7.31610E-02, 7.33471E-02, 7.34814E-02, 7.34717E-02, &
                                7.35072E-02, 7.34939E-02, 7.35202E-02, 7.33249E-02, 7.31713E-02, 7.35462E-02, 7.36920E-02, &
                                7.23677E-02, 7.25023E-02, 7.24258E-02, 7.20766E-02, 7.18284E-02, 7.32757E-02, 7.31645E-02, &
                                7.33277E-02, 7.36128E-02, 7.33752E-02, 7.28965E-02, 7.24924E-02, 7.23307E-02, 7.21050E-02, &
                                7.12620E-02, 7.10903E-02, 7.12714E-02, 7.08012E-02, 7.03752E-02, 7.00350E-02, 6.98639E-02, &
                                6.90690E-02, 6.87621E-02, 6.52080E-02, 6.65184E-02, 6.60038E-02, 6.47615E-02, 6.44831E-02, &
                                6.37206E-02, 6.24102E-02, 6.18698E-02, 6.06320E-02, 5.83498E-02, 5.67028E-02, 5.51232E-02, &
                                5.48645E-02, 5.12340E-02, 4.85581E-02, 4.85010E-02, 4.79220E-02, 4.44058E-02, 4.48718E-02, &
                                4.29373E-02, 4.15242E-02, 3.81744E-02, 3.16342E-02, 2.99615E-02, 2.92740E-02, 2.67484E-02, &
                                1.76904E-02, 1.40049E-02, 1.46224E-02, 1.39993E-02, 1.19574E-02, 1.06386E-02, 1.00980E-02, &
                                8.63808E-03, 6.52736E-03, 4.99410E-03, 4.39350E-03, 2.21676E-03, 1.33812E-03, 1.12320E-03, &
                                5.59000E-04, 3.60000E-04, 2.98080E-04, 7.46294E-05 /)
      else if (nbands == 18 .and. nfrqpts == 38) then
        !file_name = 'INPUT/esf_sw_input_data_n38b18'
        solflxbandref(1:nbands)=(/ 11.8139,  40.2847, 233.3027, 241.2290, 221.9901, 138.1064, 181.7117, 105.2140,  73.1235, &
                                  51.6579,  48.8855,   4.6404,   5.5280,   1.9658,   1.6787,   3.5985,   0.9767,   1.1128 /)
        nfreqpts(1:nbands)=(/ 1, 8, 9, 5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        endwvnbands(1:nbands)=(/  2500,  4200,  8200, 11500, 14600, 16700, 20000, 22300, 24600, 27500, 32400, 33300, 34500, &
                                35300, 36500, 40000, 43300, 57600 /)
        FIRSTRAYBAND=4
        NIRBANDS=5
        powph2o(1:NH2OBANDS)=(/ 0.00, 0.80, 0.73, 0.90, 0.00, 0.00, 0.00, 0.00, 0.00 /)
        p0h2o(1:NH2OBANDS)=(/ 1013.25,  200.00,  800.00,  600.00, 1013.25, 1013.25, 1013.25, 1013.25, 1013.25 /)
        c1co2(1:NH2OBANDS)=(/ 6.2E+01, 2.8E+01, 4.6E-03, 2.0E+00, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1co2str(1:NH2OBANDS)=(/ 7.0E-03, 2.8E-03, 3.1E-04, 7.5E-01, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2co2(1:NH2OBANDS)=(/ 7.6E-04, 3.2E-02, 3.0E+00, 9.6E+05, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2co2str(1:NH2OBANDS)=(/ 2.4E-06, 5.6E-04, 1.3E-01, 1.0E+03, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3co2(1:NH2OBANDS)=(/ 1.6E-04, 3.1E-04, 3.4E-01, 4.9E-02, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3co2str(1:NH2OBANDS)=(/ 3.9E-01, 4.6E-01, 5.2E-01, 9.3E-04, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 9.8E-05, 5.9E-05, 1.3E-03, 9.8E+00, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 9.5E-02, 2.5E-02, 9.7E-02, 2.3E-06, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1o2strschrun=2.7E-02
        c2o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 4.1E+04, 4.8E+03, 1.4E+03, 9.5E+05, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 1.5E+04, 7.0E+02, 1.9E+02, 2.2E+04, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2o2strschrun=1.2E-01
        c3o2(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 3.4E-01, 2.4E-01, 2.8E-01, 3.4E-04, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3o2str(1:NH2OBANDS)=(/ 1.0E-99, 1.0E-99, 4.2E-03, 4.1E-04, 4.2E-03, 5.1E-01, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3o2strschrun=1.5E-01
        c1ch4(1:NH2OBANDS)=(/ 2.1E-02, 3.4E-02, 3.9E+00, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1ch4str(1:NH2OBANDS)=(/ 3.5E-03, 9.7E-03, 3.7E-03, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2ch4(1:NH2OBANDS)=(/ 9.6E-03, 4.5E-03, 2.8E+00, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2ch4str(1:NH2OBANDS)=(/ 1.2E-03, 2.0E-03, 2.0E-01, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3ch4(1:NH2OBANDS)=(/ 3.1E-01, 4.5E-01, 9.6E-04, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3ch4str(1:NH2OBANDS)=(/ 4.8E-01, 4.0E-01, 7.1E-01, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1n2o(1:NH2OBANDS)=(/ 9.7E-01, 7.2E-01, 6.6E+02, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c1n2ostr(1:NH2OBANDS)=(/ 4.5E-02, 2.0E-02, 5.8E-04, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2n2o(1:NH2OBANDS)=(/ 3.9E-03, 7.8E-01, 9.9E+04, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c2n2ostr(1:NH2OBANDS)=(/ 2.5E-05, 5.2E-03, 3.7E-01, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3n2o(1:NH2OBANDS)=(/ 1.6E-02, 2.3E-02, 8.5E-03, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        c3n2ostr(1:NH2OBANDS)=(/ 6.0E-01, 7.9E-01, 7.1E-01, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /)
        wtfreq(1:nfrqpts)=(/ 1.000000E+00, 1.371226E-02, 4.076332E-02, 1.079868E-01, 1.985192E-01, 2.061013E-01, 2.034842E-01, &
                             1.485116E-01, 8.092132E-02, 1.055215E-03, 9.130769E-04, 1.033041E-02, 6.189962E-02, 1.476831E-01, &
                             1.448796E-01, 9.362167E-02, 2.026028E-01, 3.370145E-01, 2.106552E-02, 1.087997E-01, 1.191990E-01, &
                             3.284672E-01, 4.224686E-01, 9.700000E-01, 3.000000E-02, 1.000000E+00, 1.000000E+00, 1.000000E+00, &
                             1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00, &
                             1.000000E+00, 1.000000E+00, 1.000000E+00 /)
        kh2o(1:nfrqpts)=(/ 0.000000E+00, 8.000000E+03, 5.405314E+02, 6.798037E+01, 1.055235E+01, 8.784034E-01, 1.116962E-01, &
                           1.967001E-02, 0.000000E+00, 1.131744E+04, 6.019127E+03, 7.207771E+02, 5.483035E+01, 6.228471E+00, &
                           6.829058E-01, 8.266925E-02, 1.699428E-02, 0.000000E+00, 2.764583E+01, 1.154150E+00, 2.736499E-01, &
                           4.159987E-02, 0.000000E+00, 1.040000E-02, 4.800000E-01, 1.944250E-03, 2.065500E-03, 3.627360E-05, &
                           6.176610E-05, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                           0.000000E+00, 0.000000E+00, 0.000000E+00 /)
        ko3(1:nfrqpts)=(/ 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
                          0.000000E+00, 0.000000E+00, 5.041992E+00, 5.041992E+00, 3.745368E+01, 4.021207E+01, 7.062558E+00, &
                          9.620166E-01, 0.000000E+00, 1.470320E+02, 2.723337E+03, 9.135737E+03, 2.706951E+04, 5.277161E+04, &
                          1.177364E+05, 1.035882E+05, 2.475921E+04 /)
        strterm(1:nfrqpts)=(/ .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                              .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                              .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                              .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
        nintsolar=6064
        allocate ( nwvnsolar (nintsolar) )
        allocate ( solint    (nintsolar) )
        nwvnsolar(1:1400)=(/                                                               &
                     100,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 /)
        nwvnsolar(1401:2800)=(/                                                            &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 /)
        nwvnsolar(2801:4200)=(/                                                            &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1,   1, &
                       1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   1,   1, &
                       1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1, &
                       1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   1, &
                       1,   2,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1, &
                       1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   1,   2,   1,   1, &
                       1,   1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   2,   1,   1, &
                       1,   1,   1,   2,   1,   1,   1,   1,   1,   1,   2,   1,   1,   1, &
                       1,   1,   2,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   1, &
                       2,   1,   1,   1,   1,   1,   2,   1,   1,   1,   1,   2,   1,   1, &
                       1,   1,   1,   2,   1,   1,   1,   1,   2,   1,   1,   1,   1,   2, &
                       1,   1,   1,   1,   2,   1,   1,   1,   1,   2,   1,   1,   1,   1, &
                       2,   1,   1,   1,   1,   2,   1,   1,   1,   1,   2,   1,   1,   1, &
                       2,   1,   1,   1,   1,   2,   1,   1,   1,   2,   1,   1,   1,   1, &
                       2,   1,   1,   1,   2,   1,   1,   1,   2,   1,   1,   1,   1,   2, &
                       1,   1,   1,   2,   1,   1,   1,   2,   1,   1,   1,   2,   1,   1, &
                       1,   2,   1,   1,   1,   2,   1,   1,   2,   1,   1,   1,   2,   1, &
                       1,   1,   2,   1,   1,   1,   2,   1,   1,   2,   1,   1,   1,   2, &
                       1,   1,   2,   1,   1,   1,   2,   1,   1,   2,   1,   1,   1,   2, &
                       1,   1,   2,   1,   1,   1,   2,   1,   1,   2,   1,   1,   2,   1, &
                       1,   1,   2,   1,   1,   2,   1,   1,   2,   1,   1,   2,   1,   1, &
                       2,   1,   1,   2,   1,   1,   1,   2,   1,   1,   2,   1,   1,   2, &
                       1,   1,   2,   1,   1,   2,   1,   1,   2,   1,   2,   1,   1,   2, &
                       1,   1,   2,   1,   1,   2,   1,   1,   2,   1,   1,   2,   1,   1, &
                       2,   1,   2,   1,   1,   2,   1,   1,   2,   1,   2,   1,   1,   2, &
                       1,   1,   2,   1,   2,   1,   1,   2,   1,   1,   2,   1,   2,   1, &
                       1,   2,   1,   2,   1,   1,   2,   1,   2,   1,   1,   2,   1,   2, &
                       1,   1,   2,   1,   2,   1,   1,   2,   1,   2,   1,   1,   2,   1, &
                       2,   1,   2,   1,   1,   2,   1,   2,   1,   2,   1,   1,   2,   1, &
                       2,   1,   2,   1,   1,   2,   1,   2,   1,   2,   1,   2,   1,   1, &
                       2,   1,   2,   1,   2,   1,   2,   1,   1,   2,   1,   2,   1,   2, &
                       1,   2,   1,   2,   1,   2,   1,   2,   1,   1,   2,   1,   2,   1, &
                       2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1, &
                       2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1, &
                       2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   2,   1,   2, &
                       1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   2,   1, &
                       2,   1,   2,   1,   2,   1,   2,   1,   2,   2,   1,   2,   1,   2, &
                       1,   2,   2,   1,   2,   1,   2,   1,   2,   2,   1,   2,   1,   2, &
                       1,   2,   2,   1,   2,   1,   2,   2,   1,   2,   1,   2,   2,   1, &
                       2,   1,   2,   2,   1,   2,   1,   2,   2,   1,   2,   1,   2,   2, &
                       1,   2,   1,   2,   2,   1,   2,   2,   1,   2,   1,   2,   2,   1, &
                       2,   2,   1,   2,   2,   1,   2,   1,   2,   2,   1,   2,   2,   1, &
                       2,   2,   1,   2,   2,   1,   2,   2,   1,   2,   2,   1,   2,   2, &
                       1,   2,   2,   1,   2,   2,   1,   2,   2,   1,   2,   2,   1,   2, &
                       2,   2,   1,   2,   2,   1,   2,   2,   1,   2,   2,   2,   1,   2, &
                       2,   1,   2,   2,   2,   1,   2,   2,   1,   2,   2,   2,   1,   2, &
                       2,   2,   1,   2,   2,   2,   1,   2,   2,   1,   2,   2,   2,   1, &
                       2,   2,   2,   1,   2,   2,   2,   2,   1,   2,   2,   2,   1,   2, &
                       2,   2,   1,   2,   2,   2,   2,   1,   2,   2,   2,   2,   1,   2, &
                       2,   2,   1,   2,   2,   2,   2,   1,   2,   2,   2,   2,   2,   1, &
                       2,   2,   2,   2,   1,   2,   2,   2,   2,   2,   1,   2,   2,   2, &
                       2,   2,   1,   2,   2,   2,   2,   2,   2,   1,   2,   2,   2,   2, &
                       2,   1,   2,   2,   2,   2,   2,   2,   2,   1,   2,   2,   2,   2, &
                       2,   2,   2,   1,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   2,   2, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   2, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   3, &
                       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, &
                       3,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   3,   2,   2, &
                       2,   2,   2,   2,   2,   2,   2,   3,   2,   2,   2,   2,   2,   2, &
                       3,   2,   2,   2,   2,   2,   2,   2,   3,   2,   2,   2,   2,   2, &
                       3,   2,   2,   2,   2,   2,   3,   2,   2,   2,   2,   2,   3,   2, &
                       2,   2,   2,   3,   2,   2,   2,   2,   3,   2,   2,   2,   2,   3, &
                       2,   2,   2,   3,   2,   2,   2,   3,   2,   2,   2,   2,   3,   2, &
                       2,   2,   3,   2,   2,   2,   3,   2,   2,   3,   2,   2,   2,   3, &
                       2,   2,   2,   3,   2,   2,   3,   2,   2,   3,   2,   2,   2,   3, &
                       2,   2,   3,   2,   2,   3,   2,   2,   3,   2,   2,   3,   2,   2, &
                       3,   2,   2,   3,   2,   2,   3,   2,   2,   3,   2,   2,   3,   2, &
                       3,   2,   2,   3,   2,   3,   2,   2,   3,   2,   2,   3,   2,   3, &
                       2,   2,   3,   2,   3,   2,   3,   2,   2,   3,   2,   3,   2,   3 /)
        nwvnsolar(4201:5600)=(/                                                            &
                       2,   2,   3,   2,   3,   2,   3,   2,   3,   2,   3,   2,   2,   3, &
                       2,   3,   2,   3,   2,   3,   2,   3,   2,   3,   2,   3,   2,   3, &
                       2,   3,   2,   3,   2,   3,   2,   3,   2,   3,   2,   3,   3,   2, &
                       3,   2,   3,   2,   3,   2,   3,   3,   2,   3,   2,   3,   2,   3, &
                       3,   2,   3,   2,   3,   2,   3,   3,   2,   3,   2,   3,   3,   2, &
                       3,   3,   2,   3,   2,   3,   3,   2,   3,   3,   2,   3,   3,   2, &
                       3,   2,   3,   3,   2,   3,   3,   2,   3,   3,   3,   2,   3,   3, &
                       2,   3,   3,   2,   3,   3,   2,   3,   3,   3,   2,   3,   3,   3, &
                       2,   3,   3,   2,   3,   3,   3,   2,   3,   3,   3,   3,   2,   3, &
                       3,   3,   2,   3,   3,   3,   3,   2,   3,   3,   3,   3,   2,   3, &
                       3,   3,   3,   2,   3,   3,   3,   3,   3,   2,   3,   3,   3,   3, &
                       3,   3,   2,   3,   3,   3,   3,   3,   3,   3,   2,   3,   3,   3, &
                       3,   3,   3,   3,   3,   3,   2,   3,   3,   3,   3,   3,   3,   3, &
                       3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   2,   3, &
                       3,   3,   3,   3,   3,   3,   3,   3,   3,   4,   3,   3,   3,   3, &
                       3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3,   3, &
                       4,   3,   3,   3,   3,   3,   3,   3,   3,   3,   4,   3,   3,   3, &
                       3,   3,   3,   3,   4,   3,   3,   3,   3,   3,   4,   3,   3,   3, &
                       3,   3,   4,   3,   3,   3,   3,   4,   3,   3,   3,   3,   4,   3, &
                       3,   3,   4,   3,   3,   3,   4,   3,   3,   3,   4,   3,   3,   4, &
                       3,   3,   3,   4,   3,   3,   4,   3,   3,   4,   3,   3,   4,   3, &
                       3,   4,   3,   3,   4,   3,   3,   4,   3,   3,   4,   3,   4,   3, &
                       3,   4,   3,   4,   3,   3,   4,   3,   4,   3,   3,   4,   3,   4, &
                       3,   4,   3,   4,   3,   4,   3,   3,   4,   3,   4,   3,   4,   3, &
                       4,   3,   4,   3,   4,   4,   3,   4,   3,   4,   3,   4,   3,   4, &
                       3,   4,   4,   3,   4,   3,   4,   4,   3,   4,   3,   4,   4,   3, &
                       4,   3,   4,   4,   3,   4,   4,   3,   4,   4,   3,   4,   4,   3, &
                       4,   4,   3,   4,   4,   3,   4,   4,   4,   3,   4,   4,   4,   3, &
                       4,   4,   4,   3,   4,   4,   4,   3,   4,   4,   4,   4,   3,   4, &
                       4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   3,   4,   4,   4, &
                       4,   4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   4,   4,   4, &
                       4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, &
                       4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, &
                       5,   4,   4,   4,   4,   4,   4,   4,   4,   5,   4,   4,   4,   4, &
                       4,   5,   4,   4,   4,   4,   4,   5,   4,   4,   4,   5,   4,   4, &
                       4,   5,   4,   4,   4,   5,   4,   4,   4,   5,   4,   4,   5,   4, &
                       4,   5,   4,   4,   5,   4,   4,   5,   4,   4,   5,   4,   5,   4, &
                       4,   5,   4,   5,   4,   5,   4,   4,   5,   4,   5,   4,   5,   4, &
                       5,   4,   5,   4,   5,   4,   5,   4,   5,   4,   5,   4,   5,   5, &
                       4,   5,   4,   5,   4,   5,   5,   4,   5,   5,   4,   5,   4,   5, &
                       5,   4,   5,   5,   4,   5,   5,   5,   4,   5,   5,   4,   5,   5, &
                       5,   4,   5,   5,   5,   5,   4,   5,   5,   5,   5,   4,   5,   5, &
                       5,   5,   5,   4,   5,   5,   5,   5,   5,   5,   5,   5,   5,   4, &
                       5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5,   5, &
                       5,   5,   5,   5,   5,   6,   5,   5,   5,   5,   5,   5,   5,   5, &
                       6,   5,   5,   5,   5,   5,   6,   5,   5,   5,   5,   6,   5,   5, &
                       5,   6,   5,   5,   5,   6,   5,   5,   5,   6,   5,   5,   6,   5, &
                       5,   6,   5,   6,   5,   5,   6,   5,   6,   5,   5,   6,   5,   6, &
                       5,   6,   5,   6,   5,   6,   5,   6,   5,   6,   5,   6,   5,   6, &
                       6,   5,   6,   5,   6,   6,   5,   6,   5,   6,   6,   5,   6,   6, &
                       6,   5,   6,   6,   5,   6,   6,   6,   5,   6,   6,   6,   6,   5, &
                       6,   6,   6,   6,   6,   5,   6,   6,   6,   6,   6,   6,   6,   6, &
                       6,   6,   5,   6,   6,   6,   6,   6,   6,   6,   7,   6,   6,   6, &
                       6,   6,   6,   6,   6,   6,   6,   7,   6,   6,   6,   6,   6,   7, &
                       6,   6,   6,   7,   6,   6,   6,   7,   6,   6,   6,   7,   6,   6, &
                       7,   6,   7,   6,   6,   7,   6,   7,   6,   6,   7,   6,   7,   6, &
                       7,   6,   7,   6,   7,   6,   7,   7,   6,   7,   6,   7,   7,   6, &
                       7,   6,   7,   7,   6,   7,   7,   7,   6,   7,   7,   7,   6,   7, &
                       7,   7,   7,   6,   7,   7,   7,   7,   7,   7,   7,   6,   7,   7, &
                       7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   8,   7,   7, &
                       7,   7,   7,   7,   7,   8,   7,   7,   7,   8,   7,   7,   7,   8, &
                       7,   7,   7,   8,   7,   8,   7,   7,   8,   7,   7,   8,   7,   8, &
                       7,   8,   7,   8,   7,   8,   7,   8,   8,   7,   8,   7,   8,   8, &
                       7,   8,   8,   7,   8,   8,   7,   8,   8,   8,   8,   7,   8,   8, &
                       8,   8,   8,   7,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, &
                       8,   8,   8,   8,   8,   9,   8,   8,   8,   8,   8,   9,   8,   8, &
                       8,   9,   8,   8,   8,   9,   8,   8,   9,   8,   9,   8,   8,   9, &
                       8,   9,   8,   9,   8,   9,   8,   9,   9,   8,   9,   9,   8,   9, &
                       9,   8,   9,   9,   8,   9,   9,   9,   9,   8,   9,   9,   9,   9, &
                       9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, &
                       9,  10,   9,   9,   9,  10,   9,   9,   9,  10,   9,   9,  10,   9, &
                      10,   9,   9,  10,   9,  10,   9,  10,   9,  10,  10,   9,  10,  10, &
                       9,  10,  10,   9,  10,  10,  10,   9,  10,  10,  10,  10,  10,  10, &
                      10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  11, &
                      10,  10,  10,  11,  10,  10,  11,  10,  10,  11,  10,  11,  10,  11, &
                      10,  11,  10,  11,  10,  11,  11,  10,  11,  11,  11,  10,  11,  11, &
                      11,  11,  10,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  12, &
                      11,  11,  11,  11,  12,  11,  11,  11,  12,  11,  12,  11,  11,  12, &
                      11,  12,  11,  12,  12,  11,  12,  11,  12,  12,  12,  11,  12,  12, &
                      12,  12,  12,  11,  12,  12,  12,  12,  12,  13,  12,  12,  12,  12, &
                      12,  13,  12,  12,  12,  13,  12,  13,  12,  12,  13,  12,  13,  13, &
                      12,  13,  12,  13,  13,  13,  12,  13,  13,  13,  13,  13,  13,  13, &
                      13,  13,  13,  13,  13,  13,  13,  13,  14,  13,  13,  14,  13,  13, &
                      14,  13,  14,  13,  14,  13,  14,  14,  13,  14,  14,  13,  14,  14, &
                      14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  15,  14,  14,  14, &
                      15,  14,  15,  14,  14,  15,  15,  14,  15,  14,  15,  15,  15,  14, &
                      15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  16,  15,  15,  15, &
                      16,  15,  16,  15,  16,  15,  16,  15,  16,  16,  15,  16,  16,  16, &
                      16,  16,  16,  16,  16,  16,  16,  16,  16,  17,  16,  16,  17,  16, &
                      17,  16,  17,  16,  17,  17,  16,  17,  17,  17,  17,  17,  17,  17, &
                      17,  17,  17,  17,  18,  17,  17,  18,  17,  18,  17,  18,  18,  17, &
                      18,  18,  18,  17,  18,  18,  18,  18,  19,  18,  18,  18,  18,  19, &
                      18,  19,  18,  19,  18,  19,  19,  19,  18,  19,  19,  19,  19,  19, &
                      19,  20,  19,  19,  19,  20,  19,  20,  19,  20,  20,  19,  20,  20, &
                      20,  20,  20,  20,  20,  20,  20,  21,  20,  20,  21,  20,  21,  21, &
                      20,  21,  21,  21,  21,  21,  21,  21,  21,  21,  13,   9,  21,  22, &
                      21,  22,  21,  22,  22,  21,  22,  22,  22,  22,  23,  22,  22,  22, &
                      23,  22,  23,  22,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23, &
                      24,  23,  24,  23,  24,  24,  24,  23,  24,  24,  25,  24,  24,  24, &
                      25,  24,  25,  24,  25,  25,  25,  25,  25,  25,  25,  25,  26,  25 /)
        nwvnsolar(5601:6064)=(/                                                            &
                      26,  25,  26,  26,  26,  26,  26,  26,  26,  26,  26,  27,  26,  27, &
                      27,  27,  26,  27,  27,  28,  27,  27,  28,  27,  28,  27,  28,  20, &
                       8,  28,  28,  28,  28,  29,  28,  29,  28,  29,  29,  29,  29,  29, &
                      29,  29,  30,  29,  30,  30,  30,  30,  30,  30,  30,  30,  31,  30, &
                      31,  31,  31,  31,  31,  31,  31,  32,  31,  32,  32,  32,  32,  32, &
                      32,  32,  33,  32,  33,  33,  33,  33,  33,  33,  34,  33,  34,  34, &
                      34,  34,  34,  34,  35,  34,  35,  35,  35,  35,  35,  35,  36,  35, &
                      36,  36,  36,  36,  36,  37,  36,  37,  37,  37,  37,  37,  38,  37, &
                      38,  38,  38,  38,  38,  39,  38,  39,  39,  39,  39,  39,  40,  40, &
                      40,  20,  20,  40,  40,  40,  41,  41,  41,  41,  41,  42,  42,  41, &
                      42,  42,  43,  42,  43,  43,  43,  43,  44,  43,  44,  44,  44,  44, &
                      45,  45,  44,  45,  46,  45,  46,  46,  46,  46,  46,  47,  47,  47, &
                      47,  48,  47,  48,  48,  49,  48,  49,  49,  49,  49,  50,   4,  46, &
                      50,  50,  51,  50,  51,  52,  51,  52,  52,  52,  52,  53,  52,  54, &
                      53,  53,  54,  54,  55,  54,  55,  55,  56,  55,  56,  56,  57,  56, &
                      57,  58,  57,  58,  58,  58,  59,  59,  59,  60,  59,  60,  61,  60, &
                      61,  62,  61,  62,  62,  63,  63,  63,  63,  64,  64,  64,  65,  65, &
                      66,  65,  67,  66,  67,  67,  67,  68,  68,  69,  69,  69,  70,  70, &
                      70,  71,  71,  71,  72,  72,  73,  73,  74,  73,  75,  74,  75,  66, &
                      10,  76,  76,  77,  77,  78,  78,  78,  79,  79,  80,  80,  81,  81, &
                      82,  82,  82,  84,  83,  84,  85,  85,  85,  86,  87,  87,  87,  88, &
                      89,  89,  90,  90,  90,  92,  92,  92,  93,  93,  94,  95,  95,  96, &
                      97,  97,  97,  99,  99,  99, 100, 101, 101, 103, 102, 104, 104,  90, &
                      14, 106, 106, 107, 107, 108, 109, 110, 110,  23,  88, 112, 113, 113, &
                     114, 115, 116, 117, 117, 118,  77,  42, 120, 120, 122, 122, 123, 124, &
                      27,  98, 125, 127, 128, 128, 130, 130, 131, 132,  71,  63, 134, 135, &
                     136, 137, 138, 140, 140, 141, 143, 143, 145, 145, 147, 148, 149, 150, &
                     152, 152, 154, 155, 156, 158, 159,  80,  80, 161, 163, 163, 166, 166, &
                     168, 170, 170, 173, 173, 175, 177, 178, 179, 181, 183, 184, 186, 104, &
                      83, 189, 191, 193, 194, 195, 198, 199, 201, 203, 205, 207, 208, 211, &
                     212, 214, 217, 218, 220, 223, 224, 227, 229, 231, 234, 235, 238, 241, &
                     242, 245, 248, 250, 252, 255, 258, 260, 263, 266, 269, 271, 274, 277, &
                     280, 283, 286, 289, 292, 295, 299, 302, 305, 309, 312, 316, 319, 323, &
                     326, 294 /)
        solint(1:700)=(/                                                                                                   &
                    0.00000E+00,   0.27147E-04,   0.27587E-04,   0.28238E-04,   0.28696E-04,   0.29369E-04,   0.29847E-04, &
                    0.30337E-04,   0.30836E-04,   0.31550E-04,   0.31868E-04,   0.32601E-04,   0.33141E-04,   0.33491E-04, &
                    0.34255E-04,   0.34625E-04,   0.35408E-04,   0.36001E-04,   0.36403E-04,   0.37016E-04,   0.37841E-04, &
                    0.38273E-04,   0.38918E-04,   0.39572E-04,   0.40237E-04,   0.40911E-04,   0.41199E-04,   0.42092E-04, &
                    0.42599E-04,   0.43316E-04,   0.43842E-04,   0.44777E-04,   0.44930E-04,   0.46282E-04,   0.46452E-04, &
                    0.47033E-04,   0.47624E-04,   0.48616E-04,   0.49623E-04,   0.49848E-04,   0.50871E-04,   0.51513E-04, &
                    0.52163E-04,   0.52825E-04,   0.53497E-04,   0.54572E-04,   0.54870E-04,   0.55572E-04,   0.56283E-04, &
                    0.57400E-04,   0.57738E-04,   0.58876E-04,   0.59630E-04,   0.60001E-04,   0.60773E-04,   0.61948E-04, &
                    0.62740E-04,   0.63154E-04,   0.63967E-04,   0.65184E-04,   0.65626E-04,   0.66469E-04,   0.67326E-04, &
                    0.68189E-04,   0.69063E-04,   0.69951E-04,   0.70454E-04,   0.71361E-04,   0.72274E-04,   0.73199E-04, &
                    0.73744E-04,   0.75081E-04,   0.75646E-04,   0.76617E-04,   0.77203E-04,   0.78186E-04,   0.79184E-04, &
                    0.79802E-04,   0.81206E-04,   0.81844E-04,   0.82490E-04,   0.83537E-04,   0.84596E-04,   0.85273E-04, &
                    0.86350E-04,   0.87439E-04,   0.88540E-04,   0.89254E-04,   0.90376E-04,   0.91111E-04,   0.91856E-04, &
                    0.92617E-04,   0.94167E-04,   0.94946E-04,   0.95737E-04,   0.96928E-04,   0.97744E-04,   0.98957E-04, &
                    0.10018E-03,   0.10102E-03,   0.10187E-03,   0.10274E-03,   0.10400E-03,   0.10488E-03,   0.10616E-03, &
                    0.10708E-03,   0.10798E-03,   0.10891E-03,   0.10985E-03,   0.11080E-03,   0.11214E-03,   0.11350E-03, &
                    0.11408E-03,   0.11508E-03,   0.11646E-03,   0.11747E-03,   0.11888E-03,   0.11952E-03,   0.12055E-03, &
                    0.12161E-03,   0.12304E-03,   0.12411E-03,   0.12519E-03,   0.12628E-03,   0.12777E-03,   0.12888E-03, &
                    0.12999E-03,   0.13073E-03,   0.13187E-03,   0.13302E-03,   0.13418E-03,   0.13535E-03,   0.13653E-03, &
                    0.13772E-03,   0.13891E-03,   0.14012E-03,   0.14134E-03,   0.14257E-03,   0.14381E-03,   0.14467E-03, &
                    0.14593E-03,   0.14720E-03,   0.14887E-03,   0.14978E-03,   0.15068E-03,   0.15238E-03,   0.15332E-03, &
                    0.15463E-03,   0.15558E-03,   0.15732E-03,   0.15830E-03,   0.15967E-03,   0.16066E-03,   0.16244E-03, &
                    0.16346E-03,   0.16486E-03,   0.16628E-03,   0.16771E-03,   0.16876E-03,   0.17061E-03,   0.17089E-03, &
                    0.17276E-03,   0.17385E-03,   0.17495E-03,   0.17606E-03,   0.17796E-03,   0.17987E-03,   0.18024E-03, &
                    0.18140E-03,   0.18332E-03,   0.18450E-03,   0.18569E-03,   0.18688E-03,   0.18885E-03,   0.19007E-03, &
                    0.19128E-03,   0.19252E-03,   0.19377E-03,   0.19580E-03,   0.19629E-03,   0.19834E-03,   0.19962E-03, &
                    0.20168E-03,   0.20220E-03,   0.20430E-03,   0.20561E-03,   0.20696E-03,   0.20907E-03,   0.21042E-03, &
                    0.21101E-03,   0.21315E-03,   0.21455E-03,   0.21592E-03,   0.21733E-03,   0.21874E-03,   0.21941E-03, &
                    0.22160E-03,   0.22305E-03,   0.22451E-03,   0.22596E-03,   0.22745E-03,   0.22890E-03,   0.23040E-03, &
                    0.23192E-03,   0.23342E-03,   0.23572E-03,   0.23647E-03,   0.23803E-03,   0.23956E-03,   0.24115E-03, &
                    0.24271E-03,   0.24508E-03,   0.24590E-03,   0.24748E-03,   0.24910E-03,   0.25072E-03,   0.25234E-03, &
                    0.25399E-03,   0.25487E-03,   0.25732E-03,   0.25822E-03,   0.25989E-03,   0.26159E-03,   0.26328E-03, &
                    0.26501E-03,   0.26674E-03,   0.26769E-03,   0.27021E-03,   0.27195E-03,   0.27295E-03,   0.27472E-03, &
                    0.27652E-03,   0.27755E-03,   0.27933E-03,   0.28114E-03,   0.28295E-03,   0.28479E-03,   0.28662E-03, &
                    0.28848E-03,   0.29033E-03,   0.29143E-03,   0.29330E-03,   0.29520E-03,   0.29709E-03,   0.29824E-03, &
                    0.30015E-03,   0.30209E-03,   0.30406E-03,   0.30597E-03,   0.30719E-03,   0.30839E-03,   0.31113E-03, &
                    0.31235E-03,   0.31355E-03,   0.31555E-03,   0.31759E-03,   0.31962E-03,   0.32167E-03,   0.32294E-03, &
                    0.32501E-03,   0.32629E-03,   0.32837E-03,   0.32971E-03,   0.33257E-03,   0.33388E-03,   0.33526E-03, &
                    0.33740E-03,   0.33951E-03,   0.34166E-03,   0.34306E-03,   0.34522E-03,   0.34663E-03,   0.34880E-03, &
                    0.35023E-03,   0.35246E-03,   0.35389E-03,   0.35611E-03,   0.35761E-03,   0.35984E-03,   0.36134E-03, &
                    0.36358E-03,   0.36585E-03,   0.36738E-03,   0.36965E-03,   0.37042E-03,   0.37353E-03,   0.37507E-03, &
                    0.37664E-03,   0.37901E-03,   0.38058E-03,   0.38295E-03,   0.38459E-03,   0.38619E-03,   0.38858E-03, &
                    0.39018E-03,   0.39188E-03,   0.39424E-03,   0.39670E-03,   0.39759E-03,   0.40004E-03,   0.40253E-03, &
                    0.40421E-03,   0.40592E-03,   0.40767E-03,   0.41021E-03,   0.41272E-03,   0.41442E-03,   0.41622E-03, &
                    0.41798E-03,   0.42053E-03,   0.42235E-03,   0.42413E-03,   0.42594E-03,   0.42854E-03,   0.43111E-03, &
                    0.43224E-03,   0.43487E-03,   0.43669E-03,   0.43861E-03,   0.44126E-03,   0.44309E-03,   0.44579E-03, &
                    0.44768E-03,   0.44960E-03,   0.45155E-03,   0.45353E-03,   0.45623E-03,   0.45819E-03,   0.46019E-03, &
                    0.46214E-03,   0.46412E-03,   0.46612E-03,   0.46892E-03,   0.47014E-03,   0.47224E-03,   0.47496E-03, &
                    0.47703E-03,   0.47913E-03,   0.48119E-03,   0.48404E-03,   0.48615E-03,   0.48820E-03,   0.49037E-03, &
                    0.49248E-03,   0.49453E-03,   0.49670E-03,   0.49882E-03,   0.50105E-03,   0.50323E-03,   0.50534E-03, &
                    0.50834E-03,   0.51052E-03,   0.51195E-03,   0.51494E-03,   0.51720E-03,   0.51940E-03,   0.52163E-03, &
                    0.52388E-03,   0.52541E-03,   0.52849E-03,   0.53075E-03,   0.53303E-03,   0.53458E-03,   0.53768E-03, &
                    0.53994E-03,   0.54158E-03,   0.54390E-03,   0.54627E-03,   0.54932E-03,   0.55098E-03,   0.55332E-03, &
                    0.55580E-03,   0.55821E-03,   0.56054E-03,   0.56300E-03,   0.56539E-03,   0.56715E-03,   0.56960E-03, &
                    0.57198E-03,   0.57448E-03,   0.57702E-03,   0.57948E-03,   0.58197E-03,   0.58371E-03,   0.58702E-03, &
                    0.58872E-03,   0.59133E-03,   0.59385E-03,   0.59565E-03,   0.59825E-03,   0.60087E-03,   0.60341E-03, &
                    0.60608E-03,   0.60791E-03,   0.61053E-03,   0.61318E-03,   0.61587E-03,   0.61770E-03,   0.62121E-03, &
                    0.62309E-03,   0.62501E-03,   0.62771E-03,   0.63045E-03,   0.63309E-03,   0.63513E-03,   0.63859E-03, &
                    0.64056E-03,   0.64332E-03,   0.64535E-03,   0.64818E-03,   0.65090E-03,   0.65378E-03,   0.65577E-03, &
                    0.65855E-03,   0.66060E-03,   0.66270E-03,   0.66621E-03,   0.66837E-03,   0.67195E-03,   0.67402E-03, &
                    0.67612E-03,   0.67826E-03,   0.68043E-03,   0.68414E-03,   0.68623E-03,   0.68834E-03,   0.69216E-03, &
                    0.69432E-03,   0.69805E-03,   0.70014E-03,   0.70239E-03,   0.70468E-03,   0.70838E-03,   0.71057E-03, &
                    0.71280E-03,   0.71506E-03,   0.71735E-03,   0.71968E-03,   0.72357E-03,   0.72580E-03,   0.72806E-03, &
                    0.73202E-03,   0.73433E-03,   0.73667E-03,   0.73904E-03,   0.74281E-03,   0.74524E-03,   0.74754E-03, &
                    0.75003E-03,   0.75391E-03,   0.75477E-03,   0.75871E-03,   0.75962E-03,   0.76362E-03,   0.76595E-03, &
                    0.77001E-03,   0.77240E-03,   0.77347E-03,   0.77745E-03,   0.77993E-03,   0.78397E-03,   0.78499E-03, &
                    0.78738E-03,   0.79150E-03,   0.79413E-03,   0.79661E-03,   0.79913E-03,   0.80166E-03,   0.80424E-03, &
                    0.80836E-03,   0.81099E-03,   0.81364E-03,   0.81615E-03,   0.81886E-03,   0.82295E-03,   0.82420E-03, &
                    0.82835E-03,   0.83100E-03,   0.83216E-03,   0.83639E-03,   0.83761E-03,   0.84170E-03,   0.84450E-03, &
                    0.84866E-03,   0.84999E-03,   0.85268E-03,   0.85540E-03,   0.85967E-03,   0.86245E-03,   0.86373E-03, &
                    0.86809E-03,   0.87096E-03,   0.87365E-03,   0.87658E-03,   0.87933E-03,   0.88211E-03,   0.88491E-03, &
                    0.88796E-03,   0.89083E-03,   0.89525E-03,   0.89796E-03,   0.89940E-03,   0.90390E-03,   0.90518E-03, &
                    0.90821E-03,   0.91258E-03,   0.91546E-03,   0.91860E-03,   0.92002E-03,   0.92452E-03,   0.92601E-03, &
                    0.92881E-03,   0.93187E-03,   0.93648E-03,   0.93786E-03,   0.94102E-03,   0.94397E-03,   0.94718E-03, &
                    0.95019E-03,   0.95322E-03,   0.95629E-03,   0.95938E-03,   0.96250E-03,   0.96568E-03,   0.96865E-03, &
                    0.97187E-03,   0.97513E-03,   0.97816E-03,   0.98126E-03,   0.98463E-03,   0.98778E-03,   0.99095E-03, &
                    0.99414E-03,   0.99735E-03,   0.10006E-02,   0.10024E-02,   0.10057E-02,   0.10087E-02,   0.10121E-02, &
                    0.10155E-02,   0.10171E-02,   0.10218E-02,   0.10253E-02,   0.10286E-02,   0.10303E-02,   0.10351E-02, &
                    0.10369E-02,   0.10403E-02,   0.10437E-02,   0.10471E-02,   0.10502E-02,   0.10537E-02,   0.10572E-02, &
                    0.10604E-02,   0.10622E-02,   0.10658E-02,   0.10691E-02,   0.10709E-02,   0.10743E-02,   0.10778E-02, &
                    0.10812E-02,   0.10847E-02,   0.10882E-02,   0.10902E-02,   0.10953E-02,   0.10986E-02,   0.11022E-02, &
                    0.11040E-02,   0.11077E-02,   0.11112E-02,   0.11149E-02,   0.11169E-02,   0.11204E-02,   0.11240E-02, &
                    0.11261E-02,   0.11312E-02,   0.11334E-02,   0.11371E-02,   0.11408E-02,   0.11427E-02,   0.11465E-02, &
                    0.11503E-02,   0.11539E-02,   0.11562E-02,   0.11598E-02,   0.11634E-02,   0.11674E-02,   0.11711E-02, &
                    0.11733E-02,   0.11770E-02,   0.11807E-02,   0.11845E-02,   0.11869E-02,   0.11907E-02,   0.11945E-02, &
                    0.11969E-02,   0.12005E-02,   0.12045E-02,   0.12081E-02,   0.12106E-02,   0.12142E-02,   0.12183E-02, &
                    0.12220E-02,   0.12258E-02,   0.12281E-02,   0.12308E-02,   0.12361E-02,   0.12385E-02,   0.12424E-02, &
                    0.12449E-02,   0.12485E-02,   0.12510E-02,   0.12550E-02,   0.12590E-02,   0.12628E-02,   0.12654E-02, &
                    0.12695E-02,   0.12733E-02,   0.12757E-02,   0.12799E-02,   0.12838E-02,   0.12863E-02,   0.12907E-02, &
                    0.12932E-02,   0.12973E-02,   0.13013E-02,   0.13054E-02,   0.13081E-02,   0.13122E-02,   0.13149E-02, &
                    0.13188E-02,   0.13215E-02,   0.13273E-02,   0.13297E-02,   0.13325E-02,   0.13369E-02,   0.13409E-02, &
                    0.13434E-02,   0.13479E-02,   0.13505E-02,   0.13546E-02,   0.13576E-02,   0.13618E-02,   0.13660E-02, &
                    0.13687E-02,   0.13730E-02,   0.13757E-02,   0.13800E-02,   0.13828E-02,   0.13857E-02,   0.13900E-02, &
                    0.13940E-02,   0.13970E-02,   0.14015E-02,   0.14040E-02,   0.14086E-02,   0.14127E-02,   0.14158E-02, &
                    0.14185E-02,   0.14231E-02,   0.14258E-02,   0.14301E-02,   0.14334E-02,   0.14392E-02,   0.14421E-02, &
                    0.14451E-02,   0.14480E-02,   0.14540E-02,   0.14570E-02,   0.14601E-02,   0.14646E-02,   0.14678E-02, &
                    0.14709E-02,   0.14751E-02,   0.14798E-02,   0.14830E-02,   0.14858E-02,   0.14906E-02,   0.14954E-02, &
                    0.14982E-02,   0.15011E-02,   0.15061E-02,   0.15090E-02,   0.15140E-02,   0.15169E-02,   0.15200E-02, &
                    0.15246E-02,   0.15281E-02,   0.15328E-02,   0.15359E-02,   0.15405E-02,   0.15437E-02,   0.15469E-02, &
                    0.15516E-02,   0.15549E-02,   0.15592E-02,   0.15625E-02,   0.15658E-02,   0.15707E-02,   0.15736E-02 /)
        solint(701:1400)=(/                                                                                                 &
                    0.15785E-02,   0.15820E-02,   0.15850E-02,   0.15900E-02,   0.15930E-02,   0.15981E-02,   0.15998E-02, &
                    0.16044E-02,   0.16082E-02,   0.16129E-02,   0.16177E-02,   0.16200E-02,   0.16248E-02,   0.16282E-02, &
                    0.16331E-02,   0.16365E-02,   0.16400E-02,   0.16434E-02,   0.16485E-02,   0.16520E-02,   0.16571E-02, &
                    0.16607E-02,   0.16643E-02,   0.16674E-02,   0.16711E-02,   0.16763E-02,   0.16795E-02,   0.16848E-02, &
                    0.16870E-02,   0.16918E-02,   0.16956E-02,   0.17004E-02,   0.17043E-02,   0.17076E-02,   0.17130E-02, &
                    0.17164E-02,   0.17198E-02,   0.17238E-02,   0.17273E-02,   0.17323E-02,   0.17358E-02,   0.17409E-02, &
                    0.17436E-02,   0.17472E-02,   0.17523E-02,   0.17559E-02,   0.17611E-02,   0.17649E-02,   0.17672E-02, &
                    0.17735E-02,   0.17758E-02,   0.17812E-02,   0.17836E-02,   0.17890E-02,   0.17924E-02,   0.17964E-02, &
                    0.18004E-02,   0.18054E-02,   0.18095E-02,   0.18136E-02,   0.18172E-02,   0.18228E-02,   0.18264E-02, &
                    0.18291E-02,   0.18342E-02,   0.18394E-02,   0.18422E-02,   0.18459E-02,   0.18512E-02,   0.18556E-02, &
                    0.18594E-02,   0.18632E-02,   0.18669E-02,   0.18707E-02,   0.18746E-02,   0.18785E-02,   0.18825E-02, &
                    0.18880E-02,   0.18920E-02,   0.18975E-02,   0.19001E-02,   0.19057E-02,   0.19083E-02,   0.19125E-02, &
                    0.19181E-02,   0.19223E-02,   0.19259E-02,   0.19302E-02,   0.19346E-02,   0.19383E-02,   0.19427E-02, &
                    0.19471E-02,   0.19509E-02,   0.19569E-02,   0.19592E-02,   0.19652E-02,   0.19691E-02,   0.19723E-02, &
                    0.19778E-02,   0.19824E-02,   0.19864E-02,   0.19904E-02,   0.19951E-02,   0.19992E-02,   0.20033E-02, &
                    0.20074E-02,   0.20115E-02,   0.20164E-02,   0.20207E-02,   0.20264E-02,   0.20306E-02,   0.20334E-02, &
                    0.20376E-02,   0.20419E-02,   0.20462E-02,   0.20521E-02,   0.20565E-02,   0.20594E-02,   0.20653E-02, &
                    0.20684E-02,   0.20737E-02,   0.20768E-02,   0.20829E-02,   0.20860E-02,   0.20900E-02,   0.20947E-02, &
                    0.20994E-02,   0.21041E-02,   0.21081E-02,   0.21129E-02,   0.21170E-02,   0.21219E-02,   0.21260E-02, &
                    0.21309E-02,   0.21351E-02,   0.21401E-02,   0.21443E-02,   0.21494E-02,   0.21537E-02,   0.21564E-02, &
                    0.21615E-02,   0.21659E-02,   0.21703E-02,   0.21755E-02,   0.21799E-02,   0.21844E-02,   0.21888E-02, &
                    0.21933E-02,   0.21978E-02,   0.22017E-02,   0.22078E-02,   0.22109E-02,   0.22156E-02,   0.22203E-02, &
                    0.22250E-02,   0.22298E-02,   0.22331E-02,   0.22379E-02,   0.22427E-02,   0.22482E-02,   0.22517E-02, &
                    0.22566E-02,   0.22616E-02,   0.22651E-02,   0.22716E-02,   0.22759E-02,   0.22795E-02,   0.22846E-02, &
                    0.22874E-02,   0.22940E-02,   0.22977E-02,   0.23020E-02,   0.23073E-02,   0.23126E-02,   0.23155E-02, &
                    0.23209E-02,   0.23254E-02,   0.23308E-02,   0.23353E-02,   0.23398E-02,   0.23452E-02,   0.23484E-02, &
                    0.23539E-02,   0.23570E-02,   0.23632E-02,   0.23673E-02,   0.23720E-02,   0.23768E-02,   0.23824E-02, &
                    0.23872E-02,   0.23905E-02,   0.23953E-02,   0.24001E-02,   0.24044E-02,   0.24093E-02,   0.24142E-02, &
                    0.24177E-02,   0.24241E-02,   0.24291E-02,   0.24327E-02,   0.24378E-02,   0.24428E-02,   0.24465E-02, &
                    0.24517E-02,   0.24554E-02,   0.24606E-02,   0.24644E-02,   0.24711E-02,   0.24739E-02,   0.24791E-02, &
                    0.24844E-02,   0.24883E-02,   0.24937E-02,   0.24982E-02,   0.25037E-02,   0.25076E-02,   0.25122E-02, &
                    0.25178E-02,   0.25219E-02,   0.25266E-02,   0.25309E-02,   0.25366E-02,   0.25399E-02,   0.25457E-02, &
                    0.25520E-02,   0.25579E-02,   0.25613E-02,   0.25657E-02,   0.25721E-02,   0.25736E-02,   0.25800E-02, &
                    0.25846E-02,   0.25881E-02,   0.25947E-02,   0.25995E-02,   0.26031E-02,   0.26097E-02,   0.26144E-02, &
                    0.26181E-02,   0.26248E-02,   0.26297E-02,   0.26335E-02,   0.26373E-02,   0.26441E-02,   0.26509E-02, &
                    0.26528E-02,   0.26597E-02,   0.26637E-02,   0.26676E-02,   0.26745E-02,   0.26756E-02,   0.26826E-02, &
                    0.26896E-02,   0.26937E-02,   0.26979E-02,   0.27021E-02,   0.27063E-02,   0.27135E-02,   0.27177E-02, &
                    0.27220E-02,   0.27263E-02,   0.27336E-02,   0.27381E-02,   0.27425E-02,   0.27499E-02,   0.27532E-02, &
                    0.27577E-02,   0.27623E-02,   0.27669E-02,   0.27732E-02,   0.27778E-02,   0.27825E-02,   0.27901E-02, &
                    0.27936E-02,   0.27983E-02,   0.28030E-02,   0.28096E-02,   0.28114E-02,   0.28191E-02,   0.28227E-02, &
                    0.28275E-02,   0.28312E-02,   0.28389E-02,   0.28426E-02,   0.28505E-02,   0.28513E-02,   0.28593E-02, &
                    0.28631E-02,   0.28680E-02,   0.28719E-02,   0.28800E-02,   0.28840E-02,   0.28881E-02,   0.28963E-02, &
                    0.28973E-02,   0.29055E-02,   0.29096E-02,   0.29167E-02,   0.29221E-02,   0.29233E-02,   0.29275E-02, &
                    0.29345E-02,   0.29400E-02,   0.29442E-02,   0.29515E-02,   0.29559E-02,   0.29602E-02,   0.29661E-02, &
                    0.29706E-02,   0.29781E-02,   0.29797E-02,   0.29842E-02,   0.29918E-02,   0.29993E-02,   0.30037E-02, &
                    0.30080E-02,   0.30125E-02,   0.30171E-02,   0.30218E-02,   0.30267E-02,   0.30346E-02,   0.30365E-02, &
                    0.30442E-02,   0.30489E-02,   0.30538E-02,   0.30588E-02,   0.30639E-02,   0.30719E-02,   0.30770E-02, &
                    0.30790E-02,   0.30856E-02,   0.30907E-02,   0.30988E-02,   0.31010E-02,   0.31091E-02,   0.31100E-02, &
                    0.31182E-02,   0.31235E-02,   0.31290E-02,   0.31361E-02,   0.31388E-02,   0.31473E-02,   0.31513E-02, &
                    0.31567E-02,   0.31592E-02,   0.31662E-02,   0.31716E-02,   0.31772E-02,   0.31845E-02,   0.31872E-02, &
                    0.31943E-02,   0.32000E-02,   0.32042E-02,   0.32100E-02,   0.32144E-02,   0.32204E-02,   0.32248E-02, &
                    0.32307E-02,   0.32381E-02,   0.32410E-02,   0.32483E-02,   0.32543E-02,   0.32588E-02,   0.32649E-02, &
                    0.32724E-02,   0.32770E-02,   0.32801E-02,   0.32848E-02,   0.32895E-02,   0.32957E-02,   0.33034E-02, &
                    0.33053E-02,   0.33146E-02,   0.33194E-02,   0.33243E-02,   0.33321E-02,   0.33356E-02,   0.33405E-02, &
                    0.33454E-02,   0.33503E-02,   0.33553E-02,   0.33617E-02,   0.33696E-02,   0.33747E-02,   0.33798E-02, &
                    0.33849E-02,   0.33898E-02,   0.33948E-02,   0.33999E-02,   0.34067E-02,   0.34120E-02,   0.34174E-02, &
                    0.34227E-02,   0.34280E-02,   0.34335E-02,   0.34390E-02,   0.34444E-02,   0.34499E-02,   0.34584E-02, &
                    0.34610E-02,   0.34665E-02,   0.34750E-02,   0.34806E-02,   0.34816E-02,   0.34874E-02,   0.34931E-02, &
                    0.34990E-02,   0.35048E-02,   0.35105E-02,   0.35190E-02,   0.35219E-02,   0.35291E-02,   0.35352E-02, &
                    0.35382E-02,   0.35468E-02,   0.35525E-02,   0.35567E-02,   0.35626E-02,   0.35687E-02,   0.35748E-02, &
                    0.35793E-02,   0.35825E-02,   0.35917E-02,   0.35961E-02,   0.36023E-02,   0.36057E-02,   0.36102E-02, &
                    0.36192E-02,   0.36254E-02,   0.36300E-02,   0.36364E-02,   0.36399E-02,   0.36474E-02,   0.36511E-02, &
                    0.36589E-02,   0.36656E-02,   0.36674E-02,   0.36741E-02,   0.36790E-02,   0.36889E-02,   0.36941E-02, &
                    0.37011E-02,   0.37032E-02,   0.37100E-02,   0.37150E-02,   0.37218E-02,   0.37268E-02,   0.37337E-02, &
                    0.37388E-02,   0.37458E-02,   0.37509E-02,   0.37563E-02,   0.37637E-02,   0.37691E-02,   0.37733E-02, &
                    0.37786E-02,   0.37870E-02,   0.37944E-02,   0.37999E-02,   0.38024E-02,   0.38109E-02,   0.38183E-02, &
                    0.38239E-02,   0.38294E-02,   0.38339E-02,   0.38395E-02,   0.38451E-02,   0.38507E-02,   0.38583E-02, &
                    0.38610E-02,   0.38696E-02,   0.38753E-02,   0.38811E-02,   0.38859E-02,   0.38947E-02,   0.38975E-02, &
                    0.39063E-02,   0.39121E-02,   0.39151E-02,   0.39210E-02,   0.39270E-02,   0.39320E-02,   0.39409E-02, &
                    0.39440E-02,   0.39530E-02,   0.39590E-02,   0.39619E-02,   0.39678E-02,   0.39738E-02,   0.39800E-02, &
                    0.39863E-02,   0.39897E-02,   0.39989E-02,   0.40052E-02,   0.40086E-02,   0.40150E-02,   0.40215E-02, &
                    0.40279E-02,   0.40342E-02,   0.40406E-02,   0.40471E-02,   0.40537E-02,   0.40574E-02,   0.40619E-02, &
                    0.40685E-02,   0.40752E-02,   0.40790E-02,   0.40856E-02,   0.40952E-02,   0.40990E-02,   0.41068E-02, &
                    0.41109E-02,   0.41180E-02,   0.41219E-02,   0.41287E-02,   0.41334E-02,   0.41433E-02,   0.41472E-02, &
                    0.41541E-02,   0.41588E-02,   0.41620E-02,   0.41708E-02,   0.41745E-02,   0.41825E-02,   0.41898E-02, &
                    0.41940E-02,   0.41981E-02,   0.42026E-02,   0.42082E-02,   0.42139E-02,   0.42184E-02,   0.42226E-02, &
                    0.42304E-02,   0.42350E-02,   0.42424E-02,   0.42475E-02,   0.42549E-02,   0.42625E-02,   0.42650E-02, &
                    0.42728E-02,   0.42812E-02,   0.42860E-02,   0.42886E-02,   0.42963E-02,   0.43044E-02,   0.43107E-02, &
                    0.43197E-02,   0.43257E-02,   0.43337E-02,   0.43363E-02,   0.43439E-02,   0.43494E-02,   0.43580E-02, &
                    0.43651E-02,   0.43695E-02,   0.43783E-02,   0.43844E-02,   0.43929E-02,   0.43989E-02,   0.44042E-02, &
                    0.44073E-02,   0.44134E-02,   0.44218E-02,   0.44307E-02,   0.44360E-02,   0.44419E-02,   0.44477E-02, &
                    0.44531E-02,   0.44592E-02,   0.44653E-02,   0.44710E-02,   0.44772E-02,   0.44834E-02,   0.44918E-02, &
                    0.44950E-02,   0.45015E-02,   0.45106E-02,   0.45142E-02,   0.45205E-02,   0.45268E-02,   0.45356E-02, &
                    0.45391E-02,   0.45485E-02,   0.45521E-02,   0.45611E-02,   0.45677E-02,   0.45743E-02,   0.45781E-02, &
                    0.45848E-02,   0.45911E-02,   0.45979E-02,   0.46046E-02,   0.46085E-02,   0.46153E-02,   0.46222E-02, &
                    0.46286E-02,   0.46383E-02,   0.46423E-02,   0.46493E-02,   0.46564E-02,   0.46635E-02,   0.46676E-02, &
                    0.46745E-02,   0.46815E-02,   0.46855E-02,   0.46922E-02,   0.46993E-02,   0.47065E-02,   0.47137E-02, &
                    0.47180E-02,   0.47251E-02,   0.47319E-02,   0.47355E-02,   0.47449E-02,   0.47490E-02,   0.47563E-02, &
                    0.47610E-02,   0.47685E-02,   0.47759E-02,   0.47804E-02,   0.47879E-02,   0.47925E-02,   0.48030E-02, &
                    0.48078E-02,   0.48125E-02,   0.48203E-02,   0.48279E-02,   0.48326E-02,   0.48401E-02,   0.48447E-02, &
                    0.48523E-02,   0.48601E-02,   0.48649E-02,   0.48698E-02,   0.48779E-02,   0.48828E-02,   0.48909E-02, &
                    0.48995E-02,   0.49053E-02,   0.49107E-02,   0.49187E-02,   0.49239E-02,   0.49289E-02,   0.49369E-02, &
                    0.49451E-02,   0.49504E-02,   0.49558E-02,   0.49640E-02,   0.49693E-02,   0.49746E-02,   0.49830E-02, &
                    0.49913E-02,   0.49968E-02,   0.50055E-02,   0.50085E-02,   0.50171E-02,   0.50256E-02,   0.50282E-02, &
                    0.50338E-02,   0.50423E-02,   0.50508E-02,   0.50565E-02,   0.50623E-02,   0.50707E-02,   0.50731E-02, &
                    0.50813E-02,   0.50868E-02,   0.50925E-02,   0.50980E-02,   0.51066E-02,   0.51124E-02,   0.51185E-02, &
                    0.51245E-02,   0.51332E-02,   0.51392E-02,   0.51451E-02,   0.51540E-02,   0.51602E-02,   0.51634E-02, &
                    0.51723E-02,   0.51754E-02,   0.51843E-02,   0.51906E-02,   0.51971E-02,   0.52034E-02,   0.52129E-02, &
                    0.52194E-02,   0.52259E-02,   0.52328E-02,   0.52397E-02,   0.52471E-02,   0.52539E-02,   0.52609E-02, &
                    0.52678E-02,   0.52748E-02,   0.52812E-02,   0.52878E-02,   0.52942E-02,   0.53009E-02,   0.53109E-02 /)
        solint(1401:2100)=(/                                                                                               &
                    0.53149E-02,   0.53219E-02,   0.53285E-02,   0.53354E-02,   0.53422E-02,   0.53493E-02,   0.53567E-02, &
                    0.53636E-02,   0.53706E-02,   0.53774E-02,   0.53814E-02,   0.53915E-02,   0.53955E-02,   0.54054E-02, &
                    0.54093E-02,   0.54161E-02,   0.54234E-02,   0.54305E-02,   0.54347E-02,   0.54421E-02,   0.54493E-02, &
                    0.54564E-02,   0.54640E-02,   0.54712E-02,   0.54783E-02,   0.54856E-02,   0.54928E-02,   0.54974E-02, &
                    0.55049E-02,   0.55086E-02,   0.55145E-02,   0.55207E-02,   0.55273E-02,   0.55339E-02,   0.55413E-02, &
                    0.55459E-02,   0.55538E-02,   0.55584E-02,   0.55658E-02,   0.55733E-02,   0.55780E-02,   0.55889E-02, &
                    0.55934E-02,   0.56009E-02,   0.56057E-02,   0.56141E-02,   0.56221E-02,   0.56271E-02,   0.56347E-02, &
                    0.56425E-02,   0.56476E-02,   0.56562E-02,   0.56611E-02,   0.56717E-02,   0.56772E-02,   0.56863E-02, &
                    0.56896E-02,   0.56986E-02,   0.57110E-02,   0.57160E-02,   0.57205E-02,   0.57280E-02,   0.57334E-02, &
                    0.57419E-02,   0.57473E-02,   0.57525E-02,   0.57607E-02,   0.57688E-02,   0.57742E-02,   0.57834E-02, &
                    0.57886E-02,   0.57958E-02,   0.58002E-02,   0.58081E-02,   0.58162E-02,   0.58186E-02,   0.58264E-02, &
                    0.58313E-02,   0.58394E-02,   0.58452E-02,   0.58538E-02,   0.58592E-02,   0.58677E-02,   0.58732E-02, &
                    0.58787E-02,   0.58879E-02,   0.58972E-02,   0.59002E-02,   0.59086E-02,   0.59136E-02,   0.59212E-02, &
                    0.59296E-02,   0.59356E-02,   0.59410E-02,   0.59424E-02,   0.59495E-02,   0.59551E-02,   0.59648E-02, &
                    0.59715E-02,   0.59812E-02,   0.59878E-02,   0.59942E-02,   0.60036E-02,   0.60103E-02,   0.60168E-02, &
                    0.60224E-02,   0.60291E-02,   0.60348E-02,   0.60438E-02,   0.60509E-02,   0.60608E-02,   0.60670E-02, &
                    0.60730E-02,   0.60782E-02,   0.60851E-02,   0.60932E-02,   0.61034E-02,   0.61101E-02,   0.61168E-02, &
                    0.61225E-02,   0.61302E-02,   0.61411E-02,   0.61484E-02,   0.61548E-02,   0.61566E-02,   0.61645E-02, &
                    0.61704E-02,   0.61769E-02,   0.61828E-02,   0.61892E-02,   0.61940E-02,   0.61995E-02,   0.62098E-02, &
                    0.62166E-02,   0.62211E-02,   0.62276E-02,   0.62364E-02,   0.62417E-02,   0.62518E-02,   0.62568E-02, &
                    0.62629E-02,   0.62698E-02,   0.62789E-02,   0.62829E-02,   0.62900E-02,   0.62985E-02,   0.63058E-02, &
                    0.63115E-02,   0.63182E-02,   0.63255E-02,   0.63321E-02,   0.63405E-02,   0.63490E-02,   0.63554E-02, &
                    0.63628E-02,   0.63709E-02,   0.63781E-02,   0.63865E-02,   0.63941E-02,   0.64002E-02,   0.64071E-02, &
                    0.64121E-02,   0.64196E-02,   0.64282E-02,   0.64344E-02,   0.64409E-02,   0.64473E-02,   0.64504E-02, &
                    0.64602E-02,   0.64628E-02,   0.64697E-02,   0.64769E-02,   0.64829E-02,   0.64913E-02,   0.64993E-02, &
                    0.65088E-02,   0.65161E-02,   0.65203E-02,   0.65269E-02,   0.65327E-02,   0.65392E-02,   0.65414E-02, &
                    0.65473E-02,   0.65525E-02,   0.65554E-02,   0.65644E-02,   0.65707E-02,   0.65779E-02,   0.65818E-02, &
                    0.65901E-02,   0.65952E-02,   0.66032E-02,   0.66142E-02,   0.66202E-02,   0.66298E-02,   0.66333E-02, &
                    0.66404E-02,   0.66485E-02,   0.66585E-02,   0.66658E-02,   0.66715E-02,   0.66777E-02,   0.66827E-02, &
                    0.66910E-02,   0.67007E-02,   0.67141E-02,   0.67212E-02,   0.67264E-02,   0.67354E-02,   0.67433E-02, &
                    0.67507E-02,   0.67600E-02,   0.67684E-02,   0.67716E-02,   0.67803E-02,   0.67894E-02,   0.67960E-02, &
                    0.68055E-02,   0.68127E-02,   0.68166E-02,   0.68236E-02,   0.68308E-02,   0.68416E-02,   0.68461E-02, &
                    0.68518E-02,   0.68612E-02,   0.68654E-02,   0.68738E-02,   0.68825E-02,   0.68903E-02,   0.68927E-02, &
                    0.68991E-02,   0.69069E-02,   0.69163E-02,   0.69258E-02,   0.69337E-02,   0.69385E-02,   0.69446E-02, &
                    0.69461E-02,   0.69608E-02,   0.69670E-02,   0.69718E-02,   0.69765E-02,   0.69833E-02,   0.69883E-02, &
                    0.69995E-02,   0.70072E-02,   0.70135E-02,   0.70179E-02,   0.70225E-02,   0.70317E-02,   0.70382E-02, &
                    0.70469E-02,   0.70519E-02,   0.70612E-02,   0.70658E-02,   0.70721E-02,   0.70796E-02,   0.70907E-02, &
                    0.70957E-02,   0.71024E-02,   0.71064E-02,   0.71166E-02,   0.71208E-02,   0.71305E-02,   0.71359E-02, &
                    0.71433E-02,   0.71484E-02,   0.71551E-02,   0.71669E-02,   0.71735E-02,   0.71781E-02,   0.71883E-02, &
                    0.71931E-02,   0.71984E-02,   0.72056E-02,   0.72111E-02,   0.72196E-02,   0.72237E-02,   0.72283E-02, &
                    0.72349E-02,   0.72412E-02,   0.72499E-02,   0.72573E-02,   0.72651E-02,   0.72732E-02,   0.72795E-02, &
                    0.72872E-02,   0.72959E-02,   0.73030E-02,   0.73106E-02,   0.73166E-02,   0.73225E-02,   0.73299E-02, &
                    0.73364E-02,   0.73413E-02,   0.73494E-02,   0.73537E-02,   0.73607E-02,   0.73715E-02,   0.73833E-02, &
                    0.73903E-02,   0.73933E-02,   0.74039E-02,   0.74122E-02,   0.74191E-02,   0.74245E-02,   0.74296E-02, &
                    0.74340E-02,   0.74427E-02,   0.74515E-02,   0.74563E-02,   0.74628E-02,   0.74691E-02,   0.74770E-02, &
                    0.74827E-02,   0.74867E-02,   0.74973E-02,   0.75018E-02,   0.75095E-02,   0.75141E-02,   0.75224E-02, &
                    0.75312E-02,   0.75396E-02,   0.75450E-02,   0.75501E-02,   0.75573E-02,   0.75674E-02,   0.75724E-02, &
                    0.75807E-02,   0.75898E-02,   0.75924E-02,   0.75962E-02,   0.76047E-02,   0.76141E-02,   0.76204E-02, &
                    0.76259E-02,   0.76329E-02,   0.76366E-02,   0.76479E-02,   0.76591E-02,   0.76656E-02,   0.76736E-02, &
                    0.76840E-02,   0.76890E-02,   0.76982E-02,   0.77057E-02,   0.77118E-02,   0.77191E-02,   0.77265E-02, &
                    0.77337E-02,   0.77394E-02,   0.77486E-02,   0.77577E-02,   0.77614E-02,   0.77669E-02,   0.77783E-02, &
                    0.77861E-02,   0.77925E-02,   0.78027E-02,   0.78071E-02,   0.78184E-02,   0.78220E-02,   0.78319E-02, &
                    0.78402E-02,   0.78456E-02,   0.78517E-02,   0.78603E-02,   0.78667E-02,   0.78756E-02,   0.78810E-02, &
                    0.78881E-02,   0.78936E-02,   0.79047E-02,   0.79086E-02,   0.79208E-02,   0.79244E-02,   0.79311E-02, &
                    0.79386E-02,   0.79496E-02,   0.79588E-02,   0.79641E-02,   0.79717E-02,   0.79743E-02,   0.79863E-02, &
                    0.79948E-02,   0.79971E-02,   0.80058E-02,   0.80120E-02,   0.80199E-02,   0.80250E-02,   0.80375E-02, &
                    0.80445E-02,   0.80502E-02,   0.80558E-02,   0.80644E-02,   0.80744E-02,   0.80840E-02,   0.80865E-02, &
                    0.80936E-02,   0.80995E-02,   0.81090E-02,   0.81153E-02,   0.81241E-02,   0.81311E-02,   0.81397E-02, &
                    0.81413E-02,   0.81517E-02,   0.81596E-02,   0.81665E-02,   0.81716E-02,   0.81790E-02,   0.81868E-02, &
                    0.81985E-02,   0.82050E-02,   0.82096E-02,   0.82143E-02,   0.82238E-02,   0.82286E-02,   0.82411E-02, &
                    0.82524E-02,   0.82568E-02,   0.82676E-02,   0.82792E-02,   0.82847E-02,   0.82923E-02,   0.83021E-02, &
                    0.83091E-02,   0.83148E-02,   0.83247E-02,   0.83363E-02,   0.83375E-02,   0.83476E-02,   0.83526E-02, &
                    0.83611E-02,   0.83720E-02,   0.83769E-02,   0.83806E-02,   0.83833E-02,   0.83941E-02,   0.83999E-02, &
                    0.84077E-02,   0.84132E-02,   0.84205E-02,   0.84267E-02,   0.84319E-02,   0.84373E-02,   0.84414E-02, &
                    0.84469E-02,   0.84533E-02,   0.84641E-02,   0.84690E-02,   0.84779E-02,   0.84826E-02,   0.84876E-02, &
                    0.84968E-02,   0.85076E-02,   0.85136E-02,   0.85172E-02,   0.85206E-02,   0.85348E-02,   0.85459E-02, &
                    0.85543E-02,   0.85670E-02,   0.85765E-02,   0.85844E-02,   0.85941E-02,   0.86020E-02,   0.86102E-02, &
                    0.86193E-02,   0.86279E-02,   0.86355E-02,   0.86482E-02,   0.86636E-02,   0.86678E-02,   0.86770E-02, &
                    0.86791E-02,   0.86921E-02,   0.86974E-02,   0.87057E-02,   0.87128E-02,   0.87209E-02,   0.87316E-02, &
                    0.87384E-02,   0.87417E-02,   0.87458E-02,   0.87614E-02,   0.87652E-02,   0.87693E-02,   0.87709E-02, &
                    0.87860E-02,   0.87907E-02,   0.87960E-02,   0.88037E-02,   0.88114E-02,   0.88134E-02,   0.88247E-02, &
                    0.88279E-02,   0.88358E-02,   0.88440E-02,   0.88531E-02,   0.88601E-02,   0.88744E-02,   0.88830E-02, &
                    0.88882E-02,   0.88891E-02,   0.89023E-02,   0.89128E-02,   0.89258E-02,   0.89304E-02,   0.89358E-02, &
                    0.89417E-02,   0.89534E-02,   0.89634E-02,   0.89689E-02,   0.89742E-02,   0.89887E-02,   0.89950E-02, &
                    0.90005E-02,   0.90085E-02,   0.90170E-02,   0.90212E-02,   0.90266E-02,   0.90346E-02,   0.90402E-02, &
                    0.90552E-02,   0.90595E-02,   0.90664E-02,   0.90755E-02,   0.90857E-02,   0.90986E-02,   0.91044E-02, &
                    0.91125E-02,   0.91203E-02,   0.91320E-02,   0.91394E-02,   0.91475E-02,   0.91602E-02,   0.91663E-02, &
                    0.91714E-02,   0.91733E-02,   0.91816E-02,   0.91923E-02,   0.92004E-02,   0.92065E-02,   0.92155E-02, &
                    0.92302E-02,   0.92370E-02,   0.92406E-02,   0.92527E-02,   0.92573E-02,   0.92663E-02,   0.92755E-02, &
                    0.92823E-02,   0.92851E-02,   0.93070E-02,   0.93127E-02,   0.93150E-02,   0.93242E-02,   0.93362E-02, &
                    0.93419E-02,   0.93477E-02,   0.93490E-02,   0.93649E-02,   0.93792E-02,   0.93914E-02,   0.93982E-02, &
                    0.94045E-02,   0.94051E-02,   0.94218E-02,   0.94191E-02,   0.94183E-02,   0.94364E-02,   0.94451E-02, &
                    0.94421E-02,   0.94529E-02,   0.94635E-02,   0.94661E-02,   0.94770E-02,   0.94782E-02,   0.94802E-02, &
                    0.94851E-02,   0.94939E-02,   0.94915E-02,   0.95052E-02,   0.95093E-02,   0.95141E-02,   0.95286E-02, &
                    0.95334E-02,   0.95420E-02,   0.95464E-02,   0.95609E-02,   0.95657E-02,   0.95743E-02,   0.95839E-02, &
                    0.95970E-02,   0.96083E-02,   0.96090E-02,   0.96199E-02,   0.96380E-02,   0.96380E-02,   0.96383E-02, &
                    0.96569E-02,   0.96673E-02,   0.96706E-02,   0.96865E-02,   0.96968E-02,   0.96972E-02,   0.97040E-02, &
                    0.97115E-02,   0.97163E-02,   0.97207E-02,   0.97287E-02,   0.97318E-02,   0.97448E-02,   0.97534E-02, &
                    0.97557E-02,   0.97604E-02,   0.97762E-02,   0.97807E-02,   0.97885E-02,   0.97921E-02,   0.98036E-02, &
                    0.98133E-02,   0.98159E-02,   0.98279E-02,   0.98348E-02,   0.98328E-02,   0.98405E-02,   0.98515E-02, &
                    0.98516E-02,   0.98508E-02,   0.98708E-02,   0.98826E-02,   0.98840E-02,   0.98850E-02,   0.98878E-02, &
                    0.98753E-02,   0.98741E-02,   0.98711E-02,   0.98741E-02,   0.98874E-02,   0.98922E-02,   0.98944E-02, &
                    0.99041E-02,   0.99150E-02,   0.99195E-02,   0.99251E-02,   0.99405E-02,   0.99459E-02,   0.99549E-02, &
                    0.99527E-02,   0.99693E-02,   0.99762E-02,   0.99712E-02,   0.99842E-02,   0.99966E-02,   0.99976E-02, &
                    0.10001E-01,   0.10019E-01,   0.10028E-01,   0.10047E-01,   0.10063E-01,   0.10074E-01,   0.10085E-01, &
                    0.10090E-01,   0.10097E-01,   0.10093E-01,   0.10088E-01,   0.10093E-01,   0.10097E-01,   0.10101E-01, &
                    0.10111E-01,   0.10119E-01,   0.10129E-01,   0.10136E-01,   0.10144E-01,   0.10156E-01,   0.10164E-01, &
                    0.10173E-01,   0.10180E-01,   0.10191E-01,   0.10203E-01,   0.10217E-01,   0.10224E-01,   0.10224E-01, &
                    0.10234E-01,   0.10243E-01,   0.10251E-01,   0.10268E-01,   0.10274E-01,   0.10286E-01,   0.10303E-01, &
                    0.10318E-01,   0.10330E-01,   0.10335E-01,   0.10338E-01,   0.10331E-01,   0.10328E-01,   0.10329E-01 /)
        solint(2101:2800)=(/                                                                                               &
                    0.10345E-01,   0.10358E-01,   0.10376E-01,   0.10384E-01,   0.10392E-01,   0.10406E-01,   0.10412E-01, &
                    0.10412E-01,   0.10416E-01,   0.10431E-01,   0.10452E-01,   0.10459E-01,   0.10473E-01,   0.10480E-01, &
                    0.10491E-01,   0.10502E-01,   0.10516E-01,   0.10525E-01,   0.10532E-01,   0.10544E-01,   0.10558E-01, &
                    0.10573E-01,   0.10592E-01,   0.10605E-01,   0.10607E-01,   0.10609E-01,   0.10605E-01,   0.10607E-01, &
                    0.10610E-01,   0.10627E-01,   0.10636E-01,   0.10650E-01,   0.10666E-01,   0.10676E-01,   0.10687E-01, &
                    0.10690E-01,   0.10704E-01,   0.10716E-01,   0.10729E-01,   0.10739E-01,   0.10753E-01,   0.10758E-01, &
                    0.10772E-01,   0.10780E-01,   0.10790E-01,   0.10805E-01,   0.10821E-01,   0.10833E-01,   0.10848E-01, &
                    0.10859E-01,   0.10872E-01,   0.10888E-01,   0.10902E-01,   0.10907E-01,   0.10913E-01,   0.10908E-01, &
                    0.10896E-01,   0.10898E-01,   0.10902E-01,   0.10913E-01,   0.10925E-01,   0.10939E-01,   0.10955E-01, &
                    0.10975E-01,   0.10986E-01,   0.11002E-01,   0.11012E-01,   0.11031E-01,   0.11047E-01,   0.11065E-01, &
                    0.11081E-01,   0.11092E-01,   0.11103E-01,   0.11117E-01,   0.11131E-01,   0.11150E-01,   0.11168E-01, &
                    0.11183E-01,   0.11198E-01,   0.11209E-01,   0.11236E-01,   0.11261E-01,   0.11281E-01,   0.11291E-01, &
                    0.11292E-01,   0.11287E-01,   0.11285E-01,   0.11286E-01,   0.11305E-01,   0.11320E-01,   0.11336E-01, &
                    0.11353E-01,   0.11373E-01,   0.11391E-01,   0.11406E-01,   0.11422E-01,   0.11445E-01,   0.11449E-01, &
                    0.11474E-01,   0.11495E-01,   0.11512E-01,   0.11531E-01,   0.11549E-01,   0.11574E-01,   0.11596E-01, &
                    0.11613E-01,   0.11629E-01,   0.11642E-01,   0.11663E-01,   0.11687E-01,   0.11721E-01,   0.11742E-01, &
                    0.11757E-01,   0.11762E-01,   0.11766E-01,   0.11752E-01,   0.11746E-01,   0.11754E-01,   0.11775E-01, &
                    0.11796E-01,   0.11816E-01,   0.11837E-01,   0.11848E-01,   0.11868E-01,   0.11886E-01,   0.11904E-01, &
                    0.11925E-01,   0.11947E-01,   0.11958E-01,   0.11979E-01,   0.12000E-01,   0.12022E-01,   0.12036E-01, &
                    0.12057E-01,   0.12081E-01,   0.12093E-01,   0.12116E-01,   0.12141E-01,   0.12169E-01,   0.12203E-01, &
                    0.12218E-01,   0.12262E-01,   0.12279E-01,   0.12286E-01,   0.12289E-01,   0.12298E-01,   0.12303E-01, &
                    0.12315E-01,   0.12328E-01,   0.12335E-01,   0.12352E-01,   0.12352E-01,   0.12364E-01,   0.12376E-01, &
                    0.12388E-01,   0.12399E-01,   0.12407E-01,   0.12418E-01,   0.12430E-01,   0.12440E-01,   0.12452E-01, &
                    0.12459E-01,   0.12465E-01,   0.12473E-01,   0.12479E-01,   0.12483E-01,   0.12479E-01,   0.12485E-01, &
                    0.12499E-01,   0.12505E-01,   0.12519E-01,   0.12520E-01,   0.12535E-01,   0.12543E-01,   0.12548E-01, &
                    0.12554E-01,   0.12566E-01,   0.12575E-01,   0.12584E-01,   0.12595E-01,   0.12596E-01,   0.12607E-01, &
                    0.12616E-01,   0.12626E-01,   0.12626E-01,   0.12637E-01,   0.12648E-01,   0.12661E-01,   0.12663E-01, &
                    0.12681E-01,   0.12702E-01,   0.12723E-01,   0.12733E-01,   0.12749E-01,   0.12766E-01,   0.12777E-01, &
                    0.12782E-01,   0.12793E-01,   0.12800E-01,   0.12817E-01,   0.12816E-01,   0.12826E-01,   0.12841E-01, &
                    0.12851E-01,   0.12864E-01,   0.12866E-01,   0.12884E-01,   0.12885E-01,   0.12906E-01,   0.12911E-01, &
                    0.12932E-01,   0.12934E-01,   0.12954E-01,   0.12957E-01,   0.12970E-01,   0.12978E-01,   0.12991E-01, &
                    0.13003E-01,   0.13004E-01,   0.13016E-01,   0.13026E-01,   0.13031E-01,   0.13043E-01,   0.13066E-01, &
                    0.13077E-01,   0.13086E-01,   0.13103E-01,   0.13113E-01,   0.13119E-01,   0.13131E-01,   0.13138E-01, &
                    0.13144E-01,   0.13155E-01,   0.13165E-01,   0.13173E-01,   0.13173E-01,   0.13186E-01,   0.13198E-01, &
                    0.13205E-01,   0.13217E-01,   0.13226E-01,   0.13234E-01,   0.13248E-01,   0.13256E-01,   0.13262E-01, &
                    0.13277E-01,   0.13286E-01,   0.13290E-01,   0.13295E-01,   0.13300E-01,   0.13307E-01,   0.13319E-01, &
                    0.13325E-01,   0.13331E-01,   0.13333E-01,   0.13314E-01,   0.13305E-01,   0.13300E-01,   0.13310E-01, &
                    0.13325E-01,   0.13336E-01,   0.13334E-01,   0.13345E-01,   0.13356E-01,   0.13372E-01,   0.13373E-01, &
                    0.13386E-01,   0.13387E-01,   0.13401E-01,   0.13418E-01,   0.13426E-01,   0.13438E-01,   0.13459E-01, &
                    0.13464E-01,   0.13481E-01,   0.13488E-01,   0.13492E-01,   0.13497E-01,   0.13500E-01,   0.13517E-01, &
                    0.13559E-01,   0.13577E-01,   0.13607E-01,   0.13626E-01,   0.13625E-01,   0.13652E-01,   0.13659E-01, &
                    0.13667E-01,   0.13675E-01,   0.13679E-01,   0.13688E-01,   0.13692E-01,   0.13693E-01,   0.13705E-01, &
                    0.13705E-01,   0.13719E-01,   0.13729E-01,   0.13728E-01,   0.13742E-01,   0.13750E-01,   0.13763E-01, &
                    0.13781E-01,   0.13794E-01,   0.13817E-01,   0.13822E-01,   0.13839E-01,   0.13850E-01,   0.13859E-01, &
                    0.13869E-01,   0.13876E-01,   0.13885E-01,   0.13902E-01,   0.13907E-01,   0.13928E-01,   0.13931E-01, &
                    0.13948E-01,   0.13958E-01,   0.13971E-01,   0.13987E-01,   0.13990E-01,   0.14001E-01,   0.14016E-01, &
                    0.14028E-01,   0.14044E-01,   0.14057E-01,   0.14076E-01,   0.14074E-01,   0.14093E-01,   0.14097E-01, &
                    0.14113E-01,   0.14121E-01,   0.14133E-01,   0.14145E-01,   0.14160E-01,   0.14172E-01,   0.14176E-01, &
                    0.14190E-01,   0.14199E-01,   0.14214E-01,   0.14229E-01,   0.14246E-01,   0.14262E-01,   0.14268E-01, &
                    0.14291E-01,   0.14296E-01,   0.14317E-01,   0.14322E-01,   0.14335E-01,   0.14349E-01,   0.14355E-01, &
                    0.14363E-01,   0.14373E-01,   0.14381E-01,   0.14393E-01,   0.14398E-01,   0.14413E-01,   0.14426E-01, &
                    0.14429E-01,   0.14440E-01,   0.14449E-01,   0.14458E-01,   0.14464E-01,   0.14467E-01,   0.14474E-01, &
                    0.14487E-01,   0.14497E-01,   0.14510E-01,   0.14520E-01,   0.14527E-01,   0.14542E-01,   0.14546E-01, &
                    0.14553E-01,   0.14569E-01,   0.14571E-01,   0.14580E-01,   0.14594E-01,   0.14596E-01,   0.14609E-01, &
                    0.14620E-01,   0.14634E-01,   0.14638E-01,   0.14657E-01,   0.14662E-01,   0.14680E-01,   0.14684E-01, &
                    0.14701E-01,   0.14714E-01,   0.14736E-01,   0.14747E-01,   0.14761E-01,   0.14765E-01,   0.14774E-01, &
                    0.14787E-01,   0.14796E-01,   0.14805E-01,   0.14808E-01,   0.14821E-01,   0.14831E-01,   0.14847E-01, &
                    0.14858E-01,   0.14877E-01,   0.14890E-01,   0.14902E-01,   0.14912E-01,   0.14928E-01,   0.14938E-01, &
                    0.14947E-01,   0.14953E-01,   0.14955E-01,   0.14964E-01,   0.14976E-01,   0.14984E-01,   0.14993E-01, &
                    0.14994E-01,   0.15001E-01,   0.15012E-01,   0.15022E-01,   0.15030E-01,   0.15039E-01,   0.15052E-01, &
                    0.15061E-01,   0.15070E-01,   0.15085E-01,   0.15094E-01,   0.15090E-01,   0.15102E-01,   0.15104E-01, &
                    0.15111E-01,   0.15120E-01,   0.15123E-01,   0.15142E-01,   0.15147E-01,   0.15165E-01,   0.15175E-01, &
                    0.15180E-01,   0.15189E-01,   0.15190E-01,   0.15204E-01,   0.15208E-01,   0.15225E-01,   0.15235E-01, &
                    0.15248E-01,   0.15253E-01,   0.15258E-01,   0.15274E-01,   0.15284E-01,   0.15293E-01,   0.15308E-01, &
                    0.15316E-01,   0.15325E-01,   0.15330E-01,   0.15343E-01,   0.15338E-01,   0.15356E-01,   0.15359E-01, &
                    0.15368E-01,   0.15366E-01,   0.15365E-01,   0.15355E-01,   0.15350E-01,   0.15342E-01,   0.15337E-01, &
                    0.15328E-01,   0.15328E-01,   0.15328E-01,   0.15343E-01,   0.15342E-01,   0.15357E-01,   0.15368E-01, &
                    0.15375E-01,   0.15402E-01,   0.15409E-01,   0.15416E-01,   0.15436E-01,   0.15450E-01,   0.15471E-01, &
                    0.15481E-01,   0.15492E-01,   0.15514E-01,   0.15525E-01,   0.15540E-01,   0.15548E-01,   0.15561E-01, &
                    0.15598E-01,   0.15615E-01,   0.15648E-01,   0.15675E-01,   0.15695E-01,   0.15725E-01,   0.15735E-01, &
                    0.15728E-01,   0.15726E-01,   0.15720E-01,   0.15725E-01,   0.15735E-01,   0.15735E-01,   0.15750E-01, &
                    0.15760E-01,   0.15777E-01,   0.15784E-01,   0.15797E-01,   0.15808E-01,   0.15815E-01,   0.15814E-01, &
                    0.15827E-01,   0.15831E-01,   0.15845E-01,   0.15858E-01,   0.15864E-01,   0.15899E-01,   0.15910E-01, &
                    0.15919E-01,   0.15940E-01,   0.15955E-01,   0.15961E-01,   0.16009E-01,   0.16026E-01,   0.16052E-01, &
                    0.16083E-01,   0.16090E-01,   0.16103E-01,   0.16109E-01,   0.16121E-01,   0.16131E-01,   0.16131E-01, &
                    0.16136E-01,   0.16148E-01,   0.16167E-01,   0.16181E-01,   0.16210E-01,   0.16221E-01,   0.16238E-01, &
                    0.16256E-01,   0.16268E-01,   0.16274E-01,   0.16286E-01,   0.16300E-01,   0.16314E-01,   0.16332E-01, &
                    0.16340E-01,   0.16353E-01,   0.16365E-01,   0.16371E-01,   0.16382E-01,   0.16400E-01,   0.16416E-01, &
                    0.16427E-01,   0.16449E-01,   0.16459E-01,   0.16472E-01,   0.16505E-01,   0.16515E-01,   0.16527E-01, &
                    0.16543E-01,   0.16551E-01,   0.16565E-01,   0.16585E-01,   0.16595E-01,   0.16608E-01,   0.16615E-01, &
                    0.16636E-01,   0.16638E-01,   0.16651E-01,   0.16665E-01,   0.16676E-01,   0.16682E-01,   0.16700E-01, &
                    0.16700E-01,   0.16717E-01,   0.16728E-01,   0.16734E-01,   0.16741E-01,   0.16753E-01,   0.16767E-01, &
                    0.16779E-01,   0.16781E-01,   0.16797E-01,   0.16802E-01,   0.16812E-01,   0.16817E-01,   0.16832E-01, &
                    0.16831E-01,   0.16842E-01,   0.16856E-01,   0.16864E-01,   0.16876E-01,   0.16890E-01,   0.16895E-01, &
                    0.16904E-01,   0.16907E-01,   0.16922E-01,   0.16931E-01,   0.16946E-01,   0.16966E-01,   0.16974E-01, &
                    0.16987E-01,   0.16992E-01,   0.17000E-01,   0.17005E-01,   0.17015E-01,   0.17022E-01,   0.17035E-01, &
                    0.17039E-01,   0.17049E-01,   0.17074E-01,   0.17085E-01,   0.17101E-01,   0.17105E-01,   0.17127E-01, &
                    0.17135E-01,   0.17143E-01,   0.17159E-01,   0.17173E-01,   0.17182E-01,   0.17190E-01,   0.17209E-01, &
                    0.17215E-01,   0.17220E-01,   0.17226E-01,   0.17242E-01,   0.17256E-01,   0.17261E-01,   0.17279E-01, &
                    0.17300E-01,   0.17317E-01,   0.17323E-01,   0.17338E-01,   0.17353E-01,   0.17361E-01,   0.17374E-01, &
                    0.17376E-01,   0.17396E-01,   0.17404E-01,   0.17412E-01,   0.17419E-01,   0.17435E-01,   0.17448E-01, &
                    0.17452E-01,   0.17469E-01,   0.17493E-01,   0.17504E-01,   0.17514E-01,   0.17531E-01,   0.17548E-01, &
                    0.17550E-01,   0.17564E-01,   0.17573E-01,   0.17592E-01,   0.17599E-01,   0.17611E-01,   0.17614E-01, &
                    0.17628E-01,   0.17645E-01,   0.17653E-01,   0.17662E-01,   0.17671E-01,   0.17687E-01,   0.17701E-01, &
                    0.17704E-01,   0.17719E-01,   0.17723E-01,   0.17746E-01,   0.17755E-01,   0.17758E-01,   0.17770E-01, &
                    0.17783E-01,   0.17789E-01,   0.17795E-01,   0.17803E-01,   0.17817E-01,   0.17835E-01,   0.17842E-01, &
                    0.17844E-01,   0.17857E-01,   0.17878E-01,   0.17885E-01,   0.17901E-01,   0.17911E-01,   0.17920E-01, &
                    0.17937E-01,   0.17946E-01,   0.17955E-01,   0.17965E-01,   0.17979E-01,   0.17991E-01,   0.18006E-01, &
                    0.18010E-01,   0.18025E-01,   0.18029E-01,   0.18051E-01,   0.18066E-01,   0.18072E-01,   0.18089E-01, &
                    0.18099E-01,   0.18120E-01,   0.18130E-01,   0.18135E-01,   0.18151E-01,   0.18161E-01,   0.18180E-01 /)
        solint(2801:3500)=(/                                                                                               &
                    0.18189E-01,   0.18193E-01,   0.18207E-01,   0.18215E-01,   0.18223E-01,   0.18230E-01,   0.18236E-01, &
                    0.18240E-01,   0.18256E-01,   0.18265E-01,   0.18284E-01,   0.18289E-01,   0.18305E-01,   0.18313E-01, &
                    0.18313E-01,   0.18318E-01,   0.18319E-01,   0.18333E-01,   0.18338E-01,   0.18352E-01,   0.18355E-01, &
                    0.18366E-01,   0.18373E-01,   0.18388E-01,   0.18393E-01,   0.18404E-01,   0.18414E-01,   0.18415E-01, &
                    0.18415E-01,   0.18419E-01,   0.18425E-01,   0.18423E-01,   0.18429E-01,   0.18442E-01,   0.18447E-01, &
                    0.18458E-01,   0.18474E-01,   0.18480E-01,   0.18491E-01,   0.18505E-01,   0.18547E-01,   0.18557E-01, &
                    0.18574E-01,   0.18582E-01,   0.18584E-01,   0.18593E-01,   0.18621E-01,   0.18633E-01,   0.18643E-01, &
                    0.18653E-01,   0.18663E-01,   0.18668E-01,   0.18689E-01,   0.18726E-01,   0.18740E-01,   0.18757E-01, &
                    0.18777E-01,   0.18796E-01,   0.18807E-01,   0.18813E-01,   0.18830E-01,   0.18835E-01,   0.18845E-01, &
                    0.18856E-01,   0.18867E-01,   0.18877E-01,   0.18887E-01,   0.18911E-01,   0.18919E-01,   0.18930E-01, &
                    0.18940E-01,   0.18950E-01,   0.18960E-01,   0.18965E-01,   0.18992E-01,   0.18998E-01,   0.19006E-01, &
                    0.19025E-01,   0.19039E-01,   0.19050E-01,   0.19065E-01,   0.19076E-01,   0.19100E-01,   0.19115E-01, &
                    0.19124E-01,   0.19136E-01,   0.19147E-01,   0.19158E-01,   0.19168E-01,   0.19172E-01,   0.19194E-01, &
                    0.19206E-01,   0.19218E-01,   0.19230E-01,   0.19238E-01,   0.19250E-01,   0.19261E-01,   0.19274E-01, &
                    0.19282E-01,   0.19258E-01,   0.19250E-01,   0.19253E-01,   0.19255E-01,   0.19269E-01,   0.19275E-01, &
                    0.19283E-01,   0.19303E-01,   0.19312E-01,   0.19331E-01,   0.19343E-01,   0.19351E-01,   0.19362E-01, &
                    0.19373E-01,   0.19384E-01,   0.19391E-01,   0.19407E-01,   0.19409E-01,   0.19409E-01,   0.19414E-01, &
                    0.19404E-01,   0.19389E-01,   0.19370E-01,   0.19365E-01,   0.19374E-01,   0.19401E-01,   0.19428E-01, &
                    0.19448E-01,   0.19460E-01,   0.19474E-01,   0.19493E-01,   0.19504E-01,   0.19510E-01,   0.19528E-01, &
                    0.19531E-01,   0.19547E-01,   0.19550E-01,   0.19559E-01,   0.19570E-01,   0.19584E-01,   0.19594E-01, &
                    0.19614E-01,   0.19631E-01,   0.19652E-01,   0.19676E-01,   0.19704E-01,   0.19746E-01,   0.19786E-01, &
                    0.19813E-01,   0.19838E-01,   0.19858E-01,   0.19877E-01,   0.19891E-01,   0.19904E-01,   0.19930E-01, &
                    0.19941E-01,   0.19950E-01,   0.19966E-01,   0.19983E-01,   0.19998E-01,   0.20007E-01,   0.20017E-01, &
                    0.20038E-01,   0.20047E-01,   0.20056E-01,   0.20075E-01,   0.20087E-01,   0.20099E-01,   0.20119E-01, &
                    0.20140E-01,   0.20145E-01,   0.20153E-01,   0.20154E-01,   0.20161E-01,   0.20172E-01,   0.20181E-01, &
                    0.20197E-01,   0.20211E-01,   0.20221E-01,   0.20236E-01,   0.20243E-01,   0.20254E-01,   0.20261E-01, &
                    0.20279E-01,   0.20282E-01,   0.20300E-01,   0.20321E-01,   0.20336E-01,   0.20349E-01,   0.20360E-01, &
                    0.20360E-01,   0.20376E-01,   0.20384E-01,   0.20397E-01,   0.20406E-01,   0.20423E-01,   0.20437E-01, &
                    0.20458E-01,   0.20477E-01,   0.20484E-01,   0.20487E-01,   0.20504E-01,   0.20508E-01,   0.20527E-01, &
                    0.20536E-01,   0.20547E-01,   0.20567E-01,   0.20585E-01,   0.20596E-01,   0.20605E-01,   0.20633E-01, &
                    0.20637E-01,   0.20642E-01,   0.20645E-01,   0.20661E-01,   0.20678E-01,   0.20689E-01,   0.20697E-01, &
                    0.20711E-01,   0.20720E-01,   0.20734E-01,   0.20749E-01,   0.20758E-01,   0.20770E-01,   0.20790E-01, &
                    0.20809E-01,   0.20819E-01,   0.20828E-01,   0.20843E-01,   0.20857E-01,   0.20865E-01,   0.20871E-01, &
                    0.20880E-01,   0.20889E-01,   0.20891E-01,   0.20910E-01,   0.20922E-01,   0.20938E-01,   0.20962E-01, &
                    0.20976E-01,   0.20991E-01,   0.20999E-01,   0.21009E-01,   0.21010E-01,   0.21018E-01,   0.21024E-01, &
                    0.21030E-01,   0.21043E-01,   0.21055E-01,   0.21061E-01,   0.21069E-01,   0.21085E-01,   0.21094E-01, &
                    0.21102E-01,   0.21122E-01,   0.21127E-01,   0.21144E-01,   0.21159E-01,   0.21180E-01,   0.21193E-01, &
                    0.21207E-01,   0.21203E-01,   0.21190E-01,   0.21195E-01,   0.21210E-01,   0.21214E-01,   0.21228E-01, &
                    0.21246E-01,   0.21266E-01,   0.21284E-01,   0.21302E-01,   0.21315E-01,   0.21327E-01,   0.21343E-01, &
                    0.21355E-01,   0.21371E-01,   0.21381E-01,   0.21390E-01,   0.21406E-01,   0.21417E-01,   0.21433E-01, &
                    0.21443E-01,   0.21451E-01,   0.21464E-01,   0.21472E-01,   0.21483E-01,   0.21513E-01,   0.21542E-01, &
                    0.21558E-01,   0.21573E-01,   0.21584E-01,   0.21597E-01,   0.21612E-01,   0.21619E-01,   0.21627E-01, &
                    0.21642E-01,   0.21650E-01,   0.21654E-01,   0.21639E-01,   0.21633E-01,   0.21632E-01,   0.21638E-01, &
                    0.21650E-01,   0.21662E-01,   0.21673E-01,   0.21689E-01,   0.21704E-01,   0.21710E-01,   0.21728E-01, &
                    0.21741E-01,   0.21759E-01,   0.21768E-01,   0.21786E-01,   0.21804E-01,   0.21810E-01,   0.21819E-01, &
                    0.21821E-01,   0.21835E-01,   0.21843E-01,   0.21861E-01,   0.21879E-01,   0.21896E-01,   0.21931E-01, &
                    0.21964E-01,   0.21994E-01,   0.22018E-01,   0.22023E-01,   0.22034E-01,   0.22045E-01,   0.22048E-01, &
                    0.22062E-01,   0.22071E-01,   0.22079E-01,   0.22092E-01,   0.22096E-01,   0.22112E-01,   0.22123E-01, &
                    0.22140E-01,   0.22146E-01,   0.22169E-01,   0.22191E-01,   0.22208E-01,   0.22225E-01,   0.22232E-01, &
                    0.22248E-01,   0.22255E-01,   0.22262E-01,   0.22278E-01,   0.22289E-01,   0.22288E-01,   0.22300E-01, &
                    0.22320E-01,   0.22333E-01,   0.22348E-01,   0.22369E-01,   0.22381E-01,   0.22394E-01,   0.22408E-01, &
                    0.22425E-01,   0.22427E-01,   0.22441E-01,   0.22453E-01,   0.22466E-01,   0.22479E-01,   0.22492E-01, &
                    0.22509E-01,   0.22520E-01,   0.22525E-01,   0.22542E-01,   0.22555E-01,   0.22567E-01,   0.22580E-01, &
                    0.22583E-01,   0.22595E-01,   0.22577E-01,   0.22553E-01,   0.22534E-01,   0.22528E-01,   0.22531E-01, &
                    0.22541E-01,   0.22549E-01,   0.22557E-01,   0.22561E-01,   0.22579E-01,   0.22593E-01,   0.22600E-01, &
                    0.22616E-01,   0.22625E-01,   0.22636E-01,   0.22649E-01,   0.22667E-01,   0.22676E-01,   0.22695E-01, &
                    0.22718E-01,   0.22728E-01,   0.22753E-01,   0.22776E-01,   0.22821E-01,   0.22869E-01,   0.22919E-01, &
                    0.22949E-01,   0.22975E-01,   0.22997E-01,   0.23014E-01,   0.23043E-01,   0.23064E-01,   0.23076E-01, &
                    0.23089E-01,   0.23103E-01,   0.23119E-01,   0.23136E-01,   0.23152E-01,   0.23167E-01,   0.23186E-01, &
                    0.23187E-01,   0.23193E-01,   0.23205E-01,   0.23218E-01,   0.23225E-01,   0.23234E-01,   0.23249E-01, &
                    0.23264E-01,   0.23280E-01,   0.23294E-01,   0.23305E-01,   0.23316E-01,   0.23327E-01,   0.23338E-01, &
                    0.23349E-01,   0.23372E-01,   0.23388E-01,   0.23402E-01,   0.23412E-01,   0.23430E-01,   0.23439E-01, &
                    0.23456E-01,   0.23466E-01,   0.23488E-01,   0.23514E-01,   0.23528E-01,   0.23551E-01,   0.23574E-01, &
                    0.23590E-01,   0.23604E-01,   0.23618E-01,   0.23627E-01,   0.23648E-01,   0.23660E-01,   0.23682E-01, &
                    0.23695E-01,   0.23712E-01,   0.23727E-01,   0.23741E-01,   0.23754E-01,   0.23766E-01,   0.23779E-01, &
                    0.23788E-01,   0.23805E-01,   0.23829E-01,   0.23840E-01,   0.23852E-01,   0.23859E-01,   0.23875E-01, &
                    0.23886E-01,   0.23893E-01,   0.23909E-01,   0.23926E-01,   0.23935E-01,   0.23941E-01,   0.23947E-01, &
                    0.23964E-01,   0.23966E-01,   0.23975E-01,   0.23989E-01,   0.24003E-01,   0.24014E-01,   0.24030E-01, &
                    0.24050E-01,   0.24060E-01,   0.24069E-01,   0.24088E-01,   0.24099E-01,   0.24116E-01,   0.24133E-01, &
                    0.24153E-01,   0.24163E-01,   0.24171E-01,   0.24181E-01,   0.24196E-01,   0.24213E-01,   0.24222E-01, &
                    0.24231E-01,   0.24243E-01,   0.24262E-01,   0.24275E-01,   0.24292E-01,   0.24308E-01,   0.24321E-01, &
                    0.24345E-01,   0.24353E-01,   0.24366E-01,   0.24376E-01,   0.24392E-01,   0.24410E-01,   0.24427E-01, &
                    0.24438E-01,   0.24455E-01,   0.24479E-01,   0.24501E-01,   0.24507E-01,   0.24524E-01,   0.24546E-01, &
                    0.24576E-01,   0.24595E-01,   0.24611E-01,   0.24634E-01,   0.24643E-01,   0.24650E-01,   0.24669E-01, &
                    0.24680E-01,   0.24698E-01,   0.24713E-01,   0.24732E-01,   0.24751E-01,   0.24763E-01,   0.24778E-01, &
                    0.24784E-01,   0.24797E-01,   0.24804E-01,   0.24814E-01,   0.24835E-01,   0.24855E-01,   0.24866E-01, &
                    0.24882E-01,   0.24897E-01,   0.24907E-01,   0.24921E-01,   0.24936E-01,   0.24947E-01,   0.24954E-01, &
                    0.24964E-01,   0.24952E-01,   0.24949E-01,   0.24931E-01,   0.24918E-01,   0.24915E-01,   0.24918E-01, &
                    0.24929E-01,   0.24933E-01,   0.24955E-01,   0.24968E-01,   0.24982E-01,   0.25000E-01,   0.25017E-01, &
                    0.25029E-01,   0.25045E-01,   0.25058E-01,   0.25079E-01,   0.25096E-01,   0.25119E-01,   0.25146E-01, &
                    0.25167E-01,   0.25209E-01,   0.25248E-01,   0.25284E-01,   0.25320E-01,   0.25346E-01,   0.25377E-01, &
                    0.25408E-01,   0.25430E-01,   0.25447E-01,   0.25465E-01,   0.25474E-01,   0.25489E-01,   0.25496E-01, &
                    0.25503E-01,   0.25514E-01,   0.25534E-01,   0.25549E-01,   0.25568E-01,   0.25583E-01,   0.25598E-01, &
                    0.25620E-01,   0.25642E-01,   0.25656E-01,   0.25685E-01,   0.25698E-01,   0.25717E-01,   0.25733E-01, &
                    0.25746E-01,   0.25759E-01,   0.25775E-01,   0.25790E-01,   0.25806E-01,   0.25830E-01,   0.25858E-01, &
                    0.25883E-01,   0.25897E-01,   0.25917E-01,   0.25928E-01,   0.25938E-01,   0.25949E-01,   0.25966E-01, &
                    0.25984E-01,   0.25992E-01,   0.26012E-01,   0.26032E-01,   0.26043E-01,   0.26065E-01,   0.26082E-01, &
                    0.26091E-01,   0.26112E-01,   0.26127E-01,   0.26135E-01,   0.26154E-01,   0.26173E-01,   0.26183E-01, &
                    0.26199E-01,   0.26214E-01,   0.26225E-01,   0.26249E-01,   0.26267E-01,   0.26283E-01,   0.26303E-01, &
                    0.26317E-01,   0.26332E-01,   0.26353E-01,   0.26369E-01,   0.26385E-01,   0.26400E-01,   0.26418E-01, &
                    0.26429E-01,   0.26451E-01,   0.26465E-01,   0.26463E-01,   0.26477E-01,   0.26494E-01,   0.26515E-01, &
                    0.26534E-01,   0.26537E-01,   0.26550E-01,   0.26556E-01,   0.26563E-01,   0.26570E-01,   0.26576E-01, &
                    0.26586E-01,   0.26599E-01,   0.26613E-01,   0.26628E-01,   0.26647E-01,   0.26657E-01,   0.26666E-01, &
                    0.26679E-01,   0.26693E-01,   0.26706E-01,   0.26715E-01,   0.26704E-01,   0.26693E-01,   0.26689E-01, &
                    0.26683E-01,   0.26682E-01,   0.26695E-01,   0.26709E-01,   0.26738E-01,   0.26751E-01,   0.26772E-01, &
                    0.26787E-01,   0.26801E-01,   0.26822E-01,   0.26840E-01,   0.26868E-01,   0.26900E-01,   0.26926E-01, &
                    0.26959E-01,   0.26984E-01,   0.27033E-01,   0.27076E-01,   0.27104E-01,   0.27134E-01,   0.27174E-01, &
                    0.27205E-01,   0.27231E-01,   0.27252E-01,   0.27269E-01,   0.27286E-01,   0.27306E-01,   0.27326E-01, &
                    0.27350E-01,   0.27361E-01,   0.27366E-01,   0.27367E-01,   0.27380E-01,   0.27389E-01,   0.27405E-01, &
                    0.27421E-01,   0.27438E-01,   0.27474E-01,   0.27497E-01,   0.27508E-01,   0.27519E-01,   0.27532E-01 /)
        solint(3501:4200)=(/                                                                                               &
                    0.27549E-01,   0.27569E-01,   0.27588E-01,   0.27598E-01,   0.27606E-01,   0.27617E-01,   0.27636E-01, &
                    0.27675E-01,   0.27690E-01,   0.27702E-01,   0.27713E-01,   0.27730E-01,   0.27753E-01,   0.27777E-01, &
                    0.27792E-01,   0.27808E-01,   0.27829E-01,   0.27860E-01,   0.27875E-01,   0.27886E-01,   0.27898E-01, &
                    0.27897E-01,   0.27926E-01,   0.27952E-01,   0.27961E-01,   0.27979E-01,   0.27996E-01,   0.28030E-01, &
                    0.28052E-01,   0.28065E-01,   0.28086E-01,   0.28096E-01,   0.28115E-01,   0.28139E-01,   0.28158E-01, &
                    0.28167E-01,   0.28187E-01,   0.28203E-01,   0.28227E-01,   0.28231E-01,   0.28236E-01,   0.28240E-01, &
                    0.28256E-01,   0.28271E-01,   0.28295E-01,   0.28310E-01,   0.28332E-01,   0.28354E-01,   0.28369E-01, &
                    0.28363E-01,   0.28346E-01,   0.28313E-01,   0.28306E-01,   0.28306E-01,   0.28309E-01,   0.28340E-01, &
                    0.28364E-01,   0.28411E-01,   0.28463E-01,   0.28508E-01,   0.28532E-01,   0.28509E-01,   0.28503E-01, &
                    0.28537E-01,   0.28568E-01,   0.28590E-01,   0.28620E-01,   0.28691E-01,   0.28753E-01,   0.28828E-01, &
                    0.28879E-01,   0.28929E-01,   0.28969E-01,   0.28985E-01,   0.28999E-01,   0.29008E-01,   0.28964E-01, &
                    0.28946E-01,   0.28999E-01,   0.29086E-01,   0.29128E-01,   0.29154E-01,   0.29180E-01,   0.29213E-01, &
                    0.29238E-01,   0.29266E-01,   0.29291E-01,   0.29299E-01,   0.29299E-01,   0.29248E-01,   0.29109E-01, &
                    0.28985E-01,   0.28942E-01,   0.28945E-01,   0.28995E-01,   0.29013E-01,   0.29021E-01,   0.29021E-01, &
                    0.29022E-01,   0.29029E-01,   0.29044E-01,   0.29055E-01,   0.29074E-01,   0.29092E-01,   0.29114E-01, &
                    0.29152E-01,   0.29219E-01,   0.29372E-01,   0.29576E-01,   0.29671E-01,   0.29750E-01,   0.29805E-01, &
                    0.29832E-01,   0.29852E-01,   0.29880E-01,   0.29903E-01,   0.29916E-01,   0.29921E-01,   0.29928E-01, &
                    0.29931E-01,   0.29924E-01,   0.29932E-01,   0.29939E-01,   0.29929E-01,   0.29934E-01,   0.29942E-01, &
                    0.29946E-01,   0.29947E-01,   0.29955E-01,   0.29963E-01,   0.29929E-01,   0.29851E-01,   0.29817E-01, &
                    0.29806E-01,   0.29796E-01,   0.29783E-01,   0.29801E-01,   0.29835E-01,   0.29860E-01,   0.29884E-01, &
                    0.29924E-01,   0.29953E-01,   0.29964E-01,   0.29979E-01,   0.30002E-01,   0.30028E-01,   0.30035E-01, &
                    0.30115E-01,   0.30203E-01,   0.30222E-01,   0.30251E-01,   0.30283E-01,   0.30300E-01,   0.30311E-01, &
                    0.30313E-01,   0.30322E-01,   0.30316E-01,   0.30308E-01,   0.30319E-01,   0.30313E-01,   0.30303E-01, &
                    0.30294E-01,   0.30321E-01,   0.30343E-01,   0.30332E-01,   0.30350E-01,   0.30386E-01,   0.30405E-01, &
                    0.30432E-01,   0.30447E-01,   0.30448E-01,   0.30465E-01,   0.30495E-01,   0.30536E-01,   0.30534E-01, &
                    0.30506E-01,   0.30527E-01,   0.30553E-01,   0.30563E-01,   0.30575E-01,   0.30611E-01,   0.30654E-01, &
                    0.30689E-01,   0.30723E-01,   0.30754E-01,   0.30793E-01,   0.30847E-01,   0.30887E-01,   0.30909E-01, &
                    0.30939E-01,   0.30987E-01,   0.31076E-01,   0.31138E-01,   0.31176E-01,   0.31218E-01,   0.31256E-01, &
                    0.31300E-01,   0.31328E-01,   0.31358E-01,   0.31357E-01,   0.31373E-01,   0.31387E-01,   0.31403E-01, &
                    0.31408E-01,   0.31421E-01,   0.31433E-01,   0.31431E-01,   0.31433E-01,   0.31434E-01,   0.31450E-01, &
                    0.31462E-01,   0.31493E-01,   0.31518E-01,   0.31539E-01,   0.31572E-01,   0.31627E-01,   0.31666E-01, &
                    0.31693E-01,   0.31725E-01,   0.31755E-01,   0.31780E-01,   0.31806E-01,   0.31846E-01,   0.31880E-01, &
                    0.31903E-01,   0.31942E-01,   0.31953E-01,   0.31967E-01,   0.31981E-01,   0.31999E-01,   0.31941E-01, &
                    0.31900E-01,   0.31904E-01,   0.31925E-01,   0.31901E-01,   0.31856E-01,   0.31839E-01,   0.31819E-01, &
                    0.31798E-01,   0.31815E-01,   0.31864E-01,   0.31891E-01,   0.31919E-01,   0.31931E-01,   0.31920E-01, &
                    0.31969E-01,   0.32061E-01,   0.32049E-01,   0.32039E-01,   0.32082E-01,   0.32182E-01,   0.32252E-01, &
                    0.32318E-01,   0.32399E-01,   0.32449E-01,   0.32456E-01,   0.32472E-01,   0.32498E-01,   0.32509E-01, &
                    0.32552E-01,   0.32523E-01,   0.32485E-01,   0.32521E-01,   0.32563E-01,   0.32572E-01,   0.32564E-01, &
                    0.32531E-01,   0.32525E-01,   0.32525E-01,   0.32509E-01,   0.32484E-01,   0.32451E-01,   0.32464E-01, &
                    0.32489E-01,   0.32534E-01,   0.32619E-01,   0.32708E-01,   0.32773E-01,   0.32825E-01,   0.32851E-01, &
                    0.32875E-01,   0.32949E-01,   0.32980E-01,   0.32997E-01,   0.33052E-01,   0.33132E-01,   0.33206E-01, &
                    0.33245E-01,   0.33256E-01,   0.33272E-01,   0.33279E-01,   0.33293E-01,   0.33295E-01,   0.33300E-01, &
                    0.33317E-01,   0.33364E-01,   0.33385E-01,   0.33427E-01,   0.33462E-01,   0.33469E-01,   0.33469E-01, &
                    0.33399E-01,   0.33318E-01,   0.33252E-01,   0.33257E-01,   0.33303E-01,   0.33345E-01,   0.33391E-01, &
                    0.33446E-01,   0.33473E-01,   0.33516E-01,   0.33535E-01,   0.33573E-01,   0.33610E-01,   0.33640E-01, &
                    0.33714E-01,   0.33861E-01,   0.33955E-01,   0.34030E-01,   0.34066E-01,   0.34060E-01,   0.34075E-01, &
                    0.34080E-01,   0.34099E-01,   0.34102E-01,   0.34102E-01,   0.34104E-01,   0.34131E-01,   0.34123E-01, &
                    0.34135E-01,   0.34114E-01,   0.34075E-01,   0.34038E-01,   0.33978E-01,   0.33900E-01,   0.33840E-01, &
                    0.33889E-01,   0.33931E-01,   0.33955E-01,   0.34015E-01,   0.34052E-01,   0.34103E-01,   0.34126E-01, &
                    0.34125E-01,   0.34074E-01,   0.34155E-01,   0.34242E-01,   0.34355E-01,   0.34449E-01,   0.34547E-01, &
                    0.34572E-01,   0.34589E-01,   0.34614E-01,   0.34619E-01,   0.34623E-01,   0.34640E-01,   0.34659E-01, &
                    0.34718E-01,   0.34806E-01,   0.34840E-01,   0.34814E-01,   0.34807E-01,   0.34812E-01,   0.34801E-01, &
                    0.34768E-01,   0.34691E-01,   0.34652E-01,   0.34676E-01,   0.34720E-01,   0.34750E-01,   0.34795E-01, &
                    0.34827E-01,   0.34867E-01,   0.34911E-01,   0.34977E-01,   0.35062E-01,   0.35122E-01,   0.35198E-01, &
                    0.35290E-01,   0.35379E-01,   0.35414E-01,   0.35430E-01,   0.35447E-01,   0.35416E-01,   0.35433E-01, &
                    0.35450E-01,   0.35447E-01,   0.35464E-01,   0.35472E-01,   0.35478E-01,   0.35487E-01,   0.35482E-01, &
                    0.35455E-01,   0.35481E-01,   0.35474E-01,   0.35394E-01,   0.35419E-01,   0.35491E-01,   0.35539E-01, &
                    0.35572E-01,   0.35597E-01,   0.35631E-01,   0.35667E-01,   0.35699E-01,   0.35730E-01,   0.35808E-01, &
                    0.35913E-01,   0.35995E-01,   0.36111E-01,   0.36217E-01,   0.36243E-01,   0.36256E-01,   0.36290E-01, &
                    0.36331E-01,   0.36378E-01,   0.36409E-01,   0.36440E-01,   0.36489E-01,   0.36488E-01,   0.36491E-01, &
                    0.36507E-01,   0.36460E-01,   0.36440E-01,   0.36439E-01,   0.36382E-01,   0.36367E-01,   0.36393E-01, &
                    0.36426E-01,   0.36458E-01,   0.36494E-01,   0.36530E-01,   0.36567E-01,   0.36606E-01,   0.36634E-01, &
                    0.36702E-01,   0.36761E-01,   0.36792E-01,   0.36755E-01,   0.36817E-01,   0.36844E-01,   0.36870E-01, &
                    0.36909E-01,   0.36917E-01,   0.36931E-01,   0.36930E-01,   0.36967E-01,   0.36992E-01,   0.37021E-01, &
                    0.37044E-01,   0.37105E-01,   0.37245E-01,   0.37286E-01,   0.37293E-01,   0.37277E-01,   0.37289E-01, &
                    0.37296E-01,   0.37253E-01,   0.37267E-01,   0.37287E-01,   0.37314E-01,   0.37354E-01,   0.37377E-01, &
                    0.37397E-01,   0.37410E-01,   0.37379E-01,   0.37380E-01,   0.37435E-01,   0.37462E-01,   0.37513E-01, &
                    0.37616E-01,   0.37650E-01,   0.37643E-01,   0.37671E-01,   0.37714E-01,   0.37757E-01,   0.37806E-01, &
                    0.37835E-01,   0.37931E-01,   0.37985E-01,   0.37997E-01,   0.37975E-01,   0.38001E-01,   0.38049E-01, &
                    0.38111E-01,   0.38182E-01,   0.38230E-01,   0.38257E-01,   0.38300E-01,   0.38352E-01,   0.38405E-01, &
                    0.38425E-01,   0.38464E-01,   0.38496E-01,   0.38543E-01,   0.38580E-01,   0.38623E-01,   0.38657E-01, &
                    0.38683E-01,   0.38722E-01,   0.38759E-01,   0.38799E-01,   0.38818E-01,   0.38817E-01,   0.38854E-01, &
                    0.38915E-01,   0.39016E-01,   0.39047E-01,   0.39090E-01,   0.39133E-01,   0.39141E-01,   0.39111E-01, &
                    0.39082E-01,   0.38966E-01,   0.38914E-01,   0.38997E-01,   0.39056E-01,   0.39085E-01,   0.39110E-01, &
                    0.39134E-01,   0.39171E-01,   0.39200E-01,   0.39226E-01,   0.39335E-01,   0.39404E-01,   0.39563E-01, &
                    0.39712E-01,   0.39751E-01,   0.39777E-01,   0.39799E-01,   0.39798E-01,   0.39737E-01,   0.39688E-01, &
                    0.39727E-01,   0.39757E-01,   0.39782E-01,   0.39816E-01,   0.39804E-01,   0.39776E-01,   0.39810E-01, &
                    0.39850E-01,   0.39789E-01,   0.39794E-01,   0.39894E-01,   0.39956E-01,   0.39981E-01,   0.39982E-01, &
                    0.40003E-01,   0.40010E-01,   0.40050E-01,   0.40089E-01,   0.40048E-01,   0.39937E-01,   0.39809E-01, &
                    0.39412E-01,   0.39229E-01,   0.39175E-01,   0.39156E-01,   0.39143E-01,   0.39173E-01,   0.39215E-01, &
                    0.39263E-01,   0.39331E-01,   0.39441E-01,   0.39618E-01,   0.39978E-01,   0.40424E-01,   0.40637E-01, &
                    0.40735E-01,   0.40810E-01,   0.40880E-01,   0.40921E-01,   0.40965E-01,   0.40993E-01,   0.40978E-01, &
                    0.40966E-01,   0.40987E-01,   0.41003E-01,   0.41036E-01,   0.41092E-01,   0.41107E-01,   0.41113E-01, &
                    0.41149E-01,   0.41174E-01,   0.41114E-01,   0.41109E-01,   0.41169E-01,   0.41217E-01,   0.41243E-01, &
                    0.41288E-01,   0.41291E-01,   0.41304E-01,   0.41338E-01,   0.41340E-01,   0.41312E-01,   0.41352E-01, &
                    0.41401E-01,   0.41391E-01,   0.41389E-01,   0.41388E-01,   0.41372E-01,   0.41353E-01,   0.41377E-01, &
                    0.41348E-01,   0.41318E-01,   0.41382E-01,   0.41424E-01,   0.41487E-01,   0.41546E-01,   0.41592E-01, &
                    0.41605E-01,   0.41572E-01,   0.41649E-01,   0.41685E-01,   0.41751E-01,   0.41856E-01,   0.41891E-01, &
                    0.41884E-01,   0.41890E-01,   0.41916E-01,   0.41923E-01,   0.41953E-01,   0.42041E-01,   0.42036E-01, &
                    0.42030E-01,   0.41970E-01,   0.41934E-01,   0.41921E-01,   0.41842E-01,   0.41835E-01,   0.41866E-01, &
                    0.41914E-01,   0.41957E-01,   0.42003E-01,   0.42067E-01,   0.42118E-01,   0.42182E-01,   0.42191E-01, &
                    0.42244E-01,   0.42383E-01,   0.42425E-01,   0.42409E-01,   0.42436E-01,   0.42481E-01,   0.42525E-01, &
                    0.42554E-01,   0.42620E-01,   0.42694E-01,   0.42808E-01,   0.42828E-01,   0.42794E-01,   0.42857E-01, &
                    0.42934E-01,   0.42978E-01,   0.43014E-01,   0.43069E-01,   0.43089E-01,   0.42936E-01,   0.42959E-01, &
                    0.43010E-01,   0.42988E-01,   0.43141E-01,   0.43170E-01,   0.43225E-01,   0.43269E-01,   0.43315E-01, &
                    0.43350E-01,   0.43466E-01,   0.43602E-01,   0.43609E-01,   0.43738E-01,   0.43850E-01,   0.43907E-01, &
                    0.43949E-01,   0.44019E-01,   0.44076E-01,   0.44139E-01,   0.44216E-01,   0.44312E-01,   0.44384E-01, &
                    0.44458E-01,   0.44509E-01,   0.44424E-01,   0.44375E-01,   0.44454E-01,   0.44450E-01,   0.44391E-01, &
                    0.44450E-01,   0.44510E-01,   0.44574E-01,   0.44538E-01,   0.44499E-01,   0.44548E-01,   0.44793E-01, &
                    0.44844E-01,   0.44893E-01,   0.45038E-01,   0.45082E-01,   0.45095E-01,   0.45070E-01,   0.45027E-01 /)
        solint(4201:4900)=(/                                                                                               &
                    0.45171E-01,   0.45277E-01,   0.45323E-01,   0.45368E-01,   0.45409E-01,   0.45472E-01,   0.45503E-01, &
                    0.45549E-01,   0.45636E-01,   0.45819E-01,   0.45905E-01,   0.45928E-01,   0.45956E-01,   0.46006E-01, &
                    0.45963E-01,   0.45977E-01,   0.45978E-01,   0.45942E-01,   0.45968E-01,   0.46019E-01,   0.46052E-01, &
                    0.46102E-01,   0.46028E-01,   0.46001E-01,   0.46105E-01,   0.46170E-01,   0.46151E-01,   0.46166E-01, &
                    0.46184E-01,   0.46241E-01,   0.46299E-01,   0.46371E-01,   0.46534E-01,   0.46547E-01,   0.46486E-01, &
                    0.46461E-01,   0.46585E-01,   0.46723E-01,   0.46779E-01,   0.46735E-01,   0.46658E-01,   0.46636E-01, &
                    0.46456E-01,   0.46608E-01,   0.46761E-01,   0.46868E-01,   0.46873E-01,   0.46724E-01,   0.46695E-01, &
                    0.46812E-01,   0.46960E-01,   0.47043E-01,   0.47255E-01,   0.47277E-01,   0.47216E-01,   0.47219E-01, &
                    0.47199E-01,   0.47319E-01,   0.47318E-01,   0.47332E-01,   0.47365E-01,   0.47431E-01,   0.47481E-01, &
                    0.47599E-01,   0.47722E-01,   0.47870E-01,   0.48120E-01,   0.48300E-01,   0.48373E-01,   0.48133E-01, &
                    0.48130E-01,   0.48093E-01,   0.48083E-01,   0.47989E-01,   0.47683E-01,   0.47028E-01,   0.46587E-01, &
                    0.46366E-01,   0.46566E-01,   0.46768E-01,   0.46930E-01,   0.46819E-01,   0.46989E-01,   0.47279E-01, &
                    0.47939E-01,   0.48518E-01,   0.48944E-01,   0.48999E-01,   0.48963E-01,   0.48860E-01,   0.48970E-01, &
                    0.48904E-01,   0.48785E-01,   0.48738E-01,   0.48764E-01,   0.48799E-01,   0.48844E-01,   0.48928E-01, &
                    0.49052E-01,   0.49156E-01,   0.49156E-01,   0.49217E-01,   0.49248E-01,   0.49275E-01,   0.49245E-01, &
                    0.49283E-01,   0.49362E-01,   0.49408E-01,   0.49342E-01,   0.49352E-01,   0.49495E-01,   0.49609E-01, &
                    0.49579E-01,   0.49625E-01,   0.49698E-01,   0.49714E-01,   0.49758E-01,   0.49755E-01,   0.49951E-01, &
                    0.50069E-01,   0.50096E-01,   0.50075E-01,   0.50123E-01,   0.50062E-01,   0.49839E-01,   0.49752E-01, &
                    0.49851E-01,   0.49700E-01,   0.49562E-01,   0.49578E-01,   0.49660E-01,   0.49714E-01,   0.49915E-01, &
                    0.50304E-01,   0.50463E-01,   0.50515E-01,   0.50760E-01,   0.51024E-01,   0.51121E-01,   0.51373E-01, &
                    0.51438E-01,   0.51396E-01,   0.51341E-01,   0.51341E-01,   0.50967E-01,   0.49877E-01,   0.49411E-01, &
                    0.49273E-01,   0.49065E-01,   0.49150E-01,   0.49275E-01,   0.49443E-01,   0.49670E-01,   0.50263E-01, &
                    0.51235E-01,   0.51616E-01,   0.51776E-01,   0.51966E-01,   0.51944E-01,   0.51847E-01,   0.51757E-01, &
                    0.51723E-01,   0.51773E-01,   0.51960E-01,   0.51932E-01,   0.51816E-01,   0.51861E-01,   0.51861E-01, &
                    0.51865E-01,   0.51986E-01,   0.52112E-01,   0.52131E-01,   0.52096E-01,   0.52019E-01,   0.51895E-01, &
                    0.51920E-01,   0.51998E-01,   0.52012E-01,   0.51778E-01,   0.51747E-01,   0.51804E-01,   0.51853E-01, &
                    0.52230E-01,   0.52179E-01,   0.52234E-01,   0.52314E-01,   0.52489E-01,   0.52718E-01,   0.52614E-01, &
                    0.52759E-01,   0.52860E-01,   0.53064E-01,   0.53018E-01,   0.53126E-01,   0.53250E-01,   0.53314E-01, &
                    0.53495E-01,   0.53636E-01,   0.53464E-01,   0.53353E-01,   0.53408E-01,   0.53319E-01,   0.52973E-01, &
                    0.52280E-01,   0.52001E-01,   0.51880E-01,   0.51982E-01,   0.52229E-01,   0.52414E-01,   0.52595E-01, &
                    0.52993E-01,   0.53746E-01,   0.54286E-01,   0.54496E-01,   0.54608E-01,   0.54661E-01,   0.54704E-01, &
                    0.54732E-01,   0.54584E-01,   0.54463E-01,   0.54362E-01,   0.54389E-01,   0.54444E-01,   0.54441E-01, &
                    0.54376E-01,   0.54377E-01,   0.54484E-01,   0.54585E-01,   0.54655E-01,   0.54650E-01,   0.54676E-01, &
                    0.54731E-01,   0.54833E-01,   0.54940E-01,   0.55040E-01,   0.55119E-01,   0.55087E-01,   0.55107E-01, &
                    0.55073E-01,   0.55055E-01,   0.55082E-01,   0.55052E-01,   0.54992E-01,   0.54964E-01,   0.54958E-01, &
                    0.55002E-01,   0.55099E-01,   0.55245E-01,   0.55317E-01,   0.55429E-01,   0.55481E-01,   0.55582E-01, &
                    0.55777E-01,   0.55837E-01,   0.55883E-01,   0.55919E-01,   0.55943E-01,   0.55880E-01,   0.55922E-01, &
                    0.55798E-01,   0.55766E-01,   0.55853E-01,   0.55917E-01,   0.55896E-01,   0.56001E-01,   0.56130E-01, &
                    0.56219E-01,   0.56465E-01,   0.56585E-01,   0.56544E-01,   0.56480E-01,   0.56485E-01,   0.56488E-01, &
                    0.56602E-01,   0.56517E-01,   0.56484E-01,   0.56494E-01,   0.56692E-01,   0.56910E-01,   0.57078E-01, &
                    0.56932E-01,   0.57059E-01,   0.57207E-01,   0.57234E-01,   0.56832E-01,   0.56288E-01,   0.56114E-01, &
                    0.55968E-01,   0.55691E-01,   0.55748E-01,   0.55988E-01,   0.56423E-01,   0.57177E-01,   0.57578E-01, &
                    0.57929E-01,   0.58278E-01,   0.58632E-01,   0.58601E-01,   0.58221E-01,   0.58178E-01,   0.58109E-01, &
                    0.58160E-01,   0.58151E-01,   0.58197E-01,   0.58242E-01,   0.58493E-01,   0.58769E-01,   0.58869E-01, &
                    0.59055E-01,   0.59124E-01,   0.58531E-01,   0.58616E-01,   0.58593E-01,   0.58437E-01,   0.58316E-01, &
                    0.58281E-01,   0.58319E-01,   0.58853E-01,   0.58884E-01,   0.58898E-01,   0.58911E-01,   0.59041E-01, &
                    0.59078E-01,   0.59170E-01,   0.59166E-01,   0.59320E-01,   0.59535E-01,   0.59706E-01,   0.59867E-01, &
                    0.60036E-01,   0.60094E-01,   0.60167E-01,   0.59995E-01,   0.60047E-01,   0.60004E-01,   0.59981E-01, &
                    0.59976E-01,   0.59991E-01,   0.59889E-01,   0.59861E-01,   0.59301E-01,   0.58907E-01,   0.58837E-01, &
                    0.58858E-01,   0.58982E-01,   0.59089E-01,   0.59210E-01,   0.59884E-01,   0.60538E-01,   0.60456E-01, &
                    0.60562E-01,   0.60663E-01,   0.60870E-01,   0.60901E-01,   0.60847E-01,   0.60785E-01,   0.61044E-01, &
                    0.61062E-01,   0.60963E-01,   0.60877E-01,   0.61189E-01,   0.61299E-01,   0.61291E-01,   0.61368E-01, &
                    0.61407E-01,   0.61458E-01,   0.61493E-01,   0.61390E-01,   0.61395E-01,   0.61475E-01,   0.61465E-01, &
                    0.61417E-01,   0.61386E-01,   0.61366E-01,   0.61505E-01,   0.61561E-01,   0.61582E-01,   0.61544E-01, &
                    0.61488E-01,   0.61542E-01,   0.61277E-01,   0.61019E-01,   0.60842E-01,   0.60474E-01,   0.60665E-01, &
                    0.60608E-01,   0.60979E-01,   0.61469E-01,   0.61631E-01,   0.61994E-01,   0.62573E-01,   0.62869E-01, &
                    0.62912E-01,   0.62980E-01,   0.63125E-01,   0.63294E-01,   0.63248E-01,   0.63143E-01,   0.63285E-01, &
                    0.63265E-01,   0.62987E-01,   0.62462E-01,   0.62388E-01,   0.62455E-01,   0.62519E-01,   0.62150E-01, &
                    0.61882E-01,   0.62172E-01,   0.62682E-01,   0.62676E-01,   0.62506E-01,   0.62447E-01,   0.62474E-01, &
                    0.62726E-01,   0.62911E-01,   0.62711E-01,   0.62960E-01,   0.63034E-01,   0.63666E-01,   0.63694E-01, &
                    0.63558E-01,   0.63414E-01,   0.63624E-01,   0.63696E-01,   0.63955E-01,   0.63709E-01,   0.63647E-01, &
                    0.64125E-01,   0.64209E-01,   0.64450E-01,   0.64185E-01,   0.63682E-01,   0.63459E-01,   0.62517E-01, &
                    0.62204E-01,   0.62087E-01,   0.62174E-01,   0.62878E-01,   0.63278E-01,   0.64234E-01,   0.64397E-01, &
                    0.64780E-01,   0.64916E-01,   0.64805E-01,   0.64753E-01,   0.63862E-01,   0.64033E-01,   0.63274E-01, &
                    0.63217E-01,   0.63163E-01,   0.63427E-01,   0.64591E-01,   0.64463E-01,   0.64945E-01,   0.64918E-01, &
                    0.64813E-01,   0.64590E-01,   0.64517E-01,   0.64488E-01,   0.64566E-01,   0.64833E-01,   0.64982E-01, &
                    0.64974E-01,   0.64918E-01,   0.64992E-01,   0.65014E-01,   0.64689E-01,   0.64476E-01,   0.64354E-01, &
                    0.64407E-01,   0.64574E-01,   0.64496E-01,   0.64696E-01,   0.64821E-01,   0.64823E-01,   0.64840E-01, &
                    0.64891E-01,   0.65077E-01,   0.65299E-01,   0.65305E-01,   0.65558E-01,   0.65445E-01,   0.65011E-01, &
                    0.64936E-01,   0.64701E-01,   0.64868E-01,   0.64984E-01,   0.65087E-01,   0.65266E-01,   0.65270E-01, &
                    0.65536E-01,   0.65552E-01,   0.64920E-01,   0.65039E-01,   0.65319E-01,   0.65346E-01,   0.65377E-01, &
                    0.65158E-01,   0.65753E-01,   0.65719E-01,   0.65421E-01,   0.65255E-01,   0.65230E-01,   0.65404E-01, &
                    0.65528E-01,   0.65513E-01,   0.65653E-01,   0.65887E-01,   0.65717E-01,   0.65776E-01,   0.65605E-01, &
                    0.65943E-01,   0.66117E-01,   0.66127E-01,   0.66167E-01,   0.66391E-01,   0.66162E-01,   0.65209E-01, &
                    0.64473E-01,   0.63842E-01,   0.63808E-01,   0.63995E-01,   0.64721E-01,   0.65719E-01,   0.66050E-01, &
                    0.67073E-01,   0.66649E-01,   0.66536E-01,   0.66516E-01,   0.66682E-01,   0.66617E-01,   0.66731E-01, &
                    0.66940E-01,   0.66888E-01,   0.65288E-01,   0.65123E-01,   0.65110E-01,   0.65183E-01,   0.65140E-01, &
                    0.65392E-01,   0.66488E-01,   0.66499E-01,   0.66460E-01,   0.66452E-01,   0.66385E-01,   0.66356E-01, &
                    0.66295E-01,   0.66177E-01,   0.65770E-01,   0.65781E-01,   0.65805E-01,   0.65831E-01,   0.65904E-01, &
                    0.66060E-01,   0.66356E-01,   0.66537E-01,   0.66556E-01,   0.66619E-01,   0.66413E-01,   0.66193E-01, &
                    0.66142E-01,   0.66075E-01,   0.66319E-01,   0.66336E-01,   0.66604E-01,   0.66650E-01,   0.66610E-01, &
                    0.66508E-01,   0.66234E-01,   0.66047E-01,   0.66067E-01,   0.65830E-01,   0.65912E-01,   0.66258E-01, &
                    0.66427E-01,   0.66442E-01,   0.66666E-01,   0.66961E-01,   0.66632E-01,   0.66530E-01,   0.66102E-01, &
                    0.66357E-01,   0.65469E-01,   0.65874E-01,   0.66000E-01,   0.66553E-01,   0.66527E-01,   0.67320E-01, &
                    0.67339E-01,   0.67198E-01,   0.67274E-01,   0.67059E-01,   0.67061E-01,   0.66622E-01,   0.66714E-01, &
                    0.66363E-01,   0.66318E-01,   0.66136E-01,   0.66407E-01,   0.66478E-01,   0.66566E-01,   0.66537E-01, &
                    0.66774E-01,   0.66850E-01,   0.66919E-01,   0.67232E-01,   0.67380E-01,   0.67459E-01,   0.67432E-01, &
                    0.67310E-01,   0.67105E-01,   0.66797E-01,   0.66722E-01,   0.66883E-01,   0.67019E-01,   0.67140E-01, &
                    0.67445E-01,   0.67448E-01,   0.67310E-01,   0.67184E-01,   0.67172E-01,   0.67156E-01,   0.66943E-01, &
                    0.66893E-01,   0.67092E-01,   0.67179E-01,   0.67339E-01,   0.67511E-01,   0.67508E-01,   0.67458E-01, &
                    0.67356E-01,   0.67372E-01,   0.67437E-01,   0.67425E-01,   0.67459E-01,   0.67474E-01,   0.67368E-01, &
                    0.67443E-01,   0.67439E-01,   0.67397E-01,   0.67305E-01,   0.67286E-01,   0.67420E-01,   0.67396E-01, &
                    0.67476E-01,   0.67648E-01,   0.67484E-01,   0.67345E-01,   0.67186E-01,   0.67321E-01,   0.67400E-01, &
                    0.67588E-01,   0.67539E-01,   0.67508E-01,   0.67297E-01,   0.67185E-01,   0.67270E-01,   0.67481E-01, &
                    0.67600E-01,   0.67682E-01,   0.67571E-01,   0.67667E-01,   0.67591E-01,   0.67466E-01,   0.67660E-01, &
                    0.67663E-01,   0.67308E-01,   0.67270E-01,   0.67331E-01,   0.67330E-01,   0.67546E-01,   0.67699E-01, &
                    0.67800E-01,   0.67676E-01,   0.67672E-01,   0.67633E-01,   0.67648E-01,   0.67842E-01,   0.67975E-01, &
                    0.67962E-01,   0.67982E-01,   0.68051E-01,   0.68109E-01,   0.67763E-01,   0.67777E-01,   0.67820E-01, &
                    0.67968E-01,   0.68014E-01,   0.68490E-01,   0.68638E-01,   0.68617E-01,   0.68395E-01,   0.68470E-01, &
                    0.68305E-01,   0.67888E-01,   0.68350E-01,   0.68261E-01,   0.68427E-01,   0.68920E-01,   0.69051E-01 /)
        solint(4901:5600)=(/                                                                                               &
                    0.69175E-01,   0.69212E-01,   0.69225E-01,   0.69266E-01,   0.69421E-01,   0.69318E-01,   0.69471E-01, &
                    0.68946E-01,   0.68746E-01,   0.68937E-01,   0.68344E-01,   0.68487E-01,   0.69155E-01,   0.69254E-01, &
                    0.69702E-01,   0.69996E-01,   0.70036E-01,   0.70084E-01,   0.70002E-01,   0.69779E-01,   0.69650E-01, &
                    0.69514E-01,   0.69555E-01,   0.69801E-01,   0.70111E-01,   0.70093E-01,   0.70149E-01,   0.70236E-01, &
                    0.70218E-01,   0.70255E-01,   0.70384E-01,   0.70393E-01,   0.70349E-01,   0.70378E-01,   0.70470E-01, &
                    0.70534E-01,   0.70769E-01,   0.71093E-01,   0.71036E-01,   0.70920E-01,   0.69986E-01,   0.67813E-01, &
                    0.67497E-01,   0.67727E-01,   0.68620E-01,   0.71148E-01,   0.71840E-01,   0.71993E-01,   0.71992E-01, &
                    0.71778E-01,   0.71620E-01,   0.71486E-01,   0.71460E-01,   0.71514E-01,   0.71383E-01,   0.71055E-01, &
                    0.70799E-01,   0.70843E-01,   0.70819E-01,   0.70950E-01,   0.71205E-01,   0.70947E-01,   0.71093E-01, &
                    0.71309E-01,   0.71187E-01,   0.71418E-01,   0.71291E-01,   0.71238E-01,   0.71331E-01,   0.71309E-01, &
                    0.71410E-01,   0.71280E-01,   0.71291E-01,   0.71368E-01,   0.71537E-01,   0.71767E-01,   0.71688E-01, &
                    0.71570E-01,   0.71544E-01,   0.71255E-01,   0.71072E-01,   0.71048E-01,   0.71005E-01,   0.71106E-01, &
                    0.71184E-01,   0.71281E-01,   0.71299E-01,   0.71566E-01,   0.71613E-01,   0.71539E-01,   0.71552E-01, &
                    0.71451E-01,   0.71659E-01,   0.71720E-01,   0.71850E-01,   0.72054E-01,   0.71604E-01,   0.71648E-01, &
                    0.71597E-01,   0.71988E-01,   0.72078E-01,   0.72005E-01,   0.72010E-01,   0.72001E-01,   0.72018E-01, &
                    0.72067E-01,   0.72089E-01,   0.72187E-01,   0.72346E-01,   0.72364E-01,   0.72368E-01,   0.72360E-01, &
                    0.71969E-01,   0.71755E-01,   0.71093E-01,   0.70455E-01,   0.71217E-01,   0.71313E-01,   0.72783E-01, &
                    0.72725E-01,   0.71632E-01,   0.71845E-01,   0.71913E-01,   0.73055E-01,   0.71816E-01,   0.71274E-01, &
                    0.70536E-01,   0.71050E-01,   0.71801E-01,   0.72580E-01,   0.72614E-01,   0.72840E-01,   0.72875E-01, &
                    0.72646E-01,   0.71938E-01,   0.71131E-01,   0.71100E-01,   0.71432E-01,   0.72576E-01,   0.72330E-01, &
                    0.71205E-01,   0.71201E-01,   0.71204E-01,   0.71509E-01,   0.72526E-01,   0.72196E-01,   0.71900E-01, &
                    0.71593E-01,   0.71300E-01,   0.71607E-01,   0.71530E-01,   0.72323E-01,   0.72543E-01,   0.72452E-01, &
                    0.71873E-01,   0.71871E-01,   0.71699E-01,   0.71846E-01,   0.71903E-01,   0.71533E-01,   0.71707E-01, &
                    0.72045E-01,   0.71303E-01,   0.71522E-01,   0.70844E-01,   0.71482E-01,   0.71910E-01,   0.72503E-01, &
                    0.72452E-01,   0.72492E-01,   0.72409E-01,   0.72232E-01,   0.72157E-01,   0.71983E-01,   0.72109E-01, &
                    0.72117E-01,   0.72226E-01,   0.72353E-01,   0.72412E-01,   0.72019E-01,   0.71926E-01,   0.71710E-01, &
                    0.71824E-01,   0.71372E-01,   0.71467E-01,   0.71358E-01,   0.71644E-01,   0.72099E-01,   0.72334E-01, &
                    0.72648E-01,   0.72023E-01,   0.72220E-01,   0.71975E-01,   0.72275E-01,   0.71776E-01,   0.71849E-01, &
                    0.71832E-01,   0.72300E-01,   0.71441E-01,   0.71543E-01,   0.71670E-01,   0.72636E-01,   0.72646E-01, &
                    0.72123E-01,   0.71660E-01,   0.71618E-01,   0.71890E-01,   0.72162E-01,   0.72269E-01,   0.72331E-01, &
                    0.71959E-01,   0.71832E-01,   0.71683E-01,   0.71966E-01,   0.72221E-01,   0.72297E-01,   0.72148E-01, &
                    0.72133E-01,   0.72214E-01,   0.72219E-01,   0.72329E-01,   0.71922E-01,   0.71866E-01,   0.71968E-01, &
                    0.72557E-01,   0.71919E-01,   0.71913E-01,   0.71247E-01,   0.71858E-01,   0.70243E-01,   0.69955E-01, &
                    0.70564E-01,   0.72247E-01,   0.73463E-01,   0.73325E-01,   0.72874E-01,   0.71544E-01,   0.71512E-01, &
                    0.72242E-01,   0.73142E-01,   0.72683E-01,   0.71368E-01,   0.70683E-01,   0.71446E-01,   0.72268E-01, &
                    0.72294E-01,   0.72295E-01,   0.72554E-01,   0.72321E-01,   0.72141E-01,   0.71531E-01,   0.72011E-01, &
                    0.72033E-01,   0.72594E-01,   0.71391E-01,   0.70715E-01,   0.70760E-01,   0.71473E-01,   0.72309E-01, &
                    0.72353E-01,   0.72852E-01,   0.72536E-01,   0.72275E-01,   0.71939E-01,   0.71950E-01,   0.71035E-01, &
                    0.71663E-01,   0.71748E-01,   0.72658E-01,   0.72632E-01,   0.72313E-01,   0.72290E-01,   0.72210E-01, &
                    0.72488E-01,   0.72577E-01,   0.72576E-01,   0.72465E-01,   0.72251E-01,   0.71764E-01,   0.72069E-01, &
                    0.72606E-01,   0.72685E-01,   0.72746E-01,   0.72846E-01,   0.72726E-01,   0.72739E-01,   0.72376E-01, &
                    0.72619E-01,   0.72805E-01,   0.72707E-01,   0.72492E-01,   0.72352E-01,   0.72794E-01,   0.73301E-01, &
                    0.73128E-01,   0.73102E-01,   0.73222E-01,   0.73496E-01,   0.73461E-01,   0.73461E-01,   0.73501E-01, &
                    0.73566E-01,   0.73045E-01,   0.72911E-01,   0.72932E-01,   0.72590E-01,   0.72883E-01,   0.73071E-01, &
                    0.73321E-01,   0.73064E-01,   0.73545E-01,   0.73689E-01,   0.73893E-01,   0.73839E-01,   0.74039E-01, &
                    0.74060E-01,   0.74144E-01,   0.73990E-01,   0.72597E-01,   0.71207E-01,   0.71460E-01,   0.73742E-01, &
                    0.74372E-01,   0.74619E-01,   0.74405E-01,   0.74227E-01,   0.74262E-01,   0.74307E-01,   0.74212E-01, &
                    0.74113E-01,   0.74083E-01,   0.74166E-01,   0.73911E-01,   0.73946E-01,   0.73436E-01,   0.73343E-01, &
                    0.73517E-01,   0.73859E-01,   0.73770E-01,   0.73927E-01,   0.74049E-01,   0.73943E-01,   0.74219E-01, &
                    0.74380E-01,   0.74295E-01,   0.73489E-01,   0.73870E-01,   0.73834E-01,   0.74069E-01,   0.74032E-01, &
                    0.73759E-01,   0.73815E-01,   0.74560E-01,   0.74393E-01,   0.74495E-01,   0.74513E-01,   0.74126E-01, &
                    0.74359E-01,   0.73781E-01,   0.73787E-01,   0.73961E-01,   0.73851E-01,   0.74395E-01,   0.74345E-01, &
                    0.74581E-01,   0.74847E-01,   0.74961E-01,   0.74498E-01,   0.74135E-01,   0.72222E-01,   0.72131E-01, &
                    0.74360E-01,   0.74944E-01,   0.74801E-01,   0.74865E-01,   0.75185E-01,   0.75232E-01,   0.74834E-01, &
                    0.74675E-01,   0.74244E-01,   0.73569E-01,   0.74018E-01,   0.72881E-01,   0.72202E-01,   0.73530E-01, &
                    0.74825E-01,   0.74882E-01,   0.75019E-01,   0.74573E-01,   0.74163E-01,   0.74713E-01,   0.75021E-01, &
                    0.74760E-01,   0.74710E-01,   0.74725E-01,   0.74544E-01,   0.74971E-01,   0.75113E-01,   0.73476E-01, &
                    0.73168E-01,   0.73927E-01,   0.72193E-01,   0.71939E-01,   0.73629E-01,   0.74278E-01,   0.75348E-01, &
                    0.75060E-01,   0.74636E-01,   0.75086E-01,   0.75055E-01,   0.75049E-01,   0.75722E-01,   0.74583E-01, &
                    0.73703E-01,   0.73620E-01,   0.73105E-01,   0.73194E-01,   0.74412E-01,   0.74745E-01,   0.74912E-01, &
                    0.75757E-01,   0.75169E-01,   0.72655E-01,   0.71581E-01,   0.72831E-01,   0.74274E-01,   0.74967E-01, &
                    0.74621E-01,   0.75144E-01,   0.74383E-01,   0.74151E-01,   0.74091E-01,   0.73800E-01,   0.74829E-01, &
                    0.75002E-01,   0.74911E-01,   0.74570E-01,   0.73200E-01,   0.72388E-01,   0.74047E-01,   0.74529E-01, &
                    0.73658E-01,   0.74079E-01,   0.72736E-01,   0.72669E-01,   0.74512E-01,   0.74620E-01,   0.74727E-01, &
                    0.74332E-01,   0.72973E-01,   0.73043E-01,   0.73290E-01,   0.74355E-01,   0.72544E-01,   0.73602E-01, &
                    0.73610E-01,   0.72645E-01,   0.71952E-01,   0.65476E-01,   0.68412E-01,   0.74010E-01,   0.74573E-01, &
                    0.75493E-01,   0.72978E-01,   0.72640E-01,   0.72470E-01,   0.72087E-01,   0.71714E-01,   0.71746E-01, &
                    0.69763E-01,   0.58431E-01,   0.64482E-01,   0.71151E-01,   0.70370E-01,   0.69733E-01,   0.67774E-01, &
                    0.71922E-01,   0.72319E-01,   0.73310E-01,   0.73922E-01,   0.72693E-01,   0.70893E-01,   0.72183E-01, &
                    0.72952E-01,   0.73913E-01,   0.73129E-01,   0.72284E-01,   0.72894E-01,   0.73416E-01,   0.73237E-01, &
                    0.71710E-01,   0.71001E-01,   0.72183E-01,   0.73523E-01,   0.72812E-01,   0.73183E-01,   0.73301E-01, &
                    0.73576E-01,   0.73540E-01,   0.73082E-01,   0.73039E-01,   0.73052E-01,   0.73049E-01,   0.71530E-01, &
                    0.72092E-01,   0.71624E-01,   0.72272E-01,   0.73146E-01,   0.73639E-01,   0.73728E-01,   0.73706E-01, &
                    0.74104E-01,   0.73869E-01,   0.73412E-01,   0.73004E-01,   0.72076E-01,   0.72461E-01,   0.72658E-01, &
                    0.72775E-01,   0.72266E-01,   0.72375E-01,   0.72620E-01,   0.72482E-01,   0.72713E-01,   0.72013E-01, &
                    0.71802E-01,   0.72567E-01,   0.72492E-01,   0.72418E-01,   0.72000E-01,   0.70657E-01,   0.70150E-01, &
                    0.71019E-01,   0.71672E-01,   0.72423E-01,   0.72404E-01,   0.72059E-01,   0.71837E-01,   0.72015E-01, &
                    0.71781E-01,   0.71402E-01,   0.70892E-01,   0.71586E-01,   0.71774E-01,   0.71173E-01,   0.71309E-01, &
                    0.70987E-01,   0.70777E-01,   0.71198E-01,   0.71716E-01,   0.72236E-01,   0.72758E-01,   0.73210E-01, &
                    0.73127E-01,   0.73218E-01,   0.73052E-01,   0.73288E-01,   0.73183E-01,   0.71480E-01,   0.72282E-01, &
                    0.73780E-01,   0.73806E-01,   0.72658E-01,   0.72465E-01,   0.73012E-01,   0.73703E-01,   0.72497E-01, &
                    0.73434E-01,   0.73123E-01,   0.72747E-01,   0.72801E-01,   0.72423E-01,   0.71935E-01,   0.71619E-01, &
                    0.71506E-01,   0.71546E-01,   0.72211E-01,   0.72471E-01,   0.72064E-01,   0.72184E-01,   0.71267E-01, &
                    0.71840E-01,   0.70645E-01,   0.69437E-01,   0.69803E-01,   0.70042E-01,   0.69792E-01,   0.71242E-01, &
                    0.71183E-01,   0.71234E-01,   0.71836E-01,   0.71652E-01,   0.71489E-01,   0.71378E-01,   0.71384E-01, &
                    0.69931E-01,   0.70079E-01,   0.71244E-01,   0.71224E-01,   0.70872E-01,   0.71317E-01,   0.71336E-01, &
                    0.70880E-01,   0.70808E-01,   0.70614E-01,   0.70324E-01,   0.70138E-01,   0.70907E-01,   0.70834E-01, &
                    0.70372E-01,   0.70476E-01,   0.70460E-01,   0.70154E-01,   0.70268E-01,   0.71081E-01,   0.70647E-01, &
                    0.70960E-01,   0.71441E-01,   0.71200E-01,   0.71905E-01,   0.72212E-01,   0.71317E-01,   0.70813E-01, &
                    0.70354E-01,   0.70653E-01,   0.71468E-01,   0.71115E-01,   0.70101E-01,   0.70905E-01,   0.70530E-01, &
                    0.70437E-01,   0.70669E-01,   0.70388E-01,   0.70477E-01,   0.70682E-01,   0.70399E-01,   0.69875E-01, &
                    0.69603E-01,   0.69458E-01,   0.69329E-01,   0.69197E-01,   0.69197E-01,   0.69080E-01,   0.68971E-01, &
                    0.68868E-01,   0.68777E-01,   0.68691E-01,   0.68613E-01,   0.68545E-01,   0.68483E-01,   0.68430E-01, &
                    0.68384E-01,   0.68348E-01,   0.68314E-01,   0.68361E-01,   0.68527E-01,   0.68621E-01,   0.68632E-01, &
                    0.68035E-01,   0.68053E-01,   0.68277E-01,   0.67940E-01,   0.67596E-01,   0.68249E-01,   0.68318E-01, &
                    0.67763E-01,   0.67030E-01,   0.66638E-01,   0.63414E-01,   0.57160E-01,   0.62048E-01,   0.66267E-01, &
                    0.67164E-01,   0.67553E-01,   0.68379E-01,   0.67865E-01,   0.65383E-01,   0.67219E-01,   0.67243E-01, &
                    0.66608E-01,   0.66360E-01,   0.66775E-01,   0.66515E-01,   0.66577E-01,   0.65718E-01,   0.65831E-01, &
                    0.66189E-01,   0.67376E-01,   0.67249E-01,   0.66171E-01,   0.66185E-01,   0.66052E-01,   0.65672E-01, &
                    0.65875E-01,   0.64911E-01,   0.65692E-01,   0.65847E-01,   0.66445E-01,   0.66257E-01,   0.65297E-01 /)
        solint(5601:6064)=(/                                                                                               &
                    0.63888E-01,   0.63959E-01,   0.64389E-01,   0.64594E-01,   0.64842E-01,   0.65137E-01,   0.64398E-01, &
                    0.63851E-01,   0.62765E-01,   0.61654E-01,   0.62679E-01,   0.63033E-01,   0.62931E-01,   0.63582E-01, &
                    0.63730E-01,   0.63476E-01,   0.63306E-01,   0.63505E-01,   0.63959E-01,   0.64074E-01,   0.64166E-01, &
                    0.64383E-01,   0.64003E-01,   0.62708E-01,   0.62297E-01,   0.62899E-01,   0.63151E-01,   0.62999E-01, &
                    0.62999E-01,   0.63420E-01,   0.63340E-01,   0.62747E-01,   0.62873E-01,   0.62513E-01,   0.62450E-01, &
                    0.62158E-01,   0.60570E-01,   0.57640E-01,   0.60053E-01,   0.62285E-01,   0.61289E-01,   0.61255E-01, &
                    0.62220E-01,   0.62522E-01,   0.62305E-01,   0.61536E-01,   0.61089E-01,   0.60395E-01,   0.59789E-01, &
                    0.60456E-01,   0.60191E-01,   0.59795E-01,   0.60538E-01,   0.60902E-01,   0.60459E-01,   0.58054E-01, &
                    0.57812E-01,   0.58252E-01,   0.58255E-01,   0.58580E-01,   0.57414E-01,   0.57234E-01,   0.57600E-01, &
                    0.57303E-01,   0.56879E-01,   0.57391E-01,   0.56732E-01,   0.55240E-01,   0.55853E-01,   0.56442E-01, &
                    0.56997E-01,   0.58037E-01,   0.58181E-01,   0.57451E-01,   0.57474E-01,   0.57416E-01,   0.57069E-01, &
                    0.57216E-01,   0.56395E-01,   0.55983E-01,   0.56498E-01,   0.56340E-01,   0.56065E-01,   0.55435E-01, &
                    0.54830E-01,   0.54062E-01,   0.52774E-01,   0.53703E-01,   0.54803E-01,   0.54302E-01,   0.54164E-01, &
                    0.55276E-01,   0.53552E-01,   0.53253E-01,   0.53871E-01,   0.54331E-01,   0.54087E-01,   0.53184E-01, &
                    0.52529E-01,   0.49230E-01,   0.47958E-01,   0.52251E-01,   0.52701E-01,   0.51019E-01,   0.50007E-01, &
                    0.50902E-01,   0.49849E-01,   0.48810E-01,   0.45479E-01,   0.43706E-01,   0.46365E-01,   0.49678E-01, &
                    0.48818E-01,   0.49315E-01,   0.50406E-01,   0.50824E-01,   0.49712E-01,   0.50138E-01,   0.50817E-01, &
                    0.51453E-01,   0.52086E-01,   0.51160E-01,   0.49221E-01,   0.49296E-01,   0.48257E-01,   0.47138E-01, &
                    0.47974E-01,   0.48628E-01,   0.48628E-01,   0.48581E-01,   0.49269E-01,   0.49535E-01,   0.49592E-01, &
                    0.49451E-01,   0.48173E-01,   0.47480E-01,   0.47471E-01,   0.49074E-01,   0.48256E-01,   0.46249E-01, &
                    0.44735E-01,   0.41228E-01,   0.43540E-01,   0.47739E-01,   0.48190E-01,   0.48699E-01,   0.49029E-01, &
                    0.48543E-01,   0.48662E-01,   0.48577E-01,   0.48474E-01,   0.47121E-01,   0.47657E-01,   0.47568E-01, &
                    0.46437E-01,   0.46933E-01,   0.45421E-01,   0.44219E-01,   0.44759E-01,   0.44981E-01,   0.44084E-01, &
                    0.44216E-01,   0.44459E-01,   0.44350E-01,   0.45247E-01,   0.45649E-01,   0.45107E-01,   0.44301E-01, &
                    0.43517E-01,   0.43684E-01,   0.44677E-01,   0.44556E-01,   0.43274E-01,   0.42882E-01,   0.40908E-01, &
                    0.41882E-01,   0.43839E-01,   0.43538E-01,   0.41575E-01,   0.41066E-01,   0.40475E-01,   0.40474E-01, &
                    0.38148E-01,   0.38310E-01,   0.38879E-01,   0.38760E-01,   0.38567E-01,   0.37341E-01,   0.35895E-01, &
                    0.34316E-01,   0.32868E-01,   0.34125E-01,   0.35133E-01,   0.33191E-01,   0.30816E-01,   0.31409E-01, &
                    0.31767E-01,   0.28699E-01,   0.23112E-01,   0.27611E-01,   0.30352E-01,   0.30307E-01,   0.31412E-01, &
                    0.31486E-01,   0.31646E-01,   0.30833E-01,   0.31149E-01,   0.32194E-01,   0.31526E-01,   0.30709E-01, &
                    0.30621E-01,   0.30575E-01,   0.31135E-01,   0.30680E-01,   0.30545E-01,   0.30369E-01,   0.29199E-01, &
                    0.28040E-01,   0.26889E-01,   0.28965E-01,   0.28663E-01,   0.28212E-01,   0.27768E-01,   0.28459E-01, &
                    0.29222E-01,   0.29843E-01,   0.30411E-01,   0.29885E-01,   0.27386E-01,   0.27211E-01,   0.25269E-01, &
                    0.13285E-01,   0.11758E-01,   0.22054E-01,   0.17043E-01,   0.73968E-02,   0.16821E-01,   0.22183E-01, &
                    0.19022E-01,   0.18566E-01,   0.14526E-01,   0.15072E-01,   0.15408E-01,   0.15175E-01,   0.14392E-01, &
                    0.99069E-02,   0.11866E-01,   0.16473E-01,   0.18092E-01,   0.15415E-01,   0.20105E-01,   0.18700E-01, &
                    0.15488E-01,   0.15771E-01,   0.13081E-01,   0.12802E-01,   0.14800E-01,   0.17861E-01,   0.15407E-01, &
                    0.17355E-01,   0.14717E-01,   0.16464E-01,   0.17038E-01,   0.15985E-01,   0.13897E-01,   0.12862E-01, &
                    0.12862E-01,   0.13602E-01,   0.11869E-01,   0.13361E-01,   0.14152E-01,   0.82699E-02,   0.11041E-01, &
                    0.12134E-01,   0.13538E-01,   0.14407E-01,   0.13513E-01,   0.11182E-01,   0.12454E-01,   0.13927E-01, &
                    0.10528E-01,   0.11656E-01,   0.11160E-01,   0.11228E-01,   0.12069E-01,   0.88372E-02,   0.12214E-01, &
                    0.12126E-01,   0.10885E-01,   0.12386E-01,   0.11203E-01,   0.10733E-01,   0.99908E-02,   0.90809E-02, &
                    0.11448E-01,   0.11204E-01,   0.10463E-01,   0.10898E-01,   0.11215E-01,   0.11302E-01,   0.12220E-01, &
                    0.99253E-02,   0.10656E-01,   0.10958E-01,   0.95999E-02,   0.85605E-02,   0.71318E-02,   0.76340E-02, &
                    0.72718E-02,   0.91489E-02,   0.73963E-02,   0.66991E-02,   0.86940E-02,   0.62377E-02,   0.68144E-02, &
                    0.66127E-02,   0.71753E-02,   0.66052E-02,   0.75957E-02,   0.63924E-02,   0.48096E-02,   0.62842E-02, &
                    0.62842E-02,   0.62181E-02,   0.53423E-02,   0.57619E-02,   0.56355E-02,   0.60386E-02,   0.42697E-02, &
                    0.44619E-02,   0.35440E-02,   0.47516E-02,   0.47517E-02,   0.37358E-02,   0.48814E-02,   0.44973E-02, &
                    0.49263E-02,   0.43174E-02,   0.48302E-02,   0.42550E-02,   0.50523E-02,   0.52832E-02,   0.39073E-02, &
                    0.39073E-02,   0.26046E-02,   0.30602E-02,   0.28091E-02,   0.10740E-02,   0.21046E-02,   0.26190E-02, &
                    0.23908E-02,   0.23909E-02,   0.17465E-02,   0.66757E-03,   0.59927E-03,   0.12490E-02,   0.19163E-02, &
                    0.19264E-02,   0.12955E-02,   0.90463E-03,   0.16344E-02,   0.16344E-02,   0.14468E-02,   0.17197E-02, &
                    0.20542E-02,   0.16730E-02,   0.17959E-02,   0.18217E-02,   0.17703E-02,   0.18683E-02,   0.18664E-02, &
                    0.11944E-02,   0.74471E-03,   0.59185E-03,   0.58978E-03,   0.72372E-03,   0.88406E-03,   0.82886E-03, &
                    0.66739E-03,   0.53119E-03,   0.36814E-03,   0.33366E-03,   0.25543E-03,   0.29257E-03,   0.37853E-03, &
                    0.35740E-03,   0.35742E-03,   0.27860E-03,   0.35558E-03,   0.29886E-03,   0.30039E-03,   0.38542E-03, &
                    0.38950E-03,   0.43804E-03,   0.28093E-03,   0.23514E-03,   0.26178E-03,   0.22044E-03,   0.31071E-03, &
                    0.25772E-03,   0.31379E-03,   0.19480E-03,   0.24348E-03,   0.29402E-03,   0.26131E-03,   0.29778E-03, &
                    0.29779E-03,   0.25090E-03,   0.29276E-03,   0.19937E-03,   0.20413E-03,   0.27195E-03,   0.30658E-03, &
                    0.33691E-03,   0.25220E-03,   0.17542E-03,   0.24820E-03,   0.23730E-03,   0.22725E-03,   0.16054E-03, &
                    0.16003E-03,   0.16747E-03,   0.20442E-03,   0.14658E-03,   0.14908E-03,   0.16289E-03,   0.13042E-03, &
                    0.97990E-04,   0.67688E-04,   0.59302E-04,   0.50648E-04,   0.48708E-04,   0.46505E-04,   0.41213E-04, &
                    0.35264E-04,   0.35420E-04,   0.31670E-04,   0.28626E-04,   0.26351E-04,   0.26013E-04,   0.25091E-04, &
                    0.22587E-04,   0.22693E-04,   0.13498E-04,   0.17758E-04,   0.16565E-04,   0.14991E-04,   0.14364E-04, &
                    0.12701E-04,   0.11813E-04,   0.10163E-04,   0.86087E-05,   0.75787E-05,   0.86740E-05,   0.80798E-05, &
                    0.87940E-05,   0.68895E-05,   0.55891E-05,   0.56224E-05,   0.49614E-05,   0.41010E-05,   0.37850E-05, &
                    0.30495E-05,   0.24522E-05 /)
      else
        write(6,*) "ERROR: input file for desired bands and frqs is not available"
        stop
      endif
 
      if (Solar_spect%tot_wvnums /=       &
                                 endwvnbands(Solar_spect%nbands)) then
        write(6,*) "ERROR: inconsistency between highest solar spectrum wavenumber"
        write(6,*) "       in esfsw_parameters_mod and in esfsw_sriver input file"
        stop
      endif

!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.34 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_340) then
          Solar_spect%w340_band_indx = ni
          Solar_spect%w340_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.38 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_380) then
          Solar_spect%w380_band_indx = ni
          Solar_spect%w380_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to blue light
!    (0.44 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_440) then
          Solar_spect%w440_band_indx = ni
          Solar_spect%w440_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to red light
!    (0.67 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_670) then
          Solar_spect%w670_band_indx = ni
          Solar_spect%w670_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to visible light
!    (0.55 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > vis_wvnum) then
          Solar_spect%visible_band_indx = ni
          Solar_spect%visible_band_indx_iz = .true.
          exit
        endif
      end do

!---------------------------------------------------------------------
!    define the band index corresponding to 870nm
!    (0.87 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_870) then
          Solar_spect%eight70_band_indx = ni
          Solar_spect%eight70_band_indx_iz = .true.
          exit
        endif
      end do
      
!---------------------------------------------------------------------
!    define the band index corresponding to near infra red band
!    (1.00 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > one_micron_wvnum) then
          Solar_spect%one_micron_indx = ni
          Solar_spect%one_micron_indx_iz = .true.
          exit
        endif
      end do
 
!---------------------------------------------------------------------
!    define the band index corresponding to               
!    (1.61 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > onepsix_micron_wvnum) then
          onepsix_band_indx = ni
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the wavenumber one solar fluxes.                       
!----------------------------------------------------------------- --
      nw1 = 1
      nw2 = nw1 + nwvnsolar(1) - 1
      do nw = nw1,nw2
        Solar_spect%solarfluxtoa(nw) = solint(1)
      end do
      do ni = 2,nintsolar
        nw1 = nw1 + nwvnsolar(ni-1)
        nw2 = nw1 + nwvnsolar(ni) - 1
        do nw = nw1,nw2
          Solar_spect%solarfluxtoa(nw) = solint(ni)
        end do
      end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      deallocate  (solint    )
      deallocate  (nwvnsolar )
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if ( .not. Rad_control%using_solar_timeseries_data) then
        Solar_spect%solflxband = solflxbandref
      endif
      Solar_spect%solflxbandref = solflxbandref
      Solar_spect%endwvnbands = endwvnbands

!--------------------------------------------------------------------
!    override the input file value of firstrayband if the nml control
!    variables indicates rayleigh effects are to be considered in
!    all bands.
!--------------------------------------------------------------------
      if (do_rayleigh_all_bands)  firstrayband = 1

!----------------------------------------------------------------------
!    convert some values to mks to match model units
!--------------------------------------------------------------------
      p0h2o(1:NH2OBANDS) = 1.0E-2/p0h2o(1:NH2OBANDS)  ! invert, and convert mb to mks
      kh2o(1:nfrqpts) = kh2o(1:nfrqpts) *1.0E-01   ! cgs to mks
      ko3(1:nfrqpts)  = ko3(1:nfrqpts) *1.0E-01    ! cgs to mks

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      do n=1,NH2OBANDS
        if (c1co2(n) /= input_flag .and.  &
            c2co2(n) /= input_flag .and.  &
            c3co2(n) /= input_flag ) then
          c4co2(n) = c1co2(n) * c2co2(n) ** c3co2(n)
          c4co2str(n) = c1co2str(n) * c2co2str(n) ** c3co2str(n)
          totco2max(n) = ( (1.0/c1co2(n) ) + c2co2(n) ** c3co2(n) ) ** &
                       (1.0/c3co2(n) ) - c2co2(n)
          if (nbands == 18) then
            if ( n /= 4) then
              totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
            else 
              totco2strmax(n) = HUGE (c4o2strschrun) 
            endif
          else
              totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
          endif
        else
          c4co2(n) = 0.0                              
          c4co2str(n) = 0.0
          totco2max(n) = 0.0                                            
          totco2strmax(n) = 0.0
        endif
        if (c1o2(n) /= input_flag .and.   &
            c2o2(n) /= input_flag .and.   &
            c3o2(n) /= input_flag ) then
          c4o2(n) = c1o2(n) * c2o2(n) ** c3o2(n)
          c4o2str(n) = c1o2str(n) * c2o2str(n) ** c3o2str(n)
          toto2max(n) = ( (1.0/c1o2(n) ) + c2o2(n) ** c3o2(n) ) ** &
                          (1.0/c3o2(n) ) - c2o2(n)
          if (nbands == 18) then
            if ( n /= 4) then
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                c2o2str(n)
            else
              toto2strmax(n) = HUGE (c4o2strschrun) 
            endif
          else
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                c2o2str(n)
          endif
        else
          c4o2(n) = 0.0                              
          c4o2str(n) = 0.0
          toto2max(n) = 0.0                                            
          toto2strmax(n) = 0.0
        endif
    if (do_ch4_sw_effects) then
        if (c1ch4(n) /= input_flag .and.  &
            c2ch4(n) /= input_flag .and.  &
            c3ch4(n) /= input_flag ) then
          c4ch4(n) = c1ch4(n) * c2ch4(n) ** c3ch4(n)
          c4ch4str(n) = c1ch4str(n) * c2ch4str(n) ** c3ch4str(n)
          totch4max(n) = ( (1.0/c1ch4(n) ) + c2ch4(n) ** c3ch4(n) ) ** &
                           (1.0/c3ch4(n) ) - c2ch4(n)
          totch4strmax(n) = ( (1.0/c1ch4str(n) ) + c2ch4str(n) ** &
                            c3ch4str(n) ) ** (1.0/c3ch4str(n) ) - &
                            c2ch4str(n)
        else
          c4ch4(n) = 0.
          c4ch4str(n) = 0.
          totch4max(n) = 0.
          totch4strmax(n) = 0.     
        endif
     endif
     if (do_n2o_sw_effects) then
        if (c1n2o(n) /= input_flag .and.  &
            c2n2o(n) /= input_flag .and.  &
            c3n2o(n) /= input_flag ) then
          c4n2o(n) = c1n2o(n) * c2n2o(n) ** c3n2o(n)
          c4n2ostr(n) = c1n2ostr(n) * c2n2ostr(n) ** c3n2ostr(n)
          totn2omax(n) = ( (1.0/c1n2o(n) ) + c2n2o(n) ** c3n2o(n) ) ** &
                           (1.0/c3n2o(n) ) - c2n2o(n)
          totn2ostrmax(n) = ( (1.0/c1n2ostr(n) ) + c2n2ostr(n) ** &
                               c3n2ostr(n) ) ** (1.0/c3n2ostr(n) ) - &
                               c2n2ostr(n)
        else
          c4n2o(n) = 0.
          c4n2ostr(n) = 0.
          totn2omax(n) = 0.
          totn2ostrmax(n) = 0.     
        endif
     endif
      end do

      c4o2strschrun = c1o2strschrun * c2o2strschrun ** c3o2strschrun
      toto2strmaxschrun = ( (1.0/c1o2strschrun) + c2o2strschrun ** &
                             c3o2strschrun) ** (1.0/c3o2strschrun) - &
                             c2o2strschrun

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if ( .not. Rad_control%using_solar_timeseries_data) then
      solflxtotal = 0.0
      do nband = 1,NBANDS
        solflxtotal = solflxtotal + Solar_spect%solflxbandref(nband)
      end do
      endif
 
!----------------------------------------------------------------------
!    define the wavenumbers to evaluate rayleigh optical depth.      
!-------------------------------------------------------------------
      endwvnbands(0) = 0
      do nband = 1,NBANDS
        freqnu(nband) = 0.5 * ( endwvnbands(nband-1) +   &
                                endwvnbands(nband) )
      end do
 
!---------------------------------------------------------------------
!    define quantities used to determine the rayleigh optical depth. 
!    notes: refquanray is the quantity which multiplies pressure /  
!           temperature to yield the molecular density.                 
!           betaddensitymol is the quantity which multiples the       
!           molecular density to yield the rayleigh scattering      
!           coefficient.                                           
!           1.39E-02 is the depolorization factor.              
!-----------------------------------------------------------------
      refquanray = densmolref * temprefray / pressrefray 
      corrfac = ( 6.0E+00 + 3.0E+00 * 1.39E-02 )/( 6.0E+00 - 7.0E+00 * &
                  1.39E-02 )
      gamma = 1.39E-02 / ( 2.0E+00 - 1.39E-02 )
      f1 = 7.5E-01 / ( gamma * 2.0E+00 + 1.0E+00 )
      f2 = gamma * 3.0E+00 + 1.0E+00 
      f3 = 1.0E+00 - gamma 
      pinteg = 2.0E+00 * PI * ( 2.0E+00 * f1 * f2 * ( 1.0E+00 + f3 / &
               f2 / 3.0E+00 ) )
      twopiesq = 2.0E+00 *  PI ** 2
      densmolrefsqt3 = 3.0E+00 * densmolref ** 2
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do nband = 1,NBANDS
        wavelength = 1.0E+04 / freqnu(nband)
        freqsq = 1.0 / ( wavelength ) ** 2
        ristdm1 = ( 6.4328E+03 + 2.94981E+06 / ( 1.46E+02 - freqsq ) + &
                    2.554E+04 / ( 4.1E+01 - freqsq ) ) * 1.0E-08
        ri = ristdm1 + 1
        betaddensitymol(nband) = twopiesq*( ri ** 2 - 1.0E+00 ) ** 2 * &
                                 corrfac / ( densmolrefsqt3 *  &
                                 wavelength ** 4 ) * pinteg * convfac
      end do
 
      gausswt(1) = 1.0
      
!---------------------------------------------------------------------
!    define the gaussian angles for evaluation of the diffuse beam.  
!--------------------------------------------------------------------
      do i = 1,nstreams
        cosangstr(i) = ( ptstr(i) + 1. ) * 5.0E-01
      end do

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.
!--------------------------------------------------------------------
! 101  format( 12f10.4 )
! 102  format( 32i4 )
! 103  format( 20i6 )
! 104  format( 12f10.2 )
! 105  format( 1p,16e8.1 )
! 106  format( 1p,3e16.6,l16 )
! 107  format( i5,1p,e14.5 )
 
!---------------------------------------------------------------------


end subroutine esfsw_driver_init 



!#################################################################
! <SUBROUTINE NAME="swresf">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call swresf(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,  &
                   Aerosol, Aerosol_props, Astro, Cldrad_props,      &
                   Cld_spec, including_volcanoes, Sw_output,         &
                   Aerosol_diags, r, including_aerosols,             &
                   naerosol_optical) 

!----------------------------------------------------------------------
!    swresf uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

integer,                       intent(in)    :: is, ie, js, je
type(atmos_input_type),        intent(in)    :: Atmos_input
type(surface_type),            intent(in)    :: Surface
type(radiative_gases_type),    intent(in)    :: Rad_gases   
type(aerosol_type),            intent(in)    :: Aerosol     
type(aerosol_properties_type), intent(in)    :: Aerosol_props
type(astronomy_type),          intent(in)    :: Astro
type(cldrad_properties_type),  intent(in)    :: Cldrad_props
type(cld_specification_type),  intent(in)    :: Cld_spec      
logical,                       intent(in)    :: including_volcanoes
type(sw_output_type),          intent(inout) :: Sw_output   
type(aerosol_diagnostics_type),intent(inout) :: Aerosol_diags
real,dimension(:,:,:,:),       intent(inout) :: r
logical,                       intent(in)    :: including_aerosols  
integer,                       intent(in)    :: naerosol_optical


!-------------------------------------------------------------------
!  intent(in) variables:
!
!    is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Surface        surface_type structure, contains variables 
!                     defining the surface characteristics
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Aerosol_props  aerosol radiative property input data for the 
!                     radiation package
!      Astro          astronomy_type structure
!      Cldrad_props   cloud radiative property input data for the 
!                     radiation package
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 
      logical, dimension (size(Atmos_input%temp,1), &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)-1)  :: &
                             cloud

      logical, dimension (size(Atmos_input%temp,1), &
                          size(Atmos_input%temp,2))  ::   &
                              cloud_in_column,   daylight

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1, &
                        Solar_spect%nbands)  :: &
            aeroasymfac,     aerosctopdep,   aeroextopdep, &
            rayopdep

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1, &
                       nfrqpts,NSOLWG)  ::     &
                                                    gasopdep

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1)  :: &
            cloudfrac,      cldfrac_band,    cldfrac,         &
            deltaz,         cloud_deltaz,                     &
            cldext,         cldsct,          cldasymm,         &
            cloudasymfac,   cloudextopdep,                     &
            cloudsctopdep,                   deltap,           &
            densitymol,     &
            fclr,            fovc,             &
            gocpdp,         gstrclr,          &
            gstrovc,        omegastrclr,      &
            omegastrovc,    rlayerdif,       rlayerdir,        &
            rlayerdifclr,   rlayerdifovc,    &
            rlayerdirclr,   rlayerdirovc,                     &
            sctopdepclr,    sctopdepovc,      &
            taustrclr,      taustrovc,        &
            tlayerde,       tlayerdeclr,      &
            tlayerdeovc,    tlayerdif,                        &
            tlayerdifclr,   tlayerdifovc,     &
            tlayerdir,      tlayerdirclr,     &
            tlayerdirovc

      real :: gclr, govc, extopdepclr, ssalbclr, extopdepovc
      real :: ssalbovc
           
      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3))  :: &
              reflectance,  transmittance,  tr_dir            

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3), &
                       Rad_control%nzens)  :: &
            dfswbandclr,     fswbandclr,    ufswbandclr, &
            dfswband,        fswband,       ufswband,  &
            sumtrclr,        sumreclr,      sumtr_dir_clr,  & 
            sumtr,           sumre,         sumtr_dir

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       Rad_control%nzens)  :: sumtr_dir_up

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3))  :: &
            press,           pflux,            pflux_mks, &
            temp,                                       &
            reflectanceclr,  transmittanceclr, tr_dirclr

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1 , &
                       Rad_control%nzens)  :: &
                hswbandclr

      real :: hswband

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2))    ::  &
            sfcalb_dir,   sfcalb_dif,  wtfac_p,    &
            fracday_p,    solarflux_p

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),   &
                                         NSOLWG) ::         &
                 cosangsolar_p

      integer :: j, i, k, ng, np, nband, nf, ns, nz
      integer :: nzens

      integer :: nprofile, nprofiles
      real    :: profiles_inverse

      integer :: ix, jx, kx, israd, jsrad, ierad, jerad, ksrad, kerad
      real    :: ssolar  
      real    :: solflxtotal_local

!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        write(6,*) "ERROR: esfsw_driver_mode has not been initialized"
        stop
      endif

!---------------------------------------------------------------------
!    define the solar_constant appropriate at Rad_time when solar
!    input is varying with time.
!---------------------------------------------------------------------
      if (Rad_control%using_solar_timeseries_data) then
        solflxtotal_local = 0.0
        do nband = 1,NBANDS
          solflxtotal_local = solflxtotal_local + Solar_spect%solflxband(nband)
        end do
      else
        solflxtotal_local = solflxtotal
      endif

!---------------------------------------------------------------------
!
      fswband = 0. ! MJT suggestion
      cldfrac = Cld_spec%camtsw
      
!---------------------------------------------------------------------
!  convert to cgs and then back to mks for consistency with previous 
!---------------------------------------------------------------------
      press(:,:,:) = 0.1*(10.0*Atmos_input%press(:,:,:))
      pflux(:,:,:) =     (10.0*Atmos_input%pflux(:,:,:))
      deltaz(:,:,:) = Atmos_input%deltaz(:,:,:)
      temp  (:,:,:) = Atmos_input%temp  (:,:,:)
      cloud_deltaz = Atmos_input%clouddeltaz

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      ix = size(temp,1)
      jx = size(temp,2)
      kx = size(temp,3) - 1
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = ix
      jerad = jx
      kerad = kx

!----------------------------------------------------------------------c
!    define a flag indicating columns in which there is sunshine during
!    this radiation time step. define a flag indicating points with both
!    sunlight and cloud.      
!----------------------------------------------------------------------c
      daylight = ( Astro%fracday /= 0.0 )
      cloud_in_column = .false.

      call compute_aerosol_optical_props     &
                (Atmos_input, Aerosol, Aerosol_props, &
                 including_volcanoes, Aerosol_diags, r,  &
                 including_aerosols, naerosol_optical, &
                 daylight, aeroextopdep, aerosctopdep, aeroasymfac) 

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ssolar = Sw_control%solar_constant*Astro%rrsun
      
!----------------------------------------------------------------------
!    define a flag indicating points with both sunlight and cloud. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
      if (.not. Cldrad_control%do_stochastic_clouds) then
        do k = KSRAD,KERAD
          cloud(:,:,k) = daylight.and.(cldfrac(:,:,k) > 0.0)
        end do
        cloud_in_column = any(cloud,3)
        where (cloud)
          cloudfrac = cldfrac
        elsewhere
          cloudfrac = 0.
        end where
      endif


!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
      pflux_mks = pflux*1.0E-1

      do k = KSRAD+1,KERAD+1
        deltap(:,:,k-1) = pflux_mks(:,:,k) - pflux_mks(:,:,k-1)
        gocpdp(:,:,k-1) = radcon_mks/deltap(:,:,k-1)
      end do
 
      call compute_gas_props (Atmos_input, Rad_gases, Astro,   &
                              daylight, gasopdep)
 
!---------------------------------------------------------------------
!    define the molecular density for use in calculating the           
!    rayleigh optical depth (deltaz is in meters).                     
!--------------------------------------------------------------------
      do k = KSRAD,KERAD
        densitymol(:,:,k) = refquanray * press(:,:,k) / temp(:,:,k)
      end do
 
! assumption is that there is 1 cloud profile for each sw band
      if (do_ica_calcs) then
        nprofiles = nbands
        profiles_inverse = 1.0/nprofiles
      else
        nprofiles = 1
        profiles_inverse = 1.0
      endif

!--------------------------------------------------------------------
!    define the rayleigh optical depths.                                
!---------------------------------------------------------------------
      do nband = 1, NBANDS
        rayopdep(:,:,:,nband) = betaddensitymol(nband)*  &
                                        densitymol(:,:,:)*deltaz(:,:,:)
      end do   ! (nband loop)

      nzens = Rad_control%nzens

!--------------------------------------------------------------------
 
      do nprofile=1, nprofiles
        if (do_ica_calcs) then
          cldfrac_band(:,:,:) = Cld_spec%camtsw_band(:,:,:,nprofile)
        endif

!----------------------------------------------------------------------c
!    np is a counter for the pseudo-monochromatic frequency point 
!    number.
!----------------------------------------------------------------------c
        np = 0

        if (Rad_control%do_totcld_forcing) then 
          dfswbandclr(:,:,:,1:nzens) = 0.0
          fswbandclr(:,:,:,1:nzens) = 0.0
          hswbandclr(:,:,:,1:nzens) = 0.0
          ufswbandclr(:,:,:,1:nzens) = 0.0
        endif
        reflectanceclr = 0.0
        transmittanceclr = 0.0
 
!----------------------------------------------------------------------c
!    begin band loop                                                   
!----------------------------------------------------------------------c
        do nband = 1,NBANDS
 
          sumtr(:,:,:,1:nzens) = 0.0
          sumtr_dir(:,:,:,1:nzens) = 0.0
          sumtr_dir_up(:,:,1:nzens) = 0.0
          sumre(:,:,:,1:nzens) = 0.0
          if (Rad_control%do_totcld_forcing) then
            sumtrclr(:,:,:,1:nzens) = 0.0
            sumreclr(:,:,:,1:nzens) = 0.0
            sumtr_dir_clr(:,:,:,1:nzens) = 0.0
          endif
          if (Cldrad_control%do_stochastic_clouds) then
            if (.not. do_ica_calcs) then
              cldfrac_band(:,:,:) = Cld_spec%camtsw_band(:,:,:,nband)
            endif
 
!----------------------------------------------------------------------
!    if stochastic clouds are activated (cloud fields differ with sw
!    parameterization band), define a flag indicating points with both
!    sunlight and cloud for the current sw parameterization band. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                cloud_in_column(i,j) = .false.
                if (daylight(i,j)) then
                  do k = KSRAD,KERAD
                    if (cldfrac_band(i,j,k) > 0.0)  then
                      cloud_in_column(i,j) = .true.
                      cloud(i,j,k) = .true.
                      cloudfrac(i,j,k) = cldfrac_band(i,j,k)
                    else
                      cloud(i,j,k) = .false.
                      cloudfrac(i,j,k) = 0.0
                    endif
                  end do
                else
                  do k = KSRAD,KERAD
                    cloud(i,j,k) = .false.
                    cloudfrac(i,j,k) = 0.0
                  end do
                endif
              end do
            end do
          endif

!---------------------------------------------------------------------
!    obtain cloud properties from the Cldrad_props input variable.
!--------------------------------------------------------------------
          cldext(:,:,:) = Cldrad_props%cldext(:,:,:,nband,nprofile)
          cldsct(:,:,:) = Cldrad_props%cldsct(:,:,:,nband,nprofile)
          cldasymm(:,:,:) = Cldrad_props%cldasymm(:,:,:,nband,nprofile)

          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                  if ( cloud(i,j,k) ) then
                    cloudextopdep(i,j,k) = 1.0E-03*cldext(i,j,k) * cloud_deltaz(i,j,k)
                    cloudsctopdep(i,j,k) = 1.0E-03*cldsct(i,j,k) * cloud_deltaz(i,j,k)
                    cloudasymfac(i,j,k) = cldasymm(i,j,k)
                  end if
              end do
            end do
          end do

!-----------------------------------------------------------------
!    define clear sky arrays
!-----------------------------------------------------------------
          if (nband >= firstrayband ) then
            do k=ksrad,kerad
              do j = JSRAD,JERAD
                do i = ISRAD,IERAD
                  if ( daylight(i,j) ) then
                    sctopdepclr(i,j,k) = rayopdep(i,j,k,nband) + aerosctopdep(i,j,k,nband)
                    gclr = aeroasymfac(i,j,k,nband)*aerosctopdep(i,j,k,nband)/sctopdepclr(i,j,k)
                    fclr(i,j,k) = aeroasymfac(i,j,k,nband)*gclr
                    gstrclr(i,j,k) = ( gclr  - fclr(i,j,k) )/( 1.0 - fclr(i,j,k) )
                  end if
                end do
              end do
            end do
          endif

!-----------------------------------------------------------------
!    define cloudy sky arrays
!-----------------------------------------------------------------
          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                if (cloud(i,j,k)) then
                  sctopdepovc(i,j,k) = rayopdep(i,j,k,nband) +     &
                                   aerosctopdep(i,j,k,nband) +     &
                                        cloudsctopdep(i,j,k) 
                  govc = ( ( cloudasymfac(i,j,k) *          &
                              cloudsctopdep(i,j,k) ) +      &
                         ( aeroasymfac(i,j,k,nband) *       &
                         aerosctopdep(i,j,k,nband)))/       &
                                          sctopdepovc(i,j,k)
                  fovc(i,j,k) = ( ( cloudasymfac(i,j,k) ** 2 *     &
                                       cloudsctopdep(i,j,k) ) +    &
                             ( aeroasymfac(i,j,k,nband) ** 2 *     &
                                 aerosctopdep(i,j,k,nband) ))/     &
                                          sctopdepovc(i,j,k)
                  gstrovc(i,j,k) = ( govc - fovc(i,j,k))/  &
                                         ( 1.0 - fovc(i,j,k) )
                end if
              end do
            end do
          end do
          
!---------------------------------------------------------------------
!    begin frequency points in the band loop.                          
!--------------------------------------------------------------------
          do nf = 1,nfreqpts(nband)
            np = np + 1
 
!---------------------------------------------------------------------
!    begin gaussian angle loop (ng > 1 only when lswg = true).        
!--------------------------------------------------------------------
            do ng = 1,NSOLWG
 
!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------
              if (nband >= firstrayband )  then
                do k=ksrad,kerad
                  do j = JSRAD,JERAD
                    do i = ISRAD,IERAD
                      if ( daylight(i,j) ) then
                        extopdepclr = gasopdep(i,j,k,np,ng) +   &
                                     rayopdep(i,j,k,nband) +    &
                                     aeroextopdep(i,j,k,nband)
                        ssalbclr = sctopdepclr(i,j,k)/          &
                                    extopdepclr
                        taustrclr(i,j,k) = extopdepclr*( 1.0 -         &
                                        ssalbclr*fclr(i,j,k) )
                        omegastrclr(i,j,k) =                           &
                               ssalbclr*((1.0 - fclr(i,j,k))/          &
                               (1.0 -  ssalbclr*fclr(i,j,k)))
                      end if
                    end do
                  end do
                end do
              endif

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                         
!--------------------------------------------------------------------
              do k = KSRAD,KERAD
                do j = JSRAD,JERAD
                  do i = ISRAD,IERAD
                    if ( cloud(i,j,k) ) then
                      extopdepovc = gasopdep(i,j,k,np,ng) +    &
                                      rayopdep(i,j,k,nband) +  &
                                   aeroextopdep(i,j,k,nband) + &
                                          cloudextopdep(i,j,k)
                      ssalbovc = sctopdepovc(i,j,k) /    &
                                         extopdepovc
                      taustrovc(i,j,k) = extopdepovc*( 1.0 - &
                                      ssalbovc*fovc(i,j,k) )
                      omegastrovc(i,j,k) = ssalbovc*( ( 1.0 - &
                                             fovc(i,j,k) )/( 1.0 -   &
                                                 ssalbovc *   &
                                                 fovc(i,j,k) ) )
                    end if
                  end do
                end do
              end do

!---------------------------------------------------------------------
!    do calculations for all desired zenith angles.
!---------------------------------------------------------------------
              do nz = 1, nzens
                if (Rad_control%hires_coszen) then
                  cosangsolar_p(:,:,ng) = Astro%cosz_p(:,:,nz)   
                else
                  cosangsolar_p(:,:,ng) = Astro%cosz(:,:)          
                endif

                where (cosangsolar_p(:,:,:) == 0.0)   &
                                            cosangsolar_p(:,:,:) = 1.0

!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                       
!---------------------------------------------------------------------
                if (nband >= firstrayband )  then

!---------------------------------------------------------------------
!    do diffuse calculation only for first zenith angle -- it is 
!    independent of zenith angle.
!---------------------------------------------------------------------
                  if (nz == 1) then
                    call deledd     &
                        (ix, jx, kx, taustrclr, omegastrclr, &
                         gstrclr, cosangsolar_p(:,:,ng), ng, daylight, &
                         rlayerdirclr, tlayerdirclr, tlayerdeclr, &
                         rlayerdif=rlayerdifclr, tlayerdif=tlayerdifclr)
                  else
                    call deledd     &
                        (ix, jx, kx, taustrclr, omegastrclr, &
                         gstrclr, cosangsolar_p(:,:,ng), ng, daylight,&
                         rlayerdirclr, tlayerdirclr, tlayerdeclr)
                  endif

!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1.
!---------------------------------------------------------------------
                  if (ng /= 1) then
                    tlayerdifclr = 0.0       
                    if (NSTREAMS == 1) then
                      tlayerdifclr(:,:,:) =   &
                             exp(-gasopdep(:,:,:,np,ng)/cosangstr(1))
                    else
                      do ns = 1,NSTREAMS
                        tlayerdifclr(:,:,:) =    &
                             tlayerdifclr(:,:,:) +  & 
                                   exp( -gasopdep(:,:,:,np,ng)/&
                                         cosangstr(ns) )*wtstr(ns)* &
                                         cosangstr(ns)
                      end do
                    endif
                    rlayerdifclr = 0.0
                  endif
  
!---------------------------------------------------------------------
!    initialize the layer reflection and transmission arrays for the   
!    non-rayleigh-scattering case.                            
!-------------------------------------------------------------------
                else
!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1, and all bands see rayleigh
!    scattering.
!---------------------------------------------------------------------
                  tlayerdifclr = 0.0       
                  if (NSTREAMS == 1) then
                    tlayerdifclr(:,:,:) =     &
                            exp( -gasopdep(:,:,:,np,ng)/cosangstr(1))
                  else
                    do ns = 1,NSTREAMS
                      tlayerdifclr(:,:,:) =   &
                         tlayerdifclr(:,:,:) +    &
                              exp(-gasopdep(:,:,:,np,ng)/&
                                cosangstr(ns))*wtstr(ns)*cosangstr(ns)
                    end do
                  endif
                  rlayerdirclr(:,:,:) = 0.0
                  do k=KSRAD,KERAD
                    tlayerdirclr(:,:,k) =      &
                                  exp( -gasopdep(:,:,k,np,ng) /   &
                                                 cosangsolar_p(:,:,ng) )
                  end do  
                  tlayerdeclr(:,:,:) = tlayerdirclr(:,:,:)
                  rlayerdifclr = 0.0
                endif

!---------------------------------------------------------------------
!    overcast sky mode                                                  
!    note: in this mode, the delta-eddington method is performed only 
!    for spatial columns containing a cloud and experiencing sunlight. 
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!    calculate the reflection and transmission in the scattering layers 
!    using the delta-eddington method.                                  
!-------------------------------------------------------------------
                if (nz == 1) then
                  call deledd      &
                       (ix, jx, kx, taustrovc, omegastrovc, gstrovc, &
                        cosangsolar_p(:,:,ng), ng, daylight, &
                        rlayerdirovc, tlayerdirovc, tlayerdeovc, &
                        rlayerdif=rlayerdifovc, tlayerdif=tlayerdifovc,&
                        cloud=cloud)
                else
                  call deledd      &
                       (ix, jx, kx, taustrovc, omegastrovc, gstrovc, &
                        cosangsolar_p(:,:,ng), ng, daylight, &
                        rlayerdirovc, tlayerdirovc, tlayerdeovc,   & 
                        cloud=cloud)
                endif
                if (ng /= 1) then
                  tlayerdifovc(:,:,:) = tlayerdifclr(:,:,:)
                  rlayerdifovc(:,:,:) = rlayerdifclr(:,:,:)
                endif
 
!-------------------------------------------------------------------- 
!    weight the reflection and transmission arrays for clear and        
!    overcast sky conditions by the cloud fraction, to calculate the    
!    resultant values.                                                  
!---------------------------------------------------------------------- 
                do k=KSRAD,KERAD
                  do j = JSRAD,JERAD
                    do i = ISRAD,IERAD
                      if ( cloud(i,j,k) ) then
                        rlayerdir(i,j,k) = cloudfrac(i,j,k)*   &
                                           rlayerdirovc(i,j,k) +  &
                                           (1.0 - cloudfrac(i,j,k) )*  &
                                           rlayerdirclr(i,j,k)
                        rlayerdif(i,j,k) = cloudfrac(i,j,k) *  &
                                           rlayerdifovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           rlayerdifclr(i,j,k)
                        tlayerdir(i,j,k) = cloudfrac(i,j,k) *   &
                                           tlayerdirovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdirclr(i,j,k)
                        tlayerdif(i,j,k) = cloudfrac(i,j,k) *   &
                                          tlayerdifovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdifclr(i,j,k)
                        tlayerde(i,j,k) =  cloudfrac(i,j,k) *   &
                                           tlayerdeovc(i,j,k) +  &
                                           (1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdeclr(i,j,k)
                      else if (daylight(i,j)) then
                        rlayerdir(i,j,k) = rlayerdirclr(i,j,k)
                        tlayerdir(i,j,k) = tlayerdirclr(i,j,k)
                        rlayerdif(i,j,k) = rlayerdifclr(i,j,k)
                        tlayerdif(i,j,k) = tlayerdifclr(i,j,k)
                        tlayerde (i,j,k) = tlayerdeclr (i,j,k)
                      else
                        rlayerdir(i,j,k) = 0. ! MJT suggestion
                        tlayerdir(i,j,k) = 0. ! MJT suggestion
                        rlayerdif(i,j,k) = 0. ! MJT suggestion
                        tlayerdif(i,j,k) = 0. ! MJT suggestion
                        tlayerde (i,j,k) = 0. ! MJT suggestion
                      end if
                    end do
                  end do
                end do
 
!---------------------------------------------------------------------
!    define the surface albedo (infrared value for infrared bands,      
!    visible value for the remaining bands).                            
!----------------------------------------------------------------------c
                if (nband <= NIRBANDS ) then
                  sfcalb_dir(:,:) = Surface%asfc_nir_dir(:,:)
                  sfcalb_dif(:,:) = Surface%asfc_nir_dif(:,:)
                else
                  sfcalb_dir(:,:) = Surface%asfc_vis_dir(:,:)
                  sfcalb_dif(:,:) = Surface%asfc_vis_dif(:,:)
                end if
 
!-------------------------------------------------------------------- 
!    calculate the reflection and transmission at flux levels from the  
!    direct and diffuse values of reflection and transmission in the  
!    corresponding layers using the adding method.                      
!---------------------------------------------------------------------
                call adding         &
                    (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif, &
                     tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,    &
                     daylight, reflectance, transmittance, tr_dir)    

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                if (Rad_control%do_totcld_forcing) then
                  call adding       &
                      (ix, jx,  kx, rlayerdirclr, tlayerdirclr,   &
                       rlayerdifclr, tlayerdifclr, tlayerdeclr,   &
                       sfcalb_dir,  sfcalb_dif, cloud_in_column,  &
                       reflectanceclr, transmittanceclr, tr_dirclr)
                endif

!---------------------------------------------------------------------- 
!    weight and sum the reflectance and transmittance to calculate the 
!    band values.                                                     
!-------------------------------------------------------------------
                wtfac_p(:,:) = wtfreq(np)*gausswt(ng)*cosangsolar_p(:,:,ng)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                do k = KSRAD,KERAD+1
                  do j = JSRAD,JERAD
                    do i = ISRAD,IERAD
                      if ( daylight(i,j) ) then
                        sumtr(i,j,k,nz) = sumtr(i,j,k,nz) +    &
                                      transmittance(i,j,k)*wtfac_p(i,j)
                        sumtr_dir(i,j,k,nz) = sumtr_dir(i,j,k,nz) +  &
                                        tr_dir(i,j,k)*wtfac_p(i,j)
                        sumre(i,j,k,nz) = sumre(i,j,k,nz) +     &
                                       reflectance(i,j,k)*wtfac_p(i,j)
                      end if
                    end do
                  end do
                end do
                where ( daylight )
                  sumtr_dir_up(:,:,nz) = sumtr_dir_up(:,:,nz) + &
                    tr_dir(:,:,KERAD+1)*sfcalb_dir(:,:)*wtfac_p(:,:)
                end where

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                if (Rad_control%do_totcld_forcing) then
                  do k = KSRAD,KERAD+1
                    do j = JSRAD,JERAD
                      do i = ISRAD,IERAD
                        if (cloud_in_column(i,j)) then
                          sumtrclr(i,j,k,nz) =    &
                                sumtrclr(i,j,k,nz) +   &
                                transmittanceclr(i,j,k)* wtfac_p(i,j) 
                          sumtr_dir_clr(i,j,k,nz) = &
                                sumtr_dir_clr(i,j,k,nz) +  &
                                tr_dirclr(i,j,k)*wtfac_p(i,j)
                          sumreclr(i,j,k,nz) = sumreclr(i,j,k,nz) +   &
                                reflectanceclr(i,j,k)*wtfac_p(i,j)
                        else if ( daylight(i,j) ) then
                          sumtrclr(i,j,k,nz) = sumtrclr(i,j,k,nz) +   &
                                transmittance(i,j,k)*wtfac_p(i,j)
                          sumtr_dir_clr(i,j,k,nz) =    &
                             sumtr_dir_clr(i,j,k,nz) + tr_dir(i,j,k)*&
                                                        wtfac_p(i,j)
                          sumreclr(i,j,k,nz) = sumreclr(i,j,k,nz) +   &
                                        reflectance(i,j,k)*wtfac_p(i,j)
                        end if
                      end do
                    end do
                  end do
                endif
              end do  ! end of nz loop
            end do    ! end of gaussian loop
          end do      ! end of frequency points in the band loop
 
!----------------------------------------------------------------------
!    normalize the solar flux in the band to the appropriate value for  
!    the given total solar insolation.                                 
!---------------------------------------------------------------------
          do nz = 1,nzens
            if (Rad_control%hires_coszen) then
              fracday_p(:,:) = Astro%fracday_p(:,:,nz)
            else
              fracday_p(:,:) = Astro%fracday(:,:)
            endif
            solarflux_p(:,:) = fracday_p(:,:)*  &
                                   Solar_spect%solflxband(nband)*  &
                                                      ssolar/solflxtotal_local
 
            if (nband == Solar_spect%visible_band_indx) then
              Sw_output%bdy_flx(:,:,1,nz) =   &
                  Sw_output%bdy_flx(:,:,1,nz) + sumre(:,:,1,nz)*   &
                                                        solarflux_p(:,:)
              Sw_output%bdy_flx(:,:,3,nz) =    &
                  Sw_output%bdy_flx(:,:,3,nz) + sumtr(:,:,KERAD+1,nz)*&
                                                solarflux_p(:,:) -  &
                                                sumre(:,:,KERAD+1,nz)* &
                                                solarflux_p(:,:) 
            endif
            if (nband == onepsix_band_indx) then
               Sw_output%bdy_flx(:,:,2,nz) =     &
                   Sw_output%bdy_flx(:,:,2,nz) + sumre(:,:,1,nz)*  &
                                                        solarflux_p(:,:)
               Sw_output%bdy_flx(:,:,4,nz) =    &
                   Sw_output%bdy_flx(:,:,4,nz) + sumtr(:,:,KERAD+1,nz)*&
                                                 solarflux_p(:,:) - &
                                                 sumre(:,:,KERAD+1,nz)*&
                                                 solarflux_p(:,:)
            endif
          
!-------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!--------------------------------------------------------------------
            if (do_esfsw_band_diagnostics) then
              do k = KSRAD,KERAD+1
                dfswband(:,:,k,nz) = sumtr(:,:,k,nz)*solarflux_p(:,:) 
                ufswband(:,:,k,nz) = sumre(:,:,k,nz)*solarflux_p(:,:)
              end do
            endif
 
!----------------------------------------------------------------------
!    sum the band fluxes and heating rates to calculate the total       
!    spectral values.                                                  
!------------------------------------------------------------------
            do k = KSRAD,KERAD+1
              do j = JSRAD,JERAD
                do i = ISRAD,IERAD
                  if ( daylight(i,j) ) then
                    Sw_output%dfsw (i,j,k,nz) =   &
                       Sw_output%dfsw(i,j,k,nz) + sumtr(i,j,k,nz)*&
                                                      solarflux_p(i,j)
                    Sw_output%ufsw (i,j,k,nz) =   &
                       Sw_output%ufsw(i,j,k,nz) + sumre(i,j,k,nz)*  &
                                                       solarflux_p(i,j)
                    fswband(i,j,k,nz) = ((sumre(i,j,k,nz)*  &
                                          solarflux_p(i,j)) - &
                                         (sumtr(i,j,k,nz)*    &
                                                      solarflux_p(i,j)))
                    Sw_output%fsw(i,j,k,nz) = Sw_output%fsw(i,j,k,nz) +&
                                                      fswband(i,j,k,nz)
                  end if
                end do
              end do
            end do
 
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                if  ( daylight(i,j) ) then
                  Sw_output%dfsw_dir_sfc(i,j,nz) =   &
                            Sw_output%dfsw_dir_sfc(i,j,nz) +   &
                              sumtr_dir(i,j,KERAD+1,nz)*solarflux_p(i,j)
                  Sw_output%ufsw_dir_sfc(i,j,nz) =   &
                            Sw_output%ufsw_dir_sfc(i,j,nz) +   &
                              sumtr_dir_up(i,j,nz)*solarflux_p(i,j)
                  Sw_output%ufsw_dif_sfc(i,j,nz) =   &
                             Sw_output%ufsw_dif_sfc(i,j,nz) +   &
                               sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                end if
              end do
            end do

            if (nband > NIRBANDS) then
              do j = JSRAD,JERAD
                do i = ISRAD,IERAD
                  if ( daylight(i,j) ) then
                    Sw_output%dfsw_vis_sfc(i,j,nz) =   &
                          Sw_output%dfsw_vis_sfc(i,j,nz) +   &
                                sumtr(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc(i,j,nz) =   &
                           Sw_output%ufsw_vis_sfc(i,j,nz) +   &
                                 sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%dfsw_vis_sfc_dir(i,j,nz) =   &
                            Sw_output%dfsw_vis_sfc_dir(i,j,nz) +   &
                              sumtr_dir(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc_dir(i,j,nz) =   &
                            Sw_output%ufsw_vis_sfc_dir(i,j,nz) +   &
                              sumtr_dir_up(i,j,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc_dif(i,j,nz) =   &
                             Sw_output%ufsw_vis_sfc_dif(i,j,nz) +   &
                                 sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                  end if
                end do
              end do
            endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
            do k = KSRAD,KERAD
              do j = JSRAD,JERAD
                do i = ISRAD,IERAD
                  if ( daylight(i,j) ) then
                    hswband = (fswband(i,j,k+1,nz) -    &
                            fswband(i,j,k,nz) )*gocpdp(i,j,k)
                    Sw_output%hsw(i,j,k,nz) =    &
                          Sw_output%hsw(i,j,k,nz) + hswband
                  end if
                end do
              end do
            end do

!----------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!---------------------------------------------------------------------
            if (nprofile == 1) then  ! clr sky need be done only for 
                                     ! first cloud profile
              if (Rad_control%do_totcld_forcing) then
                where ( daylight )
                  Sw_output%dfsw_dir_sfc_clr(:,:,nz) =   &
                     Sw_output%dfsw_dir_sfc_clr(:,:,nz) +   &
                      sumtr_dir_clr(:,:,KERAD+1,nz)*solarflux_p(:,:)
                end where
                if (nband > NIRBANDS) then
                  where ( daylight )
                    Sw_output%dfsw_vis_sfc_clr(:,:,nz) =   &
                       Sw_output%dfsw_vis_sfc_clr(:,:,nz) +   &
                          sumtrclr(:,:,KERAD+1,nz)*solarflux_p(:,:)
                  end where
                endif
                if (nband == Solar_spect%visible_band_indx) then
                  Sw_output%bdy_flx_clr(:,:,1,nz) =      &
                           sumreclr(:,:,1,nz)*solarflux_p(:,:)
                  Sw_output%bdy_flx_clr(:,:,3,nz) =    &
                          sumtrclr(:,:,KERAD+1,nz)*solarflux_p(:,:) - &
                          sumreclr(:,:,KERAD+1,nz)*solarflux_p(:,:) 
                endif
                if (nband == onepsix_band_indx) then
                  Sw_output%bdy_flx_clr(:,:,2,nz) =    &
                          sumreclr(:,:,1,nz)*solarflux_p(:,:)
                  Sw_output%bdy_flx_clr(:,:,4,nz) =    &
                         sumtrclr(:,:,KERAD+1,nz)*solarflux_p(:,:)  -  &
                         sumreclr(:,:,KERAD+1,nz)*solarflux_p(:,:) 
                endif
          
                if (do_esfsw_band_diagnostics) then
                  do k = KSRAD,KERAD+1
                    dfswbandclr(:,:,k,nz) =     &
                               sumtrclr(:,:,k,nz)*solarflux_p(:,:)
                    ufswbandclr(:,:,k,nz) =    &
                               sumreclr(:,:,k,nz)*solarflux_p(:,:)
                  end do
                endif

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total     
!    spectral values.                                                 
!----------------------------------------------------------------------c
                do k = KSRAD,KERAD+1
                  Sw_output%dfswcf(:,:,k,nz) =    &
                          Sw_output%dfswcf(:,:,k,nz) +  &
                                sumtrclr(:,:,k,nz)*solarflux_p(:,:)
                  Sw_output%ufswcf(:,:,k,nz) =      &
                          Sw_output%ufswcf(:,:,k,nz) +  &
                                sumreclr(:,:,k,nz)*solarflux_p(:,:)
                  fswbandclr(:,:,k,nz) =    &
                        (sumreclr(:,:,k,nz)*solarflux_p(:,:)) - &
                        (sumtrclr(:,:,k,nz)*solarflux_p(:,:))
                  Sw_output%fswcf(:,:,k,nz) =    &
                                Sw_output%fswcf(:,:,k,nz) +    &
                                              fswbandclr(:,:,k,nz)
                end do

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total    
!    spectral values.                                               
!----------------------------------------------------------------------c
                do k = KSRAD,KERAD
                  hswbandclr(:,:,k,nz) =    &
                            (fswbandclr(:,:,k+1,nz) -      &
                                fswbandclr(:,:,k,nz))*gocpdp(:,:,k)
                  Sw_output%hswcf(:,:,k,nz) =   &
                               Sw_output%hswcf(:,:,k,nz) +    &
                                               hswbandclr(:,:,k,nz)
                end do
              endif
            endif ! (nprofile == 1)
          end do  ! (nz loop)
        end do    ! (end of band loop)
      end do      ! (end of nprofile loop)

!----------------------------------------------------------------------
!    if the ica calculation was being done, the fluxes and heating rates
!    which have been summed over nprofiles cloud profiles must be 
!    averaged.
!------------------------------------------------------------------
      if (do_ica_calcs) then
        Sw_output%dfsw_dir_sfc (:,:,:) = Sw_output%dfsw_dir_sfc(:,:,:)*profiles_inverse
        Sw_output%ufsw_dif_sfc (:,:,:) = Sw_output%ufsw_dif_sfc(:,:,:)*profiles_inverse
        Sw_output%dfsw_vis_sfc (:,:,:) = Sw_output%dfsw_vis_sfc(:,:,:)*profiles_inverse
        Sw_output%ufsw_vis_sfc (:,:,:) = Sw_output%ufsw_vis_sfc(:,:,:)*profiles_inverse
        Sw_output%dfsw_vis_sfc_dir (:,:,:) = Sw_output%dfsw_vis_sfc_dir(:,:,:)*profiles_inverse
        Sw_output%ufsw_vis_sfc_dif (:,:,:) = Sw_output%ufsw_vis_sfc_dif(:,:,:)*profiles_inverse
        Sw_output%bdy_flx(:,:,:,:) = Sw_output%bdy_flx(:,:,:,:)*profiles_inverse
        Sw_output%dfsw (:,:,:,:) = Sw_output%dfsw(:,:,:,:)*profiles_inverse
        Sw_output%ufsw (:,:,:,:) = Sw_output%ufsw(:,:,:,:)*profiles_inverse
        Sw_output%fsw(:,:,:,:) = Sw_output%fsw(:,:,:,:)*profiles_inverse
        Sw_output%hsw(:,:,:,:) = Sw_output%hsw(:,:,:,:)*profiles_inverse
      endif

      do nz=1,nzens
        where ( daylight )
          Sw_output%dfsw_dif_sfc(:,:,nz ) =   &
                          Sw_output%dfsw(:,:,KERAD+1,nz ) -   &
                                  Sw_output%dfsw_dir_sfc(:,:,nz )
        end where
      end do

      if (Rad_control%do_totcld_forcing) then
        do nz=1,nzens
          where ( daylight )
            Sw_output%dfsw_dif_sfc_clr(:,:,nz ) =   &
                              Sw_output%dfswcf(:,:,KERAD+1,nz) -   &
                              Sw_output%dfsw_dir_sfc_clr(:,:,nz)
          end where
        end do
      endif

      do nz=1,nzens
        where ( daylight )
          Sw_output%dfsw_vis_sfc_dif(:,:,nz) =   &
                            Sw_output%dfsw_vis_sfc(:,:,nz) -   &
                            Sw_output%dfsw_vis_sfc_dir(:,:,nz) 
        end where
      end do

!--------------------------------------------------------------------
!    convert sw fluxes to cgs and then back to  mks units.
!---------------------------------------------------------------------
      Sw_output%fsw(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%fsw(:,:,:,:))
      Sw_output%dfsw(:,:,:,:) =    &
                            1.0E-03*(1.0E+03*Sw_output%dfsw(:,:,:,:))
      Sw_output%ufsw(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufsw(:,:,:,:))
      if (Rad_control%do_totcld_forcing) then
        Sw_output%fswcf(:,:,:,:) =   &
                            1.0E-03*(1.0E+03*Sw_output%fswcf(:,:,:,:))
        Sw_output%dfswcf(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%dfswcf(:,:,:,:))
        Sw_output%ufswcf(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufswcf(:,:,:,:))
      endif

!---------------------------------------------------------------------


end subroutine swresf



!####################################################################

subroutine esfsw_driver_end

!---------------------------------------------------------------------
!    esfsw_driver_end is the destructor for esfsw_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        write(6,*) "ERROR: esfsw_driver_mod has not been initialized"
        stop
      endif

!--------------------------------------------------------------------
!    close out the modules that this module initialized.
!--------------------------------------------------------------------
      call esfsw_parameters_end

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine esfsw_driver_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#################################################################
! <SUBROUTINE NAME="compute_aerosol_optical_props">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comput(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine compute_aerosol_optical_props    &
          (Atmos_input, Aerosol, Aerosol_props, including_volcanoes,  &
           Aerosol_diags, r, including_aerosols, naerosol_optical, &
           daylight, aeroextopdep, aerosctopdep, aeroasymfac)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

type(atmos_input_type),        intent(in)    :: Atmos_input
type(aerosol_type),            intent(in)    :: Aerosol     
type(aerosol_properties_type), intent(in)    :: Aerosol_props
logical,                       intent(in)    :: including_volcanoes
type(aerosol_diagnostics_type),intent(inout) :: Aerosol_diags
real,dimension(:,:,:,:),       intent(inout) :: r
logical,                       intent(in)    :: including_aerosols  
integer,                       intent(in)    :: naerosol_optical
logical,dimension(:,:),        intent(in)    :: daylight
real, dimension(:,:,:,:),      intent(out)   :: aeroextopdep, &
                                                aerosctopdep, &
                                                aeroasymfac 


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Aerosol_props  aerosol radiative property input data for the 
!                     radiation package
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 
      real, dimension (size(Atmos_input%temp,1))  :: &
                      arprod, asymm,   arprod2
      real, dimension (size(atmos_input%temp,1),size(atmos_input%temp,2),size(atmos_input%temp,3)-1)  :: &
                      sum_g_omega_tau, sum_ext,      sum_sct

      integer, dimension (size(Atmos_input%temp,1),size(Atmos_input%temp,2),size (Atmos_input%temp, 3)-1 ) ::   &
                      opt_index_v3

      real, dimension (naerosol_optical)   ::            &
                      aerext,          aerssalb,       aerasymm
      
      real        :: aerext_i, aerssalb_i, aerasymm_i
      integer     :: j, i, k, nband, nsc, irh
      integer     :: israd, jsrad, ierad, jerad, ksrad, kerad
      integer     :: opt_index_v4, opt_index_v5, opt_index_v6, opt_index_v7, opt_index_v8, opt_index_v9, opt_index_v10
      integer     :: naerosoltypes_used


!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------

      if (including_volcanoes) then
        write(6,*) "ERROR: including_volcanoes is not supported"
        stop
      end if

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size(Atmos_input%temp,1)
      jerad = size(Atmos_input%temp,2)
      kerad = size(Atmos_input%temp,3) - 1

      naerosoltypes_used = size(Aerosol%aerosol,4)

      if (including_aerosols) then

        if (Rad_control%using_im_bcsul) then
          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD

!-------------------------------------------------------------------
!    define the local variables for the band values of aerosol and 
!    cloud single scattering parameters.                               
!    note: the unit for the aerosol extinction is kilometer**(-1).     
!--------------------------------------------------------------------
                irh = MIN(100, MAX( 0,     &
                    NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                opt_index_v3(i,j,k) = &
                            Aerosol_props%sulfate_index (irh, &
                                          Aerosol_props%ivol(i,j,k))
              end do
            end do
          end do
        else
          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                irh = MIN(100, MAX( 0,     &
                    NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                opt_index_v3(i,j,k) = &
                            Aerosol_props%sulfate_index (irh, 0)
              end do
            end do
          end do
        end if ! (Rad_control%using_im_bcsul) ..else..
              
        do nband = 1,Solar_spect%nbands
          aerext(:) = Aerosol_props%aerextband(nband,:)
          aerssalb(:) = Aerosol_props%aerssalbband(nband,:)
          aerasymm(:) = Aerosol_props%aerasymmband(nband,:)

!---------------------------------------------------------------------
!    calculate scattering properties for all aerosol constituents 
!    combined.
!---------------------------------------------------------------------
          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                sum_g_omega_tau(i,j,k) = 0.0
                sum_ext(i,j,k) = 0.
                sum_sct(i,j,k) = 0.
              end do
            end do
          end do
          do nsc = 1,NAEROSOLTYPES_USED
            do k = KSRAD,KERAD
              do j = JSRAD,JERAD
                if (Aerosol_props%optical_index(nsc) > 0) then
                  aerext_i =     &
                          aerext(Aerosol_props%optical_index(nsc))
                  aerssalb_i =     &
                          aerssalb(Aerosol_props%optical_index(nsc))
                  aerasymm_i =     &
                          aerasymm(Aerosol_props%optical_index(nsc))
                  do i = ISRAD,IERAD
                    arprod(i) =    &
                          aerext_i*(1.e3*Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +     &
                                    aerasymm_i*(aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                                    Aerosol_props%sulfate_flag) then
                  do i = ISRAD,IERAD
                    aerext_i = aerext(opt_index_v3(i,j,k))
                    aerssalb_i = aerssalb(opt_index_v3(i,j,k))
                    aerasymm_i = aerasymm(opt_index_v3(i,j,k))
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) +    &
                                 aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                                 aerasymm_i*                &
                                 (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                                   Aerosol_props%bc_flag) then
                  do i = ISRAD,IERAD
                    aerext_i = aerext(opt_index_v3(i,j,k))
                    aerssalb_i = aerssalb(opt_index_v3(i,j,k))
                    aerasymm_i = aerasymm(opt_index_v3(i,j,k))
                    arprod(i) = aerext_i *    &
                                  (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) +       &
                                  aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) + &
                                 aerasymm_i *              &
                                 (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                             Aerosol_props%omphilic_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v4 =    &
                            Aerosol_props%omphilic_index( irh )
                    aerext_i = aerext(opt_index_v4)
                    aerssalb_i = aerssalb(opt_index_v4)
                    aerasymm_i = aerasymm(opt_index_v4)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +   &
                                aerasymm_i*                  &
                                (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                                Aerosol_props%bcphilic_flag) then
                  if (Rad_control%using_im_bcsul) then
                    do i = ISRAD,IERAD
                      aerext_i = aerext(opt_index_v3(i,j,k))
                      aerssalb_i = aerssalb(opt_index_v3(i,j,k))
                      aerasymm_i = aerasymm(opt_index_v3(i,j,k))
                      arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(i) = aerssalb_i*arprod(i)
                      asymm(i)   = aerasymm_i
                      sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                      sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                      sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                                aerasymm_i *                  &
                                (aerssalb_i*arprod(i))
                    end do
                  else  ! (using_im_bcsul)
                    do i = ISRAD,IERAD
                      irh = MIN(100, MAX( 0,     &
                          NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                      opt_index_v5 =    &
                                  Aerosol_props%bcphilic_index( irh)
                      aerext_i = aerext(opt_index_v5)
                      aerssalb_i = aerssalb(opt_index_v5)
                      aerasymm_i = aerasymm(opt_index_v5)
                      arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(i) = aerssalb_i*arprod(i)
                      asymm(i)   = aerasymm_i
                      sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                      sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                      sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                               aerasymm_i *                   &
                               (aerssalb_i*arprod(i))
                    end do
                  endif  !(using_im_bcsul)
                else if (Aerosol_props%optical_index(nsc) == &
                              Aerosol_props%seasalt1_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v6 =    &
                                Aerosol_props%seasalt1_index( irh )
                    aerext_i = aerext(opt_index_v6)
                    aerssalb_i = aerssalb(opt_index_v6)
                    aerasymm_i = aerasymm(opt_index_v6)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                     aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                                aerasymm_i *                &
                                (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                        Aerosol_props%seasalt2_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v7 =    &
                            Aerosol_props%seasalt2_index( irh )
                    aerext_i = aerext(opt_index_v7)
                    aerssalb_i = aerssalb(opt_index_v7)
                    aerasymm_i = aerasymm(opt_index_v7)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) +  arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) + &
                               aerasymm_i *                &
                               (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                               Aerosol_props%seasalt3_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v8 =    &
                            Aerosol_props%seasalt3_index( irh )
                    aerext_i = aerext(opt_index_v8)
                    aerssalb_i = aerssalb(opt_index_v8)
                    aerasymm_i = aerasymm(opt_index_v8)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                                aerasymm_i *                &
                                (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                             Aerosol_props%seasalt4_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v9 =    &
                            Aerosol_props%seasalt4_index( irh )
                    aerext_i = aerext(opt_index_v9)
                    aerssalb_i = aerssalb(opt_index_v9)
                    aerasymm_i = aerasymm(opt_index_v9)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                   aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) +  &
                                 aerasymm_i*                &
                                 (aerssalb_i*arprod(i))
                  end do
                else if (Aerosol_props%optical_index(nsc) == &
                           Aerosol_props%seasalt5_flag) then
                  do i = ISRAD,IERAD
                    irh = MIN(100, MAX( 0,     &
                        NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                    opt_index_v10 =    &
                            Aerosol_props%seasalt5_index( irh )
                    aerext_i = aerext(opt_index_v10)
                    aerssalb_i = aerssalb(opt_index_v10)
                    aerasymm_i = aerasymm(opt_index_v10)
                    arprod(i) = aerext_i *    &
                                 (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                    arprod2(i) = aerssalb_i*arprod(i)
                    asymm(i)   = aerasymm_i
                    sum_ext(i,j,k) = sum_ext(i,j,k) + arprod(i)
                    sum_sct(i,j,k) = sum_sct(i,j,k) + &
                                  aerssalb_i*arprod(i)
                    sum_g_omega_tau(i,j,k) = sum_g_omega_tau(i,j,k) + &
                               aerasymm_i *                &
                               (aerssalb_i*arprod(i))
                  end do
                endif

                if (Sw_control%do_cmip_diagnostics) then
                  if (nband == Solar_spect%visible_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,1) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,1) =    &
                                            arprod(1) - arprod2(1)
                      Aerosol_diags%asymdep(i,j,k,nsc,1) = asymm(1)
                    end do
                  endif
                  if (nband == Solar_spect%eight70_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,6) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,6) =    &
                                            arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,6) = asymm(i)
                    end do
                  endif
                  if (nband == Solar_spect%one_micron_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,2) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,2) =    &
                                               arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,2) = asymm(i)
                    end do
                  endif
                  if (nband == Solar_spect%w340_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,7) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,7) =    &
                                                arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,7) = asymm(i)
                    end do
                  endif
                  if (nband == Solar_spect%w380_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,8) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,8) =    &
                                                arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,8) = asymm(i)
                    end do
                  endif
                  if (nband == Solar_spect%w440_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,9) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,9) =    &
                                               arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,9) = asymm(i)
                    end do
                  endif
                  if (nband == Solar_spect%w670_band_indx) then
                    do i = ISRAD,IERAD
                      Aerosol_diags%extopdep(i,j,k,nsc,10) = arprod(i)
                      Aerosol_diags%absopdep(i,j,k,nsc,10) =    &
                                                arprod(i) - arprod2(i)
                      Aerosol_diags%asymdep(i,j,k,nsc,10) = asymm(i)
                    end do
                  endif
                endif

              end do
            end do
          end do ! (nsc)

!----------------------------------------------------------------------
!    add the effects of volcanic aerosols, if they are to be included.
!    include generation of diagnostics in the visible (0.55 micron) and
!    nir band (1.0 micron).
!----------------------------------------------------------------------
!            if (including_volcanoes) then
!              do k = KSRAD,KERAD
!                do j = JSRAD,JERAD
!                  do i = ISRAD,IERAD
!                    deltaz = Atmos_input%deltaz(i,j,k)              
!                    sum_ext(i,j,k) = sum_ext(i,j,k) +    &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                    sum_sct(i,j,k) = sum_sct(i,j,k) +    &
!                                 Aerosol_props%sw_ssa(i,j,k,nband)*  &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                    sum_g_omega_tau(i,j,k) =   &
!                                 sum_g_omega_tau(i,j,k) +&
!                                 Aerosol_props%sw_asy(i,j,k,nband)* &
!                                 Aerosol_props%sw_ssa(i,j,k,nband)*  &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                    if (Sw_control%do_cmip_diagnostics) then
!                      if (nband == Solar_spect%visible_band_indx) then
!                           Aerosol_diags%extopdep_vlcno(i,j,k,1) =   &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                           Aerosol_diags%absopdep_vlcno(i,j,k,1) =   &
!                            (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
!                                Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                deltaz
!                      endif
!                      if (nband == Solar_spect%eight70_band_indx) then
!                           Aerosol_diags%extopdep_vlcno(i,j,k,3) =   &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                           Aerosol_diags%absopdep_vlcno(i,j,k,3) =   &
!                            (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
!                                Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                deltaz
!                      endif
!                      if (nband == Solar_spect%one_micron_indx) then
!                           Aerosol_diags%extopdep_vlcno(i,j,k,2) =   &
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                  deltaz
!                            Aerosol_diags%absopdep_vlcno(i,j,k,2) =   &
!                             (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
!                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
!                                 deltaz
!                       endif
!                    endif
!                  end do
!                end do
!              end do
!            endif   ! (including_volcanoes)
!
!----------------------------------------------------------------------
          do k = KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                aeroextopdep(i,j,k,nband) = sum_ext(i,j,k) 
                aerosctopdep(i,j,k,nband) = sum_sct(i,j,k) 
                aeroasymfac(i,j,k,nband) = sum_g_omega_tau(i,j,k) / &
                                              (sum_sct(i,j,k) + 1.0e-30 )
              end do
            end do
          end do

        end do ! (nband)

        !if (including_volcanoes) then
        !  if (do_coupled_stratozone ) then
        !    nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
        !    if (nextinct  /= NO_TRACER) &
        !          r(i,j,:,nextinct) = Aerosol_props%sw_ext(i,j,:,4)
        !  endif  
        !endif  
          
      else  ! (if not including_aerosols)
          
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (daylight(i,j) .or. Sw_control%do_cmip_diagnostics) then
              do nband = 1,Solar_spect%nbands                
                do k = KSRAD,KERAD
                  aeroextopdep(i,j,k,nband) = 0.0                    
                  aerosctopdep(i,j,k,nband) = 0.0                  
                  aeroasymfac(i,j,k,nband) = 0.0                 
                end do
              end do ! (nband)
            endif  ! (daylight or cmip_diagnostics)              
            !if (including_volcanoes) then
            !  if (do_coupled_stratozone ) then
            !    nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
            !    if (nextinct  /= NO_TRACER) &
            !          r(i,j,:,nextinct) = Aerosol_props%sw_ext(i,j,:,4)
            !  endif  
            !endif  
          end do ! (i loop)
        end do   ! (j loop)

      endif ! (including_aerosols)
    
!---------------------------------------------------------------------
!

!---------------------------------------------------------------------


end subroutine compute_aerosol_optical_props

!#################################################################
! <SUBROUTINE NAME="compute_gas_props">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comput(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine compute_gas_props (Atmos_input, Rad_gases, Astro,   &
                              daylight, gasopdep)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

type(atmos_input_type),        intent(in)    :: Atmos_input
type(radiative_gases_type),    intent(in)    :: Rad_gases   
type(astronomy_type),          intent(in)    :: Astro
logical, dimension(:,:),       intent(in)    :: daylight
real, dimension(:,:,:,:,:),    intent(out)   :: gasopdep              


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Astro          astronomy_type structure
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 

      real, dimension (size(Atmos_input%temp,1),size(Atmos_input%temp,2),size(Atmos_input%temp,3)-1)  :: &
                   deltaz,     qo3,         rh2o,                    &
                   efftauo2,   efftauco2,   efftauch4,   efftaun2o, &
                   wh2ostr,    wo3,         wo2,         quenchfac, &
                   delpdig,     deltap,      tco2,                  &
                   tch4,       tn2o,        to2,         wh2o
      
      real, dimension(size(atmos_input%temp,1),size(Atmos_input%temp,2)) :: &
                   opdep

           
      real, dimension (size(Atmos_input%temp,1),size(Atmos_input%temp,2),size(Atmos_input%temp,3))  :: &
            alphaco2,        alphaco2str,    alphao2,          &
            alphao2str,      alphach4,       alphach4str,      &
            alphan2o,        alphan2ostr,    scale,            &
            scalestr,        totco2,         totco2str,        &
            toto2,           toto2str,       totch4,           &
            totch4str,       totn2o,         totn2ostr,        &
            press,           pflux,          pflux_mks,        &
            temp,            z

      real, dimension (size(Atmos_input%temp,1),size(Atmos_input%temp,2)) :: cosangsolar
      real, dimension (size(Atmos_input%temp,1),size(Atmos_input%temp,2)) :: denom
      real :: wtquench
      real :: rrvco2 
      real :: rrvch4, rrvn2o

      integer  :: j, i, k, ng, nband, kq
      integer  :: np, nf
      integer  :: israd, jsrad, ierad, jerad, ksrad, kerad


!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size(Atmos_input%temp,1)
      jerad = size(Atmos_input%temp,2)
      kerad = size(Atmos_input%temp,3) - 1

!---------------------------------------------------------------------
!    initialize local variables.                                        
!------------------------------------------------------------------
      alphaco2   (:,:,1) = 0.0
      alphaco2str(:,:,1) = 0.0
      alphao2    (:,:,1) = 0.0
      alphao2str (:,:,1) = 0.0
      alphach4   (:,:,1) = 0.0
      alphach4str(:,:,1) = 0.0
      alphan2o   (:,:,1) = 0.0
      alphan2ostr(:,:,1) = 0.0

      rrvco2 = Rad_gases%rrvco2
      rrvch4 = Rad_gases%rrvch4
      rrvn2o = Rad_gases%rrvn2o

!---------------------------------------------------------------------
!  convert to cgs and then back to mks for consistency with previous 
!---------------------------------------------------------------------
      do k = KSRAD,KERAD+1
        press(:,:,k) = 0.1*(10.0*Atmos_input%press(:,:,k))
        pflux(:,:,k) =     (10.0*Atmos_input%pflux(:,:,k))
        temp (:,:,k) = Atmos_input%temp  (:,:,k)
      end do
      do k = KSRAD,KERAD
        rh2o  (:,:,k) = Atmos_input%rh2o(:,:,k)
        qo3   (:,:,k) = Rad_gases%qo3(:,:,k)
        deltaz(:,:,k) = Atmos_input%deltaz(:,:,k)
      end do

!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
      pflux_mks(:,:,:) = pflux(:,:,:)*1.0E-1

      do k = KSRAD+1,KERAD+1
        deltap(:,:,k-1) = pflux_mks(:,:,k) - pflux_mks(:,:,k-1)
        delpdig(:,:,k-1) = deltap(:,:,k-1)/ GRAV
        scalestr(:,:,k) = pflux_mks(:,:,k) 
        scale(:,:,k) = scalestr(:,:,k)*pflux_mks(:,:,k)/pstd_mks
      end do
 

      do k = KSRAD,KERAD
        wh2ostr(:,:,k) = rh2o(:,:,k)*delpdig(:,:,k)
        wo3(:,:,k)     = qo3(:,:,k)*delpdig(:,:,k)
        wo2(:,:,k) = o2mixrat*(WTMO2/WTMAIR)*delpdig(:,:,k)
      end do
  
!---------------------------------------------------------------------
!    if quenching factor effects are desired, calculate the height above
!    the surface of the model flux levels.
!---------------------------------------------------------------------
      if (do_quench) then
        z(:,:,KERAD+1) = 0.0
        do k = KERAD,KSRAD,-1
          z(:,:,k) = z(:,:,k+1) + deltaz(:,:,k)
        end do
          
!---------------------------------------------------------------------
!    define the quenching factor for each grid point.
!---------------------------------------------------------------------
        do k = KSRAD,KERAD
          do j = jsrad,jerad
            do i = israd,ierad
              if (z(i,j,k) < co2_quenchfac_height(1) ) then
                quenchfac(i,j,k) = 1.0
              else if (z(i,j,k) > co2_quenchfac_height(30) ) then 
                quenchfac(i,j,k) = 0.037
              else
                do kq = 1,29
                  if (z(i,j,k) > co2_quenchfac_height(kq) .and. &
                        z(i,j,k) <= co2_quenchfac_height(kq+1)) then
                    wtquench = (z(i,j,k) - co2_quenchfac_height(kq))/ &
                                 (co2_quenchfac_height(kq+1) - &
                                  co2_quenchfac_height(kq))
                    quenchfac(i,j,k) = (1. - wtquench)*   &
                                           co2_quenchfac(kq) +   &
                                      wtquench*co2_quenchfac(kq+1)
                    exit
                  endif
                end do
              endif
            end do
          end do
        end do
      else
        quenchfac(:,:,:) = 1.0
      endif !(do_quench)


      do ng = 1,NSOLWG
        cosangsolar  = Astro%cosz(:,:)
        where ( cosangsolar==0.0 )
          cosangsolar = 1.0  
        end where

!----------------------------------------------------------------------c
!    define the scaled and unscaled co2 and o2 pathlengths in 
!    centimeter-atm, and the unscaled h2o and o3 amounts in   
!    kgrams/meter**2. 
!    cm-atm needed as units because of c2co2 having those units.
!----------------------------------------------------------------------c
        denom = 1.0/(GRAV*rhoair*cosangsolar*2.0)
        do k = KSRAD+1,KERAD+1
          totco2(:,:,k) = 1.0E+02*rrvco2*scale(:,:,k)*denom       
          totco2str(:,:,k) = 2.0E+02*rrvco2*scalestr(:,:,k)*denom     
          toto2(:,:,k) = 1.0E+02*o2mixrat*scale(:,:,k)*denom      
          toto2str(:,:,k) = 2.0E+02*o2mixrat*scalestr(:,:,k)*denom    
        end do
        if (do_ch4_sw_effects) then
          do k = KSRAD+1,KERAD+1    
            totch4(:,:,k) = 1.0E+02*rrvch4*scale(:,:,k)*denom     
            totch4str(:,:,k) = 2.0E+02*rrvch4*scalestr(:,:,k)*denom     
          end do
        endif
        if (do_n2o_sw_effects) then
          do k = KSRAD+1,KERAD+1      
            totn2o(:,:,k) = 1.0E+02*rrvn2o*scale(:,:,k)*denom     
            totn2ostr(:,:,k) = 2.0E+02*rrvn2o*scalestr(:,:,k)*denom      
          end do
        endif

        np = 0
        do nband = 1, NBANDS

!-------------------------------------------------------------------
!    define the h2o scaled gas amounts in kgrams/meter**2            
!---------------------------------------------------------------------
          if (nband <= nh2obands) then
            do k = KSRAD,KERAD
              wh2o(:,:,k) = rh2o(:,:,k)*delpdig(:,:,k)*   &
                  exp(powph2o(nband)*alog(press(:,:,k)*p0h2o(nband)))
            end do
 
!---------------------------------------------------------------------
!    calculate the "effective" co2, o2, ch4 and n2o gas optical depths 
!    for the appropriate absorbing bands.                               
!    note: for large zenith angles, alpha can exceed 1. In this case,a
!    the optical depths are set to the previous layer values.          
!-------------------------------------------------------------------
            if ( c1co2(nband).ne.1.0E-99 ) then
              do k = KSRAD+1,KERAD+1
                do j = jsrad,jerad
                  do i = israd,ierad
                    if (totco2(i,j,k) < totco2max(nband) .and.  &
                       totco2str(i,j,k) < totco2strmax(nband))  then
                      alphaco2(i,j,k) =     &
                           c1co2(nband)*exp(c3co2(nband)* &
                               alog((totco2(i,j,k) + c2co2(nband))))  -  &
                                                        c4co2(nband)
                      alphaco2str(i,j,k) = &
                        c1co2str(nband)*exp(c3co2str(nband)*  &
                          alog((totco2str(i,j,k) + c2co2str(nband)))) - &
                                                  c4co2str(nband)
                      tco2(i,j,k-1) =      &
                         (1.0 - alphaco2(i,j,k))*   &
                                        (1.0 - alphaco2str(i,j,k))/ &
                         ((1.0 - alphaco2(i,j,k-1))*    &
                                        (1.0 - alphaco2str(i,j,k-1)))
                      efftauco2(i,j,k-1) = -cosangsolar(i,j)*alog( tco2(i,j,k-1))
                    else if (k > KSRAD+1) then
                      efftauco2(i,j,k-1) = efftauco2(i,j,k-2)
                    else
                      efftauco2(i,j,k-1) = 0.0
                    end if
                  end do
                end do
              end do
            else    !( c1co2(nband).ne.1.0E-99 ) 
              efftauco2(:,:,:) = 0.0
            end if  !( c1co2(nband).ne.1.0E-99 ) 

            if (do_ch4_sw_effects) then
              if (c1ch4(nband).ne.1.0E-99 ) then
                do k = KSRAD+1,KERAD+1
                  do j = jsrad,jerad
                    do i = israd,ierad
                      if (totch4(i,j,k) < totch4max(nband) .and.  &
                          totch4str(i,j,k) < totch4strmax(nband))  then
                        alphach4(i,j,k) =    &
                            c1ch4(nband)*exp(c3ch4(nband)*&
                            alog((totch4(i,j,k) + c2ch4(nband))))  -   &
                                                         c4ch4(nband)
                        alphach4str(i,j,k) = &
                          c1ch4str(nband)*exp(c3ch4str(nband)*  &
                           alog((totch4str(i,j,k) + c2ch4str(nband)))) - &
                                                      c4ch4str(nband)
                        tch4(i,j,k-1) = &
                                (1.0 - alphach4(i,j,k))*    &
                                         (1.0 - alphach4str(i,j,k))/ &
                                 ((1.0 - alphach4(i,j,k-1))*   &
                                         (1.0 - alphach4str(i,j,k-1)))
                        efftauch4(i,j,k-1) = -cosangsolar(i,j)*alog(tch4(i,j,k-1))
                      else if (k > KSRAD+1) then
                        efftauch4(i,j,k-1) = efftauch4(i,j,k-2)
                      else
                        efftauch4(i,j,k-1) = 0.0
                      end if
                    end do
                  end do
                end do
              else    !( c1ch4(nband).ne.1.0E-99 )
                efftauch4(:,:,:) = 0.0
              end if  !( c1ch4(nband).ne.1.0E-99 )
            else    !do_ch4 = .false.
              efftauch4(:,:,:) = 0.0
            end if

            if (do_n2o_sw_effects) then
              if ( c1n2o(nband).ne.1.0E-99 ) then
                do k = KSRAD+1,KERAD+1
                  do j = jsrad,jerad
                    do i = israd,ierad
                      if (totn2o(i,j,k) < totn2omax(nband) .and.  &
                             totn2ostr(i,j,k) < totn2ostrmax(nband)) then
                        alphan2o(i,j,k) = &
                             c1n2o(nband)*exp(c3n2o(nband)* &
                                alog((totn2o(i,j,k) +c2n2o(nband)))) -  &
                                                       c4n2o(nband)
                        alphan2ostr(i,j,k) = &
                         c1n2ostr(nband)*exp(c3n2ostr(nband)*  &
                         alog((totn2ostr(i,j,k) + c2n2ostr(nband)))) -  &
                                                    c4n2ostr(nband)
                        tn2o(i,j,k-1) = &
                                 (1.0 - alphan2o(i,j,k)) *  &
                                           (1.0 - alphan2ostr(i,j,k))/ &
                                 (( 1.0 - alphan2o(i,j,k-1)) *  &
                                         (1.0 - alphan2ostr(i,j,k-1)))
                        efftaun2o(i,j,k-1) = -cosangsolar(i,j)*alog(tn2o(i,j,k-1))
                      else if (k > KSRAD+1) then
                        efftaun2o(i,j,k-1) = efftaun2o(i,j,k-2)
                      else
                        efftaun2o(i,j,k-1) = 0.0
                      end if
                    end do
                  end do
                end do
              else    !( c1n2o(nband).ne.1.0E-99 )
                efftaun2o(:,:,:) = 0.0
              end if  !( c1n2o(nband).ne.1.0E-99 )
            else  !do_n2o = .false.
              efftaun2o(:,:,:) = 0.0
            end if

            if ( c1o2(nband).ne.1.0E-99 ) then
              do k = KSRAD+1,KERAD+1
                do j = jsrad,jerad
                  do i = israd,ierad
                    if (toto2(i,j,k) .lt. toto2max(nband) .and.   &
                        toto2str(i,j,k) .lt. toto2strmax(nband)) then
                      alphao2(i,j,k) = c1o2(nband)*exp( c3o2(nband)* &
                                   alog((toto2(i,j,k) + c2o2(nband)))) - &
                                                          c4o2(nband)
                      alphao2str(i,j,k) = &
                          c1o2str(nband)*exp(c3o2str(nband)*  &
                              alog((toto2str(i,j,k) + c2o2str(nband)))) &
                                                      - c4o2str(nband)
                      to2(i,j,k-1) = &
                             (1.0 - alphao2(i,j,k))*  &
                                  (1.0 - alphao2str(i,j,k) )/ &
                              ((1.0 - alphao2(i,j,k-1)) *  &
                                          (1.0 - alphao2str(i,j,k-1)))
                      efftauo2(i,j,k-1) = -cosangsolar(i,j)*alog(to2(i,j,k-1))
                    else if ( k>KSRAD+1 ) then
                      efftauo2(i,j,k-1) = efftauo2(i,j,k-2)
                    else
                      efftauo2(i,j,k-1) = 0.0
                    end if
                  end do
                end do
              end do
            else   !  ( c1o2(nband).ne.1.0E-99 ) 
              efftauo2(:,:,:) = 0.0
            end if  !  ( c1o2(nband).ne.1.0E-99 ) 
              
          end if  ! (nband <= nh2obands)
 
!---------------------------------------------------------------------
!    calculate the "effective" o2 gas optical depths for the Schuman- 
!    Runge band.                                                        
!-------------------------------------------------------------------
          if ( nband.EQ.NBANDS ) then
            if (do_herzberg) then  
              do k = KSRAD+1,KERAD+1
                do j = jsrad,jerad
                  do i = israd,ierad
                    if ( toto2str(i,j,k)<toto2strmaxschrun) then
                      alphao2str(i,j,k) =  &
                           c1o2strschrun*exp( c3o2strschrun*&
                              alog((toto2str(i,j,k) + c2o2strschrun))) - &
                                                         c4o2strschrun
                      to2(i,j,k-1) = &
                            (1.0 - alphao2str(i,j,k))/(1.0 - alphao2str(i,j,k-1)) 
                      efftauo2(i,j,k-1) =  -cosangsolar(i,j)*alog(to2(i,j,k-1) )
                      efftauo2(i,j,k-1) = efftauo2(i,j,k-1) +     &
                                                 wo2(i,j,k-1)*herzberg_fac
                    else if (k>KSRAD+1) then
                      efftauo2(i,j,k-1) = efftauo2(i,j,k-2)
                    else
                      efftauo2(i,j,k-1) = 0.0
                    end if
                  end do
                end do
              end do
            else
              do k = KSRAD+1,KERAD+1
                do j = jsrad,jerad
                  do i = israd,ierad
                    if ( toto2str(i,j,k)<toto2strmaxschrun) then
                      alphao2str(i,j,k) =  &
                           c1o2strschrun*exp( c3o2strschrun*&
                              alog((toto2str(i,j,k) + c2o2strschrun))) - &
                                                         c4o2strschrun
                      to2(i,j,k-1) = &
                            (1.0 - alphao2str(i,j,k))/(1.0 - alphao2str(i,j,k-1)) 
                      efftauo2(i,j,k-1) =  -cosangsolar(i,j)*alog(to2(i,j,k-1) )
                    else if (k>KSRAD+1) then
                      efftauo2(i,j,k-1) = efftauo2(i,j,k-2)
                    else
                      efftauo2(i,j,k-1) = 0.0
                    end if
                  end do
                end do
              end do
            end if
          end if

          do nf =1,nfreqpts(nband)
            np = np + 1

            do k = ksrad,kerad
!---------------------------------------------------------------------
!    define the h2o + o3 gas optical depths.                           
!--------------------------------------------------------------------
              if (strterm(np)) then
                opdep(:,:) = kh2o(np)*wh2ostr(:,:,k) + ko3(np)*wo3(:,:,k)
              else
                opdep(:,:) = kh2o(np)*wh2o(:,:,k) + ko3(np)*wo3(:,:,k)
              end if

              where ( daylight )
                gasopdep(:,:,k,np,ng) =    &
                         opdep(:,:) + quenchfac(:,:,k)*efftauco2(:,:,k) +   &
                            efftauo2(:,:,k) + efftauch4(:,:,k) + efftaun2o(:,:,k)
              end where
              
            end do
          end do  ! (nf loop)
          
        end do   ! (nband loop)
      end do  ! (ng loop)

!---------------------------------------------------------------------



end subroutine compute_gas_props


!#####################################################################
!<SUBROUTINE NAME="adding">
! <OVERVIEW>
!  Subroutine that implements doubling and adding technique to combine
!  multiple atmospheric layers to calculate solar fluxes
! </OVERVIEW>
! <DESCRIPTION>
!  This subroutine implements the standard doubling and adding
!  technique to combine reflectance and transmittance of multiple 
!  atmospheric layers to compute solar flux and heating rate.
! </DESCRIPTION>
! <TEMPLATE>
!  call adding ( ix, jx, kx, &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,  &
!                tlayerde, sfcalb, calc_flag, reflectance,   &
!                transmittance)
! </TEMPLATE>
! <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="jx" TYPE="integer">
!  jx is the current latitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="kx" TYPE="integer">
!  ix is the current vertical index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="rlayerdir" TYPE="real">
!  layer reflectivity to direct incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to direct incident beam
! </IN>
! <IN NAME="rlayerdif" TYPE="real">
!  layer reflectivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerde" TYPE="real">
!  layer diffuse transmissivity to direct incident beam
! </IN>
! <IN NAME="sfcalb" TYPE="real">
!  surface albedo
! </IN>
! <IN NAME="calcflag" TYPE="integer">
!  flag to indicate columns where adding is to be done
! </IN>
! <OUT NAME="reflectance" TYPE="real">
!  diffuse reflectance at a level
! </OUT>
! <OUT NAME="transmittance" TYPE="real">
!  diffuse transmittance at a level
! </OUT>
!</SUBROUTINE>
!

subroutine adding (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif,   &
                   tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,  &
                   calc_flag, reflectance, transmittance, tr_dir)
 
!-------------------------------------------------------------------
!    adding calculates the reflection and transmission at flux levels 
!    from the direct and diffuse values of reflection and transmission
!    in the corresponding layers using the adding method.           
!    references:                                                        
!    bowen, m.m., and v. ramaswamy, effects of changes in radiatively
!        active species upon the lower stratospheric temperatures.,    
!        j. geophys. res., 18909-18921, 1994.                         
!--------------------------------------------------------------------

integer, intent(in)                    :: ix, jx, kx
real, dimension(:,:,:),   intent(in)   :: rlayerdir, rlayerdif, &
                                          tlayerdir, tlayerdif, & 
                                          tlayerde
real, dimension (:,:),    intent(in)   :: sfcalb_dir, sfcalb_dif
logical, dimension (:,:), intent(in)   :: calc_flag
real, dimension(:,:,:),   intent(out)  :: reflectance, transmittance, &
                                          tr_dir

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx        dimensions of current physics window            
!    rlayerdir       layer reflectivity to a direct incident beam      
!    tlayerdir       layer transmissivity to a direct incident beam   
!    rlayerdif       layer reflectivity to a diffuse incident beam  
!    tlayerdif       layer transmissivity to a diffuse incident beam  
!    tlayerde        layer transmissivity (non-scattered) to the direct 
!                    incident beam                                 
!    sfcalb_dir      surface albedo, direct beam 
!    sfcalb_dif      surface albedo, diffuse beam
!    calc_flag       flag to indicate columns where adding is to be 
!                    done. calculations not done in "dark" columns and 
!                    on clr sky pass in columns without any clouds.
!
!  intent(out) variables:
!
!    reflectance     reflectance of the scattered radiation at a level 
!    transmittance   transmittance of the scattered radiation at a level
!    tr_dir
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:
 
      real, dimension (size(calc_flag,1),size(calc_flag,2),lbound(rlayerdir,3):ubound(rlayerdir,3)+1) ::  &
                            raddupdif2, raddupdir2

      real, dimension (size(calc_flag,1),size(calc_flag,2),lbound(rlayerdir,3):ubound(rlayerdir,3)  ) ::  &
                                      radddowndif2,  tadddowndir2

      real, dimension (size(calc_flag,1),size(calc_flag,2)) ::      &
              raddupdif2p, raddupdir2p, tlevel2p, radddowndifm,     &
              tadddowndirm
      real :: dm1tl2, dm32, dm3r2, dm3r1p2, alpp2, dm2tl2, rdm2tl2
      integer     ::  k, i, j

!-------------------------------------------------------------------
!   local variables:
!
!      raddupdif2
!      raddupdir2
!      tlevel2
!      radddowndif2
!      tadddowndir2
!      dm1tl2
!      dm2tl2
!      rdm2tl2
!      dm32
!      dm3r2
!      dm3r1p2
!      alpp2
!      raddupdif2
!      raddupdir2p
!      tlevel2p
!      radddowndifm
!      tadddowndirm
!      i,j,k
!
!--------------------------------------------------------------------

!----------------------------------------------------------------------c
!    initialization for the surface layer.                           
!----------------------------------------------------------------------c
 
!------------------------------------------------------------------ 
!    add the inhomogeneous layers upward from the surface to the top of
!    the atmosphere.                                                  
!    radiation incident from above for diffuse beam, reflection of  
!    direct beam and conversion to diffuse.                           
!--------------------------------------------------------------------
      raddupdif2p = sfcalb_dif(:,:)
      raddupdir2p = sfcalb_dir(:,:)
      do k = kx, 1,-1
        do j = 1,size(calc_flag,2)
          do i = 1,size(calc_flag,1)
            dm2tl2    = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*     &
                        raddupdif2p(i,j) )
            rdm2tl2    = dm2tl2*raddupdif2p(i,j)     
            raddupdif2(i,j,k) = rlayerdif(i,j,k) + tlayerdif(i,j,k)*  &
                                rdm2tl2    
            raddupdir2(i,j,k) = rlayerdir(i,j,k) + tlayerde(i,j,k)*   &
                            raddupdir2p(i,j)* dm2tl2 +                     &     
                            (tlayerdir(i,j,k) - tlayerde(i,j,k))*     &
                            rdm2tl2   
            raddupdir2p(i,j) = raddupdir2(i,j,k)
            raddupdif2p(i,j) = raddupdif2(i,j,k)
          end do
        end do
      end do
 
!---------------------------------------------------------------------
!    define the direct transmittance. add the inhomogeneous layers 
!    downward from the second layer to the surface. radiation incident
!    from below for diffuse beam, transmission of direct beam and 
!    conversion to diffuse.                             
!-------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    initialization for the first scattering layer.                   
!-------------------------------------------------------------------
      tlevel2p         = tlayerde(:,:,1)
      radddowndifm    =  rlayerdif(:,:,1)
      tadddowndirm    =  tlayerdir(:,:,1)
      do k= 2,kx    
        do j = 1,size(calc_flag,2)
          do i = 1,size(calc_flag,1)
            dm1tl2 = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*  &
                     radddowndifm(i,j))
            radddowndif2(i,j,k) = rlayerdif(i,j,k) + radddowndifm(i,j)* &
                                  tlayerdif(i,j,k)*dm1tl2      
            tadddowndir2(i,j,k) = tlevel2p(i,j)*(tlayerdir(i,j,k) + &
                                  rlayerdir(i,j,k)*radddowndifm(i,j)* &
                                  dm1tl2) + (tadddowndirm(i,j) -  &
                                  tlevel2p(i,j))*dm1tl2           
    
!---------------------------------------------------------------------
!    add downward to calculate the resultant reflectances and           
!    transmittances at flux levels.                                    
!------------------------------------------------------------------
            dm32  = 1.0/(1.0 - raddupdif2(i,j,k)*radddowndifm(i,j))
            dm3r2 = dm32*radddowndifm(i,j)      
            dm3r1p2 = 1.0 + raddupdif2(i,j,k)*dm3r2   
            alpp2 = (tadddowndirm(i,j) - tlevel2p(i,j))*dm32   
            if ( calc_flag(i,j) ) then
              transmittance(i,j,k) = (tlevel2p(i,j)*(1.0 + raddupdir2(i,j,k)* &
                                     dm3r2) + alpp2)
              tr_dir(i,j,k) = tlevel2p(i,j)
              reflectance(i,j,k) = (tlevel2p(i,j)*raddupdir2(i,j,k)*   &
                                    dm3r1p2 + alpp2*   &
                                    raddupdif2(i,j,k))
            end if
            tlevel2p(i,j) = tlevel2p(i,j)*tlayerde (i,j,k) 
            radddowndifm(i,j) = radddowndif2(i,j,k)
            tadddowndirm(i,j) = tadddowndir2(i,j,k)
          end do
        end do
      end do
      do j = 1,size(calc_flag,2)
        do i = 1,size(calc_flag,1)
!! CORRECT ???
!         dm32  = 1.0/(1.0 - sfcalb(i,j)*radddowndifm(i,j))
          dm32          = 1.0/(1.0 - sfcalb_dif(i,j)*   &
                          radddowndifm(i,j)       )
          dm3r2 = dm32*radddowndifm(i,j)       
!! CORRECT ???
!         dm3r1p2 = 1.0 + sfcalb(i,j)*dm3r2         
          dm3r1p2          = 1.0 + sfcalb_dif(i,j) * dm3r2
          alpp2 = (tadddowndirm(i,j) - tlevel2p(i,j))*dm32
          if ( calc_flag(i,j) ) then
            transmittance(i,j,kx+1) = (tlevel2p(i,j)*(1.0 +   &
!! CORRECT ???
!                                 sfcalb(i,j)* &
!12-08-03:  CHANGE THIS TO _dir as per SMF  sfcalb_dif(i,j)* &
                                      sfcalb_dir(i,j)* &
                                      dm3r2) + alpp2)
            tr_dir(i,j,kx+1) = tlevel2p(i,j)
            reflectance(i,j,kx+1) = (tlevel2p(i,j)*  &
!! CORRECT ???
!                                   sfcalb(i,j)*   &
                                    sfcalb_dir(i,j)* &
                                    dm3r1p2 + alpp2* &
!! CORRECT ???
!                                   sfcalb(i,j) )
                                    sfcalb_dif(i,j))  
            reflectance(i,j,1) = raddupdir2p(i,j)         
            transmittance(i,j,1) = 1.0
            tr_dir(i,j,1) = 1.0
          end if
        end do
      end do

!------------------------------------------------------------------


end subroutine adding 



!####################################################################
! <SUBROUTINE NAME="deledd">
!  <OVERVIEW>
!   Subroutine that calculates reflectivity and transmissivity in a
!   scattering layer using delta-eddington method
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes layer optical depth, single scattering abledo,
!   and asymmetry parameter, using delta-eddington method, to calculate
!   direct/diffuse reflectivity/transmissivity to direct/diffuse incident
!   radiation. The approximation uses the strong forward scattering of
!   aerosol particles.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deledd (ix, jx, kx,  &
!                taustr, omegastr, gstr, cosang, ng , daylight,  &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,   &
!                tlayerde,  cloud)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   jx is the current latitudinal index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   ix is the current vertical index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="taustr" TYPE="real">
!   the scaled optical depth, true optical depth normalized using
!   delta-eddington approximation
!  </IN>
!  <IN NAME="omegastr" TYPE="real">
!   the scaled single-scattering albedo
!  </IN>
!  <IN NAME="gstr" TYPE="real">
!   the scaled asymmetry factor
!  </IN>
!  <IN NAME="cosang" TYPE="real">
!   cosine of the solar zenith angle
!  </IN>
!  <IN NAME="ng" TYPE="real">
!   the number of gaussian angles to compute the diurnally    
!   averaged solar radiation (=1 unless lswg = true)
!  </IN>
!  <IN NAME="cloud" TYPE="real">
!   flag for existence of a cloud (used only in 'ovc' mode)
!  </IN>
!  <OUT NAME="rlayerdir" TYPE="real">
!   layer reflectivity to direct incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to direct incident beam
!  </OUT>
!  <OUT NAME="rlayerdif" TYPE="real">
!   layer reflectivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerde" TYPE="real">
!   layer diffuse transmissivity to direct incident beam
!  </OUT>
! </SUBROUTINE>
!
subroutine deledd (ix, jx, kx, taustr, omegastr, gstr, cosang, ng, &
                   daylight, rlayerdir, tlayerdir, tlayerde,   &
                   rlayerdif, tlayerdif, cloud)

!---------------------------------------------------------------------- 
!    deledd calculates the reflection and transmission in the 
!    scattering layers using the delta-eddington method.         
!    references:                                                   
!      joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington
!      approximation for radiative flux transfer.,j. atmos. sci.,33,  
!      2452-2459, 1976.                                              
!-------------------------------------------------------------------

integer,                   intent(in)              :: ix, jx, kx
real, dimension(:,:,:),    intent(inout)           :: taustr, omegastr
real, dimension(:,:,:),    intent(in)              :: gstr
real, dimension(:,:),    intent(in)                ::  cosang
integer,                   intent(in)              :: ng
logical, dimension(:,:),   intent(in)              :: daylight
real, dimension(:,:,:),    intent(out)             :: rlayerdir,   &
                                                      tlayerdir,   &
                                                      tlayerde
real, dimension(:,:,:),    intent(inout), optional :: rlayerdif,   &
                                                      tlayerdif
logical, dimension(:,:,:), intent(in), optional    :: cloud         

!----------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx
!    gstr        the scaled asymmetry factor                       
!    cosang      the cosine of the solar zenith angle    
!    ng          the number of gaussian angles to compute the diurnally 
!                averaged solar radiation (=1 unless lswg = true)       
!    daylight
!
!  intent(inout) variables:
!
!    taustr      the scaled extinction optical depth                    
!    omegastr    the scaled single-scattering albedo               
!
!  intent(out) variables:
!
!    rlayerdir   the layer reflectivity to a direct incident beam      
!    tlayerdir   the layer transmissivity to a direct incident beam   
!    rlayerdif   the layer reflectivity to a diffuse incident beam   
!    tlayerdif   the layer transmissivity to a diffuse incident beam
!    tlayerde    the layer transmissivity (non-scattered) to the direct 
!                incident beam                                       
!
! intent(in),optional:
!
!    cloud       flag for existence of a cloud (used only in 'ovc' mode)
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real        :: qq(7), rr(5), ss(8), tt(8), ww(4)
      real        :: rsum, tsum
      real        :: tmp
      real, parameter :: onedi3 = 1.0/3.0           
      real, parameter :: twodi3 = 2.0/3.0             
      integer     :: k, ns, j, nn, ntot

      real,    dimension(ix)                  ::   &
                                          gstr2, taustr2, omegastr2, &
                                           cosangzk2, rlayerdir2,    &
                                           tlayerde2, tlayerdir2, &
                                           sumr, sumt


!----------------------------------------------------------------------
!  local variables:
!
!      qq
!      rr
!      ss
!      tt
!      ww
!      rsum
!      tsum
!      alpha
!      onedi3
!      twodi3
!      i,j,k
!      ns
!      nn
!      ntot
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=1,kx         
        do j=1,jx         

!---------------------------------------------------------------------
!    overcast sky mode. note: in this mode, the delta-eddington method
!    is performed only for spatial points containing a cloud.   
!-------------------------------------------------------------------
          if (present(cloud)) then
            nn = count(cloud(1:ix,j,k))
            gstr2(1:nn) = pack( gstr(1:ix,j,k), cloud(1:ix,j,k) )
            taustr2(1:nn) = pack( taustr(1:ix,j,k), cloud(1:ix,j,k) )
            omegastr2(1:nn) = pack( omegastr(1:ix,j,k), cloud(1:ix,j,k) )
            cosangzk2(1:nn) = pack( cosang(1:ix,j), cloud(1:ix,j,k) )
!----------------------------------------------------------------------
!    note: the following are done to avoid the conservative scattering 
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                      
!----------------------------------------------------------------------c
            omegastr2(1:nn) = min( omegastr2(1:nn), 9.9999999E-01 )
            taustr2(1:nn) = min( taustr2(1:nn), 1.0E+02 )

!----------------------------------------------------------------------c
!    clear sky mode. note: in this mode, the delta-eddington method is 
!    performed for all spatial points.                 
!----------------------------------------------------------------------c
          else
            nn = count(daylight(1:ix,j))
            gstr2(1:nn) = pack( gstr(1:ix,j,k), daylight(1:ix,j) )
            taustr2(1:nn) = pack( taustr(1:ix,j,k), daylight(1:ix,j) )
            omegastr2(1:nn) = pack( omegastr(1:ix,j,k), daylight(1:ix,j) )
            cosangzk2(1:nn) = pack( cosang(1:ix,j), daylight(1:ix,j) )            
!----------------------------------------------------------------------c
!    note: the following are done to avoid the conservative scattering  
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                    
!----------------------------------------------------------------------c
            omegastr2(1:nn) = min( omegastr2(1:nn), 9.9999999E-01 )
            taustr2(1:nn) = min( taustr2(1:nn), 1.0E+02 )
          endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          ntot = nn

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          if (present (tlayerdif) .and. present(rlayerdif) .and. ng==1 ) then

            if ( nstreams==1 ) then

              do nn=1,ntot      

!----------------------------------------------------------------------c
!    direct quantities                                            
!----------------------------------------------------------------------c
                ww(1) = omegastr2(nn)
                ww(2) = gstr2(nn)
                ww(3) = taustr2(nn)
                ww(4) = cosangzk2(nn)

                qq(1) = 3.0 * ( 1.0 - ww(1) )
                qq(2) = 1.0 - ww(1) * ww(2)
                qq(3) = qq(1)/qq(2)
                qq(4) = sqrt( qq(1) * qq(2) )
                qq(5) = sqrt (qq(3))
                qq(6) = 1.0 + twodi3 * qq(5)         
                qq(7) = 1.0 - twodi3 * qq(5)       

                rr(1) = 1./qq(6)
                rr(2) = qq(7)*rr(1)
                rr(3) = exp( -ww(3) * qq(4) )
                rr(4) = 1.0/rr(3)
                rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) * rr(2) )

                tmp   = 1.0 - (qq(4)*ww(4)) ** 2
                if ( abs(tmp)<1.e-20 ) tmp=1.e-20 ! MJT suggestion
                ss(1) = 0.75 * ww(1)/tmp
                ss(2) = ss(1)*ww(4)*( 1.0 + ww(2)*qq(1)*onedi3)
                ss(3) = ss(1)*(1.0 + ww(2)*qq(1)*ww(4)** 2 )
                ss(4) = ss(2) - twodi3*ss(3)     
                ss(5) = ss(2) + twodi3*ss(3)     
                ss(6) = exp( -ww(3) / ww(4) )
                ss(7) = (ss(4)*ss(6) - ss(5)*rr(3)*rr(2))*rr(5)
                ss(8) = (ss(5) - qq(7)*ss(7))*rr(1)
                
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
                rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - ss(4)
                tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                                   qq(7) * rr(4) * ss(7) -  &
                                   ss(5) * ss(6) ) + ss(6) )
                tlayerde2(nn) = ss(6)

!----------------------------------------------------------------------c
!    diffuse quantities                                       
!    notes: the number of streams for the diffuse beam is fixed at 4.   
!    this calculation is done only for ng=1.                 
!----------------------------------------------------------------------c
 
                tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                        cosangstr(1) ) ** 2 )
                tt(2) = tt(1) * cosangstr(1) * ( 1.0 +  &
                        ww(2)        * qq(1) * onedi3 )
                tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                        cosangstr(1) ** 2 )
                tt(4) = tt(2) - twodi3 * tt(3)
                tt(5) = tt(2) + twodi3 * tt(3)
                tt(6) = exp( -ww(3)          / cosangstr(1) )
                tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                        rr(3) * rr(2)   ) * rr(5)
                tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)
                sumr(nn) = (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))
                sumt(nn) = ( (rr(3)*qq(6)*tt(8) +    &
                           qq(7)*rr(4)*tt(7) -   &
                           tt(5)*tt(6)) + tt(6))
              end do  ! ntot loop
              
            else    
                
              do nn=1,ntot      

!----------------------------------------------------------------------c
!    direct quantities                                            
!----------------------------------------------------------------------c
                ww(1) = omegastr2(nn)
                ww(2) = gstr2(nn)
                ww(3) = taustr2(nn)
                ww(4) = cosangzk2(nn)

                qq(1) = 3.0 * ( 1.0 - ww(1) )
                qq(2) = 1.0 - ww(1) * ww(2)
                qq(3) = qq(1)/qq(2)
                qq(4) = sqrt( qq(1) * qq(2) )
                qq(5) = sqrt (qq(3))
                qq(6) = 1.0 + twodi3 * qq(5)         
                qq(7) = 1.0 - twodi3 * qq(5)       

                rr(1) = 1./qq(6)
                rr(2) = qq(7)*rr(1)
                rr(3) = exp( -ww(3)          * qq(4) )
                rr(4) = 1.0/rr(3)
                rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) * rr(2) )

                tmp   = 1.0 - (qq(4)*ww(4)) ** 2
                if ( abs(tmp)<1.e-20 ) tmp=1.e-20 ! MJT suggestion
                ss(1) = 0.75 * ww(1)/tmp
                ss(2) = ss(1)*ww(4)*( 1.0 + ww(2)*qq(1)*onedi3)
                ss(3) = ss(1)*(1.0 + ww(2)*qq(1)*ww(4)** 2 )
                ss(4) = ss(2) - twodi3*ss(3)     
                ss(5) = ss(2) + twodi3*ss(3)     
                ss(6) = exp( -ww(3) / ww(4) )
                ss(7) = (ss(4)*ss(6) - ss(5)*rr(3)*rr(2))*rr(5)
                ss(8) = (ss(5) - qq(7)*ss(7))*rr(1)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
                rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - ss(4)
                tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                                   qq(7) * rr(4) * ss(7) -  &
                                   ss(5) * ss(6) ) + ss(6) )
                tlayerde2(nn) = ss(6)

!----------------------------------------------------------------------c
!    diffuse quantities                                       
!    notes: the number of streams for the diffuse beam is fixed at 4.   
!    this calculation is done only for ng=1.                 
!----------------------------------------------------------------------c
                rsum = 0.0
                tsum = 0.0
                do ns = 1,NSTREAMS
                  tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                          cosangstr(ns) ) ** 2 )
                  tt(2) = tt(1) * cosangstr(ns) * ( 1.0 +  &
                          ww(2)        * qq(1) * onedi3 )
                  tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                          cosangstr(ns) ** 2 )
                  tt(4) = tt(2) - twodi3 * tt(3)
                  tt(5) = tt(2) + twodi3 * tt(3)
                  tt(6) = exp( -ww(3)          / cosangstr(ns) )
                  tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                          rr(3) * rr(2)   ) * rr(5)
                  tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)
                  rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))* &
                         wtstr(ns)*cosangstr(ns)
                  tsum = tsum + ((rr(3)*qq(6)*tt(8) +   &
                                  qq(7)*rr(4)*tt(7) -   &
                                  tt(5)*tt(6)) + tt(6))*  &
                                  wtstr(ns)*cosangstr(ns)
                end do
                sumr(nn) = rsum
                sumt(nn) = tsum
              
              end do  ! ntot loop
              
            end if
             
          else
              
            do nn=1,ntot      

!----------------------------------------------------------------------c
!    direct quantities                                            
!----------------------------------------------------------------------c
              ww(1) = omegastr2(nn)
              ww(2) = gstr2(nn)
              ww(3) = taustr2(nn)
              ww(4) = cosangzk2(nn)

              qq(1) = 3.0 * ( 1.0 - ww(1) )
              qq(2) = 1.0 - ww(1) * ww(2)
              qq(3) = qq(1)/qq(2)
              qq(4) = sqrt( qq(1) * qq(2) )
              qq(5) = sqrt( qq(3) )
              qq(6) = 1.0 + twodi3 * qq(5)         
              qq(7) = 1.0 - twodi3 * qq(5)       

              rr(1) = 1./qq(6)
              rr(2) = qq(7)*rr(1)
              rr(3) = exp( -ww(3) * qq(4) )
              rr(4) = 1.0/rr(3)
              rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) * rr(2) )

              tmp   = 1.0 - (qq(4)*ww(4)) ** 2
              if ( abs(tmp)<1.e-20 ) tmp=1.e-20 ! MJT suggestion
              ss(1) = 0.75 * ww(1)/tmp
              ss(2) = ss(1)*ww(4)*( 1.0 + ww(2)*qq(1)*onedi3)
              ss(3) = ss(1)*(1.0 + ww(2)*qq(1)*ww(4)** 2 )
              ss(4) = ss(2) - twodi3*ss(3)     
              ss(5) = ss(2) + twodi3*ss(3)     
              ss(6) = exp( -ww(3) / ww(4) )
              ss(7) = (ss(4)*ss(6) - ss(5)*rr(3)*rr(2))*rr(5)
              ss(8) = (ss(5) - qq(7)*ss(7))*rr(1)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
              rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - ss(4)
              tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                                 qq(7) * rr(4) * ss(7) -  &
                                 ss(5) * ss(6) ) + ss(6) )
              tlayerde2(nn) = ss(6)

            end do  ! ntot loop

          end if ! present(tlayerdiff)..
!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------
          if ( present(cloud) ) then
            rlayerdir(1:ix,j,k) = unpack( rlayerdir2(1:ntot), cloud(1:ix,j,k), 0. )
            tlayerdir(1:ix,j,k) = unpack( tlayerdir2(1:ntot), cloud(1:ix,j,k), 0. )
            tlayerde(1:ix,j,k) = unpack( tlayerde2(1:ntot), cloud(1:ix,j,k), 0. )
            if ( present(tlayerdif) .and. ng==1 ) then
              rlayerdif(1:ix,j,k) = unpack( sumr(1:ntot), cloud(1:ix,j,k), 0. )
              tlayerdif(1:ix,j,k) = unpack( sumt(1:ntot), cloud(1:ix,j,k), 0. )
            end if
          else
            rlayerdir(1:ix,j,k) = unpack( rlayerdir2(1:ntot), daylight(1:ix,j), 0. )
            tlayerdir(1:ix,j,k) = unpack( tlayerdir2(1:ntot), daylight(1:ix,j), 0. )
            tlayerde(1:ix,j,k) = unpack( tlayerde2(1:ntot), daylight(1:ix,j), 0. )
            if ( present(tlayerdif) .and. ng==1 ) then
              rlayerdif(1:ix,j,k) = unpack( sumr(1:ntot), daylight(1:ix,j), 0. )
              tlayerdif(1:ix,j,k) = unpack( sumt(1:ntot), daylight(1:ix,j), 0. )
            end if
          end if

        end do
      end do

!---------------------------------------------------------------------
 
end subroutine deledd



!#####################################################################


                   end module esfsw_driver_mod
