! Interface for SEA-ESF radiation scheme (from GFDL AM3) with CCAM.

! Interface developed by Marcus Thatcher and Maciej Golebiewski

! - This routine assumes that only one month at a time is integrated in RCM mode

module seaesfrad_m

use rad_utilities_mod, only : atmos_input_type,surface_type,astronomy_type,aerosol_type,           &
                              aerosol_properties_type,radiative_gases_type,cldrad_properties_type, &
                              cld_specification_type,lw_output_type,sw_output_type,                &
                              aerosol_diagnostics_type,time_type,microphysics_type,                &
                              microrad_properties_type,lw_diagnostics_type,lw_table_type,          &
                              Sw_control,Lw_control, Rad_control,Cldrad_control,Lw_parameters,     &
                              thickavg
use esfsw_driver_mod, only : swresf,esfsw_driver_init
use sealw99_mod, only : sealw99,sealw99_init, sealw99_time_vary
use esfsw_parameters_mod, only : Solar_spect,esfsw_parameters_init

private
public seaesfrad

real, parameter :: cp       = 1004.64     ! Specific heat of dry air at const P
real, parameter :: grav     = 9.80616     ! Acceleration due to gravity
real, parameter :: stefbo   = 5.67e-8     ! Stefan-Boltzmann constant
real, parameter :: rdry     = 287.04      ! Gas constant for dry air
real, parameter :: rhow     = 1000.       ! Density of water
real, parameter :: lv       = 2.5104e6
real, parameter :: lf       = 3.36e5
real, parameter :: ls       = lv+lf
real, parameter :: pi       = 3.1415927   ! pi
real, parameter :: csolar   = 1365        ! Solar constant in W/m^2
real, parameter :: siglow   = 0.68        ! sigma level for top of low cloud (diagnostic)
real, parameter :: sigmid   = 0.44        ! sigma level for top of medium cloud (diagnostic)
real, parameter :: ratco2mw = 1.519449738 ! conversion factor for CO2 diagnostic
integer, parameter :: naermodels         = 900
integer, parameter :: N_AEROSOL_BANDS_FR = 8
integer, parameter :: N_AEROSOL_BANDS_CO = 1
integer, parameter :: N_AEROSOL_BANDS_CN = 1
integer, parameter :: N_AEROSOL_BANDS    = N_AEROSOL_BANDS_FR+N_AEROSOL_BANDS_CO
integer, parameter :: nfields            = 10
logical, parameter :: do_totcld_forcing  = .true.
logical, parameter :: include_volcanoes  = .false.
logical, save :: do_aerosol_forcing

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CCAM interface with GFDL SEA-ESF radiation
!

subroutine seaesfrad(imax,odcalc)

use aerointerface                                   ! Aerosol interface
use aerosolldr                                      ! LDR prognostic aerosols
use arrays_m                                        ! Atmosphere dyamics prognostic arrays
use ateb                                            ! Urban
use cc_mpi                                          ! CC MPI routines
use cfrac_m                                         ! Cloud fraction
use extraout_m                                      ! Additional diagnostics
use histave_m, only : alb_ave,fbeam_ave             ! Time average arrays
use infile                                          ! Input file routines
use latlong_m                                       ! Lat/lon coordinates
use microphys_rad_mod                               ! SEA/ESF microphysics
use mlo                                             ! Ocean physics and prognostic arrays
use nharrs_m                                        ! Non-hydrostatic atmosphere arrays
use nsibd_m                                         ! Land-surface arrays
use ozoneread                                       ! Ozone input routines
use pbl_m                                           ! Boundary layer arrays
use raddiag_m                                       ! Radiation diagnostic
use radisw_m, only : rrco2,rrvco2,rrvch4,rrvn2o, &  ! GHG data
    rrvf11,rrvf12,rrvf113,rrvf22
use sigs_m                                          ! Atmosphere sigma levels
use soil_m                                          ! Soil and surface data
use soilsnow_m                                      ! Soil, snow and surface data
use work3f_m                                        ! Grid work arrays
use zenith_m                                        ! Astronomy routines

implicit none

include 'parm.h'
include 'newmpar.h'
include 'kuocom.h'

logical, intent(in) :: odcalc  ! True for full radiation calculation
integer, intent(in) :: imax
integer jyear,jmonth,jday,jhour,jmin
integer k,ksigtop,mins
integer i,j,iq,istart,iend,kr,nr
integer ierr,ktop,kbot
integer, save :: nlow,nmid
real, dimension(:), allocatable, save :: sgamp
real, dimension(:,:), allocatable, save :: rtt
real, dimension(imax,kl) :: duo3n,rhoa
real, dimension(imax,kl) :: p2,cd2,dumcf,dumql,dumqf,dumt,dz
real, dimension(imax) :: qsat,coszro2,taudar2,coszro,taudar,mx
real, dimension(imax) :: sg,sint,sout,sgdn,rg,rt,rgdn,sgdnvis,sgdnnir
real, dimension(imax) :: soutclr,sgclr,rtclr,rgclr,sga
real, dimension(imax) :: sgvis,sgdnvisdir,sgdnvisdif,sgdnnirdir,sgdnnirdif
real, dimension(imax) :: dzrho,dumfbeam,tv,tnhs
real, dimension(imax) :: cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif
real, dimension(kl+1) :: sigh
real(kind=8), dimension(:,:), allocatable, save :: pref
real r1,dlt,alp,slag,dhr,fjd
real ttbg,ar1,exp_ar1,ar2,exp_ar2,ar3,snr
real dnsnow,snrat,dtau,alvo,aliro,fage,cczen,fzen,fzenm
real alvd,alv,alird,alir
real f1,f2,cosz,delta
logical, save :: first = .true.

type(time_type), save ::                    Rad_time
type(atmos_input_type), save ::             Atmos_input
type(surface_type), save ::                 Surface     
type(astronomy_type), save ::               Astro
type(aerosol_type), save ::                 Aerosol
type(aerosol_properties_type), save ::      Aerosol_props
type(radiative_gases_type), save ::         Rad_gases
type(cldrad_properties_type), save ::       Cldrad_props
type(cld_specification_type), save ::       Cld_spec
type(microphysics_type), save ::            Cloud_microphysics
type(microrad_properties_type), save ::     Lscrad_props
type(lw_output_type), dimension(1), save :: Lw_output
type(sw_output_type), dimension(1), save :: Sw_output
type(aerosol_diagnostics_type), save ::     Aerosol_diags
type(lw_table_type), save ::                Lw_tables
real(kind=8), dimension(:,:,:,:), allocatable, save :: r

call START_LOG(radmisc_begin)

if ( nmaxpr==1 .and. myid==0 ) then
  write(6,*) "seaesfrad: Starting SEA-ESF radiation"
end if

! Aerosol flag
do_aerosol_forcing = abs(iaero)>=2

! set-up half levels ------------------------------------------------
sigh(1:kl) = sigmh(1:kl)
sigh(kl+1) = 0.

! astronomy ---------------------------------------------------------
! Set up number of minutes from beginning of year
call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
fjd = float(mod(mins,525600))/1440. ! restrict to 365 day calendar

! Calculate sun position
call solargh(fjd,bpyear,r1,dlt,alp,slag)

! Prepare SEA-ESF arrays --------------------------------------------
Rad_time%days   =int(fjd)
Rad_time%seconds=mod(mins,1440)
Rad_time%ticks  =0

! Initialisation ----------------------------------------------------
if ( first ) then
  first = .false.

  if ( myid==0 ) write(6,*) "Initalising SEA-ESF radiation"
  allocate(sgamp(ifull),rtt(ifull,kl))

  ! initialise co2
  call co2_read(sig,jyear)
  rrco2=rrvco2*ratco2mw

  ! initialise ozone
  if ( amipo3 ) then
    if ( myid==0 ) write(6,*) 'AMIP2 ozone input'
    call o3read_amip
  else
    call o3_read(sig,jyear,jmonth)
  end if

  ! set-up standard pressure levels
  allocate(pref(kl+1,2))
  pref(kl+1,1)=101325.
  pref(kl+1,2)=81060. !=0.8*pref(kl+1,1)
  do k=1,kl
    kr=kl+1-k
    pref(kr,:)=sig(k)*pref(kl+1,:)
  end do  
  
  Cldrad_control%do_strat_clouds_iz      =.true.
  Cldrad_control%do_sw_micro_iz          =.true.
  Cldrad_control%do_lw_micro_iz          =.true.
  Cldrad_control%do_sw_micro             =.true.
  Cldrad_control%do_lw_micro             =.true.
  Cldrad_control%do_ica_calcs            =.false. ! must change allocations below if true
  Cldrad_control%do_no_clouds            =.false.
  Cldrad_control%do_donner_deep_clouds   =.false.
  Cldrad_control%do_stochastic_clouds    =.false.
  Cldrad_control%using_fu2007            =.false.
  Sw_control%solar_constant              =csolar
  Sw_control%do_cmip_diagnostics         =do_aerosol_forcing ! Need for aerosol optical depths
  Lw_control%do_lwcldemiss               =.true.
  Lw_control%do_o3_iz                    =.true.
  Lw_control%do_co2_iz                   =.true.
  Lw_control%do_ch4_iz                   =.true.
  Lw_control%do_n2o_iz                   =.true.
  Lw_control%do_o3                       =.true.
  Lw_control%do_co2                      =.true.
  Lw_control%do_ch4                      =rrvch4>0.
  Lw_control%do_n2o                      =rrvch4>0.
  Lw_control%do_h2o                      =.true.
  Lw_control%do_cfc                      =rrvch4>0.
  Rad_control%using_solar_timeseries_data=.false.
  Rad_control%do_totcld_forcing          =do_totcld_forcing
  Rad_control%rad_time_step              =kountr*dt
  Rad_control%rad_time_step_iz           =.true.
  Rad_control%do_aerosol                 =do_aerosol_forcing
  Rad_control%do_swaerosol_forcing       =do_aerosol_forcing
  Rad_control%do_lwaerosol_forcing       =do_aerosol_forcing
  Rad_control%hires_coszen               =.false.
  Rad_control%hires_coszen_iz            =.true.
  Rad_control%nzens                      =1
  Rad_control%nzens_iz                   =.true.
  Rad_control%using_im_bcsul             =.false.
  Rad_control%using_im_bcsul_iz          =.true.
  Astro%rrsun                            =1./(r1*r1)

  call sealw99_init(pref, Lw_tables)
  call esfsw_parameters_init
  call esfsw_driver_init
  call microphys_rad_init

  deallocate(pref)
  
  allocate ( Atmos_input%press(imax, 1, kl+1) )
  allocate ( Atmos_input%phalf(imax, 1, kl+1) )
  allocate ( Atmos_input%temp(imax, 1, kl+1) )
  allocate ( Atmos_input%rh2o(imax, 1, kl) )
  allocate ( Atmos_input%rel_hum(imax, 1, kl) )
  allocate ( Atmos_input%clouddeltaz(imax, 1,kl) )
  allocate ( Atmos_input%deltaz(imax, 1, kl) )
  allocate ( Atmos_input%pflux(imax, 1, kl+1) )
  allocate ( Atmos_input%tflux(imax, 1, kl+1) )
  allocate ( Atmos_input%psfc(imax, 1) )
  allocate ( Atmos_input%tsfc(imax, 1) )
  !if (use_co2_tracer_field) then
  !  allocate ( Atmos_input%tracer_co2(imax, 1, kl) )
  !endif
  
  allocate( Rad_gases%qo3(imax, 1, kl) )

  allocate( Cloud_microphysics%size_rain(imax, 1, kl) )
  allocate( Cloud_microphysics%size_drop(imax, 1, kl) )
  allocate( Cloud_microphysics%size_ice(imax, 1, kl) )
  allocate( Cloud_microphysics%size_snow(imax, 1, kl) )
  allocate( Cloud_microphysics%conc_drop(imax, 1, kl) )
  allocate( Cloud_microphysics%conc_ice(imax, 1, kl) )
  allocate( Cloud_microphysics%conc_rain(imax, 1, kl) )
  allocate( Cloud_microphysics%conc_snow(imax, 1, kl) )

  allocate( Cldrad_props%cldext(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props%cldsct(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props%cldemiss(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props%emmxolw(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props%emrndlw(imax, 1, kl, Cldrad_control%nlwcldb,1) )

  allocate( Lscrad_props%cldext(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props%cldsct(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb) )

  allocate( Cld_spec%camtsw(imax, 1, kl) )
  allocate( Cld_spec%cmxolw(imax, 1, kl) )
  allocate( Cld_spec%crndlw(imax, 1, kl) )
  
  allocate( Surface%asfc_vis_dir(imax, 1) )
  allocate( Surface%asfc_nir_dir(imax, 1) )
  allocate( Surface%asfc_vis_dif(imax, 1) )
  allocate( Surface%asfc_nir_dif(imax, 1) )

  allocate( Astro%cosz(imax, 1) )
  allocate( Astro%fracday(imax, 1) )

  allocate( Lw_output(1)%heatra(imax, 1, kl) )
  allocate( Lw_output(1)%flxnet(imax, 1, kl+1) )
  allocate( Lw_output(1)%bdy_flx(imax, 1, 4) )
  if (do_totcld_forcing) then
    allocate ( Lw_output(1)%heatracf(imax, 1, kl) )
    allocate ( Lw_output(1)%flxnetcf(imax, 1, kl+1) )
    allocate ( Lw_output(1)%bdy_flx_clr(imax, 1, 4) )
  endif

  allocate( Sw_output(1)%dfsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(1)%ufsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(1)%dfsw_dir_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%dfsw_dif_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%ufsw_dir_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%ufsw_dif_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%fsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(1)%hsw(imax, 1, kl, Rad_control%nzens) )
  allocate( Sw_output(1)%dfsw_vis_sfc(imax, 1, Rad_control%nzens ) )
  allocate( Sw_output(1)%ufsw_vis_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%dfsw_vis_sfc_dir(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%dfsw_vis_sfc_dif(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%ufsw_vis_sfc_dir(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%ufsw_vis_sfc_dif(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(1)%bdy_flx(imax, 1, 4, Rad_control%nzens) )
  if (do_totcld_forcing) then
    allocate( Sw_output(1)%dfswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(1)%ufswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(1)%fswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(1)%hswcf(imax, 1, kl, Rad_control%nzens) )
    allocate( Sw_output(1)%dfsw_dir_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(1)%dfsw_dif_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(1)%dfsw_vis_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(1)%bdy_flx_clr(imax, 1, 4, Rad_control%nzens) )
  endif

  if (do_aerosol_forcing) then
    if ( Rad_control%using_im_bcsul ) then
      allocate( Aerosol_props%sulfate_index(0:100, 0:100) )
    else
      allocate( Aerosol_props%sulfate_index(0:100, 0:0) )
    end if
    allocate( Aerosol_props%omphilic_index(0:100) )
    allocate( Aerosol_props%bcphilic_index(0:100) )
    allocate( Aerosol_props%seasalt1_index(0:100) )
    allocate( Aerosol_props%seasalt2_index(0:100) )
    allocate( Aerosol_props%seasalt3_index(0:100) )
    allocate( Aerosol_props%seasalt4_index(0:100) )
    allocate( Aerosol_props%seasalt5_index(0:100) )
    allocate( Aerosol_props%optical_index(nfields) )
    allocate( Aerosol%aerosol(imax, 1, kl, nfields) )
    allocate( Atmos_input%aerosolrelhum(imax, 1, kl) )
    allocate( Aerosol_props%ivol(imax, 1, kl) )
    allocate( Aerosol_props%aerextband(Solar_spect%nbands, naermodels) )
    allocate( Aerosol_props%aerssalbband(Solar_spect%nbands, naermodels) )
    allocate( Aerosol_props%aerasymmband(Solar_spect%nbands, naermodels) )
    allocate( Aerosol_props%aerssalbbandlw(N_AEROSOL_BANDS, naermodels) )
    allocate( Aerosol_props%aerextbandlw(N_AEROSOL_BANDS, naermodels) )
    allocate( Aerosol_props%aerssalbbandlw_cn(N_AEROSOL_BANDS, naermodels) )
    allocate( Aerosol_props%aerextbandlw_cn(N_AEROSOL_BANDS, naermodels) )
    allocate( Aerosol_diags%extopdep(imax, 1, kl, nfields, 10) )
    allocate( Aerosol_diags%absopdep(imax, 1, kl, nfields, 10) )
    allocate( Aerosol_diags%asymdep(imax, 1, kl, nfields, 10) )
    
    Aerosol_props%sulfate_flag =0
    Aerosol_props%omphilic_flag=-1
    Aerosol_props%bcphilic_flag=-2
    Aerosol_props%seasalt1_flag=-3
    Aerosol_props%seasalt2_flag=-4
    Aerosol_props%seasalt3_flag=-5
    Aerosol_props%seasalt4_flag=-6
    Aerosol_props%seasalt5_flag=-7
    Aerosol_props%bc_flag      =-8
    Lw_parameters%n_lwaerosol_bands=N_AEROSOL_BANDS
    Aerosol_props%optical_index(1)=Aerosol_props%sulfate_flag     ! so4
    if ( Rad_control%using_im_bcsul ) then
      Aerosol_props%optical_index(2)=Aerosol_props%bc_flag        ! so4/bc mixture
      Aerosol_props%optical_index(3)=Aerosol_props%bc_flag        ! so4/bc mixture
    else
      Aerosol_props%optical_index(2)=719                          ! black carbon (hydrophobic)
      Aerosol_props%optical_index(3)=Aerosol_props%bcphilic_flag  ! black carbon (hydrophillic)
    end if
    Aerosol_props%optical_index(4)=720                            ! organic carbon (hydrophobic)
    Aerosol_props%optical_index(5)=Aerosol_props%omphilic_flag    ! organic carbon (hydrophillic)
    Aerosol_props%optical_index(6)=897                            ! dust_0.7  (using 0.73)
    Aerosol_props%optical_index(7)=898                            ! dust_1.4  (using 1.4)
    Aerosol_props%optical_index(8)=899                            ! dust_2.4  (using 2.4)
    Aerosol_props%optical_index(9)=900                            ! dust_4.5  (using 4.5)
    Aerosol_props%optical_index(10)=705                           ! sea_salt (film drop + jet drop)
    ! GFDL bins dust1=0.1-0.5, dust2=0.5-1, dust3=1-2.5, dust4=2.5-5, dust5=5-10
    ! GFDL bins salt1=0.1-0.5, salt2=0.5-1, salt3=1-2.5, salt4=2.5-5, dust5=5-10

    Aerosol_props%sulfate_index( 0: 13, 0)=(/  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
    Aerosol_props%sulfate_index(14: 27, 0)=(/  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
    Aerosol_props%sulfate_index(28: 41, 0)=(/  1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3 /)
    Aerosol_props%sulfate_index(42: 55, 0)=(/  3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6 /)
    Aerosol_props%sulfate_index(56: 69, 0)=(/  6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9 /)
    Aerosol_props%sulfate_index(70: 83, 0)=(/  9, 9, 9,10,10,10,10,10,11,11,11,11,12,12 /)
    Aerosol_props%sulfate_index(84: 97, 0)=(/ 13,13,14,14,15,15,16,17,18,19,20,21,22,23 /)
    Aerosol_props%sulfate_index(98:100, 0)=(/ 23,23,23 /)
    if ( Rad_control%using_im_bcsul ) then
      Aerosol_props%sulfate_index(:, 1)=Aerosol_props%sulfate_index(:, 0)      
      Aerosol_props%sulfate_index(:, 2)=Aerosol_props%sulfate_index(:, 1) + 26 ! 98%
      Aerosol_props%sulfate_index(:, 3)=Aerosol_props%sulfate_index(:, 2)
      Aerosol_props%sulfate_index(:, 4)=Aerosol_props%sulfate_index(:, 3) + 26 ! 96%
      Aerosol_props%sulfate_index(:, 5)=Aerosol_props%sulfate_index(:, 4)
      Aerosol_props%sulfate_index(:, 6)=Aerosol_props%sulfate_index(:, 5) + 26 ! 94%
      Aerosol_props%sulfate_index(:, 7)=Aerosol_props%sulfate_index(:, 6)
      Aerosol_props%sulfate_index(:, 8)=Aerosol_props%sulfate_index(:, 7) + 26 ! 92%
      Aerosol_props%sulfate_index(:, 9)=Aerosol_props%sulfate_index(:, 8)
      Aerosol_props%sulfate_index(:,10)=Aerosol_props%sulfate_index(:, 9) + 26 ! 90%
      Aerosol_props%sulfate_index(:,11)=Aerosol_props%sulfate_index(:,10)
      Aerosol_props%sulfate_index(:,12)=Aerosol_props%sulfate_index(:,11) + 26 ! 88%
      Aerosol_props%sulfate_index(:,13)=Aerosol_props%sulfate_index(:,12)
      Aerosol_props%sulfate_index(:,14)=Aerosol_props%sulfate_index(:,13) + 26 ! 86%
      Aerosol_props%sulfate_index(:,15)=Aerosol_props%sulfate_index(:,14)
      Aerosol_props%sulfate_index(:,16)=Aerosol_props%sulfate_index(:,15) + 26 ! 84%
      Aerosol_props%sulfate_index(:,17)=Aerosol_props%sulfate_index(:,16)
      Aerosol_props%sulfate_index(:,18)=Aerosol_props%sulfate_index(:,17) + 26 ! 82%
      Aerosol_props%sulfate_index(:,19)=Aerosol_props%sulfate_index(:,18)
      Aerosol_props%sulfate_index(:,20)=Aerosol_props%sulfate_index(:,19) + 26 ! 80%
      Aerosol_props%sulfate_index(:,21)=Aerosol_props%sulfate_index(:,20)
      Aerosol_props%sulfate_index(:,22)=Aerosol_props%sulfate_index(:,21)
      Aerosol_props%sulfate_index(:,23)=Aerosol_props%sulfate_index(:,22) + 26
      Aerosol_props%sulfate_index(:,24)=Aerosol_props%sulfate_index(:,23)
      Aerosol_props%sulfate_index(:,25)=Aerosol_props%sulfate_index(:,24)      ! 75%
      Aerosol_props%sulfate_index(:,26)=Aerosol_props%sulfate_index(:,25)
      Aerosol_props%sulfate_index(:,27)=Aerosol_props%sulfate_index(:,26)
      Aerosol_props%sulfate_index(:,28)=Aerosol_props%sulfate_index(:,27) + 26
      Aerosol_props%sulfate_index(:,29)=Aerosol_props%sulfate_index(:,28)
      Aerosol_props%sulfate_index(:,30)=Aerosol_props%sulfate_index(:,29)      ! 70%
      Aerosol_props%sulfate_index(:,31)=Aerosol_props%sulfate_index(:,30)
      Aerosol_props%sulfate_index(:,32)=Aerosol_props%sulfate_index(:,31)
      Aerosol_props%sulfate_index(:,33)=Aerosol_props%sulfate_index(:,32) + 26
      Aerosol_props%sulfate_index(:,34)=Aerosol_props%sulfate_index(:,33)
      Aerosol_props%sulfate_index(:,35)=Aerosol_props%sulfate_index(:,34)      ! 65%
      Aerosol_props%sulfate_index(:,36)=Aerosol_props%sulfate_index(:,35)
      Aerosol_props%sulfate_index(:,37)=Aerosol_props%sulfate_index(:,36)
      Aerosol_props%sulfate_index(:,38)=Aerosol_props%sulfate_index(:,37) + 26
      Aerosol_props%sulfate_index(:,39)=Aerosol_props%sulfate_index(:,38)
      Aerosol_props%sulfate_index(:,40)=Aerosol_props%sulfate_index(:,39)      ! 60%
      Aerosol_props%sulfate_index(:,41)=Aerosol_props%sulfate_index(:,40)
      Aerosol_props%sulfate_index(:,42)=Aerosol_props%sulfate_index(:,41)
      Aerosol_props%sulfate_index(:,43)=Aerosol_props%sulfate_index(:,42) + 26
      Aerosol_props%sulfate_index(:,44)=Aerosol_props%sulfate_index(:,43)
      Aerosol_props%sulfate_index(:,45)=Aerosol_props%sulfate_index(:,44)      ! 55%
      Aerosol_props%sulfate_index(:,46)=Aerosol_props%sulfate_index(:,45)
      Aerosol_props%sulfate_index(:,47)=Aerosol_props%sulfate_index(:,46)
      Aerosol_props%sulfate_index(:,48)=Aerosol_props%sulfate_index(:,47) + 26
      Aerosol_props%sulfate_index(:,49)=Aerosol_props%sulfate_index(:,48)
      Aerosol_props%sulfate_index(:,50)=Aerosol_props%sulfate_index(:,49)      ! 50%
      Aerosol_props%sulfate_index(:,51)=Aerosol_props%sulfate_index(:,50)
      Aerosol_props%sulfate_index(:,52)=Aerosol_props%sulfate_index(:,51)
      Aerosol_props%sulfate_index(:,53)=Aerosol_props%sulfate_index(:,52) + 26
      Aerosol_props%sulfate_index(:,54)=Aerosol_props%sulfate_index(:,53)
      Aerosol_props%sulfate_index(:,55)=Aerosol_props%sulfate_index(:,54)      ! 45%
      Aerosol_props%sulfate_index(:,56)=Aerosol_props%sulfate_index(:,55)
      Aerosol_props%sulfate_index(:,57)=Aerosol_props%sulfate_index(:,56)
      Aerosol_props%sulfate_index(:,58)=Aerosol_props%sulfate_index(:,57) + 26
      Aerosol_props%sulfate_index(:,59)=Aerosol_props%sulfate_index(:,58)
      Aerosol_props%sulfate_index(:,60)=Aerosol_props%sulfate_index(:,59)      ! 40%
      Aerosol_props%sulfate_index(:,61)=Aerosol_props%sulfate_index(:,60)
      Aerosol_props%sulfate_index(:,62)=Aerosol_props%sulfate_index(:,61)
      Aerosol_props%sulfate_index(:,63)=Aerosol_props%sulfate_index(:,62) + 26
      Aerosol_props%sulfate_index(:,64)=Aerosol_props%sulfate_index(:,63)
      Aerosol_props%sulfate_index(:,65)=Aerosol_props%sulfate_index(:,64)      ! 35%
      Aerosol_props%sulfate_index(:,66)=Aerosol_props%sulfate_index(:,65)
      Aerosol_props%sulfate_index(:,67)=Aerosol_props%sulfate_index(:,66)
      Aerosol_props%sulfate_index(:,68)=Aerosol_props%sulfate_index(:,67) + 26
      Aerosol_props%sulfate_index(:,69)=Aerosol_props%sulfate_index(:,68)
      Aerosol_props%sulfate_index(:,70)=Aerosol_props%sulfate_index(:,69)      ! 30%
      Aerosol_props%sulfate_index(:,71)=Aerosol_props%sulfate_index(:,70)
      Aerosol_props%sulfate_index(:,72)=Aerosol_props%sulfate_index(:,71)
      Aerosol_props%sulfate_index(:,73)=Aerosol_props%sulfate_index(:,72) + 26
      Aerosol_props%sulfate_index(:,74)=Aerosol_props%sulfate_index(:,73)
      Aerosol_props%sulfate_index(:,75)=Aerosol_props%sulfate_index(:,74)      ! 25%
      Aerosol_props%sulfate_index(:,76)=Aerosol_props%sulfate_index(:,75)
      Aerosol_props%sulfate_index(:,77)=Aerosol_props%sulfate_index(:,76)
      Aerosol_props%sulfate_index(:,78)=Aerosol_props%sulfate_index(:,77) + 26
      Aerosol_props%sulfate_index(:,79)=Aerosol_props%sulfate_index(:,78)
      Aerosol_props%sulfate_index(:,80)=Aerosol_props%sulfate_index(:,79)      ! 20%
      Aerosol_props%sulfate_index(:,81)=Aerosol_props%sulfate_index(:,80)
      Aerosol_props%sulfate_index(:,82)=Aerosol_props%sulfate_index(:,81)
      Aerosol_props%sulfate_index(:,83)=Aerosol_props%sulfate_index(:,82) + 26
      Aerosol_props%sulfate_index(:,84)=Aerosol_props%sulfate_index(:,83)
      Aerosol_props%sulfate_index(:,85)=Aerosol_props%sulfate_index(:,84)      ! 15%
      Aerosol_props%sulfate_index(:,86)=Aerosol_props%sulfate_index(:,85)
      Aerosol_props%sulfate_index(:,87)=Aerosol_props%sulfate_index(:,86)
      Aerosol_props%sulfate_index(:,88)=Aerosol_props%sulfate_index(:,87) + 26
      Aerosol_props%sulfate_index(:,89)=Aerosol_props%sulfate_index(:,88)
      Aerosol_props%sulfate_index(:,90)=Aerosol_props%sulfate_index(:,89)      ! 10%
      Aerosol_props%sulfate_index(:,91)=Aerosol_props%sulfate_index(:,90)
      Aerosol_props%sulfate_index(:,92)=Aerosol_props%sulfate_index(:,91)
      Aerosol_props%sulfate_index(:,93)=Aerosol_props%sulfate_index(:,92) + 26
      Aerosol_props%sulfate_index(:,94)=Aerosol_props%sulfate_index(:,93)
      Aerosol_props%sulfate_index(:,95)=Aerosol_props%sulfate_index(:,94)      ! 5%
      Aerosol_props%sulfate_index(:,96)=Aerosol_props%sulfate_index(:,95)
      Aerosol_props%sulfate_index(:,97)=Aerosol_props%sulfate_index(:,96)
      Aerosol_props%sulfate_index(:,98)=Aerosol_props%sulfate_index(:,97) + 26
      Aerosol_props%sulfate_index(:,99)=Aerosol_props%sulfate_index(:,98)
      Aerosol_props%sulfate_index(:,100)=Aerosol_props%sulfate_index(:,99)     ! 0%
    end if
    Aerosol_props%bcphilic_index( 0: 13)=(/ 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722 /)
    Aerosol_props%bcphilic_index(14: 27)=(/ 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722, 722 /)
    Aerosol_props%bcphilic_index(28: 41)=(/ 722, 722, 722, 722, 722, 723, 723, 723, 723, 723, 724, 724, 724, 724 /)
    Aerosol_props%bcphilic_index(42: 55)=(/ 724, 725, 725, 725, 725, 725, 726, 726, 726, 726, 726, 727, 727, 727 /)
    Aerosol_props%bcphilic_index(56: 69)=(/ 727, 727, 728, 728, 728, 728, 728, 729, 729, 729, 729, 729, 730, 730 /)
    Aerosol_props%bcphilic_index(70: 83)=(/ 730, 730, 730, 731, 731, 731, 731, 731, 732, 732, 732, 732, 733, 733 /)
    Aerosol_props%bcphilic_index(84: 97)=(/ 734, 734, 735, 735, 736, 736, 737, 738, 739, 740, 741, 742, 743, 744 /)
    Aerosol_props%bcphilic_index(98:100)=(/ 744, 744, 744 /)
    Aerosol_props%omphilic_index(:) = Aerosol_props%bcphilic_index(:) + 25
    Aerosol_props%seasalt1_index(:) = Aerosol_props%omphilic_index(:) + 25
    Aerosol_props%seasalt2_index(:) = Aerosol_props%seasalt1_index(:) + 25
    Aerosol_props%seasalt3_index(:) = Aerosol_props%seasalt2_index(:) + 25
    Aerosol_props%seasalt4_index(:) = Aerosol_props%seasalt3_index(:) + 25
    Aerosol_props%seasalt5_index(:) = Aerosol_props%seasalt4_index(:) + 25
    call loadaerooptical(Aerosol_props)
    
    Aerosol_diags%extopdep=0.
    Aerosol_diags%absopdep=0.
    Aerosol_diags%asymdep=0.

    if (include_volcanoes) then
      write(6,*) "ERROR: Prescribed aerosol properties for"
      write(6,*) "volcanoes is currently unsupported"
      call ccmpi_abort(-1)
      !allocate(Aerosol_props%sw_ext(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%sw_ssa(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%sw_asy(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%lw_ext(imax,1,kl,N_AEROSOL_BANDS))      
      !allocate(Aerosol_diags%lw_extopdep_vlcno(imax,1,kl+1,2))
      !allocate(Aerosol_diags%lw_absopdep_vlcno(imax,1,kl+1,2))
    end if

  end if

  ! assign GHG concentrations
  Rad_gases%rrvco2  = real(rrvco2 ,8)
  Rad_gases%rrvch4  = real(rrvch4 ,8)
  Rad_gases%rrvn2o  = real(rrvn2o ,8)
  Rad_gases%rrvf11  = real(rrvf11 ,8)
  Rad_gases%rrvf12  = real(rrvf12 ,8)
  Rad_gases%rrvf113 = real(rrvf113,8)
  Rad_gases%rrvf22  = real(rrvf22 ,8)
  call sealw99_time_vary(Rad_time, Rad_gases)
  
  ! define diagnostic cloud levels
  f1 = 1.
  f2 = 1.
  do k = 1,kl-1
    if ( abs(sigmh(k+1)-siglow)<f1 ) then
      f1   = abs(sigmh(k+1)-siglow)
      nlow = k
    end if
    if ( abs(sigmh(k+1)-sigmid)<f2 ) then
      f2   = abs(sigmh(k+1)-sigmid)
      nmid = k
    end if
  end do

  ! initialise VIS fraction of SW radiation
  swrsave = 0.5
  
end if  ! (first)

if ( nmaxpr==1 .and. myid==0 ) then
  write(6,*) "seaesfrad: Prepare SEA-ESF arrays"
end if

! error checking
if ( ldr==0 ) then
  write(6,*) "ERROR: SEA-ESF radiation requires ldr/=0"
  call ccmpi_abort(-1)
end if

if ( mod(ifull,imax)/=0 ) then
  ! imax should be automatically set-up in globpe.f
  ! so an error here should indicate a bug in globpe.f
  write(6,*) 'nproc,il,jl,ifull,imax ',nproc,il,jl,ifull,imax
  write(6,*) 'illegal setting of imax in rdparm'
  call ccmpi_abort(-1)
endif

! main loop ---------------------------------------------------------
do j = 1,jl,imax/il
  istart = 1+(j-1)*il
  iend   = istart+imax-1
  
  if ( nmaxpr==1 .and. myid==0 ) then
    write(6,*) "seaesfrad: Main SEA-ESF loop for istart,iend ",istart,iend
  end if

  ! Calculate zenith angle for the solarfit calculation.
  ! This call averages zenith angle just over this time step.
  dhr = dt/3600.
  call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro2,taudar2)
  call atebccangle(istart,imax,coszro2(1:imax),rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt,sin(dlt))

  ! Call radiation --------------------------------------------------
  if ( odcalc ) then     ! Do the calculation
  
    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Update radiation"
    end if


    ! Average the zenith angle over the time (hours) between radiation
    ! calculations
    dhr = kountr*dt/3600.0
    call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro,taudar)
    
    ! Set up ozone for this time and row
    if (amipo3) then
      call o3set_amip( rlatt(istart:iend), imax, mins,sigh, ps(istart:iend), duo3n )
    else
      ! note levels are inverted
      call o3set( imax, istart, mins, duo3n, sig, ps(istart:iend) )
    end if
    Rad_gases%qo3(:,1,:)=max(1.e-10_8,real(duo3n,8))

    ! Set-up albedo
    ! Land albedo ---------------------------------------------------
    if (nsib==6.or.nsib==7) then
      ! CABLE version
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albvisdir(istart:iend) ! from cable (inc snow)
        cirrf_dir(1:imax) = albnirdir(istart:iend) ! from cable (inc snow)
        cuvrf_dif(1:imax) = albvisdif(istart:iend) ! from cable (inc snow)
        cirrf_dif(1:imax) = albnirdif(istart:iend) ! from cable (inc snow)
      end where
    else
      ! nsib=3 version (calculate snow)
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albvissav(istart:iend) ! from albfile (indata.f)
        cirrf_dir(1:imax) = albnirsav(istart:iend) ! from albnirfile (indata.f)
        cuvrf_dif(1:imax) = cuvrf_dir(1:imax)      ! assume DIR and DIF are the same
        cirrf_dif(1:imax) = cirrf_dir(1:imax)      ! assume DIR and DIF are the same
      end where
      ! The following snow calculation should be done by sib3 (sflux.f)
      do i=1,imax
        iq=i+(j-1)*il
        if (land(iq)) then
          if (snowd(iq)>0.) then
            dnsnow=min(1.,.1*max(0.,snowd(iq)-osnowd(iq)))
            ttbg=real(isflag(iq))*tggsn(iq,1) + real(1-isflag(iq))*tgg(iq,1)
            ttbg=min(ttbg,273.1)
            ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
            exp_ar1=exp(ar1)                  ! e.g. exp(0 to -4)
            ar2 = 10.*ar1                     ! freezing of melt water
            exp_ar2=exp(ar2)                  ! e.g. exp(0 to -40)
            snr=snowd(iq)/max(ssdnn(iq),100.)
            if(isoilm(iq)==9)then   ! fixes for Arctic & Antarctic
              ar3=.001
              dnsnow=max(dnsnow,.0015)
              snrat=min(1.,snr/(snr+.001))
            else
              ar3=.3               ! accumulation of dirt
              snrat=min(1.,snr/(snr+.02))
            endif
            dtau=1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
            if(snowd(iq)<= 1.)then
              snage(iq)=0.
            else
              snage(iq)=max(0.,(snage(iq) + dtau)*(1.-dnsnow))
            endif
            alvo = 0.95         !alb. for vis. on a new snow
            aliro = 0.65        !alb. for near-infr. on a new snow
            fage = 1.-1./(1.+snage(iq))  !age factor
            cczen=max(.17365, coszro(i))
            fzen=( 1.+1./2.)/(1.+2.*2.*cczen) -1./2.
            if( cczen > 0.5 ) fzen = 0.
            fzenm = max ( fzen, 0. )
            alvd = alvo * (1.0-0.2*fage)
            alv = .4 * fzenm * (1.-alvd) + alvd
            alird = aliro*(1.-.5*fage)
            alir = .4 * fzenm * (1.0-alird) + alird
            cuvrf_dir(i)=(1.-snrat)*cuvrf_dir(i) + snrat*alv
            cirrf_dir(i)=(1.-snrat)*cirrf_dir(i) + snrat*alir
            cuvrf_dif(i)=cuvrf_dir(i) ! assume DIR and DIF are the same
            cirrf_dif(i)=cirrf_dir(i) ! assume DIR and DIF are the same
          end if
        end if
      end do
    end if

    ! Water/Ice albedo --------------------------------------------
    if (nmlo==0) then
      ! NCAR CCMS3.0 scheme (Briegleb et al, 1986,
      ! J. Clim. and Appl. Met., v. 27, 214-226)
      where (.not.land(istart:iend).and.coszro(1:imax)>=0.)
        cuvrf_dir(1:imax)=0.026/(coszro(1:imax)**1.7+0.065)                  &
          +0.15*(coszro(1:imax)-0.1)*(coszro(1:imax)-0.5)*(coszro(1:imax)-1.)
      elsewhere (.not.land(istart:iend))
        cuvrf_dir(1:imax)=0.3925 ! coszen=0 value of above expression
      end where
      where (.not.land(istart:iend))
        cuvrf_dif(1:imax)=0.06
        cirrf_dir(1:imax)=cuvrf_dir(1:imax)
        cirrf_dif(1:imax)=0.06
        cuvrf_dir(1:imax)=0.85*fracice(istart:iend)+(1.-fracice(istart:iend))*cuvrf_dir(1:imax)
        cuvrf_dif(1:imax)=0.85*fracice(istart:iend)+(1.-fracice(istart:iend))*cuvrf_dif(1:imax)
        cirrf_dir(1:imax)=0.45*fracice(istart:iend)+(1.-fracice(istart:iend))*cirrf_dir(1:imax)
        cirrf_dif(1:imax)=0.45*fracice(istart:iend)+(1.-fracice(istart:iend))*cirrf_dif(1:imax)
      end where
    elseif (abs(nmlo)<=9) then
      ! MLO albedo ----------------------------------------------------
      call mloalb4(istart,imax,coszro,cuvrf_dir,cuvrf_dif,cirrf_dir,cirrf_dif,0)
    else
      ! PCOM
      write(6,*) "ERROR: This PCOM option for SEA-ESF radiation is not currently supported"
      call ccmpi_abort(-1)
    end if

    ! Urban albedo --------------------------------------------------
    call atebalb1(istart,imax,cuvrf_dir(1:imax),0,split=1)
    call atebalb1(istart,imax,cirrf_dir(1:imax),0,split=1)
    call atebalb1(istart,imax,cuvrf_dif(1:imax),0,split=2)
    call atebalb1(istart,imax,cirrf_dif(1:imax),0,split=2)

    ! Aerosols -------------------------------------------------------
    tnhs=phi_nh(istart:iend,1)/bet(1)
    tv=t(istart:iend,1)*(1.+0.61*qg(istart:iend,1)-qlrad(istart:iend,1)-qfrad(istart:iend,1))
    rhoa(:,1)=ps(istart:iend)*sig(1)/(rdry*tv) !density of air
    dz(:,1)=-rdry*dsig(1)*(tv+tnhs)/(grav*sig(1))
    do k=2,kl
      ! representing non-hydrostatic term as a correction to air temperature
      tnhs=(phi_nh(istart:iend,k)-phi_nh(istart:iend,k-1)-betm(k)*tnhs)/bet(k)
      tv=t(istart:iend,k)*(1.+0.61*qg(istart:iend,k)-qlrad(istart:iend,k)-qfrad(istart:iend,k))
      rhoa(:,k)=ps(istart:iend)*sig(k)/(rdry*tv) !density of air
      dz(:,k)=-rdry*dsig(k)*(tv+tnhs)/(grav*sig(k))
    end do
    select case (abs(iaero))
      case(0)
        ! no aerosols
      case(1)
        ! aerosols are read in (direct effect only)
        do i=1,imax
          iq=i+(j-1)*il
          cosz = max ( coszro(i), 1.e-4)
          delta = coszro(i)*0.29*8.*so4t(iq)*((1.-0.25*(cuvrf_dir(i)+cuvrf_dif(i)+cirrf_dir(i)+cirrf_dif(i)))/cosz)**2
          cuvrf_dir(i)=min(0.99, delta+cuvrf_dir(i)) ! still broadband
          cirrf_dir(i)=min(0.99, delta+cirrf_dir(i)) ! still broadband
          cuvrf_dif(i)=min(0.99, delta+cuvrf_dif(i)) ! still broadband
          cirrf_dif(i)=min(0.99, delta+cirrf_dif(i)) ! still broadband
        end do ! i=1,imax
      case(2)
        ! prognostic aerosols
        ! convert to units kg / m^2
        do k=1,kl
          kr=kl+1-k
          dzrho=rhoa(:,k)*dz(:,k)
          ! Factor of 132.14/32.06 converts from sulfur to ammmonium sulfate
          Aerosol%aerosol(:,1,kr,1) =real((132.14/32.06)*xtg(istart:iend,k,3)*dzrho,8) ! so4
          Aerosol%aerosol(:,1,kr,2) =real(xtg(istart:iend,k,4)*dzrho,8)                ! bc hydrophobic
          Aerosol%aerosol(:,1,kr,3) =real(xtg(istart:iend,k,5)*dzrho,8)                ! bc hydrophilic
          Aerosol%aerosol(:,1,kr,4) =real(xtg(istart:iend,k,6)*dzrho,8)                ! oc hydrophobic
          Aerosol%aerosol(:,1,kr,5) =real(xtg(istart:iend,k,7)*dzrho,8)                ! oc hydrophilic
          Aerosol%aerosol(:,1,kr,6) =real(xtg(istart:iend,k,8)*dzrho,8)                ! dust 0.8
          Aerosol%aerosol(:,1,kr,7) =real(xtg(istart:iend,k,9)*dzrho,8)                ! dust 1.0
          Aerosol%aerosol(:,1,kr,8) =real(xtg(istart:iend,k,10)*dzrho,8)               ! dust 2.0
          Aerosol%aerosol(:,1,kr,9) =real(xtg(istart:iend,k,11)*dzrho,8)               ! dust 4.0
          !Aerosol%aerosol(:,1,kr,10)=real((2.64e-18*ssn(istart:iend,k,1)  & ! Small film sea salt (0.035)
          !                                +1.38e-15*ssn(istart:iend,k,2)) & ! Large jet sea salt (0.35)
          !                           /rhoa(:,k)*dzrho,8)   
          Aerosol%aerosol(:,1,kr,10)=real((5.3e-17*ssn(istart:iend,k,1)  & ! Small film sea salt (0.1)
                                          +9.1e-15*ssn(istart:iend,k,2)) & ! Large jet sea salt (0.5)
                                         *dzrho/rhoa(:,k),8)                
        end do
        Aerosol%aerosol=max(Aerosol%aerosol,0._8)
        
        if ( Rad_control%using_im_bcsul ) then
          Aerosol_props%ivol(:,1,:)=100-nint(100.*Aerosol%aerosol(:,1,:,1)/ &
              max(Aerosol%aerosol(:,1,:,1)+Aerosol%aerosol(:,1,:,2)+Aerosol%aerosol(:,1,:,3),1.E-30_8))
          Aerosol_props%ivol(:,1,:)=max(min(Aerosol_props%ivol(:,1,:),100),0)
        else
          Aerosol_props%ivol=0 ! no mixing of bc with so4
        end if

      case DEFAULT
        write(6,*) "ERROR: unknown iaero option ",iaero
        call ccmpi_abort(-1)
    end select

    ! define droplet size distribution ------------------------------
    call aerodrop(istart,imax,kl,cd2,rhoa,land(istart:iend),rlatt(istart:iend))
    
    ! Cloud fraction diagnostics ------------------------------------
    cloudlo(istart:iend)=0.
    cloudmi(istart:iend)=0.
    cloudhi(istart:iend)=0.
    ! Diagnose low, middle and high clouds
    if (nmr>0) then
      ! max-rnd cloud overlap
      mx=0.
      do k=1,nlow
        mx=max(mx,cfrac(istart:iend,k))
        where (cfrac(istart:iend,k)==0.)
          cloudlo(istart:iend)=cloudlo(istart:iend)+mx*(1.-cloudlo(istart:iend))
          mx=0.
        end where
      end do
      cloudlo(istart:iend)=cloudlo(istart:iend)+mx*(1.-cloudlo(istart:iend))
      mx=0.
      do k=nlow+1,nmid
        mx=max(mx,cfrac(istart:iend,k))
        where (cfrac(istart:iend,k)==0.)
          cloudmi(istart:iend)=cloudmi(istart:iend)+mx*(1.-cloudmi(istart:iend))
          mx=0.
        end where
      end do
      cloudmi(istart:iend)=cloudmi(istart:iend)+mx*(1.-cloudmi(istart:iend))
      mx=0.
      do k=nmid+1,kl-1
        mx=max(mx,cfrac(istart:iend,k))
        where (cfrac(istart:iend,k)==0.)
          cloudhi(istart:iend)=cloudhi(istart:iend)+mx*(1.-cloudhi(istart:iend))
          mx=0.
        end where
      end do
      cloudhi(istart:iend)=cloudhi(istart:iend)+mx*(1.-cloudhi(istart:iend))  
    else
      ! rnd cloud overlap
      do k=1,nlow
        cloudlo(istart:iend)=cloudlo(istart:iend)+cfrac(istart:iend,k)*(1.-cloudlo(istart:iend))
      enddo
      do k=nlow+1,nmid
        cloudmi(istart:iend)=cloudmi(istart:iend)+cfrac(istart:iend,k)*(1.-cloudmi(istart:iend))
      enddo
      do k=nmid+1,kl-1
        cloudhi(istart:iend)=cloudhi(istart:iend)+cfrac(istart:iend,k)*(1.-cloudhi(istart:iend))
      enddo
    end if

    ! Prepare SEA-ESF arrays ----------------------------------------
    do k=1,kl
      kr=kl+1-k
      dumt(:,k)=t(istart:iend,k)
      tv=dumt(:,k)*(1.+0.61*qg(istart:iend,k)-qlrad(istart:iend,k)-qfrad(istart:iend,k))
      p2(:,k)=ps(istart:iend)*sig(k)
      call getqsat(imax,qsat,dumt(:,k),p2(:,k))
      Atmos_input%deltaz(:,1,kr)  = real(dz(:,k),8)
      Atmos_input%rh2o(:,1,kr)    = max(real(qg(istart:iend,k),8),2.E-7_8)
      Atmos_input%temp(:,1,kr)    = min(max(real(dumt(:,k),8),100._8),370._8)    
      Atmos_input%press(:,1,kr)   = real(p2(:,k),8)
      Atmos_input%rel_hum(:,1,kr) = min(real(qg(istart:iend,k)/qsat,8),1._8)
    end do
    Atmos_input%temp(:,1,kl+1)  = min(max(real(tss(istart:iend),8),100._8),370._8)
    Atmos_input%press(:,1,kl+1) = real(ps(istart:iend),8)
    Atmos_input%pflux(:,1,1  )  = 0._8
    Atmos_input%tflux(:,1,1  )  = Atmos_input%temp(:,1,1)
    do k=1,kl-1
      kr=kl+1-k
      Atmos_input%pflux(:,1,kr) = real(rathb(k)*p2(:,k)+ratha(k)*p2(:,k+1),8)
      Atmos_input%tflux(:,1,kr) = min(max(real(rathb(k)*dumt(:,k)+ratha(k)*dumt(:,k+1),8),100._8),370._8)
    end do
    Atmos_input%pflux(:,1,kl+1) = real(ps(istart:iend),8)
    Atmos_input%tflux(:,1,kl+1) = min(max(real(tss(istart:iend),8),100._8),370._8)
    Atmos_input%clouddeltaz     = Atmos_input%deltaz

    Atmos_input%psfc(:,1)    = real(ps(istart:iend),8)
    Atmos_input%tsfc(:,1)    = min(max(real(tss(istart:iend),8),100._8),370._8)
    Atmos_input%phalf(:,1,1) = 0._8
    do k=1,kl-1
      kr=kl+1-k
      Atmos_input%phalf(:,1,kr) = real(rathb(k)*p2(:,k)+ratha(k)*p2(:,k+1),8)
    end do
    Atmos_input%phalf(:,1,kl+1) = real(ps(istart:iend),8)
    if (do_aerosol_forcing) then
      Atmos_input%aerosolrelhum = Atmos_input%rel_hum
    end if
    
    ! cloud overlap
    if (nmr>0) then
      do i=1,imax ! maximum-random overlap
        iq=i+istart-1
        Cld_spec%cmxolw(i,1,:)=0._8
        k=1
        do while (k<kl)
          ktop=k
          if (cfrac(iq,k)>0.) then
            kbot=k ! found bottom of cloud
            do while (ktop<kl.and.cfrac(iq,ktop+1)>0.)
              ktop=ktop+1 ! search for top of cloud
            end do
            if (ktop>kbot) then ! if multi-layer cloud, calculate common max overlap fraction
              Cld_spec%cmxolw(i,1,kl+1-ktop:kl+1-kbot)=real(minval(cfrac(iq,kbot:ktop)),8)
            end if
          end if
          k=ktop+1
        end do
        do k=1,kl
          kr=kl+1-k
          Cld_spec%camtsw(i,1,kr)=real(cfrac(iq,k),8)                         ! Max+Rnd overlap clouds for SW
          Cld_spec%cmxolw(i,1,kr)=min(Cld_spec%cmxolw(i,1,kr),0.999_8)        ! Max overlap for LW
          Cld_spec%crndlw(i,1,kr)=real(cfrac(iq,k),8)-Cld_spec%cmxolw(i,1,kr) ! Rnd overlap for LW
        end do
      end do
    else
      do i=1,imax ! random overlap
        iq=i+istart-1
        do k=1,kl
          kr=kl+1-k
          Cld_spec%camtsw(i,1,kr)=real(cfrac(iq,k),8) ! Max+Rnd overlap clouds for SW
          Cld_spec%crndlw(i,1,kr)=real(cfrac(iq,k),8) ! Rnd overlap for LW
          Cld_spec%cmxolw(i,1,kr)=0._8
        end do
      end do
    end if

    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Calculate microphysics properties"
    end if

    ! cloud microphysics for radiation
    ! cfrac, qlrad and qfrad also include convective cloud as well as qfg and qlg
    dumcf=cfrac(istart:iend,:)
    dumql=qlrad(istart:iend,:)
    dumqf=qfrad(istart:iend,:)
    call cloud3(Cloud_microphysics%size_drop,Cloud_microphysics%size_ice,       &
                Cloud_microphysics%conc_drop,Cloud_microphysics%conc_ice,       &
                dumcf,dumql,dumqf,p2,dumt,cd2,imax,kl)
    Cloud_microphysics%size_drop=max(Cloud_microphysics%size_drop,1.e-20_8)
    Cloud_microphysics%size_ice =max(Cloud_microphysics%size_ice, 1.e-20_8)                
    Cloud_microphysics%size_rain=1.e-20_8
    Cloud_microphysics%conc_rain=0._8
    Cloud_microphysics%size_snow=1.e-20_8
    Cloud_microphysics%conc_snow=0._8
    
    Lscrad_props%cldext   = 0._8
    Lscrad_props%cldsct   = 0._8
    Lscrad_props%cldasymm = 1._8
    Lscrad_props%abscoeff = 0._8
    call microphys_lw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    call microphys_sw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    Cldrad_props%cldsct(:,:,:,:,1)  =Lscrad_props%cldsct(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldext(:,:,:,:,1)  =Lscrad_props%cldext(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldasymm(:,:,:,:,1)=Lscrad_props%cldasymm(:,:,:,:) ! Large scale cloud properties only
    Cldrad_props%abscoeff(:,:,:,:,1)=Lscrad_props%abscoeff(:,:,:,:) ! Large scale cloud properties only
    
    call lwemiss_calc(Atmos_input%clouddeltaz,Cldrad_props%abscoeff,Cldrad_props%cldemiss)
    Cldrad_props%emmxolw = Cldrad_props%cldemiss
    Cldrad_props%emrndlw = Cldrad_props%cldemiss

    Surface%asfc_vis_dir(:,1)=real(cuvrf_dir(:),8)
    Surface%asfc_nir_dir(:,1)=real(cirrf_dir(:),8)
    Surface%asfc_vis_dif(:,1)=real(cuvrf_dif(:),8)
    Surface%asfc_nir_dif(:,1)=real(cirrf_dif(:),8)
   
    Astro%cosz(:,1)   =max(real(coszro,8),0._8)
    Astro%fracday(:,1)=real(taudar,8)

    call END_LOG(radmisc_end)

    call START_LOG(radlw_begin)
    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Longwave radiation"
    end if
    call longwave_driver (1, imax, 1, 1, Rad_time, Atmos_input,  &
                          Rad_gases, Aerosol, Aerosol_props,     &
                          Cldrad_props, Cld_spec, Aerosol_diags, &
                          Lw_output)
    call END_LOG(radlw_end)

    call START_LOG(radsw_begin)
    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Shortwave radiation"
    end if
    call shortwave_driver (1, imax, 1, 1, Atmos_input, Surface,      &
                           Astro, Aerosol, Aerosol_props, Rad_gases, &
                           Cldrad_props, Cld_spec, Sw_output,        &
                           Aerosol_diags, r)
    call END_LOG(radsw_end)
    
    call START_LOG(radmisc_begin)

    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Process SEA-ESF output"
    end if

    ! store shortwave and fbeam data --------------------------------
    sgdn    = real(Sw_output(1)%dfsw(:,1,kl+1,1))
    sgdnvis = real(Sw_output(1)%dfsw_vis_sfc(:,1,1))
    sgdnnir = sgdn-sgdnvis
    sg      = sgdn-real(Sw_output(1)%ufsw(:,1,kl+1,1))
    sgvis   = sgdnvis-real(Sw_output(1)%ufsw_vis_sfc(:,1,1))
    !sgvisdir = Sw_output(1)%dfsw_vis_sfc_dir(:,1,1)
    !sgvisdif = Sw_output(1)%dfsw_vis_sfc_dif(:,1,1)-Sw_output(1)%ufsw_vis_sfc_dif(:,1,1)
    !sgnirdir = Sw_output(1)%dfsw_dir_sfc(:,1,1)-sgvisdir
    !sgnirdif = Sw_output(1)%dfsw_dif_sfc(:,1,1)-Sw_output(1)%ufsw_dif_sfc(:,1,1)-sgvisdif
    !sgdir    = Sw_output(1)%dfsw_dir_sfc(:,1,1)
    !sgdif    = Sw_output(1)%dfsw_dif_sfc(:,1,1)-Sw_output(1)%ufsw_dif_sfc(:,1,1)
    sgdnvisdir = real(Sw_output(1)%dfsw_vis_sfc_dir(:,1,1))
    sgdnvisdif = real(Sw_output(1)%dfsw_vis_sfc_dif(:,1,1))
    sgdnnirdir = real(Sw_output(1)%dfsw_dir_sfc(:,1,1))-sgdnvisdir
    sgdnnirdif = real(Sw_output(1)%dfsw_dif_sfc(:,1,1))-sgdnvisdif
    
    swrsave(istart:iend)=sgdnvis/max(sgdn,0.01)
    fbeamvis(istart:iend)=sgdnvisdir/max(sgdnvis,0.01)
    fbeamnir(istart:iend)=sgdnnirdir/max(sgdnnir,0.01)
    
    ! Store albedo data ---------------------------------------------
    albvisnir(istart:iend,1) = real(Surface%asfc_vis_dir(:,1))*fbeamvis(istart:iend) &
                              +real(Surface%asfc_vis_dif(:,1))*(1.-fbeamvis(istart:iend))
    albvisnir(istart:iend,2) = real(Surface%asfc_nir_dir(:,1))*fbeamnir(istart:iend) &
                              +real(Surface%asfc_nir_dif(:,1))*(1.-fbeamnir(istart:iend))
    
    ! longwave output -----------------------------------------------
    rg(1:imax) = real(Lw_output(1)%flxnet(:,1,kl+1))          ! longwave at surface
    rt(1:imax) = real(Lw_output(1)%flxnet(:,1,1))             ! longwave at top
    ! rg is net upwards = stefbo T^4 - Rdown
    rgdn(1:imax) = stefbo*tss(istart:iend)**4 - rg(1:imax)

    ! shortwave output ----------------------------------------------
    sint(1:imax) = real(Sw_output(1)%dfsw(:,1,1,1))   ! solar in top
    sout(1:imax) = real(Sw_output(1)%ufsw(:,1,1,1))   ! solar out top
    !sgdn(1:imax) = sg(1:imax) / ( 1. - swrsave(istart:iend)*albvisnir(istart:iend,1) &
    !              -(1.-swrsave(istart:iend))*albvisnir(istart:iend,2) )

    ! Clear sky calculation -----------------------------------------
    if (do_totcld_forcing) then
      soutclr(1:imax) = real(Sw_output(1)%ufswcf(:,1,1,1))      ! solar out top
      sgclr(1:imax)   = -real(Sw_output(1)%fswcf(:,1,kl+1,1))   ! solar absorbed at the surface
      rtclr(1:imax)   = real(Lw_output(1)%flxnetcf(:,1,1))      ! clr sky lw at top
      rgclr(1:imax)   = real(Lw_output(1)%flxnetcf(:,1,kl+1))   ! clear sky longwave at surface
    else
      soutclr(1:imax) = 0.
      sgclr(1:imax)   = 0.
      rtclr(1:imax)   = 0.
      rgclr(1:imax)   = 0.
    end if

    ! heating rate --------------------------------------------------
    do k=1,kl
      ! total heating rate (convert deg K/day to deg K/sec)
      rtt(istart:iend,kl+1-k)=-real(Sw_output(1)%hsw(:,1,k,1)+Lw_output(1)%heatra(:,1,k))/86400.
    end do
    
    ! aerosol optical depths ----------------------------------------
    if ( do_aerosol_forcing ) then
      opticaldepth(istart:iend,:,:)=0.
      ! Sulfate
      do k=1,kl
        opticaldepth(istart:iend,3,1)=opticaldepth(istart:iend,3,1)+real(Aerosol_diags%extopdep(1:imax,1,k,1,1)) ! Visible
        opticaldepth(istart:iend,3,2)=opticaldepth(istart:iend,3,2)+real(Aerosol_diags%extopdep(1:imax,1,k,1,2)) ! Near IR
        opticaldepth(istart:iend,3,3)=opticaldepth(istart:iend,3,3)+real(Aerosol_diags%extopdep(1:imax,1,k,1,3)) ! Longwave
      end do
      ! BC
      do nr=2,3
        do k=1,kl          
          opticaldepth(istart:iend,5,1)=opticaldepth(istart:iend,5,1)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,5,2)=opticaldepth(istart:iend,5,2)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,5,3)=opticaldepth(istart:iend,5,3)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! OC
      do nr=4,5
        do k=1,kl    
          opticaldepth(istart:iend,6,1)=opticaldepth(istart:iend,6,1)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,6,2)=opticaldepth(istart:iend,6,2)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,6,3)=opticaldepth(istart:iend,6,3)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! Small dust
      do k=1,kl
        opticaldepth(istart:iend,1,1)=opticaldepth(istart:iend,1,1)+real(Aerosol_diags%extopdep(1:imax,1,k,6,1)) ! Visible
        opticaldepth(istart:iend,1,2)=opticaldepth(istart:iend,1,2)+real(Aerosol_diags%extopdep(1:imax,1,k,6,2)) ! Near IR
        opticaldepth(istart:iend,1,3)=opticaldepth(istart:iend,1,3)+real(Aerosol_diags%extopdep(1:imax,1,k,6,3)) ! Longwave
      end do
      ! Large dust
      do nr=7,9
        do k=1,kl
          opticaldepth(istart:iend,2,1)=opticaldepth(istart:iend,2,1)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,2,2)=opticaldepth(istart:iend,2,2)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,2,3)=opticaldepth(istart:iend,2,3)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! Seasalt
      do k=1,kl
        opticaldepth(istart:iend,7,1)=opticaldepth(istart:iend,7,1)+real(Aerosol_diags%extopdep(1:imax,1,k,10,1)) ! Visible
        opticaldepth(istart:iend,7,2)=opticaldepth(istart:iend,7,2)+real(Aerosol_diags%extopdep(1:imax,1,k,10,2)) ! Near IR
        opticaldepth(istart:iend,7,3)=opticaldepth(istart:iend,7,3)+real(Aerosol_diags%extopdep(1:imax,1,k,10,3)) ! Longwave
      end do
      ! Aerosol
      do nr=1,nfields
        do k=1,kl
          opticaldepth(istart:iend,4,1)=opticaldepth(istart:iend,4,1)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,4,2)=opticaldepth(istart:iend,4,2)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,4,3)=opticaldepth(istart:iend,4,3)+real(Aerosol_diags%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
    end if

    ! Calculate the amplitude of the diurnal cycle of solar radiation
    ! at the surface (using the value for the middle of the radiation
    ! step) and use this value to get solar radiation at other times.
    ! Use the zenith angle and daylight fraction calculated in zenith
    ! to remove these factors.
    where (coszro(1:imax)*taudar(1:imax)<=1.E-5)
      ! The sun isn't up at all over the radiation period so no 
      ! fitting need be done.
      sga(1:imax)=0.
    elsewhere
      sga(1:imax)=sg(1:imax)/(coszro(1:imax)*taudar(1:imax))
    end where

    ! Save things for non-radiation time steps ----------------------
    sgsave(istart:iend)   = sg(1:imax)   ! repeated after solarfit
    sgamp(istart:iend)    = sga(1:imax)
    ! Save the value excluding Ts^4 part.  This is allowed to change.
    rgsave(istart:iend)   = rg(1:imax)-stefbo*tss(istart:iend)**4
    sintsave(istart:iend) = sint(1:imax) 
    rtsave(istart:iend)   = rt(1:imax) 
    rtclsave(istart:iend) = rtclr(1:imax)  
    sgclsave(istart:iend) = sgclr(1:imax)

    ! cloud amounts for saving --------------------------------------
    cloudtot(istart:iend)=1.-(1.-cloudlo(istart:iend))*(1.-cloudmi(istart:iend))*(1.-cloudhi(istart:iend))

    if (nmaxpr==1.and.myid==0) then
      write(6,*) "seaesfrad: Calculate averages"
    end if

    ! Use explicit indexing rather than array notation so that we can run
    ! over the end of the first index
    if(ktau>0)then ! averages not added at time zero
      if(j==1)koundiag=koundiag+1  
      sint_ave(istart:iend) = sint_ave(istart:iend) + sint(1:imax)
      sot_ave(istart:iend)  = sot_ave(istart:iend)  + sout(1:imax)
      soc_ave(istart:iend)  = soc_ave(istart:iend)  + soutclr(1:imax)
      rtu_ave(istart:iend)  = rtu_ave(istart:iend)  + rt(1:imax)
      rtc_ave(istart:iend)  = rtc_ave(istart:iend)  + rtclr(1:imax)
      rgn_ave(istart:iend)  = rgn_ave(istart:iend)  + rg(1:imax)
      rgc_ave(istart:iend)  = rgc_ave(istart:iend)  + rgclr(1:imax)
      rgdn_ave(istart:iend) = rgdn_ave(istart:iend) + rgdn(1:imax)
      sgdn_ave(istart:iend) = sgdn_ave(istart:iend) + sgdn(1:imax)
      sgc_ave(istart:iend)  = sgc_ave(istart:iend)  + sgclr(1:imax)
      cld_ave(istart:iend)  = cld_ave(istart:iend)  + cloudtot(istart:iend)
      cll_ave(istart:iend)  = cll_ave(istart:iend)  + cloudlo(istart:iend)
      clm_ave(istart:iend)  = clm_ave(istart:iend)  + cloudmi(istart:iend)
      clh_ave(istart:iend)  = clh_ave(istart:iend)  + cloudhi(istart:iend)
      alb_ave(istart:iend)  = alb_ave(istart:iend)+swrsave(istart:iend)*albvisnir(istart:iend,1) &
                             +(1.-swrsave(istart:iend))*albvisnir(istart:iend,2)
      fbeam_ave(istart:iend)= fbeam_ave(istart:iend)+fbeamvis(istart:iend)*swrsave(istart:iend) &
                             +fbeamnir(istart:iend)*(1.-swrsave(istart:iend))
    endif   ! (ktau>0)
    
    ! Store fraction of direct radiation in urban scheme
    dumfbeam=fbeamvis(istart:iend)*swrsave(istart:iend)+fbeamnir(istart:iend)*(1.-swrsave(istart:iend))
    call atebfbeam(istart,imax,dumfbeam,0)

  end if  ! odcalc

  if (nmaxpr==1.and.myid==0) then
    write(6,*) "seaesfrad: Solarfit"
  end if

  ! Calculate the solar using the saved amplitude.
  sg(1:imax) = sgamp(istart:iend)*coszro2(1:imax)*taudar2(1:imax)
  if(ktau>0)then ! averages not added at time zero
    sgn_ave(istart:iend)  = sgn_ave(istart:iend)  + sg(1:imax)
    where (sg(1:imax)/ ( 1. - swrsave(istart:iend)*albvisnir(istart:iend,1) &
           -(1.-swrsave(istart:iend))*albvisnir(istart:iend,2) )>120.)
      sunhours(istart:iend)=sunhours(istart:iend)+86400.
    end where
  endif  ! (ktau>0)
      
  ! Set up the CC model radiation fields
  ! slwa is negative net radiational htg at ground
  ! Note that this does not include the upward LW radiation from the surface.
  ! That is included in sflux.f
  sgsave(istart:iend) = sg(1:imax)   ! this is the repeat after solarfit
  slwa(istart:iend) = -sgsave(istart:iend)+rgsave(istart:iend)

end do  ! Row loop (j)  j=1,jl,imax/il

! Calculate net radiational cooling of atmosphere (K/s)
t(1:ifull,:)=t(1:ifull,:)-dt*rtt(1:ifull,:)

if (nmaxpr==1.and.myid==0) then
  write(6,*) "seaesfrad: Finishing SEA-ESF radiation"
end if

call END_LOG(radmisc_end)

return
end subroutine seaesfrad

subroutine longwave_driver (is, ie, js, je, Rad_time, Atmos_input, &
                            Rad_gases, Aerosol, Aerosol_props,     &
                            Cldrad_props, Cld_spec, Aerosol_diags, &
                            Lw_output)

implicit none

!--------------------------------------------------------------------
!    longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!--------------------------------------------------------------------

integer,                            intent(in)     :: is, ie, js, je
type(time_type),                    intent(in)     :: Rad_time
type(atmos_input_type),             intent(in)     :: Atmos_input  
type(radiative_gases_type),         intent(inout)  :: Rad_gases   
type(aerosol_type),                 intent(in)     :: Aerosol     
type(aerosol_properties_type),      intent(inout)  :: Aerosol_props
type(aerosol_diagnostics_type),     intent(inout)  :: Aerosol_diags
type(cldrad_properties_type),       intent(in)     :: Cldrad_props
type(cld_specification_type),       intent(in)     :: Cld_spec     
type(lw_output_type), dimension(:), intent(inout)  :: Lw_output
type(lw_diagnostics_type), save                    :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Rad_time       time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol 
!                     fields that are seen by the longwave radiation 
!                     package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(inout) variables:
!
!      Aerosol_props  aerosol_properties_type variable containing the 
!                     aerosol radiative properties needed by the rad-
!                     iation package 
!      Lw_output      lw_output_type variable containing longwave 
!                     radiation output data 
!      Lw_diagnostics lw_diagnostics_type variable containing diagnostic
!                     longwave output used by the radiation diagnostics
!                     module
!  
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
        call sealw99 (is, ie, js, je, Rad_time, Atmos_input,           &
                      Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
                      Cld_spec, Aerosol_diags, Lw_output(1),           &
                      Lw_diagnostics, do_aerosol_forcing)

!---------------------------------------------------------------------

end subroutine longwave_driver

subroutine shortwave_driver (is, ie, js, je, Atmos_input, Surface,     &
                             Astro, Aerosol, Aerosol_props, Rad_gases, &
                             Cldrad_props,  Cld_spec, Sw_output,       &
                             Aerosol_diags, r) 

!---------------------------------------------------------------------
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!---------------------------------------------------------------------

implicit none

integer,                         intent(in)    :: is, ie, js, je
type(atmos_input_type),          intent(in)    :: Atmos_input     
type(surface_type),              intent(in)    :: Surface     
type(astronomy_type),            intent(in)    :: Astro           
type(radiative_gases_type),      intent(in)    :: Rad_gases   
type(aerosol_type),              intent(in)    :: Aerosol     
type(aerosol_properties_type),   intent(inout) :: Aerosol_props
type(cldrad_properties_type),    intent(in)    :: Cldrad_props
type(cld_specification_type),    intent(in)    :: Cld_spec
type(sw_output_type), dimension(:), intent(inout) :: Sw_output
type(aerosol_diagnostics_type), intent(inout)  :: Aerosol_diags
real(kind=8), dimension(:,:,:,:),        intent(inout) :: r

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Surface        surface_type variable containing the surface input
!                     fields needed by the radiation package
!      Astro          astronomy_type variable containing the astronom-
!                     ical input fields needed by the radiation package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol input
!                     data needed by the radiation package
!      Aerosol_props  aerosol_properties_type variable containing the
!                     aerosol radiative properties input data needed by
!                     the radiation package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(out) variables:
!
!      Sw_output      sw_output_type variable containing shortwave 
!                     radiation output data 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      type(sw_output_type) :: Sw_output_ad, Sw_output_std
      integer  :: naerosol_optical
      integer  :: i, j       
      logical  :: skipswrad

!---------------------------------------------------------------------
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      with_clouds  are clouds to be considered in determining
!                   the sw fluxes and heating rates ?
!      i,j          do-loop indices
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    allocate and initialize fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      Sw_output(1)%fsw   (:,:,:,:)    = 0.0_8
      Sw_output(1)%dfsw  (:,:,:,:)    = 0.0_8
      Sw_output(1)%ufsw  (:,:,:,:)    = 0.0_8
      Sw_output(1)%hsw   (:,:,:,:)    = 0.0_8
      Sw_output(1)%dfsw_dir_sfc     = 0.0_8
      Sw_output(1)%dfsw_dif_sfc     = 0.0_8
      Sw_output(1)%ufsw_dif_sfc     = 0.0_8
      Sw_output(1)%dfsw_vis_sfc     = 0._8
      Sw_output(1)%ufsw_vis_sfc     = 0._8
      Sw_output(1)%dfsw_vis_sfc_dir = 0._8
      Sw_output(1)%dfsw_vis_sfc_dif = 0._8
      Sw_output(1)%ufsw_vis_sfc_dif = 0._8
      Sw_output(1)%bdy_flx(:,:,:,:)   = 0.0_8       

!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then
        Sw_output(1)%fswcf (:,:,:,:) = 0.0_8
        Sw_output(1)%dfswcf(:,:,:,:) = 0.0_8
        Sw_output(1)%ufswcf(:,:,:,:) = 0.0_8
        Sw_output(1)%hswcf (:,:,:,:) = 0.0_8
        Sw_output(1)%dfsw_dir_sfc_clr = 0.0_8
        Sw_output(1)%dfsw_dif_sfc_clr  = 0.0_8
        Sw_output(1)%bdy_flx_clr (:,:,:,:) = 0.0_8
      endif

!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
          if (do_aerosol_forcing) then
            naerosol_optical = size (Aerosol_props%aerextband,2)
          else
            naerosol_optical = 0  
          endif 
          call swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,    &
                       Aerosol, Aerosol_props, Astro, Cldrad_props,        &
                       Cld_spec, include_volcanoes,                        &
                       Sw_output(1), Aerosol_diags, r,                     &
                       do_aerosol_forcing, naerosol_optical)
!--------------------------------------------------------------------

end subroutine shortwave_driver

subroutine getqsat(ilen,qsatout,temp,psin)

use estab

implicit none

integer, intent(in) :: ilen
integer iq
real, dimension(ilen), intent(in) :: temp,psin
real, dimension(ilen), intent(out) :: qsatout
real, dimension(ilen) :: esatf

do iq=1,ilen
  esatf(iq)=establ(temp(iq))
end do
qsatout=0.622*esatf/max(psin-esatf,0.1)

return
end subroutine getqsat

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlg,qfg,prf,ttg,cdrop,imax,kl)

implicit none

include 'parm.h'

integer, intent(in) :: imax,kl
integer iq,k,kr
real, dimension(imax,kl), intent(in) :: cfrac,qlg,qfg,prf,ttg
real, dimension(imax,kl), intent(in) :: cdrop
real(kind=8), dimension(imax,kl), intent(out) :: Rdrop,Rice,conl,coni
real, dimension(imax,kl) :: reffl,reffi,Wliq,rhoa,cfl,cfi
real, dimension(imax,kl) :: eps,rk,Wice
real, parameter :: scale_factor = 1. ! account for the plane-parallel homo-
                                     ! genous cloud bias  (e.g. Cahalan effect)
logical, parameter :: do_brenguier = .false. ! Adjust effective radius for vertically
                                             ! stratified cloud

rhoa=prf/(rdry*ttg)
cfl=cfrac*qlg/max(qlg+qfg,1.E-10)
cfi=max(cfrac-cfl,0.)

! Reffl is the effective radius calculated following
! Martin etal 1994, JAS 51, 1823-1842
where ( qlg>1.E-10 .and. cfl>0. .and. cfrac>0. )
  Wliq=rhoa*qlg/cfrac !kg/m^3
  ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
  !eps = 1.-0.7*exp(-0.008e-6*cdrop)  !upper bound
  eps = 1.-0.7*exp(-0.003e-6*cdrop)   !mid range
  !eps = 1.-0.7*exp(-0.001e-6*cdrop)  !lower bound
  rk  = (1.+eps**2)/(1.+2.*eps**2)**2
  
  ! k_ratio = rk**(-1./3.)  
  ! GFDL        k_ratio (land) 1.143 (water) 1.077
  ! mid range   k_ratio (land) 1.393 (water) 1.203
  ! lower bound k_ratio (land) 1.203 (water) 1.050

  ! Martin et al 1994
  reffl=(3.*(rhoa*qlg/cfl)/(4.*pi*rhow*rk*cdrop))**(1./3.)
elsewhere
  reffl=0.
  Wliq=0.
end where

! (GFDL NOTES)
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly 
!    stratified cloud in liquid specific humidity, the cloud top 
!    effective radius is greater than the effective radius of the 
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!    single layer liquid or mixed phase clouds.
if (do_brenguier) then
  if (nmr>0) then
    where (cfrac(:,2)==0.)
      !reffl(:,1)=reffl(:,1)*1.134
      reffl(:,1)=reffl(:,1)*1.2599
    end where
    do k=2,kl-1
      where (cfrac(:,k-1)==0..and.cfrac(:,k+1)==0.)
        !reffl(:,k)=reffl(:,k)*1.134
        reffl(:,k)=reffl(:,k)*1.2599
      end where
    end do
    where (cfrac(:,kl-1)==0.)
      !reffl(:,kl)=reffl(:,kl)*1.134
      reffl(:,kl)=reffl(:,kl)*1.2599
    end where 
  else
    !reffl=reffl*1.134
    reffl=reffl*1.2599
  end if
end if

!Lohmann et al.(1999)
!where ( qfg>1.E-10 .and. cfi>0. .and. cfrac>0. )
!  Wice=rhoa*qfg/cfrac !kg/m**3
!  reffi=min(150.e-6,3.73e-4*(rhoa*qfg/cfi)**0.216) 
!elsewhere
!  Wice=0.
!  reffi=0.
!end where

!Donner et al (1997)
do k=1,kl
  do iq=1,imax
    if (qfg(iq,k)>1.E-10.and.cfi(iq,k)>0..and.cfrac(iq,k)>0.) then
      Wice(iq,k)=rhoa(iq,k)*qfg(iq,k)/cfrac(iq,k) ! kg/m**3
      if (ttg(iq,k)>248.16) then
        reffi(iq,k)=5.E-7*100.6
      elseif (ttg(iq,k)>243.16) then
        reffi(iq,k)=5.E-7*80.8
      elseif (ttg(iq,k)>238.16) then
        reffi(iq,k)=5.E-7*93.5
      elseif (ttg(iq,k)>233.16) then
        reffi(iq,k)=5.E-7*63.9
      elseif (ttg(iq,k)>228.16) then
        reffi(iq,k)=5.E-7*42.5
      elseif (ttg(iq,k)>223.16) then
        reffi(iq,k)=5.E-7*39.9
      elseif (ttg(iq,k)>218.16) then
        reffi(iq,k)=5.E-7*21.6
      else
        reffi(iq,k)=5.E-7*20.2
      end if
    else
      reffi(iq,k)=0.
      Wice(iq,k)=0.
    end if
  end do
end do

do k=1,kl
  kr=kl+1-k
  Rdrop(:,kr)=real(2.E6*reffl(:,k),8) ! convert to diameter and microns
  Rice(:,kr) =real(2.E6*reffi(:,k),8)
  conl(:,kr) =real(1000.*scale_factor*Wliq(:,k),8) !g/m^3
  coni(:,kr) =real(1000.*scale_factor*Wice(:,k),8)
end do

Rdrop=min(max(Rdrop,8.4_8),33.2_8) ! constrain diameter to acceptable range (see microphys_rad.f90)
Rice=min(max(Rice,18.6_8),130.2_8)

return
end subroutine cloud3

subroutine loadaerooptical(Aerosol_props)

use cc_mpi

implicit none

include 'filnames.h'

integer :: n, nmodel, unit, num_wavenumbers, num_input_categories
integer :: noptical, nivl3, nband, nw, ierr, na, ni
integer, dimension(:), allocatable, save :: endaerwvnsf
integer, dimension(:), allocatable, save :: nivl1aero,nivl2aero
real(kind=8) :: sumsol3
real(kind=8), dimension(:,:), allocatable, save :: aeroextivl, aerossalbivl, aeroasymmivl
real(kind=8), dimension(:,:), allocatable, save :: sflwwts, sflwwts_cn
real(kind=8), dimension(:,:), allocatable, save :: solivlaero
real(kind=8), dimension(:), allocatable, save :: aeroext_in, aerossalb_in, aeroasymm_in
character(len=64), dimension(naermodels) :: aerosol_optical_names
character(len=64) :: name_in
character(len=110) :: filename
type(aerosol_properties_type), intent(inout) :: Aerosol_props

aerosol_optical_names( 1: 4)=(/ "sulfate_30%_100%", "sulfate_35%_100%", "sulfate_40%_100%", "sulfate_45%_100%" /)
aerosol_optical_names( 5: 8)=(/ "sulfate_50%_100%", "sulfate_55%_100%", "sulfate_60%_100%", "sulfate_65%_100%" /)
aerosol_optical_names( 9:12)=(/ "sulfate_70%_100%", "sulfate_75%_100%", "sulfate_80%_100%", "sulfate_82%_100%" /)
aerosol_optical_names(13:16)=(/ "sulfate_84%_100%", "sulfate_86%_100%", "sulfate_88%_100%", "sulfate_90%_100%" /)
aerosol_optical_names(17:20)=(/ "sulfate_91%_100%", "sulfate_92%_100%", "sulfate_93%_100%", "sulfate_94%_100%" /)
aerosol_optical_names(21:24)=(/ "sulfate_95%_100%", "sulfate_96%_100%", "sulfate_97%_100%", "sulfate_98%_100%" /)
aerosol_optical_names(25)=      "sulfate_99%_100%"
aerosol_optical_names(26)=      "sulfate_100%_100%"
aerosol_optical_names(27:30)=(/ "sulfate_30%_98%", "sulfate_35%_98%", "sulfate_40%_98%", "sulfate_45%_98%" /)
aerosol_optical_names(31:34)=(/ "sulfate_50%_98%", "sulfate_55%_98%", "sulfate_60%_98%", "sulfate_65%_98%" /)
aerosol_optical_names(35:38)=(/ "sulfate_70%_98%", "sulfate_75%_98%", "sulfate_80%_98%", "sulfate_82%_98%" /)
aerosol_optical_names(39:42)=(/ "sulfate_84%_98%", "sulfate_86%_98%", "sulfate_88%_98%", "sulfate_90%_98%" /)
aerosol_optical_names(43:46)=(/ "sulfate_91%_98%", "sulfate_92%_98%", "sulfate_93%_98%", "sulfate_94%_98%" /)
aerosol_optical_names(47:50)=(/ "sulfate_95%_98%", "sulfate_96%_98%", "sulfate_97%_98%", "sulfate_98%_98%" /)
aerosol_optical_names(51)=      "sulfate_99%_98%"
aerosol_optical_names(52)=      "sulfate_100%_98%"
aerosol_optical_names(53:56)=(/ "sulfate_30%_96%", "sulfate_35%_96%", "sulfate_40%_96%", "sulfate_45%_96%" /)
aerosol_optical_names(57:60)=(/ "sulfate_50%_96%", "sulfate_55%_96%", "sulfate_60%_96%", "sulfate_65%_96%" /)
aerosol_optical_names(61:64)=(/ "sulfate_70%_96%", "sulfate_75%_96%", "sulfate_80%_96%", "sulfate_82%_96%" /)
aerosol_optical_names(65:68)=(/ "sulfate_84%_96%", "sulfate_86%_96%", "sulfate_88%_96%", "sulfate_90%_96%" /)
aerosol_optical_names(69:72)=(/ "sulfate_91%_96%", "sulfate_92%_96%", "sulfate_93%_96%", "sulfate_94%_96%" /)
aerosol_optical_names(73:76)=(/ "sulfate_95%_96%", "sulfate_96%_96%", "sulfate_97%_96%", "sulfate_98%_96%" /)
aerosol_optical_names(77)=      "sulfate_99%_96%"
aerosol_optical_names(78)=      "sulfate_100%_96%"
aerosol_optical_names(79: 82)=(/ "sulfate_30%_94%", "sulfate_35%_94%", "sulfate_40%_94%", "sulfate_45%_94%" /)
aerosol_optical_names(83: 86)=(/ "sulfate_50%_94%", "sulfate_55%_94%", "sulfate_60%_94%", "sulfate_65%_94%" /)
aerosol_optical_names(87: 90)=(/ "sulfate_70%_94%", "sulfate_75%_94%", "sulfate_80%_94%", "sulfate_82%_94%" /)
aerosol_optical_names(91: 94)=(/ "sulfate_84%_94%", "sulfate_86%_94%", "sulfate_88%_94%", "sulfate_90%_94%" /)
aerosol_optical_names(95: 98)=(/ "sulfate_91%_94%", "sulfate_92%_94%", "sulfate_93%_94%", "sulfate_94%_94%" /)
aerosol_optical_names(99:102)=(/ "sulfate_95%_94%", "sulfate_96%_94%", "sulfate_97%_94%", "sulfate_98%_94%" /)
aerosol_optical_names(103)=      "sulfate_99%_94%"
aerosol_optical_names(104)=      "sulfate_100%_94%"
aerosol_optical_names(105:108)=(/ "sulfate_30%_92%", "sulfate_35%_92%", "sulfate_40%_92%", "sulfate_45%_92%" /)
aerosol_optical_names(109:112)=(/ "sulfate_50%_92%", "sulfate_55%_92%", "sulfate_60%_92%", "sulfate_65%_92%" /)
aerosol_optical_names(113:116)=(/ "sulfate_70%_92%", "sulfate_75%_92%", "sulfate_80%_92%", "sulfate_82%_92%" /)
aerosol_optical_names(117:120)=(/ "sulfate_84%_92%", "sulfate_86%_92%", "sulfate_88%_92%", "sulfate_90%_92%" /)
aerosol_optical_names(121:124)=(/ "sulfate_91%_92%", "sulfate_92%_92%", "sulfate_93%_92%", "sulfate_94%_92%" /)
aerosol_optical_names(125:128)=(/ "sulfate_95%_92%", "sulfate_96%_92%", "sulfate_97%_92%", "sulfate_98%_92%" /)
aerosol_optical_names(129)=       "sulfate_99%_92%"
aerosol_optical_names(130)=       "sulfate_100%_92%"
aerosol_optical_names(131:134)=(/ "sulfate_30%_90%", "sulfate_35%_90%", "sulfate_40%_90%", "sulfate_45%_90%" /)
aerosol_optical_names(135:138)=(/ "sulfate_50%_90%", "sulfate_55%_90%", "sulfate_60%_90%", "sulfate_65%_90%" /)
aerosol_optical_names(139:142)=(/ "sulfate_70%_90%", "sulfate_75%_90%", "sulfate_80%_90%", "sulfate_82%_90%" /)
aerosol_optical_names(143:146)=(/ "sulfate_84%_90%", "sulfate_86%_90%", "sulfate_88%_90%", "sulfate_90%_90%" /)
aerosol_optical_names(147:150)=(/ "sulfate_91%_90%", "sulfate_92%_90%", "sulfate_93%_90%", "sulfate_94%_90%" /)
aerosol_optical_names(151:154)=(/ "sulfate_95%_90%", "sulfate_96%_90%", "sulfate_97%_90%", "sulfate_98%_90%" /)
aerosol_optical_names(155)=       "sulfate_99%_90%"
aerosol_optical_names(156)=       "sulfate_100%_90%"
aerosol_optical_names(157:160)=(/ "sulfate_30%_88%", "sulfate_35%_88%", "sulfate_40%_88%", "sulfate_45%_88%" /)
aerosol_optical_names(161:164)=(/ "sulfate_50%_88%", "sulfate_55%_88%", "sulfate_60%_88%", "sulfate_65%_88%" /)
aerosol_optical_names(165:168)=(/ "sulfate_70%_88%", "sulfate_75%_88%", "sulfate_80%_88%", "sulfate_82%_88%" /)
aerosol_optical_names(169:172)=(/ "sulfate_84%_88%", "sulfate_86%_88%", "sulfate_88%_88%", "sulfate_90%_88%" /)
aerosol_optical_names(173:176)=(/ "sulfate_91%_88%", "sulfate_92%_88%", "sulfate_93%_88%", "sulfate_94%_88%" /)
aerosol_optical_names(177:180)=(/ "sulfate_95%_88%", "sulfate_96%_88%", "sulfate_97%_88%", "sulfate_98%_88%" /)
aerosol_optical_names(181)=       "sulfate_99%_88%"
aerosol_optical_names(182)=       "sulfate_100%_88%"
aerosol_optical_names(183:186)=(/ "sulfate_30%_86%", "sulfate_35%_86%", "sulfate_40%_86%", "sulfate_45%_86%" /)
aerosol_optical_names(187:190)=(/ "sulfate_50%_86%", "sulfate_55%_86%", "sulfate_60%_86%", "sulfate_65%_86%" /)
aerosol_optical_names(191:194)=(/ "sulfate_70%_86%", "sulfate_75%_86%", "sulfate_80%_86%", "sulfate_82%_86%" /)
aerosol_optical_names(195:198)=(/ "sulfate_84%_86%", "sulfate_86%_86%", "sulfate_88%_86%", "sulfate_90%_86%" /)
aerosol_optical_names(199:202)=(/ "sulfate_91%_86%", "sulfate_92%_86%", "sulfate_93%_86%", "sulfate_94%_86%" /)
aerosol_optical_names(203:206)=(/ "sulfate_95%_86%", "sulfate_96%_86%", "sulfate_97%_86%", "sulfate_98%_86%" /)
aerosol_optical_names(207)=       "sulfate_99%_86%"
aerosol_optical_names(208)=       "sulfate_100%_86%"
aerosol_optical_names(209:212)=(/ "sulfate_30%_84%", "sulfate_35%_84%", "sulfate_40%_84%", "sulfate_45%_84%" /)
aerosol_optical_names(213:216)=(/ "sulfate_50%_84%", "sulfate_55%_84%", "sulfate_60%_84%", "sulfate_65%_84%" /)
aerosol_optical_names(217:220)=(/ "sulfate_70%_84%", "sulfate_75%_84%", "sulfate_80%_84%", "sulfate_82%_84%" /)
aerosol_optical_names(221:224)=(/ "sulfate_84%_84%", "sulfate_86%_84%", "sulfate_88%_84%", "sulfate_90%_84%" /)
aerosol_optical_names(225:228)=(/ "sulfate_91%_84%", "sulfate_92%_84%", "sulfate_93%_84%", "sulfate_94%_84%" /)
aerosol_optical_names(229:232)=(/ "sulfate_95%_84%", "sulfate_96%_84%", "sulfate_97%_84%", "sulfate_98%_84%" /)
aerosol_optical_names(233)=       "sulfate_99%_84%"
aerosol_optical_names(234)=       "sulfate_100%_84%"
aerosol_optical_names(235:238)=(/ "sulfate_30%_82%", "sulfate_35%_82%", "sulfate_40%_82%", "sulfate_45%_82%" /)
aerosol_optical_names(239:242)=(/ "sulfate_50%_82%", "sulfate_55%_82%", "sulfate_60%_82%", "sulfate_65%_82%" /)
aerosol_optical_names(243:246)=(/ "sulfate_70%_82%", "sulfate_75%_82%", "sulfate_80%_82%", "sulfate_82%_82%" /)
aerosol_optical_names(247:250)=(/ "sulfate_84%_82%", "sulfate_86%_82%", "sulfate_88%_82%", "sulfate_90%_82%" /)
aerosol_optical_names(251:254)=(/ "sulfate_91%_82%", "sulfate_92%_82%", "sulfate_93%_82%", "sulfate_94%_82%" /)
aerosol_optical_names(255:258)=(/ "sulfate_95%_82%", "sulfate_96%_82%", "sulfate_97%_82%", "sulfate_98%_82%" /)
aerosol_optical_names(259)=       "sulfate_99%_82%"
aerosol_optical_names(260)=       "sulfate_100%_82%"
aerosol_optical_names(261:264)=(/ "sulfate_30%_80%", "sulfate_35%_80%", "sulfate_40%_80%", "sulfate_45%_80%" /)
aerosol_optical_names(265:268)=(/ "sulfate_50%_80%", "sulfate_55%_80%", "sulfate_60%_80%", "sulfate_65%_80%" /)
aerosol_optical_names(269:272)=(/ "sulfate_70%_80%", "sulfate_75%_80%", "sulfate_80%_80%", "sulfate_82%_80%" /)
aerosol_optical_names(273:276)=(/ "sulfate_84%_80%", "sulfate_86%_80%", "sulfate_88%_80%", "sulfate_90%_80%" /)
aerosol_optical_names(277:280)=(/ "sulfate_91%_80%", "sulfate_92%_80%", "sulfate_93%_80%", "sulfate_94%_80%" /)
aerosol_optical_names(281:284)=(/ "sulfate_95%_80%", "sulfate_96%_80%", "sulfate_97%_80%", "sulfate_98%_80%" /)
aerosol_optical_names(285)=       "sulfate_99%_80%"
aerosol_optical_names(286)=       "sulfate_100%_80%"
aerosol_optical_names(287:290)=(/ "sulfate_30%_75%", "sulfate_35%_75%", "sulfate_40%_75%", "sulfate_45%_75%" /)
aerosol_optical_names(291:294)=(/ "sulfate_50%_75%", "sulfate_55%_75%", "sulfate_60%_75%", "sulfate_65%_75%" /)
aerosol_optical_names(295:298)=(/ "sulfate_70%_75%", "sulfate_75%_75%", "sulfate_80%_75%", "sulfate_82%_75%" /)
aerosol_optical_names(299:302)=(/ "sulfate_84%_75%", "sulfate_86%_75%", "sulfate_88%_75%", "sulfate_90%_75%" /)
aerosol_optical_names(303:306)=(/ "sulfate_91%_75%", "sulfate_92%_75%", "sulfate_93%_75%", "sulfate_94%_75%" /)
aerosol_optical_names(307:310)=(/ "sulfate_95%_75%", "sulfate_96%_75%", "sulfate_97%_75%", "sulfate_98%_75%" /)
aerosol_optical_names(311)=       "sulfate_99%_75%"
aerosol_optical_names(312)=       "sulfate_100%_75%"
aerosol_optical_names(313:316)=(/ "sulfate_30%_70%", "sulfate_35%_70%", "sulfate_40%_70%", "sulfate_45%_70%" /)
aerosol_optical_names(317:320)=(/ "sulfate_50%_70%", "sulfate_55%_70%", "sulfate_60%_70%", "sulfate_65%_70%" /)
aerosol_optical_names(321:324)=(/ "sulfate_70%_70%", "sulfate_75%_70%", "sulfate_80%_70%", "sulfate_82%_70%" /)
aerosol_optical_names(325:328)=(/ "sulfate_84%_70%", "sulfate_86%_70%", "sulfate_88%_70%", "sulfate_90%_70%" /)
aerosol_optical_names(329:332)=(/ "sulfate_91%_70%", "sulfate_92%_70%", "sulfate_93%_70%", "sulfate_94%_70%" /)
aerosol_optical_names(333:336)=(/ "sulfate_95%_70%", "sulfate_96%_70%", "sulfate_97%_70%", "sulfate_98%_70%" /)
aerosol_optical_names(337)=       "sulfate_99%_70%"
aerosol_optical_names(338)=       "sulfate_100%_70%"
aerosol_optical_names(339:342)=(/ "sulfate_30%_65%", "sulfate_35%_65%", "sulfate_40%_65%", "sulfate_45%_65%" /)
aerosol_optical_names(343:346)=(/ "sulfate_50%_65%", "sulfate_55%_65%", "sulfate_60%_65%", "sulfate_65%_65%" /)
aerosol_optical_names(347:350)=(/ "sulfate_70%_65%", "sulfate_75%_65%", "sulfate_80%_65%", "sulfate_82%_65%" /)
aerosol_optical_names(351:354)=(/ "sulfate_84%_65%", "sulfate_86%_65%", "sulfate_88%_65%", "sulfate_90%_65%" /)
aerosol_optical_names(355:358)=(/ "sulfate_91%_65%", "sulfate_92%_65%", "sulfate_93%_65%", "sulfate_94%_65%" /)
aerosol_optical_names(359:362)=(/ "sulfate_95%_65%", "sulfate_96%_65%", "sulfate_97%_65%", "sulfate_98%_65%" /)
aerosol_optical_names(363)=       "sulfate_99%_65%"
aerosol_optical_names(364)=       "sulfate_100%_65%"
aerosol_optical_names(365:368)=(/ "sulfate_30%_60%", "sulfate_35%_60%", "sulfate_40%_60%", "sulfate_45%_60%" /)
aerosol_optical_names(369:372)=(/ "sulfate_50%_60%", "sulfate_55%_60%", "sulfate_60%_60%", "sulfate_65%_60%" /)
aerosol_optical_names(373:376)=(/ "sulfate_70%_60%", "sulfate_75%_60%", "sulfate_80%_60%", "sulfate_82%_60%" /)
aerosol_optical_names(377:380)=(/ "sulfate_84%_60%", "sulfate_86%_60%", "sulfate_88%_60%", "sulfate_90%_60%" /)
aerosol_optical_names(381:384)=(/ "sulfate_91%_60%", "sulfate_92%_60%", "sulfate_93%_60%", "sulfate_94%_60%" /)
aerosol_optical_names(385:388)=(/ "sulfate_95%_60%", "sulfate_96%_60%", "sulfate_97%_60%", "sulfate_98%_60%" /)
aerosol_optical_names(389)=       "sulfate_99%_60%"
aerosol_optical_names(390)=       "sulfate_100%_60%"
aerosol_optical_names(391:394)=(/ "sulfate_30%_55%", "sulfate_35%_55%", "sulfate_40%_55%", "sulfate_45%_55%" /)
aerosol_optical_names(395:398)=(/ "sulfate_50%_55%", "sulfate_55%_55%", "sulfate_60%_55%", "sulfate_65%_55%" /)
aerosol_optical_names(399:402)=(/ "sulfate_70%_55%", "sulfate_75%_55%", "sulfate_80%_55%", "sulfate_82%_55%" /)
aerosol_optical_names(403:406)=(/ "sulfate_84%_55%", "sulfate_86%_55%", "sulfate_88%_55%", "sulfate_90%_55%" /)
aerosol_optical_names(407:410)=(/ "sulfate_91%_55%", "sulfate_92%_55%", "sulfate_93%_55%", "sulfate_94%_55%" /)
aerosol_optical_names(411:414)=(/ "sulfate_95%_55%", "sulfate_96%_55%", "sulfate_97%_55%", "sulfate_98%_55%" /)
aerosol_optical_names(415)=       "sulfate_99%_55%"
aerosol_optical_names(416)=       "sulfate_100%_55%"
aerosol_optical_names(417:420)=(/ "sulfate_30%_50%", "sulfate_35%_50%", "sulfate_40%_50%", "sulfate_45%_50%" /)
aerosol_optical_names(421:424)=(/ "sulfate_50%_50%", "sulfate_55%_50%", "sulfate_60%_50%", "sulfate_65%_50%" /)
aerosol_optical_names(425:428)=(/ "sulfate_70%_50%", "sulfate_75%_50%", "sulfate_80%_50%", "sulfate_82%_50%" /)
aerosol_optical_names(429:432)=(/ "sulfate_84%_50%", "sulfate_86%_50%", "sulfate_88%_50%", "sulfate_90%_50%" /)
aerosol_optical_names(433:436)=(/ "sulfate_91%_50%", "sulfate_92%_50%", "sulfate_93%_50%", "sulfate_94%_50%" /)
aerosol_optical_names(437:440)=(/ "sulfate_95%_50%", "sulfate_96%_50%", "sulfate_97%_50%", "sulfate_98%_50%" /)
aerosol_optical_names(441)=       "sulfate_99%_50%"
aerosol_optical_names(442)=       "sulfate_100%_50%"
aerosol_optical_names(443:446)=(/ "sulfate_30%_45%", "sulfate_35%_45%", "sulfate_40%_45%", "sulfate_45%_45%" /)
aerosol_optical_names(447:450)=(/ "sulfate_50%_45%", "sulfate_55%_45%", "sulfate_60%_45%", "sulfate_65%_45%" /)
aerosol_optical_names(451:454)=(/ "sulfate_70%_45%", "sulfate_75%_45%", "sulfate_80%_45%", "sulfate_82%_45%" /)
aerosol_optical_names(455:458)=(/ "sulfate_84%_45%", "sulfate_86%_45%", "sulfate_88%_45%", "sulfate_90%_45%" /)
aerosol_optical_names(459:462)=(/ "sulfate_91%_45%", "sulfate_92%_45%", "sulfate_93%_45%", "sulfate_94%_45%" /)
aerosol_optical_names(463:466)=(/ "sulfate_95%_45%", "sulfate_96%_45%", "sulfate_97%_45%", "sulfate_98%_45%" /)
aerosol_optical_names(467)=       "sulfate_99%_45%"
aerosol_optical_names(468)=       "sulfate_100%_45%"
aerosol_optical_names(469:472)=(/ "sulfate_30%_40%", "sulfate_35%_40%", "sulfate_40%_40%", "sulfate_45%_40%" /)
aerosol_optical_names(473:476)=(/ "sulfate_50%_40%", "sulfate_55%_40%", "sulfate_60%_40%", "sulfate_65%_40%" /)
aerosol_optical_names(477:480)=(/ "sulfate_70%_40%", "sulfate_75%_40%", "sulfate_80%_40%", "sulfate_82%_40%" /)
aerosol_optical_names(481:484)=(/ "sulfate_84%_40%", "sulfate_86%_40%", "sulfate_88%_40%", "sulfate_90%_40%" /)
aerosol_optical_names(485:488)=(/ "sulfate_91%_40%", "sulfate_92%_40%", "sulfate_93%_40%", "sulfate_94%_40%" /)
aerosol_optical_names(489:492)=(/ "sulfate_95%_40%", "sulfate_96%_40%", "sulfate_97%_40%", "sulfate_98%_40%" /)
aerosol_optical_names(493)=       "sulfate_99%_40%"
aerosol_optical_names(494)=       "sulfate_100%_40%"
aerosol_optical_names(495:498)=(/ "sulfate_30%_35%", "sulfate_35%_35%", "sulfate_40%_35%", "sulfate_45%_35%" /)
aerosol_optical_names(499:502)=(/ "sulfate_50%_35%", "sulfate_55%_35%", "sulfate_60%_35%", "sulfate_65%_35%" /)
aerosol_optical_names(503:506)=(/ "sulfate_70%_35%", "sulfate_75%_35%", "sulfate_80%_35%", "sulfate_82%_35%" /)
aerosol_optical_names(507:510)=(/ "sulfate_84%_35%", "sulfate_86%_35%", "sulfate_88%_35%", "sulfate_90%_35%" /)
aerosol_optical_names(511:514)=(/ "sulfate_91%_35%", "sulfate_92%_35%", "sulfate_93%_35%", "sulfate_94%_35%" /)
aerosol_optical_names(515:518)=(/ "sulfate_95%_35%", "sulfate_96%_35%", "sulfate_97%_35%", "sulfate_98%_35%" /)
aerosol_optical_names(519)=       "sulfate_99%_35%"
aerosol_optical_names(520)=       "sulfate_100%_35%"
aerosol_optical_names(521:524)=(/ "sulfate_30%_30%", "sulfate_35%_30%", "sulfate_40%_30%", "sulfate_45%_30%" /)
aerosol_optical_names(525:528)=(/ "sulfate_50%_30%", "sulfate_55%_30%", "sulfate_60%_30%", "sulfate_65%_30%" /)
aerosol_optical_names(529:532)=(/ "sulfate_70%_30%", "sulfate_75%_30%", "sulfate_80%_30%", "sulfate_82%_30%" /)
aerosol_optical_names(533:536)=(/ "sulfate_84%_30%", "sulfate_86%_30%", "sulfate_88%_30%", "sulfate_90%_30%" /)
aerosol_optical_names(537:540)=(/ "sulfate_91%_30%", "sulfate_92%_30%", "sulfate_93%_30%", "sulfate_94%_30%" /)
aerosol_optical_names(541:544)=(/ "sulfate_95%_30%", "sulfate_96%_30%", "sulfate_97%_30%", "sulfate_98%_30%" /)
aerosol_optical_names(545)=       "sulfate_99%_30%"
aerosol_optical_names(546)=       "sulfate_100%_30%"
aerosol_optical_names(547:550)=(/ "sulfate_30%_25%", "sulfate_35%_25%", "sulfate_40%_25%", "sulfate_45%_25%" /)
aerosol_optical_names(551:554)=(/ "sulfate_50%_25%", "sulfate_55%_25%", "sulfate_60%_25%", "sulfate_65%_25%" /)
aerosol_optical_names(555:558)=(/ "sulfate_70%_25%", "sulfate_75%_25%", "sulfate_80%_25%", "sulfate_82%_25%" /)
aerosol_optical_names(559:562)=(/ "sulfate_84%_25%", "sulfate_86%_25%", "sulfate_88%_25%", "sulfate_90%_25%" /)
aerosol_optical_names(563:566)=(/ "sulfate_91%_25%", "sulfate_92%_25%", "sulfate_93%_25%", "sulfate_94%_25%" /)
aerosol_optical_names(567:570)=(/ "sulfate_95%_25%", "sulfate_96%_25%", "sulfate_97%_25%", "sulfate_98%_25%" /)
aerosol_optical_names(571)=       "sulfate_99%_25%"
aerosol_optical_names(572)=       "sulfate_100%_25%"
aerosol_optical_names(573:576)=(/ "sulfate_30%_20%", "sulfate_35%_20%", "sulfate_40%_20%", "sulfate_45%_20%" /)
aerosol_optical_names(577:580)=(/ "sulfate_50%_20%", "sulfate_55%_20%", "sulfate_60%_20%", "sulfate_65%_20%" /)
aerosol_optical_names(581:584)=(/ "sulfate_70%_20%", "sulfate_75%_20%", "sulfate_80%_20%", "sulfate_82%_20%" /)
aerosol_optical_names(585:588)=(/ "sulfate_84%_20%", "sulfate_86%_20%", "sulfate_88%_20%", "sulfate_90%_20%" /)
aerosol_optical_names(589:592)=(/ "sulfate_91%_20%", "sulfate_92%_20%", "sulfate_93%_20%", "sulfate_94%_20%" /)
aerosol_optical_names(593:596)=(/ "sulfate_95%_20%", "sulfate_96%_20%", "sulfate_97%_20%", "sulfate_98%_20%" /)
aerosol_optical_names(597)=       "sulfate_99%_20%"
aerosol_optical_names(598)=       "sulfate_100%_20%"
aerosol_optical_names(599:602)=(/ "sulfate_30%_15%", "sulfate_35%_15%", "sulfate_40%_15%", "sulfate_45%_15%" /)
aerosol_optical_names(603:606)=(/ "sulfate_50%_15%", "sulfate_55%_15%", "sulfate_60%_15%", "sulfate_65%_15%" /)
aerosol_optical_names(607:610)=(/ "sulfate_70%_15%", "sulfate_75%_15%", "sulfate_80%_15%", "sulfate_82%_15%" /)
aerosol_optical_names(611:614)=(/ "sulfate_84%_15%", "sulfate_86%_15%", "sulfate_88%_15%", "sulfate_90%_15%" /)
aerosol_optical_names(615:618)=(/ "sulfate_91%_15%", "sulfate_92%_15%", "sulfate_93%_15%", "sulfate_94%_15%" /)
aerosol_optical_names(619:622)=(/ "sulfate_95%_15%", "sulfate_96%_15%", "sulfate_97%_15%", "sulfate_98%_15%" /)
aerosol_optical_names(623)=       "sulfate_99%_15%"
aerosol_optical_names(624)=       "sulfate_100%_15%"
aerosol_optical_names(625:628)=(/ "sulfate_30%_10%", "sulfate_35%_10%", "sulfate_40%_10%", "sulfate_45%_10%" /)
aerosol_optical_names(629:632)=(/ "sulfate_50%_10%", "sulfate_55%_10%", "sulfate_60%_10%", "sulfate_65%_10%" /)
aerosol_optical_names(633:636)=(/ "sulfate_70%_10%", "sulfate_75%_10%", "sulfate_80%_10%", "sulfate_82%_10%" /)
aerosol_optical_names(637:640)=(/ "sulfate_84%_10%", "sulfate_86%_10%", "sulfate_88%_10%", "sulfate_90%_10%" /)
aerosol_optical_names(641:644)=(/ "sulfate_91%_10%", "sulfate_92%_10%", "sulfate_93%_10%", "sulfate_94%_10%" /)
aerosol_optical_names(645:648)=(/ "sulfate_95%_10%", "sulfate_96%_10%", "sulfate_97%_10%", "sulfate_98%_10%" /)
aerosol_optical_names(649)=       "sulfate_99%_10%"
aerosol_optical_names(650)=       "sulfate_100%_10%"
aerosol_optical_names(651:654)=(/ "sulfate_30%_5%", "sulfate_35%_5%", "sulfate_40%_5%", "sulfate_45%_5%" /)
aerosol_optical_names(655:658)=(/ "sulfate_50%_5%", "sulfate_55%_5%", "sulfate_60%_5%", "sulfate_65%_5%" /)
aerosol_optical_names(659:662)=(/ "sulfate_70%_5%", "sulfate_75%_5%", "sulfate_80%_5%", "sulfate_82%_5%" /)
aerosol_optical_names(663:666)=(/ "sulfate_84%_5%", "sulfate_86%_5%", "sulfate_88%_5%", "sulfate_90%_5%" /)
aerosol_optical_names(667:670)=(/ "sulfate_91%_5%", "sulfate_92%_5%", "sulfate_93%_5%", "sulfate_94%_5%" /)
aerosol_optical_names(671:674)=(/ "sulfate_95%_5%", "sulfate_96%_5%", "sulfate_97%_5%", "sulfate_98%_5%" /)
aerosol_optical_names(675)=       "sulfate_99%_5%"
aerosol_optical_names(676)=       "sulfate_100%_5%"
aerosol_optical_names(677:680)=(/ "sulfate_30%_0%", "sulfate_35%_0%", "sulfate_40%_0%", "sulfate_45%_0%" /)
aerosol_optical_names(681:684)=(/ "sulfate_50%_0%", "sulfate_55%_0%", "sulfate_60%_0%", "sulfate_65%_0%" /)
aerosol_optical_names(685:688)=(/ "sulfate_70%_0%", "sulfate_75%_0%", "sulfate_80%_0%", "sulfate_82%_0%" /)
aerosol_optical_names(689:692)=(/ "sulfate_84%_0%", "sulfate_86%_0%", "sulfate_88%_0%", "sulfate_90%_0%" /)
aerosol_optical_names(693:696)=(/ "sulfate_91%_0%", "sulfate_92%_0%", "sulfate_93%_0%", "sulfate_94%_0%" /)
aerosol_optical_names(697:700)=(/ "sulfate_95%_0%", "sulfate_96%_0%", "sulfate_97%_0%", "sulfate_98%_0%" /)
aerosol_optical_names(701)=       "sulfate_99%_0%"
aerosol_optical_names(702)=       "sulfate_100%_0%"
aerosol_optical_names(703)=       "organic_carbon"
aerosol_optical_names(704)=       "soot"
aerosol_optical_names(705:708)=(/ "sea_salt",    "dust_0.1",    "dust_0.2",    "dust_0.4" /)
aerosol_optical_names(709:712)=(/ "dust_0.8",    "dust_1.0",    "dust_2.0",    "dust_4.0" /)
aerosol_optical_names(713)   =    "dust_8.0"
aerosol_optical_names(714:717)=(/ "dust1",       "dust2",       "dust3",       "dust4"    /)
aerosol_optical_names(718)=       "dust5"
aerosol_optical_names(719:720)=(/ "bcphobic",    "omphobic" /)
aerosol_optical_names(721)=       "bcdry"
aerosol_optical_names(722:725)=(/ "bcphilic_30%", "bcphilic_35%", "bcphilic_40%", "bcphilic_45%" /)
aerosol_optical_names(726:729)=(/ "bcphilic_50%", "bcphilic_55%", "bcphilic_60%", "bcphilic_65%" /)
aerosol_optical_names(730:733)=(/ "bcphilic_70%", "bcphilic_75%", "bcphilic_80%", "bcphilic_82%" /)
aerosol_optical_names(734:737)=(/ "bcphilic_84%", "bcphilic_86%", "bcphilic_88%", "bcphilic_90%" /)
aerosol_optical_names(738:741)=(/ "bcphilic_91%", "bcphilic_92%", "bcphilic_93%", "bcphilic_94%" /)
aerosol_optical_names(742:745)=(/ "bcphilic_95%", "bcphilic_96%", "bcphilic_97%", "bcphilic_98%" /)
aerosol_optical_names(746)=       "bcphilic_99%"
aerosol_optical_names(747:750)=(/ "omphilic_30%", "omphilic_35%", "omphilic_40%", "omphilic_45%" /)
aerosol_optical_names(751:754)=(/ "omphilic_50%", "omphilic_55%", "omphilic_60%", "omphilic_65%" /)
aerosol_optical_names(755:758)=(/ "omphilic_70%", "omphilic_75%", "omphilic_80%", "omphilic_82%" /)
aerosol_optical_names(759:762)=(/ "omphilic_84%", "omphilic_86%", "omphilic_88%", "omphilic_90%" /)
aerosol_optical_names(763:766)=(/ "omphilic_91%", "omphilic_92%", "omphilic_93%", "omphilic_94%" /)
aerosol_optical_names(767:770)=(/ "omphilic_95%", "omphilic_96%", "omphilic_97%", "omphilic_98%" /)
aerosol_optical_names(771)=       "omphilic_99%"
aerosol_optical_names(772:775)=(/ "seasalt1_30%", "seasalt1_35%", "seasalt1_40%", "seasalt1_45%" /)
aerosol_optical_names(776:779)=(/ "seasalt1_50%", "seasalt1_55%", "seasalt1_60%", "seasalt1_65%" /)
aerosol_optical_names(780:783)=(/ "seasalt1_70%", "seasalt1_75%", "seasalt1_80%", "seasalt1_82%" /)
aerosol_optical_names(784:787)=(/ "seasalt1_84%", "seasalt1_86%", "seasalt1_88%", "seasalt1_90%" /)
aerosol_optical_names(788:791)=(/ "seasalt1_91%", "seasalt1_92%", "seasalt1_93%", "seasalt1_94%" /)
aerosol_optical_names(792:795)=(/ "seasalt1_95%", "seasalt1_96%", "seasalt1_97%", "seasalt1_98%" /)
aerosol_optical_names(796)=       "seasalt1_99%"
aerosol_optical_names(797:800)=(/ "seasalt2_30%", "seasalt2_35%", "seasalt2_40%", "seasalt2_45%" /)
aerosol_optical_names(801:804)=(/ "seasalt2_50%", "seasalt2_55%", "seasalt2_60%", "seasalt2_65%" /)
aerosol_optical_names(805:808)=(/ "seasalt2_70%", "seasalt2_75%", "seasalt2_80%", "seasalt2_82%" /)
aerosol_optical_names(809:812)=(/ "seasalt2_84%", "seasalt2_86%", "seasalt2_88%", "seasalt2_90%" /)
aerosol_optical_names(813:816)=(/ "seasalt2_91%", "seasalt2_92%", "seasalt2_93%", "seasalt2_94%" /)
aerosol_optical_names(817:820)=(/ "seasalt2_95%", "seasalt2_96%", "seasalt2_97%", "seasalt2_98%" /)
aerosol_optical_names(821)=       "seasalt2_99%"
aerosol_optical_names(822:825)=(/ "seasalt3_30%", "seasalt3_35%", "seasalt3_40%", "seasalt3_45%" /)
aerosol_optical_names(826:829)=(/ "seasalt3_50%", "seasalt3_55%", "seasalt3_60%", "seasalt3_65%" /)
aerosol_optical_names(830:833)=(/ "seasalt3_70%", "seasalt3_75%", "seasalt3_80%", "seasalt3_82%" /)
aerosol_optical_names(834:837)=(/ "seasalt3_84%", "seasalt3_86%", "seasalt3_88%", "seasalt3_90%" /)
aerosol_optical_names(838:841)=(/ "seasalt3_91%", "seasalt3_92%", "seasalt3_93%", "seasalt3_94%" /)
aerosol_optical_names(842:845)=(/ "seasalt3_95%", "seasalt3_96%", "seasalt3_97%", "seasalt3_98%" /)
aerosol_optical_names(846)=       "seasalt3_99%"
aerosol_optical_names(847:850)=(/ "seasalt4_30%", "seasalt4_35%", "seasalt4_40%", "seasalt4_45%" /)
aerosol_optical_names(851:854)=(/ "seasalt4_50%", "seasalt4_55%", "seasalt4_60%", "seasalt4_65%" /)
aerosol_optical_names(855:858)=(/ "seasalt4_70%", "seasalt4_75%", "seasalt4_80%", "seasalt4_82%" /)
aerosol_optical_names(859:862)=(/ "seasalt4_84%", "seasalt4_86%", "seasalt4_88%", "seasalt4_90%" /)
aerosol_optical_names(863:866)=(/ "seasalt4_91%", "seasalt4_92%", "seasalt4_93%", "seasalt4_94%" /)
aerosol_optical_names(867:870)=(/ "seasalt4_95%", "seasalt4_96%", "seasalt4_97%", "seasalt4_98%" /)
aerosol_optical_names(871)=       "seasalt4_99%"
aerosol_optical_names(872:875)=(/ "seasalt5_30%", "seasalt5_35%", "seasalt5_40%", "seasalt5_45%" /)
aerosol_optical_names(876:879)=(/ "seasalt5_50%", "seasalt5_55%", "seasalt5_60%", "seasalt5_65%" /)
aerosol_optical_names(880:883)=(/ "seasalt5_70%", "seasalt5_75%", "seasalt5_80%", "seasalt5_82%" /)
aerosol_optical_names(884:887)=(/ "seasalt5_84%", "seasalt5_86%", "seasalt5_88%", "seasalt5_90%" /)
aerosol_optical_names(888:891)=(/ "seasalt5_91%", "seasalt5_92%", "seasalt5_93%", "seasalt5_94%" /)
aerosol_optical_names(892:895)=(/ "seasalt5_95%", "seasalt5_96%", "seasalt5_97%", "seasalt5_98%" /)
aerosol_optical_names(896)=       "seasalt5_99%"

! shortwave optical models

if (myid==0) then
  filename=trim(cnsdir) // '/Ginoux_Reddy_2005'
  unit=16
  open(unit,file=filename,iostat=ierr,status='old')
  if (ierr/=0) then
    write(6,*) "ERROR: Cannot open ",trim(filename)
    call ccmpi_abort(-1)
  end if
  write(6,*) "Loading aerosol optical properties"

  !----------------------------------------------------------------------
  !    read the dimension information contained in the input file.
  !----------------------------------------------------------------------
  read ( unit,* ) num_wavenumbers
  read ( unit,* ) num_input_categories
  
  !----------------------------------------------------------------------
  !    read wavenumber limits for aerosol parameterization bands from 
  !    the input file.
  !----------------------------------------------------------------------
  call ccmpi_bcast(num_wavenumbers,0,comm_world)
  allocate ( endaerwvnsf(num_wavenumbers) )
  allocate ( aeroextivl  (num_wavenumbers, naermodels) )
  allocate ( aerossalbivl(num_wavenumbers, naermodels) )
  allocate ( aeroasymmivl (num_wavenumbers, naermodels) )
  allocate ( aeroext_in  (num_wavenumbers ) )
  allocate ( aerossalb_in(num_wavenumbers ) )
  allocate ( aeroasymm_in(num_wavenumbers ) )
  allocate ( nivl1aero (Solar_spect%nbands) )
  allocate ( nivl2aero (Solar_spect%nbands) )
  allocate ( solivlaero(Solar_spect%nbands, num_wavenumbers))
  allocate (sflwwts(N_AEROSOL_BANDS, num_wavenumbers))
  allocate (sflwwts_cn(N_AEROSOL_BANDS_CN, num_wavenumbers))           
  
  read (unit,* )
  read (unit,* ) endaerwvnsf
 
  !----------------------------------------------------------------------
  !    match the names of optical property categories from input file with
  !    those specified in the namelist, and store the following data
  !    appropriately.
  !----------------------------------------------------------------------
  do n=1,num_input_categories
    read( unit,* ) name_in
    read( unit,* )
    read( unit,* ) aeroext_in
    read( unit,* )
    read( unit,* ) aerossalb_in
    read( unit,* )
    read( unit,* ) aeroasymm_in
    do noptical=1,naermodels-4
      if (aerosol_optical_names(noptical) == name_in) then
        write(6,*) "Loading optical model for ",trim(name_in)
        aeroextivl(:,noptical)   = aeroext_in
        aerossalbivl(:,noptical) = aerossalb_in
        aeroasymmivl(:,noptical) = aeroasymm_in
        exit
      endif
    end do
  end do
  ! Dust_0.73
  aeroextivl(:,897)   = 0.175*aeroextivl(:,708)  +0.825*aeroextivl(:,709)
  aerossalbivl(:,897) = 0.175*aerossalbivl(:,708)+0.825*aerossalbivl(:,709)
  aeroasymmivl(:,897) = 0.175*aeroasymmivl(:,708)+0.825*aeroasymmivl(:,709)
  ! Dust 1.4
  aeroextivl(:,898)   = 0.6*aeroextivl(:,710)  +0.4*aeroextivl(:,711)
  aerossalbivl(:,898) = 0.6*aerossalbivl(:,710)+0.4*aerossalbivl(:,711)
  aeroasymmivl(:,898) = 0.6*aeroasymmivl(:,710)+0.4*aeroasymmivl(:,711)
  ! Dust 2.4
  aeroextivl(:,899)   = 0.8*aeroextivl(:,711)  +0.2*aeroextivl(:,712)
  aerossalbivl(:,899) = 0.8*aerossalbivl(:,711)+0.2*aerossalbivl(:,712)
  aeroasymmivl(:,899) = 0.8*aeroasymmivl(:,711)+0.2*aeroasymmivl(:,712)
  ! Dust 4.5
  aeroextivl(:,900)   = 0.875*aeroextivl(:,712)  +0.125*aeroextivl(:,713)
  aerossalbivl(:,900) = 0.875*aerossalbivl(:,712)+0.125*aerossalbivl(:,713)
  aeroasymmivl(:,900) = 0.875*aeroasymmivl(:,712)+0.125*aeroasymmivl(:,713)  

  close(unit)
  deallocate (aeroasymm_in,aerossalb_in,aeroext_in)

else
  call ccmpi_bcast(num_wavenumbers,0,comm_world)
  allocate ( endaerwvnsf(num_wavenumbers) )
  allocate ( aeroextivl  (num_wavenumbers, naermodels) )
  allocate ( aerossalbivl(num_wavenumbers, naermodels) )
  allocate ( aeroasymmivl(num_wavenumbers, naermodels) )
  allocate ( nivl1aero (Solar_spect%nbands) )
  allocate ( nivl2aero (Solar_spect%nbands) )
  allocate ( solivlaero(Solar_spect%nbands, num_wavenumbers) )
  allocate (sflwwts(N_AEROSOL_BANDS, num_wavenumbers) )
  allocate (sflwwts_cn(N_AEROSOL_BANDS_CN, num_wavenumbers) )  
end if

call ccmpi_bcast(endaerwvnsf,0,comm_world)
call ccmpi_bcastr8(aeroextivl,0,comm_world)
call ccmpi_bcastr8(aerossalbivl,0,comm_world)
call ccmpi_bcastr8(aeroasymmivl,0,comm_world)

!---------------------------------------------------------------------
!    define the solar weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the solar
!    spectral intervals and so determine the single-scattering proper-
!    ties on the solar spectral intervals.
!--------------------------------------------------------------------
nivl3 = 1
sumsol3 = 0.0_8
nband = 1
solivlaero(:,:) = 0.0_8
nivl1aero(1) = 1
do nw = 1,Solar_spect%endwvnbands(Solar_spect%nbands)
  sumsol3 = sumsol3 + Solar_spect%solarfluxtoa(nw)
  if (nw == endaerwvnsf(nivl3) ) then
    solivlaero(nband,nivl3) = sumsol3
    sumsol3 = 0.0_8
  end if
  if ( nw == Solar_spect%endwvnbands(nband) ) then
    if ( nw /= endaerwvnsf(nivl3) ) then
      solivlaero(nband,nivl3) = sumsol3 
      sumsol3 = 0.0_8
    end if
    nivl2aero(nband) = nivl3
    nband = nband + 1
    if ( nband <= Solar_spect%nbands ) then
      if ( nw == endaerwvnsf(nivl3) ) then
        nivl1aero(nband) = nivl3 + 1
      else
        nivl1aero(nband) = nivl3
      end if
    end if
  end if
  if ( nw == endaerwvnsf(nivl3) ) nivl3 = nivl3 + 1
end do

Aerosol_props%aerextband=0._8
Aerosol_props%aerssalbband=0._8
Aerosol_props%aerasymmband=0._8

do nmodel=1,naermodels
  call thickavg (nivl1aero, nivl2aero, num_wavenumbers,    &
                 Solar_spect%nbands, aeroextivl(:,nmodel), &
                 aerossalbivl(:,nmodel),                   &
                 aeroasymmivl(:,nmodel), solivlaero,       &
                 Solar_spect%solflxbandref,                & 
                 Aerosol_props%aerextband(:,nmodel),       &
                 Aerosol_props%aerssalbband(:,nmodel),     &
                 Aerosol_props%aerasymmband(:,nmodel))
end do

! longwave optical models
call lw_aerosol_interaction(num_wavenumbers,sflwwts,sflwwts_cn,endaerwvnsf)

!    the units of extinction coefficient (aeroextivl) are m**2/gm.
!    to make the lw band extinction coefficient (aerextbandlw) have
!    units (m**2/Kg) consistent with the units in FMS models, one
!    must multiply by 1000. this is done below.
  
Aerosol_props%aerextbandlw=0._8
Aerosol_props%aerssalbbandlw=0._8
Aerosol_props%aerextbandlw_cn=0._8
Aerosol_props%aerssalbbandlw_cn=0._8

do nw=1,naermodels    
  do na=1,N_AEROSOL_BANDS  
    do ni=1,num_wavenumbers 
      Aerosol_props%aerextbandlw(na,nw) =    &
                      Aerosol_props%aerextbandlw(na,nw) + &
                      aeroextivl(ni,nw)*sflwwts(na,ni)*  &
                      1.0E+03_8
      Aerosol_props%aerssalbbandlw(na,nw) =     &
                      Aerosol_props%aerssalbbandlw(na,nw) +&
                      aerossalbivl(ni,nw)*sflwwts(na,ni)
    end do
  end do
end do
do nw=1,naermodels    
  do na=1,N_AEROSOL_BANDS_CN
    do ni=1,num_wavenumbers 
      Aerosol_props%aerextbandlw_cn(na,nw) =    &
                      Aerosol_props%aerextbandlw_cn(na,nw) + &
                      aeroextivl(ni,nw)*sflwwts_cn(na,ni)*  &
                      1.0E+03_8
      Aerosol_props%aerssalbbandlw_cn(na,nw) =     &
                      Aerosol_props%aerssalbbandlw_cn(na,nw) +&
                      aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
    end do
  end do
end do

deallocate (sflwwts_cn,sflwwts)
deallocate (solivlaero,nivl2aero,nivl1aero) 
deallocate (aeroasymmivl,aerossalbivl,aeroextivl)
deallocate (endaerwvnsf)

return
end subroutine loadaerooptical

subroutine lw_aerosol_interaction(num_wavenumbers,sflwwts,sflwwts_cn,endaerwvnsf)      

use longwave_params_mod

implicit none

integer, intent(in) :: num_wavenumbers
integer, dimension(num_wavenumbers), intent(in) :: endaerwvnsf
real(kind=8), dimension(N_AEROSOL_BANDS, num_wavenumbers), intent(out) :: sflwwts
real(kind=8), dimension(N_AEROSOL_BANDS_CN, num_wavenumbers), intent(out) :: sflwwts_cn

!----------------------------------------------------------------------
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    the following arrays define the wavenumber ranges for the separate
!    aerosol emissivity bands in the model infrared parameterization. 
!    these may be changed only by the keeper of the radiation code.
!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code.
!
!      aerbandlo_fr      low wavenumber limit for the non-continuum 
!                        aerosol emissivity bands
!      aerbandhi_fr      high wavenumber limit for the non-continuum
!                        aerosol emissivity bands
!      istartaerband_fr  starting wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      iendaerband_fr    ending wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      aerbandlo_co      low wavenumber limit for the continuum 
!                        aerosol emissivity bands
!      aerbandhi_co      high wavenumber limit for the continuum
!                        aerosol emissivity bands
!      istartaerband_co  starting wavenumber index for the continuum
!                        aerosol emissivity bands
!      iendaerband_co    ending wavenumber index for the continuum
!                        aerosol emissivity bands
!      aerbandlo         low wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      aerbandhi         high wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      istartaerband     starting wavenumber index for the entire set of
!                        aerosol emissivity bands
!      iendaerband       ending wavenumber index for the entire set of
!                        aerosol emissivity bands
!
!----------------------------------------------------------------------
      real(kind=8), dimension (N_AEROSOL_BANDS_FR)     :: aerbandlo_fr =  &
      (/ 560.0_8, 630.0_8, 700.0_8, 800.0_8, 900.0_8,  990.0_8, 1070.0_8, 1200.0_8 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_FR)     :: aerbandhi_fr =  &
      (/ 630.0_8, 700.0_8, 800.0_8, 900.0_8, 990.0_8, 1070.0_8, 1200.0_8, 1400.0_8 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CO)     :: aerbandlo_co =  &
      (/ 560.0_8 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CO)     :: aerbandhi_co =  &
      (/ 800.0_8 /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: istartaerband_co =  &
      (/ 57  /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: iendaerband_co =  &
      (/ 80  /)
      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandlo_cn =  &
      (/ 800.0_8 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandhi_cn =  &
      (/ 1200.0_8 /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: istartaerband_cn =  &
      (/ 81  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: iendaerband_cn =  &
      (/ 120 /)

      real(kind=8),    dimension(N_AEROSOL_BANDS)      :: aerbandlo, aerbandhi
      integer, dimension(N_AEROSOL_BANDS)      :: istartaerband,    &
                                                  iendaerband

!---------------------------------------------------------------------
!    the following arrays define how the ir aerosol band structure 
!    relates to the aerosol parameterization bands.
!
!      nivl1aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl2aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl1aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl2aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl1aer(n)       aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber for the 
!                        ir aerosol emissivity band n
!      nivl2aer(n)       aerosol parameterization band index corres-
!                        ponding to the highest wavenumber for the 
!                        ir aerosol emissivity band n
!      planckaerband(n)  planck function summed over each lw param-
!                        eterization band that is contained in the 
!                        ir aerosol emissivity band n
!
!---------------------------------------------------------------------
      integer, dimension (N_AEROSOL_BANDS_FR)  :: nivl1aer_fr,   &
                                                  nivl2aer_fr
      integer, dimension (N_AEROSOL_BANDS_CO)  :: nivl1aer_co,   &
                                                  nivl2aer_co
      integer, dimension (N_AEROSOL_BANDS_CN)  :: nivl1aer_cn,   &
                                                  nivl2aer_cn
      integer, dimension (N_AEROSOL_BANDS)     :: nivl1aer, nivl2aer
      real(kind=8),    dimension (N_AEROSOL_BANDS)     :: planckaerband
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN)  :: planckaerband_cn

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
!        aerextbandlw_fr(n,na)  band averaged extinction coefficient
!                               for non-continuum aerosol emissivity 
!                               band n and aerosol properties type na
!        aerssalbbandlw_fr(n,na)
!                               band averaged single-scattering
!                               coefficient for non-continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        aerextbandlw_co(n,na)  band averaged extinction coefficient
!                               for the continuum aerosol emissivity
!                               band n and aerosol properties type na
!        aerssalbbandlw_co(n,na)
!                               band averaged single-scattering
!                               coefficient for continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        planckivlaer_fr(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity non-
!                               continuum band n and aerosol parameter-
!                               ization band ni
!        planckivlaer_co(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity continuum 
!                               band n and aerosol parameterization 
!                               band ni
!        sflwwts_fr(n,ni)       band weights for the aerosol emissivity
!                               non-continuum band n and the aerosol 
!                               parameterization band ni 
!        sflwwts_co(n,ni)       band weights for the aerosol emissivity
!                               continuum band n and the aerosol 
!                               parameterization band ni 
!        planckivlaer(n,ni)     planck function over the spectral range
!                               common to aerosol emissivity band n and
!                               aerosol parameterization band ni
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real(kind=8),    dimension (N_AEROSOL_BANDS_FR, naermodels) ::   &
                                                  aerextbandlw_fr, &
                                                  aerssalbbandlw_fr
      real(kind=8),    dimension (N_AEROSOL_BANDS_CO, naermodels) ::   &
                                                  aerextbandlw_co, &
                                                  aerssalbbandlw_co
      real(kind=8),    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real(kind=8),    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN, num_wavenumbers) :: &
                                                  planckivlaer_cn   
      real(kind=8),    dimension (N_AEROSOL_BANDS, num_wavenumbers)  ::    &
                                                  planckivlaer
      integer, dimension (num_wavenumbers)    ::  iendsfbands

!---------------------------------------------------------------------
!    variables associated with the planck function calculation.
!    the planck function is defined for each of the NBLW longwave 
!    parameterization bands.
!---------------------------------------------------------------------
      real(kind=8), dimension(NBLW)  :: c1, centnb, sc, src1nb, x, x1
      real(kind=8)                   :: del, xtemv, sumplanck

!---------------------------------------------------------------------
!    miscellaneous variables:

     logical         :: do_band1   !  should we do special calculation 
                                   !  for band 1 ?
     integer         :: ib, nw, nivl, nband, n, ni 
                                   !  do-loop indices and counters

!--------------------------------------------------------------------
!    define arrays containing the characteristics of all the ir aerosol
!    emissivity bands, both continuum and non-continuum.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        aerbandlo(n)     = aerbandlo_fr(n)
        aerbandhi(n)     = aerbandhi_fr(n)
        istartaerband(n) = istartaerband_fr(n)
        iendaerband(n)   = iendaerband_fr(n)
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        aerbandlo(n)     = aerbandlo_co     (n - N_AEROSOL_BANDS_FR)
        aerbandhi(n)     = aerbandhi_co     (n - N_AEROSOL_BANDS_FR)
        istartaerband(n) = istartaerband_co (n - N_AEROSOL_BANDS_FR)
        iendaerband(n)   = iendaerband_co   (n - N_AEROSOL_BANDS_FR)
      end do

!---------------------------------------------------------------------
!    define the number of aerosol ir bands to be used in other modules.
!    set the initialization flag to .true.
!---------------------------------------------------------------------
      Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS
      Lw_parameters%n_lwaerosol_bands_iz = .true.

!--------------------------------------------------------------------
!    define the ending aerosol band index for each of the aerosol
!    parameterization bands.
!--------------------------------------------------------------------
      iendsfbands(:) = INT((endaerwvnsf(:) + 0.01_8)/10.0_8)

!--------------------------------------------------------------------
!    compute the planck function at 10C over each of the longwave
!    parameterization bands to be used as the weighting function. 
!--------------------------------------------------------------------
      do n=1,NBLW 
        del  = 10.0E+00_8
        xtemv = 283.15_8
        centnb(n) = 5.0_8 + real(n - 1,8)*del
        c1(n)     = (3.7412E-05_8)*centnb(n)**3
        x(n)      = 1.4387E+00_8*centnb(n)/xtemv
        x1(n)     = EXP(x(n))
        sc(n)     = c1(n)/(x1(n) - 1.0E+00_8)
        src1nb(n) = del*sc(n)
      end do
 
!--------------------------------------------------------------------
!    sum the weighting function calculated over the longwave param-
!    eterization bands that are contained in each of the aerosol 
!    emissivity bands. 
!--------------------------------------------------------------------
      planckaerband(:) = 0.0E+00_8
      do n = 1,N_AEROSOL_BANDS
        do ib = istartaerband(n),iendaerband(n)
          planckaerband(n) = planckaerband(n) + src1nb(ib)
        end do
      end do
      planckaerband_cn(:) = 0.0E+00_8
      do n = 1,N_AEROSOL_BANDS_CN
        do ib = istartaerband_cn(n),iendaerband_cn(n)
          planckaerband_cn(n) = planckaerband_cn(n) + src1nb(ib)
        end do
      end do
 
!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the non-
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0_8
      nband = 1
      planckivlaer_fr(:,:) = 0.0_8
      nivl1aer_fr(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_fr(nband,nivl) = sumplanck
          sumplanck = 0.0_8
        end if
        if ( nw == iendaerband_fr(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck 
            sumplanck = 0.0_8
          end if
          nivl2aer_fr(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_FR ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_fr(nband) = nivl + 1
            else
              nivl1aer_fr(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.   &
              iendsfbands(nivl-1) >= istartaerband_fr(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_fr(1)) then
            nivl1aer_fr(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if (nw >= iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0_8
      nband = 1
      planckivlaer_co(:,:) = 0.0_8
      nivl1aer_co(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_co(nband,nivl) = sumplanck
          sumplanck = 0.0_8
        end if
        if ( nw == iendaerband_co(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck 
            sumplanck = 0.0_8
          end if
          nivl2aer_co(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CO ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_co(nband) = nivl + 1
            else
              nivl1aer_co(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_co(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_co(1)) then
            nivl1aer_co(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_co(N_AEROSOL_BANDS_CO) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0_8
      nband = 1
      planckivlaer_cn(:,:) = 0.0_8
      nivl1aer_cn(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_cn(nband,nivl) = sumplanck
          sumplanck = 0.0_8
        end if
        if ( nw == iendaerband_cn(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_cn(nband,nivl) = sumplanck 
            sumplanck = 0.0_8
          end if
          nivl2aer_cn(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CN ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_cn(nband) = nivl + 1
            else
              nivl1aer_cn(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_cn(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_cn(1)) then
            nivl1aer_cn(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_cn(N_AEROSOL_BANDS_CN) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the planck-function-weighted band weights for the aerosol
!    parameterization bands onto the non-continuum and continuum ir 
!    aerosol emissivity bands.
!--------------------------------------------------------------------
      sflwwts_fr(:,:) = 0.0E+00_8
      do n=1,N_AEROSOL_BANDS_FR
        do ni=nivl1aer_fr(n),nivl2aer_fr(n)
          sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
        end do
      end do
      sflwwts_co(:,:) = 0.0E+00_8
      do n=1,N_AEROSOL_BANDS_CO
        do ni=nivl1aer_co(n),nivl2aer_co(n)
          sflwwts_co(n,ni) = planckivlaer_co(n,ni)/     &
                             planckaerband(N_AEROSOL_BANDS_FR+n)
        end do
      end do
      sflwwts_cn(:,:) = 0.0E+00_8
      do n=1,N_AEROSOL_BANDS_CN
        do ni=nivl1aer_cn(n),nivl2aer_cn(n)
          sflwwts_cn(n,ni) = planckivlaer_cn(n,ni)/     &
                             planckaerband_cn(n)
        end do
      end do

!--------------------------------------------------------------------
!    consolidate the continuum and non-continuum weights into an
!    array covering all ir aerosol emissivity bands.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_fr(n,ni)
        end do
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        end do
      end do

!----------------------------------------------------------------------

end subroutine lw_aerosol_interaction


end module seaesfrad_m
