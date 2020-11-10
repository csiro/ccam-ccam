! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
! Interface for SEA-ESF radiation scheme (from GFDL AM3) with CCAM.

! Interface developed by Marcus Thatcher, Maciej Golebiewski and Paul Ryan

! This routine assumes that only one month at a time is integrated in RCM mode

module seaesfrad_m

use rad_utilities_mod, only : atmos_input_type,surface_type,astronomy_type,aerosol_type,           &
                              aerosol_properties_type,radiative_gases_type,cldrad_properties_type, &
                              cld_specification_type,lw_output_type,sw_output_type,                &
                              aerosol_diagnostics_type,time_type,microphysics_type,                &
                              microrad_properties_type,lw_diagnostics_type,lw_table_type,          &
                              Sw_control,Lw_control, Rad_control,Cldrad_control,Lw_parameters,     &
                              thickavg,lw_clouds_type,optical_path_type,gas_tf_type
use esfsw_driver_mod, only : swresf,esfsw_driver_init
use sealw99_mod, only : sealw99,sealw99_init, sealw99_time_vary
use esfsw_parameters_mod, only : Solar_spect,esfsw_parameters_init,sw_resolution,sw_diff_streams
use microphys_rad_mod, only : microphys_rad_init,microphys_sw_driver,microphys_lw_driver,          &
                              lwemiss_calc,lwem_form

private
public seaesfrad_settime, seaesfrad, seaesfrad_init, sw_resolution, sw_diff_streams, liqradmethod, iceradmethod
public carbonradmethod, so4radmethod, dustradmethod, seasaltradmethod, lwem_form
public csolar

real, parameter :: rhow     = 1000.            ! Density of water (kg/m^3)
real, save      :: csolar   = 1365.            ! Solar constant in W/m^2
real, parameter :: siglow   = 0.68             ! sigma level for top of low cloud (diagnostic)
real, parameter :: sigmid   = 0.44             ! sigma level for top of medium cloud (diagnostic)
real, parameter :: ratco2mw = 1.519449738      ! conversion factor for CO2 diagnostic
integer, parameter :: naermodels         = 218 ! number of aerosol optical models
integer, parameter :: N_AEROSOL_BANDS_FR = 8
integer, parameter :: N_AEROSOL_BANDS_CO = 1
integer, parameter :: N_AEROSOL_BANDS_CN = 1
integer, parameter :: N_AEROSOL_BANDS    = N_AEROSOL_BANDS_FR + N_AEROSOL_BANDS_CO
integer, parameter :: nfields            = 11 ! number of aerosol fields for radiation
integer, save :: liqradmethod = 0     ! Method for calculating radius of liquid droplets
                                      ! (0=Martin Mid/NBer, 1=Martin Low/NBer, 2=Martin High/NBer
                                      !  3=Martin Mid/Ber,  4=Martin Low/Ber,  5=Martin High/Ber )
integer, save :: iceradmethod = 1     ! Method for calculating radius of ice droplets
                                      ! (0=Lohmann, 1=Donner smooth, 2=Fu, 3=Donner orig)
integer, save :: so4radmethod     = 0 ! Method for SO4 direct effects      (0=on, -1=off)
integer, save :: carbonradmethod  = 0 ! Method for carbon direct effects   (0=phobic/phillic, 1=generic, -1=off)
integer, save :: dustradmethod    = 0 ! Method for dust direct effects     (0=on, -1=off)
integer, save :: seasaltradmethod = 0 ! Method for sea-salt direct effects (0=on, -1=off)
logical, parameter :: do_totcld_forcing  = .true.
logical, parameter :: include_volcanoes  = .false.
logical, save :: do_aerosol_forcing ! =.true. when abs(iaero)>=2

integer, save :: nlow, nmid
real(kind=8), dimension(:,:), allocatable, save :: pref

integer, save :: mins
real, save :: r1, dlt, alp, slag, fjd

type(time_type), save ::                                           Rad_time
type(aerosol_properties_type), save ::                             Aerosol_props
type(lw_table_type), save ::                                       Lw_tables
type(atmos_input_type), dimension(:), allocatable, save ::         Atmos_input
type(surface_type), dimension(:), allocatable, save ::             Surface     
type(astronomy_type), dimension(:), allocatable, save ::           Astro
type(aerosol_type), dimension(:), allocatable, save ::             Aerosol
type(radiative_gases_type), dimension(:), allocatable, save ::     Rad_gases
type(cldrad_properties_type), dimension(:), allocatable, save ::   Cldrad_props
type(cld_specification_type), dimension(:), allocatable, save ::   Cld_spec
type(microphysics_type), dimension(:), allocatable, save ::        Cloud_microphysics
type(microrad_properties_type), dimension(:), allocatable, save :: Lscrad_props
type(lw_output_type), dimension(:), allocatable, save ::           Lw_output
type(sw_output_type), dimension(:), allocatable, save ::           Sw_output
type(aerosol_diagnostics_type), dimension(:), allocatable, save :: Aerosol_diags
type(lw_diagnostics_type), dimension(:), allocatable, save ::      Lw_diagnostics
type(lw_clouds_type), dimension(:), allocatable, save ::           Lw_clouds
type(optical_path_type), dimension(:), allocatable, save ::        Optical
type(gas_tf_type), dimension(:), allocatable, save ::              Gas_tf

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation step for GFDL SEA-ESF radiation (single thread only)
!

subroutine seaesfrad_settime

use infile                                          ! Input file routines
use nharrs_m                                        ! Non-hydrostatic atmosphere arrays
use parm_m                                          ! Model configuration
use raddiag_m                                       ! Radiation diagnostic
use zenith_m                                        ! Astronomy routines

implicit none

integer jyear, jmonth, jday, jhour, jmin

! True for full radiation calculation
odcalc = mod(ktau,kountr)==0 .or. (ktau==1.and.((.not.lrestart_radiation).or.always_mspeca))

! astronomy ---------------------------------------------------------
! Set up number of minutes from beginning of year
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
fjd = float(mod(mins, 525600))/1440. ! restrict to 365 day calendar
! Calculate sun position
call solargh(fjd,bpyear,r1,dlt,alp,slag)

! Prepare SEA-ESF arrays --------------------------------------------
Rad_time%days    = int(fjd)
Rad_time%seconds = mod(mins, 1440)*60
Rad_time%ticks   = 0

return
end subroutine seaesfrad_settime
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CCAM interface with GFDL SEA-ESF radiation
!

subroutine seaesfrad(koundiag)

use aerointerface                                   ! Aerosol interface
use aerosolldr                                      ! LDR prognostic aerosols
use arrays_m                                        ! Atmosphere dyamics prognostic arrays
use ateb, only : atebccangle,atebalb1,atebfbeam     ! Urban
use cc_mpi                                          ! CC MPI routines
use cc_omp                                          ! CC OpenMP routines
use cfrac_m                                         ! Cloud fraction
use const_phys                                      ! Physical constants
use extraout_m                                      ! Additional diagnostics
use estab                                           ! Liquid saturation function
use infile                                          ! Input file routines
use latlong_m                                       ! Lat/lon coordinates
use mlo, only : mloalb4                             ! Ocean physics and prognostic arrays
use newmpar_m                                       ! Grid parameters
use nharrs_m                                        ! Non-hydrostatic atmosphere arrays
use nsibd_m                                         ! Land-surface arrays
use ozoneread                                       ! Ozone input routines
use parm_m                                          ! Model configuration
use pbl_m                                           ! Boundary layer arrays
use raddiag_m                                       ! Radiation diagnostic
use sigs_m                                          ! Atmosphere sigma levels
use soil_m                                          ! Soil and surface data
use soilsnow_m                                      ! Soil, snow and surface data
use work3f_m                                        ! Grid work arrays
use zenith_m                                        ! Astronomy routines

implicit none

include 'kuocom.h'                                  ! Convection parameters

integer, intent(inout) :: koundiag
integer k
integer i, iq, istart, iend, kr, nr, iq_tile
integer ktop, kbot, mythread
real, dimension(imax,kl) :: duo3n, rhoa
real, dimension(imax,kl) :: p2, cd2, dumcf, dumql, dumqf, dumt, dz
real, dimension(imax) :: coszro2, taudar2, coszro, taudar, mx
real, dimension(imax) :: sgdnvis, sgdnnir
real, dimension(imax) :: sgvis, sgdnvisdir, sgdnvisdif, sgdnnirdir, sgdnnirdif
real, dimension(imax) :: dzrho, dumtss, alb
real, dimension(imax) :: cuvrf_dir, cirrf_dir, cuvrf_dif, cirrf_dif, fbeam, sgn_save
real, dimension(imax) :: sgn
real, dimension(kl+1) :: sigh
real, dimension(kl) :: diag_temp
real dhr, cosz, delta

real(kind=8), dimension(1,1,1,1) :: r

if ( diag .and. mydiag ) then
  diag_temp(:) = t(idjd,:)
  write(6,*) "tdiag ",diag_temp
  diag_temp(:) = qg(idjd,:)
  write(6,*) "qgdiag ",diag_temp
  diag_temp(:) = qlrad(idjd,:)
  write(6,*) "qlraddiag ",diag_temp
  diag_temp(:) = qfrad(idjd,:)
  write(6,*) "qfraddiag ",diag_temp
  if ( abs(iaero)>=2 ) then
    diag_temp(:) = xtg(idjd,:,3)  
    write(6,*) "SO4diag ",diag_temp
    diag_temp(:) = xtg(idjd,:,4)
    write(6,*) "BCphobdiag ",diag_temp
    diag_temp(:) = xtg(idjd,:,5)
    write(6,*) "BCphildiag ",diag_temp
    diag_temp(:) = xtg(idjd,:,6)
    write(6,*) "OCphobdiag ",diag_temp
    diag_temp(:) = xtg(idjd,:,7)
    write(6,*) "OCphildiag ",diag_temp
    diag_temp(:) = xtg(idjd,:,8)
    write(6,*) "dust0.8diag ",diag_temp
    diag_temp(:) = xtg(idjd,:,9)
    write(6,*) "dust1.0diag ",diag_temp
    diag_temp(:) = xtg(idjd,:,10)
    write(6,*) "dust2.0diag ",diag_temp
    diag_temp(:) = xtg(idjd,:,11)
    write(6,*) "dust4.0diag ",diag_temp
    diag_temp(:) = xtg(idjd,:,12)
    write(6,*) "saltfilmdiag ",diag_temp
    diag_temp(:) = xtg(idjd,:,13)
    write(6,*) "saltjetdiag  ",diag_temp
  end if
end if


! main loop ---------------------------------------------------------
!$omp do schedule(static) private(istart,iend,mythread,dhr,coszro2,taudar2,coszro,taudar), &
!$omp private(sigh,duo3n,cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif,rhoa,dz,i,iq,cosz),       &
!$omp private(delta,k,kr,dzrho,cd2,mx,dumt,p2,ktop,kbot,dumcf),                            &
!$omp private(dumql,dumqf,sgdnvis,sgdnnir,sgvis,sgdnvisdir,sgdnvisdif,sgdnnirdir),         &
!$omp private(sgdnnirdif,rt,fbeam,sgn_save,sgn)
do iq_tile = 1,ifull,imax
  istart = iq_tile
  iend   = istart + imax - 1
  
  call ccomp_mythread(mythread)
  
  ! Calculate zenith angle for the solarfit calculation.
  dhr = dt/3600.
  call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro2,taudar2)
  call atebccangle(istart,imax,coszro2,rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt,sin(dlt))

  ! Set-up albedo
  ! Land albedo ---------------------------------------------------
  if ( nsib==6 .or. nsib==7 ) then
    ! CABLE version
    where ( land(istart:iend) )
      cuvrf_dir(1:imax) = albvisdir(istart:iend) ! from cable (inc snow)
      cirrf_dir(1:imax) = albnirdir(istart:iend) ! from cable (inc snow)
      cuvrf_dif(1:imax) = albvisdif(istart:iend) ! from cable (inc snow)
      cirrf_dif(1:imax) = albnirdif(istart:iend) ! from cable (inc snow)
    end where
  else
    ! nsib=3 version (calculate snow)
    where ( land(istart:iend) )
      cuvrf_dir(1:imax) = albvissav(istart:iend) ! from albfile (indata.f)
      cirrf_dir(1:imax) = albnirsav(istart:iend) ! from albnirfile (indata.f)
      cuvrf_dif(1:imax) = cuvrf_dir(1:imax)      ! assume DIR and DIF are the same
      cirrf_dif(1:imax) = cirrf_dir(1:imax)      ! assume DIR and DIF are the same
    end where
    call calc_snow_albedo(coszro2,cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif,iq_tile)
  end if

  ! Water/Ice albedo --------------------------------------------
  if ( nmlo==0 ) then ! prescribed SSTs
    ! NCAR CCMS3.0 scheme (Briegleb et al, 1986,
    ! J. Clim. and Appl. Met., v. 27, 214-226)
    where ( .not.land(istart:iend) .and. coszro2>=1.e-10 )
      cuvrf_dir(1:imax) = 0.026/(coszro2**1.7+0.065)                 &
        + 0.15*(coszro2-0.1)*(coszro2-0.5)*(coszro2-1.)
    elsewhere ( .not.land(istart:iend) )
      cuvrf_dir(1:imax) = 0.3925 ! coszen=0 value of above expression
    end where
    where ( .not.land(istart:iend) )
      cuvrf_dif(1:imax) = 0.06
      cirrf_dir(1:imax) = cuvrf_dir(1:imax)
      cirrf_dif(1:imax) = 0.06
      cuvrf_dir(1:imax) = 0.85*fracice(istart:iend) + (1.-fracice(istart:iend))*cuvrf_dir(1:imax)
      cuvrf_dif(1:imax) = 0.85*fracice(istart:iend) + (1.-fracice(istart:iend))*cuvrf_dif(1:imax)
      cirrf_dir(1:imax) = 0.45*fracice(istart:iend) + (1.-fracice(istart:iend))*cirrf_dir(1:imax)
      cirrf_dif(1:imax) = 0.45*fracice(istart:iend) + (1.-fracice(istart:iend))*cirrf_dif(1:imax)
    end where
  elseif (abs(nmlo)<=9) then
    ! MLO albedo ----------------------------------------------------
    call mloalb4(istart,imax,coszro2,cuvrf_dir,cuvrf_dif,cirrf_dir,cirrf_dif,0)
  else
    write(6,*) "ERROR: PCOM albedo is required in seaesfrad"
    call ccmpi_abort(-1)
  end if

  ! Urban albedo --------------------------------------------------
  call atebalb1(istart,imax,cuvrf_dir,0,split=1) ! direct
  call atebalb1(istart,imax,cirrf_dir,0,split=1) ! direct
  call atebalb1(istart,imax,cuvrf_dif,0,split=2) ! diffuse
  call atebalb1(istart,imax,cirrf_dif,0,split=2) ! diffuse
   
  ! Call radiation --------------------------------------------------
  if ( odcalc ) then     ! Do the calculation

    ! Average the zenith angle over the time (hours) between radiation
    ! calculations
    dhr = real(kountr)*dt/3600.
    call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro,taudar)

    ! Recalculate albedo for backwards compatibility
    if ( always_mspeca ) then
      if ( .not.(nsib==6.or.nsib==7) ) then
        ! nsib=3 version (calculate snow)
        where ( land(istart:iend) )
          cuvrf_dir(1:imax) = albvissav(istart:iend) ! from albfile (indata.f)
          cirrf_dir(1:imax) = albnirsav(istart:iend) ! from albnirfile (indata.f)
          cuvrf_dif(1:imax) = cuvrf_dir(1:imax)      ! assume DIR and DIF are the same
          cirrf_dif(1:imax) = cirrf_dir(1:imax)      ! assume DIR and DIF are the same
        end where
        call calc_snow_albedo(coszro,cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif,iq_tile)
      end if  
      if ( nmlo==0 ) then ! prescribed SSTs
        ! NCAR CCMS3.0 scheme (Briegleb et al, 1986,
        ! J. Clim. and Appl. Met., v. 27, 214-226)
        where ( .not.land(istart:iend) .and. coszro>=0. )
          cuvrf_dir(1:imax) = 0.026/(coszro**1.7+0.065)                 &
            + 0.15*(coszro-0.1)*(coszro-0.5)*(coszro-1.)
        elsewhere ( .not.land(istart:iend) )
         cuvrf_dir(1:imax) = 0.3925 ! coszen=0 value of above expression
        end where
        where ( .not.land(istart:iend) )
          cuvrf_dif(1:imax) = 0.06
          cirrf_dir(1:imax) = cuvrf_dir(1:imax)
          cirrf_dif(1:imax) = 0.06
          cuvrf_dir(1:imax) = 0.85*fracice(istart:iend) + (1.-fracice(istart:iend))*cuvrf_dir(1:imax)
          cuvrf_dif(1:imax) = 0.85*fracice(istart:iend) + (1.-fracice(istart:iend))*cuvrf_dif(1:imax)
          cirrf_dir(1:imax) = 0.45*fracice(istart:iend) + (1.-fracice(istart:iend))*cirrf_dir(1:imax)
          cirrf_dif(1:imax) = 0.45*fracice(istart:iend) + (1.-fracice(istart:iend))*cirrf_dif(1:imax)
        end where
      end if
    end if ! always_mspeca
    
    ! Set up ozone for this time and row
    if ( amipo3 ) then
      sigh(1:kl) = sigmh(1:kl)
      sigh(kl+1) = 0.
      call o3set_amip(rlatt(istart:iend),imax,mins,sigh,ps(istart:iend),duo3n)
    else
      ! note levels are inverted
      call o3set(imax,istart,mins,duo3n,sig,ps(istart:iend))
    end if
    Rad_gases(mythread)%qo3(:,1,:) = real( max( 1.e-10, duo3n ), 8 )

    ! Aerosols -------------------------------------------------------
    do k = 1,kl
      rhoa(:,k) = ps(istart:iend)*sig(k)/(rdry*t(istart:iend,k)) !density of air
      dz(:,k) = -rdry*dsig(k)*t(istart:iend,k)/(grav*sig(k))
      dz(:,k) = min( max(dz(:,k), 1.), 2.e4 )
    end do
    select case (abs(iaero))
      case(0)
        ! no aerosols
      case(1)
        ! aerosols are read in (direct effect only)
        do i = 1,imax
          iq = i + iq_tile - 1
          cosz = max( coszro(i), 1.e-4 )
          delta = coszro(i)*0.29*8.*so4t(iq)*((1.-0.25*(cuvrf_dir(i)+cuvrf_dif(i)+cirrf_dir(i)+cirrf_dif(i)))/cosz)**2
          cuvrf_dir(i) = min(0.99, delta+cuvrf_dir(i)) ! still broadband
          cirrf_dir(i) = min(0.99, delta+cirrf_dir(i)) ! still broadband
          cuvrf_dif(i) = min(0.99, delta+cuvrf_dif(i)) ! still broadband
          cirrf_dif(i) = min(0.99, delta+cirrf_dif(i)) ! still broadband
        end do ! i=1,imax
      case(2,3)
        ! prognostic aerosols
        ! convert to units kg / m^2
        Aerosol(mythread)%aerosol(:,:,:,:) = 0._8
        if ( so4radmethod/=-1 ) then
          do k = 1,kl
            kr = kl + 1 - k
            dzrho = rhoa(:,k)*dz(:,k)
            ! Factor of 132.14/32.06 converts from sulfur to ammmonium sulfate
            Aerosol(mythread)%aerosol(:,1,kr,1) = real((132.14/32.06)*xtg(istart:iend,k,3)*dzrho, 8) ! so4
          end do
        end if
        if ( carbonradmethod/=-1 ) then
          do k = 1,kl
            kr = kl + 1 - k
            dzrho = rhoa(:,k)*dz(:,k)
            Aerosol(mythread)%aerosol(:,1,kr,2) = real(xtg(istart:iend,k,4)*dzrho, 8)        ! bc hydrophobic
            Aerosol(mythread)%aerosol(:,1,kr,3) = real(xtg(istart:iend,k,5)*dzrho, 8)        ! bc hydrophilic
            Aerosol(mythread)%aerosol(:,1,kr,4) = real(xtg(istart:iend,k,6)*dzrho, 8)        ! oc hydrophobic
            Aerosol(mythread)%aerosol(:,1,kr,5) = real(xtg(istart:iend,k,7)*dzrho, 8)        ! oc hydrophilic
          end do
        end if
        if ( dustradmethod/=-1 ) then
          do k = 1,kl
            kr = kl + 1 - k
            dzrho = rhoa(:,k)*dz(:,k)
            Aerosol(mythread)%aerosol(:,1,kr,6) = real(xtg(istart:iend,k,8)*dzrho, 8)        ! dust 0.8
            Aerosol(mythread)%aerosol(:,1,kr,7) = real(xtg(istart:iend,k,9)*dzrho, 8)        ! dust 1.0
            Aerosol(mythread)%aerosol(:,1,kr,8) = real(xtg(istart:iend,k,10)*dzrho, 8)       ! dust 2.0
            Aerosol(mythread)%aerosol(:,1,kr,9) = real(xtg(istart:iend,k,11)*dzrho, 8)       ! dust 4.0
          end do
        end if
        if ( seasaltradmethod/=-1 ) then
          do k = 1,kl
            kr = kl + 1 - k
            !Aerosol(mythread)%aerosol(:,1,kr,10) = real(xtg(istart:iend,k,12)*dzrho         & ! Small film sea salt (0.1)
            !                                           +xtg(istart:iend,k,13)*dzrho,8)        ! Large jet sea salt (0.5)
            Aerosol(mythread)%aerosol(:,1,kr,10) = real(xtg(istart:iend,k,12)*dzrho,8)        ! Small film sea salt (0.1)
            Aerosol(mythread)%aerosol(:,1,kr,11) = real(xtg(istart:iend,k,13)*dzrho,8)        ! Large jet sea salt (0.5)
          end do
        end if
        !Aerosol(mythread)%aerosol=min(max(Aerosol(mythread)%aerosol, 0._8), 2.e-4_8)
        Aerosol(mythread)%aerosol=max(Aerosol(mythread)%aerosol, 0._8)
        
        !if ( Rad_control%using_im_bcsul ) then
        !  Aerosol_props(mythread)%ivol(:,1,:) = 100-nint(100.*Aerosol(mythread)%aerosol(:,1,:,1)/ &
        !      max(Aerosol(mythread)%aerosol(:,1,:,1)+Aerosol(mythread)%aerosol(:,1,:,2)+Aerosol(mythread)%aerosol(:,1,:,3), 1.E-30_8))
        !  Aerosol_props(mythread)%ivol(:,1,:) = max(min(Aerosol_props(mythread)%ivol(:,1,:), 100), 0)
        !else
          Aerosol_props%ivol(:,:,:) = 0 ! no mixing of bc with so4
        !end if
    end select

    ! define droplet size distribution ------------------------------
    call aerodrop(istart,cd2,rhoa)
    
    ! Cloud fraction diagnostics ------------------------------------
    cloudlo(istart:iend) = 0.
    cloudmi(istart:iend) = 0.
    cloudhi(istart:iend) = 0.
    ! Diagnose low, middle and high clouds
    if ( nmr>=1 ) then
      ! max-rnd cloud overlap
      mx = 0.
      do k = 1,nlow
        mx = max(mx, cfrac(istart:iend,k))
        where ( cfrac(istart:iend,k)<1.e-4 )
          cloudlo(istart:iend) = cloudlo(istart:iend) + mx*(1.-cloudlo(istart:iend))
          mx = 0.
        end where
      end do
      cloudlo(istart:iend) = cloudlo(istart:iend) + mx*(1.-cloudlo(istart:iend))
      mx = 0.
      do k = nlow+1,nmid
        mx = max(mx, cfrac(istart:iend,k))
        where ( cfrac(istart:iend,k)<1.e-4 )
          cloudmi(istart:iend) = cloudmi(istart:iend) + mx*(1.-cloudmi(istart:iend))
          mx = 0.
        end where
      end do
      cloudmi(istart:iend) = cloudmi(istart:iend) + mx*(1.-cloudmi(istart:iend))
      mx = 0.
      do k = nmid+1,kl-1
        mx = max(mx, cfrac(istart:iend,k))
        where ( cfrac(istart:iend,k)<1.e-4 )
          cloudhi(istart:iend) = cloudhi(istart:iend) + mx*(1.-cloudhi(istart:iend))
          mx = 0.
        end where
      end do
      cloudhi(istart:iend) = cloudhi(istart:iend) + mx*(1.-cloudhi(istart:iend))  
    else
      ! rnd cloud overlap
      do k = 1,nlow
        cloudlo(istart:iend) = cloudlo(istart:iend) + cfrac(istart:iend,k)*(1.-cloudlo(istart:iend))
      end do
      do k = nlow+1,nmid
        cloudmi(istart:iend) = cloudmi(istart:iend) + cfrac(istart:iend,k)*(1.-cloudmi(istart:iend))
      end do
      do k = nmid+1,kl-1
        cloudhi(istart:iend) = cloudhi(istart:iend) + cfrac(istart:iend,k)*(1.-cloudhi(istart:iend))
      end do
    end if

    ! Prepare SEA-ESF arrays ----------------------------------------
    dumtss = min( max( tss(istart:iend), 100.), 365.)
    do k = 1,kl
      kr = kl + 1 - k
      dumt(:,k) = min( max( t(istart:iend,k), 100.), 365.)
      p2(:,k) = ps(istart:iend)*sig(k)
      Atmos_input(mythread)%deltaz(:,1,kr)  = real(dz(:,k), 8)
      Atmos_input(mythread)%rh2o(:,1,kr)    = max(real(qg(istart:iend,k), 8), 2.E-7_8)
      Atmos_input(mythread)%temp(:,1,kr)    = real(dumt(:,k),8)    
      Atmos_input(mythread)%press(:,1,kr)   = real(p2(:,k), 8)
      Atmos_input(mythread)%rel_hum(:,1,kr) = min(real(qg(istart:iend,k)/qsat(p2(:,k),dumt(:,k),imax), 8), 1._8)
    end do
    Atmos_input(mythread)%temp(:,1,kl+1)  = real(dumtss, 8)
    Atmos_input(mythread)%press(:,1,kl+1) = real(ps(istart:iend), 8)
    Atmos_input(mythread)%pflux(:,1,1  )  = 0._8
    Atmos_input(mythread)%tflux(:,1,1  )  = Atmos_input(mythread)%temp(:,1,1)
    do k = 1,kl-1
      kr = kl + 1 - k
      Atmos_input(mythread)%pflux(:,1,kr) = real(rathb(k)*p2(:,k)+ratha(k)*p2(:,k+1), 8)
      Atmos_input(mythread)%tflux(:,1,kr) = real(rathb(k)*dumt(:,k)+ratha(k)*dumt(:,k+1), 8)
    end do
    Atmos_input(mythread)%pflux(:,1,kl+1) = real(ps(istart:iend), 8)
    Atmos_input(mythread)%tflux(:,1,kl+1) = real(dumtss, 8)
    Atmos_input(mythread)%clouddeltaz     = Atmos_input(mythread)%deltaz

    Atmos_input(mythread)%psfc(:,1)    = real(ps(istart:iend), 8)
    Atmos_input(mythread)%tsfc(:,1)    = real(dumtss, 8)
    Atmos_input(mythread)%phalf(:,1,1) = 0._8
    do k = 1,kl-1
      kr = kl + 1 - k
      Atmos_input(mythread)%phalf(:,1,kr) = real(rathb(k)*p2(:,k)+ratha(k)*p2(:,k+1), 8)
    end do
    Atmos_input(mythread)%phalf(:,1,kl+1) = real(ps(istart:iend), 8)
    if ( do_aerosol_forcing ) then
      Atmos_input(mythread)%aerosolrelhum = Atmos_input(mythread)%rel_hum
    end if
    
    ! cloud overlap
    if ( nmr>=1 ) then
      do i = 1,imax ! maximum-random overlap
        iq = i + istart - 1
        Cld_spec(mythread)%cmxolw(i,1,:) = 0._8
        k = 1
        do while ( k<kl )
          ktop = k
          if ( cfrac(iq,k)>0. ) then
            kbot = k ! found bottom of cloud
            do while ( ktop<kl .and. cfrac(iq,min(ktop+1, kl))>0. )
              ktop = ktop + 1 ! search for top of cloud
            end do
            if ( ktop>kbot ) then ! if multi-layer cloud, calculate common max overlap fraction
              Cld_spec(mythread)%cmxolw(i,1,kl+1-ktop:kl+1-kbot) = real(minval(cfrac(iq,kbot:ktop)), 8)
            end if
          end if
          k = ktop + 1
        end do
      end do  
      do k = 1,kl
        kr = kl + 1 - k
        do i = 1,imax
          iq = i + istart - 1  
          ! Max+Rnd overlap clouds for SW
          Cld_spec(mythread)%camtsw(i,1,kr) = real(cfrac(iq,k), 8)                                   
          ! Max overlap for LW
          Cld_spec(mythread)%cmxolw(i,1,kr) = min(Cld_spec(mythread)%cmxolw(i,1,kr), 0.999_8)        
          ! Rnd overlap for LW
          Cld_spec(mythread)%crndlw(i,1,kr) = real(cfrac(iq,k), 8)-Cld_spec(mythread)%cmxolw(i,1,kr) 
        end do
      end do
    else
      do k = 1,kl
        kr = kl + 1 - k  
        do i = 1,imax ! random overlap
          iq = i + istart - 1
          Cld_spec(mythread)%camtsw(i,1,kr) = real(cfrac(iq,k), 8) ! Max+Rnd overlap clouds for SW
          Cld_spec(mythread)%crndlw(i,1,kr) = real(cfrac(iq,k), 8) ! Rnd overlap for LW
          Cld_spec(mythread)%cmxolw(i,1,kr) = 0._8
        end do
      end do
    end if

    ! cloud microphysics for radiation
    ! cfrac, qlrad and qfrad also include convective cloud as well as qfg and qlg
    dumcf = cfrac(istart:iend,:)
    dumql = qlrad(istart:iend,:)
    dumqf = qfrad(istart:iend,:)
    call cloud3(Cloud_microphysics(mythread)%size_drop,Cloud_microphysics(mythread)%size_ice,       &
                Cloud_microphysics(mythread)%conc_drop,Cloud_microphysics(mythread)%conc_ice,       &
                dumcf,dumql,dumqf,p2,dumt,cd2,imax,kl)
    Cloud_microphysics(mythread)%size_drop = max(Cloud_microphysics(mythread)%size_drop, 1.e-20_8)
    Cloud_microphysics(mythread)%size_ice  = max(Cloud_microphysics(mythread)%size_ice,  1.e-20_8)                
    Cloud_microphysics(mythread)%size_rain = 1.e-20_8
    Cloud_microphysics(mythread)%conc_rain = 0._8
    Cloud_microphysics(mythread)%size_snow = 1.e-20_8
    Cloud_microphysics(mythread)%conc_snow = 0._8
    
    Lscrad_props(mythread)%cldext   = 0._8
    Lscrad_props(mythread)%cldsct   = 0._8
    Lscrad_props(mythread)%cldasymm = 1._8
    Lscrad_props(mythread)%abscoeff = 0._8
    call microphys_lw_driver(1, imax, 1, 1, Cloud_microphysics(mythread), Micro_rad_props=Lscrad_props(mythread))
    call microphys_sw_driver(1, imax, 1, 1, Cloud_microphysics(mythread), Micro_rad_props=Lscrad_props(mythread))
    Cldrad_props(mythread)%cldsct(:,:,:,:,1)   = Lscrad_props(mythread)%cldsct(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props(mythread)%cldext(:,:,:,:,1)   = Lscrad_props(mythread)%cldext(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props(mythread)%cldasymm(:,:,:,:,1) = Lscrad_props(mythread)%cldasymm(:,:,:,:) ! Large scale cloud properties only
    Cldrad_props(mythread)%abscoeff(:,:,:,:,1) = Lscrad_props(mythread)%abscoeff(:,:,:,:) ! Large scale cloud properties only
    
    call lwemiss_calc(Atmos_input(mythread)%clouddeltaz,Cldrad_props(mythread)%abscoeff, &
                      Cldrad_props(mythread)%cldemiss)
    Cldrad_props(mythread)%emmxolw = Cldrad_props(mythread)%cldemiss
    Cldrad_props(mythread)%emrndlw = Cldrad_props(mythread)%cldemiss

    Surface(mythread)%asfc_vis_dir(:,1) = real(cuvrf_dir(:), 8)
    Surface(mythread)%asfc_nir_dir(:,1) = real(cirrf_dir(:), 8)
    Surface(mythread)%asfc_vis_dif(:,1) = real(cuvrf_dif(:), 8)
    Surface(mythread)%asfc_nir_dif(:,1) = real(cirrf_dif(:), 8)
   
    Astro(mythread)%cosz(:,1)    = max(real(coszro, 8), 0._8)
    Astro(mythread)%fracday(:,1) = real(taudar, 8)

    
    call longwave_driver (1, imax, 1, 1, Rad_time, Atmos_input(mythread),         &
                          Rad_gases(mythread), Aerosol(mythread), Aerosol_props,  &
                          Cldrad_props(mythread), Cld_spec(mythread),             &
                          Aerosol_diags(mythread), Lw_output(mythread),           &
                          Lw_diagnostics(mythread), Lw_clouds(mythread),          &
                          Optical(mythread), Gas_tf(mythread))

    
    call shortwave_driver (1, imax, 1, 1, Atmos_input(mythread), Surface(mythread),  &
                           Astro(mythread), Aerosol(mythread), Aerosol_props,        &
                           Rad_gases(mythread), Cldrad_props(mythread),              &
                           Cld_spec(mythread), Sw_output(mythread:mythread),         &
                           Aerosol_diags(mythread), r)

    
    ! store shortwave and fbeam data --------------------------------
    sgdn(istart:iend) = real(Sw_output(mythread)%dfsw(:,1,kl+1,1))
    sgdnvis           = real(Sw_output(mythread)%dfsw_vis_sfc(:,1,1))
    sgdnnir           = sgdn(istart:iend) - sgdnvis
    ! sg = +Snet = Sdown - Sup
    sgn        = sgdn(istart:iend) - real(Sw_output(mythread)%ufsw(:,1,kl+1,1))
    sgvis      = sgdnvis - real(Sw_output(mythread)%ufsw_vis_sfc(:,1,1))
    !sgvisdir  = Sw_output(mythread)%dfsw_vis_sfc_dir(:,1,1)-Sw_output(mythread)%ufsw_vis_sfc_dir(:,1,1)
    !sgvisdif  = Sw_output(mythread)%dfsw_vis_sfc_dif(:,1,1)-Sw_output(mythread)%ufsw_vis_sfc_dif(:,1,1)
    !sgnirdir  = Sw_output(mythread)%dfsw_dir_sfc(:,1,1)-sgvisdir
    !sgnirdif  = Sw_output(mythread)%dfsw_dif_sfc(:,1,1)-Sw_output(mythread)%ufsw_dif_sfc(:,1,1)-sgvisdif
    !sgdir     = Sw_output(mythread)%dfsw_dir_sfc(:,1,1)
    !sgdif     = Sw_output(mythread)%dfsw_dif_sfc(:,1,1)-Sw_output(mythread)%ufsw_dif_sfc(:,1,1)
    sgdnvisdir = real(Sw_output(mythread)%dfsw_vis_sfc_dir(:,1,1))
    sgdnvisdif = real(Sw_output(mythread)%dfsw_vis_sfc_dif(:,1,1))
    sgdnnirdir = real(Sw_output(mythread)%dfsw_dir_sfc(:,1,1)) - sgdnvisdir
    sgdnnirdif = real(Sw_output(mythread)%dfsw_dif_sfc(:,1,1)) - sgdnvisdif
    where ( sgdn(istart:iend)<0.001 )
      swrsave(istart:iend)  = 0.5
    elsewhere  
      swrsave(istart:iend)  = sgdnvis/sgdn(istart:iend)
    end where
    where ( sgdnvis<0.001 )
      fbeamvis(istart:iend) = 0.
    elsewhere  
      fbeamvis(istart:iend) = sgdnvisdir/sgdnvis
    end where
    where ( sgdnnir<0.001 )
      fbeamnir(istart:iend) = 0.
    elsewhere  
      fbeamnir(istart:iend) = sgdnnirdir/sgdnnir
    end where
    where ( coszro<=1.e-5 )
      dni(istart:iend) = 0.
    elsewhere
      dni(istart:iend) = (sgdnvisdir+sgdnnirdir)/coszro
    end where
        
    ! longwave output -----------------------------------------------
    rgn(istart:iend) = real(Lw_output(mythread)%flxnet(:,1,kl+1))          ! longwave at surface
    rt(istart:iend) = real(Lw_output(mythread)%flxnet(:,1,1))              ! longwave at top
    ! rg = -Rnet = Rup - Rdown = stefbo T^4 - Rdown
    rgdn(istart:iend) = stefbo*tss(istart:iend)**4 - rgn(istart:iend)

    ! shortwave output ----------------------------------------------
    sint(istart:iend) = real(Sw_output(mythread)%dfsw(:,1,1,1))   ! solar in top
    sout(istart:iend) = real(Sw_output(mythread)%ufsw(:,1,1,1))   ! solar out top

    ! Clear sky calculation -----------------------------------------
    if ( do_totcld_forcing ) then
      soutclr(istart:iend) = real(Sw_output(mythread)%ufswcf(:,1,1,1))      ! solar out top
      sgclr(istart:iend)   = -real(Sw_output(mythread)%fswcf(:,1,kl+1,1))   ! solar absorbed at the surface
      rtclr(istart:iend)   = real(Lw_output(mythread)%flxnetcf(:,1,1))      ! clr sky lw at top
      rgclr(istart:iend)   = real(Lw_output(mythread)%flxnetcf(:,1,kl+1))   ! clear sky longwave at surface
    else
      soutclr(istart:iend) = 0.
      sgclr(istart:iend)   = 0.
      rtclr(istart:iend)   = 0.
      rgclr(istart:iend)   = 0.
    end if

    ! heating rate --------------------------------------------------
    do k = 1,kl
      ! total heating rate (convert deg K/day to deg K/sec)
      sw_tend(istart:iend,kl+1-k) = -real(Sw_output(mythread)%hsw(:,1,k,1)/86400._8)
      lw_tend(istart:iend,kl+1-k) = -real(Lw_output(mythread)%heatra(:,1,k)/86400._8)
    end do
    
#ifdef seaesfdebug
    if ( any(Sw_output(mythread)%hsw(:,1,:,1)/=Sw_output(mythread)%hsw(:,1,:,1)) ) then
      write(6,*) "ERROR: NaN detected in hsw for seaesfrad on myid=",myid
      call ccmpi_abort(-1)
    end if
    if ( any(Sw_output(mythread)%hsw(:,1,:,1)>8640000.) .or. any(Sw_output(mythread)%hsw(:,1,:,1)<-8640000.) ) then
      write(6,*) "ERROR: hsw is out-of-range in seaesfrad on myid=",myid  
      write(6,*) "minval,maxval ",minval(Sw_output(mythread)%hsw(:,1,:,1)),maxval(Sw_output(mythread)%hsw(:,1,:,1))
      write(6,*) "minloc,maxloc ",minloc(Sw_output(mythread)%hsw(:,1,:,1)),maxloc(Sw_output(mythread)%hsw(:,1,:,1))
      call ccmpi_abort(-1)
    end if
    if ( any(Lw_output(mythread)%heatra(:,1,:)/=Lw_output(mythread)%heatra(:,1,:)) ) then
      write(6,*) "ERROR: NaN detected in heatra for seaesfrad on myid=",myid
      call ccmpi_abort(-1)
    end if
    if ( any(Lw_output(mythread)%heatra(:,1,:)>8640000.) .or. any(Lw_output(mythread)%heatra(:,1,:)<-8640000.) ) then
      write(6,*) "ERROR: heatra is out-of-range in seaesfrad on myid=",myid  
      write(6,*) "minval,maxval ",minval(Lw_output(mythread)%heatra(:,1,:)),maxval(Lw_output(mythread)%heatra(:,1,:))
      write(6,*) "minloc,maxloc ",minloc(Lw_output(mythread)%heatra(:,1,:)),maxloc(Lw_output(mythread)%heatra(:,1,:))
      call ccmpi_abort(-1)
    end if
#endif
    
    ! aerosol optical depths ----------------------------------------
    if ( do_aerosol_forcing ) then
      opticaldepth(istart:iend,:,:) = 0.
      ! Sulfate
      do k = 1,kl
        opticaldepth(istart:iend,3,1) = opticaldepth(istart:iend,3,1) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,1,1)) ! Visible
        opticaldepth(istart:iend,3,2) = opticaldepth(istart:iend,3,2) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,1,2)) ! Near IR
        opticaldepth(istart:iend,3,3) = opticaldepth(istart:iend,3,3) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,1,3)) ! Longwave
      end do
      ! BC
      do nr = 2,3
        do k = 1,kl          
          opticaldepth(istart:iend,5,1) = opticaldepth(istart:iend,5,1) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,5,2) = opticaldepth(istart:iend,5,2) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,5,3) = opticaldepth(istart:iend,5,3) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! OC
      do nr = 4,5
        do k = 1,kl    
          opticaldepth(istart:iend,6,1) = opticaldepth(istart:iend,6,1) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,6,2) = opticaldepth(istart:iend,6,2) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,6,3) = opticaldepth(istart:iend,6,3) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! Small dust
      do k = 1,kl
        opticaldepth(istart:iend,1,1) = opticaldepth(istart:iend,1,1) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,6,1)) ! Visible
        opticaldepth(istart:iend,1,2) = opticaldepth(istart:iend,1,2) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,6,2)) ! Near IR
        opticaldepth(istart:iend,1,3) = opticaldepth(istart:iend,1,3) &
                                      + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,6,3)) ! Longwave
      end do
      ! Large dust
      do nr = 7,9
        do k = 1,kl
          opticaldepth(istart:iend,2,1) = opticaldepth(istart:iend,2,1) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,2,2) = opticaldepth(istart:iend,2,2) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,2,3) = opticaldepth(istart:iend,2,3) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
      ! Seasalt
      do nr = 10,11
        do k = 1,kl
          opticaldepth(istart:iend,7,1) = opticaldepth(istart:iend,7,1) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,7,2) = opticaldepth(istart:iend,7,2) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,7,3) = opticaldepth(istart:iend,7,3) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do  
      ! Aerosol
      do nr = 1,nfields
        do k = 1,kl
          opticaldepth(istart:iend,4,1) = opticaldepth(istart:iend,4,1) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,1)) ! Visible
          opticaldepth(istart:iend,4,2) = opticaldepth(istart:iend,4,2) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,2)) ! Near IR
          opticaldepth(istart:iend,4,3) = opticaldepth(istart:iend,4,3) &
                                        + real(Aerosol_diags(mythread)%extopdep(1:imax,1,k,nr,3)) ! Longwave
        end do
      end do
    end if

    if ( always_mspeca ) then
      where ( coszro*taudar<=1.E-5 )
        sgdn_amp(istart:iend) = 0.
        sgn_amp(istart:iend)  = 0.
        dni_amp(istart:iend)  = 0.
        sint_amp(istart:iend) = 0.
        sout_amp(istart:iend) = 0.
        soutclr_amp(istart:iend) = 0.
        sgclr_amp(istart:iend)   = 0.
      elsewhere
        sgdn_amp(istart:iend) = sgdn(istart:iend)
        sgn_amp(istart:iend)  = sgn/(coszro*taudar)
        dni_amp(istart:iend)  = dni(istart:iend)/taudar
        sint_amp(istart:iend) = sint(istart:iend)
        sout_amp(istart:iend) = sout(istart:iend)
        soutclr_amp(istart:iend) = soutclr(istart:iend)
        sgclr_amp(istart:iend)   = sgclr(istart:iend)
      end where
      do k = 1,kl
        where ( coszro*taudar<=1.E-5 )
          sw_tend_amp(istart:iend,k)  = 0.
        elsewhere
          sw_tend_amp(istart:iend,k) = sw_tend(istart:iend,k)
        end where   
      end do
      sgsave(istart:iend) = sgn   ! Net solar radiation (after solar fit)
      if ( ktau>0 ) then
        if ( istart==1 ) koundiag = koundiag + 1  
        cloudtot(istart:iend) = 1. - (1.-cloudlo(istart:iend))*(1.-cloudmi(istart:iend))*(1.-cloudhi(istart:iend))
        sint_ave(istart:iend) = sint_ave(istart:iend) + sint(istart:iend)
        sot_ave(istart:iend)  = sot_ave(istart:iend)  + sout(istart:iend)
        soc_ave(istart:iend)  = soc_ave(istart:iend)  + soutclr(istart:iend)
        rtu_ave(istart:iend)  = rtu_ave(istart:iend)  + rt(istart:iend)
        rtc_ave(istart:iend)  = rtc_ave(istart:iend)  + rtclr(istart:iend)
        rgn_ave(istart:iend)  = rgn_ave(istart:iend)  + rgn(istart:iend)
        rgc_ave(istart:iend)  = rgc_ave(istart:iend)  + rgclr(istart:iend)
        rgdn_ave(istart:iend) = rgdn_ave(istart:iend) + rgdn(istart:iend)
        sgc_ave(istart:iend)  = sgc_ave(istart:iend)  + sgclr(istart:iend)
        cld_ave(istart:iend)  = cld_ave(istart:iend)  + cloudtot(istart:iend)
        cll_ave(istart:iend)  = cll_ave(istart:iend)  + cloudlo(istart:iend)
        clm_ave(istart:iend)  = clm_ave(istart:iend)  + cloudmi(istart:iend)
        clh_ave(istart:iend)  = clh_ave(istart:iend)  + cloudhi(istart:iend)
      end if
    else
      ! Calculate the amplitude of the diurnal cycle of solar radiation
      ! at the surface (using the value for the middle of the radiation
      ! step) and use this value to get solar radiation at other times.
      ! Use the zenith angle and daylight fraction calculated in zenith
      ! to remove these factors.
      where ( coszro*taudar<=1.E-5 )
        ! The sun is not up at all over the radiation period so no 
        ! fitting need be done.
        sgdn_amp(istart:iend) = 0.
        sgn_amp(istart:iend)  = 0.
        dni_amp(istart:iend)  = 0.
        sint_amp(istart:iend) = 0.
        sout_amp(istart:iend) = 0.
        soutclr_amp(istart:iend) = 0.
        sgclr_amp(istart:iend)   = 0.
      elsewhere
        sgdn_amp(istart:iend) = sgdn(istart:iend)/(coszro*taudar)
        sgn_amp(istart:iend)  = sgn/(coszro*taudar)
        dni_amp(istart:iend)  = dni(istart:iend)/taudar
        sint_amp(istart:iend) = sint(istart:iend)/(coszro*taudar)
        sout_amp(istart:iend) = sout(istart:iend)/(coszro*taudar)
        soutclr_amp(istart:iend) = soutclr(istart:iend)/(coszro*taudar)
        sgclr_amp(istart:iend)   = sgclr(istart:iend)/(coszro*taudar)
      end where
      do k = 1,kl
        where ( coszro*taudar<=1.E-5 )
          sw_tend_amp(istart:iend,k)  = 0.
        elsewhere
          sw_tend_amp(istart:iend,k) = sw_tend(istart:iend,k)/(coszro*taudar)
        end where   
      end do
    end if  

    ! Save things for non-radiation time steps ----------------------
    ! Save the value excluding Ts^4 part.  This is allowed to change.
    rgsave(istart:iend) = rgn(istart:iend) - stefbo*tss(istart:iend)**4
    slwa(istart:iend) = -sgsave(istart:iend) + rgsave(istart:iend)
   
  end if  ! odcalc

  ! Store albedo data ---------------------------------------------
  albvisnir(istart:iend,1) = cuvrf_dir*fbeamvis(istart:iend)      &
                           + cuvrf_dif*(1.-fbeamvis(istart:iend))
  albvisnir(istart:iend,2) = cirrf_dir*fbeamnir(istart:iend)      &
                           + cirrf_dif*(1.-fbeamnir(istart:iend))
  
  ! Store fraction of direct radiation in urban scheme
  fbeam = fbeamvis(istart:iend)*swrsave(istart:iend) &
        + fbeamnir(istart:iend)*(1.-swrsave(istart:iend))
  call atebfbeam(istart,imax,fbeam,0)
  
  ! cloud amounts for saving --------------------------------------
  cloudtot(istart:iend) = 1. - (1.-cloudlo(istart:iend))*(1.-cloudmi(istart:iend))*(1.-cloudhi(istart:iend))
   
  ! Calculate the solar using the saved amplitude.
  if ( always_mspeca ) then
    sgdn(istart:iend) = sgdn_amp(istart:iend)  
    sgn = sgn_amp(istart:iend)*coszro2*taudar2
    dni(istart:iend) = dni_amp(istart:iend)*taudar2
    sint(istart:iend) = sint_amp(istart:iend)
    sout(istart:iend) = sout_amp(istart:iend)
    soutclr(istart:iend) = soutclr_amp(istart:iend)
    sgclr(istart:iend)   = sgclr_amp(istart:iend)
    sgsave(istart:iend) = sgn   ! Net solar radiation (after solar fit)
    slwa(istart:iend) = -sgsave(istart:iend) + rgsave(istart:iend)
    do k = 1,kl
      sw_tend(istart:iend,k) = sw_tend_amp(istart:iend,k)
    end do
  else
    sgdn(istart:iend) = sgdn_amp(istart:iend)*coszro2*taudar2
    alb = swrsave(istart:iend)*albvisnir(istart:iend,1)     &
       + (1.-swrsave(istart:iend))*albvisnir(istart:iend,2)
    sgn = sgdn(istart:iend)*(1.-alb)
    sgn_save = sgn_amp(istart:iend)*coszro2*taudar2
    dni(istart:iend)  = dni_amp(istart:iend)*taudar2
    sint(istart:iend) = sint_amp(istart:iend)*coszro2*taudar2
    sout(istart:iend) = sout_amp(istart:iend)*coszro2*taudar2
    sout(istart:iend) = sout(istart:iend) + sgn_save - sgn ! correct for changing albedo
    soutclr(istart:iend) = soutclr_amp(istart:iend)*coszro2*taudar2
    sgclr(istart:iend)   = sgclr_amp(istart:iend)*coszro2*taudar2
    ! Set up the CC model radiation fields
    ! slwa is negative net radiational htg at ground
    ! Note that this does not include the upward LW radiation from the surface.
    ! That is included in sflux.f90
    sgsave(istart:iend) = sgn   ! Net solar radiation (after solar fit)
    slwa(istart:iend) = -sgsave(istart:iend) + rgsave(istart:iend)
    ! Update tendencies
    do k = 1,kl
      sw_tend(istart:iend,k) = sw_tend_amp(istart:iend,k)*coszro2(1:imax)*taudar2(1:imax)
    end do
  end if  

end do  ! iq_tile = 1,ifull,imax
!$omp end do nowait


if ( diag .and. mydiag ) then
  diag_temp(:) = t(idjd,:)
  write(6,*) "tdiag ",diag_temp
  diag_temp(:) = qg(idjd,:)
  write(6,*) "qgdiag ",diag_temp
end if

return
end subroutine seaesfrad

subroutine longwave_driver (is, ie, js, je, Rad_time, Atmos_input, &
                            Rad_gases, Aerosol, Aerosol_props,     &
                            Cldrad_props, Cld_spec, Aerosol_diags, &
                            Lw_output, Lw_diagnostics, Lw_clouds,  &
                            Optical, Gas_tf)

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
type(lw_output_type), dimension(1), intent(inout)  :: Lw_output
type(lw_diagnostics_type),          intent(inout)  :: Lw_diagnostics
type(lw_clouds_type),               intent(inout)  :: Lw_clouds
type(optical_path_type),            intent(inout)  :: Optical
type(gas_tf_type),                  intent(inout)  :: Gas_tf

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
 
LW_output(1)%bdy_flx(:,:,:) = 0._8

call sealw99 (is, ie, js, je, Rad_time, Atmos_input,           &
              Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
              Cld_spec, Aerosol_diags, Lw_output(1),           &
              Lw_diagnostics, do_aerosol_forcing, Lw_clouds,   &
              Optical, Gas_tf)

return
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

integer,                         intent(in)       :: is, ie, js, je
type(atmos_input_type),          intent(in)       :: Atmos_input     
type(surface_type),              intent(in)       :: Surface     
type(astronomy_type),            intent(in)       :: Astro           
type(radiative_gases_type),      intent(in)       :: Rad_gases   
type(aerosol_type),              intent(in)       :: Aerosol     
type(aerosol_properties_type),   intent(inout)    :: Aerosol_props
type(cldrad_properties_type),    intent(in)       :: Cldrad_props
type(cld_specification_type),    intent(in)       :: Cld_spec
type(sw_output_type), dimension(1), intent(inout) :: Sw_output
type(aerosol_diagnostics_type), intent(inout)     :: Aerosol_diags
real(kind=8), dimension(:,:,:,:), intent(inout)   :: r
integer :: naerosol_optical

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

Sw_output(1)%fsw (:,:,:,:)    = 0.0_8
Sw_output(1)%dfsw(:,:,:,:)    = 0.0_8
Sw_output(1)%ufsw(:,:,:,:)    = 0.0_8
Sw_output(1)%hsw (:,:,:,:)    = 0.0_8
Sw_output(1)%dfsw_dir_sfc     = 0.0_8
Sw_output(1)%dfsw_dif_sfc     = 0.0_8
Sw_output(1)%ufsw_dir_sfc     = 0.0_8
Sw_output(1)%ufsw_dif_sfc     = 0.0_8
Sw_output(1)%dfsw_vis_sfc     = 0._8
Sw_output(1)%ufsw_vis_sfc     = 0._8
Sw_output(1)%dfsw_vis_sfc_dir = 0._8
Sw_output(1)%dfsw_vis_sfc_dif = 0._8
SW_output(1)%dfsw_vis_sfc_clr = 0._8
Sw_output(1)%ufsw_vis_sfc_dir = 0._8
Sw_output(1)%ufsw_vis_sfc_dif = 0._8
Sw_output(1)%bdy_flx(:,:,:,:) = 0.0_8       

if (Rad_control%do_totcld_forcing) then
  Sw_output(1)%fswcf (:,:,:,:) = 0.0_8
  Sw_output(1)%dfswcf(:,:,:,:) = 0.0_8
  Sw_output(1)%ufswcf(:,:,:,:) = 0.0_8
  Sw_output(1)%hswcf (:,:,:,:) = 0.0_8
  Sw_output(1)%dfsw_dir_sfc_clr = 0.0_8
  Sw_output(1)%dfsw_dif_sfc_clr  = 0.0_8
  Sw_output(1)%bdy_flx_clr (:,:,:,:) = 0.0_8
end if

if (do_aerosol_forcing) then
  naerosol_optical = size(Aerosol_props%aerextband,2)
else
  naerosol_optical = 0  
end if 

call swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,    &
             Aerosol, Aerosol_props, Astro, Cldrad_props,        &
             Cld_spec, include_volcanoes,                        &
             Sw_output(1), Aerosol_diags, r,                     &
             do_aerosol_forcing, naerosol_optical)

return
end subroutine shortwave_driver

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlrad,qfrad,prf,ttg,cdrop,imax,kl)

use cc_mpi              ! CC MPI routines
use const_phys          ! Physical constants
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: imax, kl
integer iq, k, kr
real, dimension(imax,kl), intent(in) :: cfrac, qlrad, qfrad, prf, ttg
real, dimension(imax,kl), intent(in) :: cdrop
real(kind=8), dimension(imax,kl), intent(out) :: Rdrop, Rice, conl, coni
real, dimension(imax,kl) :: reffl, reffi, Wliq, rhoa
real, dimension(imax,kl) :: eps, rk, Wice, basesize, xwgt
real :: liqparm
! parameters
logical :: do_brenguier   ! Adjust effective radius for vertically stratified cloud
real :: scale_factor      ! account for the plane-parallel homogenous cloud bias  (e.g. Cahalan effect)

scale_factor = 1.
liqparm = -0.003e-6
do_brenguier = .false.

rhoa(:,:) = prf(:,:)/(rdry*ttg(:,:))

select case(liqradmethod)
  case(0)
    liqparm = -0.003e-6 ! mid range
    do_brenguier = .false.
    scale_factor = 1.
  case(1)
    liqparm = -0.001e-6 ! lower bound
    do_brenguier = .false.
    scale_factor = 1.
  case(2)
    liqparm = -0.008e-6 ! upper bound
    do_brenguier = .false.
    scale_factor = 1.
  case(3)
    liqparm = -0.003e-6 ! mid range
    do_brenguier = .true.
    scale_factor = 0.85
  case(4)
    liqparm = -0.001e-6 ! lower bound
    do_brenguier = .true.
    scale_factor = 0.85
  case(5)
    liqparm = -0.008e-6 ! upper bound
    do_brenguier = .true.
    scale_factor = 0.85
  case default
    write(6,*) "ERROR: Unknown option for liqradmethod ",liqradmethod
    call ccmpi_abort(-1)
end select

! Reffl is the effective radius calculated following
! Martin etal 1994, JAS 51, 1823-1842
where ( qlrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
  Wliq(:,:) = rhoa(:,:)*qlrad(:,:)/cfrac(:,:) !kg/m^3
  ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
  eps(:,:) = 1. - 0.7*exp(liqparm*cdrop(:,:))
  rk(:,:)  = (1.+eps(:,:)**2)/(1.+2.*eps(:,:)**2)**2

  ! k_ratio = rk**(-1./3.)  
  ! GFDL        k_ratio (land) 1.143 (water) 1.077
  ! mid range   k_ratio (land) 1.393 (water) 1.203
  ! lower bound k_ratio (land) 1.203 (water) 1.050

  ! Martin et al 1994
  reffl(:,:) = (3.*Wliq(:,:)/(4.*pi*rhow*rk(:,:)*cdrop(:,:)))**(1./3.)
elsewhere
  reffl(:,:) = 0.
  Wliq(:,:) = 0.
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
if ( do_brenguier ) then
  if ( nmr>=1 ) then
    ! Max-Rnd overlap
    where ( cfrac(:,2)<1.e-4 )
      reffl(:,1) = reffl(:,1)*1.134
    end where
    do k = 2,kl-1
      where ( cfrac(:,k-1)<1.e-4 .and. cfrac(:,k+1)<1.e-4 )
        reffl(:,k) = reffl(:,k)*1.134
      end where
    end do
    where ( cfrac(:,kl-1)<1.e-4 )
      reffl(:,kl) = reffl(:,kl)*1.134
    end where 
  else
    ! Rnd overlap
    reffl(:,:) = reffl(:,:)*1.134
  end if
end if

select case(iceradmethod)
  case(0)
    !Lohmann et al.(1999)
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = 0.5*min(150.e-6, 3.73e-4*Wice(:,:)**0.216) 
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where
    
  case(1)
    !Donner et al (1997)
    ! linear interpolation by MJT
    where ( ttg(:,:)>250.66 )
      basesize(:,:) = 100.6
    elsewhere ( ttg(:,:)>245.66 )
      xwgt(:,:) = (ttg(:,:)-245.66)/5.
      basesize(:,:) = 80.8*(1.-xwgt(:,:)) + xwgt(:,:)*100.6
    elsewhere ( ttg(:,:)>240.66 )
      xwgt(:,:) = (ttg(:,:)-240.66)/5.
      basesize(:,:) = 93.5*(1.-xwgt(:,:)) + xwgt(:,:)*80.6
    elsewhere ( ttg(:,:)>235.66 )
      xwgt(:,:) = (ttg(:,:)-235.66)/5.
      basesize(:,:) = 63.6*(1.-xwgt(:,:)) + xwgt(:,:)*93.6
    elsewhere ( ttg(:,:)>230.66 )
      xwgt(:,:) = (ttg(:,:)-230.66)/5.
      basesize(:,:) = 42.5*(1.-xwgt(:,:)) + xwgt(:,:)*63.6
    elsewhere ( ttg(:,:)>225.66 )
      xwgt(:,:) = (ttg(:,:)-225.66)/5.
      basesize(:,:) = 39.9*(1.-xwgt(:,:)) + xwgt(:,:)*42.5
    elsewhere ( ttg(:,:)>220.66 )
      xwgt(:,:) = (ttg(:,:)-220.66)/5.
      basesize(:,:) = 21.6*(1.-xwgt(:,:)) + xwgt(:,:)*39.9
    elsewhere ( ttg(:,:)>215.66 )
      xwgt(:,:) = (ttg(:,:)-215.66)/5.
      basesize(:,:) = 20.2*(1.-xwgt(:,:)) + xwgt(:,:)*21.6
    elsewhere
      basesize(:,:) = 20.2
    end where
    where ( qfrad(:,:)>1.e-8 .and. cfrac(:,:)>1.e-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) ! kg/m**3
      reffi(:,:) = 5.e-7*basesize(:,:)
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where
   
  case(2)
    ! Fu 2007
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = 5.E-7*(47.05+0.6624*(ttg(:,:)-273.16)+0.001741*(ttg(:,:)-273.16)**2)
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where

  case(3)
    do k = 1,kl
      do iq = 1,imax
        if ( qfrad(iq,k)>1.E-8 .and. cfrac(iq,k)>1.E-4 ) then
          Wice(iq,k) = rhoa(iq,k)*qfrad(iq,k)/cfrac(iq,k) ! kg/m**3
          if ( ttg(iq,k)>248.16 ) then
            reffi(iq,k) = 5.E-7*100.6
          elseif ( ttg(iq,k)>243.16 ) then
            reffi(iq,k) = 5.E-7*80.8
          elseif ( ttg(iq,k)>238.16 ) then
            reffi(iq,k) = 5.E-7*93.5
          elseif ( ttg(iq,k)>233.16 ) then
            reffi(iq,k) = 5.E-7*63.9
          elseif ( ttg(iq,k)>228.16 ) then
            reffi(iq,k) = 5.E-7*42.5
          elseif ( ttg(iq,k)>223.16 ) then
            reffi(iq,k) = 5.E-7*39.9
          elseif ( ttg(iq,k)>218.16 ) then
            reffi(iq,k) = 5.E-7*21.6
          else
            reffi(iq,k) = 5.E-7*20.2
          end if
        else
          reffi(iq,k) = 0.
          Wice(iq,k) = 0.
        end if
      end do
    end do

end select 
  
do k = 1,kl
  kr = kl + 1 - k
  Rdrop(:,kr) = real(2.E6*reffl(:,k), 8) ! convert to diameter and microns
  Rice(:,kr)  = real(2.E6*reffi(:,k), 8)
  conl(:,kr)  = real(1000.*scale_factor*Wliq(:,k), 8) !g/m^3
  coni(:,kr)  = real(1000.*scale_factor*Wice(:,k), 8)
end do

Rdrop(:,:) = min(max(Rdrop(:,:), 8.4_8), 33.2_8) ! constrain diameter to acceptable range (see microphys_rad.f90)
Rice(:,:) = min(max(Rice(:,:), 18.6_8), 130.2_8)

return
end subroutine cloud3

subroutine calc_snow_albedo(coszro,cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif,iq_tile)

use cc_omp                                          ! CC OpenMP routines
use nsibd_m                                         ! Land-surface arrays
use parm_m                                          ! Model configuration
use soil_m                                          ! Soil and surface data
use soilsnow_m                                      ! Soil, snow and surface data

implicit none

integer, intent(in) :: iq_tile
integer i, iq
real, dimension(imax) :: coszro
real, dimension(imax), intent(inout) :: cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif
real alvo, aliro, dnsnow, ttbg, ar1, exp_ar1, ar2, exp_ar2, snr, ar3
real snrat, dtau, fage, cczen, fzen, fzenm, alvd, alv, alird, alir

! The following snow calculation should be done by sib3 (sflux.f)
alvo = 0.95         !alb. for vis. on a new snow
aliro = 0.65        !alb. for near-infr. on a new snow      
do i = 1,imax
  iq = i + iq_tile - 1
  if ( land(iq) .and. snowd(iq)>0. ) then
    dnsnow = min(1., .1*max(0.,snowd(iq)-osnowd(iq)))
    ttbg = real(isflag(iq))*tggsn(iq,1) + real(1-isflag(iq))*tgg(iq,1)
    ttbg = min(ttbg, 273.1)
    ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
    exp_ar1 = exp(ar1)                ! e.g. exp(0 to -4)
    ar2 = 10.*ar1                     ! freezing of melt water
    exp_ar2 = exp(ar2)                ! e.g. exp(0 to -40)
    snr = snowd(iq)/max(ssdnn(iq), 100.)
    if ( isoilm(iq)==9 ) then   ! fixes for Arctic & Antarctic
      ar3 = .001
      dnsnow = max(dnsnow, .0015)
      snrat = min(1., snr/(snr+.001))
    else
      ar3 = .3             ! accumulation of dirt
      snrat = min(1., snr/(snr+.02))
    end if
    dtau = 1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
    if ( snowd(iq)<=1. )then
      snage(iq) = 0.
    else
      snage(iq) = max(0., (snage(iq)+dtau)*(1.-dnsnow))
    end if
    fage = 1. - 1./(1.+snage(iq))  !age factor
    cczen = max(.17365, coszro(i))
    fzen = (1.+1./2.)/(1.+2.*2.*cczen) - 1./2.
    if ( cczen>0.5 ) fzen = 0.
    fzenm = max( fzen, 0. )
    alvd = alvo * (1.-0.2*fage)
    alv = .4 * fzenm * (1.-alvd) + alvd
    alird = aliro*(1.-.5*fage)
    alir = .4 * fzenm * (1.-alird) + alird
    cuvrf_dir(i) = (1.-snrat)*cuvrf_dir(i) + snrat*alv
    cirrf_dir(i) = (1.-snrat)*cirrf_dir(i) + snrat*alir
    cuvrf_dif(i) = cuvrf_dir(i) ! assume DIR and DIF are the same
    cirrf_dif(i) = cirrf_dir(i) ! assume DIR and DIF are the same
  end if
end do

return
end subroutine calc_snow_albedo
    
subroutine seaesfrad_init

use cc_mpi
use cc_omp
use co2_read_m
use extraout_m
use infile
use newmpar_m
use ozoneread
use parm_m
use radisw_m, only : rrco2,rrvco2,rrvch4,rrvn2o, &
    rrvf11,rrvf12,rrvf113,rrvf22
use sigs_m
use zenith_m

implicit none

include 'kuocom.h'

integer k, kr, mythread
integer jyear, jmonth, jday, jhour, jmin, mins
real f1, f2, r1, fjd, dlt, alp, slag

call START_LOG(radinit_begin)

if ( myid==0 ) write(6,*) "Initalise SEA-ESF radiation"

! Aerosol flag
do_aerosol_forcing = abs(iaero)>=2

! Set up number of minutes from beginning of year
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
fjd = float(mod(mins, 525600))/1440. ! restrict to 365 day calendar
! Calculate sun position
call solargh(fjd,bpyear,r1,dlt,alp,slag)

! initialise co2
call co2_read(sig,jyear,csolar)
rrco2 = rrvco2*ratco2mw

! initialise ozone
if ( amipo3 ) then
  if ( myid==0 ) write(6,*) "AMIP2 ozone input"
  call o3read_amip
else
  call o3_read(sig,jyear,jmonth)
end if

! set-up standard pressure levels
allocate( pref(kl+1,2) )
pref(kl+1,1) = 101325.
pref(kl+1,2) = 81060. !=0.8*pref(kl+1,1)
do k = 1,kl
  kr = kl + 1 - k
  pref(kr,:) = sig(k)*pref(kl+1,:)
end do  
  
Cldrad_control%do_strat_clouds_iz       = .true.
Cldrad_control%do_sw_micro_iz           = .true.
Cldrad_control%do_lw_micro_iz           = .true.
Cldrad_control%do_sw_micro              = .true.
Cldrad_control%do_lw_micro              = .true.
Cldrad_control%do_ica_calcs             = .false. ! must change allocations below if true
Cldrad_control%do_no_clouds             = .false.
Cldrad_control%do_donner_deep_clouds    = .false.
Cldrad_control%do_stochastic_clouds     = .false.
Cldrad_control%using_fu2007             = .false.
Sw_control%solar_constant               = csolar
Sw_control%do_cmip_diagnostics          = do_aerosol_forcing ! Need for aerosol optical depths
Lw_control%do_lwcldemiss                = .true.
Lw_control%do_o3_iz                     = .true.
Lw_control%do_co2_iz                    = .true.
Lw_control%do_ch4_iz                    = .true.
Lw_control%do_n2o_iz                    = .true.
Lw_control%do_o3                        = .true.
Lw_control%do_co2                       = .true.
Lw_control%do_ch4                       = rrvch4>0.
Lw_control%do_n2o                       = rrvch4>0.
Lw_control%do_h2o                       = .true.
Lw_control%do_cfc                       = rrvch4>0.
Rad_control%using_solar_timeseries_data = .false.
Rad_control%do_totcld_forcing           = do_totcld_forcing
Rad_control%rad_time_step               = nint(real(kountr)*dt)
Rad_control%rad_time_step_iz            = .true.
Rad_control%do_aerosol                  = do_aerosol_forcing
Rad_control%do_swaerosol_forcing        = do_aerosol_forcing
Rad_control%do_lwaerosol_forcing        = do_aerosol_forcing
Rad_control%hires_coszen                = .false.
Rad_control%hires_coszen_iz             = .true.
Rad_control%nzens                       = 1
Rad_control%nzens_iz                    = .true.
Rad_control%using_im_bcsul              = .false.
Rad_control%using_im_bcsul_iz           = .true.

if ( myid==0 ) write(6,*) "Reading SEA-ESF model datasets" 
call sealw99_init(pref,Lw_tables)
call esfsw_parameters_init
call esfsw_driver_init
call microphys_rad_init

deallocate( pref )
  
allocate( Atmos_input(0:maxthreads-1) )
allocate( Rad_gases(0:maxthreads-1) )
allocate( Cloud_microphysics(0:maxthreads-1) )
allocate( Cldrad_props(0:maxthreads-1) )
allocate( Lscrad_props(0:maxthreads-1) )
allocate( Cld_spec(0:maxthreads-1) )
allocate( Surface(0:maxthreads-1) )
allocate( Astro(0:maxthreads-1) )
allocate( Lw_output(0:maxthreads-1) )
allocate( Sw_output(0:maxthreads-1) )
allocate( Aerosol(0:maxthreads-1) )
allocate( Aerosol_diags(0:maxthreads-1) )
allocate( Lw_diagnostics(0:maxthreads-1) ) ! working arrays
allocate( Lw_clouds(0:maxthreads-1) )      ! working arrays
allocate( Optical(0:maxthreads-1) )        ! working arrays
allocate( Gas_tf(0:maxthreads-1) )         ! working arrays

do mythread = 0,maxthreads-1

  allocate ( Atmos_input(mythread)%press(imax, 1, kl+1) )
  allocate ( Atmos_input(mythread)%phalf(imax, 1, kl+1) )
  allocate ( Atmos_input(mythread)%temp(imax, 1, kl+1) )
  allocate ( Atmos_input(mythread)%rh2o(imax, 1, kl) )
  allocate ( Atmos_input(mythread)%rel_hum(imax, 1, kl) )
  allocate ( Atmos_input(mythread)%clouddeltaz(imax, 1,kl) )
  allocate ( Atmos_input(mythread)%deltaz(imax, 1, kl) )
  allocate ( Atmos_input(mythread)%pflux(imax, 1, kl+1) )
  allocate ( Atmos_input(mythread)%tflux(imax, 1, kl+1) )
  allocate ( Atmos_input(mythread)%psfc(imax, 1) )
  allocate ( Atmos_input(mythread)%tsfc(imax, 1) )
  if ( do_aerosol_forcing ) then
    allocate( Atmos_input(mythread)%aerosolrelhum(imax, 1, kl) )
  end if  
  !if (use_co2_tracer_field) then
  !  allocate ( Atmos_input(mythread)%tracer_co2(imax, 1, kl) )
  !endif
  
  allocate( Rad_gases(mythread)%qo3(imax, 1, kl) )

  allocate( Cloud_microphysics(mythread)%size_rain(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%size_drop(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%size_ice(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%size_snow(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%conc_drop(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%conc_ice(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%conc_rain(imax, 1, kl) )
  allocate( Cloud_microphysics(mythread)%conc_snow(imax, 1, kl) )

  allocate( Cldrad_props(mythread)%cldext(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props(mythread)%cldsct(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props(mythread)%cldasymm(imax, 1, kl, Solar_spect%nbands, 1) )
  allocate( Cldrad_props(mythread)%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props(mythread)%cldemiss(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props(mythread)%emmxolw(imax, 1, kl, Cldrad_control%nlwcldb,1) )
  allocate( Cldrad_props(mythread)%emrndlw(imax, 1, kl, Cldrad_control%nlwcldb,1) )

  allocate( Lscrad_props(mythread)%cldext(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props(mythread)%cldsct(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props(mythread)%cldasymm(imax, 1, kl, Solar_spect%nbands) )
  allocate( Lscrad_props(mythread)%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb) )

  allocate( Cld_spec(mythread)%camtsw(imax, 1, kl) )
  allocate( Cld_spec(mythread)%cmxolw(imax, 1, kl) )
  allocate( Cld_spec(mythread)%crndlw(imax, 1, kl) )
  
  allocate( Surface(mythread)%asfc_vis_dir(imax, 1) )
  allocate( Surface(mythread)%asfc_nir_dir(imax, 1) )
  allocate( Surface(mythread)%asfc_vis_dif(imax, 1) )
  allocate( Surface(mythread)%asfc_nir_dif(imax, 1) )

  allocate( Astro(mythread)%cosz(imax, 1) )
  allocate( Astro(mythread)%fracday(imax, 1) )
  Astro(mythread)%rrsun = 1./(r1*r1)

  allocate( Lw_output(mythread)%heatra(imax, 1, kl) )
  allocate( Lw_output(mythread)%flxnet(imax, 1, kl+1) )
  allocate( Lw_output(mythread)%bdy_flx(imax, 1, 4) )
  if ( do_totcld_forcing ) then
    allocate( Lw_output(mythread)%heatracf(imax, 1, kl) )
    allocate( Lw_output(mythread)%flxnetcf(imax, 1, kl+1) )
    allocate( Lw_output(mythread)%bdy_flx_clr(imax, 1, 4) )
  endif

  allocate( Sw_output(mythread)%dfsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%ufsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%dfsw_dir_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%dfsw_dif_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%ufsw_dir_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%ufsw_dif_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%fsw(imax, 1, kl+1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%hsw(imax, 1, kl, Rad_control%nzens) )
  allocate( Sw_output(mythread)%dfsw_vis_sfc(imax, 1, Rad_control%nzens ) )
  allocate( Sw_output(mythread)%ufsw_vis_sfc(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%dfsw_vis_sfc_dir(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%dfsw_vis_sfc_dif(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%ufsw_vis_sfc_dir(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%ufsw_vis_sfc_dif(imax, 1, Rad_control%nzens) )
  allocate( Sw_output(mythread)%bdy_flx(imax, 1, 4, Rad_control%nzens) )
  if ( do_totcld_forcing ) then
    allocate( Sw_output(mythread)%dfswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%ufswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%fswcf(imax, 1, kl+1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%hswcf(imax, 1, kl, Rad_control%nzens) )
    allocate( Sw_output(mythread)%dfsw_dir_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%dfsw_dif_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%dfsw_vis_sfc_clr(imax, 1, Rad_control%nzens) )
    allocate( Sw_output(mythread)%bdy_flx_clr(imax, 1, 4, Rad_control%nzens) )
  end if
  
  if ( do_aerosol_forcing ) then
    allocate( Aerosol(mythread)%aerosol(imax, 1, kl, nfields) )
    
    allocate( Aerosol_diags(mythread)%extopdep(imax, 1, kl, nfields, 10) )
    allocate( Aerosol_diags(mythread)%absopdep(imax, 1, kl, nfields, 10) )
    allocate( Aerosol_diags(mythread)%asymdep(imax, 1, kl, nfields, 10) )
    Aerosol_diags(mythread)%extopdep = 0._8
    Aerosol_diags(mythread)%absopdep = 0._8
    Aerosol_diags(mythread)%asymdep  = 0._8
  end if

  ! assign GHG concentrations
  Rad_gases(mythread)%rrvco2  = real(rrvco2 ,8)
  Rad_gases(mythread)%rrvch4  = real(rrvch4 ,8)
  Rad_gases(mythread)%rrvn2o  = real(rrvn2o ,8)
  Rad_gases(mythread)%rrvf11  = real(rrvf11 ,8)
  Rad_gases(mythread)%rrvf12  = real(rrvf12 ,8)
  Rad_gases(mythread)%rrvf113 = real(rrvf113,8)
  Rad_gases(mythread)%rrvf22  = real(rrvf22 ,8)
  Rad_gases(mythread)%time_varying_co2 = .false.
  Rad_gases(mythread)%time_varying_ch4 = .false.
  Rad_gases(mythread)%time_varying_n2o = .false.
  Rad_gases(mythread)%time_varying_f11 = .false.
  Rad_gases(mythread)%time_varying_f12 = .false.
  Rad_gases(mythread)%time_varying_f113 = .false.
  Rad_gases(mythread)%time_varying_f22 = .false.
  Rad_gases(mythread)%use_model_supplied_co2 = .false.
  if ( mythread==0 ) then
    call sealw99_time_vary(Rad_time, Rad_gases(mythread)) ! sets values in gas_tf.f90
  end if    
  
end do ! mythread = 0,maxthreads-1
  
if ( do_aerosol_forcing ) then
  !if ( Rad_control%using_im_bcsul ) then
  !  allocate( Aerosol_props%sulfate_index(0:100, 0:100) )
  !else
    allocate( Aerosol_props%sulfate_index(0:100, 0:0) )
  !end if
  allocate( Aerosol_props%omphilic_index(0:100) )
  allocate( Aerosol_props%bcphilic_index(0:100) )
  allocate( Aerosol_props%seasalt1_index(0:100) )
  allocate( Aerosol_props%seasalt2_index(0:100) )
  allocate( Aerosol_props%seasalt3_index(0:100) )
  allocate( Aerosol_props%seasalt4_index(0:100) )
  allocate( Aerosol_props%seasalt5_index(0:100) )
  allocate( Aerosol_props%optical_index(nfields) )
  allocate( Aerosol_props%ivol(imax, 1, kl) ) ! currently set to zero so we can ignore threads for now
  allocate( Aerosol_props%aerextband(Solar_spect%nbands, naermodels) )
  allocate( Aerosol_props%aerssalbband(Solar_spect%nbands, naermodels) )
  allocate( Aerosol_props%aerasymmband(Solar_spect%nbands, naermodels) )
  allocate( Aerosol_props%aerssalbbandlw(N_AEROSOL_BANDS, naermodels) )
  allocate( Aerosol_props%aerextbandlw(N_AEROSOL_BANDS, naermodels) )
  allocate( Aerosol_props%aerssalbbandlw_cn(N_AEROSOL_BANDS, naermodels) )
  allocate( Aerosol_props%aerextbandlw_cn(N_AEROSOL_BANDS, naermodels) )
    
  Aerosol_props%sulfate_flag  = 0
  Aerosol_props%omphilic_flag = -1
  Aerosol_props%bcphilic_flag = -2
  Aerosol_props%seasalt1_flag = -3
  Aerosol_props%seasalt2_flag = -4
  Aerosol_props%seasalt3_flag = -5
  Aerosol_props%seasalt4_flag = -6
  Aerosol_props%seasalt5_flag = -7
  Aerosol_props%bc_flag       = -8
  Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS
  Aerosol_props%optical_index(1) = Aerosol_props%sulfate_flag        ! so4
  select case( carbonradmethod )
    case(0,-1)    
      !if ( Rad_control%using_im_bcsul ) then
      !  Aerosol_props%optical_index(2) = Aerosol_props%bc_flag      ! so4/bc mixture
      !  Aerosol_props%optical_index(3) = Aerosol_props%bc_flag      ! so4/bc mixture
      !else
        Aerosol_props%optical_index(2) = 2                           ! black carbon (hydrophobic)
        Aerosol_props%optical_index(3) = Aerosol_props%bcphilic_flag ! black carbon (hydrophillic)
        Aerosol_props%optical_index(4) = 3                           ! organic carbon (hydrophobic)
        Aerosol_props%optical_index(5) = Aerosol_props%omphilic_flag ! organic carbon (hydrophillic)
      !end if    
    case(1)
      Aerosol_props%optical_index(2) = naermodels - 4                ! black carbon (soot)
      Aerosol_props%optical_index(3) = naermodels - 4                ! black carbon (soot)
      Aerosol_props%optical_index(4) = naermodels - 5                ! organic carbon (organic cabron)
      Aerosol_props%optical_index(5) = naermodels - 5                ! organic carbon (organic carbon)
    case default
      write(6,*) "ERROR: Invalid carbonradmethod ",carbonradmethod
      call ccmpi_abort(-1)
  end select
  Aerosol_props%optical_index(6) = naermodels - 3                    ! dust_0.7  (using 0.73)
  Aerosol_props%optical_index(7) = naermodels - 2                    ! dust_1.4  (using 1.4)
  Aerosol_props%optical_index(8) = naermodels - 1                    ! dust_2.4  (using 2.4)
  Aerosol_props%optical_index(9) = naermodels                        ! dust_4.5  (using 4.5)
  !Aerosol_props%optical_index(10) = 1                               ! sea_salt (film drop + jet drop)
  Aerosol_props%optical_index(10) = Aerosol_props%seasalt1_flag      ! sea salt film drop (0.1)
  Aerosol_props%optical_index(11) = Aerosol_props%seasalt2_flag      ! sea salt jet drop (0.5)
  ! GFDL bins dust1=0.1-0.5, dust2=0.5-1, dust3=1-2.5, dust4=2.5-5, dust5=5-10
  ! GFDL bins salt1=0.1-0.5, salt2=0.5-1, salt3=1-2.5, salt4=2.5-5, dust5=5-10

  Aerosol_props%sulfate_index( 0: 13, 0) = (/ 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62 /)
  Aerosol_props%sulfate_index(14: 27, 0) = (/ 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62 /)
  Aerosol_props%sulfate_index(28: 41, 0) = (/ 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 64, 64, 64, 64 /)
  Aerosol_props%sulfate_index(42: 55, 0) = (/ 64, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 67, 67, 67 /)
  Aerosol_props%sulfate_index(56: 69, 0) = (/ 67, 67, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 70, 70 /)
  Aerosol_props%sulfate_index(70: 83, 0) = (/ 70, 70, 70, 71, 71, 71, 71, 71, 72, 72, 72, 72, 73, 73 /)
  Aerosol_props%sulfate_index(84: 97, 0) = (/ 74, 74, 75, 75, 76, 76, 77, 78, 79, 80, 81, 82, 83, 84 /)
  Aerosol_props%sulfate_index(98:100, 0) = (/ 84, 84, 84 /)
  !if ( Rad_control%using_im_bcsul ) then
  !  Aerosol_props%sulfate_index(:, 1)=Aerosol_props%sulfate_index(:, 0)      
  !  Aerosol_props%sulfate_index(:, 2)=Aerosol_props%sulfate_index(:, 1) + 26 ! 98%
  !  Aerosol_props%sulfate_index(:, 3)=Aerosol_props%sulfate_index(:, 2)
  !  Aerosol_props%sulfate_index(:, 4)=Aerosol_props%sulfate_index(:, 3) + 26 ! 96%
  !  Aerosol_props%sulfate_index(:, 5)=Aerosol_props%sulfate_index(:, 4)
  !  Aerosol_props%sulfate_index(:, 6)=Aerosol_props%sulfate_index(:, 5) + 26 ! 94%
  !  Aerosol_props%sulfate_index(:, 7)=Aerosol_props%sulfate_index(:, 6)
  !  Aerosol_props%sulfate_index(:, 8)=Aerosol_props%sulfate_index(:, 7) + 26 ! 92%
  !  Aerosol_props%sulfate_index(:, 9)=Aerosol_props%sulfate_index(:, 8)
  !  Aerosol_props%sulfate_index(:,10)=Aerosol_props%sulfate_index(:, 9) + 26 ! 90%
  !  Aerosol_props%sulfate_index(:,11)=Aerosol_props%sulfate_index(:,10)
  !  Aerosol_props%sulfate_index(:,12)=Aerosol_props%sulfate_index(:,11) + 26 ! 88%
  !  Aerosol_props%sulfate_index(:,13)=Aerosol_props%sulfate_index(:,12)
  !  Aerosol_props%sulfate_index(:,14)=Aerosol_props%sulfate_index(:,13) + 26 ! 86%
  !  Aerosol_props%sulfate_index(:,15)=Aerosol_props%sulfate_index(:,14)
  !  Aerosol_props%sulfate_index(:,16)=Aerosol_props%sulfate_index(:,15) + 26 ! 84%
  !  Aerosol_props%sulfate_index(:,17)=Aerosol_props%sulfate_index(:,16)
  !  Aerosol_props%sulfate_index(:,18)=Aerosol_props%sulfate_index(:,17) + 26 ! 82%
  !  Aerosol_props%sulfate_index(:,19)=Aerosol_props%sulfate_index(:,18)
  !  Aerosol_props%sulfate_index(:,20)=Aerosol_props%sulfate_index(:,19) + 26 ! 80%
  !  Aerosol_props%sulfate_index(:,21)=Aerosol_props%sulfate_index(:,20)
  !  Aerosol_props%sulfate_index(:,22)=Aerosol_props%sulfate_index(:,21)
  !  Aerosol_props%sulfate_index(:,23)=Aerosol_props%sulfate_index(:,22) + 26
  !  Aerosol_props%sulfate_index(:,24)=Aerosol_props%sulfate_index(:,23)
  !  Aerosol_props%sulfate_index(:,25)=Aerosol_props%sulfate_index(:,24)      ! 75%
  !  Aerosol_props%sulfate_index(:,26)=Aerosol_props%sulfate_index(:,25)
  !  Aerosol_props%sulfate_index(:,27)=Aerosol_props%sulfate_index(:,26)
  !  Aerosol_props%sulfate_index(:,28)=Aerosol_props%sulfate_index(:,27) + 26
  !  Aerosol_props%sulfate_index(:,29)=Aerosol_props%sulfate_index(:,28)
  !  Aerosol_props%sulfate_index(:,30)=Aerosol_props%sulfate_index(:,29)      ! 70%
  !  Aerosol_props%sulfate_index(:,31)=Aerosol_props%sulfate_index(:,30)
  !  Aerosol_props%sulfate_index(:,32)=Aerosol_props%sulfate_index(:,31)
  !  Aerosol_props%sulfate_index(:,33)=Aerosol_props%sulfate_index(:,32) + 26
  !  Aerosol_props%sulfate_index(:,34)=Aerosol_props%sulfate_index(:,33)
  !  Aerosol_props%sulfate_index(:,35)=Aerosol_props%sulfate_index(:,34)      ! 65%
  !  Aerosol_props%sulfate_index(:,36)=Aerosol_props%sulfate_index(:,35)
  !  Aerosol_props%sulfate_index(:,37)=Aerosol_props%sulfate_index(:,36)
  !  Aerosol_props%sulfate_index(:,38)=Aerosol_props%sulfate_index(:,37) + 26
  !  Aerosol_props%sulfate_index(:,39)=Aerosol_props%sulfate_index(:,38)
  !  Aerosol_props%sulfate_index(:,40)=Aerosol_props%sulfate_index(:,39)      ! 60%
  !  Aerosol_props%sulfate_index(:,41)=Aerosol_props%sulfate_index(:,40)
  !  Aerosol_props%sulfate_index(:,42)=Aerosol_props%sulfate_index(:,41)
  !  Aerosol_props%sulfate_index(:,43)=Aerosol_props%sulfate_index(:,42) + 26
  !  Aerosol_props%sulfate_index(:,44)=Aerosol_props%sulfate_index(:,43)
  !  Aerosol_props%sulfate_index(:,45)=Aerosol_props%sulfate_index(:,44)      ! 55%
  !  Aerosol_props%sulfate_index(:,46)=Aerosol_props%sulfate_index(:,45)
  !  Aerosol_props%sulfate_index(:,47)=Aerosol_props%sulfate_index(:,46)
  !  Aerosol_props%sulfate_index(:,48)=Aerosol_props%sulfate_index(:,47) + 26
  !  Aerosol_props%sulfate_index(:,49)=Aerosol_props%sulfate_index(:,48)
  !  Aerosol_props%sulfate_index(:,50)=Aerosol_props%sulfate_index(:,49)      ! 50%
  !  Aerosol_props%sulfate_index(:,51)=Aerosol_props%sulfate_index(:,50)
  !  Aerosol_props%sulfate_index(:,52)=Aerosol_props%sulfate_index(:,51)
  !  Aerosol_props%sulfate_index(:,53)=Aerosol_props%sulfate_index(:,52) + 26
  !  Aerosol_props%sulfate_index(:,54)=Aerosol_props%sulfate_index(:,53)
  !  Aerosol_props%sulfate_index(:,55)=Aerosol_props%sulfate_index(:,54)      ! 45%
  !  Aerosol_props%sulfate_index(:,56)=Aerosol_props%sulfate_index(:,55)
  !  Aerosol_props%sulfate_index(:,57)=Aerosol_props%sulfate_index(:,56)
  !  Aerosol_props%sulfate_index(:,58)=Aerosol_props%sulfate_index(:,57) + 26
  !  Aerosol_props%sulfate_index(:,59)=Aerosol_props%sulfate_index(:,58)
  !  Aerosol_props%sulfate_index(:,60)=Aerosol_props%sulfate_index(:,59)      ! 40%
  !  Aerosol_props%sulfate_index(:,61)=Aerosol_props%sulfate_index(:,60)
  !  Aerosol_props%sulfate_index(:,62)=Aerosol_props%sulfate_index(:,61)
  !  Aerosol_props%sulfate_index(:,63)=Aerosol_props%sulfate_index(:,62) + 26
  !  Aerosol_props%sulfate_index(:,64)=Aerosol_props%sulfate_index(:,63)
  !  Aerosol_props%sulfate_index(:,65)=Aerosol_props%sulfate_index(:,64)      ! 35%
  !  Aerosol_props%sulfate_index(:,66)=Aerosol_props%sulfate_index(:,65)
  !  Aerosol_props%sulfate_index(:,67)=Aerosol_props%sulfate_index(:,66)
  !  Aerosol_props%sulfate_index(:,68)=Aerosol_props%sulfate_index(:,67) + 26
  !  Aerosol_props%sulfate_index(:,69)=Aerosol_props%sulfate_index(:,68)
  !  Aerosol_props%sulfate_index(:,70)=Aerosol_props%sulfate_index(:,69)      ! 30%
  !  Aerosol_props%sulfate_index(:,71)=Aerosol_props%sulfate_index(:,70)
  !  Aerosol_props%sulfate_index(:,72)=Aerosol_props%sulfate_index(:,71)
  !  Aerosol_props%sulfate_index(:,73)=Aerosol_props%sulfate_index(:,72) + 26
  !  Aerosol_props%sulfate_index(:,74)=Aerosol_props%sulfate_index(:,73)
  !  Aerosol_props%sulfate_index(:,75)=Aerosol_props%sulfate_index(:,74)      ! 25%
  !  Aerosol_props%sulfate_index(:,76)=Aerosol_props%sulfate_index(:,75)
  !  Aerosol_props%sulfate_index(:,77)=Aerosol_props%sulfate_index(:,76)
  !  Aerosol_props%sulfate_index(:,78)=Aerosol_props%sulfate_index(:,77) + 26
  !  Aerosol_props%sulfate_index(:,79)=Aerosol_props%sulfate_index(:,78)
  !  Aerosol_props%sulfate_index(:,80)=Aerosol_props%sulfate_index(:,79)      ! 20%
  !  Aerosol_props%sulfate_index(:,81)=Aerosol_props%sulfate_index(:,80)
  !  Aerosol_props%sulfate_index(:,82)=Aerosol_props%sulfate_index(:,81)
  !  Aerosol_props%sulfate_index(:,83)=Aerosol_props%sulfate_index(:,82) + 26
  !  Aerosol_props%sulfate_index(:,84)=Aerosol_props%sulfate_index(:,83)
  !  Aerosol_props%sulfate_index(:,85)=Aerosol_props%sulfate_index(:,84)      ! 15%
  !  Aerosol_props%sulfate_index(:,86)=Aerosol_props%sulfate_index(:,85)
  !  Aerosol_props%sulfate_index(:,87)=Aerosol_props%sulfate_index(:,86)
  !  Aerosol_props%sulfate_index(:,88)=Aerosol_props%sulfate_index(:,87) + 26
  !  Aerosol_props%sulfate_index(:,89)=Aerosol_props%sulfate_index(:,88)
  !  Aerosol_props%sulfate_index(:,90)=Aerosol_props%sulfate_index(:,89)      ! 10%
  !  Aerosol_props%sulfate_index(:,91)=Aerosol_props%sulfate_index(:,90)
  !  Aerosol_props%sulfate_index(:,92)=Aerosol_props%sulfate_index(:,91)
  !  Aerosol_props%sulfate_index(:,93)=Aerosol_props%sulfate_index(:,92) + 26
  !  Aerosol_props%sulfate_index(:,94)=Aerosol_props%sulfate_index(:,93)
  !  Aerosol_props%sulfate_index(:,95)=Aerosol_props%sulfate_index(:,94)      ! 5%
  !  Aerosol_props%sulfate_index(:,96)=Aerosol_props%sulfate_index(:,95)
  !  Aerosol_props%sulfate_index(:,97)=Aerosol_props%sulfate_index(:,96)
  !  Aerosol_props%sulfate_index(:,98)=Aerosol_props%sulfate_index(:,97) + 26
  !  Aerosol_props%sulfate_index(:,99)=Aerosol_props%sulfate_index(:,98)
  !  Aerosol_props%sulfate_index(:,100)=Aerosol_props%sulfate_index(:,99)     ! 0%
  !end if
  Aerosol_props%bcphilic_index( 0: 13)=(/  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12 /)
  Aerosol_props%bcphilic_index(14: 27)=(/  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12 /)
  Aerosol_props%bcphilic_index(28: 41)=(/  12,  12,  12,  12,  12,  13,  13,  13,  13,  13,  14,  14,  14,  14 /)
  Aerosol_props%bcphilic_index(42: 55)=(/  14,  15,  15,  15,  15,  15,  16,  16,  16,  16,  16,  17,  17,  17 /)
  Aerosol_props%bcphilic_index(56: 69)=(/  17,  17,  18,  18,  18,  18,  18,  19,  19,  19,  19,  19,  20,  20 /)
  Aerosol_props%bcphilic_index(70: 83)=(/  20,  20,  20,  21,  21,  21,  21,  21,  22,  22,  22,  22,  23,  23 /)
  Aerosol_props%bcphilic_index(84: 97)=(/  24,  24,  25,  25,  26,  26,  27,  28,  29,  30,  31,  32,  33,  34 /)
  Aerosol_props%bcphilic_index(98:100)=(/  34,  34,  34 /)
  Aerosol_props%omphilic_index(:) = Aerosol_props%bcphilic_index(:) + 25
  Aerosol_props%seasalt1_index( 0: 13)=(/  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88 /)
  Aerosol_props%seasalt1_index(14: 27)=(/  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88,  88 /)
  Aerosol_props%seasalt1_index(28: 41)=(/  88,  88,  88,  88,  88,  89,  89,  89,  89,  89,  90,  90,  90,  90 /)
  Aerosol_props%seasalt1_index(42: 55)=(/  90,  91,  91,  91,  91,  91,  92,  92,  92,  92,  92,  93,  93,  93 /)
  Aerosol_props%seasalt1_index(56: 69)=(/  93,  93,  94,  94,  94,  94,  94,  95,  95,  95,  95,  95,  96,  96 /)
  Aerosol_props%seasalt1_index(70: 83)=(/  96,  96,  96,  97,  97,  97,  97,  97,  98,  98,  98,  98,  99,  99 /)
  Aerosol_props%seasalt1_index(84: 97)=(/ 100, 100, 101, 101, 102, 102, 103, 104, 105, 106, 107, 108, 109, 110 /)
  Aerosol_props%seasalt1_index(98:100)=(/ 111, 111, 111 /)
  Aerosol_props%seasalt2_index(:) = Aerosol_props%seasalt1_index(:) + 25
  Aerosol_props%seasalt3_index(:) = Aerosol_props%seasalt2_index(:) + 25
  Aerosol_props%seasalt4_index(:) = Aerosol_props%seasalt3_index(:) + 25
  Aerosol_props%seasalt5_index(:) = Aerosol_props%seasalt4_index(:) + 25
  call loadaerooptical(Aerosol_props)
    
  if ( include_volcanoes ) then
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
swrsave(:) = 0.5

! error checking
!if ( ldr==0 ) then
!  write(6,*) "ERROR: SEA-ESF radiation requires ldr/=0"
!  call ccmpi_abort(-1)
!end if
  
call END_LOG(radinit_end)

return
end subroutine seaesfrad_init

subroutine loadaerooptical(Aerosol_props)

use aerosolldr
use cc_mpi
use filnames_m
use infile

implicit none

integer :: n, nmodel, unit, num_wavenumbers, num_input_categories
integer :: noptical, nivl3, nband, nw, ierr, na, ni
integer, dimension(:), allocatable, save :: nivl1aero, nivl2aero
integer, dimension(:), allocatable, save :: endaerwvnsf
real(kind=8), dimension(:,:), allocatable, save :: aeroextivl, aerossalbivl, aeroasymmivl
real(kind=8), dimension(:,:), allocatable, save :: sflwwts, sflwwts_cn
real(kind=8), dimension(:,:), allocatable, save :: solivlaero
real(kind=8), dimension(:), allocatable, save :: aeroext_in, aerossalb_in, aeroasymm_in
real(kind=8) :: sumsol3, frac
real(kind=8), dimension(:,:), allocatable, save :: dum_aero
logical, dimension(naermodels-4) :: optical_check
character(len=64), dimension(naermodels) :: aerosol_optical_names
character(len=64) :: name_in
character(len=110) :: filename
type(aerosol_properties_type), intent(inout) :: Aerosol_props

aerosol_optical_names(1)=       "sea_salt"
aerosol_optical_names(2:3)=  (/ "bcphobic",    "omphobic" /)
aerosol_optical_names(4:5)=  (/ "dust_0.1",    "dust_0.2" /)
aerosol_optical_names(6:7)=  (/ "dust_0.4",    "dust_0.8" /)
aerosol_optical_names(8:9)=  (/ "dust_1.0",    "dust_2.0" /)
aerosol_optical_names(10:11)=(/ "dust_4.0",    "dust_8.0" /)
aerosol_optical_names(12:13)=(/ "bcphilic_30%", "bcphilic_35%" /)
aerosol_optical_names(14:15)=(/ "bcphilic_40%", "bcphilic_45%" /)
aerosol_optical_names(16:17)=(/ "bcphilic_50%", "bcphilic_55%" /)
aerosol_optical_names(18:19)=(/ "bcphilic_60%", "bcphilic_65%" /)
aerosol_optical_names(20:21)=(/ "bcphilic_70%", "bcphilic_75%" /)
aerosol_optical_names(22:23)=(/ "bcphilic_80%", "bcphilic_82%" /)
aerosol_optical_names(24:25)=(/ "bcphilic_84%", "bcphilic_86%" /)
aerosol_optical_names(26:27)=(/ "bcphilic_88%", "bcphilic_90%" /)
aerosol_optical_names(28:29)=(/ "bcphilic_91%", "bcphilic_92%" /)
aerosol_optical_names(30:31)=(/ "bcphilic_93%", "bcphilic_94%" /)
aerosol_optical_names(32:33)=(/ "bcphilic_95%", "bcphilic_96%" /)
aerosol_optical_names(34:35)=(/ "bcphilic_97%", "bcphilic_98%" /)
aerosol_optical_names(36)=      "bcphilic_99%"
aerosol_optical_names(37:38)=(/ "omphilic_30%", "omphilic_35%" /)
aerosol_optical_names(39:40)=(/ "omphilic_40%", "omphilic_45%" /)
aerosol_optical_names(41:42)=(/ "omphilic_50%", "omphilic_55%" /)
aerosol_optical_names(43:44)=(/ "omphilic_60%", "omphilic_65%" /)
aerosol_optical_names(45:46)=(/ "omphilic_70%", "omphilic_75%" /)
aerosol_optical_names(47:48)=(/ "omphilic_80%", "omphilic_82%" /)
aerosol_optical_names(49:50)=(/ "omphilic_84%", "omphilic_86%" /)
aerosol_optical_names(51:52)=(/ "omphilic_88%", "omphilic_90%" /)
aerosol_optical_names(53:54)=(/ "omphilic_91%", "omphilic_92%" /)
aerosol_optical_names(55:56)=(/ "omphilic_93%", "omphilic_94%" /)
aerosol_optical_names(57:58)=(/ "omphilic_95%", "omphilic_96%" /)
aerosol_optical_names(59:60)=(/ "omphilic_97%", "omphilic_98%" /)
aerosol_optical_names(61)=      "omphilic_99%"
aerosol_optical_names(62:63)=(/ "sulfate_30%_100%", "sulfate_35%_100%" /)
aerosol_optical_names(64:65)=(/ "sulfate_40%_100%", "sulfate_45%_100%" /)
aerosol_optical_names(66:67)=(/ "sulfate_50%_100%", "sulfate_55%_100%" /)
aerosol_optical_names(68:69)=(/ "sulfate_60%_100%", "sulfate_65%_100%" /)
aerosol_optical_names(70:71)=(/ "sulfate_70%_100%", "sulfate_75%_100%" /)
aerosol_optical_names(72:73)=(/ "sulfate_80%_100%", "sulfate_82%_100%" /)
aerosol_optical_names(74:75)=(/ "sulfate_84%_100%", "sulfate_86%_100%" /)
aerosol_optical_names(76:77)=(/ "sulfate_88%_100%", "sulfate_90%_100%" /)
aerosol_optical_names(78:79)=(/ "sulfate_91%_100%", "sulfate_92%_100%" /)
aerosol_optical_names(80:81)=(/ "sulfate_93%_100%", "sulfate_94%_100%" /)
aerosol_optical_names(82:83)=(/ "sulfate_95%_100%", "sulfate_96%_100%" /)
aerosol_optical_names(84:85)=(/ "sulfate_97%_100%", "sulfate_98%_100%" /)
aerosol_optical_names(86)=      "sulfate_99%_100%"
aerosol_optical_names(87)=      "sulfate_100%_100%"
aerosol_optical_names( 88: 91)=(/ "seasalt1_30%", "seasalt1_35%", "seasalt1_40%", "seasalt1_45%" /)
aerosol_optical_names( 92: 95)=(/ "seasalt1_50%", "seasalt1_55%", "seasalt1_60%", "seasalt1_65%" /)
aerosol_optical_names( 96: 99)=(/ "seasalt1_70%", "seasalt1_75%", "seasalt1_80%", "seasalt1_82%" /)
aerosol_optical_names(100:103)=(/ "seasalt1_84%", "seasalt1_86%", "seasalt1_88%", "seasalt1_90%" /)
aerosol_optical_names(104:107)=(/ "seasalt1_91%", "seasalt1_92%", "seasalt1_93%", "seasalt1_94%" /)
aerosol_optical_names(108:111)=(/ "seasalt1_95%", "seasalt1_96%", "seasalt1_97%", "seasalt1_98%" /)
aerosol_optical_names(112)=       "seasalt1_99%"
aerosol_optical_names(113:116)=(/ "seasalt2_30%", "seasalt2_35%", "seasalt2_40%", "seasalt2_45%" /)
aerosol_optical_names(117:120)=(/ "seasalt2_50%", "seasalt2_55%", "seasalt2_60%", "seasalt2_65%" /)
aerosol_optical_names(121:124)=(/ "seasalt2_70%", "seasalt2_75%", "seasalt2_80%", "seasalt2_82%" /)
aerosol_optical_names(125:128)=(/ "seasalt2_84%", "seasalt2_86%", "seasalt2_88%", "seasalt2_90%" /)
aerosol_optical_names(129:132)=(/ "seasalt2_91%", "seasalt2_92%", "seasalt2_93%", "seasalt2_94%" /)
aerosol_optical_names(133:136)=(/ "seasalt2_95%", "seasalt2_96%", "seasalt2_97%", "seasalt2_98%" /)
aerosol_optical_names(137)=       "seasalt2_99%"
aerosol_optical_names(138:141)=(/ "seasalt3_30%", "seasalt3_35%", "seasalt3_40%", "seasalt3_45%" /)
aerosol_optical_names(142:145)=(/ "seasalt3_50%", "seasalt3_55%", "seasalt3_60%", "seasalt3_65%" /)
aerosol_optical_names(146:149)=(/ "seasalt3_70%", "seasalt3_75%", "seasalt3_80%", "seasalt3_82%" /)
aerosol_optical_names(150:153)=(/ "seasalt3_84%", "seasalt3_86%", "seasalt3_88%", "seasalt3_90%" /)
aerosol_optical_names(154:157)=(/ "seasalt3_91%", "seasalt3_92%", "seasalt3_93%", "seasalt3_94%" /)
aerosol_optical_names(158:161)=(/ "seasalt3_95%", "seasalt3_96%", "seasalt3_97%", "seasalt3_98%" /)
aerosol_optical_names(162)=       "seasalt3_99%"
aerosol_optical_names(163:166)=(/ "seasalt4_30%", "seasalt4_35%", "seasalt4_40%", "seasalt4_45%" /)
aerosol_optical_names(167:170)=(/ "seasalt4_50%", "seasalt4_55%", "seasalt4_60%", "seasalt4_65%" /)
aerosol_optical_names(171:174)=(/ "seasalt4_70%", "seasalt4_75%", "seasalt4_80%", "seasalt4_82%" /)
aerosol_optical_names(175:178)=(/ "seasalt4_84%", "seasalt4_86%", "seasalt4_88%", "seasalt4_90%" /)
aerosol_optical_names(179:182)=(/ "seasalt4_91%", "seasalt4_92%", "seasalt4_93%", "seasalt4_94%" /)
aerosol_optical_names(183:186)=(/ "seasalt4_95%", "seasalt4_96%", "seasalt4_97%", "seasalt4_98%" /)
aerosol_optical_names(187)=       "seasalt4_99%"
aerosol_optical_names(188:191)=(/ "seasalt5_30%", "seasalt5_35%", "seasalt5_40%", "seasalt5_45%" /)
aerosol_optical_names(192:195)=(/ "seasalt5_50%", "seasalt5_55%", "seasalt5_60%", "seasalt5_65%" /)
aerosol_optical_names(196:199)=(/ "seasalt5_70%", "seasalt5_75%", "seasalt5_80%", "seasalt5_82%" /)
aerosol_optical_names(200:203)=(/ "seasalt5_84%", "seasalt5_86%", "seasalt5_88%", "seasalt5_90%" /)
aerosol_optical_names(204:207)=(/ "seasalt5_91%", "seasalt5_92%", "seasalt5_93%", "seasalt5_94%" /)
aerosol_optical_names(208:211)=(/ "seasalt5_95%", "seasalt5_96%", "seasalt5_97%", "seasalt5_98%" /)
aerosol_optical_names(212)=       "seasalt5_99%"
!aerosol_optical_names( 88: 91)=(/ "sulfate_30%_98%", "sulfate_35%_98%", "sulfate_40%_98%", "sulfate_45%_98%" /)
!aerosol_optical_names( 92: 95)=(/ "sulfate_50%_98%", "sulfate_55%_98%", "sulfate_60%_98%", "sulfate_65%_98%" /)
!aerosol_optical_names( 96: 99)=(/ "sulfate_70%_98%", "sulfate_75%_98%", "sulfate_80%_98%", "sulfate_82%_98%" /)
!aerosol_optical_names(100:103)=(/ "sulfate_84%_98%", "sulfate_86%_98%", "sulfate_88%_98%", "sulfate_90%_98%" /)
!aerosol_optical_names(104:107)=(/ "sulfate_91%_98%", "sulfate_92%_98%", "sulfate_93%_98%", "sulfate_94%_98%" /)
!aerosol_optical_names(108:111)=(/ "sulfate_95%_98%", "sulfate_96%_98%", "sulfate_97%_98%", "sulfate_98%_98%" /)
!aerosol_optical_names(112)=       "sulfate_99%_98%"
!aerosol_optical_names(113)=       "sulfate_100%_98%"
!aerosol_optical_names(114:117)=(/ "sulfate_30%_96%", "sulfate_35%_96%", "sulfate_40%_96%", "sulfate_45%_96%" /)
!aerosol_optical_names(118:121)=(/ "sulfate_50%_96%", "sulfate_55%_96%", "sulfate_60%_96%", "sulfate_65%_96%" /)
!aerosol_optical_names(122:125)=(/ "sulfate_70%_96%", "sulfate_75%_96%", "sulfate_80%_96%", "sulfate_82%_96%" /)
!aerosol_optical_names(126:129)=(/ "sulfate_84%_96%", "sulfate_86%_96%", "sulfate_88%_96%", "sulfate_90%_96%" /)
!aerosol_optical_names(130:133)=(/ "sulfate_91%_96%", "sulfate_92%_96%", "sulfate_93%_96%", "sulfate_94%_96%" /)
!aerosol_optical_names(134:137)=(/ "sulfate_95%_96%", "sulfate_96%_96%", "sulfate_97%_96%", "sulfate_98%_96%" /)
!aerosol_optical_names(138)=       "sulfate_99%_96%"
!aerosol_optical_names(139)=       "sulfate_100%_96%"
!aerosol_optical_names(140:143)=(/ "sulfate_30%_94%", "sulfate_35%_94%", "sulfate_40%_94%", "sulfate_45%_94%" /)
!aerosol_optical_names(144:147)=(/ "sulfate_50%_94%", "sulfate_55%_94%", "sulfate_60%_94%", "sulfate_65%_94%" /)
!aerosol_optical_names(148:151)=(/ "sulfate_70%_94%", "sulfate_75%_94%", "sulfate_80%_94%", "sulfate_82%_94%" /)
!aerosol_optical_names(152:155)=(/ "sulfate_84%_94%", "sulfate_86%_94%", "sulfate_88%_94%", "sulfate_90%_94%" /)
!aerosol_optical_names(156:159)=(/ "sulfate_91%_94%", "sulfate_92%_94%", "sulfate_93%_94%", "sulfate_94%_94%" /)
!aerosol_optical_names(160:163)=(/ "sulfate_95%_94%", "sulfate_96%_94%", "sulfate_97%_94%", "sulfate_98%_94%" /)
!aerosol_optical_names(164)=       "sulfate_99%_94%"
!aerosol_optical_names(165)=       "sulfate_100%_94%"
!aerosol_optical_names(166:169)=(/ "sulfate_30%_92%", "sulfate_35%_92%", "sulfate_40%_92%", "sulfate_45%_92%" /)
!aerosol_optical_names(170:173)=(/ "sulfate_50%_92%", "sulfate_55%_92%", "sulfate_60%_92%", "sulfate_65%_92%" /)
!aerosol_optical_names(174:177)=(/ "sulfate_70%_92%", "sulfate_75%_92%", "sulfate_80%_92%", "sulfate_82%_92%" /)
!aerosol_optical_names(178:180)=(/ "sulfate_84%_92%", "sulfate_86%_92%", "sulfate_88%_92%", "sulfate_90%_92%" /)
!aerosol_optical_names(182:185)=(/ "sulfate_91%_92%", "sulfate_92%_92%", "sulfate_93%_92%", "sulfate_94%_92%" /)
!aerosol_optical_names(186:189)=(/ "sulfate_95%_92%", "sulfate_96%_92%", "sulfate_97%_92%", "sulfate_98%_92%" /)
!aerosol_optical_names(190)=       "sulfate_99%_92%"
!aerosol_optical_names(191)=       "sulfate_100%_92%"
!aerosol_optical_names(192:195)=(/ "sulfate_30%_90%", "sulfate_35%_90%", "sulfate_40%_90%", "sulfate_45%_90%" /)
!aerosol_optical_names(196:199)=(/ "sulfate_50%_90%", "sulfate_55%_90%", "sulfate_60%_90%", "sulfate_65%_90%" /)
!aerosol_optical_names(200:203)=(/ "sulfate_70%_90%", "sulfate_75%_90%", "sulfate_80%_90%", "sulfate_82%_90%" /)
!aerosol_optical_names(204:207)=(/ "sulfate_84%_90%", "sulfate_86%_90%", "sulfate_88%_90%", "sulfate_90%_90%" /)
!aerosol_optical_names(208:211)=(/ "sulfate_91%_90%", "sulfate_92%_90%", "sulfate_93%_90%", "sulfate_94%_90%" /)
!aerosol_optical_names(212:215)=(/ "sulfate_95%_90%", "sulfate_96%_90%", "sulfate_97%_90%", "sulfate_98%_90%" /)
!aerosol_optical_names(216)=       "sulfate_99%_90%"
!aerosol_optical_names(217)=       "sulfate_100%_90%"
!aerosol_optical_names(218:221)=(/ "sulfate_30%_88%", "sulfate_35%_88%", "sulfate_40%_88%", "sulfate_45%_88%" /)
!aerosol_optical_names(222:225)=(/ "sulfate_50%_88%", "sulfate_55%_88%", "sulfate_60%_88%", "sulfate_65%_88%" /)
!aerosol_optical_names(226:229)=(/ "sulfate_70%_88%", "sulfate_75%_88%", "sulfate_80%_88%", "sulfate_82%_88%" /)
!aerosol_optical_names(230:233)=(/ "sulfate_84%_88%", "sulfate_86%_88%", "sulfate_88%_88%", "sulfate_90%_88%" /)
!aerosol_optical_names(234:237)=(/ "sulfate_91%_88%", "sulfate_92%_88%", "sulfate_93%_88%", "sulfate_94%_88%" /)
!aerosol_optical_names(238:241)=(/ "sulfate_95%_88%", "sulfate_96%_88%", "sulfate_97%_88%", "sulfate_98%_88%" /)
!aerosol_optical_names(242)=       "sulfate_99%_88%"
!aerosol_optical_names(243)=       "sulfate_100%_88%"
!aerosol_optical_names(244:247)=(/ "sulfate_30%_86%", "sulfate_35%_86%", "sulfate_40%_86%", "sulfate_45%_86%" /)
!aerosol_optical_names(248:251)=(/ "sulfate_50%_86%", "sulfate_55%_86%", "sulfate_60%_86%", "sulfate_65%_86%" /)
!aerosol_optical_names(252:255)=(/ "sulfate_70%_86%", "sulfate_75%_86%", "sulfate_80%_86%", "sulfate_82%_86%" /)
!aerosol_optical_names(256:259)=(/ "sulfate_84%_86%", "sulfate_86%_86%", "sulfate_88%_86%", "sulfate_90%_86%" /)
!aerosol_optical_names(260:263)=(/ "sulfate_91%_86%", "sulfate_92%_86%", "sulfate_93%_86%", "sulfate_94%_86%" /)
!aerosol_optical_names(264:267)=(/ "sulfate_95%_86%", "sulfate_96%_86%", "sulfate_97%_86%", "sulfate_98%_86%" /)
!aerosol_optical_names(268)=       "sulfate_99%_86%"
!aerosol_optical_names(269)=       "sulfate_100%_86%"
!aerosol_optical_names(270:273)=(/ "sulfate_30%_84%", "sulfate_35%_84%", "sulfate_40%_84%", "sulfate_45%_84%" /)
!aerosol_optical_names(274:277)=(/ "sulfate_50%_84%", "sulfate_55%_84%", "sulfate_60%_84%", "sulfate_65%_84%" /)
!aerosol_optical_names(278:281)=(/ "sulfate_70%_84%", "sulfate_75%_84%", "sulfate_80%_84%", "sulfate_82%_84%" /)
!aerosol_optical_names(282:285)=(/ "sulfate_84%_84%", "sulfate_86%_84%", "sulfate_88%_84%", "sulfate_90%_84%" /)
!aerosol_optical_names(286:289)=(/ "sulfate_91%_84%", "sulfate_92%_84%", "sulfate_93%_84%", "sulfate_94%_84%" /)
!aerosol_optical_names(290:293)=(/ "sulfate_95%_84%", "sulfate_96%_84%", "sulfate_97%_84%", "sulfate_98%_84%" /)
!aerosol_optical_names(294)=       "sulfate_99%_84%"
!aerosol_optical_names(295)=       "sulfate_100%_84%"
!aerosol_optical_names(296:299)=(/ "sulfate_30%_82%", "sulfate_35%_82%", "sulfate_40%_82%", "sulfate_45%_82%" /)
!aerosol_optical_names(300:303)=(/ "sulfate_50%_82%", "sulfate_55%_82%", "sulfate_60%_82%", "sulfate_65%_82%" /)
!aerosol_optical_names(304:307)=(/ "sulfate_70%_82%", "sulfate_75%_82%", "sulfate_80%_82%", "sulfate_82%_82%" /)
!aerosol_optical_names(308:311)=(/ "sulfate_84%_82%", "sulfate_86%_82%", "sulfate_88%_82%", "sulfate_90%_82%" /)
!aerosol_optical_names(312:315)=(/ "sulfate_91%_82%", "sulfate_92%_82%", "sulfate_93%_82%", "sulfate_94%_82%" /)
!aerosol_optical_names(316:319)=(/ "sulfate_95%_82%", "sulfate_96%_82%", "sulfate_97%_82%", "sulfate_98%_82%" /)
!aerosol_optical_names(320)=       "sulfate_99%_82%"
!aerosol_optical_names(321)=       "sulfate_100%_82%"
!aerosol_optical_names(322:325)=(/ "sulfate_30%_80%", "sulfate_35%_80%", "sulfate_40%_80%", "sulfate_45%_80%" /)
!aerosol_optical_names(326:329)=(/ "sulfate_50%_80%", "sulfate_55%_80%", "sulfate_60%_80%", "sulfate_65%_80%" /)
!aerosol_optical_names(330:333)=(/ "sulfate_70%_80%", "sulfate_75%_80%", "sulfate_80%_80%", "sulfate_82%_80%" /)
!aerosol_optical_names(334:337)=(/ "sulfate_84%_80%", "sulfate_86%_80%", "sulfate_88%_80%", "sulfate_90%_80%" /)
!aerosol_optical_names(338:341)=(/ "sulfate_91%_80%", "sulfate_92%_80%", "sulfate_93%_80%", "sulfate_94%_80%" /)
!aerosol_optical_names(342:345)=(/ "sulfate_95%_80%", "sulfate_96%_80%", "sulfate_97%_80%", "sulfate_98%_80%" /)
!aerosol_optical_names(346)=       "sulfate_99%_80%"
!aerosol_optical_names(347)=       "sulfate_100%_80%"
!aerosol_optical_names(348:351)=(/ "sulfate_30%_75%", "sulfate_35%_75%", "sulfate_40%_75%", "sulfate_45%_75%" /)
!aerosol_optical_names(352:355)=(/ "sulfate_50%_75%", "sulfate_55%_75%", "sulfate_60%_75%", "sulfate_65%_75%" /)
!aerosol_optical_names(356:359)=(/ "sulfate_70%_75%", "sulfate_75%_75%", "sulfate_80%_75%", "sulfate_82%_75%" /)
!aerosol_optical_names(360:363)=(/ "sulfate_84%_75%", "sulfate_86%_75%", "sulfate_88%_75%", "sulfate_90%_75%" /)
!aerosol_optical_names(364:367)=(/ "sulfate_91%_75%", "sulfate_92%_75%", "sulfate_93%_75%", "sulfate_94%_75%" /)
!aerosol_optical_names(368:371)=(/ "sulfate_95%_75%", "sulfate_96%_75%", "sulfate_97%_75%", "sulfate_98%_75%" /)
!aerosol_optical_names(372)=       "sulfate_99%_75%"
!aerosol_optical_names(373)=       "sulfate_100%_75%"
!aerosol_optical_names(374:377)=(/ "sulfate_30%_70%", "sulfate_35%_70%", "sulfate_40%_70%", "sulfate_45%_70%" /)
!aerosol_optical_names(378:381)=(/ "sulfate_50%_70%", "sulfate_55%_70%", "sulfate_60%_70%", "sulfate_65%_70%" /)
!aerosol_optical_names(382:385)=(/ "sulfate_70%_70%", "sulfate_75%_70%", "sulfate_80%_70%", "sulfate_82%_70%" /)
!aerosol_optical_names(386:389)=(/ "sulfate_84%_70%", "sulfate_86%_70%", "sulfate_88%_70%", "sulfate_90%_70%" /)
!aerosol_optical_names(390:393)=(/ "sulfate_91%_70%", "sulfate_92%_70%", "sulfate_93%_70%", "sulfate_94%_70%" /)
!aerosol_optical_names(394:397)=(/ "sulfate_95%_70%", "sulfate_96%_70%", "sulfate_97%_70%", "sulfate_98%_70%" /)
!aerosol_optical_names(398)=       "sulfate_99%_70%"
!aerosol_optical_names(399)=       "sulfate_100%_70%"
!aerosol_optical_names(400:403)=(/ "sulfate_30%_65%", "sulfate_35%_65%", "sulfate_40%_65%", "sulfate_45%_65%" /)
!aerosol_optical_names(404:407)=(/ "sulfate_50%_65%", "sulfate_55%_65%", "sulfate_60%_65%", "sulfate_65%_65%" /)
!aerosol_optical_names(408:411)=(/ "sulfate_70%_65%", "sulfate_75%_65%", "sulfate_80%_65%", "sulfate_82%_65%" /)
!aerosol_optical_names(412:415)=(/ "sulfate_84%_65%", "sulfate_86%_65%", "sulfate_88%_65%", "sulfate_90%_65%" /)
!aerosol_optical_names(416:419)=(/ "sulfate_91%_65%", "sulfate_92%_65%", "sulfate_93%_65%", "sulfate_94%_65%" /)
!aerosol_optical_names(420:423)=(/ "sulfate_95%_65%", "sulfate_96%_65%", "sulfate_97%_65%", "sulfate_98%_65%" /)
!aerosol_optical_names(424)=       "sulfate_99%_65%"
!aerosol_optical_names(425)=       "sulfate_100%_65%"
!aerosol_optical_names(426:429)=(/ "sulfate_30%_60%", "sulfate_35%_60%", "sulfate_40%_60%", "sulfate_45%_60%" /)
!aerosol_optical_names(430:433)=(/ "sulfate_50%_60%", "sulfate_55%_60%", "sulfate_60%_60%", "sulfate_65%_60%" /)
!aerosol_optical_names(434:437)=(/ "sulfate_70%_60%", "sulfate_75%_60%", "sulfate_80%_60%", "sulfate_82%_60%" /)
!aerosol_optical_names(438:441)=(/ "sulfate_84%_60%", "sulfate_86%_60%", "sulfate_88%_60%", "sulfate_90%_60%" /)
!aerosol_optical_names(442:445)=(/ "sulfate_91%_60%", "sulfate_92%_60%", "sulfate_93%_60%", "sulfate_94%_60%" /)
!aerosol_optical_names(446:449)=(/ "sulfate_95%_60%", "sulfate_96%_60%", "sulfate_97%_60%", "sulfate_98%_60%" /)
!aerosol_optical_names(450)=       "sulfate_99%_60%"
!aerosol_optical_names(451)=       "sulfate_100%_60%"
!aerosol_optical_names(452:455)=(/ "sulfate_30%_55%", "sulfate_35%_55%", "sulfate_40%_55%", "sulfate_45%_55%" /)
!aerosol_optical_names(456:459)=(/ "sulfate_50%_55%", "sulfate_55%_55%", "sulfate_60%_55%", "sulfate_65%_55%" /)
!aerosol_optical_names(460:463)=(/ "sulfate_70%_55%", "sulfate_75%_55%", "sulfate_80%_55%", "sulfate_82%_55%" /)
!aerosol_optical_names(464:467)=(/ "sulfate_84%_55%", "sulfate_86%_55%", "sulfate_88%_55%", "sulfate_90%_55%" /)
!aerosol_optical_names(468:471)=(/ "sulfate_91%_55%", "sulfate_92%_55%", "sulfate_93%_55%", "sulfate_94%_55%" /)
!aerosol_optical_names(472:475)=(/ "sulfate_95%_55%", "sulfate_96%_55%", "sulfate_97%_55%", "sulfate_98%_55%" /)
!aerosol_optical_names(476)=       "sulfate_99%_55%"
!aerosol_optical_names(477)=       "sulfate_100%_55%"
!aerosol_optical_names(478:481)=(/ "sulfate_30%_50%", "sulfate_35%_50%", "sulfate_40%_50%", "sulfate_45%_50%" /)
!aerosol_optical_names(482:485)=(/ "sulfate_50%_50%", "sulfate_55%_50%", "sulfate_60%_50%", "sulfate_65%_50%" /)
!aerosol_optical_names(486:489)=(/ "sulfate_70%_50%", "sulfate_75%_50%", "sulfate_80%_50%", "sulfate_82%_50%" /)
!aerosol_optical_names(490:493)=(/ "sulfate_84%_50%", "sulfate_86%_50%", "sulfate_88%_50%", "sulfate_90%_50%" /)
!aerosol_optical_names(494:497)=(/ "sulfate_91%_50%", "sulfate_92%_50%", "sulfate_93%_50%", "sulfate_94%_50%" /)
!aerosol_optical_names(498:501)=(/ "sulfate_95%_50%", "sulfate_96%_50%", "sulfate_97%_50%", "sulfate_98%_50%" /)
!aerosol_optical_names(502)=       "sulfate_99%_50%"
!aerosol_optical_names(503)=       "sulfate_100%_50%"
!aerosol_optical_names(504:507)=(/ "sulfate_30%_45%", "sulfate_35%_45%", "sulfate_40%_45%", "sulfate_45%_45%" /)
!aerosol_optical_names(508:511)=(/ "sulfate_50%_45%", "sulfate_55%_45%", "sulfate_60%_45%", "sulfate_65%_45%" /)
!aerosol_optical_names(512:515)=(/ "sulfate_70%_45%", "sulfate_75%_45%", "sulfate_80%_45%", "sulfate_82%_45%" /)
!aerosol_optical_names(516:519)=(/ "sulfate_84%_45%", "sulfate_86%_45%", "sulfate_88%_45%", "sulfate_90%_45%" /)
!aerosol_optical_names(520:523)=(/ "sulfate_91%_45%", "sulfate_92%_45%", "sulfate_93%_45%", "sulfate_94%_45%" /)
!aerosol_optical_names(524:527)=(/ "sulfate_95%_45%", "sulfate_96%_45%", "sulfate_97%_45%", "sulfate_98%_45%" /)
!aerosol_optical_names(528)=       "sulfate_99%_45%"
!aerosol_optical_names(529)=       "sulfate_100%_45%"
!aerosol_optical_names(530:533)=(/ "sulfate_30%_40%", "sulfate_35%_40%", "sulfate_40%_40%", "sulfate_45%_40%" /)
!aerosol_optical_names(534:537)=(/ "sulfate_50%_40%", "sulfate_55%_40%", "sulfate_60%_40%", "sulfate_65%_40%" /)
!aerosol_optical_names(538:541)=(/ "sulfate_70%_40%", "sulfate_75%_40%", "sulfate_80%_40%", "sulfate_82%_40%" /)
!aerosol_optical_names(542:545)=(/ "sulfate_84%_40%", "sulfate_86%_40%", "sulfate_88%_40%", "sulfate_90%_40%" /)
!aerosol_optical_names(546:549)=(/ "sulfate_91%_40%", "sulfate_92%_40%", "sulfate_93%_40%", "sulfate_94%_40%" /)
!aerosol_optical_names(550:553)=(/ "sulfate_95%_40%", "sulfate_96%_40%", "sulfate_97%_40%", "sulfate_98%_40%" /)
!aerosol_optical_names(554)=       "sulfate_99%_40%"
!aerosol_optical_names(555)=       "sulfate_100%_40%"
!aerosol_optical_names(556:559)=(/ "sulfate_30%_35%", "sulfate_35%_35%", "sulfate_40%_35%", "sulfate_45%_35%" /)
!aerosol_optical_names(560:563)=(/ "sulfate_50%_35%", "sulfate_55%_35%", "sulfate_60%_35%", "sulfate_65%_35%" /)
!aerosol_optical_names(564:567)=(/ "sulfate_70%_35%", "sulfate_75%_35%", "sulfate_80%_35%", "sulfate_82%_35%" /)
!aerosol_optical_names(568:571)=(/ "sulfate_84%_35%", "sulfate_86%_35%", "sulfate_88%_35%", "sulfate_90%_35%" /)
!aerosol_optical_names(572:575)=(/ "sulfate_91%_35%", "sulfate_92%_35%", "sulfate_93%_35%", "sulfate_94%_35%" /)
!aerosol_optical_names(576:579)=(/ "sulfate_95%_35%", "sulfate_96%_35%", "sulfate_97%_35%", "sulfate_98%_35%" /)
!aerosol_optical_names(580)=       "sulfate_99%_35%"
!aerosol_optical_names(581)=       "sulfate_100%_35%"
!aerosol_optical_names(582:585)=(/ "sulfate_30%_30%", "sulfate_35%_30%", "sulfate_40%_30%", "sulfate_45%_30%" /)
!aerosol_optical_names(586:589)=(/ "sulfate_50%_30%", "sulfate_55%_30%", "sulfate_60%_30%", "sulfate_65%_30%" /)
!aerosol_optical_names(590:593)=(/ "sulfate_70%_30%", "sulfate_75%_30%", "sulfate_80%_30%", "sulfate_82%_30%" /)
!aerosol_optical_names(594:597)=(/ "sulfate_84%_30%", "sulfate_86%_30%", "sulfate_88%_30%", "sulfate_90%_30%" /)
!aerosol_optical_names(598:601)=(/ "sulfate_91%_30%", "sulfate_92%_30%", "sulfate_93%_30%", "sulfate_94%_30%" /)
!aerosol_optical_names(602:605)=(/ "sulfate_95%_30%", "sulfate_96%_30%", "sulfate_97%_30%", "sulfate_98%_30%" /)
!aerosol_optical_names(606)=       "sulfate_99%_30%"
!aerosol_optical_names(607)=       "sulfate_100%_30%"
!aerosol_optical_names(608:611)=(/ "sulfate_30%_25%", "sulfate_35%_25%", "sulfate_40%_25%", "sulfate_45%_25%" /)
!aerosol_optical_names(612:615)=(/ "sulfate_50%_25%", "sulfate_55%_25%", "sulfate_60%_25%", "sulfate_65%_25%" /)
!aerosol_optical_names(616:619)=(/ "sulfate_70%_25%", "sulfate_75%_25%", "sulfate_80%_25%", "sulfate_82%_25%" /)
!aerosol_optical_names(620:623)=(/ "sulfate_84%_25%", "sulfate_86%_25%", "sulfate_88%_25%", "sulfate_90%_25%" /)
!aerosol_optical_names(624:627)=(/ "sulfate_91%_25%", "sulfate_92%_25%", "sulfate_93%_25%", "sulfate_94%_25%" /)
!aerosol_optical_names(628:631)=(/ "sulfate_95%_25%", "sulfate_96%_25%", "sulfate_97%_25%", "sulfate_98%_25%" /)
!aerosol_optical_names(632)=       "sulfate_99%_25%"
!aerosol_optical_names(633)=       "sulfate_100%_25%"
!aerosol_optical_names(634:637)=(/ "sulfate_30%_20%", "sulfate_35%_20%", "sulfate_40%_20%", "sulfate_45%_20%" /)
!aerosol_optical_names(638:641)=(/ "sulfate_50%_20%", "sulfate_55%_20%", "sulfate_60%_20%", "sulfate_65%_20%" /)
!aerosol_optical_names(642:645)=(/ "sulfate_70%_20%", "sulfate_75%_20%", "sulfate_80%_20%", "sulfate_82%_20%" /)
!aerosol_optical_names(646:649)=(/ "sulfate_84%_20%", "sulfate_86%_20%", "sulfate_88%_20%", "sulfate_90%_20%" /)
!aerosol_optical_names(650:653)=(/ "sulfate_91%_20%", "sulfate_92%_20%", "sulfate_93%_20%", "sulfate_94%_20%" /)
!aerosol_optical_names(654:657)=(/ "sulfate_95%_20%", "sulfate_96%_20%", "sulfate_97%_20%", "sulfate_98%_20%" /)
!aerosol_optical_names(658)=       "sulfate_99%_20%"
!aerosol_optical_names(659)=       "sulfate_100%_20%"
!aerosol_optical_names(660:663)=(/ "sulfate_30%_15%", "sulfate_35%_15%", "sulfate_40%_15%", "sulfate_45%_15%" /)
!aerosol_optical_names(664:667)=(/ "sulfate_50%_15%", "sulfate_55%_15%", "sulfate_60%_15%", "sulfate_65%_15%" /)
!aerosol_optical_names(668:671)=(/ "sulfate_70%_15%", "sulfate_75%_15%", "sulfate_80%_15%", "sulfate_82%_15%" /)
!aerosol_optical_names(672:675)=(/ "sulfate_84%_15%", "sulfate_86%_15%", "sulfate_88%_15%", "sulfate_90%_15%" /)
!aerosol_optical_names(676:679)=(/ "sulfate_91%_15%", "sulfate_92%_15%", "sulfate_93%_15%", "sulfate_94%_15%" /)
!aerosol_optical_names(680:683)=(/ "sulfate_95%_15%", "sulfate_96%_15%", "sulfate_97%_15%", "sulfate_98%_15%" /)
!aerosol_optical_names(684)=       "sulfate_99%_15%"
!aerosol_optical_names(685)=       "sulfate_100%_15%"
!aerosol_optical_names(686:689)=(/ "sulfate_30%_10%", "sulfate_35%_10%", "sulfate_40%_10%", "sulfate_45%_10%" /)
!aerosol_optical_names(690:693)=(/ "sulfate_50%_10%", "sulfate_55%_10%", "sulfate_60%_10%", "sulfate_65%_10%" /)
!aerosol_optical_names(694:697)=(/ "sulfate_70%_10%", "sulfate_75%_10%", "sulfate_80%_10%", "sulfate_82%_10%" /)
!aerosol_optical_names(698:701)=(/ "sulfate_84%_10%", "sulfate_86%_10%", "sulfate_88%_10%", "sulfate_90%_10%" /)
!aerosol_optical_names(702:705)=(/ "sulfate_91%_10%", "sulfate_92%_10%", "sulfate_93%_10%", "sulfate_94%_10%" /)
!aerosol_optical_names(706:709)=(/ "sulfate_95%_10%", "sulfate_96%_10%", "sulfate_97%_10%", "sulfate_98%_10%" /)
!aerosol_optical_names(710)=       "sulfate_99%_10%"
!aerosol_optical_names(711)=       "sulfate_100%_10%"
!aerosol_optical_names(712:715)=(/ "sulfate_30%_5%", "sulfate_35%_5%", "sulfate_40%_5%", "sulfate_45%_5%" /)
!aerosol_optical_names(716:719)=(/ "sulfate_50%_5%", "sulfate_55%_5%", "sulfate_60%_5%", "sulfate_65%_5%" /)
!aerosol_optical_names(720:723)=(/ "sulfate_70%_5%", "sulfate_75%_5%", "sulfate_80%_5%", "sulfate_82%_5%" /)
!aerosol_optical_names(724:727)=(/ "sulfate_84%_5%", "sulfate_86%_5%", "sulfate_88%_5%", "sulfate_90%_5%" /)
!aerosol_optical_names(728:731)=(/ "sulfate_91%_5%", "sulfate_92%_5%", "sulfate_93%_5%", "sulfate_94%_5%" /)
!aerosol_optical_names(732:735)=(/ "sulfate_95%_5%", "sulfate_96%_5%", "sulfate_97%_5%", "sulfate_98%_5%" /)
!aerosol_optical_names(736)=       "sulfate_99%_5%"
!aerosol_optical_names(737)=       "sulfate_100%_5%"
!aerosol_optical_names(738:741)=(/ "sulfate_30%_0%", "sulfate_35%_0%", "sulfate_40%_0%", "sulfate_45%_0%" /)
!aerosol_optical_names(742:745)=(/ "sulfate_50%_0%", "sulfate_55%_0%", "sulfate_60%_0%", "sulfate_65%_0%" /)
!aerosol_optical_names(746:751)=(/ "sulfate_70%_0%", "sulfate_75%_0%", "sulfate_80%_0%", "sulfate_82%_0%" /)
!aerosol_optical_names(750:753)=(/ "sulfate_84%_0%", "sulfate_86%_0%", "sulfate_88%_0%", "sulfate_90%_0%" /)
!aerosol_optical_names(754:757)=(/ "sulfate_91%_0%", "sulfate_92%_0%", "sulfate_93%_0%", "sulfate_94%_0%" /)
!aerosol_optical_names(758:761)=(/ "sulfate_95%_0%", "sulfate_96%_0%", "sulfate_97%_0%", "sulfate_98%_0%" /)
!aerosol_optical_names(762)=       "sulfate_99%_0%"
!aerosol_optical_names(763)=       "sulfate_100%_0%"
!aerosol_optical_names(766:769)=(/ "dust1",       "dust2",       "dust3",       "dust4"    /)
!aerosol_optical_names(770)=       "dust5"
!aerosol_optical_names(771)=       "bcdry"
aerosol_optical_names(naermodels-5)="organic_carbon"
aerosol_optical_names(naermodels-4)="soot"
aerosol_optical_names(naermodels-3)="dust_0.73"
aerosol_optical_names(naermodels-2)="dust_1.4"
aerosol_optical_names(naermodels-1)="dust_2.4"
aerosol_optical_names(naermodels)  ="dust_4.5"

! shortwave optical models

if ( myid==0 ) then
  filename = trim(cnsdir)//'/Ginoux_Reddy_2005'
  unit = 16
  open(unit,file=filename,iostat=ierr,status='old')
  if ( ierr/=0 ) then
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
  allocate( endaerwvnsf(num_wavenumbers) )
  allocate( aeroextivl  (num_wavenumbers, naermodels) )
  allocate( aerossalbivl(num_wavenumbers, naermodels) )
  allocate( aeroasymmivl(num_wavenumbers, naermodels) )
  allocate( aeroext_in  (num_wavenumbers ) )
  allocate( aerossalb_in(num_wavenumbers ) )
  allocate( aeroasymm_in(num_wavenumbers ) )
  allocate( nivl1aero (Solar_spect%nbands) )
  allocate( nivl2aero (Solar_spect%nbands) )
  allocate( solivlaero(Solar_spect%nbands, num_wavenumbers) )
  allocate( sflwwts(N_AEROSOL_BANDS, num_wavenumbers) )
  allocate( sflwwts_cn(N_AEROSOL_BANDS_CN, num_wavenumbers) )           
  allocate( dum_aero(num_wavenumbers, 3*naermodels) )

  read(unit,* )
  read(unit,* ) endaerwvnsf
 
  !----------------------------------------------------------------------
  !    match the names of optical property categories from input file with
  !    those specified in the namelist, and store the following data
  !    appropriately.
  !----------------------------------------------------------------------
  optical_check(:) = .false.
  do n = 1,num_input_categories
    read( unit,* ) name_in
    read( unit,* )
    read( unit,* ) aeroext_in
    read( unit,* )
    read( unit,* ) aerossalb_in
    read( unit,* )
    read( unit,* ) aeroasymm_in
    do noptical = 1,naermodels-4
      if ( trim(aerosol_optical_names(noptical))==trim(name_in) ) then
        write(6,*) "Loading optical model for ",trim(name_in)
        optical_check(noptical) = .true.
        aeroextivl(:,noptical)   = aeroext_in
        aerossalbivl(:,noptical) = aerossalb_in
        aeroasymmivl(:,noptical) = aeroasymm_in
        exit
      endif
    end do
  end do
  if ( .not.all(optical_check) ) then
    write(6,*) "ERROR: Cannot idenify optical model"
    do noptical = 1,naermodels-4
      if ( .not.optical_check(noptical) ) then
        write(6,*) "ERROR: Cannot idenify ",trim(aerosol_optical_names(noptical))
      end if    
    end do
    call ccmpi_abort(-1)
  end if
  ! Dust_0.73
  frac = (real(dustreff(1),8) - 0.4E-6_8)/0.4E-6_8
  aeroextivl(:,naermodels-3)   = (1.-frac)*aeroextivl(:,6)   + frac*aeroextivl(:,7)
  aerossalbivl(:,naermodels-3) = (1.-frac)*aerossalbivl(:,6) + frac*aerossalbivl(:,7)
  aeroasymmivl(:,naermodels-3) = (1.-frac)*aeroasymmivl(:,6) + frac*aeroasymmivl(:,7)
  ! Dust 1.4
  frac = (real(dustreff(2),8) - 1.E-6_8)/1.E-6_8
  aeroextivl(:,naermodels-2)   = (1.-frac)*aeroextivl(:,8)   + frac*aeroextivl(:,9)
  aerossalbivl(:,naermodels-2) = (1.-frac)*aerossalbivl(:,8) + frac*aerossalbivl(:,9)
  aeroasymmivl(:,naermodels-2) = (1.-frac)*aeroasymmivl(:,8) + frac*aeroasymmivl(:,9)
  ! Dust 2.4
  frac = (real(dustreff(3),8) - 2.E-6_8)/2.E-6_8
  aeroextivl(:,naermodels-1)   = (1.-frac)*aeroextivl(:,9)   + frac*aeroextivl(:,10)
  aerossalbivl(:,naermodels-1) = (1.-frac)*aerossalbivl(:,9) + frac*aerossalbivl(:,10)
  aeroasymmivl(:,naermodels-1) = (1.-frac)*aeroasymmivl(:,9) + frac*aeroasymmivl(:,10)
  ! Dust 4.5
  frac = (real(dustreff(4),8) - 4.E-6_8)/4.E-6_8
  aeroextivl(:,naermodels)   = (1.-frac)*aeroextivl(:,10)   + frac*aeroextivl(:,11)
  aerossalbivl(:,naermodels) = (1.-frac)*aerossalbivl(:,10) + frac*aerossalbivl(:,11)
  aeroasymmivl(:,naermodels) = (1.-frac)*aeroasymmivl(:,10) + frac*aeroasymmivl(:,11)  

  close(unit)
  deallocate( aeroasymm_in, aerossalb_in, aeroext_in )
  
  dum_aero(:,1:naermodels)                = aeroextivl(:,:)
  dum_aero(:,naermodels+1:2*naermodels)   = aerossalbivl(:,:)
  dum_aero(:,2*naermodels+1:3*naermodels) = aeroasymmivl(:,:)

else
  call ccmpi_bcast(num_wavenumbers,0,comm_world)
  allocate ( endaerwvnsf(num_wavenumbers) )
  allocate ( aeroextivl  (num_wavenumbers, naermodels) )
  allocate ( aerossalbivl(num_wavenumbers, naermodels) )
  allocate ( aeroasymmivl(num_wavenumbers, naermodels) )
  allocate ( nivl1aero (Solar_spect%nbands) )
  allocate ( nivl2aero (Solar_spect%nbands) )
  allocate ( solivlaero(Solar_spect%nbands, num_wavenumbers) )
  allocate ( sflwwts(N_AEROSOL_BANDS, num_wavenumbers) )
  allocate ( sflwwts_cn(N_AEROSOL_BANDS_CN, num_wavenumbers) )  
  allocate ( dum_aero(num_wavenumbers, 3*naermodels) )
end if

call ccmpi_bcast(endaerwvnsf,0,comm_world)
call ccmpi_bcastr8(dum_aero,0,comm_world)

aeroextivl(:,:)   = dum_aero(:,1:naermodels)
aerossalbivl(:,:) = dum_aero(:,naermodels+1:2*naermodels)
aeroasymmivl(:,:) = dum_aero(:,2*naermodels+1:3*naermodels)
deallocate( dum_aero )

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
  if ( nw==endaerwvnsf(nivl3) ) then
    solivlaero(nband,nivl3) = sumsol3
    sumsol3 = 0.0_8
  end if
  if ( nw==Solar_spect%endwvnbands(nband) ) then
    if ( nw/=endaerwvnsf(nivl3) ) then
      solivlaero(nband,nivl3) = sumsol3 
      sumsol3 = 0.0_8
    end if
    nivl2aero(nband) = nivl3
    nband = nband + 1
    if ( nband<=Solar_spect%nbands ) then
      if ( nw==endaerwvnsf(nivl3) ) then
        nivl1aero(nband) = nivl3 + 1
      else
        nivl1aero(nband) = nivl3
      end if
    end if
  end if
  if ( nw==endaerwvnsf(nivl3) ) nivl3 = nivl3 + 1
end do

Aerosol_props%aerextband   = 0._8
Aerosol_props%aerssalbband = 0._8
Aerosol_props%aerasymmband = 0._8

do nmodel = 1,naermodels
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
  
Aerosol_props%aerextbandlw      = 0._8
Aerosol_props%aerssalbbandlw    = 0._8
Aerosol_props%aerextbandlw_cn   = 0._8
Aerosol_props%aerssalbbandlw_cn = 0._8

do nw=1,naermodels    
  do na=1,N_AEROSOL_BANDS  
    do ni=1,num_wavenumbers 
      Aerosol_props%aerextbandlw(na,nw) =                   &
                      Aerosol_props%aerextbandlw(na,nw) +   &
                      aeroextivl(ni,nw)*sflwwts(na,ni)*     &
                      1.0E+03_8
      Aerosol_props%aerssalbbandlw(na,nw) =                 &
                      Aerosol_props%aerssalbbandlw(na,nw) + &
                      aerossalbivl(ni,nw)*sflwwts(na,ni)
    end do
  end do
end do
do nw=1,naermodels    
  do na=1,N_AEROSOL_BANDS_CN
    do ni=1,num_wavenumbers 
      Aerosol_props%aerextbandlw_cn(na,nw) =                   &
                      Aerosol_props%aerextbandlw_cn(na,nw) +   &
                      aeroextivl(ni,nw)*sflwwts_cn(na,ni)*     &
                      1.0E+03_8
      Aerosol_props%aerssalbbandlw_cn(na,nw) =                 &
                      Aerosol_props%aerssalbbandlw_cn(na,nw) + &
                      aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
    end do
  end do
end do

deallocate ( sflwwts_cn, sflwwts )
deallocate ( solivlaero, nivl2aero, nivl1aero ) 
deallocate ( aeroasymmivl, aerossalbivl, aeroextivl )
deallocate ( endaerwvnsf )

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
!      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandlo_cn =  &
!      (/ 800.0_8 /)
!
!      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandhi_cn =  &
!      (/ 1200.0_8 /)

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
      real(kind=8),    dimension (N_AEROSOL_BANDS)     :: planckaerband
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN)  :: planckaerband_cn

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
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
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real(kind=8),    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real(kind=8),    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN, num_wavenumbers) :: &
                                                  planckivlaer_cn   
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

     nivl2aer_cn=0
     nivl2aer_co=0
     
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
