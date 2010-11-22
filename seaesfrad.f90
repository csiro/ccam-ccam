! Interface for SEA-ESF radiation scheme from GFDL with CCAM.

module seaesfrad_m

use rad_utilities_mod, only: atmos_input_type,surface_type,astronomy_type,aerosol_type, &
                             aerosol_properties_type,radiative_gases_type,cldrad_properties_type, &
                             cld_specification_type,lw_output_type,sw_output_type, &
                             aerosol_diagnostics_type,time_type,microphysics_type, &
                             microrad_properties_type,lw_diagnostics_type,lw_table_type, &
                             Sw_control,Lw_control, Rad_control,Cldrad_control,Lw_parameters
use esfsw_driver_mod, only : swresf,esfsw_driver_init
use sealw99_mod, only : sealw99,sealw99_init
use esfsw_parameters_mod, only:  Solar_spect,esfsw_parameters_init

private
public seaesfrad

real, parameter :: cp       = 1004.64    ! Specific heat of dry air at const P
real, parameter :: grav     = 9.80616    ! Acceleration of gravity
real, parameter :: stefbo   = 5.67e-8    ! Stefan-Boltzmann constant
real, parameter :: rdry     = 287.04     ! gas constant for dry air
real, parameter :: rhow     = 1000.      ! Density of water
real, parameter :: pi       = 3.1415927  ! pi
real, parameter :: csolar   = 1365       ! Solar constant in W/m^2
real, parameter :: siglow   =.68         ! sigma level for top of low cloud
real, parameter :: sigmid   =.44         ! sigma level for top of medium cloud
real, parameter :: ratco2mw =1.519449738
real, parameter :: cong     = cp/grav
logical, parameter :: do_totcld_forcing             = .true.
logical, parameter :: calculate_volcanic_sw_heating = .false.

logical, save :: do_aerosol_forcing

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CCAM interface
!

subroutine seaesfrad(odcalc,iaero)

use zenith_m
use microphys_rad_mod, only: microphys_sw_driver,microphys_lw_driver,lwemiss_calc,microphys_rad_init

use ateb
use cc_mpi
use cable_ccam, only : CABLE
use latlong_m
use mlo
use ozoneread

implicit none

include 'parm.h'
include 'newmpar.h'
integer, parameter :: imax=il*nrows_rad
include 'arrays.h'
include 'cparams.h'
include 'dates.h'
include 'extraout.h'
include 'kuocom.h'
include 'liqwpar.h'
include 'nsibd.h'
include 'pbl.h'
include 'raddiag.h'
include 'sigs.h'
include 'soil.h'
include 'soilsnow.h'

logical, intent(in) :: odcalc  ! True for full radiation calculation
integer, intent(in) :: iaero
integer, dimension(12) :: ndoy   ! days from beginning of year (1st Jan is 0)
integer jyear,jmonth,jday,jhour,jmin
integer k,ksigtop,mstart,mins
integer i,j,iq,istart,iend,kr
integer, save :: nlow,nmid
real, dimension(ifull), save :: sgamp
real, dimension(ifull,kl), save :: rtt
real, dimension(imax) :: qsat,coszro2,taudar2,coszro,taudar
real, dimension(imax) :: sg,sint,sout,sgdn,rg,rt,rgdn
real, dimension(imax) :: soutclr,sgclr,rtclr,rgclr,sga
real, dimension(imax) :: sgvis,sgdnvisdir,sgdnvisdif,sgdnnirdir,sgdnnirdif
real, dimension(imax,kl) :: duo3n
real, dimension(imax) :: cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif
real, dimension(imax,kl) :: p2,cd2
real, dimension(kl+1) :: sigh
real(kind=8), dimension(kl+1,2) :: pref
real dduo3n,ddo3n2,ddo3n3,ddo3n4
real qccon,qlrad,qfrad,cfrac
real r1,dlt,alp,slag
real dhr,fjd,bpyear
real ttbg,ar1,exp_ar1,ar2,exp_ar2,ar3,snr
real dnsnow,snrat,dtau,alvo,aliro,fage,cczen,fzen,fzenm
real alvd,alv,alird,alir
real rrvco2,rrvch4,rrvn2o,rrvf11,rrvf12,rrvf113,rrvf22,ssolar,rrco2
real f1,f2,cosz,delta
logical maxover,newcld
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
real(kind=8), dimension(:,:,:,:), allocatable :: r

common/cfrac/cfrac(ifull,kl)
common/work3f/qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl) ! ditto
common /radisw2/ rrco2, ssolar, rrvco2,rrvch4,rrvn2o,rrvf11,rrvf12,rrvf113,rrvf22
common /o3dat/ dduo3n(37,kl),ddo3n2(37,kl),ddo3n3(37,kl),ddo3n4(37,kl)

data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

! Aerosol flag
do_aerosol_forcing=abs(iaero).gt.1

! set-up half levels ------------------------------------------------
sigh(1:kl) = sigmh(1:kl)
sigh(kl+1) = 0.

! set-up standard pressure levels -----------------------------------
pref(kl+1,1)=101325.
pref(kl+1,2)=81060.
do k=1,kl
  kr=kl+1-k
  pref(kr,:)=sig(k)*pref(kl+1,:)
end do

! astronomy ---------------------------------------------------------
! Set up number of minutes from beginning of year
! For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
jyear=kdate/10000
jmonth=(kdate-jyear*10000)/100
jday=kdate-jyear*10000-jmonth*100
jhour=ktime/100
jmin=ktime-jhour*100
mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y
! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

! Set number of years before present for orbital parameters.
! Allowed values are 0, 6000 and 21000.
bpyear = 0.
if(nhstest<0)then  ! aquaplanet test
  fjd = 79.+mod(mins,1440)/1440.       ! set to 21 March +frac of day
else
  fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
endif

! Calculate sun position
call solargh(fjd,bpyear,r1,dlt,alp,slag)
ssolar = csolar / (r1*r1)

! Initialisation ----------------------------------------------------
if ( first ) then
  first = .false.

  ! initialise co2
  call co2_read(sig,jyear)
  rrco2=rrvco2*ratco2mw

  ! initialise ozone
  if(amipo3)then
    write(6,*) 'AMIP2 ozone input'
    call o3read_amip
  else
    call o3_read(sig,jyear,jmonth)
    call resetd(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*kl)
  end if
  
  Cldrad_control%do_strat_clouds_iz      =.true.
  Cldrad_control%do_sw_micro_iz          =.true.
  Cldrad_control%do_lw_micro_iz          =.true.
  Cldrad_control%do_sw_micro             =.true.
  Cldrad_control%do_lw_micro             =.true.
  Cldrad_control%do_ica_calcs            =.false. ! must change allocations below if true
  Cldrad_control%do_no_clouds            =.false.
  Cldrad_control%do_donner_deep_clouds   =.false.
  Cldrad_control%do_stochastic_clouds    =.false.  
  Sw_control%solar_constant              =csolar
  Sw_control%do_cmip_diagnostics         =.false.
  Lw_control%do_lwcldemiss               =.true.
  Lw_control%do_o3_iz                    =.true.
  Lw_control%do_co2_iz                   =.true.
  Lw_control%do_ch4_iz                   =.true.
  Lw_control%do_n2o_iz                   =.true.
  Lw_control%do_o3                       =.true.
  Lw_control%do_co2                      =.true.
  Lw_control%do_ch4                      =rrvch4.gt.0.
  Lw_control%do_n2o                      =rrvch4.gt.0.
  Lw_control%do_h2o                      =.true.
  Lw_control%do_cfc                      =rrvch4.gt.0.
  Rad_control%using_solar_timeseries_data=.false.
  Rad_control%do_totcld_forcing          =do_totcld_forcing
  Rad_control%rad_time_step              =kountr*dt
  Rad_control%rad_time_step_iz           =.true.
  Rad_control%do_aerosol                 =.false.
  Rad_control%do_swaerosol_forcing       =.false.
  Rad_control%do_lwaerosol_forcing       =.false.
  Astro%rrsun                            =1./(r1*r1)

  call sealw99_init(pref, Lw_tables)
  call esfsw_parameters_init
  call esfsw_driver_init
  call microphys_rad_init

  allocate ( Atmos_input%press(imax, 1, kl+1) )
  allocate ( Atmos_input%phalf(imax, 1, kl+1) )
  allocate ( Atmos_input%temp (imax, 1, kl+1) )
  allocate ( Atmos_input%rh2o (imax, 1, kl  ) )
  allocate ( Atmos_input%rel_hum(imax, 1, kl  ) )
  !allocate ( Atmos_input%cloudtemp(imax, 1,kl  ) )
  !allocate ( Atmos_input%cloudvapor(imax, 1,kl  ) )
  allocate ( Atmos_input%clouddeltaz(imax, 1,kl  ) )
  !allocate ( Atmos_input%aerosoltemp(imax, 1,kl  ) )
  !allocate ( Atmos_input%aerosolpress(imax, 1,kl+1) )
  !allocate ( Atmos_input%aerosolvapor(imax, 1,kl  ) )
  !allocate ( Atmos_input%aerosolrelhum(imax, 1,kl  ) )
  allocate ( Atmos_input%deltaz(imax, 1, kl ) )
  allocate ( Atmos_input%pflux (imax, 1, kl+1) )
  allocate ( Atmos_input%tflux (imax, 1, kl+1) )
  allocate ( Atmos_input%psfc (imax, 1 ) )
  allocate ( Atmos_input%tsfc (imax, 1 ) )
  !if (use_co2_tracer_field) then
  !  allocate ( Atmos_input%tracer_co2(imax, 1, kl ) )
  !endif
  
  allocate(Rad_gases%qo3(imax,1,kl))

  allocate(Cloud_microphysics%size_rain(imax, 1, kl))
  allocate(Cloud_microphysics%size_drop(imax, 1, kl))
  allocate(Cloud_microphysics%size_ice (imax, 1, kl))
  allocate(Cloud_microphysics%size_snow (imax, 1, kl))
  allocate(Cloud_microphysics%conc_drop(imax, 1, kl))
  allocate(Cloud_microphysics%conc_ice (imax, 1, kl))
  allocate(Cloud_microphysics%conc_rain(imax, 1, kl))
  allocate(Cloud_microphysics%conc_snow(imax, 1, kl))

  allocate (Cldrad_props%cldext  (imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%cldsct  (imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%cldemiss(imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%emmxolw (imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%emrndlw (imax, 1, kl, Cldrad_control%nlwcldb,1))

  allocate (Lscrad_props%cldext(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%cldsct(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb) )

  allocate ( Cld_spec%camtsw (imax, 1, kl ) )
  allocate ( Cld_spec%cmxolw (imax, 1, kl ) )
  allocate ( Cld_spec%crndlw (imax, 1, kl ) )
  
  allocate (Surface%asfc_vis_dir (imax, 1 ) )
  allocate (Surface%asfc_nir_dir (imax, 1 ) )
  allocate (Surface%asfc_vis_dif (imax, 1 ) )
  allocate (Surface%asfc_nir_dif (imax, 1 ) )

  allocate ( Astro%cosz   (imax, 1 ) )
  allocate ( Astro%fracday(imax, 1 ) )

  allocate (Lw_output(1)%heatra(imax,1,kl)  )
  allocate (Lw_output(1)%flxnet(imax,1,kl+1))
  allocate (Lw_output(1)%bdy_flx(imax,1,4)  )
  if (do_totcld_forcing) then
    allocate (Lw_output(1)%heatracf(imax,1,kl)  )
    allocate (Lw_output(1)%flxnetcf(imax,1,kl+1))
    allocate (Lw_output(1)%bdy_flx_clr(imax,1,4))
  endif

  allocate (Sw_output(1)%dfsw(imax,1,kl+1))
  allocate (Sw_output(1)%ufsw(imax,1,kl+1))
  allocate (Sw_output(1)%dfsw_dir_sfc(imax,1) )
  allocate (Sw_output(1)%dfsw_dif_sfc(imax,1) )
  allocate (Sw_output(1)%ufsw_dif_sfc(imax,1) )
  allocate (Sw_output(1)%fsw(imax,1,kl+1)     )
  allocate (Sw_output(1)%hsw(imax,1,kl)     )
  allocate (Sw_output(1)%dfsw_vis_sfc(imax,1)    )
  allocate (Sw_output(1)%ufsw_vis_sfc(imax,1)    )
  allocate (Sw_output(1)%dfsw_vis_sfc_dir(imax,1)    )
  allocate (Sw_output(1)%dfsw_vis_sfc_dif(imax,1)    )
  allocate (Sw_output(1)%ufsw_vis_sfc_dif(imax,1)    )
  allocate (Sw_output(1)%bdy_flx(imax,1,4))
  if (do_totcld_forcing) then
    allocate (Sw_output(1)%dfswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%ufswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%fswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%hswcf(imax,1,kl)   )
    allocate (Sw_output(1)%dfsw_dir_sfc_clr(imax,1))
    allocate (Sw_output(1)%dfsw_dif_sfc_clr(imax,1))
    allocate (Sw_output(1)%bdy_flx_clr(imax,1,4))
  endif

  ! define diagnostic cloud levels
  f1=1.
  f2=1.
  do k=1,kl-1
    if(abs(sigh(k+1)-siglow)<f1)then
      f1=abs(sigh(k+1)-siglow)
      nlow=k
    endif
    if(abs(sigh(k+1)-sigmid)<f2)then
      f2=abs(sigh(k+1)-sigmid)
      nmid=k
    endif
  enddo

  ! initialise VIS fraction of SW radiation
  swrsave=0.5
  
end if  ! (first)

! Prepare SEA-ESF arrays --------------------------------------------
Rad_time%days   =fjd
Rad_time%seconds=60*(60*jhour + jmin)
Rad_time%ticks  =0

if (ldr.eq.0) then
  write(6,*) "ERROR: SEA-ESF radiation requires ldr.ne.0"
  stop
end if

! main loop ---------------------------------------------------------
if(mod(ifull,imax).ne.0)then
  write(6,*) 'nproc,il,jl,ifull,imax ',nproc,il,jl,ifull,imax
  stop 'illegal setting of imax in rdparm'
endif

do j=1,jl,imax/il
  istart=1+(j-1)*il
  iend=istart+imax-1

  ! Calculate zenith angle for the solarfit calculation.
  ! This call averages zenith angle just over this time step.
  dhr = dt/3600.
  call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro2,taudar2)
  call atebccangle(istart,imax,coszro2(1:imax),rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt,sin(dlt))

  ! Call radiation --------------------------------------------------
  if ( odcalc ) then     ! Do the calculation

    ! Average the zenith angle over the time (hours) between radiation
    ! calculations
    dhr = kountr*dt/3600.0
    call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro,taudar)
    
    ! Set up ozone for this time and row
    if (amipo3) then
      call o3set_amip ( rlatt(istart:iend), imax, mins,sigh, ps(istart:iend), Rad_gases%qo3(:,1,:) )
      Rad_gases%qo3(:,1,:)=max(1.e-10,Rad_gases%qo3(:,1,:))    ! July 2008
    else
      call o3set(rlatt(istart:iend),rlongg(istart:iend),imax,mins,duo3n,sig,ps(istart:iend))
      ! Conversion of o3 from units of cm stp to gm/gm
      do k=1,kl
        Rad_gases%qo3(:,1,k) = max(1.e-10,duo3n(1:imax,k))
      end do
    end if

    ! Set-up albedo
    ! Land albedo ---------------------------------------------------
    if (nsib.eq.CABLE.or.nsib.eq.6.or.nsib.eq.7) then
      ! CABLE version
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albvisdir(istart:iend)    ! from cable (inc snow)
        cirrf_dir(1:imax) = albnirdir(istart:iend)    ! from cable (inc snow)
        cuvrf_dif(1:imax) = albvisdif(istart:iend)    ! from cable (inc snow)
        cirrf_dif(1:imax) = albnirdif(istart:iend)    ! from cable (inc snow)
      end where
    else
      ! nsib=3 version (calculate snow)
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albsav(istart:iend)    ! from albfile (indata.f)
        cirrf_dir(1:imax) = albnirsav(istart:iend) ! from albnirfile (indata.f)
        cuvrf_dif(1:imax) = cuvrf_dir(1:imax)      ! assume DIR and DIF are the same
        cirrf_dif(1:imax) = cirrf_dir(1:imax)      ! assume DIR and DIF are the same
      end where
      ! The following snow calculation show be done by sib3 (sflux.f)
      do i=1,imax
        iq=i+(j-1)*il
        if (land(iq)) then
          if (snowd(iq).gt.0.) then
            ttbg=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
            ttbg=min(ttbg,273.1)
            ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
            exp_ar1=exp(ar1)                  ! e.g. exp(0 to -4)
            ar2 = 10.*ar1                     ! freezing of melt water
            exp_ar2=exp(ar2)                  ! e.g. exp(0 to -40)
            snr=snowd(iq)/max(ssdnn(iq),100.)
            if(isoilm(iq).eq.9)then   ! fixes for Arctic & Antarctic
              ar3=.001
              dnsnow=max(dnsnow,.0015)
              snrat=min(1.,snr/(snr+.001))
            else
              ar3=.3               ! accumulation of dirt
              snrat=min(1.,snr/(snr+.02))
            endif
            dtau=1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
            if(snowd(iq).le. 1.)then
              snage(iq)=0.
            else
              snage(iq)=max(0.,(snage(iq) + dtau)*(1.-dnsnow))
            endif
            alvo = 0.95	        !alb. for vis. on a new snow
            aliro = 0.65        !alb. for near-infr. on a new snow
            fage = 1.-1./(1.+snage(iq))	 !age factor
            cczen=max(.17365, coszro(i))
            fzen=( 1.+1./2.)/(1.+2.*2.*cczen) -1./2.
            if( cczen .gt. 0.5 ) fzen = 0.
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
    ! NCAR CCMS3.0 scheme (Briegleb et al, 1986,
    ! J. Clim. and Appl. Met., v. 27, 214-226)     ! 
    where (.not.land(istart:iend).and.coszro(1:imax).ge.0.)
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
    
    ! MLO albedo ----------------------------------------------------
    call mloalb4(istart,imax,coszro,cuvrf_dir,cuvrf_dif,cirrf_dir,cirrf_dif,0)    

    ! Urban albedo --------------------------------------------------
    call atebalb1(istart,imax,cuvrf_dir(1:imax),0)
    call atebalb1(istart,imax,cirrf_dir(1:imax),0)
    call atebalb1(istart,imax,cuvrf_dif(1:imax),0)
    call atebalb1(istart,imax,cirrf_dif(1:imax),0)    

    ! Aerosols -------------------------------------------------------
    select case (abs(iaero))
      case(0)
        ! no aerosols
      case(1)
        ! aerosols are read in (direct effect only)
        do i=1,imax
          iq=i+(j-1)*il
          cosz = max ( coszro(i), 1.e-4)
          delta =  coszro(i)*0.29*8.*so4t(iq)*((1.-0.25*(cuvrf_dir(i)+cuvrf_dif(i)+cirrf_dir(i)+cirrf_dif(i)))/cosz)**2
          cuvrf_dir(i)=min(0.99, delta+cuvrf_dir(i)) ! still broadband
          cirrf_dir(i)=min(0.99, delta+cirrf_dir(i)) ! still broadband
          cuvrf_dif(i)=min(0.99, delta+cuvrf_dif(i)) ! still broadband
          cirrf_dif(i)=min(0.99, delta+cirrf_dif(i)) ! still broadband
        end do ! i=1,imax
      case(2)
        ! prognostic aerosols
        !Aerosol=
        !Aerosol_props=
        write(6,*) "ERROR: prognostic aerosols are not supported"
        stop
      case DEFAULT
        write(6,*) "ERROR: unknown iaero option ",iaero
        stop
    end select

    ! define droplet size (from radriv90.f) -------------------------
    if (iaero.ne.2) then
      where (land(istart:iend).and.rlatt(istart:iend)>0.)
        cd2(1:imax,1)=cdropl_nh
      else where (land(istart:iend))
        cd2(1:imax,1)=cdropl_sh
      else where (rlatt(istart:iend)>0.)
        cd2(1:imax,1)=cdrops_nh
      else where
        cd2(1:imax,1)=cdrops_sh
      end where
      do k=2,kl
        cd2(1:imax,k)=cd2(1:imax,1)
      enddo
    else
      !call cldrop(cd2,rhoa)
    end if
    
    ! Cloud fraction diagnostics ------------------------------------
    cloudlo(istart:iend)=0.
    cloudmi(istart:iend)=0.
    cloudhi(istart:iend)=0.
    ! Diagnose low, middle and high clouds
    do k=1,nlow
      cloudlo(istart:iend)=cloudlo(istart:iend)+cfrac(istart:iend,k)-cloudlo(istart:iend)*cfrac(istart:iend,k)
    enddo
    do k=nlow+1,nmid
      cloudmi(istart:iend)=cloudmi(istart:iend)+cfrac(istart:iend,k)-cloudmi(istart:iend)*cfrac(istart:iend,k)
    enddo
    do k=nmid+1,kl-1
      cloudhi(istart:iend)=cloudhi(istart:iend)+cfrac(istart:iend,k)-cloudhi(istart:iend)*cfrac(istart:iend,k)
    enddo

    ! Prepare SEA-ESF arrays ----------------------------------------
    do k=1,kl
      kr=kl+1-k
      Atmos_input%deltaz(:,1,kr) =(-dsig(k)/sig(k))*rdry*t(istart:iend,k)/grav
      Atmos_input%rh2o(:,1,kr)   =max(qg(istart:iend,k) ,2.e-7)
      Atmos_input%temp(:,1,kr)   =t(istart:iend,k)      
      Atmos_input%press(:,1,kr)  =ps(istart:iend)*sig(k)
      call getqsat(imax,qsat,t(istart:iend,k),ps(istart:iend)*sig(k))
      Atmos_input%rel_hum(:,1,kr)=min(qg(istart:iend,k)/qsat,1.)
    end do
    Atmos_input%temp(:,1,kl+1)  = tss(istart:iend)
    Atmos_input%press(:,1,kl+1) = ps(istart:iend)
    Atmos_input%pflux(:,1,1  )  = 0.
    Atmos_input%tflux(:,1,1  )  = Atmos_input%temp (:,1,1  )
    do k=2,kl
      Atmos_input%pflux(:,1,k) = 0.5*(Atmos_input%press(:,1,k-1)+Atmos_input%press(:,1,k))
      Atmos_input%tflux(:,1,k) = 0.5*(Atmos_input%temp (:,1,k-1)+Atmos_input%temp (:,1,k))
    end do
    Atmos_input%pflux(:,1,kl+1) = Atmos_input%press(:,1,kl+1)
    Atmos_input%tflux(:,1,kl+1) = Atmos_input%temp (:,1,kl+1)
    Atmos_input%clouddeltaz     = Atmos_input%deltaz        
    !Atmos_input%cloudtemp       = Atmos_input%temp ! fix
    !do k=1,kl
    !  kr = kl+1-k
    !  Atmos_input%cloudvapor(:,1,kr)=qg(istart:iend,k) ! fix
    !end do
    !Atmos_input%aerosolrelhum   =Atmos_input%rel_hum
    !Atmos_input%aerosoltemp     =Atmos_input%temp ! fix
    !Atmos_input%aerosolvapor    =Atmos_input%cloudvapor ! fix
    !Atmos_input%aerosolpress    =Atmos_input%press  ! fix    

    Atmos_input%psfc(:,1)    =ps(istart:iend)
    Atmos_input%tsfc(:,1)    =tss(istart:iend)
    do k=1,kl+1
      kr=kl+2-k
      Atmos_input%phalf(:,1,kr)=ps(istart:iend)*sigh(k)
    end do
    
    Rad_gases%rrvco2 =rrvco2
    Rad_gases%rrvch4 =rrvch4
    Rad_gases%rrvn2o =rrvn2o
    Rad_gases%rrvf11 =rrvf11
    Rad_gases%rrvf12 =rrvf12
    Rad_gases%rrvf113=rrvf113
    Rad_gases%rrvf22 =rrvf22
    
    Cld_spec%camtsw=0.
    Cld_spec%crndlw=0.
    Cld_spec%cmxolw=0.
    !Cld_spec%ncldsw=0
    !Cld_spec%nrndlw=0
    !Cld_spec%nmxolw=0
    if (nmr.eq.0) then
      do i=1,imax
        iq=i+istart-1
        newcld=.true.
        do k=1,kl
          kr=kl+1-k
          if (cfrac(iq,k).gt.0.) then
            Cld_spec%camtsw(i,1,kr)=cfrac(iq,k) ! Max+Rnd overlap clouds for SW
            !Cld_spec%ncldsw=Cld_spec%ncldsw+1
            Cld_spec%crndlw(i,1,kr)=cfrac(iq,k) ! Rnd overlap for LW
            !Cld_spec%nrndlw=Cld_spec%nrndlw+1
          end if
        end do
      end do
    else
      do i=1,imax
        iq=i+istart-1
        newcld=.true.
        do k=1,kl
          kr=kl+1-k
          if (cfrac(iq,k).gt.0.) then
            Cld_spec%camtsw(i,1,kr)=cfrac(iq,k) ! Max+Rnd overlap clouds for SW
            !Cld_spec%ncldsw=Cld_spec%ncldsw+1
            maxover=.false.
            if (k.gt.1) then
              if (cfrac(iq,k-1).gt.0.) maxover=.true.
            end if
            if (k.lt.kl) then
              if (cfrac(iq,k+1).gt.0.) maxover=.true.
            end if
            if (maxover) then
              Cld_spec%cmxolw(i,1,kr)=cfrac(iq,k) ! Max overlap for LW
              !if (newcld) Cld_spec%nmxolw=Cld_spec%nmxolw+1
              newcld=.false.
            else
              Cld_spec%crndlw(i,1,kr)=cfrac(iq,k) ! Rnd overlap for LW
              !Cld_spec%nrndlw=Cld_spec%nrndlw+1
            end if
          else
            newcld=.true.
          end if
        end do
      end do
    end if

    do k=1,kl
      p2(:,k)=ps(istart:iend)*sig(k) !
    end do
    call cloud3(Cloud_microphysics%size_drop,Cloud_microphysics%size_ice,       &
                Cloud_microphysics%conc_drop,Cloud_microphysics%conc_ice,       &
                cfrac(istart:iend,:),qlrad(istart:iend,:),qfrad(istart:iend,:), &
                p2,t(istart:iend,:),cd2,imax,kl)
    Cloud_microphysics%size_drop=max(Cloud_microphysics%size_drop,1.e-20)
    Cloud_microphysics%size_ice =max(Cloud_microphysics%size_ice,1.e-20)                
    Cloud_microphysics%size_rain=1.e-20
    Cloud_microphysics%conc_rain=0.
    Cloud_microphysics%size_snow=1.e-20
    Cloud_microphysics%conc_snow=0.

    Lscrad_props%cldext   = 0.
    Lscrad_props%cldsct   = 0.
    Lscrad_props%cldasymm = 1.
    Lscrad_props%abscoeff = 0.
    call microphys_lw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    call microphys_sw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    Cldrad_props%cldsct(:,:,:,:,1)  =Lscrad_props%cldsct(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldext(:,:,:,:,1)  =Lscrad_props%cldext(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldasymm(:,:,:,:,1)=Lscrad_props%cldasymm(:,:,:,:) ! Large scale cloud properties only
    Cldrad_props%abscoeff(:,:,:,:,1)=Lscrad_props%abscoeff(:,:,:,:) ! Large scale cloud properties only
    
    call lwemiss_calc(Atmos_input%clouddeltaz,Cldrad_props%abscoeff,Cldrad_props%cldemiss)
    Cldrad_props%emmxolw = Cldrad_props%cldemiss
    Cldrad_props%emrndlw = Cldrad_props%cldemiss
    
    Surface%asfc_vis_dir(:,1)=cuvrf_dir(:)
    Surface%asfc_nir_dir(:,1)=cirrf_dir(:)
    Surface%asfc_vis_dif(:,1)=cuvrf_dif(:)
    Surface%asfc_nir_dif(:,1)=cirrf_dif(:)
   
    Astro%cosz(:,1)   =max(coszro,0.)
    Astro%fracday(:,1)=taudar

    call longwave_driver (1, imax, 1, 1, Rad_time, Atmos_input,  &
                          Rad_gases, Aerosol, Aerosol_props,     &
                          Cldrad_props, Cld_spec, Aerosol_diags, &
                          Lw_output)

    call shortwave_driver (1, imax, 1, 1, Atmos_input, Surface,      &
                           Astro, Aerosol, Aerosol_props, Rad_gases, &
                           Cldrad_props, Cld_spec, Sw_output,        &
                           Aerosol_diags, r)

    ! store shortwave and fbeam data --------------------------------
    sg=Sw_output(1)%dfsw(:,1,kl+1)-Sw_output(1)%ufsw(:,1,kl+1)
    sgvis=Sw_output(1)%dfsw_vis_sfc(:,1)-Sw_output(1)%ufsw_vis_sfc(:,1)
    !sgvisdir=Sw_output(1)%dfsw_vis_sfc_dir(:,1)
    !sgvisdif=Sw_output(1)%dfsw_vis_sfc_dif(:,1)-Sw_output(1)%ufsw_vis_sfc_dif(:,1)
    !sgnirdir=Sw_output(1)%dfsw_dir_sfc(:,1)-sgvisdir
    !sgnirdif=Sw_output(1)%dfsw_dif_sfc(:,1)-Sw_output(1)%ufsw_dif_sfc(:,1)-sgvisdif
    !sgdir=Sw_output(1)%dfsw_dir_sfc(:,1)
    !sgdif=Sw_output(1)%dfsw_dif_sfc(:,1)-Sw_output(1)%ufsw_dif_sfc(:,1)
    sgdnvisdir=Sw_output(1)%dfsw_vis_sfc_dir(:,1)
    sgdnvisdif=Sw_output(1)%dfsw_vis_sfc_dif(:,1)
    sgdnnirdir=Sw_output(1)%dfsw_dir_sfc(:,1)-sgdnvisdir
    sgdnnirdif=Sw_output(1)%dfsw_dif_sfc(:,1)-sgdnvisdif
    
    where (sg.gt.0.1)
      swrsave(istart:iend)=sgvis/sg
    elsewhere
      swrsave(istart:iend)=0.
    end where
    where (sgdnvisdir+sgdnvisdif.gt.0.1)
      fbeamvis(istart:iend)=sgdnvisdir/(sgdnvisdir+sgdnvisdif)
    elsewhere
      fbeamvis(istart:iend)=0.
    end where
    where (sgdnnirdir+sgdnnirdif.gt.0.1)
      fbeamnir(istart:iend)=sgdnnirdir/(sgdnnirdir+sgdnnirdif)
    elsewhere
      fbeamnir(istart:iend)=0.
    end where
    
    ! Store albedo data ---------------------------------------------
    albvisnir(istart:iend,1)=cuvrf_dir(1:imax)*fbeamvis(istart:iend)+cuvrf_dif(1:imax)*(1.-fbeamvis(istart:iend))
    albvisnir(istart:iend,2)=cirrf_dir(1:imax)*fbeamnir(istart:iend)+cirrf_dif(1:imax)*(1.-fbeamnir(istart:iend))
    
    ! longwave output -----------------------------------------------
    rg(1:imax) = Lw_output(1)%flxnet(:,1,kl+1)          ! longwave at surface
    rt(1:imax) = Lw_output(1)%flxnet(:,1,1)             ! longwave at top
    ! rg is net upwards = sigma T^4 - Rdown
    rgdn(1:imax) = stefbo*tss(istart:iend)**4 - rg(1:imax)

    ! shortwave output ----------------------------------------------
    sint(1:imax) = Sw_output(1)%dfsw(:,1,1)   ! solar in top
    sout(1:imax) = Sw_output(1)%ufsw(:,1,1)   ! solar out top
    sgdn(1:imax) = sg(1:imax) / ( 1. - swrsave(istart:iend)*albvisnir(istart:iend,1) &
                  -(1.-swrsave(istart:iend))*albvisnir(istart:iend,2) ) ! MJT albedo

    ! Clear sky calculation -----------------------------------------
    if (do_totcld_forcing) then
      soutclr(1:imax) = Sw_output(1)%ufswcf(:,1,1)      ! solar out top
      sgclr(1:imax)   = -Sw_output(1)%fswcf(:,1,kl+1)   ! solar absorbed at the surface
      rtclr(1:imax)   = Lw_output(1)%flxnetcf(:,1,1)    ! clr sky lw at top
      rgclr(1:imax)   = Lw_output(1)%flxnetcf(:,1,kl+1) ! clear sky longwave at surface
    else
      soutclr(1:imax) = 0.
      sgclr(1:imax) = 0.
      rtclr(1:imax) = 0.
      rgclr(1:imax) = 0.
    end if

    ! heating rate --------------------------------------------------
    do k=1,kl
      ! total heating rate (convert deg K/day to deg K/sec)
      rtt(istart:iend,kl+1-k)=-(Sw_output(1)%hsw(:,1,k)+Lw_output(1)%heatra(:,1,k))/86400.
    end do

    ! Calculate the amplitude of the diurnal cycle of solar radiation
    ! at the surface (using the value for the middle of the radiation
    ! step) and use this value to get solar radiation at other times.
    ! Use the zenith angle and daylight fraction calculated in zenith
    ! to remove these factors.
    where (coszro(1:imax)*taudar(1:imax).le.1.E-5)
      ! The sun isn't up at all over the radiation period so no 
      ! fitting need be done.
      sga(1:imax)=0.
    else where
      sga(1:imax)=sg(1:imax)/(coszro(1:imax)*taudar(1:imax))
    end where

    ! Save things for non-radiation time steps ----------------------
    sgsave(istart:iend)   = sg(1:imax)   ! repeated after solarfit
    sgamp(istart:iend)    = sga(1:imax)
    ! Save the value excluding Ts^4 part.  This is allowed to change.
    rgsave(istart:iend)   = rg(1:imax)-stefbo*tss(istart:iend)**4  ! opposite sign to prev. darlam scam
    sintsave(istart:iend) = sint(1:imax) 
    rtsave(istart:iend)   = rt(1:imax) 
    rtclsave(istart:iend) = rtclr(1:imax)  
    sgclsave(istart:iend) = sgclr(1:imax)

    ! cloud amounts for saving --------------------------------------
    cloudtot(istart:iend)=1.-(1.-cloudlo(istart:iend))*(1.-cloudmi(istart:iend))*(1.-cloudhi(istart:iend))

    ! Use explicit indexing rather than array notation so that we can run
    ! over the end of the first index
    if(ktau>1)then ! averages not added at time zero
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
      cld_ave(istart:iend)  = cld_ave(istart:iend)  + cloudtot(istart:iend)
      cll_ave(istart:iend)  = cll_ave(istart:iend)  + cloudlo(istart:iend)
      clm_ave(istart:iend)  = clm_ave(istart:iend)  + cloudmi(istart:iend)
      clh_ave(istart:iend)  = clh_ave(istart:iend)  + cloudhi(istart:iend)
    endif   ! (ktau>1)

  end if  ! odcalc
      
  ! Calculate the solar using the saved amplitude.
  sg(1:imax) = sgamp(istart:iend)*coszro2(1:imax)*taudar2(1:imax)
  if(ktau>1)then ! averages not added at time zero
    sgn_ave(istart:iend)  = sgn_ave(istart:iend)  + sg(1:imax)
  endif  ! (ktau>1)
      
  ! Set up the CC model radiation fields
  ! slwa is negative net radiational htg at ground
  ! Note that this does not include the upward LW radiation from the surface.
  ! That is included in sflux.f
  slwa(istart:iend) = -sg(1:imax)+rgsave(istart:iend)
  sgsave(istart:iend) = sg(1:imax)   ! this is the repeat after solarfit 26/7/02

end do  ! Row loop (j)  j=1,jl,imax/il

! Calculate net radiational cooling of atmosphere (K/s)
t(1:ifull,:)=t(1:ifull,:)-dt*rtt(1:ifull,:)

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

integer,                      intent(in)     :: is, ie, js, je
type(time_type),              intent(in)     :: Rad_time
type(atmos_input_type),       intent(in)     :: Atmos_input  
type(radiative_gases_type),   intent(inout)  :: Rad_gases   
type(aerosol_type),           intent(in)     :: Aerosol     
type(aerosol_properties_type),intent(inout)  :: Aerosol_props
type(aerosol_diagnostics_type),intent(inout)  :: Aerosol_diags
type(cldrad_properties_type), intent(in)     :: Cldrad_props
type(cld_specification_type), intent(in)     :: Cld_spec     
type(lw_output_type), dimension(:),  intent(inout)  :: Lw_output
type(lw_diagnostics_type)  :: Lw_diagnostics

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

!--------------------------------------------------------------------
!   local variables

      integer  :: ix, jx, kx  ! dimensions of current physics window

!--------------------------------------------------------------------
!    call longwave_driver_alloc to allocate component arrays of a
!    lw_output_type variable.
!----------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
        call sealw99 (is, ie, js, je, Rad_time, Atmos_input,           &
                      Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
                      Cld_spec, Aerosol_diags, Lw_output(1),           &
                      Lw_diagnostics, do_aerosol_forcing)

!--------------------------------------------------------------------
!    deallocate the components of Lw_diagnostics.
!--------------------------------------------------------------------
      deallocate (Lw_diagnostics%flx1e1)
      deallocate (Lw_diagnostics%fluxn )
      deallocate (Lw_diagnostics%cts_out)
      deallocate (Lw_diagnostics%cts_outcf)
      deallocate (Lw_diagnostics%gxcts )
      deallocate (Lw_diagnostics%excts )
      deallocate (Lw_diagnostics%exctsn)
      deallocate (Lw_diagnostics%fctsg )
      deallocate (Lw_diagnostics%flx1e1f)
      if (Rad_control%do_totcld_forcing) then
        deallocate (Lw_diagnostics%fluxncf)
      endif        

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
      logical  :: skipswrad
      logical  :: with_clouds
      logical  :: calc_includes_aerosols
      integer  :: naerosol_optical
      integer  :: i, j       
      integer  :: ix, jx, kx

!---------------------------------------------------------------------
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      with_clouds  are clouds to be considered in determining
!                   the sw fluxes and heating rates ?
!      ix,jx,kx     dimensions of current physics window
!      i,j          do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    call shortwave_driver_alloc to initialize shortwave fluxes and 
!    heating rates.
!--------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1

!--------------------------------------------------------------------
!    allocate and initialize fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      Sw_output(1)%fsw   (:,:,:) = 0.0
      Sw_output(1)%dfsw  (:,:,:) = 0.0
      Sw_output(1)%ufsw  (:,:,:) = 0.0
      Sw_output(1)%hsw   (:,:,:) = 0.0
      Sw_output(1)%dfsw_dir_sfc = 0.0
      Sw_output(1)%dfsw_dif_sfc  = 0.0
      Sw_output(1)%ufsw_dif_sfc = 0.0
      Sw_output(1)%dfsw_vis_sfc = 0.
      Sw_output(1)%ufsw_vis_sfc = 0.
      Sw_output(1)%dfsw_vis_sfc_dir = 0.
      Sw_output(1)%dfsw_vis_sfc_dif = 0.
      Sw_output(1)%ufsw_vis_sfc_dif = 0.
      Sw_output(1)%bdy_flx(:,:,:) = 0.0       

!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then
        Sw_output(1)%fswcf (:,:,:) = 0.0
        Sw_output(1)%dfswcf(:,:,:) = 0.0
        Sw_output(1)%ufswcf(:,:,:) = 0.0
        Sw_output(1)%hswcf (:,:,:) = 0.0
        Sw_output(1)%dfsw_dir_sfc_clr = 0.0
        Sw_output(1)%dfsw_dif_sfc_clr  = 0.0
        Sw_output(1)%bdy_flx_clr (:,:,:) = 0.0
      endif

!--------------------------------------------------------------------
!    determine when the no-sun case exists at all points within the 
!    physics window and bypass the sw radiation calculations for that 
!    window. for do_annual_mean or do_daily_mean, only one cosz in a
!    model row need be tested, since all points in i have the same 
!    zenith angle.
!--------------------------------------------------------------------
      skipswrad = .true.
      do j=1,jx        
        do i = 1,ix         
          if (Astro%cosz(i,j) > 0.0 )  then
            skipswrad = .false.
            exit
          endif
        end do
      end do


!--------------------------------------------------------------------
!    if the sun is shining nowhere in the physics window allocate
!    output fields which will be needed later, set them to a flag
!    value and return.
!--------------------------------------------------------------------
      if (skipswrad)  then

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      else 
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
          if (do_aerosol_forcing) then
            naerosol_optical = size (Aerosol_props%aerextband,2)
          else
            naerosol_optical = 0  
          endif 
          call swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,&
                       Aerosol, Aerosol_props, Astro, Cldrad_props,  &
                       Cld_spec, calculate_volcanic_sw_heating, &
                       Sw_output(1), Aerosol_diags, r,  &
                       do_aerosol_forcing, naerosol_optical)
      endif
!--------------------------------------------------------------------

end subroutine shortwave_driver

subroutine getqsat(ifull,qsat,temp,ps)

implicit none

integer, intent(in) :: ifull
real, dimension(ifull), intent(in) :: temp,ps
real, dimension(ifull), intent(out) :: qsat
real, dimension(ifull) :: esatf
real, parameter :: latent= 2.50E6
real, parameter :: latsub= 2.83E6
real, parameter :: rv    = 461.5

where (temp.ge.273.15)
  esatf = 610.*exp(latent/rv*(1./273.15-1./min(max(temp,123.),343.)))
elsewhere
  esatf = 610.*exp(latsub/rv*(1./273.15-1./min(max(temp,123.),343.)))
endwhere
qsat = 0.622*esatf/(ps-0.378*esatf)

return
end subroutine getqsat

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlg,qfg,prf,ttg,cdrop,imax,kl)

implicit none

integer, intent(in) :: imax,kl
integer k,kr,mg
real, dimension(imax,kl), intent(in) :: cfrac,qlg,qfg,prf,ttg
real, dimension(imax,kl), intent(in) :: cdrop
real(kind=8), dimension(imax,kl), intent(out) :: Rdrop,Rice,conl,coni
real, dimension(imax,kl) :: reffl,reffi,fice,cfl,Wliq,rhoa
real, dimension(imax,kl) :: eps,rk,Wice

!--------------------------------------------------------------------
!    if liquid is present in the layer, compute the effective drop
!    radius. the following formula, recommended by (Martin et al., 
!    J. Atmos. Sci, vol 51, pp. 1823-1842) is used for liquid droplets:
!    reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = density of cloud droplets (number per cubic meter)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------

! Reffl is the effective radius at the top of the cloud (calculated following
! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
! formula for reffl. Use mid cloud value of Reff for emissivity.

where (cfrac.gt.0.)
  fice= qfg/max(qfg+qlg,1.e-12)
else where
  fice=0.
end where
cfl=cfrac*(1.-fice)
rhoa=prf/(rdry*ttg)

where (qlg.gt.1.E-8.and.cfrac.gt.0.)
  Wliq=rhoa*qlg/cfl     !kg/m^3
  ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
  eps = 1. - 0.7 * exp(-0.003e-6*cdrop) !mid range
  rk  = (1.+eps**2)/(1.+2.*eps**2)**2
  ! Martin et al 1994
  reffl=(3.*2.*Wliq/(4.*pi*rhow*rk*cdrop))**(1./3.)
elsewhere
  reffl=0.
  Wliq=0.
end where
where (qfg.gt.1.E-8.and.cfrac.gt.0.)
  Wice=rhoa*qfg/(cfrac*fice) !kg/m**3
  reffi=min(150.e-6,3.73e-4*Wice**0.216) !Lohmann et al.(1999)
elsewhere
  reffi=0.
  Wice=0.
end where

do k=1,kl
  kr=kl+1-k
  Rdrop(:,kr)=2.*reffl(:,k)*1.E6 ! convert to diameter and microns
  Rice(:,kr)=2.*reffi(:,k)*1.E6
  conl(:,kr)=1000.*Wliq(:,k)
  coni(:,kr)=1000.*Wice(:,k)
end do

where (Rdrop.gt.0.)
  Rdrop=min(max(Rdrop,8.4),33.2) ! constrain diameter to acceptable range (see microphys_rad.f90)
endwhere
where (Rice.gt.0.)
  Rice=min(max(Rice,18.6),130.2)
endwhere

return
end subroutine cloud3


end module seaesfrad_m
