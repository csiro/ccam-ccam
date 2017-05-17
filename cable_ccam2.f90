! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! CABLE interface originally developed by the CABLE group
! Subsequently modified by MJT for maxtile=5 tile mosaic and SEAESF radiation scheme
  
! - Currently all tiles have the same soil texture, but independent soil temperatures,
!   moisture, etc.
! - Data can be read as CSIRO PFT or as IGBP
! - LAI can be interpolated between timesteps using a PW-linear fit to the LAI integral
!   or LAI can be taken as constant for the month.  CASA-CNP can predict LAI and vcmax,
!   although this can take considerable time to spin-up.
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.
! - Options exist for using SLI soil model (soil_struc=1), climate feedbacks
!   (cable_climate=1), phenology based on climate (phenology_switch=1) and the POP
!   model (cable_pop=1)

! The following mappings between IGBP and CSIRO PFT were recommended by Rachel Law
    
! ivegt   IGBP type                             CSIRO PFT
! 1       Evergreen Needleleaf Forest           1.  Evergreen Needleleaf
! 2       Evergreen Broadleaf Forest            1.  Evergreen Broadleaf
! 3       Deciduous Needleaf Forest             1.  Deciduous Needleleaf
! 4       Deciduous Broadleaf Forest            1.  Deciduous Broadleaf
! 5       Mixed Forest                          1.  Deciduous Broadleaf                              when -25<lat<25
!                                               0.5 Evergreen Needleleaf    0.5 Deciduous Broadleaf  when lat<-25 or lat>25
! 6       Closed Shrublands                     0.8 Shrub                   0.2 (Grass)
! 7       Open Shrublands                       0.2 Shrub                   0.8 (Grass)
! 8       Woody Savannas                        0.6 (Grass)                 0.4 Evergreen Needleleaf when lat<-40 or lat>40
!                                               0.6 (Grass)                 0.4 Evergreen Broadleaf  when -40<lat<40
! 9       Savannas                              0.9 (Grass)                 0.1 Evergreen Needleleaf when lat<-40 or lat>40
!                                               0.9 (Grass)                 0.1 Evergreen Broadleaf  when -40<lat<40
! 10      Grasslands                            1.  (Grass)
! 11      Permanent Wetlands                    1.  Wetland
! 12      Croplands                             1.  (Crop)
! 13      Urban and Built-up                    1.  Urban
! 14      Cropland/Natural Vegetation Mosaic    1.  (Crop)
! 15      Snow and Ice                          1.  Ice
! 16      Barren or Sparsely Vegetated          1.  Barren
! 17      Water Bodies                          1.  Lakes

! where:
!   (Grass)   0.9  C3 0.1  C4 0. Tundra   40<lat<50 or -50<lat<-40
!             0.8  C3 0.2  C4 0. Tundra   30<lat<40 or -40<lat<-30
!             0.5  C3 0.5  C4 0. Tundra   25<lat<30 or -30<lat<-25
!             0.05 C3 0.95 C4 0. Tundra  -25<lat<25
!             0.   C3 0.   C4 1. Tundra   lat<-50 or lat>50

!   (Crop)    0.7 C3  0.3 C4   -30<lat<30
!             0.9 C3  0.1 C4   30<lat<40 or -40<lat<-30
!             1.  C3  0.  C4   lat<-40   or lat>40

! *** NOTE MJT SPECIAL
! The PFT evergreen broadleaf's canopy height is reduced for woody savannas for improved roughness length

! CSIRO PFT index
! 1  Evergreen Needleleaf
! 2  Evergreen Broadleaf
! 3  Deciduous Needleaf
! 4  Deciduous Broadleaf
! 5  Shrub
! 6  C3 grass
! 7  C4 grass
! 8  Tundra
! 9  C3 crop
! 10 C4 crop
! 11 Wetland
! 12 Not used
! 13 Not used
! 14 Barren
! 15 Urban
! 16 Lakes
! 17 Ice
! (18 Evergreen Broadleaf (Savanna)) - MJT defined
  
! isoilm  type
! 0       water/ocean
! 1       coarse               sand/loamy_sand
! 2       medium               clay-loam/silty-clay-loam/silt-loam
! 3       fine                 clay
! 4       coarse-medium        sandy-loam/loam
! 5       coarse-fine          sandy-clay
! 6       medium-fine          silty-clay 
! 7       coarse-medium-fine   sandy-clay-loam
! 8       organi!              peat
! 9       land ice

module cable_ccam

use cable_air_module
use cable_albedo_module
use cable_canopy_module
use cable_carbon_module
use cable_common_module
use cable_data_module
use cable_def_types_mod, cbm_ms => ms
use cable_radiation_module
use cable_roughness_module
use cable_soil_snow_module
use casa_cnp_module
use casadimension
use casaparm, xroot => froot
use casavariable
use pop_types

implicit none

private
public sib4
public loadcbmparm, cbmparm, loadtile, defaulttile, savetiledef, savetile, newcbmwb
public cablesettemp, cableinflow, cbmemiss
public proglai, progvcmax, maxtile, soil_struc, cable_pop, phenology_switch
public fwsoil_switch, cable_litter, gs_switch, cable_climate

! CABLE - CCAM options.
integer, save :: proglai            = -1         ! -1, piece-wise linear prescribed LAI, 0 PWCB prescribed LAI, 1 prognostic LAI
integer, save :: progvcmax          = 0          ! 0 prescribed, 1 prognostic vcmax (default), 2 prognostic vcmax (Walker2014)
integer, save :: soil_struc         = 0          ! 0 default, 1 SLI soil model
integer, save :: cable_pop          = 0          ! 0 POP off, 1 POP on
integer, save :: fwsoil_switch      = 0          ! 0 default, 1 non-linear, 2 Lai and Ktaul, 3 Haverd2013
integer, save :: gs_switch          = 0          ! 0 leuning, 1 medlyn
integer, save :: cable_litter       = 0          ! 0 off, 1 on
integer, save :: cable_climate      = 0          ! 0 off, 1 on
integer, save :: phenology_switch   = 0          ! 0 default, 1 climatology
integer, parameter :: tracerco2     = 0          ! 0 use radiation CO2, 1 use tracer CO2 
integer, parameter :: maxtile       = 5          ! maximum possible number of tiles
integer, parameter :: COLDEST_DAY_NHEMISPHERE = 355
integer, parameter :: COLDEST_DAY_SHEMISPHERE = 172
real, parameter :: minfrac = 0.01                ! minimum non-zero tile fraction (improves load balancing)

integer, dimension(maxtile,2), save :: pind  
integer, save :: maxnb                  ! maximum number of actual tiles
real, dimension(:), allocatable, save :: sv, vl1, vl2, vl3, vl4
logical, dimension(:,:), allocatable, save :: tmap
type (air_type), save            :: air
type (bgc_pool_type), save       :: bgc
type (met_type), save            :: met
type (balances_type), save       :: bal
type (radiation_type), save      :: rad
type (roughness_type), save      :: rough
type (soil_snow_type), save      :: ssnow
type (sum_flux_type), save       :: sum_flux
type (climate_type), save        :: climate
type (veg_parameter_type), save  :: veg
type (soil_parameter_type), save :: soil
type (canopy_type), save         :: canopy
type (casa_balance), save        :: casabal
type (casa_biome), save          :: casabiome
type (casa_flux), save           :: casaflux
type (casa_met), save            :: casamet
type (casa_pool), save           :: casapool
type (phen_variable), save       :: phen
type (pop_type), save            :: pop
type (physical_constants), save  :: c

contains
! ****************************************************************************

! CABLE-CCAM interface
subroutine sib4

use arrays_m
use carbpools_m
use const_phys
use dates_m
use estab
use extraout_m
use infile
use latlong_m
use liqwpar_m
use morepbl_m
use newmpar_m
use nharrs_m
use nsibd_m
use parm_m
use pbl_m
use permsurf_m
use prec_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use vegpar_m
use work2_m, only : qsttg,zo,zoh,zoq,theta,vmod,wetfac
use work3_m, only : ga
use zenith_m
  
implicit none

real fjd, r1, dlt, slag, dhr, alp
real, dimension(ifull) :: coszro2, taudar2, tmps, atmco2
real, dimension(ifull) :: swdwn, alb, qsttg_land
real, dimension(mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
real, dimension(mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
real, dimension(mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
real, dimension(:), allocatable, save :: tmp
real(r_2), dimension(mp) :: xKNlimiting, xkleafcold, xkleafdry
real(r_2), dimension(mp) :: xkleaf, xnplimit, xNPuptake, xklitter
real(r_2), dimension(mp) :: xksoil
integer jyear, jmonth, jday, jhour, jmin, idoy
integer k, mins, nb, j
integer is, ie, casaperiod, npercasa
integer lalloc

cansto = 0.
fwet = 0.
vlai = 0.

! abort calculation if no land points on this processor  
if ( mp<=0 ) return

! calculate time from beginning of the simulation
! jyear, jmonth, jday, jhour, jmin are for the start date of the simulation
! mins is the time from the start of the year
! mtimer is the time elapsed from the start of the simulation
dhr = dt/3600.
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
! calculate zenith angle
fjd = float(mod(mins, 525600))/1440. ! restrict to 365 day calendar
idoy = int(fjd)
call solargh(fjd,bpyear,r1,dlt,alp,slag)
call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro2,taudar2)

! calculate CO2 concentration
call setco2for(atmco2)
if ( any(atmco2<1.) ) then
  write(6,*) "ERROR: Invalid CO2 mixing ratio in cable_ccam2 ",minval(atmco2)
  stop
end if

! set meteorological forcing
! swdwn is downwelling shortwave (positive) W/m^2
albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif
albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif
alb   = swrsave*albvissav + (1.-swrsave)*albnirsav
if ( any(alb<1.e-20) ) then
  write(6,*) "ERROR: Invalid albedo in cable_ccam2 ",minval(alb)
  stop
end if
swdwn = sgsave/(1.-alb)
do nb = 1,maxnb
  is = pind(nb,1)
  ie = pind(nb,2)
  met%tk(is:ie)          = pack(theta,  tmap(:,nb))
  met%ua(is:ie)          = pack(vmod,   tmap(:,nb))
  met%ca(is:ie)          = pack(atmco2, tmap(:,nb))*1.e-6
  met%coszen(is:ie)      = pack(coszro2,tmap(:,nb))              ! use instantaneous value
  met%qv(is:ie)          = pack(qg(1:ifull,1),tmap(:,nb))        ! specific humidity in kg/kg
  met%pmb(is:ie)         = pack(ps(1:ifull),  tmap(:,nb))*0.01   ! pressure in mb at ref height
  met%precip(is:ie)      = pack(condx,  tmap(:,nb))              ! in mm not mm/sec
  met%precip_sn(is:ie)   = pack(conds+condg,  tmap(:,nb))        ! in mm not mm/sec
  ! swrsave indicates the fraction of net VIS radiation (compared to NIR)
  ! fbeamvis indicates the beam fraction of downwelling direct radiation (compared to diffuse) for VIS
  ! fbeamnir indicates the beam fraction of downwelling direct radiation (compared to diffuse) for NIR
  met%fsd(is:ie,1)       = pack(swrsave*swdwn,        tmap(:,nb))
  met%fsd(is:ie,2)       = pack((1.-swrsave)*swdwn,   tmap(:,nb))
  met%fld(is:ie)         = pack(-rgsave,              tmap(:,nb))      ! long wave down (positive) W/m^2
  rad%fbeam(is:ie,1)     = pack(fbeamvis,             tmap(:,nb))
  rad%fbeam(is:ie,2)     = pack(fbeamnir,             tmap(:,nb))
  rough%za_tq(is:ie)     = pack(bet(1)*t(1:ifull,1)+phi_nh(:,1),tmap(:,nb))/grav ! reference height
end do
met%doy         = fjd
met%hod         = (met%doy-int(met%doy))*24. + rad%longitude/15.
met%hod         = mod(met%hod, 24.)
met%tvair       = met%tk
met%tvrad       = met%tk
met%qvair       = met%qv
met%ua          = max(met%ua, c%umin)
met%coszen      = max(met%coszen, 1.e-8) 
rough%za_uv     = rough%za_tq
rad%fbeam(:,3)  = 0.            ! dummy for now

! Interpolate LAI.  Also need sigmf for LDR prognostic aerosols.
call setlai(sigmf,jyear,jmonth,jday,jhour,jmin,mp)

! Calculate vcmax
call vcmax_feedback

!--------------------------------------------------------------
! CABLE
ktau_gl          = 900
kend_gl          = 999
met%ofsd         = met%fsd(:,1) + met%fsd(:,2)
ssnow%owetfac    = ssnow%wetfac
canopy%oldcansto = canopy%cansto
call ruff_resist(veg,rough,ssnow,canopy)
call define_air(met,air)
call init_radiation(met,rad,veg,canopy)
call surface_albedo(ssnow,veg,met,rad,soil,canopy)
call define_canopy(bal,rad,rough,air,met,dt,ssnow,soil,veg,canopy,climate)
ssnow%otss_0  = ssnow%otss
ssnow%otss    = ssnow%tss
ssnow%owetfac = ssnow%wetfac
select case ( soil_struc )
  case(0)  
    call soil_snow(dt,soil,ssnow,canopy,met,bal,veg)
  case(1)
    call sli_main(999,dt,veg,soil,ssnow,met,canopy,air,rad,0)    
  case default
    write(6,*) "ERROR: Unknown option soil_struc ",soil_struc
    stop
end select
! adjust for new soil temperature
ssnow%deltss     = ssnow%tss - ssnow%otss
canopy%fhs       = canopy%fhs + ssnow%deltss*ssnow%dfh_dtg
canopy%fes       = canopy%fes + ssnow%deltss*ssnow%dfe_ddq*ssnow%ddq_dtg
canopy%fns       = canopy%fns + ssnow%deltss*ssnow%dfn_dtg
canopy%ga        = real(canopy%ga  + ssnow%deltss*canopy%dgdtg)
!canopy%fhs_cor   = canopy%fhs_cor + ssnow%deltss*ssnow%dfh_dtg
!canopy%fes_cor   = canopy%fes_cor + ssnow%deltss*ssnow%dfe_ddq*ssnow%ddq_dtg
canopy%fh        = canopy%fhv + canopy%fhs
canopy%fev       = real(canopy%fevc + canopy%fevw)
canopy%fe        = real(canopy%fev + canopy%fes)
canopy%rnet      = canopy%fns + canopy%fnv
rad%trad         = ( (1.-rad%transd)*canopy%tv**4 + rad%transd*ssnow%tss**4 )**0.25

! note that conservation is still preserved at this point
! canopy%ga    = canopy%rnet - canopy%fh - canopy%fe
! canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq*ssnow%ddq_dtg

! EK suggestion
!canopy%cdtq =  max( 0.1*canopy%cduv, canopy%cdtq )
! MJT suggestion
canopy%cdtq =  max( 0., canopy%cdtq )

!--------------------------------------------------------------
! CABLE CLIMATE
call cableclimate(idoy,jmonth,jyear)

!--------------------------------------------------------------
! CASA CNP
select case (icycle)
  case(0) ! off
    !call plantcarb(veg,bgc,met,canopy)
    !call soilcarb(soil,ssnow,veg,bgc,met,canopy)
    !call carbon_pl(dt,soil,ssnow,veg,canopy,bgc)
    !canopy%fnpp = -canopy%fpn - canopy%frp
    !canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
  case(3) ! C+N+P
    ! update casamet
    if ( cable_pop==1 ) then
      lalloc = 3  
    else
      lalloc = 0  
    end if
    ! run CASA CNP once per day
    casaperiod = nint(86400.*deltpool)
    if ( abs(mod(real(casaperiod),dt))>1.e-20 ) then
      write(6,*) "ERROR: Invalid casaperiod ",casaperiod
      write(6,*) "Period must be an integer multiple of the time-step dt ",dt
      write(6,*) "Please adjust deltpool from CABLE CASA-CNP accordingly"
      stop
    end if
    npercasa = max( nint(real(casaperiod)/dt), 1 )
    casamet%tairk = casamet%tairk + met%tk
    casamet%tsoil = casamet%tsoil + ssnow%tgg
    casamet%moist = casamet%moist + ssnow%wb
    casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dt
    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dt
    if ( mod(ktau,npercasa)==0 ) then
      casamet%tairk=casamet%tairk/real(npercasa)
      casamet%tsoil=casamet%tsoil/real(npercasa)
      casamet%moist=casamet%moist/real(npercasa)
      xKNlimiting = 1.
      if ( phenology_switch==1 .and. cable_climate==1 ) then
        call cable_phenology_clim
      else
        call phenology(idoy,veg,phen)
      end if
      call avgsoil(veg,soil,casamet)
      call casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
      if ( cable_pop/=1 ) then
        call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,lalloc)
      end if
      call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casamet,phen)
      call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool,casaflux,casamet,phen)
      call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)
      if ( cable_pop==1 ) then
        call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,lalloc)
        where (pop%pop_grid(:)%cmass_sum_old>0.001 .and. pop%pop_grid(:)%cmass_sum>0.001 )
          casaflux%frac_sapwood(POP%Iwood) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
          casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6_r_2)
          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
          where (pop%pop_grid(:)%LU==2)
            casaflux%kplant(POP%Iwood,2) =  1.0 -                   &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality +   &
              POP%pop_grid(:)%crowding_mortality                    &
              + POP%pop_grid(:)%fire_mortality )                    &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
              1.0/veg%disturbance_interval(POP%Iwood,1), 0.99_r_2), 0.0_r_2))**(1./365.)
          elsewhere
            casaflux%kplant(POP%Iwood,2) =  1.0 -                               &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality +               &
              POP%pop_grid(:)%crowding_mortality                                &
              + POP%pop_grid(:)%fire_mortality+POP%pop_grid(:)%cat_mortality  ) &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99_r_2), 0.0_r_2))**(1./365.)
          end where
          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
        elsewhere
          casaflux%frac_sapwood(POP%Iwood) = 1.
          casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6_r_2)
          casaflux%kplant(POP%Iwood,2) = 0.
          veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
        end where
      end if  
      call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
      call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
      call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
      do j = 1,mlitter
        casaflux%klitter(:,j) = casaflux%klitter(:,j)*xkNlimiting
      end do
      call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      call casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      call casa_delplant(veg,casabiome,casapool,casaflux,casamet, &
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,   &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,   &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
      casaflux%Cplant_turnover_disturbance = 0
      casaflux%Cplant_turnover_crowding = 0
      casaflux%Cplant_turnover_resource_limitation = 0
      if ( cable_pop==1 ) then
        if (.not.allocated(tmp)) allocate(tmp(size(POP%pop_grid)))
        tmp = (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality &
             +POP%pop_grid(:)%cat_mortality &
             + POP%pop_grid(:)%fire_mortality  )
        where ( tmp>1.e-12 )
           casaflux%Cplant_turnover_disturbance(POP%Iwood) =                         &
                casaflux%Cplant_turnover(POP%Iwood,2)*(POP%pop_grid(:)%cat_mortality &
                + POP%pop_grid(:)%fire_mortality  )/tmp
           casaflux%Cplant_turnover_crowding(POP%Iwood) =                                     &
                casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%crowding_mortality/tmp
           casaflux%Cplant_turnover_resource_limitation(POP%Iwood) =                          &
                casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%stress_mortality/tmp
        endwhere
      end if
      call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
      call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet,lalloc)
      call casa_cnpbal(casapool,casaflux,casabal)
      casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp*deltpool
      casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp*deltpool
      casabal%FCrmleafyear = casabal%FCrmleafyear + casaflux%Crmplant(:,leaf)*deltpool
      casabal%FCrmwoodyear = casabal%FCrmwoodyear + casaflux%Crmplant(:,wood)*deltpool
      casabal%FCrmrootyear = casabal%FCrmrootyear + casaflux%Crmplant(:,froot)*deltpool
      casabal%FCrgrowyear  = casabal%FCrgrowyear  + casaflux%Crgplant*deltpool
      casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp*deltpool
      casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil*deltpool
      casabal%FCneeyear = casabal%FCneeyear + (casaflux%Cnpp-casaflux%Crsoil)*deltpool
      casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool
      casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
      casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
      casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
      casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
      casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
      casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
      casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
      casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
      casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
      casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
      casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool  
      casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool 
      ! reset casamet for next call
      casamet%tairk = 0.
      casamet%tsoil = 0.
      casamet%moist = 0.
      casaflux%cgpp = 0.
      casaflux%crmplant(:,leaf) = 0.
    end if
    canopy%frp  = real((casaflux%crmplant(:,wood)+casaflux%crmplant(:,xroot)+casaflux%crgplant(:))/real(casaperiod,8))
    canopy%frs  = real(casaflux%Crsoil(:)/real(casaperiod,8))
    canopy%frpw = real(casaflux%crmplant(:,wood)/real(casaperiod,8))
    canopy%frpr = real(casaflux%crmplant(:,xroot)/real(casaperiod,8))
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnee = real((casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/real(casaperiod,8))
    
    sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dt
    sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dt
    sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dt
    sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dt
    sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dt
    sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dt
    sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dt
    sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dt
    sum_flux%sumrs  = sum_flux%sumrs  + canopy%frs*dt
  case default
    write(6,*) "ERROR: Unsupported carbon cycle option with ccycle=",icycle
    stop
end select  

!--------------------------------------------------------------
! POP
if ( cable_pop==1 ) then
  ! update once per year
  ! call popdriver(casaflux,casabal,veg,pop)  
  ! casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
  ! casabal%LAImax = max(casamet%glai, casabal%LAImax)
  ! casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
  ! casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/real(LOY)/1000.
else
  !casaflux%stemnpp = 0.  
end if
  
!--------------------------------------------------------------
      
! Unpack tiles into grid point averages.
! Note that albsav and albnirsav are the VIS and NIR albedo output from CABLE to
! be used by the radiadiation scheme at the next time step.  albvisnir(:,1) and
! albvisnir(:,2) are the VIS and NIR albedo used by the radiation scheme for the
! current time step.
do k = 1,ms
  where ( land(1:ifull) )
    tgg(:,k) = 0.
    wb(:,k) = 0.
    wbice(:,k) = 0.
  end where
end do
do k = 1,3
  where ( land(1:ifull) )
    tggsn(:,k) = 0.
    smass(:,k) = 0.
    ssdn(:,k) = 0.
  end where
end do
where ( land(1:ifull) )
  albvisdir = 0.
  albvisdif = 0.
  albnirdir = 0.
  albnirdif = 0.
  fg = 0.
  eg = 0.
  ga = 0.
  rnet = 0.
  epot = 0.
  tss = 0.
  zo = 0.
  zoh = 0.
  cduv = 0.
  cdtq = 0.
  ustar = 0.
  wetfac = 0.
  rsmin = 0.
  ssdnn = 0.
  snowd = 0.
  snage = 0.
  ! screen and 10m diagnostics - rhscrn calculated in sflux.f
  !tscrn = 0.
  !uscrn = 0.
  !qgscrn = 0.
  !u10 = 0.
end where
tmps = 0. ! average isflag

do nb = 1,maxnb
  is = pind(nb,1)
  ie = pind(nb,2)
  ! albedo
  albvisdir = albvisdir + unpack(sv(is:ie)*rad%reffbm(is:ie,1),tmap(:,nb),0.)
  albnirdir = albnirdir + unpack(sv(is:ie)*rad%reffbm(is:ie,2),tmap(:,nb),0.)
  albvisdif = albvisdif + unpack(sv(is:ie)*rad%reffdf(is:ie,1),tmap(:,nb),0.)
  albnirdif = albnirdif + unpack(sv(is:ie)*rad%reffdf(is:ie,2),tmap(:,nb),0.)
  ! fluxes
  fg = fg + unpack(sv(is:ie)*canopy%fh(is:ie),tmap(:,nb),0.)
  eg = eg + unpack(sv(is:ie)*canopy%fe(is:ie),tmap(:,nb),0.)
  ga = ga + unpack(sv(is:ie)*canopy%ga(is:ie),tmap(:,nb),0.)
  rnet = rnet + unpack(sv(is:ie)*canopy%rnet(is:ie),tmap(:,nb),0.)
  tss = tss + unpack(sv(is:ie)*rad%trad(is:ie)**4,tmap(:,nb),0.) ! ave longwave radiation
  ! drag and mixing
  zo   = zo   + unpack(sv(is:ie)/log(zmin/rough%z0m(is:ie))**2,tmap(:,nb),0.)
  cduv = cduv + unpack(sv(is:ie)*canopy%cduv(is:ie),tmap(:,nb),0.)
  cdtq = cdtq + unpack(sv(is:ie)*canopy%cdtq(is:ie),tmap(:,nb),0.)
  ! soil
  do k = 1,ms
    tgg(:,k)   = tgg(:,k)   + unpack(sv(is:ie)*ssnow%tgg(is:ie,k),        tmap(:,nb),0.)
    wb(:,k)    = wb(:,k)    + unpack(sv(is:ie)*real(ssnow%wb(is:ie,k)),   tmap(:,nb),0.)
    wbice(:,k) = wbice(:,k) + unpack(sv(is:ie)*real(ssnow%wbice(is:ie,k)),tmap(:,nb),0.)
  end do
  ! hydrology
  runoff = runoff + unpack(sv(is:ie)*ssnow%runoff(is:ie)*dt,tmap(:,nb),0.) ! convert mm/s to mm
  runoff_surface = runoff_surface + unpack(sv(is:ie)*ssnow%rnof1(is:ie)*dt,tmap(:,nb),0.) ! convert mm/s to mm
  fwet = fwet + unpack(sv(is:ie)*canopy%fwet(is:ie),tmap(:,nb),0.)         ! used for aerosols
  wetfac = wetfac + unpack(sv(is:ie)*ssnow%wetfac(is:ie),tmap(:,nb),0.)    ! used for aerosols
  cansto = cansto + unpack(sv(is:ie)*canopy%cansto(is:ie),tmap(:,nb),0.)   ! not used
  ! snow
  tmps = tmps + unpack(sv(is:ie)*real(ssnow%isflag(is:ie)),tmap(:,nb),0.)  ! used in radiation (for nsib==3)
  do k = 1,3
    tggsn(:,k) = tggsn(:,k) + unpack(sv(is:ie)*ssnow%tggsn(is:ie,k),tmap(:,nb),0.) ! for restart file
    smass(:,k) = smass(:,k) + unpack(sv(is:ie)*ssnow%smass(is:ie,k),tmap(:,nb),0.) ! for restart file
    ssdn(:,k)  = ssdn(:,k)  + unpack(sv(is:ie)*ssnow%ssdn(is:ie,k),tmap(:,nb),0.)  ! for restart file
  end do
  ssdnn = ssdnn + unpack(sv(is:ie)*ssnow%ssdnn(is:ie),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snage = snage + unpack(sv(is:ie)*ssnow%snage(is:ie),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snowd = snowd + unpack(sv(is:ie)*ssnow%snowd(is:ie),tmap(:,nb),0.)
  snowmelt = snowmelt + unpack(sv(is:ie)*ssnow%smelt(is:ie),tmap(:,nb),0.)
  ! diagnostic
  epot = epot + unpack(sv(is:ie)*ssnow%potev(is:ie),tmap(:,nb),0.)         ! diagnostic in history file
  vlai = vlai + unpack(sv(is:ie)*veg%vlai(is:ie),tmap(:,nb),0.)
  rsmin = rsmin + unpack(sv(is:ie)*canopy%gswx_T(is:ie),tmap(:,nb),0.)     ! diagnostic in history file
  !tscrn = tscrn + unpack(sv(pind(nb,1):pind(nb,2))*canopy%tscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
  !uscrn = uscrn + unpack(sv(pind(nb,1):pind(nb,2))*canopy%uscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
  !qgscrn = qgscrn + unpack(sv(pind(nb,1):pind(nb,2))*canopy%qscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
end do

if ( icycle==0 ) then
  !cplant = 0.
  !csoil = 0.
  !fnee = 0.
  !fpn = 0.
  !frd = 0.
  !frp = 0.
  !frpw = 0.
  !frpr = 0.
  !frs = 0.
  !do nb = 1,maxnb
  !  is = pind(nb,1)
  !  ie = pind(nb,2)
  !  do k = 1,ncp
  !    cplant(:,k) = cplant(:,k) + unpack(sv(is:ie)*bgc%cplant(is:ie,k),tmap(:,nb),0.)
  !  end do
  !  do k = 1,ncs
  !    csoil(:,k) = csoil(:,k) + unpack(sv(is:ie)*bgc%csoil(is:ie,k),   tmap(:,nb),0.)
  !  end do
  !end do
  ! carbon cycle
  !fnee = fnee + unpack(sv(is:ie)*canopy%fnee(is:ie), tmap(:,nb),0.)
  !fpn  = fpn  + unpack(sv(is:ie)*canopy%fpn(is:ie),  tmap(:,nb),0.)
  !frd  = frd  + unpack(sv(is:ie)*canopy%frday(is:ie),tmap(:,nb),0.)
  !frp  = frp  + unpack(sv(is:ie)*canopy%frp(is:ie),  tmap(:,nb),0.)
  !frpw = frpw + unpack(sv(is:ie)*canopy%frpw(is:ie), tmap(:,nb),0.)
  !frs  = frs  + unpack(sv(is:ie)*canopy%frs(is:ie),  tmap(:,nb),0.)
else
  cplant = 0.
  niplant = 0.
  pplant = 0.
  clitter = 0.
  nilitter = 0.
  plitter = 0.
  csoil = 0.
  nisoil = 0.
  psoil = 0.
  !glai = 0.
  fnee = 0.
  fpn = 0.
  frd = 0.
  frp = 0.
  frpw = 0.
  frpr = 0.
  frs = 0.
  do nb = 1,maxnb
    is = pind(nb,1)
    ie = pind(nb,2)
    do k = 1,mplant
      cplant(:,k)  = cplant(:,k)  + unpack(sv(is:ie)*real(casapool%cplant(is:ie,k)),tmap(:,nb),0.)
      niplant(:,k) = niplant(:,k) + unpack(sv(is:ie)*real(casapool%nplant(is:ie,k)),tmap(:,nb),0.)
      pplant(:,k)  = pplant(:,k)  + unpack(sv(is:ie)*real(casapool%pplant(is:ie,k)),tmap(:,nb),0.)
    end do
    do k = 1,mlitter
      clitter(:,k)  = clitter(:,k)  + unpack(sv(is:ie)*real(casapool%clitter(is:ie,k)),tmap(:,nb),0.)
      nilitter(:,k) = nilitter(:,k) + unpack(sv(is:ie)*real(casapool%nlitter(is:ie,k)),tmap(:,nb),0.)
      plitter(:,k)  = plitter(:,k)  + unpack(sv(is:ie)*real(casapool%plitter(is:ie,k)),tmap(:,nb),0.)
    end do
    do k = 1,msoil
      csoil(:,k)  = csoil(:,k)  + unpack(sv(is:ie)*real(casapool%csoil(is:ie,k)),tmap(:,nb),0.)
      nisoil(:,k) = nisoil(:,k) + unpack(sv(is:ie)*real(casapool%nsoil(is:ie,k)),tmap(:,nb),0.)
      psoil(:,k)  = psoil(:,k)  + unpack(sv(is:ie)*real(casapool%psoil(is:ie,k)),tmap(:,nb),0.)
    end do
    !glai = glai + unpack(sv(is:ie)*real(casamet%glai(is:ie)),tmap(:,nb),0.)
    ! carbon cycle
    fnee = fnee + unpack(sv(is:ie)*canopy%fnee(is:ie), tmap(:,nb),0.)
    fpn  = fpn  + unpack(sv(is:ie)*canopy%fpn(is:ie),  tmap(:,nb),0.)
    frd  = frd  + unpack(sv(is:ie)*canopy%frday(is:ie),tmap(:,nb),0.)
    frp  = frp  + unpack(sv(is:ie)*canopy%frp(is:ie),  tmap(:,nb),0.)
    frpw = frpw + unpack(sv(is:ie)*canopy%frpw(is:ie), tmap(:,nb),0.)
    frs  = frs  + unpack(sv(is:ie)*canopy%frs(is:ie),  tmap(:,nb),0.)
  end do
end if

! MJT notes - ustar, cduv, fg and eg are passed to the boundary layer turbulence scheme
! zoh, zoq and zo are passed to the scrnout diagnostics routines
! rsmin is typically used by CTM

where ( land(1:ifull) )
  zo        = zmin*exp(-1./sqrt(zo))
  zoh       = zo/7.4
  zoq       = zoh
  ustar     = sqrt(cduv)*vmod  
  cduv      = cduv*vmod           ! cduv is Cd*vmod in CCAM
  cdtq      = cdtq*vmod
  tss       = tss**0.25
  rsmin     = 1./rsmin
  ! update albedo and tss before calculating net radiation
  albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif
  albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif  
  alb       = swrsave*albvissav + (1.-swrsave)*albnirsav
  !rnet      = sgsave - rgsave - stefbo*tss**4
  isflag(:) = nint(tmps(:)) ! tmps is average isflag
end where
qsttg_land(:) = qsat(ps(1:ifull),tss(1:ifull)) ! must wait for tss to be updated first
where ( land(1:ifull) )
  qsttg(:)  = qsttg_land(:)
end where

return
end subroutine sib4

! *************************************************************************************
subroutine setco2for(atmco2)
! set co2 forcing for cable
! constant: atmospheric co2 = 360 ppm 
! host: atmospheric co2 follows that from CCAM radiation scheme
! interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)

use cc_mpi, only : myid
use newmpar_m
use parm_m
use radisw_m, only : rrvco2
use tracers_m, only : tr,ngas,tractype,tracname

implicit none

integer ico2,igas
real, dimension(ifull), intent(out) :: atmco2

ico2=0
atmco2 = 1.E6*rrvco2          ! from radiative CO2 forcings
if ( tracerco2==1 ) then
  do igas=1,ngas
    if ( trim(tractype(igas))=='online' .and. trim(tracname(igas))=='cbmnep' ) then
      ico2=igas
      exit
    end if
  end do
  if ( ico2>0 ) then
    atmco2 = tr(1:ifull,1,ico2) ! use interactive tracers
  end if
end if
if ( myid==0 .and. ktau==1 ) then
  if ( ico2==0 ) then
    write(6,*) "CABLE using prescribed CO2 from radiative forcings"
  else
    write(6,*) "CABLE using prognostic CO2 from tracer"
  end if
end if

return
end subroutine setco2for

! *************************************************************************************
subroutine cbmemiss(trsrc,mvegt,mode)
  
use newmpar_m
use parm_m

implicit none
  
integer, intent(in) :: mvegt,mode
integer nb
real, dimension(ifull), intent(out) :: trsrc
real, dimension(ifull) :: fpn,frd,frp,frs
  
if ( nsib/=6 .and. nsib/=7 ) then
  write(6,*) "ERROR: Attempted to read CABLE emissions with CABLE disabled"
  stop
end if

fpn=0.
frd=0.
frp=0.
frs=0.
do nb=1,maxnb
  where ( veg%iveg(pind(nb,1):pind(nb,2))==mvegt )
    fpn=fpn+unpack(sv(pind(nb,1):pind(nb,2))*canopy%fpn(pind(nb,1):pind(nb,2)),  tmap(:,nb),0.)
    frd=frd+unpack(sv(pind(nb,1):pind(nb,2))*canopy%frday(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
    frp=frp+unpack(sv(pind(nb,1):pind(nb,2))*canopy%frp(pind(nb,1):pind(nb,2)),  tmap(:,nb),0.)
    frs=frs+unpack(sv(pind(nb,1):pind(nb,2))*canopy%frs(pind(nb,1):pind(nb,2)),  tmap(:,nb),0.)
  end where
end do
  
select case( mode )
  case(1)
    trsrc=fpn-frd
  case(2)
    trsrc=frp+frd
  case(3)
    trsrc=frs
  case default
    write(6,*) "ERROR: Unknown mode for cbmemiss ",mode
    stop
end select
  
return
end subroutine cbmemiss

! *************************************************************************************
subroutine setlai(sigmf,jyear,jmonth,jday,jhour,jmin,mp)

use cc_mpi
use dates_m
use newmpar_m
use parm_m

implicit none
  
integer, intent(in) :: jyear,jmonth,jday,jhour,jmin,mp
integer monthstart,nb
integer, dimension(12) :: imonth
real, parameter :: vextkn = 0.4
real, dimension(ifull), intent(out) :: sigmf
real, dimension(mp) :: a0, a1, a2, aa, bb, cc, mp1, mp2, c2, c3, c4
real x

select case( proglai )
  case(-1) ! piece-wise linear interpolated LAI
    imonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    if ( leap==1 ) then
      if ( mod(jyear,4)  ==0 ) imonth(2)=29
      if ( mod(jyear,100)==0 ) imonth(2)=28
      if ( mod(jyear,400)==0 ) imonth(2)=29
    end if
    monthstart = 1440*(jday-1) + 60*jhour + jmin ! mins from start month
    x = min(max(real(mtimer+monthstart)/real(1440*imonth(jmonth)),0.),1.)
    a0(:) = 0.5*vl1(:)
    a1(:) = -vl2(:)
    a2(:) = 0.5*vl3(:)
    aa(:) = a0(:) + a1(:) + a2(:)
    bb(:) = -3.*a0(:) - 2.*a1(:) - a2(:)
    cc(:) = 2.*a0(:)
    mp1(:) = 0.25 *aa(:) + 0.5*bb(:) + cc(:) ! start of month value
    a0(:) = 0.5*vl2(:)
    a1(:) = -vl3(:)
    a2(:) = 0.5*vl4(:)
    aa(:) = a0(:) + a1(:) + a2(:)
    bb(:) = -3.*a0(:) - 2.*a1(:) - a2(:)
    cc(:) = 2.*a0(:)
    mp2(:) = 0.25 *aa(:) + 0.5*bb(:) + cc(:) ! end of month value
    if ( x<0.5 ) then
      c4(:) = 2.*vl2(:) - 0.5*mp1(:) - 0.5*mp2(:) ! mid-point value
      c2(:) = mp1(:)                              ! intercept
      c3(:) = 2.*(c4(:)-c2(:))                    ! gradient
      veg%vlai(:) = c3(:)*x + c2(:)
    else
      c4(:) = 2.*vl2(:) - 0.5*mp1(:) - 0.5*mp2(:) ! mid-point value
      c3(:) = 2.*(mp2(:)-c4(:))                   ! gradient
      c2(:) = 2.*c4(:) - mp2(:)                   ! intercept
      veg%vlai(:) = c3(:)*x + c2(:)
    end if
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.01 )
    elsewhere
      veg%vlai = 1.E-8
    end where
    
  case(0) ! PWCB interpolated LAI
    imonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    if ( leap==1 ) then
      if ( mod(jyear,4)  ==0 ) imonth(2) = 29
      if ( mod(jyear,100)==0 ) imonth(2) = 28
      if ( mod(jyear,400)==0 ) imonth(2) = 29
    end if
    monthstart = 1440*(jday-1) + 60*jhour + jmin ! mins from start month
    x = min(max(real(mtimer+monthstart)/real(1440*imonth(jmonth)),0.),1.)
    aa(:) = vl1(:)
    bb(:) = vl2(:) + vl1(:)
    cc(:) = vl3(:) + vl2(:) + vl1(:)
    a0(:) = 0.5*bb(:)
    a1(:) = 4.*bb(:) - 5.*aa(:) - cc(:)
    a2(:) = 1.5*(cc(:) + 3.*aa(:) - 3.*bb(:))
    veg%vlai = a0(:)+a1(:)*x+a2(:)*x*x     ! LAI as a function of time
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.01 )
    elsewhere
      veg%vlai = 1.E-8
    end where

  case(1) ! prognostic LAI
    if ( icycle/=3 ) then
      write(6,*) "ERROR: CASA CNP LAI is not operational"
      write(6,*) "Prognostic LAI requires ccycle=3"
      call ccmpi_abort(-1)
    end if
    veg%vlai(:) = real(casamet%glai(:))
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.01 )
    elsewhere
      veg%vlai = 1.E-8
    end where

  case default
    write(6,*) "ERROR: Unknown proglai option ",proglai
    call ccmpi_abort(-1)
end select

! diagnose greeness fraction (e.g., for aerosols)  
sigmf(:) = 0.
do nb = 1,maxnb
  sigmf(:) = sigmf(:) + unpack(sv(pind(nb,1):pind(nb,2))*(1.-exp(-vextkn*veg%vlai(pind(nb,1):pind(nb,2)))),tmap(:,nb),0.)
end do
sigmf = min( sigmf, 1. )
  
return
end subroutine setlai

! *************************************************************************************
! Calculate prognostic vcmax from CASA-CNP
subroutine vcmax_feedback

use newmpar_m, only : mxvt

implicit none

integer np, ivt
real, dimension(mp) :: ncleafx, npleafx, pleafx, nleafx
real, dimension(mxvt), parameter :: xnslope = (/ 0.8,1.,2.,1.,1.,1.,0.5,1.,0.34,1.,1.,1.,1.,1.,1.,1.,1. /)

if ( progvcmax>0 .and. icycle==3 ) then

  ! initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx(:) = 14.2

  do np = 1,mp
    ivt = veg%iveg(np)
    if ( casamet%iveg2(np)/=icewater .and. casamet%glai(np)>casabiome%glaimin(ivt) .and. &
         casapool%cplant(np,leaf)>0. ) then
      ncleafx(np) = min( casabiome%ratioNCplantmax(ivt,leaf),                  &
                    max( casabiome%ratioNCplantmin(ivt,leaf),                  &
                         casapool%nplant(np,leaf)/casapool%cplant(np,leaf) ) )
      if ( casapool%pplant(np,leaf)>0. ) then
        npleafx(np) = min(30.,max(8.,real(casapool%nplant(np,leaf)/casapool%pplant(np,leaf))))
      end if
    end if
  end do

  ! method for prognostic vcmax
  select case ( progvcmax )
    case(1)
      do np = 1,mp
        ivt = veg%iveg(np)  
        if ( casamet%glai(np)>casabiome%glaimin(ivt) ) then
          if ( ivt/=2 ) then
            veg%vcmax(np) = ( casabiome%nintercept(ivt)                            &
                    + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.e-6
          else
            if ( casapool%nplant(np,leaf)>0. .and. casapool%pplant(np,leaf)>0. ) then
              veg%vcmax(np) = ( casabiome%nintercept(ivt)        &
                   + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                   * ncleafx(np)/casabiome%sla(ivt) )*1.e-6
            else
              veg%vcmax(np) = ( casabiome%nintercept(ivt)                           &
                   + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.e-6
            end if
          end if
          veg%vcmax(np) =veg%vcmax(np)*xnslope(ivt)
        end if
      end do
  
    case(2)
      !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax – 
      !to leaf nitrogen, leaf phosphorus, and specific leaf area: 
      !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
      do np = 1,mp
        ivt = veg%iveg(np)  
        if ( ivt==7 ) then
          veg%vcmax(np) = 1.e-5 ! special for C4 grass: set here to value from  parameter file
        else
          nleafx(np) = ncleafx(np)/casabiome%sla(ivt) ! leaf N in g N m-2 leaf
          pleafx(np) = nleafx(np)/npleafx(np)         ! leaf P in g P m-2 leaf
          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))
        end if
      end do
 
    case default
      write(6,*) "ERROR: Invalid progvcmax ",progvcmax
      stop
  end select

  ! update ejmax  
  veg%ejmax = 2.*veg%vcmax

end if
  
return
end subroutine vcmax_feedback

! *************************************************************************************
! Uphade phenology depending on climate
subroutine cable_phenology_clim

implicit none

integer :: np, days, ivt
real :: gdd0
real :: phengdd5ramp
real(r_2) :: phen_tmp
real, parameter :: k_chilla = 0, k_chillb = 100, k_chillk = 0.05
real, parameter :: APHEN_MAX = 200.0, mmoisture_min=0.30

!phen%doyphase(np,1) ! DOY for greenup
!phen%doyphase(np,2) ! DOY for steady LAI
!phen%doyphase(np,3) ! DOY for leaf senescence
!phen%doyphase(np,4) ! DOY for minimal LAI season

do np = 1,mp
  ivt = veg%iveg(np)  
  phen_tmp = 0.
  
  ! evergreen pfts
  if ( ivt==31 .or. ivt==2 .or. ivt==5 ) then
    phen%doyphase(np,1) = -50
    phen%doyphase(np,2) = phen%doyphase(np,1) + 14
    phen%doyphase(np,3) = 367
    phen%doyphase(np,4) = phen%doyphase(np,3) + 14
    phen%phase(np) = 2
  end if

  ! summergreen woody pfts
  if ( ivt==3 .or. ivt==4 ) then  ! deciduous needleleaf(3) and broadleaf(4)
    ! Calculate GDD0  base value (=gdd to bud burst) for this PFT given
    !  current length of chilling period (Sykes et al 1996, Eqn 1)
    gdd0 = k_chilla + k_chillb*exp(-k_chillk*real(climate%chilldays(np)))
    phengdd5ramp = 200
    if ( climate%gdd5(np)>gdd0 .and. phen%aphen(np)<APHEN_MAX) then
      phen_tmp = min(1., (climate%gdd5(np)-gdd0)/phengdd5ramp)
    end if
  end if

  ! summergreen grass or crops
  if ( ivt>=6 .and. ivt<=10 )  then     ! grass or crops
    phengdd5ramp = 50
    phen_tmp = min(1., climate%gdd5(np)/phengdd5ramp)
  end if

  ! raingreen pfts
  if ( ivt>=6 .and. ivt<=10 ) then ! (grass or crops) need to include raingreen savanna trees here too
    if (climate%dmoist(np)<mmoisture_min) then
      phen_tmp = 0.
    end if  
  end if

  if ( (ivt==3.or.ivt==4) .or. (ivt>=6.and.ivt<=10) ) then

    if (phen_tmp>0. .and.( phen%phase(np)==3 .or. phen%phase(np)==0 )) then
      phen%phase(np) = 1 ! greenup
      phen%doyphase(np,1) = climate%doy
    elseif (phen_tmp>=1._r_2 .and. phen%phase(np)==1) then
      phen%phase(np) = 2 ! steady LAI
      phen%doyphase(np,2) = climate%doy
    elseif (phen_tmp<1._r_2 .and. phen%phase(np)==2) then
      phen%phase(np) = 3 ! senescence
      phen%doyphase(np,3) = climate%doy
    endif

    if ( phen%phase(np)==3 ) then
       days = min(climate%doy,365)-phen%doyphase(np,3)
       if ( days<0 ) days = days + 365
       if ( days>14 ) phen%phase(np) = 0          ! mimimum LAI
    endif

    ! Update annual leaf-on sum
    if ( (rad%latitude(np)>=0. .and. climate%doy==COLDEST_DAY_NHEMISPHERE) .or. &
         (rad%latitude(np)<0. .and. climate%doy==COLDEST_DAY_SHEMISPHERE) ) then
      phen%aphen(np) = 0
    end if
    phen%phen(np) = phen_tmp
    phen%aphen(np) = phen%aphen(np) + phen%phen(np)

  end if

end do  ! end loop over patches

return
end subroutine cable_phenology_clim

! *************************************************************************************
! track climate feedback into CABLE
subroutine cableclimate(idoy,imonth,iyear)

use parm_m, only : dt, leap

implicit none

real, dimension(mp) :: RhoA, PPc, EpsA, phiEq
real, dimension(mp) :: mtemp_last !, tmp
real, parameter :: Gaero  = 0.015    ! (m s-1) aerodynmaic conductance (for use in PT evap)
real, parameter :: Capp   = 29.09    ! isobaric spec heat air    [J/molA/K]
real, parameter :: SBoltz = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
real, parameter :: CoeffPT = 1.26
integer, intent(in) :: idoy, imonth,iyear
integer, save :: climate_daycount = 0  ! counter for daily averages
integer d
integer, dimension(12) :: ndoy
integer, dimension(12), parameter :: odoy=(/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /) 

if ( cable_climate==0 ) return

climate%doy = idoy

! accumulate annual evaporation and potential evaporation
!RhoA = met%pmb*100./(8.314*met%tk) ! air density [molA/m3]
!PPc  = Gaero/(Gaero+4.*SBoltz*(met%tk**3/(RhoA*Capp)))
!tmp = met%tk - 273.16 ! avoid array temporary
!EpSA  = Epsif(tmp, met%pmb)
!phiEq = canopy%rniso*(PPc*EpsA)/(PPc*EpsA+1.) ! equil ltnt heat flux  [W/m2]

! accumulate daily temperature, evap and potential evap
if ( real(climate_daycount)*dt>86400. ) then
  climate%dtemp   = 0.
  climate%dmoist  = 0.
  climate_daycount = 0    
end if
!climate%evap_PT = climate%evap_PT + max(phiEq,1.)*CoeffPT/air%rlam*dt  ! mm
!climate%aevap   = climate%aevap + met%precip ! mm
climate%dtemp  = climate%dtemp + met%tk - 273.15
climate%dmoist = climate%dmoist + canopy%fwsoil
climate_daycount = climate_daycount + 1

if ( real(climate_daycount)*dt>86400. ) then
  climate%dtemp = climate%dtemp/real(climate_daycount)
  climate%dmoist = climate%dmoist/real(climate_daycount)

  ! On first day of year ...
  if ( idoy==1 ) then
    climate%agdd5 = 0.
    climate%agdd0 = 0.
    !climate%evap_PT = 0.     ! annual PT evap [mm]
    !climate%aevap  = 0.      ! annual evap [mm]  
  end if
  ! In midwinter, reset GDD counter for summergreen phenology  
  where ( (rad%latitude>=0..and.idoy==COLDEST_DAY_NHEMISPHERE) .or. &
          (rad%latitude<0..and.idoy==COLDEST_DAY_SHEMISPHERE) )
    climate%gdd5 = 0.
    climate%gdd0 = 0.
  end where
  ! Update GDD counters and chill day count
  climate%gdd0  = climate%gdd0 + max(0.,climate%dtemp)
  climate%agdd0 = climate%agdd0 + max(0.,climate%dtemp)
  climate%gdd5  = climate%gdd5 + max(0.,climate%dtemp-5.)
  climate%agdd5 = climate%agdd5 + max(0.,climate%dtemp-5.)
  where ( climate%dtemp<5. .and. climate%chilldays<=365 )
    climate%chilldays = climate%chilldays + 1
  end where

  ! Save yesterday's mean temperature for the last month
  mtemp_last = climate%mtemp

  ! Update daily temperatures, and mean overall temperature, for last 31 days
  climate%mtemp = climate%dtemp
  climate%qtemp = climate%dtemp
  climate%mmoist = climate%dmoist
  do d = 1,30
    climate%dtemp_31(:,d) = climate%dtemp_31(:,d+1) ! why not use dtemp_91(:,61:91)?
    climate%mtemp = climate%mtemp + climate%dtemp_31(:,d)
    climate%dmoist_31(:,d) = climate%dmoist_31(:,d+1)
    climate%mmoist = climate%mmoist + climate%dmoist_31(:,d)
  end do
  do d = 1,90
    climate%dtemp_91(:,d) = climate%dtemp_91(:,d+1)
    climate%qtemp = climate%qtemp + climate%dtemp_91(:,d)
  end do
  climate%dtemp_91(:,91) = climate%dtemp
  climate%dtemp_31(:,31) = climate%dtemp
  climate%dmoist_31(:,31) = climate%dmoist
  climate%qtemp = climate%qtemp/91.   ! average temperature over the last quarter
  climate%mtemp = climate%mtemp/31.   ! average temperature over the last month
  climate%mmoist = climate%mmoist/31. ! average moisture index over the last month

  ! Reset GDD and chill day counter if mean monthly temperature falls below base
  ! temperature
  where ( mtemp_last>=5. .and. climate%mtemp<5. )
    climate%gdd5 = 0.
    climate%chilldays = 0
  end where

  ! calculate end of the month
  ndoy(1:12) = odoy(1:12)
  if ( leap==1 ) then
    if ( mod(iyear,4)  ==0 ) ndoy(3:12)=odoy(3:12)+1
    if ( mod(iyear,100)==0 ) ndoy(3:12)=odoy(3:12)
    if ( mod(iyear,400)==0 ) ndoy(3:12)=odoy(3:12)+1
  end if
  
  ! check if end of month
  if ( idoy==ndoy(imonth) ) then

    ! Update mean temperature for the last 12 months
    ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)
    !climate%atemp_mean = climate%atemp_mean*(11./12.) + climate%mtemp*(1./12.)

    ! Record minimum and maximum monthly temperatures
    if ( imonth==1 ) then
      climate%mtemp_min = climate%mtemp;
      climate%mtemp_max = climate%mtemp;
      climate%qtemp_max_last_year = climate%qtemp_max;
      climate%qtemp_max = climate%qtemp;
    end if
    climate%mtemp_min = min(climate%mtemp, climate%mtemp_min)
    climate%mtemp_max = max(climate%mtemp, climate%mtemp_max)
    climate%qtemp_max = max(climate%qtemp, climate%qtemp_max)

    ! On 31 December update records of minimum monthly temperatures for the last
    ! 20 years and find minimum monthly temperature for the last 20 years
    !if ( imonth==12 ) then
    !  climate%nyears = climate%nyears + 1
    !
    !  startyear=20-min(19,climate%nyears-1)
    !  climate%mtemp_min20=0.0
    !  climate%mtemp_max20=0.0
    !  climate%alpha_PT20 = 0.0
    ! 
    !  if (startyear<20) then
    !    do y = startyear,19
    !      climate%mtemp_min_20(:,y) = climate%mtemp_min_20(:,y+1)
    !      climate%mtemp_min20 = climate%mtemp_min20 + climate%mtemp_min_20(:,y)
    !      climate%mtemp_max_20(:,y) = climate%mtemp_max_20(:,y+1)
    !      climate%mtemp_max20 = climate%mtemp_max20 + climate%mtemp_max_20(:,y)
    !      climate%alpha_PT_20(:,y) = climate%alpha_PT_20(:,y+1)
    !      climate%alpha_PT20 = climate%alpha_PT20 + climate%alpha_PT_20(:,y)
    !    end do
    !    climate%mtemp_min20=climate%mtemp_min20/real(20-startyear)
    !    climate%mtemp_max20=climate%mtemp_max20/real(20-startyear)
    !    climate%alpha_PT20=climate%alpha_PT20/real(20-startyear)
    !  else
    !    ! only occurs when climate%nyears = 1
    !    climate%mtemp_min20 = climate%mtemp_min
    !    climate%mtemp_max20 = climate%mtemp_max
    !    climate%alpha_PT20 = climate%alpha_PT
    !  end if
    ! 
    !  climate%mtemp_min_20(:,20)=climate%mtemp_min
    !  climate%mtemp_max_20(:,20)=climate%mtemp_max
    !  climate%alpha_PT = max(climate%aevap/climate%evap_PT, 0.) ! ratio of annual evap to annual PT evap
    !  climate%alpha_PT_20(:,20)=climate%alpha_PT
    ! 
    !  call biome1_pft
    !
    !end if  ! last month of year
  end if    ! last day of month
end if      ! end of day

return
end subroutine cableclimate

function epsif(tc,pmb) result(epsif_out)

implicit none

real, dimension(mp), intent(in) :: tc, pmb
real, dimension(mp) :: epsif_out
real, dimension(mp) :: tctmp, es, desdt
real,parameter:: A = 6.106      ! Teten coefficients
real,parameter:: B = 17.27      ! Teten coefficients
real,parameter:: C = 237.3      ! Teten coefficients
real,parameter:: Rlat      = 44140.0  ! lat heat evap H2O at 20C  [J/molW]
real,parameter:: Capp      = 29.09    ! isobaric spec heat air    [J/molA/K]

TCtmp = max( min( TC, 100.), -40. )     ! constrain TC to (-40.0,100.0)
ES    = A*EXP(B*TCtmp/(C+TCtmp))        ! sat vapour pressure
dESdT = ES*B*C/(C+TCtmp)**2             ! d(sat VP)/dT: (mb/K)
Epsif_out = (Rlat/Capp)*dESdT/Pmb       ! dimensionless (ES/Pmb = molW/molA)

return
end function epsif

!subroutine biome1_pft
!
!implicit none
!
!integer :: k, npft
!integer, dimension(4) :: pft_biome1
!real :: alpha_pt_scaled
!
!! TABLE 1 , Prentice et al. J. Biogeog., 19, 117-134, 1992
!! pft_biome1: Trees (1)tropical evergreen; (2) tropical raingreen; (3) warm temp evergreen ;
!! (4) temperate summergreen; (5) cool-temperate conifer; (6) boreal evergreen conifer;
!! (7) boreal summergreen;  
!! Non-trees: (8) sclerophyll/succulent; (9) warm grass/shrub; (10) cool grass/shrub; 
!! (11) cold grass/shrub; (12) hot desert shrub; (13) cold desert shrub.
!
!do k = 1,mp
!
!   alpha_PT_scaled =  min(climate%alpha_PT20(k), 1.)
!   pft_biome1(:) = 999
!
!   IF ( climate%mtemp_min20(k)>=15.5 ) THEN
!     IF ( alpha_PT_scaled>=0.80 ) THEN
!       pft_biome1(1) = 1
!       IF ( alpha_PT_scaled<=0.85 ) THEN
!         pft_biome1(2) = 2
!       ENDIF
!     ELSEIF ( alpha_PT_scaled>=0.4 .and. alpha_PT_scaled<0.80 ) THEN
!       pft_biome1(1) = 2
!     END IF
!   END IF
!   
!   IF ( climate%mtemp_min20(k)>=5 .and. alpha_PT_scaled>=0.4 .and. &
!        pft_biome1(1)==999 ) THEN
!     pft_biome1(1) = 3
!   END IF
!   
!   IF ( climate%mtemp_min20(k)>=-15 .and. climate%mtemp_min20(k)<=15.5 .and. &
!        alpha_PT_scaled>=0.35 .and. climate%agdd5(k)>1200 .and.              &
!        pft_biome1(1)>3 ) THEN
!     pft_biome1(1) = 4
!   END IF
!   
!   IF ( climate%mtemp_min20(k)>=-19 .and. climate%mtemp_min20(k)<=5 .and. &
!        alpha_PT_scaled>=0.35 .and. climate%agdd5(k)>900 )  THEN
!     IF (pft_biome1(1)>4 ) THEN
!       pft_biome1(1) = 5
!     ELSEIF (pft_biome1(1)==4 ) THEN
!       pft_biome1(2) = 5
!     END IF
!   END IF
!   
!   IF ( climate%mtemp_min20(k)>=-35 .and. climate%mtemp_min20(k)<=-2 .and. &
!        alpha_PT_scaled>=0.35 .and. climate%agdd5(k)>350 )  THEN
!     IF ( pft_biome1(1)==999 ) THEN
!       pft_biome1(1) = 6
!     ELSEIF ( pft_biome1(2)==999 ) THEN
!       pft_biome1(2) = 6
!     ELSE
!       pft_biome1(3) = 6
!     END IF
!   END IF
!   
!   IF ( climate%mtemp_min20(k)<=5 .and. alpha_PT_scaled>=0.45 .and. &
!        climate%agdd5(k)>350 ) THEN
!     IF (pft_biome1(1)==999 ) THEN
!       pft_biome1(1) = 7
!     ELSEIF (pft_biome1(2)==999 ) THEN
!       pft_biome1(2) = 7
!     ELSEIF (pft_biome1(3)==999 ) THEN
!       pft_biome1(3) = 7  
!     ELSE
!       pft_biome1(4) = 7
!     END IF
!   END IF
!   
!   IF ( climate%mtemp_min20(k)>=5 .and. alpha_PT_scaled>=0.2 .and. &
!        pft_biome1(1)==999 ) THEN
!     pft_biome1(1) = 8
!   END IF
!   
!   IF (climate%mtemp_max20(k)>=22 .and.alpha_PT_scaled>=0.1 &
!        .and. pft_biome1(1).eq.999 ) THEN
!     pft_biome1(1) = 9
!   END IF
!   
!   IF ( climate%agdd5(k)>=500 .and.alpha_PT_scaled>=0.33 .and. &
!        pft_biome1(1)==999 ) THEN
!     pft_biome1(1) = 10
!   END IF
!   
!   IF ( climate%agdd0(k)>=100 .and. alpha_PT_scaled>=0.33 ) THEN 
!     IF ( pft_biome1(1)==999 ) THEN
!       pft_biome1(1) = 11
!     ELSEIF ( pft_biome1(1)==10 ) THEN
!       pft_biome1(2) = 11
!     END IF
!   END IF
!   
!   IF ( climate%mtemp_max20(k)>=22 .and. pft_biome1(1)==999 ) THEN
!     pft_biome1(1) = 12
!   END IF
!   
!   IF ( climate%agdd0(k)>=100 .and. pft_biome1(1)==999 ) THEN
!     pft_biome1(1) = 13
!   END IF
!   
!  ! end of evironmental constraints on pft
!  npft = count( pft_biome1(:)/=999 )
!  
!  ! MAP to Biome1 biome and CABLE pft
!  
!  ! (1) Tropical Rainforest
!  if ( pft_biome1(1)==1 .and. npft==1 ) then
!    climate%biome(k) = 1
!    climate%iveg(k) = 2
!  end if
!
!  ! (2) Tropical Seasonal forest
!  if ( pft_biome1(1)==1 .and. pft_biome1(2)==2 .and. npft==2 ) then
!    climate%biome(k) = 2
!    climate%iveg(k) = 2
!  end if
!
!  ! (3) Tropical dry forest/savanna
!  if ( pft_biome1(1)==2 .and. npft==1 ) then
!    climate%biome(k) = 3
!    climate%iveg(k) = 2  ! N.B. need to include c4 grass
!  end if
!
!  ! (4) Broad-leaved evergreen/warm mixed-forest
!  if ( pft_biome1(1)==3 .and. npft==1 ) then
!    climate%biome(k) = 4
!    climate%iveg(k) = 2
!  end if
!
!  ! (5) Temperate deciduous forest
!  if ( pft_biome1(1)==4 .and. pft_biome1(2)==5 .and. &
!       pft_biome1(3)==7 .and. npft==3 ) then
!    climate%biome(k) = 5
!    climate%iveg(k) = 4
!  end if
!
!  ! (6) Cool mixed forest
!  if ( pft_biome1(1)==4 .and. pft_biome1(2)==5 .and. &
!       pft_biome1(3)==6 .and. pft_biome1(4)==7 .and. &
!       npft==4 ) then
!    climate%biome(k) = 6
!    climate%iveg(k) = 4
!  end if
!
!  ! (7) Cool conifer forest
!  if ( pft_biome1(1)==5 .and. pft_biome1(2)==6 .and. &
!       pft_biome1(3)==7 .and. npft==3 ) then
!    climate%biome(k) = 7
!    climate%iveg(k) = 1
!  end if
!      
!  ! (8) Taiga
!  if ( pft_biome1(1)==6 .and. pft_biome1(2)==7 .and. npft==2 ) then
!    climate%biome(k) = 8
!    climate%iveg(k) = 1
!  end if
!
!  ! (9) Cold mixed forest
!  if ( pft_biome1(1)==5 .and. pft_biome1(2)==7 .and. npft==2 ) then
!    climate%biome(k) = 9
!    climate%iveg(k) = 1
!  end if
!
!  ! (10) Cold deciduous forest
!  if ( pft_biome1(k,1)==7 .and. npft==1 ) then
!    climate%biome(k) = 10
!    climate%iveg(k) = 3
!  end if
!
!  ! (11) Xerophytic woods/scrub
!  if ( pft_biome1(1)==8 .and. npft==1 ) then
!    climate%biome(k) = 11
!    climate%iveg(k) = 5
!  end if
!
!  ! (12) Warm grass/shrub
!  if ( pft_biome1(1)==9 .and. npft==1 ) then
!    climate%biome(k) = 12
!    climate%iveg(k) = 5  ! include C4 grass tile ?
!  end if
!
!  ! (13) Cool grass/shrub
!  if ( pft_biome1(1)==10 .and. pft_biome1(2)==11 .and. npft==2 ) then
!    climate%biome(k) = 13
!    climate%iveg(k) = 5  ! include C3 grass tile ?
!  end if
!
!  ! (14) Tundra
!  if ( pft_biome1(1)==11 .and. npft==1 ) then
!    climate%biome(k) = 14
!    climate%iveg(k) = 8
!  end if
!
!  ! (15) Hot desert
!  if ( pft_biome1(1)==12 .and. npft==1 ) then
!    climate%biome(k) = 15
!    climate%iveg(k) = 14
!  end if
!
!  ! (16) Semidesert
!  if ( pft_biome1(1)==13 .and. npft==1 ) then
!    climate%biome(k) = 16
!    climate%iveg(k) = 5
!  end if
!
!  ! (17) Ice/polar desert
!  if ( climate%biome(k)==999 ) then
!    climate%biome(k) = 17
!    climate%iveg(k) = 17
!  end if
!
!  ! check for DBL or NEL in SH: set to EBL instead
!  if ( (climate%iveg(k)==1.or.climate%iveg(k)==3.or.climate%iveg(k)==4) &
!     .and. rad%latitude(k)<0. ) then
!    climate%iveg(k) = 2
!  end if
!     
!end do
!
!return
!end subroutine biome1_pft

! *************************************************************************************
subroutine loadcbmparm(fveg,fvegprev,fvegnext,fvegnext2,fphen,casafile, &
                       ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,        &
                       casapoint,greenup,fall,phendoy1)

use cc_mpi
use newmpar_m
use parm_m
  
implicit none

integer n
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs,vlin,vlinprev,vlinnext,vlinnext2
real, dimension(ifull,maxtile), intent(out) :: casapoint
real cableformat
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
character(len=*), intent(in) :: fveg,fvegprev,fvegnext,fvegnext2,fphen,casafile

! read CABLE biome and LAI data
if ( myid==0 ) then
  write(6,*) "Reading tiled surface data for CABLE"
  call vegta(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,fvegprev,fveg,fvegnext,fvegnext2, &
             cableformat)
else
  call vegtb(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,fvegprev,fvegnext,fvegnext2, &
             cableformat)
end if
do n = 1,maxtile
  svs(:,n) = svs(:,n)/sum(svs,dim=2)
end do

if ( fvegprev==' ' .and. fvegnext==' ' ) then
  vlinprev = vlin
  vlinnext = vlin
end if
if ( fvegnext2==' ' ) then
  vlinnext2 = vlinnext
end if

if ( abs(cableformat-1.)<1.e-20 ) then
  if ( myid==0 ) write(6,*) "Procesing CSIRO PFTs"    
else
  if ( myid==0 ) write(6,*) "Processing IGBP and converting to CSIRO PFTs"
  call convertigbp(ivs,svs,vlin,vlinprev,vlinnext,vlinnext2)
end if

if ( ccycle==0 ) then
  casapoint(:,:) = 0.
  greenup(:,:)  = -50
  fall(:,:)     = 367
  phendoy1(:,:) = 2  
else
  call casa_readpoint(casafile,casapoint)         ! read point sources
  call casa_readphen(fphen,greenup,fall,phendoy1) ! read MODIS leaf phenology
end if

return
end subroutine loadcbmparm


subroutine convertigbp(ivs,svs,vlin,vlinprev,vlinnext,vlinnext2)

use cc_mpi
use const_phys
use latlong_m
use newmpar_m
use soil_m

implicit none

integer, dimension(ifull,maxtile), intent(inout) :: ivs
integer, dimension(1) :: pos
integer iq, n, ipos, iv
real, dimension(ifull,maxtile), intent(inout) :: svs,vlin,vlinprev,vlinnext,vlinnext2
real, dimension(18,0:3) :: newlai
real, dimension(18) :: newgrid
real fc3, fc4, ftu, fg3, fg4, clat, nsum
real xp

do iq = 1,ifull
  if ( land(iq) ) then
    newgrid(:)  = 0.
    newlai(:,:) = 0.      
    clat = rlatt(iq)*180./pi
    ! grass
    if (abs(clat)>50.5) then
      fg3=0.
      fg4=0.
    else if (abs(clat)>49.5) then
      xp=abs(clat)-49.5
      fg3=(1.-xp)*0.9
      fg4=(1.-xp)*0.1
    else if (abs(clat)>40.5) then
      fg3=0.9
      fg4=0.1
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fg3=(1.-xp)*0.8+xp*0.9
      fg4=(1.-xp)*0.2+xp*0.1
    else if (abs(clat)>30.5) then
      fg3=0.8
      fg4=0.2
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fg3=(1.-xp)*0.5+xp*0.8
      fg4=(1.-xp)*0.5+xp*0.2
    else if (abs(clat)>25.5) then
      fg3=0.5
      fg4=0.5
    else if (abs(clat)>24.5) then
      xp=abs(clat)-24.5
      fg3=(1.-xp)*0.05+xp*0.5
      fg4=(1.-xp)*0.95+xp*0.5
    else
      fg3=0.05
      fg4=0.95
    end if
    ftu=1.-fg3-fg4
    ! crops
    if (abs(clat)>40.5) then
      fc3=1.
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fc3=(1.-xp)*0.9+xp
    else if (abs(clat)>30.5) then
      fc3=0.9
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fc3=(1.-xp)*0.7+xp*0.9
    else
      fc3=0.7
    end if
    fc4=1.-fc3
    do n = 1,maxtile
      select case (ivs(iq,n))
        case (1,2,3,4,11)
          newgrid(ivs(iq,n))=newgrid(ivs(iq,n))+svs(iq,n)
          newlai(ivs(iq,n),0)=newlai(ivs(iq,n),0)+svs(iq,n)*vlinprev(iq,n)
          newlai(ivs(iq,n),1)=newlai(ivs(iq,n),1)+svs(iq,n)*vlin(iq,n)
          newlai(ivs(iq,n),2)=newlai(ivs(iq,n),2)+svs(iq,n)*vlinnext(iq,n)
          newlai(ivs(iq,n),3)=newlai(ivs(iq,n),3)+svs(iq,n)*vlinnext2(iq,n)
        case (5)
          if (abs(clat)>25.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5
            newlai(1,0)=newlai(1,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(1,1)=newlai(1,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(1,2)=newlai(1,2)+svs(iq,n)*0.5*vlinnext(iq,n)
            newlai(1,3)=newlai(1,3)+svs(iq,n)*0.5*vlinnext2(iq,n)
            newgrid(4)=newgrid(4)+svs(iq,n)*0.5
            newlai(4,0)=newlai(4,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(4,1)=newlai(4,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(4,2)=newlai(4,2)+svs(iq,n)*0.5*vlinnext(iq,n)
            newlai(4,3)=newlai(4,3)+svs(iq,n)*0.5*vlinnext2(iq,n)
          else if (abs(clat)>24.5) then
            xp=abs(clat)-24.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5*xp
            newlai(1,0)=newlai(1,0)+svs(iq,n)*0.5*vlinprev(iq,n)*xp
            newlai(1,1)=newlai(1,1)+svs(iq,n)*0.5*vlin(iq,n)*xp
            newlai(1,2)=newlai(1,2)+svs(iq,n)*0.5*vlinnext(iq,n)*xp
            newlai(1,3)=newlai(1,3)+svs(iq,n)*0.5*vlinnext2(iq,n)*xp
            newgrid(4)=newgrid(4)+svs(iq,n)*(1.-0.5*xp)
            newlai(4,0)=newlai(4,0)+svs(iq,n)*vlinprev(iq,n)*(1.-0.5*xp)
            newlai(4,1)=newlai(4,1)+svs(iq,n)*vlin(iq,n)*(1.-0.5*xp)
            newlai(4,2)=newlai(4,2)+svs(iq,n)*vlinnext(iq,n)*(1.-0.5*xp)
            newlai(4,3)=newlai(4,3)+svs(iq,n)*vlinnext2(iq,n)*(1.-0.5*xp)
          else
            newgrid(4)=newgrid(4)+svs(iq,n)
            newlai(4,0)=newlai(4,0)+svs(iq,n)*vlinprev(iq,n)
            newlai(4,1)=newlai(4,1)+svs(iq,n)*vlin(iq,n)
            newlai(4,2)=newlai(4,2)+svs(iq,n)*vlinnext(iq,n)
            newlai(4,3)=newlai(4,3)+svs(iq,n)*vlinnext2(iq,n)
          end if
        case (6)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.8
          newlai(5,0)=newlai(5,0)+svs(iq,n)*0.8*vlinprev(iq,n)
          newlai(5,1)=newlai(5,1)+svs(iq,n)*0.8*vlin(iq,n)
          newlai(5,2)=newlai(5,2)+svs(iq,n)*0.8*vlinnext(iq,n)
          newlai(5,3)=newlai(5,3)+svs(iq,n)*0.8*vlinnext2(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.2*fg3
          newlai(6,0)=newlai(6,0)+svs(iq,n)*0.2*fg3*vlinprev(iq,n)
          newlai(6,1)=newlai(6,1)+svs(iq,n)*0.2*fg3*vlin(iq,n)
          newlai(6,2)=newlai(6,2)+svs(iq,n)*0.2*fg3*vlinnext(iq,n)
          newlai(6,3)=newlai(6,3)+svs(iq,n)*0.2*fg3*vlinnext2(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.2*fg4
          newlai(7,0)=newlai(7,0)+svs(iq,n)*0.2*fg4*vlinprev(iq,n)
          newlai(7,1)=newlai(7,1)+svs(iq,n)*0.2*fg4*vlin(iq,n)
          newlai(7,2)=newlai(7,2)+svs(iq,n)*0.2*fg4*vlinnext(iq,n)
          newlai(7,3)=newlai(7,3)+svs(iq,n)*0.2*fg4*vlinnext2(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.2*ftu
          newlai(8,0)=newlai(8,0)+svs(iq,n)*0.2*ftu*vlinprev(iq,n)
          newlai(8,1)=newlai(8,1)+svs(iq,n)*0.2*ftu*vlin(iq,n)
          newlai(8,2)=newlai(8,2)+svs(iq,n)*0.2*ftu*vlinnext(iq,n)
          newlai(8,3)=newlai(8,3)+svs(iq,n)*0.2*ftu*vlinnext2(iq,n)
        case (7)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.2
          newlai(5,0)=newlai(5,0)+svs(iq,n)*0.2*vlinprev(iq,n)
          newlai(5,1)=newlai(5,1)+svs(iq,n)*0.2*vlin(iq,n)
          newlai(5,2)=newlai(5,2)+svs(iq,n)*0.2*vlinnext(iq,n)
          newlai(5,3)=newlai(5,3)+svs(iq,n)*0.2*vlinnext2(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.8*fg3
          newlai(6,0)=newlai(6,0)+svs(iq,n)*0.8*fg3*vlinprev(iq,n)
          newlai(6,1)=newlai(6,1)+svs(iq,n)*0.8*fg3*vlin(iq,n)
          newlai(6,2)=newlai(6,2)+svs(iq,n)*0.8*fg3*vlinnext(iq,n)
          newlai(6,3)=newlai(6,3)+svs(iq,n)*0.8*fg3*vlinnext2(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.8*fg4
          newlai(7,0)=newlai(7,0)+svs(iq,n)*0.8*fg4*vlinprev(iq,n)
          newlai(7,1)=newlai(7,1)+svs(iq,n)*0.8*fg4*vlin(iq,n)
          newlai(7,2)=newlai(7,2)+svs(iq,n)*0.8*fg4*vlinnext(iq,n)
          newlai(7,3)=newlai(7,3)+svs(iq,n)*0.8*fg4*vlinnext2(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.8*ftu
          newlai(8,0)=newlai(8,0)+svs(iq,n)*0.8*ftu*vlinprev(iq,n)
          newlai(8,1)=newlai(8,1)+svs(iq,n)*0.8*ftu*vlin(iq,n)
          newlai(8,2)=newlai(8,2)+svs(iq,n)*0.8*ftu*vlinnext(iq,n)
          newlai(8,3)=newlai(8,3)+svs(iq,n)*0.8*ftu*vlinnext2(iq,n)
        case (8)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4
            newlai(1,0)=newlai(1,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(1,1)=newlai(1,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(1,2)=newlai(1,2)+svs(iq,n)*0.4*vlinnext(iq,n)
            newlai(1,3)=newlai(1,3)+svs(iq,n)*0.4*vlinnext2(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4*xp
            newlai(1,0)=newlai(1,0)+svs(iq,n)*vlinprev(iq,n)*0.4*xp
            newlai(1,1)=newlai(1,1)+svs(iq,n)*vlin(iq,n)*0.4*xp
            newlai(1,2)=newlai(1,2)+svs(iq,n)*vlinnext(iq,n)*0.4*xp   
            newlai(1,3)=newlai(1,3)+svs(iq,n)*vlinnext2(iq,n)*0.4*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4*(1.-xp)
            newlai(18,0)=newlai(18,0)+svs(iq,n)*vlinprev(iq,n)*0.4*(1.-xp)
            newlai(18,1)=newlai(18,1)+svs(iq,n)*vlin(iq,n)*0.4*(1.-xp)
            newlai(18,2)=newlai(18,2)+svs(iq,n)*vlinnext(iq,n)*0.4*(1.-xp)
            newlai(18,3)=newlai(18,3)+svs(iq,n)*vlinnext2(iq,n)*0.4*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4
            newlai(18,0)=newlai(18,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(18,1)=newlai(18,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(18,2)=newlai(18,2)+svs(iq,n)*0.4*vlinnext(iq,n)
            newlai(18,3)=newlai(18,3)+svs(iq,n)*0.4*vlinnext2(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.6*fg3
          newlai(6,0)=newlai(6,0)+svs(iq,n)*0.6*fg3*vlinprev(iq,n)
          newlai(6,1)=newlai(6,1)+svs(iq,n)*0.6*fg3*vlin(iq,n)
          newlai(6,2)=newlai(6,2)+svs(iq,n)*0.6*fg3*vlinnext(iq,n)
          newlai(6,3)=newlai(6,3)+svs(iq,n)*0.6*fg3*vlinnext2(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.6*fg4
          newlai(7,0)=newlai(7,0)+svs(iq,n)*0.6*fg4*vlinprev(iq,n)
          newlai(7,1)=newlai(7,1)+svs(iq,n)*0.6*fg4*vlin(iq,n)
          newlai(7,2)=newlai(7,2)+svs(iq,n)*0.6*fg4*vlinnext(iq,n)
          newlai(7,3)=newlai(7,3)+svs(iq,n)*0.6*fg4*vlinnext2(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.6*ftu
          newlai(8,0)=newlai(8,0)+svs(iq,n)*0.6*ftu*vlinprev(iq,n)
          newlai(8,1)=newlai(8,1)+svs(iq,n)*0.6*ftu*vlin(iq,n)
          newlai(8,2)=newlai(8,2)+svs(iq,n)*0.6*ftu*vlinnext(iq,n)
          newlai(8,3)=newlai(8,3)+svs(iq,n)*0.6*ftu*vlinnext2(iq,n)
        case (9)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1
            newlai(1,0)=newlai(1,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(1,1)=newlai(1,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(1,2)=newlai(1,2)+svs(iq,n)*0.1*vlinnext(iq,n)
            newlai(1,3)=newlai(1,3)+svs(iq,n)*0.1*vlinnext2(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1*xp
            newlai(1,0)=newlai(1,0)+svs(iq,n)*vlinprev(iq,n)*0.1*xp
            newlai(1,1)=newlai(1,1)+svs(iq,n)*vlin(iq,n)*0.1*xp
            newlai(1,2)=newlai(1,2)+svs(iq,n)*vlinnext(iq,n)*0.1*xp
            newlai(1,3)=newlai(1,3)+svs(iq,n)*vlinnext2(iq,n)*0.1*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1*(1.-xp)
            newlai(18,0)=newlai(18,0)+svs(iq,n)*vlinprev(iq,n)*0.1*(1.-xp)
            newlai(18,1)=newlai(18,1)+svs(iq,n)*vlin(iq,n)*0.1*(1.-xp)
            newlai(18,2)=newlai(18,2)+svs(iq,n)*vlinnext(iq,n)*0.1*(1.-xp)
            newlai(18,3)=newlai(18,3)+svs(iq,n)*vlinnext2(iq,n)*0.1*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1
            newlai(18,0)=newlai(18,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(18,1)=newlai(18,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(18,2)=newlai(18,2)+svs(iq,n)*0.1*vlinnext(iq,n)
            newlai(18,3)=newlai(18,3)+svs(iq,n)*0.1*vlinnext2(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.9*fg3
          newlai(6,0)=newlai(6,0)+svs(iq,n)*0.9*fg3*vlinprev(iq,n)
          newlai(6,1)=newlai(6,1)+svs(iq,n)*0.9*fg3*vlin(iq,n)
          newlai(6,2)=newlai(6,2)+svs(iq,n)*0.9*fg3*vlinnext(iq,n)
          newlai(6,3)=newlai(6,3)+svs(iq,n)*0.9*fg3*vlinnext2(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.9*fg4
          newlai(7,0)=newlai(7,0)+svs(iq,n)*0.9*fg4*vlinprev(iq,n)
          newlai(7,1)=newlai(7,1)+svs(iq,n)*0.9*fg4*vlin(iq,n)
          newlai(7,2)=newlai(7,2)+svs(iq,n)*0.9*fg4*vlinnext(iq,n)
          newlai(7,3)=newlai(7,3)+svs(iq,n)*0.9*fg4*vlinnext2(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.9*ftu
          newlai(8,0)=newlai(8,0)+svs(iq,n)*0.9*ftu*vlinprev(iq,n)
          newlai(8,1)=newlai(8,1)+svs(iq,n)*0.9*ftu*vlin(iq,n)
          newlai(8,2)=newlai(8,2)+svs(iq,n)*0.9*ftu*vlinnext(iq,n)
          newlai(8,3)=newlai(8,3)+svs(iq,n)*0.9*ftu*vlinnext2(iq,n)
        case (10)
          newgrid(6)=newgrid(6)+svs(iq,n)*fg3
          newlai(6,0)=newlai(6,0)+svs(iq,n)*fg3*vlinprev(iq,n)
          newlai(6,1)=newlai(6,1)+svs(iq,n)*fg3*vlin(iq,n)
          newlai(6,2)=newlai(6,2)+svs(iq,n)*fg3*vlinnext(iq,n)
          newlai(6,3)=newlai(6,3)+svs(iq,n)*fg3*vlinnext2(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*fg4
          newlai(7,0)=newlai(7,0)+svs(iq,n)*fg4*vlinprev(iq,n)
          newlai(7,1)=newlai(7,1)+svs(iq,n)*fg4*vlin(iq,n)
          newlai(7,2)=newlai(7,2)+svs(iq,n)*fg4*vlinnext(iq,n)
          newlai(7,3)=newlai(7,3)+svs(iq,n)*fg4*vlinnext2(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*ftu
          newlai(8,0)=newlai(8,0)+svs(iq,n)*ftu*vlinprev(iq,n)
          newlai(8,1)=newlai(8,1)+svs(iq,n)*ftu*vlin(iq,n)
          newlai(8,2)=newlai(8,2)+svs(iq,n)*ftu*vlinnext(iq,n)
          newlai(8,3)=newlai(8,3)+svs(iq,n)*ftu*vlinnext2(iq,n)
        case (12,14)
          newgrid(9)=newgrid(9)+svs(iq,n)*fc3
          newlai(9,0)=newlai(9,0)+svs(iq,n)*fc3*vlinprev(iq,n)
          newlai(9,1)=newlai(9,1)+svs(iq,n)*fc3*vlin(iq,n)
          newlai(9,2)=newlai(9,2)+svs(iq,n)*fc3*vlinnext(iq,n)
          newlai(9,3)=newlai(9,3)+svs(iq,n)*fc3*vlinnext2(iq,n)
          newgrid(10)=newgrid(10)+svs(iq,n)*fc4
          newlai(10,0)=newlai(10,0)+svs(iq,n)*fc4*vlinprev(iq,n)
          newlai(10,1)=newlai(10,1)+svs(iq,n)*fc4*vlin(iq,n)
          newlai(10,2)=newlai(10,2)+svs(iq,n)*fc4*vlinnext(iq,n)
          newlai(10,3)=newlai(10,3)+svs(iq,n)*fc4*vlinnext2(iq,n)
        case (13)
          newgrid(15)=newgrid(15)+svs(iq,n)
          newlai(15,0)=newlai(15,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(15,1)=newlai(15,1)+svs(iq,n)*vlin(iq,n)
          newlai(15,2)=newlai(15,2)+svs(iq,n)*vlinnext(iq,n)
          newlai(15,3)=newlai(15,3)+svs(iq,n)*vlinnext2(iq,n)
        case (15)
          newgrid(17)=newgrid(17)+svs(iq,n)
          newlai(17,0)=newlai(17,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(17,1)=newlai(17,1)+svs(iq,n)*vlin(iq,n)
          newlai(17,2)=newlai(17,2)+svs(iq,n)*vlinnext(iq,n)
          newlai(17,3)=newlai(17,3)+svs(iq,n)*vlinnext2(iq,n)
        case (16)
          newgrid(14)=newgrid(14)+svs(iq,n)
          newlai(14,0)=newlai(14,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(14,1)=newlai(14,1)+svs(iq,n)*vlin(iq,n)
          newlai(14,2)=newlai(14,2)+svs(iq,n)*vlinnext(iq,n)
          newlai(14,3)=newlai(14,3)+svs(iq,n)*vlinnext2(iq,n)
        case (17)
          newgrid(16)=newgrid(16)+svs(iq,n)
          newlai(16,0)=newlai(16,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(16,1)=newlai(16,1)+svs(iq,n)*vlin(iq,n)
          newlai(16,2)=newlai(16,2)+svs(iq,n)*vlinnext(iq,n)
          newlai(16,3)=newlai(16,3)+svs(iq,n)*vlinnext2(iq,n)
        case DEFAULT
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq,ivs,land=",myid,iq,ivs(iq,n),land(iq)
          call ccmpi_abort(-1)
      end select
    end do
    where ( newgrid(:)>0. )
      newlai(:,0) = newlai(:,0)/newgrid(:)
      newlai(:,1) = newlai(:,1)/newgrid(:)
      newlai(:,2) = newlai(:,2)/newgrid(:)
      newlai(:,3) = newlai(:,3)/newgrid(:)
    end where
    ipos = count(newgrid(:)>0.)
    do while ( ipos>maxtile )
      pos = minloc(newgrid(:), newgrid(:)>0.)
      newgrid(pos(1)) = 0.
      nsum = sum(newgrid(:))
      newgrid(:) = newgrid(:)/nsum
      ipos = count(newgrid(:)>0.)
    end do    
    do while ( any(newgrid(:)<minfrac.and.newgrid(:)>0.) )
      pos = minloc(newgrid(:), newgrid(:)>0.)
      newgrid(pos(1)) = 0.
      nsum = sum(newgrid(:))
      newgrid(:) = newgrid(:)/nsum
    end do

    n = 0
    ivs(iq,:) = 0
    svs(iq,:) = 0.
    vlinprev(iq,:) = 0.
    vlin(iq,:)     = 0.
    vlinnext(iq,:) = 0.
    vlinnext2(iq,:) = 0.
    do iv = 1,18
      if ( newgrid(iv)>0. ) then
        n = n + 1
        ivs(iq,n)      = iv
        svs(iq,n)      = newgrid(iv)
        vlinprev(iq,n) = newlai(iv,0)
        vlin(iq,n)     = newlai(iv,1)
        vlinnext(iq,n) = newlai(iv,2)
        vlinnext2(iq,n) = newlai(iv,3)
      end if
    end do

  end if
end do

return
end subroutine convertigbp


subroutine cbmparm(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2, &
                   casapoint,greenup,fall,phendoy1)

use carbpools_m
use cc_mpi
use const_phys
use infile
use latlong_m
use newmpar_m
use nsibd_m
use parm_m
use pbl_m
use sigs_m
use soil_m
use soilsnow_m
use soilv_m
use vegpar_m
  
implicit none
  
integer, dimension(ifull,maxtile), intent(in) :: ivs
integer, dimension(271,mxvt), intent(in) :: greenup, fall, phendoy1
integer iq,n,k,ipos,iv,ilat,ivp
integer jyear,jmonth,jday,jhour,jmin,mins
integer landcount
integer, dimension(:), allocatable :: cveg
real fjd, ivmax
real, dimension(mxvt,mplant) :: ratiocnplant
real, dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real, dimension(ifull,maxtile), intent(in) :: svs,vlin,vlinprev,vlinnext,vlinnext2
real, dimension(ifull,5), intent(in) :: casapoint
real, dimension(ifull,2) :: albsoilsn
real, dimension(12,msoil) :: rationpsoil
real, dimension(mxvt) :: leafage,woodage,frootage,metage
real, dimension(mxvt) :: strage,cwdage,micage,slowage,passage
real, dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real, dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real, dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real, dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real, dimension(mxvt) :: cleaf,cwood,cfroot,cmet,cstr,ccwd,cmic,cslow,cpass,nleaf
real, dimension(mxvt) :: nwood,nfroot,nmet,nstr,ncwd,nmic,nslow,npass,xpleaf,xpwood
real, dimension(mxvt) :: xpfroot,xpmet,xpstr,xpcwd,xpmic,xpslow,xppass,clabileage
real, dimension(ifull) :: albsoil
real, dimension(12) :: xkmlabp,xpsorbmax,xfPleach

if ( myid==0 ) write(6,*) "Initialising CABLE"

if ( cbm_ms/=ms ) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  call ccmpi_abort(-1)
end if

! redefine rhos
rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

icycle=ccycle

if ( myid==0 ) write(6,*) "Define CABLE and CASA CNP arrays"

! default values (i.e., no land)  
ivegt=0
albsoilsn=0.08  
albsoil=0.08
albvisdir=0.08
albvisdif=0.08
albnirdir=0.08
albnirdif=0.08
zolnd=0.
!cplant=0.
!csoil=0.
pind=ifull+1
mvtype=mxvt
mstype=mxst

! calculate CABLE vector length
do iq = 1,ifull
  if ( land(iq) ) then
    landcount = count(svs(iq,:)>0.)
    mp = mp + landcount
    if ( landcount==0 ) then
      write(6,*) "ERROR: No CABLE tiles assigned to land point: myid,iq",myid,iq
      call ccmpi_abort(-1)
    end if
  end if
end do


!if (nmaxpr==1) then
!  !write(6,*) "myid,landtile ",myid,mp
!
!  lndtst(1)=0
!  if (mp>0) lndtst(1)=1
!  call ccmpi_reduce(lndtst(1:1),lndtst_g(1:1),"sum",0,comm_world)
!  if (myid==0) then
!    write(6,*) "Processors with land ",lndtst_g(1),nproc
!  end if
!end if

! if CABLE is present on this processor, then start allocating arrays
! Write messages here in case myid==0 has no land-points (mp==0)
if (myid==0) then
  write(6,*) "Allocating CABLE and CASA CNP arrays"
  if (icycle==0) then
    write(6,*) "Using CABLE without carbon cycle"
  else
    write(6,*) "Using CASA CNP"
  end if
end if

if (mp>0) then
  
  allocate(sv(mp))
  allocate(vl1(mp),vl2(mp),vl3(mp),vl4(mp))
  allocate(tmap(ifull,maxtile))
  call alloc_cbm_var(air, mp)
  call alloc_cbm_var(bgc, mp)
  call alloc_cbm_var(canopy, mp)
  call alloc_cbm_var(met, mp)
  call alloc_cbm_var(bal, mp)
  call alloc_cbm_var(rad, mp)
  call alloc_cbm_var(rough, mp)
  call alloc_cbm_var(soil, mp)
  call alloc_cbm_var(ssnow, mp)
  call alloc_cbm_var(sum_flux, mp)
  call alloc_cbm_var(veg, mp)
  if ( cable_climate==1 ) then
    call alloc_cbm_var(climate, mp)
  end if
  allocate(cveg(mp))

  ! Cable configuration
  cable_user%ssnow_POTEV = "P-M"
  cable_user%MetType = "default"
  cable_user%diag_soil_resp = "ON"
  cable_user%leaf_respiration = "ON"
  cable_user%run_diag_level = "NONE"
  cable_user%consistency_check = .false.
  cable_user%logworker = .false.
  cable_user%l_new_roughness_soil = .false.
  select case ( cable_climate )
    case(1)
      cable_user%call_climate = .true.
    case default
      cable_user%call_climate = .false.
  end select
  select case ( cable_pop )
    case(1)
      cable_user%call_pop = .true.
    case default
      cable_user%call_pop = .false.
  end select
  select case ( soil_struc )
    case(1)
      cable_user%soil_struc = "sli"  
    case default
      cable_user%soil_struc = "default"
  end select
  select case ( fwsoil_switch )
    case(3)
      cable_user%fwsoil_switch = "Haverd2003"
    case(2)
      cable_user%fwsoil_switch = "Lai and Ktaul 2000"  
    case(1)
      cable_user%fwsoil_switch = "non-linear extrapolation"    
    case default
      cable_user%fwsoil_switch = "standard"      
  end select
  select case ( gs_switch )  
    case(1)
      cable_user%gs_switch = "medlyn"  
    case default
      cable_user%gs_switch = "leuning"
  end select
  select case ( cable_litter )
    case(1)
      cable_user%litter = .true.  
    case default
      cable_user%litter = .false.
  end select
  select case ( phenology_switch )
    case(1)
      cable_user%phenology_switch = "climate"
      if ( cable_climate/=1 ) then
        write(6,*) "ERROR: phenology_switch=1 requires cable_climate=1"
        stop
      end if
    case default
      cable_user%phenology_switch = "MODIS"
  end select
  select case ( progvcmax )
    case(2)
      cable_user%vcmax = "Walker2014"        
    case default
      cable_user%vcmax = "standard"    
  end select
  kwidth_gl = nint(dt) ! MJT notes - what happens when the timestep is less than a second?
  if ( kwidth_gl==0) then
    write(6,*) "ERROR: Timestep too small for CABLE"
    call ccmpi_abort(-1)
  end if
  
  ! soil parameters
  soil%zse        = zse ! soil layer thickness
  soil%zshh(1)    = 0.5 * soil%zse(1)
  soil%zshh(ms+1) = 0.5 * soil%zse(ms)
  soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
 
  sv=0.
  vl1=0.
  vl2=0.
  vl3=0.
  vl4=0.
  tmap=.false.
  cveg=0

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  ! now up to maxtile=5 PFT tiles from 5 IGBP classes (need correct order for vectorisation)
  ipos = 0
  do n = 1,maxtile
    pind(n,1) = ipos + 1
    do iq = 1,ifull
      if ( land(iq) ) then
        if ( svs(iq,n)>0. ) then
          ipos = ipos + 1
          tmap(iq,n)        = .true.
          iv                = ivs(iq,n)
          sv(ipos)          = svs(iq,n)
          cveg(ipos)        = iv
          soil%isoilm(ipos) = isoilm(iq)
          vl1(ipos) = max( vlinprev(iq,n), 0.01 )
          vl2(ipos) = max( vlin(iq,n), 0.01 )
          vl3(ipos) = max( vlinnext(iq,n), 0.01 )
          vl4(ipos) = max( vlinnext2(iq,n), 0.01 )
        end if
      end if
    end do
    pind(n,2) = ipos
  end do
  
  if ( ipos/=mp ) then
    write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
    call ccmpi_abort(-1)
  end if

  ! Load CABLE biophysical arrays
  do iq = 1,ifull
    if ( land(iq) ) then
      ivmax = -1.
      do n = 1,maxtile
        if ( svs(iq,n)>ivmax ) then
          ivmax = svs(iq,n)
          ivegt(iq) = ivs(iq,n) ! diagnostic (CSIRO pft)
        end if
      end do
    end if
  end do
  
  call cable_biophysic_parm(cveg)
  
  where ( veg%iveg>=14 .and. veg%iveg<=17 )
    vl1(:) = 1.e-8  
    vl2(:) = 1.e-8
    vl3(:) = 1.e-8
    vl4(:) = 1.e-8
  end where
  
  deallocate( cveg )

  ! calculate actual max tile number
  do n = 1,maxtile
    if ( pind(n,1)<=mp ) then
      maxnb = n
    end if
  end do
  
  ! Load CABLE soil data
  soil%bch     = bch(soil%isoilm)
  soil%css     = css(soil%isoilm)
  soil%rhosoil = rhos(soil%isoilm)
  soil%cnsd    = cnsd(soil%isoilm)
  soil%hyds    = hyds(soil%isoilm)
  soil%sucs    = sucs(soil%isoilm)
  soil%hsbh    = hsbh(soil%isoilm)
  soil%sfc     = sfc(soil%isoilm)
  soil%ssat    = ssat(soil%isoilm)
  soil%swilt   = swilt(soil%isoilm)
  soil%ibp2    = ibp2(soil%isoilm)
  soil%i2bp3   = i2bp3(soil%isoilm)
  soil%pwb_min = (soil%swilt/soil%ssat)**soil%ibp2
  soil%clay    = clay(soil%isoilm)
  soil%sand    = sand(soil%isoilm)
  soil%silt    = silt(soil%isoilm)
  
  ! depeciated
  !bgc%ratecp(:) = ratecp(:)
  !bgc%ratecs(:) = ratecs(:)

  ! store bare soil albedo and define snow free albedo
  do n = 1,maxnb
    soil%albsoil(pind(n,1):pind(n,2),1) = pack(albvisnir(:,1),tmap(:,n))
    soil%albsoil(pind(n,1):pind(n,2),2) = pack(albvisnir(:,2),tmap(:,n))
  end do
  soil%albsoil(:,3) = 0.05
    
  where ( land(1:ifull) )
    albsoil(1:ifull) = 0.5*sum(albvisnir(1:ifull,:),2)
  end where
  where ( albsoil(1:ifull)<=0.14 .and. land(1:ifull) )
    !sfact=0.5 for alb <= 0.14
    albsoilsn(1:ifull,1) = (1.00/1.50)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.50)*albsoil(1:ifull)
  elsewhere ( albsoil(1:ifull)<=0.2 .and. land(1:ifull) )
    !sfact=0.62 for 0.14 < alb <= 0.20
    albsoilsn(1:ifull,1) = (1.24/1.62)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.62)*albsoil(1:ifull)
  elsewhere ( land(1:ifull) )
    !sfact=0.68 for 0.2 < alb
    albsoilsn(1:ifull,1) = (1.36/1.68)*albsoil(1:ifull)
    albsoilsn(1:ifull,2) = (2.00/1.68)*albsoil(1:ifull)
  end where

  do n = 1,maxnb
    ! MJT patch
    soil%albsoil(pind(n,1):pind(n,2),1)   =pack(albsoil,       tmap(:,n))
    soil%albsoil(pind(n,1):pind(n,2),2)   =pack(albsoil,       tmap(:,n))
    ssnow%albsoilsn(pind(n,1):pind(n,2),1)=pack(albsoilsn(:,1),tmap(:,n)) ! overwritten by CABLE
    ssnow%albsoilsn(pind(n,1):pind(n,2),2)=pack(albsoilsn(:,2),tmap(:,n)) ! overwritten by CABLE
    rad%albedo_T(pind(n,1):pind(n,2))     =pack(albsoil,       tmap(:,n))
    rad%trad(pind(n,1):pind(n,2))         =pack(tss,           tmap(:,n))
    rad%latitude(pind(n,1):pind(n,2))     =pack(rlatt,         tmap(:,n))*180./pi
    rad%longitude(pind(n,1):pind(n,2))    =pack(rlongg,        tmap(:,n))*180./pi
  end do
    
  ssnow%albsoilsn(:,3)=0.05    
  ssnow%t_snwlr=0.05
  ssnow%pudsmx=0.

  canopy%oldcansto=0.  
  canopy%ghflux=0.
  canopy%sghflux=0.
  canopy%ga=0.
  canopy%dgdtg=0.
  canopy%fhs_cor=0.
  canopy%fes_cor=0.
  canopy%us=0.01
  ssnow%wb_lake=0. ! not used when mlo.f90 is active
  ssnow%fland=1.
  ssnow%ifland=soil%isoilm
    
  ! Initialise sum flux variables
  sum_flux%sumpn=0.
  sum_flux%sumrp=0.
  sum_flux%sumrs=0.
  sum_flux%sumrd=0.
  sum_flux%sumrpw=0.
  sum_flux%sumrpr=0.
  sum_flux%dsumpn=0.
  sum_flux%dsumrp=0.
  sum_flux%dsumrs=0.
  sum_flux%dsumrd=0.
  
  bal%evap_tot=0.
  bal%precip_tot=0.
  bal%ebal_tot=0.
  bal%rnoff_tot=0.
  
  if ( cable_climate==1 ) then
    climate%doy = 0
    climate%nyears = 0    ! not used?
    climate%biome = 0     ! not used?
    climate%iveg = 0      ! not used?
    climate%chilldays = 0 ! used by phenology
    climate%gdd5 = 0.     ! used by phenology
    climate%gdd0 = 0.
    climate%agdd5 = 0.
    climate%agdd0 = 0.
    climate%dmoist = 0.
    climate%dtemp = 0.
    climate%mtemp = 0.
    climate%mmoist = 0.
    climate%qtemp = 0.
    climate%evap_PT = 0.      ! not used?
    climate%aevap = 0.        ! not used?
    climate%mtemp_max = 0.
    climate%mtemp_min = 0.
    climate%qtemp_max = 0.
    climate%qtemp_max_last_year = 0. ! used by cable_canopy and CASA-CNP
    climate%dtemp_91 = 0.
    climate%dtemp_31 = 0.
    climate%dmoist_31 = 0.
    climate%atemp_mean = 0.   ! not used?
    climate%mtemp_min_20 = 0. ! not used?
    climate%mtemp_min20 = 0.  ! not used?
    climate%mtemp_max_20 = 0. ! not used?
    climate%mtemp_max20 = 0.  ! not used?
    climate%alpha_PT_20 = 0.  ! not used?
    climate%alpha_PT20 = 0.   ! not used?
    climate%alpha_PT = 0.     ! not used?
  end if
  
  if ( icycle==0 ) then
    ! Initialise CABLE carbon pools
    !cplant=0.
    !csoil=0.
    !do n=1,maxnb
    !  do k=1,ncp
    !    cplant(:,k)=cplant(:,k)+unpack(sv(pind(n,1):pind(n,2))*tcplant(veg%iveg(pind(n,1):pind(n,2)),k),tmap(:,n),0.)
    !  end do
    !  do k=1,ncs
    !    csoil(:,k) =csoil(:,k) +unpack(sv(pind(n,1):pind(n,2))*tcsoil(veg%iveg(pind(n,1):pind(n,2)),k), tmap(:,n),0.)
    !  end do
    !end do
  else
    ! CASA CNP
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    call alloc_phenvariable(phen,mp)
    
    casamet%lat=rad%latitude
    
    do n = 1,maxnb
      casamet%isorder(pind(n,1):pind(n,2))  = nint(pack(casapoint(:,1),tmap(:,n)))
      casaflux%Nmindep(pind(n,1):pind(n,2)) = pack(casapoint(:,2),tmap(:,n))/365.*1.E-3
      casaflux%Nminfix(pind(n,1):pind(n,2)) = pack(casapoint(:,3),tmap(:,n))/365.
      casaflux%Pdep(pind(n,1):pind(n,2))    = pack(casapoint(:,4),tmap(:,n))/365.
      casaflux%Pwea(pind(n,1):pind(n,2))    = pack(casapoint(:,5),tmap(:,n))/365.
    end do

    where ( veg%iveg==9 .or. veg%iveg==10 ) ! crops
      ! P fertilizer =13 Mt P globally in 1994
      casaflux%Pdep = casaflux%Pdep + 0.7/365.
      ! N fertilizer =86 Mt N globally in 1994
      casaflux%Nminfix = casaflux%Nminfix + 4.3/365.
    end where

    if ( any( casamet%isorder<1 .or. casamet%isorder>12 ) ) then
      write(6,*) "ERROR: Invalid isorder in CASA-CNP"
      call ccmpi_abort(-1)
    end if

    leafage =(/ 2.0,1.5,1.0,1.0,1.0,0.8,0.8,1.0,     0.8,     0.8,1.0,1.0,1.0,1.0,1.0,1.0,1.0 /)
    woodage =(/ 70.,60.,80.,40.,40.,1.0,1.0,1.0,     1.0,     1.0,1.0,1.0,1.0,5.0,1.0,1.0,1.0 /)
    frootage=(/ 18.,10.,10.,10.,5.0,3.0,3.0,3.0,0.884227,0.884227,1.0,1.0,1.0,4.0,1.0,1.0,1.0 /)
    metage=0.04
    strage=0.23
    cwdage=0.824
    micage=0.137
    slowage=5.
    passage=222.22
    clabileage=0.2

    xfherbivore   =(/ 0.068,0.406,0.068,0.134,0.022,0.109,0.109,0.109,0.140,0.140,0.000,0.000,0.000,0.010,0.000,0.000,0.000 /)
    xxkleafcoldmax=(/   0.2,  0.1,  0.1,  0.6,   1.,  0.2,  0.2,  0.2,  0.3,  0.3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1 /)
    xxkleafdrymax =(/   0.1,  0.1,  0.1,   1.,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1 /)
    xratioNPleafmin =(/ 10.92308,15.95339,9.254839,12.73848,12.07217,13.51473,   14.05,12.57800,15.12262,10.,13.,10.,10., 16.2336, &
                        10.,10.,10. /)
    xratioNPleafmax =(/ 12.07288, 17.6327,10.22903,14.07938,13.34292,14.93733,15.52895,  13.902,16.71447,10.,13.,10.,10., 17.9424, &
                        10.,10.,10. /)
    xratioNPwoodmin =(/ 20.30167,15.89425,17.48344,19.08018,22.46035,     15.,     15.,   15.96,   20.52,15.,15.,15.,15., 17.5275, &
                        15.,15.,15. /)
    xratioNPwoodmax =(/ 22.43869,17.56733, 19.3238,21.08862, 24.8246,     15.,     15.,   17.64,   20.52,15.,15.,15.,15., 19.3725, &
                        15.,15.,15. /)
    xratioNPfrootmin=(/ 20.29341,15.87155,17.39767, 19.0601,22.49363,15.63498,16.08255,14.49241,22.69109,15.,15.,15.,15.,22.13268, &
                        15.,15.,15. /)
    xratioNPfrootmax=(/ 22.42955,17.54224,  19.229,21.06643,24.86138,17.28077,17.77545,16.01793,25.07962,15.,15.,15.,15.,24.46244, &
                        15.,15.,15. /)
    xfNminloss=0.05
    xfNminleach=0.05
    xnfixrate=(/ 0.08,2.6,0.21,1.64,0.37,0.95,0.95,0.95,4.,4.,0.,0.,0.,0.35,0.,0.,0. /)
    xnsoilmin=1000.
    
    ratiocnplant(:,1)=(/  49.8, 23.1, 59.3, 31.4, 37.6, 34.8,  44., 49.2, 21.6, 25., 30., 30., 30., 50., 40., 40., 40. /)
    ratiocnplant(:,2)=(/ 238.1,134.9,243.8,156.2,142.1, 150., 150.,147.3, 150.,125.,150.,150.,150.,150.,150.,135.,150. /)
    ratiocnplant(:,3)=(/  73.7, 61.2,  75., 63.2, 67.1, 64.5, 62.7,  69., 60.7, 71., 71., 71., 71., 71., 71., 71., 71. /)
    ratiocnsoil(:,1)=8.
    ratiocnsoil(:,2)=(/ 16.1,12.8,24.8, 30.,19.3,13.1,13.1,13.1,13.2,13.2,13.1,13.1,13.1,26.8, 20., 20., 20. /)
    ratiocnsoil(:,3)=(/ 16.1,12.8,24.8, 30.,19.3,13.1,13.1,13.1,13.2,13.2,13.1,13.1,13.1,26.8, 20., 20., 20. /)
    ratiocnsoilmin(:,1)=3.
    ratiocnsoilmin(:,2)=12.
    ratiocnsoilmin(:,3)=7.
    ratiocnsoilmax(:,1)=15.
    ratiocnsoilmax(:,2)=30.
    ratiocnsoilmax(:,3)=15.
     
    ! Initial values for CNP pools over 3*plant, 3*litter and 3*soil (=27 pools in total)
    cleaf  =(/ 384.6037,    273.,96.59814,150.2638,     88.,137.1714,137.1714,137.1714,    160.,    160.,0.,0.,0.,      0.,0.,0., &
              0. /)
    cwood  =(/ 7865.396,  11451.,5683.402,10833.74,    372.,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,      0.,0.,0., &
              0. /)
    cfroot =(/     250.,   2586.,    220.,    220.,    140.,    263.,    263.,    263.,    240.,    240.,0.,0.,0.,      0.,0.,0., &
              0. /)
    cmet   =(/ 6.577021,44.63457,7.127119,10.97797,3.229374,28.57245,28.57245,28.57245,28.57245,28.57245,0.,0.,0.,1.457746,0.,0., &
              0. /)
    cstr   =(/ 209.1728,433.7626,277.7733,312.5492,39.44449,50.91091,50.91091,50.91091,50.91091,50.91091,0.,0.,0.,4.956338,0.,0., &
              0. /)
    ccwd   =(/ 606.0255,1150.765,776.7331,888.5864,111.5864,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,28.44085,0.,0., &
              0. /) 
    cmic   =(/  528.664,11.37765,597.0785,405.5554,168.0451,425.6431,425.6431,425.6431,512.4247,512.4247,0.,0.,0.,57.77585,0.,0., &
              0. /)
    cslow  =(/ 13795.94,311.8092,16121.12,11153.25,4465.478,5694.437,5694.437,5694.437,6855.438,6855.438,0.,0.,0.,1325.052,0.,0., &
              0. /)
    cpass  =(/ 4425.396,13201.81,5081.802,5041.192,1386.477, 4179.92, 4179.92, 4179.92,5032.137,5032.137,0.,0.,0.,517.1719,0.,0., &
              0. /)
    nleaf  =(/ 7.541249,     9.9,1.609969,3.756594,2.933333,4.572381,4.572381,4.572381,5.333333,5.333333,0.,0.,0.,     0.5,0.,0., &
              0. /)
    nwood  =(/ 31.46159,    102.,22.73361,80.24989,2.755555,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,0.125926,0.,0., &
              0. /)
    nfroot =(/ 6.097561,     38.,5.365854,5.365854,3.414634,6.414634,6.414634,6.414634,5.853659,5.853659,0.,0.,0.,1.536585,0.,0., &
              0. /)
    nmet   =(/ 0.064481, 0.74391,0.059393,0.137225,0.053823,0.476208,0.476208,0.476208,0.476208,0.476208,0.,0.,0.,0.018222,0.,0., &
              0. /)
    nstr   =(/ 1.394485,2.891751,1.851822,2.083661,0.262963,0.339406,0.339406,0.339406,0.339406,0.339406,0.,0.,0.,0.033042,0.,0., &
              0. /)
    ncwd   =(/ 2.424102,8.524183,3.106932,6.581996,0.826566,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,0.210673,0.,0., &
              0. /)
    nmic   =(/  52.8664,1.137765,59.70785,40.55554,16.80451,42.56431,42.56431,42.56431,51.24247,51.24247,0.,0.,0.,5.777585,0.,0., &
               0. /)
    nslow  =(/ 919.7293,20.78728,1074.741,743.5501,297.6985,379.6291,379.6291,379.6291,457.0292,457.0292,0.,0.,0.,88.33682,0.,0., &
               0. /)
    npass  =(/ 295.0264,880.1209,338.7868,336.0795, 92.4318,278.6613,278.6613,278.6613,335.4758,335.4758,0.,0.,0.,34.47813,0.,0., &
               0. /)
    xpleaf =(/ 0.191648,   0.415,0.115988,0.135453,0.022821, 0.15125, 0.15125, 0.15125, 0.15125, 0.15125,0.,0.,0.,   0.007,0.,0., &
               0. /)
    xpwood =(/ 0.953979,    5.88, 0.64438,2.424778,      0.,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,      0.,0.,0., &
               0. /)
    xpfroot=(/ 0.076659,    1.95,0.080548,0.141097,0.037083, 0.15125, 0.15125, 0.15125, 0.15125, 0.15125,0.,0.,0., 0.00875,0.,0., &
               0. /)
    xpmet  =(/ 0.004385,0.029756,0.004751,0.007319,0.002153,0.019048,0.019048,0.019048,0.019048,0.019048,0.,0.,0.,0.000972,0.,0., &
               0. /)
    xpstr  =(/ 0.069724,0.144588,0.092591,0.104183,0.013148, 0.01697, 0.01697, 0.01697, 0.01697, 0.01697,0.,0.,0.,0.001652,0.,0., &
               0. /)
    xpcwd  =(/ 0.101004,0.191794,0.129456,0.148095,0.018598,      0.,      0.,      0.,      0.,      0.,0.,0.,0.,      0.,0.,0., &
               0. /)
    xpmic  =(/ 6.872632, 0.14791,7.762021, 5.27222,2.184586,5.533361,5.533361,5.533361,6.661522,6.661522,0.,0.,0.,0.751086,0.,0., &
               0. /)
    xpslow =(/ 119.5648,2.702347,139.7164,96.66152,38.70081,49.35178,49.35178,49.35178, 59.4138, 59.4138,0.,0.,0.,11.48379,0.,0., &
               0. /)
    xppass =(/ 38.35343,114.4157,44.04228,43.69033,12.01613,36.22598,36.22598,36.22598,43.61185,43.61185,0.,0.,0.,4.482157,0.,0., &
               0. /)
    xplab  =(/   26.737,  19.947,  29.107,  30.509,  23.206,  25.538,  25.538,  25.538,  27.729,  27.729,0.,0.,0.,  21.038,0.,0., &
                 0.103 /)
    xpsorb =(/   126.73,  92.263, 134.639, 132.012,  173.47, 186.207, 186.207, 186.207, 155.518, 155.518,0.,0.,0.,  255.79,0.,0., &
                 1.176 /)
    xpocc  =(/  138.571, 120.374,  138.22, 148.083, 114.496, 145.163, 145.163, 145.163, 158.884, 158.884,0.,0.,0., 108.897,0.,0., &
                0.688 /)
 
    xkmlabp  =(/ 74.5408, 68.1584,  77.952,64.41918,64.41918,70.5856, 64.5888,54.1692, 9.7704, 28.29,  63.963,  32.402 /)
    xpsorbmax=(/ 745.408,788.0815,1110.816, 744.847, 744.847,816.146,746.8081,722.256,293.112,311.19,373.1175,615.6381 /)
    xfPleach =0.0005
    rationpsoil(:,1)=4.
    rationpsoil(:,2)=(/ 5.,5.,5.,15.,5.,5.,5.,5.,7.,7.,7.,7. /)
    rationpsoil(:,3)=(/ 5.,5.,5.,15.,5.,5.,5.,5.,7.,7.,7.,7. /)
 
    casabiome%ivt2     =(/   3,  3,  3,  3,  2,  1,  1,  2,  1,  1,  0,  0,  0,  1,  0,  0,  0 /)
    casabiome%kroot    =(/ 5.5,3.9,5.5,3.9,2.0,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,2.0,2.0,5.5,5.5 /)
    casabiome%rootdepth=(/ 1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.5,0.5 /)
    casabiome%kuptake  =(/ 2.0,1.9,2.0,2.0,1.8,2.0,2.0,2.0,1.6,1.6,1.6,1.8,1.8,1.8,1.8,1.8,1.8 /)
    casabiome%krootlen =(/ 14.87805,14.38596,14.02597,18.94737,32.30769,84.,84.,84.,120.5,120.5, &
                           0.,0.,0.,30.76923,0.,0.,0. /)
    casabiome%kminN=2
    casabiome%kuplabP=0.5
    casabiome%fracnpptoP(:,leaf) =(/ 0.25,0.20,0.40,0.35,0.35,0.35,0.35,0.50,0.50,0.50,0.50,0.50,0.50,0.25,0.50,0.60,0.50 /)
    casabiome%fracnpptoP(:,wood) =(/ 0.40,0.35,0.30,0.25,0.25,0.00,0.00,0.10,0.00,0.00,0.00,0.00,0.00,0.25,0.00,0.40,0.00 /)
    casabiome%fracnpptoP(:,xroot)=(/ 0.35,0.45,0.30,0.40,0.40,0.65,0.65,0.40,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.00,0.50 /) 
    casabiome%rmplant(:,leaf)    =0.1
    casabiome%rmplant(:,wood)    =(/ 2.0,1.0,1.5,0.8,0.5,0.5,0.4,1.8,2.0,1.0,1.0,1.0,1.0,2.0,1.0,1.0,1.0 /)
    casabiome%rmplant(:,xroot)   =(/ 10.,2.0,7.5,2.5,4.5,4.5,4.0,15.,25.,10.,10.,10.,10.,10.,10.,10.,10. /)
    casabiome%ftransNPtoL(:,leaf) =0.5
    casabiome%ftransNPtoL(:,wood) =0.95
    casabiome%ftransNPtoL(:,xroot)=0.9
    casabiome%fracligninplant(:,leaf) =(/ 0.25,0.20,0.20,0.20,0.20,0.10,0.10,0.10,0.10,0.10,0.15,0.15,0.15,0.15,0.15,0.25,0.10 /)
    casabiome%fracligninplant(:,wood) =0.4
    casabiome%fracligninplant(:,xroot)=(/ 0.25,0.20,0.20,0.20,0.20,0.10,0.10,0.10,0.10,0.10,0.15,0.15,0.15,0.15,0.15,0.25,0.10 /)
    casabiome%glaimax=(/ 7.,7.,7.,7.,3.,3.,3.,3.,6.,6., 5., 5., 5., 1.,6., 1.,0. /)
    casabiome%glaimin=(/ 1.,1.,.5,.5,.1,.1,.1,.1,.1,.1,.05,.05,.05,.05,0.,.05,0. /)
    phen%TKshed=(/ 268.,260.,263.15,268.15,277.15,275.15,275.15,275.15,278.15,278.15,277.15,277.15,277.15,277.15,277.15,277.15, &
                   283.15 /)
    casabiome%xkleafcoldexp=3.
    casabiome%xkleafdryexp=3.
    casabiome%ratioNCplantmin(:,leaf) =(/     0.02,    0.04,0.016667,0.028571,   0.025, 0.02631,    0.02,    0.02,    0.04, &
                                              0.04,0.033333,   0.025,   0.025,0.018182,   0.025,   0.025,   0.025 /)
    casabiome%ratioNCplantmax(:,leaf) =(/    0.024,   0.048,    0.02,0.034286,    0.03,0.031572,   0.024,   0.024,   0.048, &
                                             0.048,    0.04,    0.03,    0.03,0.022222,    0.03,    0.03,    0.03 /)
    casabiome%ratioNCplantmin(:,wood) =(/    0.004,0.006667,   0.004,0.005714,0.006667,0.006667,0.006667,0.006667,   0.008, &
                                             0.008,0.006667,0.006667,0.006667,0.006667,0.006667,0.007307,0.006667 /)
    casabiome%ratioNCplantmax(:,wood) =(/   0.0048,   0.008,  0.0048,0.006857,   0.008,   0.008,   0.008,   0.008,  0.0096, &
                                            0.0096,   0.008,   0.008,   0.008,   0.008,   0.008,0.008889,   0.008 /)
    casabiome%ratioNCplantmin(:,xroot)=(/ 0.012821,0.014706,0.012821,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085, &
                                          0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085 /)
    casabiome%ratioNCplantmax(:,xroot)=(/ 0.015385,0.017647,0.015385,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901, &
                                          0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901 /)
    casabiome%ftransPPtoL(:,leaf)=0.5
    casabiome%ftransPPtoL(:,wood)=0.95
    casabiome%ftransPPtoL(:,xroot)=0.9
    casabiome%ratioPcplantmin(:,leaf)  = 1./(xratioNPleafmax*ratioCNplant(:,leaf))
    casabiome%ratioPcplantmax(:,leaf)  = 1./(xratioNPleafmin*ratioCNplant(:,leaf))
    casabiome%ratioPcplantmin(:,wood)  = 1./(xratioNPwoodmax*ratioCNplant(:,wood))
    casabiome%ratioPcplantmax(:,wood)  = 1./(xratioNPwoodmin*ratioCNplant(:,wood))
    casabiome%ratioPcplantmin(:,xroot) = 1./(xratioNPfrootmax*ratioCNplant(:,xroot))
    casabiome%ratioPcplantmax(:,xroot) = 1./(xratioNPfrootmin*ratioCNplant(:,xroot))
    casabiome%sla                = 0.025/sqrt(leafage) ! see eqn A1 of Arora and Boer, GCB, 2005
    casabiome%fraclabile(:,leaf) = deltcasa*0.6    !1/day
    casabiome%fraclabile(:,xroot)= deltcasa*0.4    !1/day
    casabiome%fraclabile(:,wood) = 0.
    casabiome%plantrate(:,leaf)  = deltcasa/(leafage*(1.-xfherbivore))
    casabiome%plantrate(:,xroot) = deltcasa/frootage
    casabiome%plantrate(:,wood)  = deltcasa/woodage
    casabiome%litterrate(:,metb) = deltcasa/metage
    casabiome%litterrate(:,str)  = deltcasa/strage
    casabiome%litterrate(:,cwd)  = deltcasa/cwdage
    casabiome%soilrate(:,mic)    = deltcasa/micage
    casabiome%soilrate(:,slow)   = deltcasa/slowage
    casabiome%soilrate(:,pass)   = deltcasa/passage
    casabiome%xkleafcoldmax      = deltcasa*xxkleafcoldmax
    casabiome%xkleafdrymax       = deltcasa*xxkleafdrymax
    casabiome%rmplant            = casabiome%rmplant*deltcasa 
    casabiome%kclabrate          = deltcasa/clabileage

    casamet%iveg2 =casabiome%ivt2(veg%iveg)
    where (casamet%iveg2==3.or.casamet%iveg2==2)
      casamet%lnonwood = 0
      casapool%cplant(:,wood)  = cwood(veg%iveg) 
      casapool%clitter(:,cwd)  = ccwd(veg%iveg)
      casapool%nplant(:,wood)  = nwood(veg%iveg) 
      casapool%nlitter(:,cwd)  = ncwd(veg%iveg)
      casapool%pplant(:,wood)  = xpwood(veg%iveg)
      casapool%plitter(:,cwd)  = xpcwd(veg%iveg)
    elsewhere
      casamet%lnonwood = 1
      casapool%cplant(:,wood)  = 0.
      casapool%clitter(:,cwd)  = 0.
      casapool%nplant(:,wood)  = 0.
      casapool%nlitter(:,cwd)  = 0.
      casapool%pplant(:,wood)  = 0.
      casapool%plitter(:,cwd)  = 0.    
    end where
    casapool%cplant(:,leaf)     = cleaf(veg%iveg)
    casapool%cplant(:,xroot)    = cfroot(veg%iveg)
    casapool%clabile            = 0.
    casapool%clitter(:,metb)    = cmet(veg%iveg)
    casapool%clitter(:,str)     = cstr(veg%iveg)
    casapool%csoil(:,mic)       = cmic(veg%iveg)
    casapool%csoil(:,slow)      = cslow(veg%iveg)
    casapool%csoil(:,pass)      = cpass(veg%iveg)
    ! initializing glai in case not reading pool file (eg. during spin)
    casamet%glai = max(casabiome%glaimin(veg%iveg), casabiome%sla(veg%iveg)*casapool%cplant(:,leaf))
    casaflux%fNminloss   = xfNminloss(veg%iveg) 
    casaflux%fNminleach  = 10.*xfNminleach(veg%iveg)*deltcasa
    casapool%nplant(:,leaf) = nleaf(veg%iveg)
    casapool%nplant(:,xroot)= nfroot(veg%iveg)
    casapool%nlitter(:,metb)= nmet(veg%iveg)
    casapool%nlitter(:,str) = cstr(veg%iveg)*ratioNCstrfix
    casapool%nsoil(:,mic)   = nmic(veg%iveg)
    casapool%nsoil(:,slow)  = nslow(veg%iveg)
    casapool%nsoil(:,pass)  = npass(veg%iveg) 
    casapool%nsoilmin       = xnsoilmin(veg%iveg) 
    casapool%pplant(:,leaf) = xpleaf(veg%iveg)
    casapool%pplant(:,xroot)= xpfroot(veg%iveg) 
    casapool%plitter(:,metb)= xpmet(veg%iveg)
    casapool%plitter(:,str) = casapool%nlitter(:,str)/ratioNPstrfix
    casapool%psoil(:,mic)   = xpmic(veg%iveg)
    casapool%psoil(:,slow)  = xpslow(veg%iveg)
    casapool%psoil(:,pass)  = xppass(veg%iveg)
    casapool%psoillab       = xplab(veg%iveg)
    casapool%psoilsorb      = xpsorb(veg%iveg)
    casapool%psoilocc       = xpocc(veg%iveg)
    casaflux%kmlabp         = xkmlabp(casamet%isorder)
    casaflux%psorbmax       = xpsorbmax(casamet%isorder)
    casaflux%fpleach        = xfPleach(casamet%isorder)
    casapool%rationcplant   = 1./ratioCNplant(veg%iveg,:)
    casapool%ratiopcplant   = casabiome%ratioPcplantmax(veg%iveg,:)
    casapool%rationclitter  = casapool%nlitter/(casapool%clitter(:,:)+1.0e-10)
    casapool%ratiopclitter  = casapool%plitter/(casapool%clitter(:,:)+1.0e-10)
    casapool%ratioNCsoil    = 1./ratioCNsoil(veg%iveg,:)
    casapool%ratioPCsoil    = 1./(ratioCNsoil(veg%iveg,:)*ratioNPsoil(casamet%isorder,:))
    casapool%ratioNCsoilmin = 1./ratioCNsoilmax(veg%iveg,:)
    casapool%ratioNCsoilmax = 1./ratioCNsoilmin(veg%iveg,:)
    casapool%ratioNCsoilnew = casapool%ratioNCsoilmax
    casapool%Nsoil          = casapool%ratioNCsoil*casapool%Csoil
    casapool%Psoil          = casapool%ratioPCsoil*casapool%Csoil
    casapool%psoilsorb      = casaflux%psorbmax*casapool%psoillab &
                              /(casaflux%kmlabp+casapool%psoillab)
    
    do n = 1,mp
      ilat = nint((rad%latitude(n)+55.25)*2.) + 1
      ilat = min( 271, max( 1, ilat ) )
      ivp = veg%iveg(n)
      phen%aphen         = 0.
      phen%phase(n)      = phendoy1(ilat,ivp)
      phen%doyphase(n,1) = greenup(ilat,ivp)          ! DOY for greenup
      phen%doyphase(n,2) = phen%doyphase(n,1) + 14    ! DOY for steady LAI
      phen%doyphase(n,3) = fall(ilat,ivp)             ! DOY for leaf senescence
      phen%doyphase(n,4) = phen%doyphase(n,3) + 14    ! DOY for minimal LAI season
      if ( phen%doyphase(n,2) > 365 ) phen%doyphase(n,2) = phen%doyphase(n,2) - 365
      if ( phen%doyphase(n,4) > 365 ) phen%doyphase(n,4) = phen%doyphase(n,4) - 365
    end do
    
    casamet%tairk     = 0.
    casamet%tsoil     = 0.
    casamet%moist     = 0.
    casaflux%cgpp     = 0.
    casaflux%Crsoil   = 0.
    casaflux%crgplant = 0.
    casaflux%crmplant = 0.
    casaflux%clabloss = 0.

    cplant=0.
    clitter=0.
    csoil=0.
    niplant=0.
    nilitter=0.
    nisoil=0.
    pplant=0.
    plitter=0.
    psoil=0.
    !glai=0.
    do n=1,maxnb
      do k=1,mplant
        cplant(:,k) =cplant(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%cplant(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        niplant(:,k)=niplant(:,k)+unpack(sv(pind(n,1):pind(n,2))*real(casapool%nplant(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        pplant(:,k) =pplant(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%pplant(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
      end do
      do k=1,mlitter
        clitter(:,k) =clitter(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%clitter(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        nilitter(:,k)=nilitter(:,k)+unpack(sv(pind(n,1):pind(n,2))*real(casapool%nlitter(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        plitter(:,k) =plitter(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%plitter(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
      end do
      do k=1,msoil
        csoil(:,k) =csoil(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%csoil(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        nisoil(:,k)=nisoil(:,k)+unpack(sv(pind(n,1):pind(n,2))*real(casapool%nsoil(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
        psoil(:,k) =psoil(:,k) +unpack(sv(pind(n,1):pind(n,2))*real(casapool%psoil(pind(n,1):pind(n,2),k)),tmap(:,n),0.)
      end do
      !glai(:)=glai(:)+unpack(sv(pind(n,1):pind(n,2))*real(casamet%glai(pind(n,1):pind(n,2))),tmap(:,n),0.)
    end do

  end if ! icycle>0

  
  ! Calculate LAI and veg fraction diagnostics
  ! (needs to occur after CASA-CNP in case prognostic LAI is required)
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jyear,jmonth,jday,jhour,jmin,mp)
  vlai(:) = 0.
  do n = 1,maxnb
    vlai(:) = vlai(:) + unpack(sv(pind(n,1):pind(n,2))*veg%vlai(pind(n,1):pind(n,2)),tmap(:,n),0.)
  end do

  ! MJT suggestion to get an approx inital albedo (before cable is called)
  where ( land(1:ifull) )
    albvisnir(:,1) = albsoilsn(:,1)*(1.-sigmf) + 0.03*sigmf
    albvisnir(:,2) = albsoilsn(:,2)*(1.-sigmf) + 0.40*sigmf
  end where
  albvisdir(:) = albvisnir(:,1) ! To be updated by CABLE
  albvisdif(:) = albvisnir(:,1) ! To be updated by CABLE
  albnirdir(:) = albvisnir(:,2) ! To be updated by CABLE
  albnirdif(:) = albvisnir(:,2) ! To be updated by CABLE
  
else

  allocate( cveg(0) )
  call cable_biophysic_parm(cveg)
  deallocate( cveg )
  
end if
  
if (myid==0) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine cbmparm

subroutine cable_biophysic_parm(cveg)

use cc_mpi     ! CC MPI routines
use darcdf_m   ! Netcdf data
use infile     ! Input file routines

implicit none

integer k, numpft
integer, dimension(mp), intent(in) :: cveg
integer, dimension(1) :: nstart, ncount
integer, dimension(:), allocatable :: csiropft
real totdepth
real, dimension(:), allocatable :: hc, xfang, leaf_w, leaf_l, canst1
real, dimension(:), allocatable :: shelrb, extkn, vcmax, rpcoef
real, dimension(:), allocatable :: rootbeta, c4frac, vbeta
real, dimension(:), allocatable :: a1gs, d0gs, alpha, convex, cfrd
real, dimension(:), allocatable :: gswmin, conkc0, conko0, ekc, eko, g0, g1
real, dimension(:), allocatable :: zr, clitt
real, dimension(:,:), allocatable :: refl, taul, froot2

numpft = -1 ! missing flag

if ( myid==0 .and. lncveg==1 ) then
  call ccnf_inq_dimlen(ncidveg,'pft',numpft,failok=.true.)
end if
    
call ccmpi_bcast(numpft,0,comm_world)

if ( numpft<1 ) then
    
  ! default biophysical parameter tables
  if ( myid==0 ) then
    write(6,*) "Using default CABLE biophysical parameter tables"
  end if
  numpft = 18
  allocate( csiropft(numpft), hc(numpft), xfang(numpft), leaf_w(numpft), leaf_l(numpft) )
  allocate( canst1(numpft), shelrb(numpft), extkn(numpft), refl(numpft,2), taul(numpft,2) )
  allocate( vcmax(numpft), rpcoef(numpft), rootbeta(numpft), c4frac(numpft), froot2(numpft,ms) )
  allocate( vbeta(numpft) )
  allocate( a1gs(numpft), d0gs(numpft), alpha(numpft), convex(numpft), cfrd(numpft) )
  allocate( gswmin(numpft), conkc0(numpft), conko0(numpft), ekc(numpft), eko(numpft), g0(numpft), g1(numpft) )
  allocate( zr(numpft), clitt(numpft) )
  csiropft=(/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 2 /)
  hc    =(/   17.,  35.,  15.5,  20.,   0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567,  0.2, 6.017,  0.2,  0.2,  0.2,  0.2, 17. /)
  xfang =(/  0.01,  0.1,  0.01, 0.25,  0.01,  -0.3,  -0.3,  -0.3, -0.3, -0.3,  -0.3,  0.1,    0.,   0.,   0.,   0.,   0., 0.1 /)
  leaf_w=(/ 0.001, 0.05, 0.001, 0.08, 0.005,  0.01,  0.01,  0.01, 0.01, 0.01,  0.01, 0.03, 0.015, 0.00,   0.,   0.,   0., 0.05 /)
  leaf_l=(/ 0.055, 0.10, 0.040, 0.15, 0.100,  0.30,  0.30,  0.30, 0.30, 0.30,  0.30, 0.30, 0.242, 0.03, 0.03, 0.03, 0.03, 0.10 /)
  canst1=0.1
  shelrb=2.
  extkn=0.001
  refl(:,1)=(/ 0.062,0.076,0.056,0.092,0.100,0.110,0.100,0.117,0.100,0.090,0.108,0.055,0.091,0.238,0.143,0.143,0.159,0.076 /)
  refl(:,2)=(/ 0.302,0.350,0.275,0.380,0.400,0.470,0.400,0.343,0.400,0.360,0.343,0.190,0.310,0.457,0.275,0.275,0.305,0.350 /)
  taul(:,1)=(/ 0.050,0.050,0.045,0.050,0.050,0.070,0.100,0.080,0.100,0.090,0.075,0.023,0.059,0.039,0.023,0.023,0.026,0.050 /)
  taul(:,2)=(/ 0.100,0.250,0.144,0.250,0.240,0.250,0.150,0.124,0.150,0.225,0.146,0.198,0.163,0.189,0.113,0.113,0.113,0.250 /)
  vcmax=(/ 40.E-6,55.E-6,40.E-6,60.E-6,40.E-6,60.E-6,10.E-6,40.E-6,80.E-6,80.E-6,60.E-6,17.E-6,1.E-6,17.E-6,17.E-6,17.E-6, &
           17.E-6,55.E-6 /)
  rpcoef=0.0832
  rootbeta=(/ 0.943,0.962,0.966,0.961,0.964,0.943,0.943,0.943,0.961,0.961,0.943,0.975,0.961,0.961,0.961,0.961,0.961,0.962 /)
  c4frac=(/ 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0. /)
  vbeta=(/ 2., 2., 2., 2., 4., 4., 4., 4., 2., 2., 4., 4., 2., 4., 4., 4., 4., 2. /)
  a1gs=(/ 9., 9., 9., 9., 9., 9., 4., 9., 9., 4., 9., 9., 9., 9., 9., 9., 9., 9. /)
  d0gs=1500.
  alpha=(/ 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.2, 0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 /)
  convex=(/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.8, 0.01, 0.01, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 /)
  cfrd=(/ 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.015, 0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, &
          0.015, 0.015 /)
  gswmin=(/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.01, 0.04, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 /)
  conkc0=302.e-6
  conko0=256.e-3
  ekc=59430.
  eko=36000.
  g0=0.
  g1=(/ 2.346064, 4.114762, 2.346064, 4.447321, 4.694803, 5.248500, 1.616178, 2.222156, 5.789377, 1.616178, 5.248500, 5.248500, &
        0.000000, 5.248500, 5.248500, 5.248500, 5.248500, 2.346064 /)
  zr=(/ 1.8, 3., 2., 2., 2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.8, 3.1, 3., 1., 1., 1., 1., 3. /)
  clitt=(/ 20., 6., 10., 13., 2., 2., 0.3, 0.3, 0., 0., 2., 2., 0., 0., 0., 0., 0., 6. /) 
  
else

  ! user defined biophysical parameter tables  
  if ( myid==0 ) then
    write(6,*) "Using user defined CABLE biophysical parameter tables"
  end if
  allocate( csiropft(numpft), hc(numpft), xfang(numpft), leaf_w(numpft), leaf_l(numpft) )
  allocate( canst1(numpft), shelrb(numpft), extkn(numpft), refl(numpft,2), taul(numpft,2) )
  allocate( vcmax(numpft), rpcoef(numpft), rootbeta(numpft), c4frac(numpft), froot2(numpft,ms) )
  allocate( vbeta(numpft) )
  allocate( a1gs(numpft), d0gs(numpft), alpha(numpft), convex(numpft), cfrd(numpft) )
  allocate( gswmin(numpft), conkc0(numpft), conko0(numpft), ekc(numpft), eko(numpft), g0(numpft), g1(numpft) )
  allocate( zr(numpft), clitt(numpft) )

  if ( myid==0 ) then
    nstart(1) = 1
    ncount(1) = numpft
    call ccnf_get_vara(ncidveg,'csiropft',nstart,ncount,csiropft)
    call ccnf_get_vara(ncidveg,'hc',nstart,ncount,hc)
    call ccnf_get_vara(ncidveg,'xfang',nstart,ncount,xfang)
    call ccnf_get_vara(ncidveg,'leaf_w',nstart,ncount,leaf_w)
    call ccnf_get_vara(ncidveg,'leaf_l',nstart,ncount,leaf_l)
    call ccnf_get_vara(ncidveg,'canst1',nstart,ncount,canst1)
    call ccnf_get_vara(ncidveg,'shelrb',nstart,ncount,shelrb)
    call ccnf_get_vara(ncidveg,'extkn',nstart,ncount,extkn)
    call ccnf_get_vara(ncidveg,'rholeaf-vis',nstart,ncount,refl(:,1))
    call ccnf_get_vara(ncidveg,'rholeaf-nir',nstart,ncount,refl(:,2))
    call ccnf_get_vara(ncidveg,'tauleaf-vis',nstart,ncount,taul(:,1))
    call ccnf_get_vara(ncidveg,'tauleaf-nir',nstart,ncount,taul(:,2))
    call ccnf_get_vara(ncidveg,'vcmax',nstart,ncount,vcmax)
    call ccnf_get_vara(ncidveg,'rpcoef',nstart,ncount,rpcoef)
    call ccnf_get_vara(ncidveg,'rootbeta',nstart,ncount,rootbeta)
    call ccnf_get_vara(ncidveg,'c4frac',nstart,ncount,c4frac)
    call ccnf_get_vara(ncidveg,'vbeta',nstart,ncount,vbeta)
    call ccnf_get_vara(ncidveg,'a1gs',nstart,ncount,a1gs)
    call ccnf_get_vara(ncidveg,'d0gs',nstart,ncount,d0gs)
    call ccnf_get_vara(ncidveg,'alpha',nstart,ncount,alpha)
    call ccnf_get_vara(ncidveg,'convex',nstart,ncount,convex)
    call ccnf_get_vara(ncidveg,'cfrd',nstart,ncount,cfrd)
    call ccnf_get_vara(ncidveg,'gswmin',nstart,ncount,gswmin)
    call ccnf_get_vara(ncidveg,'conkc0',nstart,ncount,conkc0)
    call ccnf_get_vara(ncidveg,'conko0',nstart,ncount,conko0)
    call ccnf_get_vara(ncidveg,'ekc',nstart,ncount,ekc)
    call ccnf_get_vara(ncidveg,'eko',nstart,ncount,eko)
    call ccnf_get_vara(ncidveg,'g0',nstart,ncount,g0)
    call ccnf_get_vara(ncidveg,'g1',nstart,ncount,g1)
    call ccnf_get_vara(ncidveg,'zr',nstart,ncount,zr)
    call ccnf_get_vara(ncidveg,'clitt',nstart,ncount,clitt)
  end if
  call ccmpi_bcast(csiropft,0,comm_world)  
  call ccmpi_bcast(hc,0,comm_world)
  call ccmpi_bcast(xfang,0,comm_world)
  call ccmpi_bcast(leaf_w,0,comm_world)
  call ccmpi_bcast(leaf_l,0,comm_world)
  call ccmpi_bcast(canst1,0,comm_world)
  call ccmpi_bcast(shelrb,0,comm_world)
  call ccmpi_bcast(extkn,0,comm_world)
  call ccmpi_bcast(refl,0,comm_world)
  call ccmpi_bcast(taul,0,comm_world)
  call ccmpi_bcast(vcmax,0,comm_world)
  call ccmpi_bcast(rpcoef,0,comm_world)
  call ccmpi_bcast(rootbeta,0,comm_world)
  call ccmpi_bcast(c4frac,0,comm_world)
  call ccmpi_bcast(vbeta,0,comm_world)
  call ccmpi_bcast(a1gs,0,comm_world)
  call ccmpi_bcast(d0gs,0,comm_world)
  call ccmpi_bcast(alpha,0,comm_world)
  call ccmpi_bcast(convex,0,comm_world)
  call ccmpi_bcast(cfrd,0,comm_world)
  call ccmpi_bcast(gswmin,0,comm_world)
  call ccmpi_bcast(conkc0,0,comm_world)
  call ccmpi_bcast(conko0,0,comm_world)
  call ccmpi_bcast(ekc,0,comm_world)
  call ccmpi_bcast(eko,0,comm_world)
  call ccmpi_bcast(g0,0,comm_world)
  call ccmpi_bcast(g1,0,comm_world)
  call ccmpi_bcast(zr,0,comm_world)
  call ccmpi_bcast(clitt,0,comm_world)
  
end if

if ( mp>0 ) then

  ! froot is now calculated from soil depth and the new parameter rootbeta 
  ! according to Jackson et al. 1996, Oceologica, 108:389-411
  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + soil%zse(k)*100.
    froot2(:,k) = min(1.,1.-rootbeta(:)**totdepth)
  end do
  do k = ms-1, 2, -1
    froot2(:,k) = froot2(:,k) - froot2(:,k-1)
  end do
  froot2(:,ms) = 1.-sum(froot2(:,1:ms-1),2)
  
  ! Eva's method for ACCESS1.3
  !froot2(:,1)=0.05
  !froot2(:,2)=0.20
  !froot2(:,3)=0.20
  !froot2(:,4)=0.20
  !froot2(:,5)=0.20
  !froot2(:,6)=0.15

  if ( maxval(cveg)>numpft .or. minval(cveg)<1 ) then
    write(6,*) "ERROR: Invalid range of vegetation classes for CABLE"
    write(6,*) "cveg min,max           = ",minval(cveg),maxval(cveg)
    write(6,*) "Expected range min,max = ",1,numpft
    call ccmpi_abort(-1)
  end if

  veg%meth      = 1
  veg%iveg      = csiropft(cveg)
  veg%hc        = hc(cveg)
  veg%xfang     = xfang(cveg)  
  veg%dleaf     = sqrt(max(leaf_w(cveg)*leaf_l(cveg),1.e-20))
  veg%canst1    = canst1(cveg)
  veg%shelrb    = shelrb(cveg)
  veg%extkn     = extkn(cveg)
  veg%refl(:,1) = refl(cveg,1)
  veg%refl(:,2) = refl(cveg,2)  
  veg%taul(:,1) = taul(cveg,1)
  veg%taul(:,2) = taul(cveg,2)  
  veg%vcmax     = vcmax(cveg)
  veg%ejmax     = 2.*veg%vcmax
  veg%rpcoef    = rpcoef(cveg)
  do k = 1,ms
    veg%froot(:,k)=froot2(cveg,k)
  end do
  veg%frac4     = c4frac(cveg)
  veg%xalbnir   = 1. ! not used
  veg%vbeta     = vbeta(cveg)
  veg%a1gs      = a1gs(cveg)   
  veg%d0gs      = d0gs(cveg)
  veg%g0        = g0(cveg)
  veg%g1        = g1(cveg)
  veg%alpha     = alpha(cveg)
  veg%convex    = convex(cveg) 
  veg%cfrd      = cfrd(cveg)
  veg%gswmin    = gswmin(cveg)
  veg%conkc0    = conkc0(cveg)
  veg%conko0    = conko0(cveg)
  veg%ekc       = ekc(cveg)
  veg%eko       = eko(cveg)
  veg%zr        = zr(cveg)
  veg%clitt     = clitt(cveg)
  
  ! depeciated
  !veg%tminvj    = tminvj(veg%iveg)
  !veg%tmaxvj    = tmaxvj(veg%iveg)
  !veg%rp20      = rp20(veg%iveg)
  !veg%rs20      = rs20(veg%iveg)
  !veg%vegcf     = vegcf(veg%iveg)

end if

deallocate( csiropft, hc, xfang, leaf_w, leaf_l )
deallocate( canst1, shelrb, extkn, refl, taul )
deallocate( vcmax, rpcoef, rootbeta, c4frac, froot2 )
deallocate( vbeta )
deallocate( a1gs, d0gs, alpha, convex, cfrd )
deallocate( gswmin, conkc0, conko0, ekc, eko, g0, g1 )
deallocate( zr, clitt )

return
end subroutine cable_biophysic_parm
                   
! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,fvegprev,fveg,fvegnext,fvegnext2, &
                 cableformat)
  
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parmgeom_m

implicit none
  
character(len=*), intent(in) :: fveg,fvegprev,fvegnext,fvegnext2
integer, dimension(ifull,maxtile), intent(out) :: ivs
integer, dimension(ifull_g,maxtile) :: ivsg  
integer, dimension(3) :: spos,npos
integer n,iq,ilx,jlx,iad 
integer ncidx,iernc,varid,ndims
real, dimension(ifull,maxtile), intent(out) :: svs, vlinprev, vlin, vlinnext, vlinnext2
real, dimension(ifull_g,maxtile) :: svsg, vling
real, dimension(ifull_g) :: savannafrac_g
real rlong0x,rlat0x,schmidtx,dsx,ra,rb
real cablever
real, intent(out) :: cableformat
character(len=47) header  
character(len=7) vname
real, parameter :: cableversion = 3939. ! version id for input data
logical tst

cableformat=0.

write(6,*) "Reading land-use maps for CABLE"
if ( lncveg==1 ) then
  spos(1:3) = 1
  npos(1) = il_g
  npos(2) = 6*il_g
  npos(3) = 1
  call ccnf_inq_dimlen(ncidveg,'longitude',ilx)
  call ccnf_inq_dimlen(ncidveg,'latitude',jlx)
  call ccnf_get_attg(ncidveg,'lon0',rlong0x)
  call ccnf_get_attg(ncidveg,'lat0',rlat0x)
  call ccnf_get_attg(ncidveg,'schmidt',schmidtx)
  if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
      abs(schmidtx-schmidt)>1.e-20) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  call ccnf_get_attg(ncidveg,'cableversion',cablever,ierr=iernc)
  if (iernc/=0) then
    write(6,*) "Missing version of CABLE data"
    write(6,*) "Regenerate land-use data with up-to-date version of igbpveg"
    call ccmpi_abort(-1)
  end if
  if (abs(cablever-cableversion)>1.e-20) then
    write(6,*) "Wrong version of CABLE data"
    write(6,*) "Expecting ",cableversion
    write(6,*) "Found     ",cablever
    call ccmpi_abort(-1)
  end if
  call ccnf_get_attg(ncidveg,'cableformat',cableformat,ierr=iernc)
  if ( iernc/=0 ) then
    cableformat=0.
  end if
  do n = 1,maxtile
    write(vname,"(A,I1.1)") "lai",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    write(vname,"(A,I1.1)") "vegt",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n)) 
    ivsg(:,n)=nint(svsg(:,n))
    write(vname,"(A,I1.1)") "vfrac",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n))
  end do
  if ( abs(cableformat-1.)>1.e-20 ) then
    vname="savanna"
    call ccnf_inq_varid(ncidveg,vname,varid,tst)
    if ( .not.tst ) then
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),savannafrac_g)
      do n = 1,maxtile
        where ( savannafrac_g>0.5*svsg(:,n) .and. ivsg(:,n)==2 )
          ivsg(:,n) = 18
        end where
      end do
    end if
  end if
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if ( fvegprev/=' ' .and. fvegnext/=' ' ) then
    call ccnf_open(fvegprev,ncidx,iernc)
    if ( iernc/=0 ) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do n = 1,maxtile
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinprev,vling)
    call ccnf_open(fvegnext,ncidx,iernc)
    if (iernc/=0) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    do n=1,maxtile
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinnext,vling)
  else
    vlinprev = -1.
    vlinnext = -1.    
  end if
  if ( fvegnext2/=' ' ) then
    call ccnf_open(fvegnext2,ncidx,iernc)
    if ( iernc/=0 ) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegnext2)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do n = 1,maxtile
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinnext2,vling)
  else
    vlinnext2 = -1.
  end if
  
else
  open(87,file=fveg,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
      abs(schmidtx-schmidt)>1.e-20) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  do iq = 1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if ( fvegprev/=' ' .and. fvegnext/=' ' ) then
    open(87,file=fvegprev,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinprev,vling)
    open(87,file=fvegnext,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if (ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinnext,vling)
  else
    vlinprev=-1.
    vlinnext=-1.    
  end if
  if ( fvegnext2/=' ' ) then
    open(87,file=fvegnext2,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-20.or.abs(rlat0x-rlat0)>1.e-20.or. &
        abs(schmidtx-schmidt)>1.e-20) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinnext2,vling)
  else
    vlinnext2=-1.
  end if
  
end if

call ccmpi_bcast(cableformat,0,comm_world)

return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlinprev,vlin,vlinnext,vlinnext2,fvegprev,fvegnext,fvegnext2, &
                 cableformat)
  
use cc_mpi
use newmpar_m
  
implicit none

character(len=*), intent(in) :: fvegprev,fvegnext,fvegnext2
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs, vlinprev, vlin, vlinnext, vlinnext2
real, intent(out) :: cableformat

cableformat = 0.

call ccmpi_distribute(ivs)
call ccmpi_distribute(svs)
call ccmpi_distribute(vlin)
if ( fvegprev/=' ' .and. fvegnext/=' ' ) then
  call ccmpi_distribute(vlinprev)
  call ccmpi_distribute(vlinnext)
else
  vlinprev = -1.
  vlinnext = -1.
end if    
if ( fvegnext2/=' ' ) then
  call ccmpi_distribute(vlinnext2)
else
  vlinnext2 = -1.
end if    
  
call ccmpi_bcast(cableformat,0,comm_world)

return
end subroutine vegtb

! *************************************************************************************  
! This subroutine loads CABLE tile data
subroutine defaulttile

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
integer k, n

if ( mp>0 ) then
  do n = 1,maxnb
    do k = 1,ms
      ssnow%tgg(pind(n,1):pind(n,2),k)   = pack(tgg(:,k),  tmap(:,n))
      ssnow%wb(pind(n,1):pind(n,2),k)    = pack(wb(:,k),   tmap(:,n))
      ssnow%wb(pind(n,1):pind(n,2),k)    = max( ssnow%wb(pind(n,1):pind(n,2),k), real(0.5*soil%swilt(pind(n,1):pind(n,2)),r_2) )
      ssnow%wb(pind(n,1):pind(n,2),k)    = min( ssnow%wb(pind(n,1):pind(n,2),k), real(soil%sfc(pind(n,1):pind(n,2)),r_2) )
      ssnow%wbice(pind(n,1):pind(n,2),k) = pack(wbice(:,k),tmap(:,n))
    end do
    do k = 1,3
      ssnow%tggsn(pind(n,1):pind(n,2),k)  = pack(tggsn(:,k),tmap(:,n))
      ssnow%smass(pind(n,1):pind(n,2),k)  = pack(smass(:,k),tmap(:,n))
      ssnow%ssdn(pind(n,1):pind(n,2),k)   = pack(ssdn(:,k), tmap(:,n))
      ssnow%sdepth(pind(n,1):pind(n,2),k) = pack(smass(:,k)/ssdn(:,k),tmap(:,n))
      ssnow%sconds(pind(n,1):pind(n,2),k) = 0.3
    end do      
    ssnow%ssdnn(pind(n,1):pind(n,2))  = pack(ssdnn, tmap(:,n))
    ssnow%isflag(pind(n,1):pind(n,2)) = pack(isflag,tmap(:,n))
    ssnow%snowd(pind(n,1):pind(n,2))  = pack(snowd, tmap(:,n))
    ssnow%snage(pind(n,1):pind(n,2))  = pack(snage, tmap(:,n))
  end do
  ssnow%rtsoil = 50.
  canopy%cansto = 0.
  canopy%ga = 0.
  canopy%us = 0.01
  ssnow%pudsto = 0.
  ssnow%wetfac = 0.
  ssnow%osnowd = ssnow%snowd
  if ( icycle==0 ) then
    !do n = 1,maxnb
    !  do k = 1,ncp
    !    bgc%cplant(pind(n,1):pind(n,2),k) = pack(cplant(:,k),tmap(:,n))
    !  end do
    !  do k = 1,ncs
    !    bgc%csoil(pind(n,1):pind(n,2),k) = pack(csoil(:,k),tmap(:,n))
    !  end do
    !end do
  else
    do n = 1,maxnb
      do k = 1,mplant
        casapool%cplant(pind(n,1):pind(n,2),k) = pack(cplant(:,k), tmap(:,n))
        casapool%nplant(pind(n,1):pind(n,2),k) = pack(niplant(:,k),tmap(:,n))
        casapool%pplant(pind(n,1):pind(n,2),k) = pack(pplant(:,k), tmap(:,n))
      end do
      do k = 1,mlitter
        casapool%clitter(pind(n,1):pind(n,2),k) = pack(clitter(:,k), tmap(:,n))
        casapool%nlitter(pind(n,1):pind(n,2),k) = pack(nilitter(:,k),tmap(:,n))
        casapool%plitter(pind(n,1):pind(n,2),k) = pack(plitter(:,k), tmap(:,n))
      end do
      do k = 1,msoil
        casapool%csoil(pind(n,1):pind(n,2),k) = pack(csoil(:,k), tmap(:,n))
        casapool%nsoil(pind(n,1):pind(n,2),k) = pack(nisoil(:,k),tmap(:,n))
        casapool%psoil(pind(n,1):pind(n,2),k) = pack(psoil(:,k), tmap(:,n))
      end do
      !casamet%glai(pind(n,1):pind(n,2)) = pack(glai,tmap(:,n))
    end do
  end if
  canopy%fhs_cor = 0.
  canopy%fes_cor = 0.
  
  call fixtile
  
end if

return
end subroutine defaulttile

subroutine newcbmwb

use soilsnow_m

implicit none

integer k, n

if ( mp>0 ) then
  do n = 1,maxnb
    do k = 1,ms
      ssnow%wb(pind(n,1):pind(n,2),k)    = pack(wb(:,k),   tmap(:,n))
    end do
  end do
end if

return
end subroutine newcbmwb

subroutine loadtile(usedefault)

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use nsibd_m, only : sigmf
use parm_m
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
logical, intent(in), optional :: usedefault
integer k, n, ierr, idv
integer jyear,jmonth,jday,jhour,jmin,mins
real, dimension(ifull) :: dat
real, dimension(ifull,ms) :: datms
real, dimension(ifull,3) :: dat3
real, dimension(ifull,mplant) :: datmplant
real, dimension(ifull,mlitter) :: datmlitter
real, dimension(ifull,msoil) :: datmsoil
real fjd
logical tst
logical defaultmode
character(len=26) vname
character(len=7) testname

! force CABLE to use generic input for all tiles
! if usedefault = defaultmode = .true.
defaultmode = .false.
if ( present(usedefault) ) then
  defaultmode = usedefault
end if

! check that CABLE data exists in restart file
! and communicate the result to all processors
! as not all processors are assigned an input file
ierr = 1
if ( io_in==1 .and. .not.defaultmode ) then
  if ( myid==0 .or. pfall ) then
    write(testname,'("t",I1.1,"_tgg1")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr = 0
    end if
  end if
  if ( .not.pfall ) then
    call ccmpi_bcast(ierr,0,comm_world)
  end if
end if
  
! Cannot locate tile data, use diagnostic data instead
if ( ierr/=0 ) then
  if ( myid==0 ) write(6,*) "Use gridbox averaged data to initialise CABLE"
  call defaulttile
else
  ! Located CABLE tile data
  if ( myid==0 ) write(6,*) "Use tiled data to initialise CABLE"
  do n = 1,maxtile
    write(vname,'("t",I1.1,"_tgg")') n
    call histrd4(iarchi-1,ierr,vname,il_g,ms,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        ssnow%tgg(pind(n,1):pind(n,2),k) = pack(datms(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_wb")') n
    call histrd4(iarchi-1,ierr,vname,il_g,ms,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        ssnow%wb(pind(n,1):pind(n,2),k) = pack(datms(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_wbice")') n
    call histrd4(iarchi-1,ierr,vname,il_g,ms,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        ssnow%wbice(pind(n,1):pind(n,2),k) = pack(datms(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_tggsn")') n
    call histrd4(iarchi-1,ierr,vname,il_g,3,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        ssnow%tggsn(pind(n,1):pind(n,2),k) = pack(dat3(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_smass")') n
    call histrd4(iarchi-1,ierr,vname,il_g,3,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        ssnow%smass(pind(n,1):pind(n,2),k) = pack(dat3(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_ssdn")') n
    call histrd4(iarchi-1,ierr,vname,il_g,3,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        ssnow%ssdn(pind(n,1):pind(n,2),k) = pack(dat3(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n
    call histrd4(iarchi-1,ierr,vname,il_g,3,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        ssnow%sdepth(pind(n,1):pind(n,2),k) = pack(dat3(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_sconds")') n
    call histrd4(iarchi-1,ierr,vname,il_g,3,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        ssnow%sconds(pind(n,1):pind(n,2),k) = pack(dat3(:,k),tmap(:,n))
      end do
    end if
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%ssdnn(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_sflag")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%isflag(pind(n,1):pind(n,2))=nint(pack(dat,tmap(:,n)))
    write(vname,'("t",I1.1,"_snd")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%snowd(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_osnd")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%osnowd(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_snage")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%snage(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%rtsoil(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_cansto")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) canopy%cansto(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_us")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) canopy%us(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_pudsto")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%pudsto(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_wetfac")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) ssnow%wetfac(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("t",I1.1,"_ga")') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n<=maxnb ) canopy%ga(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    if ( icycle==0 ) then
      !write(vname,'("t",I1.1,"_cplant")') n
      !call histrd4(iarchi-1,ierr,vname,il_g,ncp,datncp(:,1:ncp),ifull)
      !if ( n<=maxnb ) then
      !  do k = 1,ncp
      !    bgc%cplant(pind(n,1):pind(n,2),k) = pack(datncp(:,k),tmap(:,n))
      !  end do
      !end if
      !write(vname,'("t",I1.1,"_csoil")') n
      !call histrd4(iarchi-1,ierr,vname,il_g,ncs,datncs(:,1:ncs),ifull)
      !if ( n<=maxnb ) then
      !  do k = 1,ncs
      !    bgc%csoil(pind(n,1):pind(n,2),k) = pack(datncs(:,k),tmap(:,n))
      !  end do
      !end if
    else
      write(vname,'("t",I1.1,"_cplant")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mplant,datmplant(:,1:mplant),ifull)
      if ( n<=maxnb ) then
        do k = 1,mplant
          casapool%cplant(pind(n,1):pind(n,2),k) = pack(datmplant(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_nplant")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mplant,datmplant(:,1:mplant),ifull)
      if ( n<=maxnb ) then
        do k = 1,mplant
          casapool%nplant(pind(n,1):pind(n,2),k) = pack(datmplant(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_pplant")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mplant,datmplant(:,1:mplant),ifull)
      if ( n<=maxnb ) then
        do k = 1,mplant
          casapool%pplant(pind(n,1):pind(n,2),k) = pack(datmplant(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_clitter")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mlitter,datmlitter(:,1:mlitter),ifull)
      if ( n<=maxnb ) then
        do k = 1,mlitter
          casapool%clitter(pind(n,1):pind(n,2),k) = pack(datmlitter(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_nlitter")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mlitter,datmlitter(:,1:mlitter),ifull)
      if ( n<=maxnb ) then
        do k = 1,mlitter
          casapool%nlitter(pind(n,1):pind(n,2),k) = pack(datmlitter(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_plitter")') n
      call histrd4(iarchi-1,ierr,vname,il_g,mlitter,datmlitter(:,1:mlitter),ifull)
      if ( n<=maxnb ) then
        do k = 1,mlitter
          casapool%plitter(pind(n,1):pind(n,2),k) = pack(datmlitter(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_csoil")') n
      call histrd4(iarchi-1,ierr,vname,il_g,msoil,datmsoil(:,1:msoil),ifull)
      if ( n<=maxnb ) then
        do k = 1,msoil
          casapool%csoil(pind(n,1):pind(n,2),k) = pack(datmsoil(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_nsoil")') n
      call histrd4(iarchi-1,ierr,vname,il_g,msoil,datmsoil(:,1:msoil),ifull)
      if ( n<=maxnb ) then
        do k = 1,msoil
          casapool%nsoil(pind(n,1):pind(n,2),k) = pack(datmsoil(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_psoil")') n
      call histrd4(iarchi-1,ierr,vname,il_g,msoil,datmsoil(:,1:msoil),ifull)
      if ( n<=maxnb ) then
        do k = 1,msoil
          casapool%psoil(pind(n,1):pind(n,2),k) = pack(datmsoil(:,k),tmap(:,n))
        end do
      end if
      write(vname,'("t",I1.1,"_glai")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casamet%glai(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_phenphase")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) phen%phase(pind(n,1):pind(n,2)) = nint(pack(dat,tmap(:,n)))
      write(vname,'("t",I1.1,"_phenaphen")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) phen%aphen(pind(n,1):pind(n,2)) = nint(pack(dat,tmap(:,n)))
      write(vname,'("t",I1.1,"_clabile")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casapool%clabile(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casapool%nsoilmin(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_psoillab")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casapool%psoillab(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casapool%psoilsorb(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) casapool%psoilocc(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
    end if
    if ( cable_climate==1 ) then
      !write(vname,'("t",I1.1,"_climateevapPT")') n
      !call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      !if ( n<=maxnb ) climate%evap_PT(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      !write(vname,'("t",I1.1,"_climateaevap")') n
      !call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      !if ( n<=maxnb ) climate%aevap(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climatemtempmax")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%mtemp_max(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climatemtempmin")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%mtemp_min(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climateqtempmax")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%qtemp_max(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%qtemp_max_last_year(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climategdd0")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%gdd0(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climategdd5")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%gdd5(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climateagdd0")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%agdd0(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climateagdd5")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%agdd5(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("t",I1.1,"_climatechilldays")') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n<=maxnb ) climate%chilldays(pind(n,1):pind(n,2)) = nint(pack(dat,tmap(:,n)))
      do k = 1,91
        write(vname,'("t",I1.1,"_climatedtemp",I2.2)') 1,k ! all tiles have the same data
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n<=maxnb ) then
          climate%dtemp_91(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
          if ( k>60 ) then
            climate%dtemp_31(pind(n,1):pind(n,2),k-60) = pack(dat,tmap(:,n))
          end if
        end if
      end do
      do k = 1,31
        write(vname,'("t",I1.1,"_climatedmoist",I2.2)') 1,k
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n<=maxnb ) climate%dmoist_31(pind(n,1):pind(n,2),k) = nint(pack(dat,tmap(:,n)))
      end do
    end if
    ! CABLE correction terms
    !write(vname,'("t",I1.1,"_fhscor")') n
    !call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    !if ( n<=maxnb ) canopy%fhs_cor(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
    !write(vname,'("t",I1.1,"_fescor")') n
    !call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    !if ( n<=maxnb ) canopy%fes_cor(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
  end do
  ! albvisdir, albvisdif, albnirdir, albnirdif are used when nrad=5
  vname = 'albvisdir'
  call histrd1(iarchi-1,ierr,vname,il_g,albvisdir,ifull)
  vname = 'albvisdif'
  call histrd1(iarchi-1,ierr,vname,il_g,albvisdif,ifull)
  vname = 'albnirdir'
  call histrd1(iarchi-1,ierr,vname,il_g,albnirdir,ifull)
  vname = 'albnirdif'
  call histrd1(iarchi-1,ierr,vname,il_g,albnirdif,ifull)
  ! albvis and albnir are used when nrad=4
  vname = 'albvis'
  call histrd1(iarchi-1,ierr,vname,il_g,albvisnir(:,1),ifull)
  vname = 'albnir'
  call histrd1(iarchi-1,ierr,vname,il_g,albvisnir(:,2),ifull)
  
  call fixtile
  
end if

! Calculate LAI and veg fraction diagnostics
vlai(:) = 0.
sigmf(:) = 0.
if ( mp>0 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jyear,jmonth,jday,jhour,jmin,mp)
  do n = 1,maxnb
    vlai(:) = vlai(:) + unpack(sv(pind(n,1):pind(n,2))*veg%vlai(pind(n,1):pind(n,2)),tmap(:,n),0.)
  end do
end if

return
end subroutine loadtile

subroutine fixtile

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
integer k
real totdepth
 
! Some fixes for rounding errors
if ( mp>0 ) then

  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + soil%zse(k)*100.
  enddo

  ssnow%wb = max(ssnow%wb,0._r_2)
  ssnow%wbice = min( max(ssnow%wbice, 0._r_2), ssnow%wb )
  ssnow%smass = max(ssnow%smass,0.)
  ssnow%rtsoil = max(ssnow%rtsoil,0.)
  ssnow%snowd = max(ssnow%snowd,0.)
  ssnow%osnowd = max(ssnow%osnowd,0.)
  ssnow%wetfac = min(max(ssnow%wetfac,0.),1.)
  canopy%cansto = max(canopy%cansto,0.)

  ssnow%wbtot = 0.
  ssnow%wbtot1 = 0.
  ssnow%wbtot2 = 0.
  ssnow%tggav = 0.
  do k = 1,ms
    ssnow%wbtot = ssnow%wbtot+ssnow%wb(:,k)*1000.0*soil%zse(k)
    ssnow%tggav = ssnow%tggav+soil%zse(k)*ssnow%tgg(:,k)/(totdepth/100.)
    ssnow%gammzz(:,k) = max((1.-soil%ssat)*soil%css* soil%rhosoil                     &
        + real(ssnow%wb(:,k)-ssnow%wbice(:,k))*4.218e3* 1000.                         &
        + real(ssnow%wbice(:,k))*2.100e3*1000.*0.9,soil%css*soil%rhosoil)*soil%zse(k) &
        + (1.-ssnow%isflag)*2090.0*ssnow%snowd
  end do

  if ( icycle == 0 ) then
    !bgc%cplant = max(bgc%cplant,0.)
    !bgc%csoil = max(bgc%csoil,0.)
  else
    casapool%cplant     = max(0._r_2,casapool%cplant)
    casapool%clitter    = max(0._r_2,casapool%clitter)
    casapool%csoil      = max(0._r_2,casapool%csoil)
    casabal%cplantlast(1:mp,1:mplant)   = casapool%cplant(1:mp,1:mplant)
    casabal%clitterlast(1:mp,1:mlitter) = casapool%clitter(1:mp,1:mlitter)
    casabal%csoillast(1:mp,1:msoil)     = casapool%csoil(1:mp,1:msoil)
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal      = 0.
    casabal%FCgppyear    = 0.
    casabal%FCrpyear     = 0
    casabal%FCrmleafyear = 0.
    casabal%FCrmwoodyear = 0.
    casabal%FCrmrootyear = 0.
    casabal%FCrgrowyear  = 0.
    casabal%FCnppyear    = 0.
    casabal%FCrsyear     = 0.
    casabal%FCneeyear    = 0.
    casabal%dCdtyear     = 0.
    casapool%nplant      = max(1.e-6_r_2,casapool%nplant)
    casapool%nlitter     = max(1.e-6_r_2,casapool%nlitter)
    casapool%nsoil       = max(1.e-6_r_2,casapool%nsoil)
    casapool%nsoilmin    = max(1.e-6_r_2,casapool%nsoilmin)
    casabal%nplantlast(1:mp,1:mplant)   = casapool%nplant(1:mp,1:mplant)
    casabal%nlitterlast(1:mp,1:mlitter) = casapool%nlitter(1:mp,1:mlitter)
    casabal%nsoillast(1:mp,1:msoil)     = casapool%nsoil(1:mp,1:msoil)      
    casabal%nsoilminlast= casapool%nsoilmin
    casabal%sumnbal      = 0.
    casabal%FNdepyear    = 0.
    casabal%FNfixyear    = 0.
    casabal%FNsnetyear   = 0.
    casabal%FNupyear     = 0.
    casabal%FNleachyear  = 0.
    casabal%FNlossyear   = 0.
    casapool%pplant      = max(1.0e-7_r_2,casapool%pplant)
    casapool%plitter     = max(1.0e-7_r_2,casapool%plitter)
    casapool%psoil       = max(1.0e-7_r_2,casapool%psoil)
    casapool%Psoillab    = max(1.0e-7_r_2,casapool%psoillab)
    casapool%psoilsorb   = max(1.0e-7_r_2,casapool%psoilsorb)
    casapool%psoilocc    = max(1.0e-7_r_2,casapool%psoilocc)
    casabal%pplantlast(1:mp,1:mplant)   = casapool%pplant(1:mp,1:mplant)
    casabal%plitterlast(1:mp,1:mlitter) = casapool%plitter(1:mp,1:mlitter)
    casabal%psoillast(1:mp,1:msoil)     = casapool%psoil(1:mp,1:msoil)       
    casabal%psoillablast = casapool%psoillab
    casabal%psoilsorblast= casapool%psoilsorb
    casabal%psoilocclast = casapool%psoilocc
    casabal%sumpbal      = 0.
    casabal%FPweayear    = 0.
    casabal%FPdustyear   = 0.
    casabal%FPsnetyear   = 0.
    casabal%FPupyear     = 0.
    casabal%FPleachyear  = 0.
    casabal%FPlossyear   = 0.
  end if
end if
  
return
end subroutine fixtile

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetiledef(idnc,local,idim,isize)

use carbpools_m
use cc_mpi, only : myid
use infile
use newmpar_m
  
implicit none
  
integer, intent(in) :: idnc, isize
integer k,n
integer, dimension(4), intent(in) :: idim  
character(len=26) vname
character(len=40) lname
logical, intent(in) :: local
  
if (myid==0.or.local) then
  if (myid==0) then
    write(6,*) "Defining CABLE tile data"
  end if
  do n=1,maxtile
    do k=1,ms
      write(lname,'("Soil temperature tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
      call attrib(idnc,idim,isize,vname,lname,'K',100.,400.,0,-1)
      write(lname,'("Soil moisture tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_wb",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'m3/m3',0.,2.6,0,-1)
      write(lname,'("Soil ice tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_wbice",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'m3/m3',0.,2.6,0,-1)
    end do
    do k=1,3
      write(lname,'("Snow temperature tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
      call attrib(idnc,idim,isize,vname,lname,'K',100.,400.,0,-1)
      write(lname,'("Snow mass tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_smass",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'K',0.,650.,0,-1)
      write(lname,'("Snow density tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'kg/m3',0.,650.,0,-1)
      write(lname,'("Snow depth tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'mm',0.,6500.,0,-1)
      write(lname,'("Snow sconds tile ",I1.1," lev ",I1.1)') n,k
      write(vname,'("t",I1.1,"_sconds",I1.1)') n,k 
      call attrib(idnc,idim,isize,vname,lname,'none',0.,6.5,0,-1)
    end do
    write(lname,'("Snow ssdnn tile ",I1.1)') n
    write(vname,'("t",I1.1,"_ssdnn")') n
    call attrib(idnc,idim,isize,vname,lname,'kg/m3',0.,650.,0,-1)
    write(lname,'("Snow flag tile ",I1.1)') n
    write(vname,'("t",I1.1,"_sflag")') n
    call attrib(idnc,idim,isize,vname,lname,'mm',0.,6.5,0,-1)
    write(lname,'("Snow depth tile ",I1.1)') n
    write(vname,'("t",I1.1,"_snd")') n
    call attrib(idnc,idim,isize,vname,lname,'mm',0.,6500.,0,-1)  ! -1=long
    write(lname,'("Old snow depth tile ",I1.1)') n
    write(vname,'("t",I1.1,"_osnd")') n
    call attrib(idnc,idim,isize,vname,lname,'mm',0.,6500.,0,-1)  ! -1=long
    write(lname,'("Snow age tile ",I1.1)') n
    write(vname,'("t",I1.1,"_snage")') n
    call attrib(idnc,idim,isize,vname,lname,'none',0.,26.,0,-1)
    write(lname,'("Soil turbulent resistance tile ",I1.1)') n
    write(vname,'("t",I1.1,"_rtsoil")') n
    call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3e5,0,-1)
    write(lname,'("cansto tile ",I1.1)') n
    write(vname,'("t",I1.1,"_cansto")') n
    call attrib(idnc,idim,isize,vname,lname,'none',0.,13.,0,-1)
    write(lname,'("us tile ",I1.1)') n
    write(vname,'("t",I1.1,"_us")') n
    call attrib(idnc,idim,isize,vname,lname,'m/s',0.,13.,0,-1)    
    write(lname,'("pudsto tile ",I1.1)') n
    write(vname,'("t",I1.1,"_pudsto")') n
    call attrib(idnc,idim,isize,vname,lname,'none',0.,13.,0,-1)
    write(lname,'("wetfac tile ",I1.1)') n
    write(vname,'("t",I1.1,"_wetfac")') n
    call attrib(idnc,idim,isize,vname,lname,'none',0.,6.5,0,-1)
    write(lname,'("ga tile ",I1.1)') n
    write(vname,'("t",I1.1,"_ga")') n
    call attrib(idnc,idim,isize,vname,lname,'none',-6500.,6500.,0,-1)
    if ( icycle==0 ) then
      !write(lname,'("Carbon leaf pool tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_cplant1")') n    
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      !write(lname,'("Carbon wood pool tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_cplant2")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      !write(lname,'("Carbon root pool tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_cplant3")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      !write(lname,'("Carbon soil fast pool tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_csoil1")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      !write(lname,'("Carbon soil slow pool tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_csoil2")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
    else
      do k = 1,mplant
        write(lname,'("C leaf pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N leaf pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P leaf pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      end do
      do k = 1,mlitter
        write(lname,'("C litter pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N litter pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P litter pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      end do
      do k = 1,msoil
        write(lname,'("C soil pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N soil pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P soil pool tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
        call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      end do
      write(lname,'("Prognostic LAI tile ",I1.1)') n
      write(vname,'("t",I1.1,"_glai")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Leaf phenology phase tile ",I1.1)') n
      write(vname,'("t",I1.1,"_phenphase")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Leaf phenology rainfall ",I1.1)') n
      write(vname,'("t",I1.1,"_phenaphen")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("C labile tile ",I1.1)') n
      write(vname,'("t",I1.1,"_clabile")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("N soilmin tile ",I1.1)') n
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soillab tile ",I1.1)') n
      write(vname,'("t",I1.1,"_psoillab")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soilsorb tile ",I1.1)') n
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soilocc tile ",I1.1)') n
      write(vname,'("t",I1.1,"_psoilocc")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
    end if
    if ( cable_climate==1 ) then
      !write(lname,'("climate evapPT tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_climateevapPT")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      !write(lname,'("climate aevap tile ",I1.1)') n
      !write(vname,'("t",I1.1,"_climateaevap")') n
      !call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate mtempmax tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climatemtempmax")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate mtempmin tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climatemtempmin")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate qtempmax tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climateqtempmax")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate qtempmaxlastyear tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate gdd0 tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climategdd0")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate gdd5 tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climategdd5")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate agdd0 tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climateagdd0")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate agdd5 tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climateagdd5")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("climate chilldays tile ",I1.1)') n
      write(vname,'("t",I1.1,"_climatechilldays")') n
      call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
      if ( n==1 ) then
        do k = 1,91
          write(lname,'("climate dtemp tile ",I1.1," lev ",I2.2)') n,k
          write(vname,'("t",I1.1,"_climatedtemp",I2.2)') n,k
          call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        end do
        do k = 1,31
          write(lname,'("climate dmoist tile ",I1.1," lev ",I2.2)') n,k
          write(vname,'("t",I1.1,"_climatedmoist",I2.2)') n,k
          call attrib(idnc,idim,isize,vname,lname,'none',0.,65000.,0,-1)
        end do
      end if
    end if
    !write(lname,'("Sensible correction term ",I1.1)') n
    !write(vname,'("t",I1.1,"_fhscor")') n
    !call attrib(idnc,idim,isize,vname,lname,'W/m2',-3000.,3000.,0,-1)
    !write(lname,'("Latent correction term ",I1.1)') n
    !write(vname,'("t",I1.1,"_fescor")') n
    !call attrib(idnc,idim,isize,vname,lname,'W/m2',-3000.,3000.,0,-1)
  end do
  lname='DIR VIS albedo'
  vname='albvisdir'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
  lname='DIF VIS albedo'
  vname='albvisdif'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
  lname='DIR NIR albedo'
  vname='albnirdir'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
  lname='DIF NIR albedo'
  vname='albnirdif'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
  lname='VIS albedo'
  vname='albvis'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
  lname='NIR albedo'
  vname='albnir'
  call attrib(idnc,idim,isize,vname,lname,'none',0.,1.3,0,-1)
end if
  
return
end subroutine savetiledef

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,iarch)

use carbpools_m
use infile
use newmpar_m
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
integer, intent(in) :: idnc,iarch
integer k,n
real, dimension(ifull) :: dat
character(len=26) vname
logical, intent(in) :: local
  
do n = 1,maxtile  ! tile
  do k = 1,ms     ! soil layer
    dat=tgg(:,k)
    if (n<=maxnb) dat=unpack(ssnow%tgg(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wb(:,k)
    if (n<=maxnb) dat=unpack(real(ssnow%wb(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_wb",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wbice(:,k)
    if (n<=maxnb) dat=unpack(real(ssnow%wbice(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_wbice",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  do k = 1,3 ! snow layer
    dat=tggsn(:,k)
    if (n<=maxnb) dat=unpack(ssnow%tggsn(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=smass(:,k)
    if (n<=maxnb) dat=unpack(ssnow%smass(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_smass",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=ssdn(:,k)
    if (n<=maxnb) dat=unpack(ssnow%ssdn(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=snowd/3.
    if (n<=maxnb) dat=unpack(ssnow%sdepth(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.2
    if (n<=maxnb) dat=unpack(ssnow%sconds(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_sconds",I1.1)') n,k
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  dat=ssdnn
  if (n<=maxnb) dat=unpack(ssnow%ssdnn(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_ssdnn")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=real(isflag)
  if (n<=maxnb) dat=unpack(real(ssnow%isflag(pind(n,1):pind(n,2))),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_sflag")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snowd
  if (n<=maxnb) dat=unpack(ssnow%snowd(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_snd")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snowd
  if (n<=maxnb) dat=unpack(ssnow%osnowd(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_osnd")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snage
  if (n<=maxnb) dat=unpack(ssnow%snage(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_snage")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=100.
  if (n<=maxnb) dat=unpack(ssnow%rtsoil(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_rtsoil")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)   
  dat=0.
  if (n<=maxnb) dat=unpack(canopy%cansto(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_cansto")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.01 ! ustar
  if (n<=maxnb) dat=unpack(canopy%us(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_us")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)  
  dat=0.
  if (n<=maxnb) dat=unpack(ssnow%pudsto(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_pudsto")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (n<=maxnb) dat=unpack(ssnow%wetfac(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_wetfac")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (n<=maxnb) dat=unpack(canopy%ga(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("t",I1.1,"_ga")') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  if ( icycle==0 ) then
    !do k = 1,ncp
    !  dat=cplant(:,k)
    !  if (n<=maxnb) dat=unpack(bgc%cplant(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    !  write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
    !  call histwrt3(dat,vname,idnc,iarch,local,.true.)    
    !end do
    !do k = 1,ncs
    !  dat=csoil(:,k)
    !  if (n<=maxnb) dat=unpack(bgc%csoil(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    !  write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
    !  call histwrt3(dat,vname,idnc,iarch,local,.true.)
    !end do
  else
    ! MJT notes - possible rounding error when using CASA CNP  
    do k = 1,mplant     
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%cplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%pplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,mlitter
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%clitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nlitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%plitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,msoil
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%csoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nsoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%psoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    dat=0.
    if (n<=maxnb) dat=unpack(real(casamet%glai(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_glai")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(phen%phase(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_phenphase")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(phen%aphen(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_phenaphen")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%clabile(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_clabile")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%nsoilmin(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_nsoilmin")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoillab(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_psoillab")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoilsorb(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_psoilsorb")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoilocc(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_psoilocc")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end if
  if ( cable_climate==1 ) then
    !dat=0.
    !if (n<=maxnb) dat=unpack(climate%evap_PT(pind(n,1):pind(n,2)),tmap(:,n),dat)
    !write(vname,'("t",I1.1,"_climateevapPT")') n
    !call histwrt3(dat,vname,idnc,iarch,local,.true.)
    !dat=0.
    !if (n<=maxnb) dat=unpack(climate%aevap(pind(n,1):pind(n,2)),tmap(:,n),dat)
    !write(vname,'("t",I1.1,"_climateaevap")') n
    !call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%mtemp_max(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climatemtempmax")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%mtemp_min(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climatemtempmin")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%qtemp_max(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climateqtempmax")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%qtemp_max_last_year(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%gdd0(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climategdd0")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%gdd5(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climategdd5")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%agdd0(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climateagdd0")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(climate%agdd5(pind(n,1):pind(n,2)),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climateagdd5")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(climate%chilldays(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("t",I1.1,"_climatechilldays")') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    if ( n==1 ) then
      do k = 1,91
        dat=0.
        if (n<=maxnb) dat=unpack(climate%dtemp_91(pind(n,1):pind(n,2),k),tmap(:,n),dat)
        write(vname,'("t",I1.1,"_climatedtemp",I2.2)') n,k
        call histwrt3(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,31
        if (n<=maxnb) dat=unpack(climate%dmoist_31(pind(n,1):pind(n,2),k),tmap(:,n),dat)
        write(vname,'("t",I1.1,"_climatedmoist",I2.2)') n,k
        call histwrt3(dat,vname,idnc,iarch,local,.true.)
      end do
    end if
  end if
  !dat=0.
  !if (n<=maxnb) dat=unpack(real(canopy%fhs_cor(pind(n,1):pind(n,2))),tmap(:,n),dat)
  !write(vname,'("t",I1.1,"_fhscor")') n
  !call histwrt3(dat,vname,idnc,iarch,local,.true.)
  !dat=0.
  !if (n<=maxnb) dat=unpack(real(canopy%fes_cor(pind(n,1):pind(n,2))),tmap(:,n),dat)
  !write(vname,'("t",I1.1,"_fescor")') n
  !call histwrt3(dat,vname,idnc,iarch,local,.true.)
end do
vname='albvisdir'
call histwrt3(albvisdir,vname,idnc,iarch,local,.true.)
vname='albvisdif'
call histwrt3(albvisdif,vname,idnc,iarch,local,.true.)
vname='albnirdir'
call histwrt3(albnirdir,vname,idnc,iarch,local,.true.)
vname='albnirdif'
call histwrt3(albnirdif,vname,idnc,iarch,local,.true.)
vname='albvis'
call histwrt3(albvisnir(:,1),vname,idnc,iarch,local,.true.)
vname='albnir'
call histwrt3(albvisnir(:,2),vname,idnc,iarch,local,.true.)
  
return
end subroutine savetile 

! *************************************************************************************
! Water inflow from river routing
subroutine cableinflow(inflow)

use newmpar_m
use soil_m

implicit none

integer nb, k
real, dimension(ifull), intent(inout) :: inflow
real, dimension(ifull) :: delflow
real, dimension(mp) :: xx, ll, delxx

if ( mp<=0 ) return

do nb = 1,maxnb
  xx(pind(nb,1):pind(nb,2)) = pack( inflow(1:ifull), tmap(:,nb) )
end do
delxx(1:mp) = 0.
do k = 1,cbm_ms
  ll(1:mp) = max( soil%sfc(1:mp)-real(ssnow%wb(1:mp,k)), 0. )*1000.*soil%zse(k)
  ll(1:mp) = min( xx(1:mp), ll(1:mp) )
  ssnow%wb(1:mp,k) = ssnow%wb(1:mp,k) + ll(1:mp)/(1000.*soil%zse(k))
  delxx(1:mp) = delxx(1:mp) - ll(1:mp)
end do
delflow(1:ifull) = 0.
do nb = 1,maxnb
  delflow(1:ifull) = delflow(1:ifull) + unpack(sv(pind(nb,1):pind(nb,2))*delxx(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
end do
inflow(1:ifull) = inflow(1:ifull) + delflow(1:ifull)

return
end subroutine cableinflow

! *************************************************************************************
! Transfer grid information from CABLE internally, read N&P input from
! integral NETCDF file
subroutine casa_readpoint(casafile,casapoint)

use cc_mpi
use infile
use newmpar_m
use parmgeom_m

implicit none

integer ncstatus, ncid, varid, tilg
integer, dimension(2) :: spos, npos
real tlat, tlon, tschmidt
real, dimension(:,:), allocatable :: dumg
real, dimension(ifull,5), intent(out) :: casapoint
character(len=*), intent(in) :: casafile
logical tst

if ( myid==0 ) then
  allocate( dumg(ifull_g,5) )
  write(6,*) "Reading ",trim(casafile)
  call ccnf_open(casafile,ncid,ncstatus)
  call ncmsg('CASA_readpoint',ncstatus)
  ! check dimensions and location
  call ccnf_get_attg(ncid,'lat0',tlat)
  call ccnf_get_attg(ncid,'lon0',tlon)
  call ccnf_get_attg(ncid,'schmidt0',tschmidt)
  if ( abs(rlong0-tlon)>1.e-20 .or. abs(rlat0-tlat)>1.e-20 .or. abs(schmidt-tschmidt)>1.e-20 ) then
    write(6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    call ccmpi_abort(-1)
  end if
  call ccnf_inq_dimlen(ncid,'longitude',tilg)
  if ( tilg /= il_g ) then
    write (6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write (6,*) "il_g,tilg ",il_g,tilg
    call ccmpi_abort(-1)
  end if
  ! load casa fields
  spos(1:2) = 1
  npos(1) = il_g
  npos(2) = il_g*6
  write(6,*) "Loading soil order"
  call ccnf_inq_varid(ncid,'sorder',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  write(6,*) "Loading N deposition rate"
  call ccnf_inq_varid(ncid,'ndep',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  write(6,*) "Loading N fixation rate"
  call ccnf_inq_varid(ncid,'nfix',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  write(6,*) "Loading P dust deposition"
  call ccnf_inq_varid(ncid,'pdust',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,4))
  write(6,*) "Loading P weathering rate"
  call ccnf_inq_varid(ncid,'pweather',varid,tst)
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,5))
  call ccnf_close(ncid)
  call ccmpi_distribute(casapoint,dumg)
  deallocate(dumg)
else
  call ccmpi_distribute(casapoint)
end if

end subroutine casa_readpoint

! *************************************************************************************  
! This subroutine reads the MODIS derived leaf phenology data
subroutine casa_readphen(fphen,greenup,fall,phendoy1)

use cc_mpi
use newmpar_m

implicit none

integer, parameter :: nphen = 8 ! was 10(IGBP). changed by Q.Zhang @01/12/2011
integer ilat, ierr
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
integer, dimension(nphen) :: greenupx, fallx, xphendoy1
integer, dimension(nphen) :: ivtx
real :: xlat
character(len=*), intent(in) :: fphen

! initilize for evergreen PFTs
greenup(:,:)  = -50
fall(:,:)     = 367
phendoy1(:,:) = 2

! MJT notes - if the CSIRO PFTs are redefined, then this phenology data will be mismatched

if ( myid==0 ) then
  write(6,*) "Reading CASA leaf phenology data"
  open(87,file=fphen,status='old',iostat=ierr)
  if ( ierr/=0 ) then
    write(6,*) "ERROR: Cannot open phenfile=",trim(fphen)
    call ccmpi_abort(-1)
  end if
  read(87,*)
  read(87,*) ivtx
  do ilat = 271,1,-1
    read(87,*) xlat,greenupx,fallx,xphendoy1 
    greenup(ilat,ivtx(:))  = greenupx(:)
    fall(ilat,ivtx(:))     = fallx(:)
    phendoy1(ilat,ivtx(:)) = xphendoy1(:)
  end do
  close(87)
end if
call ccmpi_bcast(greenup, 0,comm_world)
call ccmpi_bcast(fall,    0,comm_world)
call ccmpi_bcast(phendoy1,0,comm_world)

return
end subroutine casa_readphen

! *************************************************************************************  
! This overrides the snow and surface temperatures for the SCM mode
subroutine cablesettemp(tset)

use newmpar_m

implicit none

integer k
real, dimension(ifull), intent(in) :: tset

ssnow%tss = tset(:)
do k = 1,3
  ssnow%tggsn(:,k) = tset(:)
end do
do k = 1,ms
  ssnow%tgg(:,k) = tset(:)
end do

return
end subroutine cablesettemp

end module cable_ccam

