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

! CABLE interface originally developed by the CABLE group
! Subsequently modified by MJT for maxtile=5 tile mosaic and SEAESF radiation scheme

! Thanks to Paul Ryan for OMP routines
    
! - Currently all tiles have the same soil texture, but independent soil temperatures,
!   moisture, etc.
! - LAI can be interpolated between timesteps using a PW-linear fit to the LAI integral
!   or LAI can be taken as constant for the month.  CASA-CNP can predict LAI and vcmax,
!   although this can take considerable time to spin-up.
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.
! - Options exist for using SLI soil model (soil_struc=1), climate feedbacks
!   (cable_climate=1) and the POP model (cable_pop=1)

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
! (18 Evergreen Broadleaf Savanna) - MJT defined
  
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
   
    
! The following mappings between IGBP and CSIRO PFT were recommended by Rachel Law
    
! ivegt   IGBP type                             CSIRO PFT
! 1       Evergreen Needleleaf Forest           1.  Evergreen Needleleaf
! 2       Evergreen Broadleaf Forest            1.  Evergreen Broadleaf
! 3       Deciduous Needleleaf Forest           1.  Deciduous Needleleaf
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

module cable_ccam

use cable_air_module
use cable_albedo_module
use cable_canopy_module
use cable_ccam3
use cable_ccam4
use cable_common_module
use cable_data_module
use cable_def_types_mod, cbm_ms => ms
use cable_optimise_JV_module, only : optimise_JV
use cable_radiation_module
use cable_roughness_module
use cable_soil_snow_module
use casa_cnp_module
use casadimension
use casaparm, xroot => froot
use casavariable
use phenvariable
use popmodule, only : pop_init, popstep
use pop_types
use sli_main_mod

implicit none

private
public sib4
public loadcbmparm, cbmparm, loadtile, defaulttile, savetiledef, savetile, newcbmwb
public cablesettemp, cableinflow, cbmemiss
public proglai, progvcmax, maxtile, soil_struc, cable_pop, ccycle
public fwsoil_switch, cable_litter, gs_switch, cable_climate
public smrf_switch, strf_switch, cable_gw_model
public POP_NPATCH, POP_NCOHORT, POP_AGEMAX

! CABLE biophysical options
integer, save :: soil_struc      = 0          ! 0 default, 1 SLI soil model
integer, save :: fwsoil_switch   = 0          ! 0 default, 1 non-linear, 2 Lai and Ktaul, 3 Haverd2013
integer, save :: cable_litter    = 0          ! 0 off, 1 on
integer, save :: gs_switch       = 0          ! 0 leuning, 1 medlyn
integer, save :: smrf_switch     = 4          ! 1 CASA-CNP, 2 SOLIN, 3 TRIFFID, 4 Trudinger2016(default),
                                              ! 5 DAMM (Soil Moist Respiration Function)
integer, save :: strf_switch     = 4          ! 1 CASA-CNP, 2 K1995, 3 PnET-CN, 4 LT1994(default),
                                              ! 5 DAMM (Soil Temp Respiration Function)
integer, save :: cable_gw_model  = 0          ! 0 off, 1 GW_Hydro
! CABLE biochemical options
integer, save :: ccycle          = 0          ! 0 off, 1 (C), 2 (CN), 3 (CNP)
integer, save :: proglai         = -1         ! -1, piece-wise linear prescribed LAI, 0 PWCB prescribed LAI, 1 prognostic LAI
integer, save :: progvcmax       = 0          ! 0 prescribed, 1 prognostic vcmax (standard), 2 prognostic vcmax (Walker2014)
! CABLE POP options
integer, save :: cable_pop       = 0          ! 0 off, 1 on
! CABLE climate options
integer, save :: cable_climate   = 0          ! 0 off, 1 on
! CABLE parameters
integer, parameter :: maxtile    = 9          ! maximum possible number of tiles (1-5=natural/secondary, 6-7=pasture/rangeland, 8-9=crops)
integer, parameter :: COLDEST_DAY_NHEMISPHERE = 355
integer, parameter :: COLDEST_DAY_SHEMISPHERE = 172
integer, save :: POP_NPATCH      = -1
integer, save :: POP_NLAYER      = -1
integer, save :: POP_NCOHORT     = -1
integer, save :: POP_HEIGHT_BINS = -1
integer, save :: POP_NDISTURB    = -1
integer, save :: POP_AGEMAX      = -1

integer, save :: maxnb                           ! maximum number of actual tiles
integer, save :: mp_global                       ! maximum number of land-points on this process
real, dimension(:), allocatable, target, save :: sv, vl1, vl2, vl3, vl4
integer, dimension(:,:), allocatable :: qmap

contains
! ****************************************************************************

! CABLE-CCAM interface
subroutine sib4

use cc_omp, only : ntiles, imax
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
use raddiag_m
use radisw_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use tracers_m
use vegpar_m
use work2_m, only : qsttg,zo,zoh,zoq,theta,vmod,wetfac
use work3_m, only : ga
use zenith_m

implicit none

integer, dimension(maxtile,2) :: ltind
integer, dimension(imax) :: lclimate_biome, lclimate_iveg, lclimate_gmd
integer :: lmaxnb, tile, is, ie, js, je
integer :: ico2, igas
real, dimension(imax,mlitter) :: lclitter, lnilitter, lplitter
real, dimension(imax,mplant) :: lcplant, lniplant, lpplant
real, dimension(imax,msoil) :: lcsoil, lnisoil, lpsoil
real, dimension(imax,ms) :: ltgg, lwb, lwbice
real, dimension(imax,3) :: lsmass, lssdn, ltggsn
real, dimension(imax) :: lcnbp, lcnpp, lfnee, lfpn, lfrd, lfrp, lfrpr, lfrpw, lfrs
real, dimension(imax) :: latmco2
real, dimension(imax) :: lclimate_min20, lclimate_max20, lclimate_alpha20
real, dimension(imax) :: lclimate_agdd5
real, dimension(imax) :: lclimate_dmoist_min20, lclimate_dmoist_max20
real, dimension(imax) :: lfevc,lplant_turnover,lplant_turnover_wood
logical, dimension(imax,maxtile) :: ltmap
type(air_type) :: lair
type(balances_type) :: lbal
type(canopy_type) :: lcanopy
type(casa_balance) :: lcasabal
type(casa_flux) :: lcasaflux
type(casa_met) :: lcasamet
type(casa_pool) :: lcasapool
type(climate_type) :: lclimate
type(met_type) :: lmet
type(phen_variable) :: lphen
type(pop_type) :: lpop
type(radiation_type) :: lrad
type(roughness_type) :: lrough
type(soil_parameter_type) :: lsoil
type(soil_snow_type) :: lssnow
type(sum_flux_type) :: lsum_flux
type(veg_parameter_type) :: lveg
type(climate_save_type) :: lclimate_save

!$omp do schedule(static) private(is,ie,js,je,ltind,ltmap,lmaxnb),                                             &
!$omp private(lclitter,lcnbp,lcnpp,lcplant,lcsoil,lfnee,lfpn),                                                 &
!$omp private(lfrd,lfrp,lfrpr,lfrpw,lfrs,lnilitter,lniplant,lnisoil,lplitter,lpplant,lpsoil),                  &
!$omp private(lsmass,lssdn,ltgg,ltggsn,lwb,lwbice,lair,lbal,lcanopy,lcasabal,latmco2),                         &
!$omp private(lcasaflux,lcasamet,lcasapool,lclimate,lmet,lphen,lpop,lrad,lrough,lsoil,lssnow,lsum_flux,lveg),  &
!$omp private(lclimate_biome,lclimate_iveg,lclimate_min20,lclimate_max20,lclimate_alpha20,lclimate_agdd5),     &
!$omp private(lclimate_gmd,lclimate_dmoist_min20,lclimate_dmoist_max20,lfevc,lplant_turnover),                 &
!$omp private(lplant_turnover_wood,lclimate_save)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  mp = tdata(tile)%mp
  lmaxnb = tdata(tile)%maxnb
  
  js = tdata(tile)%toffset + 1
  je = tdata(tile)%toffset + tdata(tile)%mp

  ltind = tdata(tile)%tind - tdata(tile)%toffset

  if ( mp>0 ) then
    ltmap = tdata(tile)%tmap       
    lsmass = smass(is:ie,:)
    lssdn = ssdn(is:ie,:)
    ltgg = tgg(is:ie,:)
    ltggsn = tggsn(is:ie,:)
    lwb = wb(is:ie,:)
    lwbice = wbice(is:ie,:)
    if ( ccycle/=0 ) then
      lclitter = clitter(is:ie,:)
      lcnbp = cnbp(is:ie)
      lcnpp = cnpp(is:ie)
      lcplant = cplant(is:ie,:)
      lcsoil = csoil(is:ie,:)
      lfnee = fnee(is:ie)
      lfpn = fpn(is:ie)
      lfrd = frd(is:ie)
      lfrp = frp(is:ie)
      lfrpr = frpr(is:ie)
      lfrpw = frpw(is:ie)
      lfrs = frs(is:ie)
      lnilitter = nilitter(is:ie,:)
      lniplant = niplant(is:ie,:)
      lnisoil = nisoil(is:ie,:)
      lplitter = plitter(is:ie,:)
      lpplant = pplant(is:ie,:)
      lpsoil = psoil(is:ie,:)
      if ( diaglevel_carbon>0 ) then
        lfevc = fevc(is:ie)
        lplant_turnover = plant_turnover(is:ie)
        lplant_turnover_wood = plant_turnover_wood(is:ie)
      end if
    end if 
    if ( cable_climate==1 ) then
      lclimate_biome = 999
      lclimate_iveg = 0
      lclimate_min20 = 0.
      lclimate_max20 = 0.
      lclimate_alpha20 = 0.
      lclimate_agdd5 = 0.
      lclimate_gmd = 0
      lclimate_dmoist_min20 = 0.
      lclimate_dmoist_max20 = 0.
    end if
    
    ! set co2 forcing for cable
    ! constant: atmospheric co2 = 360 ppm 
    ! host: atmospheric co2 follows that from CCAM radiation scheme
    ! interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)
    latmco2 = 1.E6*rrvco2          ! from radiative CO2 forcings
    ico2 = 0
    do igas = 1,ngas
      if ( trim(tractype(igas))=='online' .and. trim(tracname(igas))=='cbmnep' ) then
        ico2 = igas
        exit
      end if
    end do
    if ( ico2>0 ) then
      latmco2 = tr(1:imax,1,ico2) ! use interactive tracers
    end if
    
    !set pointers to pass through
    call setp(air,lair,tile)
    call setp(bal,lbal,tile)
    call setp(canopy,lcanopy,tile)
    call setp(met,lmet,tile)
    call setp(rad,lrad,tile)
    call setp(rough,lrough,tile)
    call setp(soil,lsoil,tile)
    call setp(ssnow,lssnow,tile)
    call setp(sum_flux,lsum_flux,tile)
    call setp(veg,lveg,tile)
    if ( ccycle/=0 ) then
      call setp(casabal,lcasabal,tile)
      call setp(casaflux,lcasaflux,tile)
      call setp(casamet,lcasamet,tile)
      call setp(casapool,lcasapool,tile)
      call setp(phen,lphen,tile)
      if ( cable_pop==1 ) then
        call setp(pop,lpop,tile)  
      end if    
    end if
    if ( cable_climate==1 ) then
      call setp(climate,lclimate,tile)
      call setp(climate_save,lclimate_save,tile)
    end if  

    call sib4_work(albnirdif(is:ie),albnirdir(is:ie),albnirsav(is:ie),albvisdif(is:ie),albvisdir(is:ie),                 &
                   albvissav(is:ie),cansto(is:ie),cdtq(is:ie),cduv(is:ie),lclitter,lcnbp,                                &
                   lcnpp,condg(is:ie),conds(is:ie),condx(is:ie),lcplant,lcsoil,eg(is:ie),epot(is:ie),fbeamnir(is:ie),    &
                   fbeamvis(is:ie),fg(is:ie),lfnee,lfpn,lfrd,lfrp,                                                       &
                   lfrpr,lfrpw,lfrs,fwet(is:ie),ga(is:ie),isflag(is:ie),land(is:ie),lmaxnb,mp,lnilitter,lniplant,        &
                   lnisoil,lplitter,lpplant,ps(is:ie),lpsoil,qg(is:ie,1),qsttg(is:ie),rgsave(is:ie),                     &
                   rlatt(is:ie),rlongg(is:ie),rnet(is:ie),rsmin(is:ie),runoff(is:ie),runoff_surface(is:ie),              &
                   sigmf(is:ie),lsmass,snage(is:ie),snowd(is:ie),snowmelt(is:ie),lssdn,ssdnn(is:ie),                     &
                   sv(js:je),sgdn(is:ie),swrsave(is:ie),t(is:ie,1),ltgg,ltggsn,theta(is:ie),ltind,ltmap,                 &
                   latmco2,tss(is:ie),ustar(is:ie),vlai(is:ie),vl1(js:je),vl2(js:je),vl3(js:je),vl4(js:je),vmod(is:ie),  &
                   lwb,lwbice,wetfac(is:ie),zo(is:ie),zoh(is:ie),zoq(is:ie),evspsbl(is:ie),sbl(is:ie),lair,lbal,c,       &
                   lcanopy,lcasabal,casabiome,lcasaflux,lcasamet,lcasapool,lclimate,lmet,lphen,lpop,lrad,lrough,lsoil,   &
                   lssnow,lsum_flux,lveg,lclimate_iveg,lclimate_biome,lclimate_min20,lclimate_max20,lclimate_alpha20,    &
                   lclimate_agdd5,lclimate_gmd,lclimate_dmoist_min20,lclimate_dmoist_max20,lfevc,lplant_turnover,        &
                   lplant_turnover_wood,lclimate_save,imax)

    smass(is:ie,:) = lsmass
    ssdn(is:ie,:) = lssdn
    tgg(is:ie,:) = ltgg
    tggsn(is:ie,:) = ltggsn
    wb(is:ie,:) = lwb
    wbice(is:ie,:) = lwbice
    if ( ccycle/=0 ) then
      clitter(is:ie,:) = lclitter
      cnbp(is:ie) = lcnbp
      cnpp(is:ie) = lcnpp
      cplant(is:ie,:) = lcplant
      csoil(is:ie,:) = lcsoil
      fnee(is:ie) = lfnee
      fpn(is:ie) = lfpn
      frd(is:ie) = lfrd
      frp(is:ie) = lfrp
      frpr(is:ie) = lfrpr
      frpw(is:ie) = lfrpw
      frs(is:ie) = lfrs
      nilitter(is:ie,:) = lnilitter
      niplant(is:ie,:) = lniplant
      nisoil(is:ie,:) = lnisoil
      plitter(is:ie,:) = lplitter
      pplant(is:ie,:) = lpplant
      psoil(is:ie,:) = lpsoil
      if ( diaglevel_carbon > 0 ) then
        fevc(is:ie) = lfevc
        plant_turnover(is:ie) = lplant_turnover
        plant_turnover_wood(is:ie) = lplant_turnover_wood
      end if
    end if  
    if ( cable_climate==1 ) then
      climate_biome(is:ie) = lclimate_biome
      climate_ivegt(is:ie) = lclimate_iveg
      climate_min20(is:ie) = lclimate_min20
      climate_max20(is:ie) = lclimate_max20
      climate_alpha20(is:ie) = lclimate_alpha20
      climate_agdd5(is:ie) = lclimate_agdd5
      climate_gmd(is:ie) = lclimate_gmd
      climate_dmoist_min20(is:ie) = lclimate_dmoist_min20
      climate_dmoist_max20(is:ie) = lclimate_dmoist_max20
    end if
  end if

end do
!$omp end do nowait

return
end subroutine sib4

! CABLE-CCAM interface
subroutine sib4_work(albnirdif,albnirdir,albnirsav,albvisdif,albvisdir,albvissav,cansto,cdtq,cduv,clitter,cnbp, &
                     cnpp,condg,conds,condx,cplant,csoil,eg,epot,fbeamnir,fbeamvis,fg,fnee,fpn,frd,frp,frpr,    &
                     frpw,frs,fwet,ga,isflag,land,maxnb,mp,nilitter,niplant,nisoil,plitter,pplant,ps,            &
                     psoil,qg,qsttg,rgsave,rlatt,rlongg,rnet,rsmin,runoff,runoff_surface,sigmf,smass,            &
                     snage,snowd,snowmelt,ssdn,ssdnn,sv,sgdn,swrsave,t,tgg,tggsn,theta,tind,tmap,atmco2,tss,    &
                     ustar,vlai,vl1,vl2,vl3,vl4,vmod,wb,wbice,wetfac,zo,zoh,zoq,evspsbl,sbl,air,bal,c,canopy,   &
                     casabal,casabiome,casaflux,casamet,casapool,climate,met,phen,pop,rad,rough,soil,ssnow,     &
                     sum_flux,veg,climate_ivegt,climate_biome,climate_min20,climate_max20,climate_alpha20,      &
                     climate_agdd5,climate_gmd,climate_dmoist_min20,climate_dmoist_max20,fevc,plant_turnover,   &
                     plant_turnover_wood,climate_save,imax)

use const_phys
use dates_m
use estab, only : qsat
use infile, only : getzinp
use newmpar_m
use parm_m
use sigs_m
use soil_m, only : zmin
use zenith_m, only : solargh, zenith
  
implicit none

integer, intent(in) :: imax, maxnb, mp
integer, dimension(imax), intent(inout) :: climate_ivegt, climate_biome, climate_gmd
integer jyear, jmonth, jday, jhour, jmin, idoy
integer k, mins, nb, j
integer is, ie, casaperiod, npercasa
integer lalloc
integer mp_POP
integer, dimension(maxtile,2), intent(in) :: tind
integer, dimension(imax), intent(inout) :: isflag
integer, dimension(13), parameter :: ndoy=(/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 0 /) 
real fjd, r1, dlt, slag, dhr, alp
real, dimension(imax), intent(in) :: atmco2
real, dimension(imax), intent(in) :: qg, t
real, dimension(imax,mlitter), intent(inout) :: clitter, nilitter, plitter
real, dimension(imax,mplant), intent(inout) :: cplant, niplant, pplant
real, dimension(imax,msoil), intent(inout) :: csoil, nisoil, psoil
real, dimension(imax,ms), intent(inout) :: tgg, wb, wbice
real, dimension(imax,3), intent(inout) :: smass, ssdn, tggsn
real, dimension(imax), intent(inout) :: albnirdif, albnirdir, albnirsav, albvisdif, albvisdir, albvissav
real, dimension(imax), intent(inout) :: cansto, cdtq, cduv, cnbp, cnpp, eg, epot, fg, fnee, fpn, frd
real, dimension(imax), intent(inout) :: frp, frpr, frpw, frs, fwet, ga, qsttg, rnet, rsmin, runoff
real, dimension(imax), intent(inout) :: runoff_surface, sigmf, snage, snowd, snowmelt, ssdnn, tss, ustar
real, dimension(imax), intent(inout) :: vlai, wetfac, zo, zoh, zoq, evspsbl, sbl
real, dimension(imax), intent(inout) :: climate_min20, climate_max20, climate_alpha20
real, dimension(imax), intent(inout) :: climate_agdd5
real, dimension(imax), intent(inout) :: climate_dmoist_min20, climate_dmoist_max20
real, dimension(imax), intent(inout) :: fevc, plant_turnover, plant_turnover_wood
real, dimension(imax), intent(in) :: condg, conds, condx, fbeamnir, fbeamvis, ps, rgsave, rlatt, rlongg
real, dimension(imax), intent(in) :: sgdn, swrsave, theta, vmod
real, dimension(imax) :: coszro2, taudar2, tmps, qsttg_land
real, dimension(imax) :: evspsbl_l, sbl_l
real, dimension(mp), intent(in) :: sv, vl1, vl2, vl3, vl4
real(kind=8) :: dtr8
real(kind=8), dimension(mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
real(kind=8), dimension(mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
real(kind=8), dimension(mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
real(r_2), dimension(mp) :: xKNlimiting, xkleafcold, xkleafdry
real(r_2), dimension(mp) :: xkleaf, xnplimit, xNPuptake, xklitter
real(r_2), dimension(mp) :: xksoil
logical, dimension(imax,maxtile), intent(in) :: tmap
logical, dimension(imax), intent(in) :: land
type(air_type), intent(inout) :: air
type(balances_type), intent(inout) :: bal
type(physical_constants), intent(in) :: c
type(canopy_type), intent(inout) :: canopy
type(casa_balance), intent(inout) :: casabal
type(casa_biome), intent(inout) :: casabiome
type(casa_flux), intent(inout) :: casaflux
type(casa_met), intent(inout) :: casamet
type(casa_pool), intent(inout) :: casapool
type(climate_type), intent(inout) :: climate
type(met_type), intent(inout) :: met
type(phen_variable), intent(inout) :: phen
type(pop_type), intent(inout) :: pop
type(radiation_type), intent(inout) :: rad
type(roughness_type), intent(inout) :: rough
type(soil_parameter_type), intent(inout) :: soil
type(soil_snow_type), intent(inout) :: ssnow
type(sum_flux_type), intent(inout) :: sum_flux
type(veg_parameter_type), intent(inout) :: veg
type(climate_save_type), intent(inout) :: climate_save

cansto = 0.
fwet = 0.
vlai = 0.

! calculate time from beginning of the simulation
! jyear, jmonth, jday, jhour, jmin are for the start date of the simulation
! mins is the time from the start of the year
! mtimer is the time elapsed from the start of the simulation
dtr8 = real(dt, 8)
dhr = dt/3600.
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
! calculate zenith angle
fjd = float(mod(mins, 525600))/1440. ! restrict to 365 day calendar
idoy = int(fjd)
call solargh(fjd,bpyear,r1,dlt,alp,slag)
call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,imax,coszro2,taudar2)

!! store soil temperature and moisture for budget calculation
!bwb(:,:) = ssnow%wb(:,:)
!btgg(:,:) = ssnow%tgg(:,:)

! run CASA CNP once per day
casaperiod = nint(86400._8*deltpool)
npercasa = max( nint(real(casaperiod,8)/dtr8), 1 )

! checks
if ( any(atmco2<1.) ) then
  write(6,*) "ERROR: Invalid CO2 mixing ratio in cable_ccam2 ",minval(atmco2)
  stop
end if

! set meteorological forcing
albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif ! for nrad=4
albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif ! for nrad=4
do nb = 1,maxnb
  is = tind(nb,1)
  ie = tind(nb,2)
  met%tk(is:ie)        = real(pack(theta,      tmap(:,nb)), 8)
  met%ua(is:ie)        = real(pack(vmod,       tmap(:,nb)), 8)
  met%ca(is:ie)        = real(pack(atmco2,     tmap(:,nb))*1.e-6, 8)
  met%coszen(is:ie)    = real(pack(coszro2,    tmap(:,nb)), 8)        ! use instantaneous value
  met%qv(is:ie)        = real(pack(qg(1:imax), tmap(:,nb)),8 )        ! specific humidity in kg/kg
  met%pmb(is:ie)       = real(pack(ps(1:imax), tmap(:,nb))*0.01, 8)   ! pressure in mb at ref height
  met%precip(is:ie)    = real(pack(condx,      tmap(:,nb)), 8)        ! in mm not mm/sec
  met%precip_sn(is:ie) = real(pack(conds+condg,tmap(:,nb)), 8)        ! in mm not mm/sec
  ! swrsave indicates the fraction of net VIS radiation (compared to NIR)
  ! fbeamvis indicates the beam fraction of downwelling direct radiation (compared to diffuse) for VIS
  ! fbeamnir indicates the beam fraction of downwelling direct radiation (compared to diffuse) for NIR
  met%fsd(is:ie,1)     = real(pack(swrsave*sgdn,     tmap(:,nb)), 8)
  met%fsd(is:ie,2)     = real(pack((1.-swrsave)*sgdn,tmap(:,nb)), 8)
  met%fld(is:ie)       = real(pack(-rgsave,          tmap(:,nb)), 8)      ! long wave down (positive) W/m^2
  rad%fbeam(is:ie,1)   = real(pack(fbeamvis,         tmap(:,nb)), 8)
  rad%fbeam(is:ie,2)   = real(pack(fbeamnir,         tmap(:,nb)), 8)
  rough%za_tq(is:ie)   = real(pack(bet(1)*t(1:imax), tmap(:,nb))/grav, 8) ! reference height
  rough%za_tq(is:ie)   = max( rough%za_tq(is:ie), veg%hc(is:ie)+1. )
end do
met%doy        = real(fjd, 8)
met%hod        = (met%doy-int(met%doy))*24._8 + rad%longitude/15._8
met%hod        = mod(met%hod, 24._8)
met%tvair      = met%tk
met%tvrad      = met%tk
met%qvair      = met%qv
met%ua         = max(met%ua, c%umin)
met%coszen     = max(met%coszen, 1.e-8_8) 
rough%za_uv    = rough%za_tq
rad%fbeam(:,3) = 0._8            ! dummy for now

! Interpolate LAI.  Also need sigmf for LDR prognostic aerosols.
call setlai(sigmf,jmonth,jday,jhour,jmin,mp,sv,vl1,vl2,vl3,vl4,casamet,veg,imax,tind,tmap,maxnb)

! Calculate vcmax
call vcmax_feedback(casabiome,casamet,casapool,veg,climate,ktau)

!--------------------------------------------------------------
! CABLE
canopy%fev       = 0._8
canopy%fes       = 0._8
canopy%fhv       = 0._8
canopy%fhs       = 0._8
met%ofsd         = met%fsd(:,1) + met%fsd(:,2)
ssnow%owetfac    = ssnow%wetfac
canopy%oldcansto = canopy%cansto
call ruff_resist(veg,rough,ssnow,canopy)
call define_air(met,air)
call init_radiation(met,rad,veg,canopy)
call surface_albedo(ssnow,veg,met,rad,soil,canopy)
call define_canopy(bal,rad,rough,air,met,dtr8,ssnow,soil,veg,canopy,climate)
ssnow%otss_0  = ssnow%otss
ssnow%otss    = ssnow%tss
ssnow%owetfac = ssnow%wetfac
select case ( soil_struc )
  case(0)  
    call soil_snow(dtr8,soil,ssnow,canopy,met,bal,veg)
  case(1)
    call sli_main(999,dtr8,veg,soil,ssnow,met,canopy,air,rad,0)   
  case default
    write(6,*) "ERROR: Unknown option soil_struc ",soil_struc
    stop
end select
! adjust for new soil temperature
ssnow%deltss     = ssnow%tss - ssnow%otss
if ( soil_struc==0 ) then
  canopy%fhs       = canopy%fhs + ssnow%deltss*ssnow%dfh_dtg
  canopy%fes       = canopy%fes + ssnow%deltss*ssnow%dfe_ddq*ssnow%ddq_dtg
  canopy%fns       = canopy%fns + ssnow%deltss*ssnow%dfn_dtg
  canopy%ga        = canopy%ga  + ssnow%deltss*canopy%dgdtg
  !canopy%fhs_cor  = canopy%fhs_cor + ssnow%deltss*ssnow%dfh_dtg
  !canopy%fes_cor  = canopy%fes_cor + ssnow%deltss*ssnow%dfe_ddq*ssnow%ddq_dtg
  !canopy%fns_cor  = canopy%fns_cor + ssnow%deltss*ssnow%dfn_dtg
  !canopy%ga_cor   = canopy%ga_cor + ssnow%deltss*canopy%dgdtg
  ! MJT fix
  ssnow%wb(:,ms) = ssnow%wb(:,ms) - (ssnow%deltss*ssnow%dfe_ddq*ssnow%ddq_dtg)*dtr8 &
                                   /(ssnow%cls*air%rlam*soil%zse(ms)*1000._8) 
end if
canopy%fh        = canopy%fhv + canopy%fhs
canopy%fev       = canopy%fevc + canopy%fevw
canopy%fe        = canopy%fev + canopy%fes
canopy%rnet      = canopy%fns + canopy%fnv
rad%trad         = ( (1._8-rad%transd)*canopy%tv**4 + rad%transd*ssnow%tss**4 )**0.25_8

! note that conservation is still preserved at this point
! canopy%ga    = canopy%rnet - canopy%fh - canopy%fe
! canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq*ssnow%ddq_dtg

! EK suggestion
!canopy%cdtq =  max( 0.1_8*canopy%cduv, canopy%cdtq )
! MJT suggestion
canopy%cdtq =  max( 0._8, canopy%cdtq )

! MJT suggestion
ssnow%wbice = max( ssnow%wbice, 0._8 )

! change in soil temperature
!deltgg(:) = 0._8
!do k = 1,ms
!  deltgg(:) = deltgg(:) + (ssnow%tgg(:,k)-btgg(:,k))*ssnow%gammzz(:,k)/dtr8
!end do  
!! energy balance
!tbal = canopy%rnet - canopy%fh - canopy%fe - deltgg
!write(6,*) "energy ",minval(tbal),maxval(tbal)
!! change in soil moisture
!delwb(:) = 0._8
!do k = 1,ms
!  delwb(:) = delwb(:) + (ssnow%wb(:,k)-bwb(:,k))*soil%zse(k)*1000.
!end do  
!! net water into soil
!bal%wbal = met%precip - canopy%delwc - ssnow%snowd + ssnow%osnowd - ssnow%runoff    &
!          - (canopy%fevw + canopy%fevc + canopy%fes/ssnow%cls)*dtr8/air%rlam - delwb
!write(6,*) "bal%wbal ",minval(bal%wbal),maxval(bal%wbal)

#ifdef cabledebug
if ( any( canopy%fhv/=canopy%fhv ) ) then
  write(6,*) "ERROR: NaN found in canopy%fhv after CABLE"
  stop -1
end if
if ( any( canopy%fhs/=canopy%fhs ) ) then
  write(6,*) "ERROR: NaN found in canopy%fhs after CABLE"
  stop -1
end if
if ( any( canopy%fev/=canopy%fev ) ) then
  write(6,*) "ERROR: NaN found in canopy%fev after CABLE"
  stop -1
end if
if ( any( canopy%fes/=canopy%fes ) ) then
  write(6,*) "ERROR: NaN found in canopy%fes after CABLE"
  stop -1
end if
if ( any( canopy%tv/=canopy%tv ) ) then
  write(6,*) "ERROR: NaN found in canopy%tv after CABLE"
  stop -1
end if
if ( any( ssnow%tss/=ssnow%tss ) ) then
  write(6,*) "ERROR: NaN found in ssnow%tss after CABLE"
  stop -1
end if
#endif

!--------------------------------------------------------------
! CABLE CLIMATE
call cableclimate(idoy,jmonth,ndoy,canopy,climate,met,rad,ssnow,veg, &
                  climate_save,npercasa,ktau)

!--------------------------------------------------------------
! CASA CNP
select case (ccycle)
  case(0) ! off

  case(1,2,3) ! C, C+N, C+N+P
    ! update casamet
    if ( cable_pop==1 ) then
      lalloc = 3
    else
      lalloc = 0  
    end if
    !if ( abs(mod(real(casaperiod),dt))>1.e-20 ) then
    !  write(6,*) "ERROR: Invalid casaperiod ",casaperiod
    !  write(6,*) "Period must be an integer multiple of the time-step dt ",dt
    !  write(6,*) "Please adjust deltpool from CABLE CASA-CNP accordingly"
    !  stop
    !end if
    casamet%tairk = casamet%tairk + met%tk
    casamet%tsoil = casamet%tsoil + ssnow%tgg
    casamet%moist = casamet%moist + ssnow%wb
    casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dtr8
    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dtr8
    if ( mod(ktau,npercasa)==0 ) then
      casamet%tairk=casamet%tairk/real(npercasa,8)
      casamet%tsoil=casamet%tsoil/real(npercasa,8)
      casamet%moist=casamet%moist/real(npercasa,8)
      xKNlimiting = 1._8
      if ( cable_climate==1 ) then
        call cable_phenology_clim(climate,phen,rad,veg)
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
        call casa_pop_firstcall(lalloc,casabiome,casaflux,casamet,casapool,phen,pop,soil,veg)  
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
        mp_POP = size(POP%pop_grid)  
        call casa_pop_secondcall(mp_POP,casaflux,pop)  
      end if
      call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
      call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet,lalloc)
      if ( ccycle<2 ) call casa_ndummy(casapool)
      if ( ccycle<3 ) call casa_pdummy(casapool)
      call casa_cnpbal(casapool,casaflux,casabal)
      casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp*deltpool
      casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp*deltpool
      casabal%FCrmleafyear = casabal%FCrmleafyear + casaflux%Crmplant(:,leaf)*deltpool
      casabal%FCrmwoodyear = casabal%FCrmwoodyear + casaflux%Crmplant(:,wood)*deltpool
      casabal%FCrmrootyear = casabal%FCrmrootyear + casaflux%Crmplant(:,xroot)*deltpool
      casabal%FCrgrowyear  = casabal%FCrgrowyear  + casaflux%Crgplant*deltpool
      casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp*deltpool
      casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil*deltpool
      casabal%FCneeyear = casabal%FCneeyear + (casaflux%Cnpp-casaflux%Crsoil)*deltpool
      casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool
      casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep*deltpool
      casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix*deltpool
      casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet*deltpool
      casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake*deltpool
      casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach*deltpool
      casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss*deltpool
      casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea*deltpool
      casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep*deltpool
      casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet*deltpool
      casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake*deltpool
      casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach*deltpool  
      casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss*deltpool 
      ! reset casamet for next call
      casamet%tairk = 0._8
      casamet%tsoil = 0._8
      casamet%moist = 0._8
      casaflux%cgpp = 0._8
      casaflux%crmplant(:,leaf) = 0._8
      ! update for POP
      if ( cable_pop==1 ) then
        if ( jmonth/=12 .or. idoy/=1 ) then
          ! do not update for extra day in leap year  
          casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp*casaflux%fracCalloc(:,2)*0.7_8
          casabal%LAImax = max(casamet%glai, casabal%LAImax)
          casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/365._8/1000._8
          casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/365._8/1000._8
        end if  
      else
        casaflux%stemnpp = 0._8
      end if
    end if
    canopy%frp  = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,xroot)+casaflux%crgplant(:))/real(casaperiod,8)
    canopy%frs  = casaflux%Crsoil(:)/real(casaperiod,8)
    canopy%frpw = casaflux%crmplant(:,wood)/real(casaperiod,8)
    canopy%frpr = casaflux%crmplant(:,xroot)/real(casaperiod,8)
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnpp = -canopy%fpn - canopy%frp
    if ( ccycle<=1 ) then
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    else
      if ( progvcmax>0 ) then
        canopy%fnee = canopy%fpn + canopy%frs + canopy%frp + casaflux%clabloss/real(casaperiod,8)
      else
        canopy%fnee = (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/real(casaperiod,8)
      end if
    end if
    
    sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dtr8
    sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dtr8
    sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dtr8
    sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dtr8
    sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dtr8
    sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dtr8
    sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dtr8
    sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dtr8
    sum_flux%sumrs  = sum_flux%sumrs  + canopy%frs*dtr8
  case default
    write(6,*) "ERROR: Unsupported carbon cycle option with ccycle=",ccycle
    stop
end select  

!--------------------------------------------------------------
! POP
if ( cable_pop==1 ) then
  ! update once per year
  if ( jmonth==12 .and. idoy==0 .and. mod(ktau,npercasa)==0 ) then
    call popdriver(casabal,casaflux,pop,veg)
    ! reset stemnpp, LAImax, Cleafmean and Crootmean for next call
    casaflux%stemnpp  = 0._8
    casabal%LAImax    = 0._8
    casabal%Cleafmean = 0._8
    casabal%Crootmean = 0._8
  end if
end if

!--------------------------------------------------------------
      
! Unpack tiles into grid point averages.
! Note that albsav and albnirsav are the VIS and NIR albedo output from CABLE to
! be used by the radiadiation scheme at the next time step.
do k = 1,ms
  where ( land(1:imax) )
    tgg(:,k) = 0.
    wb(:,k) = 0.
    wbice(:,k) = 0.
  end where
end do
do k = 1,3
  where ( land(1:imax) )
    tggsn(:,k) = 0.
    smass(:,k) = 0.
    ssdn(:,k) = 0.
  end where
end do
where ( land(1:imax) )
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
  zoq = 0.
  cduv = 0.
  cdtq = 0.
  ustar = 0.
  wetfac = 0.
  rsmin = 0.
  ssdnn = 0.
  snowd = 0.
  snage = 0.
end where
evspsbl_l = 0.
sbl_l = 0.
tmps = 0. ! average isflag

do nb = 1,maxnb
  is = tind(nb,1)
  ie = tind(nb,2)
  ! albedo
  albvisdir = albvisdir + unpack(sv(is:ie)*real(rad%reffbm(is:ie,1)),tmap(:,nb),0.)
  albnirdir = albnirdir + unpack(sv(is:ie)*real(rad%reffbm(is:ie,2)),tmap(:,nb),0.)
  albvisdif = albvisdif + unpack(sv(is:ie)*real(rad%reffdf(is:ie,1)),tmap(:,nb),0.)
  albnirdif = albnirdif + unpack(sv(is:ie)*real(rad%reffdf(is:ie,2)),tmap(:,nb),0.)
  ! fluxes
  fg = fg + unpack(sv(is:ie)*real(canopy%fh(is:ie)),tmap(:,nb),0.)
  eg = eg + unpack(sv(is:ie)*real(canopy%fe(is:ie)),tmap(:,nb),0.)
  ga = ga + unpack(sv(is:ie)*real(canopy%ga(is:ie)),tmap(:,nb),0.)
  rnet = rnet + unpack(sv(is:ie)*real(canopy%rnet(is:ie)),tmap(:,nb),0.)
  tss = tss + unpack(sv(is:ie)*real(rad%trad(is:ie)**4),tmap(:,nb),0.) ! ave longwave radiation
  ! drag and mixing
  zo   = zo   + unpack(sv(is:ie)/real(log(real(zmin,8)/rough%z0m(is:ie))**2),tmap(:,nb),0.)
  cduv = cduv + unpack(sv(is:ie)*real(canopy%cduv(is:ie)),tmap(:,nb),0.)
  cdtq = cdtq + unpack(sv(is:ie)*real(canopy%cdtq(is:ie)),tmap(:,nb),0.)
  ! soil
  do k = 1,ms
    tgg(:,k)   = tgg(:,k)   + unpack(sv(is:ie)*real(ssnow%tgg(is:ie,k)),  tmap(:,nb),0.)
    wb(:,k)    = wb(:,k)    + unpack(sv(is:ie)*real(ssnow%wb(is:ie,k)),   tmap(:,nb),0.)
    wbice(:,k) = wbice(:,k) + unpack(sv(is:ie)*real(ssnow%wbice(is:ie,k)),tmap(:,nb),0.)
  end do
  ! hydrology
  runoff = runoff + unpack(sv(is:ie)*real(ssnow%runoff(is:ie)*dtr8),tmap(:,nb),0.) ! convert mm/s to mm
  runoff_surface = runoff_surface + unpack(sv(is:ie)*real(ssnow%rnof1(is:ie)*dtr8),tmap(:,nb),0.) ! convert mm/s to mm
  fwet = fwet + unpack(sv(is:ie)*real(canopy%fwet(is:ie)),tmap(:,nb),0.)         ! used for aerosols
  wetfac = wetfac + unpack(sv(is:ie)*real(ssnow%wetfac(is:ie)),tmap(:,nb),0.)    ! used for aerosols
  cansto = cansto + unpack(sv(is:ie)*real(canopy%cansto(is:ie)),tmap(:,nb),0.)   ! not used
  ! snow
  tmps = tmps + unpack(sv(is:ie)*real(ssnow%isflag(is:ie)),tmap(:,nb),0.)  ! used in radiation (for nsib==3)
  do k = 1,3
    tggsn(:,k) = tggsn(:,k) + unpack(sv(is:ie)*real(ssnow%tggsn(is:ie,k)),tmap(:,nb),0.) ! for restart file
    smass(:,k) = smass(:,k) + unpack(sv(is:ie)*real(ssnow%smass(is:ie,k)),tmap(:,nb),0.) ! for restart file
    ssdn(:,k)  = ssdn(:,k)  + unpack(sv(is:ie)*real(ssnow%ssdn(is:ie,k)),tmap(:,nb),0.)  ! for restart file
  end do
  ssdnn = ssdnn + unpack(sv(is:ie)*real(ssnow%ssdnn(is:ie)),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snage = snage + unpack(sv(is:ie)*real(ssnow%snage(is:ie)),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snowd = snowd + unpack(sv(is:ie)*real(ssnow%snowd(is:ie)),tmap(:,nb),0.)
  snowmelt = snowmelt + unpack(sv(is:ie)*real(ssnow%smelt(is:ie)),tmap(:,nb),0.)
  evspsbl_l = evspsbl_l + unpack(sv(is:ie)*real((canopy%fev(is:ie)+canopy%fesp(is:ie) &
                                +canopy%fess(is:ie)/ssnow%cls(is:ie))/C%HL),tmap(:,nb),0.)
  sbl_l = sbl_l + unpack(sv(is:ie)*real(ssnow%evapsn(is:ie)/dtr8),tmap(:,nb),0.)
  ! diagnostic
  ! diagnostic
  epot = epot + unpack(sv(is:ie)*real(ssnow%potev(is:ie)),tmap(:,nb),0.)         ! diagnostic in history file
  vlai = vlai + unpack(sv(is:ie)*real(veg%vlai(is:ie)),tmap(:,nb),0.)
  rsmin = rsmin + unpack(sv(is:ie)*real(canopy%gswx_T(is:ie)),tmap(:,nb),0.)     ! diagnostic in history file
end do

if ( ccycle/=0 ) then
  cplant = 0.
  niplant = 0.
  pplant = 0.
  clitter = 0.
  nilitter = 0.
  plitter = 0.
  csoil = 0.
  nisoil = 0.
  psoil = 0.
  fnee = 0.
  fpn = 0.
  frd = 0.
  frp = 0.
  frpw = 0.
  frpr = 0.
  frs = 0.
  cnpp = 0.
  cnbp = 0.
  if ( diaglevel_carbon > 0 ) then
    fevc = 0.
    plant_turnover = 0.
    plant_turnover_wood = 0.
  end if
  do nb = 1,maxnb
    is = tind(nb,1)
    ie = tind(nb,2)
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
    fnee = fnee + unpack(sv(is:ie)*real(canopy%fnee(is:ie)),  tmap(:,nb),0.)
    fpn  = fpn  + unpack(sv(is:ie)*real(canopy%fpn(is:ie)),   tmap(:,nb),0.)
    frd  = frd  + unpack(sv(is:ie)*real(canopy%frday(is:ie)), tmap(:,nb),0.)
    frp  = frp  + unpack(sv(is:ie)*real(canopy%frp(is:ie)),   tmap(:,nb),0.)
    frpw = frpw + unpack(sv(is:ie)*real(canopy%frpw(is:ie)),  tmap(:,nb),0.)
    frpr = frpr + unpack(sv(is:ie)*real(canopy%frpr(is:ie)),  tmap(:,nb),0.)
    frs  = frs  + unpack(sv(is:ie)*real(canopy%frs(is:ie)),   tmap(:,nb),0.)
    cnpp = cnpp + unpack(sv(is:ie)*real(casaflux%cnpp(is:ie))/real(casaperiod),tmap(:,nb),0.)
    cnbp = cnbp + unpack(sv(is:ie)*real(casaflux%Crsoil(is:ie)-casaflux%cnpp(is:ie)-casapool%dClabiledt(is:ie)) &
        /real(casaperiod),tmap(:,nb),0.)
    if ( diaglevel_carbon > 0 ) then
      fevc = fevc + unpack(sv(is:ie)*real(canopy%fevc(is:ie)),  tmap(:,nb),0.)
      plant_turnover      = plant_turnover      + unpack(sv(is:ie)* &
                            real(sum(casaflux%Cplant_turnover(is:ie,:),2))/real(casaperiod),  tmap(:,nb),0.)
      plant_turnover_wood = plant_turnover_wood + unpack(sv(is:ie)* &
                            real(casaflux%Cplant_turnover(is:ie,2))/real(casaperiod),         tmap(:,nb),0.)
    end if
  end do
end if

if ( cable_climate==1 ) then
  ! just extract the first tile for now  
  climate_ivegt        = unpack(climate%iveg,tmap(:,1),0)
  climate_biome        = unpack(climate%biome,tmap(:,1),0)
  climate_min20        = real(unpack(climate%mtemp_min20,tmap(:,1),0._8),4)
  climate_max20        = real(unpack(climate%mtemp_max20,tmap(:,1),0._8),4)
  climate_alpha20      = real(unpack(climate%alpha_PT20,tmap(:,1),0._8),4)
  climate_agdd5        = real(unpack(climate%agdd5,tmap(:,1),0._8),4)
  climate_gmd          = unpack(climate%gmd,tmap(:,1),0)
  climate_dmoist_min20 = real(unpack(climate%dmoist_min20,tmap(:,1),0._8),4)
  climate_dmoist_max20 = real(unpack(climate%dmoist_max20,tmap(:,1),0._8),4)
end if

! MJT notes - ustar, cduv, fg and eg are passed to the boundary layer turbulence scheme
! zoh, zoq and zo are passed to the scrnout diagnostics routines
! rsmin is typically used by CTM

where ( land(1:imax) )
  zo        = zmin*exp(-1./sqrt(zo))
  zoh       = zo/7.4
  zoq       = zoh
  ustar     = sqrt(cduv)*vmod  
  cduv      = cduv*vmod           ! cduv is Cd*vmod in CCAM
  cdtq      = cdtq*vmod
  tss       = tss**0.25
  rsmin     = 1./rsmin
  evspsbl   = evspsbl + dt*evspsbl_l
  sbl       = sbl + dt*sbl_l
  ! update albedo and tss before calculating net radiation
  albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif ! for nrad=4
  albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif ! for nrad=4 
  isflag(:) = nint(tmps(:)) ! tmps is average isflag
end where
qsttg_land(:) = qsat(ps,tss,imax) ! must wait for tss to be updated first
where ( land(1:imax) )
  qsttg(:)  = qsttg_land(:)
end where

return
end subroutine sib4_work

! *************************************************************************************
subroutine cbmemiss(trsrc,mvegt,mode,tile,imax)
  
use newmpar_m
use parm_m

implicit none
  
integer, intent(in) :: tile,imax
integer, intent(in) :: mvegt,mode
real, dimension(imax), intent(out) :: trsrc

call cbmemiss_work(trsrc,mvegt,mode,imax,tdata(tile)%tind,tdata(tile)%tmap,tdata(tile)%maxnb)
  
return
end subroutine cbmemiss

subroutine cbmemiss_work(trsrc,mvegt,mode,imax,tind,tmap,maxnb)
  
use newmpar_m
use parm_m

implicit none
  
integer, intent(in) :: imax
integer, intent(in) :: mvegt,mode
integer nb
real, dimension(imax), intent(out) :: trsrc
integer, dimension(maxtile,2), intent(in) :: tind
logical, dimension(imax,maxtile), intent(in) :: tmap
integer, intent(in) :: maxnb
real, dimension(imax) :: fpn,frd,frp,frs
integer :: is, ie

if ( nsib/=6 .and. nsib/=7 ) then
  write(6,*) "ERROR: Attempted to read CABLE emissions with CABLE disabled"
  stop
end if

fpn=0.
frd=0.
frp=0.
frs=0.
do nb = 1,maxnb
  is = tind(nb,1)
  ie = tind(nb,2)
  where ( veg%iveg(is:ie)==mvegt )
    fpn = fpn + unpack(sv(is:ie)*real(canopy%fpn(is:ie)),  tmap(:,nb),0.)
    frd = frd + unpack(sv(is:ie)*real(canopy%frday(is:ie)),tmap(:,nb),0.)
    frp = frp + unpack(sv(is:ie)*real(canopy%frp(is:ie)),  tmap(:,nb),0.)
    frs = frs + unpack(sv(is:ie)*real(canopy%frs(is:ie)),  tmap(:,nb),0.)
  end where
end do
  
select case( mode )
  case(1)
    trsrc = fpn - frd
  case(2)
    trsrc = frp + frd
  case(3)
    trsrc = frs
  case default
    write(6,*) "ERROR: Unknown mode for cbmemiss ",mode
    stop
end select
  
return
end subroutine cbmemiss_work

! *************************************************************************************
subroutine setlai(sigmf,jmonth,jday,jhour,jmin,mp,sv,vl1,vl2,vl3,vl4,casamet,veg,imax,tind,tmap,maxnb)

use cc_mpi
use dates_m
use newmpar_m
use parm_m

implicit none
  
integer, intent(in) :: jmonth,jday,jhour,jmin,mp
integer, intent(in) :: imax
integer, optional :: maxnb
integer monthstart, nb, is, ie
integer, dimension(12), parameter :: imonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
integer, dimension(maxtile,2), intent(in), optional :: tind
real, dimension(mp), intent(in) :: sv, vl1, vl2, vl3, vl4
real, parameter :: vextkn = 0.4
real, dimension(imax), intent(out) :: sigmf
real, dimension(mp) :: a0, a1, a2, aa, bb, cc, mp1, mp2, c2, c3, c4
real, dimension(mp) :: dummy_unpack
real x
logical, dimension(imax,maxtile), intent(in), optional :: tmap
type(casa_met), intent(in) :: casamet
type(veg_parameter_type), intent(inout) :: veg

select case( proglai )
  case(-1) ! piece-wise linear interpolated LAI
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
      veg%vlai(:) = real(c3(:)*x + c2(:), 8)
    else
      c4(:) = 2.*vl2(:) - 0.5*mp1(:) - 0.5*mp2(:) ! mid-point value
      c3(:) = 2.*(mp2(:)-c4(:))                   ! gradient
      c2(:) = 2.*c4(:) - mp2(:)                   ! intercept
      veg%vlai(:) = real(c3(:)*x + c2(:), 8)
    end if
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.02_8 )
    elsewhere
      veg%vlai = 1.E-8_8
    end where
    
  case(0) ! PWCB interpolated LAI
    write(6,*) "ERROR: proglai=0 has been depreciated.  Please use proglai=-1"
    call ccmpi_abort(-1)

  case(1) ! prognostic LAI
    if ( ccycle/=2 .and. ccycle/=3 ) then
      write(6,*) "ERROR: CASA CNP LAI is not operational"
      write(6,*) "Prognostic LAI requires ccycle=2 or ccycle=3"
      call ccmpi_abort(-1)
    end if
    veg%vlai(:) = casamet%glai(:)
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.02_8 )
    elsewhere
      veg%vlai = 1.E-8_8
    end where

  case default
    write(6,*) "ERROR: Unknown proglai option ",proglai
    call ccmpi_abort(-1)
end select

! diagnose greeness fraction (e.g., for aerosols)  
sigmf(:) = 0.
if ( present(tind) .and. present(tmap) .and. present(maxnb) ) then
  !inside OpenMP region
  do nb = 1,maxnb
    is = tind(nb,1)
    ie = tind(nb,2)
    sigmf(:) = sigmf(:) + unpack(sv(is:ie)*(1.-exp(-vextkn*real(veg%vlai(is:ie)))),tmap(:,nb),0.)
  end do
else
  dummy_unpack = sv*(1.-exp(-vextkn*real(veg%vlai)))  
  call cable_unpack(dummy_unpack,sigmf)
endif
sigmf = min( sigmf, 1. )
  
return
end subroutine setlai

! *************************************************************************************
! Calculate prognostic vcmax from CASA-CNP
subroutine vcmax_feedback(casabiome,casamet,casapool,veg,climate,ktau)

use newmpar_m, only : mxvt
use parm_m, only : nperday

implicit none

type(casa_biome), intent(in) :: casabiome
type(casa_met), intent(in) :: casamet
type(casa_pool), intent(in) :: casapool
type(veg_parameter_type), intent(inout) :: veg
type(climate_type), intent(in) :: climate
integer, intent(in) :: ktau
integer np, ivt
real(kind=8) bjvref
real(kind=8), dimension(mp) :: ncleafx, npleafx, pleafx, nleafx
real(kind=8), dimension(mxvt), parameter :: xnslope = (/ 0.8,1.,2.,1.,1.,1.,0.5,1.,0.34,1.,1.,1.,1.,1.,1.,1.,1. /)

if ( progvcmax>0 .and. ccycle>=2 ) then

  ! initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx(:) = casabiome%ratioNPplantmin(veg%iveg(:),leaf)
  bjvref = 1.7_8 ! Walker 2014

  do np = 1,mp
    ivt = veg%iveg(np)
    if ( casamet%iveg2(np)/=icewater .and. casamet%glai(np)>casabiome%glaimin(ivt) .and. &
         casapool%cplant(np,leaf)>0._8 ) then

      if (icycle>1 .and. casapool%cplant(np,leaf)>0._8) then
         ncleafx(np) = min( casabiome%ratioNCplantmax(ivt,leaf),                  &
                       max( casabiome%ratioNCplantmin(ivt,leaf),                  &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf) ) )
      end if
      if ( ccycle>2 .and. casapool%pplant(np,leaf)>0._8 ) then
        npleafx(np) = min(30._8,max(8._8,casapool%nplant(np,leaf)/casapool%pplant(np,leaf)))
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
                    + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.e-6_8
          else
            if ( casapool%nplant(np,leaf)>0. .and. casapool%pplant(np,leaf)>0._8 ) then
              veg%vcmax(np) = ( casabiome%nintercept(ivt)           &
                   + casabiome%nslope(ivt)*(0.4_8+9._8/npleafx(np)) &
                   * ncleafx(np)/casabiome%sla(ivt) )*1.e-6_8
            else
              veg%vcmax(np) = ( casabiome%nintercept(ivt)                           &
                   + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.e-6_8
            end if
          end if
          veg%vcmax(np) =veg%vcmax(np)*xnslope(ivt)
        end if
        veg%ejmax = 2._8*veg%vcmax
      end do
  
    case(2)
      !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax – 
      !to leaf nitrogen, leaf phosphorus, and specific leaf area: 
      !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
      do np = 1,mp
        ivt = veg%iveg(np)  
        nleafx(np) = ncleafx(np)/casabiome%sla(ivt) ! leaf N in g N m-2 leaf
        pleafx(np) = nleafx(np)/npleafx(np)         ! leaf P in g P m-2 leaf

        if (ivt==7 .or. ivt==9  ) then
          veg%vcmax(np) = 1.0e-5_8 ! special for C4 grass: set here to value from  parameter file
          veg%ejmax(np) = 2.0_8 * veg%vcmax(np)
        elseif (ivt==1) then
          ! account here for spring recovery
          !veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))*casabiome%vcmax_scalar(ivt) &
          !               *climate%frec(np)
          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))*climate%frec(np)  
          veg%ejmax(np) =bjvref * veg%vcmax(np)
        else
          !veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))*casabiome%vcmax_scalar(ivt)
          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))
          veg%ejmax(np) =bjvref * veg%vcmax(np)
        endif
      end do

      if ( cable_user%finite_gm ) then
        ! vcmax and jmax modifications according to Sun et al. 2014 Table S3  
        where( veg%iveg==1 )
          veg%vcmax = veg%vcmax*2.2_8
          veg%ejmax = veg%vcmax*1.1_8
        elsewhere ( veg%iveg==2 )
          veg%vcmax = veg%vcmax*1.9_8
          veg%ejmax = veg%vcmax*1.2_8
        elsewhere ( veg%iveg==3 )
          veg%vcmax = veg%vcmax*1.4_8
          veg%ejmax = veg%vcmax*1.5_8
        elsewhere ( veg%iveg==4 )
          veg%vcmax = veg%vcmax*1.45_8
          veg%ejmax = veg%vcmax*1.3_8
        elsewhere ( veg%iveg==5 )
          veg%vcmax = veg%vcmax*1.7_8
          veg%ejmax = veg%vcmax*1.2_8
        elsewhere ( veg%iveg==6 .or. veg%iveg==8 .or. veg%iveg==9 )
          veg%vcmax = veg%vcmax*1.6_8
          veg%ejmax = veg%vcmax*1.2_8
        end where
      end if    
      
    case default
      write(6,*) "ERROR: Invalid progvcmax ",progvcmax
      stop
  end select
  
  if ( mod(ktau,nperday)==1 ) then    
    if ( cable_climate==1 .and. (progvcmax==2 .or. progvcmax==1) ) then
      if ( all(climate%tleaf_sun(:,1)>0.1_8) ) then  
        call optimise_JV(veg,climate,24,bjvref)
      end if
    end if  
  end if
  
end if
  
return
end subroutine vcmax_feedback

! *************************************************************************************
! Update phenology depending on climate
subroutine cable_phenology_clim(climate,phen,rad,veg)

implicit none

integer :: np, days, ivt
integer, parameter :: ndays_raingreenup = 60
real(kind=8) :: gdd0
real(kind=8) :: phengdd5ramp
real(kind=8) :: phen_tmp
real(kind=8), parameter :: k_chilla = 0._8, k_chillb = 100._8, k_chillk = 0.05_8
real(kind=8), parameter :: APHEN_MAX = 200._8
type(climate_type), intent(in) :: climate
type(phen_variable), intent(inout) :: phen
type(radiation_type), intent(in) :: rad
type(veg_parameter_type), intent(in) :: veg

!phen%doyphase(np,1) ! DOY for greenup
!phen%doyphase(np,2) ! DOY for steady LAI
!phen%doyphase(np,3) ! DOY for leaf senescence
!phen%doyphase(np,4) ! DOY for minimal LAI season

do np = 1,mp
  ivt = veg%iveg(np)  
  phen_tmp = 1.0_8
  
  ! evergreen pfts
  if ( ivt==1 .or. ivt==2 .or. ivt==5 ) then
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
    gdd0 = k_chilla + k_chillb*exp(-k_chillk*real(climate%chilldays(np),8))
    phengdd5ramp = 200._8
    if ( climate%gdd5(np)>gdd0 .and. phen%aphen(np)<APHEN_MAX) then
      phen_tmp = min(1._8, (climate%gdd5(np)-gdd0)/phengdd5ramp)
    else
      phen_tmp = 0._8  
    end if
  end if

  ! summergreen grass or crops
  if ( ivt>=6 .and. ivt<=10 )  then     ! grass or crops
    if ( climate%gdd5(np)>0.1_8 ) then  
      phengdd5ramp = 200._8
      phen_tmp = min(1._8, climate%gdd5(np)/phengdd5ramp)
    else
      phen_tmp = 0._8
    end if  
  end if
  
  ! raingreen pfts
  if ( ivt>=6 .and. ivt<=10 ) then
    if ( climate%gmd(np)>=1 .and. climate%gmd(np)<ndays_raingreenup ) then    
      phen_tmp = min( phen_tmp, 0.99_8 )
    else if ( climate%gmd(np)==0 ) then
      phen_tmp = 0._8
    else if ( climate%gmd(np)>=ndays_raingreenup ) then
      phen_tmp = min( phen_tmp, 1._8 )
    end if
  end if  

  if ( (ivt==3.or.ivt==4) .or. (ivt>=6.and.ivt<=10) ) then

    if (phen_tmp>0._8 .and. (phen_tmp<1._8)) then
      phen%phase(np) = 1 ! greenup
      phen%doyphase(np,1) = climate%doy
    elseif (phen_tmp>=1._8) then
      phen%phase(np) = 2 ! steady LAI
      phen%doyphase(np,2) = climate%doy
    elseif (phen_tmp<0.0000001_8) then
      phen%phase(np) = 3 ! senescence
      phen%doyphase(np,3) = climate%doy
    endif

    if ( phen%phase(np)==3 ) then
       days = min(climate%doy,365) - phen%doyphase(np,3)
       if ( days<0 ) days = days + 365
       if ( days>14 ) phen%phase(np) = 0          ! mimimum LAI
    endif

    ! Update annual leaf-on sum
    if ( (rad%latitude(np)>=0._8 .and. climate%doy==COLDEST_DAY_NHEMISPHERE) .or. &
         (rad%latitude(np)<0._8 .and. climate%doy==COLDEST_DAY_SHEMISPHERE) ) then
      phen%aphen(np) = 0._8
    end if
    phen%phen(np) = phen_tmp
    phen%aphen(np) = phen%aphen(np) + phen%phen(np)

  end if

end do  ! end loop over patches

return
end subroutine cable_phenology_clim

! *************************************************************************************
! POP subroutines
subroutine casa_pop_firstcall(lalloc,casabiome,casaflux,casamet,casapool,phen,pop,soil,veg)

implicit none

integer, intent(in) :: lalloc
type(casa_biome), intent(inout) :: casabiome
type(casa_flux), intent(inout) :: casaflux
type(casa_met), intent(inout) :: casamet
type(casa_pool), intent(inout) :: casapool
type(phen_variable), intent(inout) :: phen
type(pop_type), intent(inout) :: pop
type(soil_parameter_type), intent(inout) :: soil
type(veg_parameter_type), intent(inout) :: veg

call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,lalloc)

if ( POP%np>0 ) then

  where ( pop%pop_grid(:)%cmass_sum_old>0.001_8 .and. pop%pop_grid(:)%cmass_sum>0.001_8 )
    casaflux%frac_sapwood(POP%Iwood) = real( POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum, 8)
    casaflux%sapwood_area(POP%Iwood) = real( max(POP%pop_grid(:)%sapwood_area/10000._dp, 1.e-6_dp), 8)
    veg%hc(POP%Iwood) = real( POP%pop_grid(:)%height_max, 8)
    where ( pop%pop_grid(:)%LU==2 )
      casaflux%kplant(POP%Iwood,2) = real( 1._dp -            &
        (1._dp-  max( min((POP%pop_grid(:)%stress_mortality + &
        POP%pop_grid(:)%crowding_mortality                    &
        + POP%pop_grid(:)%fire_mortality )                    &
        /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
        1._dp/veg%disturbance_interval(POP%Iwood,1), 0.99_dp), 0._dp))**(1._dp/365._dp), 8)
    elsewhere
      casaflux%kplant(POP%Iwood,2) = real( 1._dp -                        &
        (1._dp-  max( min((POP%pop_grid(:)%stress_mortality +             &
        POP%pop_grid(:)%crowding_mortality                                &
        + POP%pop_grid(:)%fire_mortality+POP%pop_grid(:)%cat_mortality  ) &
        /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99_dp), 0._dp))**(1._dp/365._dp), 8)
    end where
    veg%hc(POP%Iwood) = real( POP%pop_grid(:)%height_max, 8)
  elsewhere
    casaflux%frac_sapwood(POP%Iwood) = 1._8
    casaflux%sapwood_area(POP%Iwood) = real( max(POP%pop_grid(:)%sapwood_area/10000._dp, 1.e-6_dp), 8)
    casaflux%kplant(POP%Iwood,2) = 0._8
    veg%hc(POP%Iwood) = real( POP%pop_grid(:)%height_max, 8)
  end where
  
end if

return
end subroutine casa_pop_firstcall

subroutine casa_pop_secondcall(mp_POP,casaflux,pop)

implicit none

integer, intent(in) :: mp_POP
type(casa_flux), intent(inout) :: casaflux
type(pop_type), intent(inout) :: pop
real(kind=dp), dimension(mp_POP) :: tmp

if ( mp_POP>0 ) then

  tmp = (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality &
        +POP%pop_grid(:)%cat_mortality + POP%pop_grid(:)%fire_mortality  )
  where ( tmp>1.e-12_dp )
    casaflux%Cplant_turnover_disturbance(POP%Iwood) =                         &
         real( casaflux%Cplant_turnover(POP%Iwood,2)*(POP%pop_grid(:)%cat_mortality &
         + POP%pop_grid(:)%fire_mortality  )/tmp, 8)
    casaflux%Cplant_turnover_crowding(POP%Iwood) =                                     &
         real( casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%crowding_mortality/tmp, 8)
    casaflux%Cplant_turnover_resource_limitation(POP%Iwood) =                          &
         real( casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%stress_mortality/tmp, 8)
  end where
  
end if

return
end subroutine casa_pop_secondcall

subroutine POPdriver(casabal,casaflux,pop,veg)

implicit none
  
type(casa_balance), intent(in) :: casabal
type(casa_flux), intent(in) :: casaflux
type(pop_type), intent(inout) :: pop
type(veg_parameter_type), intent(in) :: veg
integer, dimension(POP%np) :: Iw
integer(kind=i4b), dimension(POP%np,2) :: disturbance_interval_tmp
real(kind=dp), dimension(POP%np,2) :: disturbance_intensity_tmp, stemNPP_tmp
real(kind=dp), dimension(POP%np) :: LAImax_tmp, Cleafmean_tmp, Crootmean_tmp, NPPtoGPP_tmp
real(kind=dp), dimension(mp) :: NPPtoGPP

if ( cable_pop==1 .and. POP%np>0 ) then ! CALL_POP
  Iw = POP%Iwood
  where ( casabal%FCgppyear>1.e-5_8 .and. casabal%FCnppyear>1.e-5_8 )
    NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
  elsewhere
    NPPtoGPP = 0.5_dp
  end where
  stemNPP_tmp(:,1) = real(max(casaflux%stemnpp(Iw)/1000._8,0.0001_8),dp)
  stemNPP_tmp(:,2) = 0.0001_dp !PAR to be consistent with original
  disturbance_interval_tmp = int(veg%disturbance_interval(Iw,:), i4b) 
  disturbance_intensity_tmp = real(veg%disturbance_intensity(Iw,:),8)
  LAImax_tmp = real(max(casabal%LAImax(Iw),0.001_8),dp)
  Cleafmean_tmp = real(casabal%cleafmean(Iw),dp)
  Crootmean_tmp = real(casabal%Crootmean(Iw),dp)
  NPPtoGPP_tmp = NPPtoGPP(Iw)

  call POPStep(pop, stemNPP_tmp, disturbance_interval_tmp,       & 
        disturbance_intensity_tmp,                               &
        LAImax_tmp, Cleafmean_tmp, Crootmean_tmp, NPPtoGPP_tmp)

end if ! CALL_POP

return
end subroutine POPdriver

! *************************************************************************************
! track climate feedback into CABLE
subroutine cableclimate(idoy,imonth,ndoy,canopy,climate,met,rad,ssnow,veg, &
                        climate_save,npercasa,ktau)

use parm_m, only : dt, nperhr

implicit none

integer, intent(in) :: idoy, imonth, npercasa, ktau
integer d, y, nsd
integer, dimension(13), intent(in) :: ndoy
real(kind=8), dimension(mp) :: mtemp_last
real(kind=8), dimension(mp) :: RhoA, PPc, EpsA, phiEq, tmp
real(kind=8), dimension(mp) :: frec0, f1, f2
type(canopy_type), intent(in) :: canopy
type(climate_type), intent(inout) :: climate
type(met_type), intent(in) :: met
type(radiation_type), intent(inout) :: rad
type(soil_snow_type), intent(in) :: ssnow
type(veg_parameter_type), intent(in) :: veg
type(climate_save_type), intent(in) :: climate_save 

real(kind=8), parameter :: Gaero   = 0.015_8    ! (m s-1) aerodynmaic conductance (for use in PT evap)
real(kind=8), parameter :: Capp    = 29.09_8    ! isobaric spec heat air    [J/molA/K]
real(kind=8), parameter :: SBoltz  = 5.67e-8_8  ! Stefan-Boltzmann constant [W/m2/K4]
real(kind=8), parameter :: moisture_min = 0.15  ! threshold for setting "growing moisture days",
                                                ! as required for drought-deciduous phenology
real(kind=8), parameter :: CoeffPT = 1.26_8
real(kind=8), parameter :: T1 = 0._8
real(kind=8), parameter :: T2 = -3._8
real(kind=8), parameter :: T3 = -4._8
real(kind=8), parameter :: T6 = -5._8
real(kind=8), parameter :: ffrost = 0.1_8
real(kind=8), parameter :: fdorm0 = 0.15_8
real(kind=8), parameter :: gdd0_rec0 = 500.

if ( cable_climate==0 ) return

nsd = 24*5 ! number of subdirunal time-steps to be accumulated
climate%doy = idoy

! accumulate annual evaporation and potential evaporation
RhoA  = met%pmb*100._8/(8.314_8*met%tk) ! air density [molA/m3]
PPc   = Gaero/(Gaero+4._8*SBoltz*(met%tk**3/(RhoA*Capp)))
tmp   = met%tk - 273.16_8 ! avoid array temporary
EpSA  = Epsif(tmp, met%pmb)
phiEq = canopy%rniso*(PPc*EpsA)/(PPc*EpsA+1._8) ! equil ltnt heat flux  [W/m2]

climate%evap_PT = climate%evap_PT + phiEq*CoeffPT/2.5014e6*dt  ! mm
climate%aevap   = climate%aevap + met%precip ! mm
climate%dtemp   = climate%dtemp + met%tk - 273.15_8
climate%dmoist  = climate%dmoist + sum(ssnow%wb(:,:)*veg%froot(:,:),2)

if ( mod(ktau,npercasa)==1 ) climate%dtemp_min = met%tk - 273.15_8
climate%dtemp_min = min( climate%dtemp_min, met%tk - 273.15_8 )

if ( mod(ktau,nperhr)==1 ) then
  climate_save%apar_leaf_sun = 0._8
  climate_save%apar_leaf_shade = 0._8
  climate_save%dleaf_sun = 0._8
  climate_save%fwsoil = 0._8
  climate_save%dleaf_shade = 0._8
  climate_save%tleaf_sun = 0._8
  climate_save%tleaf_shade = 0._8
  climate_save%cs_sun = 0._8
  climate_save%cs_shade = 0._8
  climate_save%scalex_sun = 0._8
  climate_save%scalex_shade = 0._8
end if
climate_save%apar_leaf_sun = climate_save%apar_leaf_sun + rad%qcan(:,1,1)*4.6_8 ! umol m-2 s-1
climate_save%apar_leaf_shade = climate_save%apar_leaf_shade + rad%qcan(:,2,1)*4.6_8 ! umod m-2 s-1
climate_save%dleaf_sun = climate_save%dleaf_sun + canopy%dlf
climate_save%fwsoil = climate_save%fwsoil + canopy%fwsoil
climate_save%dleaf_shade = climate_save%dleaf_shade + canopy%dlf
climate_save%tleaf_sun = climate_save%tleaf_sun + canopy%tlf
climate_save%tleaf_shade = climate_save%tleaf_shade + canopy%tlf
climate_save%cs_sun = climate_save%cs_sun + canopy%cs_sl ! ppm
climate_save%cs_shade = climate_save%cs_shade + canopy%cs_sh ! ppm
climate_save%scalex_sun = climate_save%scalex_sun + rad%scalex(:,1)
climate_save%scalex_shade = climate_save%scalex_shade + rad%scalex(:,2)

! midday fraction of incoming visible radiation absorbed by the canopy
!if ( mod(ktau,npercasa) == npercasa/2 ) THEN
!  where (rad%fbeam(:,1)>=0.01)
!    climate%fapar_ann_max = max(1.- rad%extkbm(:,1)*canopy%vlaiw, climate%fapar_ann_max)
!  endwhere
!endif

! accumulate sub-diurnal sub- and shade-leaf net variables that are revant for calc of Anet
if ( mod(ktau,nperhr)==0 ) then
  climate%APAR_leaf_sun(:,1:nsd-1)   = climate%APAR_leaf_sun(:,2:nsd)
  climate%APAR_leaf_shade(:,1:nsd-1) = climate%APAR_leaf_shade(:,2:nsd)
  climate%Dleaf_sun(:,1:nsd-1)       = climate%Dleaf_sun(:,2:nsd)
  climate%fwsoil(:,1:nsd-1)          = climate%fwsoil(:,2:nsd)
  climate%Dleaf_shade(:,1:nsd-1)     = climate%Dleaf_shade(:,2:nsd)
  climate%Tleaf_sun(:,1:nsd-1)       = climate%Tleaf_sun(:,2:nsd)
  climate%Tleaf_shade(:,1:nsd-1)     = climate%Tleaf_shade(:,2:nsd)
  climate%cs_sun(:,1:nsd-1)          = climate%cs_sun(:,2:nsd)
  climate%cs_shade(:,1:nsd-1)        = climate%cs_shade(:,2:nsd)
  climate%scalex_sun(:,1:nsd-1)      = climate%scalex_sun(:,2:nsd)
  climate%scalex_shade(:,1:nsd-1)    = climate%scalex_shade(:,2:nsd)
  climate%APAR_leaf_sun(:,nsd)   = climate_save%apar_leaf_sun/real(nperhr,8)
  climate%APAR_leaf_shade(:,nsd) = climate_save%apar_leaf_shade/real(nperhr,8)
  climate%Dleaf_sun(:,nsd)       = climate_save%dleaf_sun/real(nperhr,8)
  climate%fwsoil(:,nsd)          = climate_save%fwsoil/real(nperhr,8)
  climate%Dleaf_shade(:,nsd)     = climate_save%dleaf_shade/real(nperhr,8)
  climate%Tleaf_sun(:,nsd)       = climate_save%tleaf_sun/real(nperhr,8)
  climate%Tleaf_shade(:,nsd)     = climate_save%tleaf_shade/real(nperhr,8)
  climate%cs_sun(:,nsd)          = climate_save%cs_sun/real(nperhr,8)
  climate%cs_shade(:,nsd)        = climate_save%cs_shade/real(nperhr,8)
  climate%scalex_sun(:,nsd)      = climate_save%scalex_sun/real(nperhr,8)
  climate%scalex_shade(:,nsd)    = climate_save%scalex_shade/real(nperhr,8)
end if
  
! accumulate daily temperature, evap and potential evap
if ( mod(ktau,npercasa)==0 .and. ktau>0 ) then

  climate%dtemp = climate%dtemp/real(npercasa,8)
  climate%dmoist = climate%dmoist/real(npercasa,8)

  ! In midwinter, reset GDD counter for summergreen phenology  
  where ( (rad%latitude>=0._8.and.idoy==COLDEST_DAY_NHEMISPHERE) .or. &
          (rad%latitude<0._8.and.idoy==COLDEST_DAY_SHEMISPHERE) )
    !climate%gdd5 = 0._8
    climate%gdd0 = 0._8
    ! reset day degree sum related to spring photosynthetic recovery
    climate%gdd0_rec = 0._8
  end where
          
  ! In minsumer, reset dormancy fraction        
  where ( (rad%latitude<=0._8 .and. idoy==COLDEST_DAY_NHEMISPHERE) .or. &
          (rad%latitude>=0._8 .and. idoy==COLDEST_DAY_SHEMISPHERE) )
    climate%fdorm = 1._8
  end where  

  where ( climate%dtemp_min>=T6 .and. climate%dtemp_min<=0._8 )
    climate%fdorm = max( climate%fdorm-ffrost*climate%dtemp_min/T6, 0._8 )
  elsewhere ( climate%dtemp_min<T6 )
    climate%fdorm = max( climate%fdorm-ffrost, 0._8 )
  end where
  
  frec0 = fdorm0 + (1._8-fdorm0)*climate%fdorm
          
  where ( climate%dtemp_min>=T1 .and. climate%dtemp>=T1 )
    f2 = 1._8
  elsewhere ( climate%dtemp_min<=T2 .and. climate%dtemp<=T2 ) 
    f2 = 0._8
  elsewhere
    f2 = min( (climate%dtemp_min-T2)/(T1-T2), (climate%dtemp-T2)/(T1-T2) )  
  end where
  
  where ( climate%dtemp_min>=T2 )
    f1 = 0._8
  elsewhere ( climate%dtemp_min<=T3 )
    f1 = 0.3_8
  elsewhere
    f1 = 0.3_8*(T2-climate%dtemp_min)/(T2-T3)
  end where  
  
  where ( climate%dtemp_min>=T2 )
    climate%gdd0_rec = max( climate%gdd0_rec+climate%dtemp*f2, 0._8 )
  elsewhere ( climate%dtemp_min<T2 )
    climate%gdd0_rec = max( climate%gdd0_rec*(1._8-f1), 0._8 )
  end where  
  
  where ( climate%gdd0_rec<=gdd0_rec0 )
    climate%frec = frec0 + (1._8-frec0)*climate%gdd0_rec/gdd0_rec0
  elsewhere
    climate%frec = 1._8
  end where  
  
  ! Update GDD counters and chill day count
  climate%gdd0  = climate%gdd0 + max(0._8,climate%dtemp)
  climate%agdd0 = climate%agdd0 + max(0._8,climate%dtemp)
  climate%gdd5  = climate%gdd5 + max(0._8,climate%dtemp-5.)
  climate%agdd5 = climate%agdd5 + max(0._8,climate%dtemp-5.)
  climate%dmoist_min = min( climate%dmoist_min, climate%dmoist )
  climate%dmoist_max = max( climate%dmoist_max, climate%dmoist )
  where ( climate%dtemp<5. )
    climate%chilldays = min(climate%chilldays + 1, 365)
  end where

  ! update GMD (growing moisture day) counter
  where ( climate%dmoist>climate%dmoist_min20+moisture_min*(climate%dmoist_max20-climate%dmoist_min20) .and. &
          climate%dmoist_max20>climate%dmoist_min20 )
    climate%gmd = climate%gmd + 1
  elsewhere
    climate%gmd = 0
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
  climate%qtemp = climate%qtemp/91._8   ! average temperature over the last quarter
  climate%mtemp = climate%mtemp/31._8   ! average temperature over the last month
  climate%mmoist = climate%mmoist/31._8 ! average moisture index over the last month

  ! Reset GDD and chill day counter if mean monthly temperature falls below base
  ! temperature
  where ( mtemp_last>=5._8 .and. climate%mtemp<5._8 )
    climate%gdd5 = 0._8
    climate%chilldays = 0._8
  end where
  
  ! check if end of month
  if ( idoy==ndoy(imonth+1) ) then

    ! Update mean temperature for the last 12 months
    ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)
    !climate%atemp_mean = climate%atemp_mean*(11./12.) + climate%mtemp*(1./12.)

    ! Record minimum and maximum monthly temperatures
    if ( imonth==1 ) then
      climate%mtemp_min = climate%mtemp
      climate%mtemp_max = climate%mtemp
      climate%qtemp_max_last_year = climate%qtemp_max
      climate%qtemp_max = climate%qtemp
    end if
    climate%mtemp_min = min(climate%mtemp, climate%mtemp_min)
    climate%mtemp_max = max(climate%mtemp, climate%mtemp_max)
    climate%qtemp_max = max(climate%qtemp, climate%qtemp_max)

    ! On 31 December update records of minimum monthly temperatures for the last
    ! 20 years and find minimum monthly temperature for the last 20 years
    if ( imonth==12 ) then

      climate%alpha_PT = max(climate%aevap/climate%evap_PT, 0._8) ! ratio of annual evap to annual PT evap  
      !climate%fapar_ann_max_last_year = climate%fapar_ann_max 
        
      do y = 1,19
        climate%mtemp_min_20(:,y) = climate%mtemp_min_20(:,y+1)
        climate%mtemp_max_20(:,y) = climate%mtemp_max_20(:,y+1)
        climate%alpha_PT_20(:,y) = climate%alpha_PT_20(:,y+1)
        climate%dmoist_min_20(:,y) = climate%dmoist_min_20(:,y+1)
        climate%dmoist_max_20(:,y) = climate%dmoist_max_20(:,y+1)
      end do
      climate%mtemp_min_20(:,20) = climate%mtemp_min
      climate%mtemp_max_20(:,20) = climate%mtemp_max
      climate%alpha_PT_20(:,20) = climate%alpha_PT
      climate%dmoist_min_20(:,20) = climate%dmoist_min
      climate%dmoist_max_20(:,20) = climate%dmoist_max

      climate%mtemp_min20 = 0._8
      climate%mtemp_max20 = 0._8
      climate%alpha_PT20 = 0._8
      climate%dmoist_min20 = 0._8
      climate%dmoist_max20 = 0._8
      climate%nyears = 0
      do y = 1,20      
        if ( .not.(all(climate%mtemp_min_20(:,y)<0.00000001_8).and.all(climate%mtemp_max_20(:,y)<0.00000001_8)) ) then
          climate%nyears = climate%nyears + 1
          climate%mtemp_min20 = climate%mtemp_min20 + climate%mtemp_min_20(:,y)
          climate%mtemp_max20 = climate%mtemp_max20 + climate%mtemp_max_20(:,y)
          climate%alpha_PT20 = climate%alpha_PT20 + climate%alpha_PT_20(:,y)
          climate%dmoist_min20 = climate%dmoist_min20 + climate%dmoist_min_20(:,y)
          climate%dmoist_max20 = climate%dmoist_max20 + climate%dmoist_max_20(:,y)
        end if
      end do
      
      if ( climate%nyears>0 ) then
        climate%mtemp_min20 = climate%mtemp_min20/real(climate%nyears,8)
        climate%mtemp_max20 = climate%mtemp_max20/real(climate%nyears,8)
        climate%alpha_PT20 = climate%alpha_PT20/real(climate%nyears,8)
        climate%dmoist_min20 = climate%dmoist_min20/real(climate%nyears,8)
        climate%dmoist_max20 = climate%dmoist_max20/real(climate%nyears,8)
        call biome1_pft(climate,rad,mp)
      end if  
 
      ! reset climate data for next year
      climate%agdd5      = 0._8
      climate%agdd0      = 0._8
      climate%evap_PT    = 0._8     ! annual PT evap [mm]
      climate%aevap      = 0._8     ! annual evap [mm]
      !climate%fapar_ann_max = 0.0 
      climate%dmoist_min = climate%dmoist
      climate%dmoist_max = climate%dmoist
      
    end if  ! last month of year
  end if    ! last day of month

  ! reset climate data for next day
  climate%dtemp   = 0._8
  climate%dmoist  = 0._8
  
end if      ! end of day

return
end subroutine cableclimate

function epsif(tc,pmb) result(epsif_out)

implicit none

real(kind=8), dimension(mp), intent(in) :: tc, pmb
real(kind=8), dimension(mp) :: epsif_out
real(kind=8), dimension(mp) :: tctmp, es, desdt
real(kind=8), parameter:: A = 6.106_8      ! Teten coefficients
real(kind=8), parameter:: B = 17.27_8      ! Teten coefficients
real(kind=8), parameter:: C = 237.3_8      ! Teten coefficients
real(kind=8), parameter:: Rlat = 44140._8  ! lat heat evap H2O at 20C  [J/molW]
real(kind=8), parameter:: Capp = 29.09_8   ! isobaric spec heat air    [J/molA/K]

TCtmp = max( min( TC, 100._8), -40._8 ) ! constrain TC to (-40.0,100.0)
ES    = A*exp(B*TCtmp/(C+TCtmp))        ! sat vapour pressure
dESdT = ES*B*C/(C+TCtmp)**2             ! d(sat VP)/dT: (mb/K)
Epsif_out = (Rlat/Capp)*dESdT/Pmb       ! dimensionless (ES/Pmb = molW/molA)

return
end function epsif

subroutine biome1_pft(climate,rad,mp)

implicit none

type(climate_type), intent(inout) :: climate
type(radiation_type), intent(in) :: rad
integer, intent(in) :: mp
integer :: iq
integer, dimension(mp) :: npft
integer, dimension(mp,4) :: pft_biome1
real(kind=8) :: alpha_pt_scaled

! TABLE 1 , Prentice et al. J. Biogeog., 19, 117-134, 1992
! pft_biome1: Trees (1)tropical evergreen; (2) tropical raingreen; (3) warm temp evergreen ;
! (4) temperate summergreen; (5) cool-temperate conifer; (6) boreal evergreen conifer;
! (7) boreal summergreen;  
! Non-trees: (8) sclerophyll/succulent; (9) warm grass/shrub; (10) cool grass/shrub; 
! (11) cold grass/shrub; (12) hot desert shrub; (13) cold desert shrub.

pft_biome1 = 999
climate%biome = 999
climate%iveg = 0

do iq = 1,mp

   alpha_PT_scaled =  min(climate%alpha_PT20(iq), 1._8)

   ! tropical evergreen (1)
   ! tropical raingreen (2)
   if ( climate%mtemp_min20(iq)>=15.5_8 ) then
     if ( alpha_PT_scaled>=0.85_8 ) then
       pft_biome1(iq,1) = 1
       if ( alpha_PT_scaled<=0.90_8 ) then
         pft_biome1(iq,2) = 2
       end if
     else if ( alpha_PT_scaled>=0.4_8 .and. alpha_pt_scaled<0.85 ) then
       pft_biome1(iq,1) = 2
     end if
   end if
   
   if ( climate%mtemp_min20(iq)>=5._8 .and. alpha_PT_scaled>=0.4_8 .and. &
        pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 3
   end if
   
   if ( climate%mtemp_min20(iq)>=-15._8 .and. climate%mtemp_min20(iq)<=15.5_8 .and. &
        alpha_PT_scaled>=0.35_8 .and. climate%agdd5(iq)>1200. .and.                  &
        pft_biome1(iq,1)>3 ) then
     pft_biome1(iq,1) = 4
   end if
   
   if ( climate%mtemp_min20(iq)>=-19._8 .and. climate%mtemp_min20(iq)<=5._8 .and. &
        alpha_PT_scaled>=0.35 .and. climate%agdd5(iq)>900. ) then
     if (pft_biome1(iq,1)>4 ) then
       pft_biome1(iq,1) = 5
     else if (pft_biome1(iq,1)==4 ) then
       pft_biome1(iq,2) = 5
     end if
   end if
   
   if ( climate%mtemp_min20(iq)>=-35._8 .and. climate%mtemp_min20(iq)<=-2._8 .and. &
        alpha_PT_scaled>=0.35_8 .and. climate%agdd5(iq)>550. ) then
     if ( pft_biome1(iq,1)==999 ) then
       pft_biome1(iq,1) = 6
     else if ( pft_biome1(iq,2)==999 ) then
       pft_biome1(iq,2) = 6
     else
       pft_biome1(iq,3) = 6
     end if
   end if
   
   if ( climate%mtemp_min20(iq)<=5._8 .and. alpha_PT_scaled>=0.35_8 .and. &
        climate%agdd5(iq)>550. ) then
     if (pft_biome1(iq,1)==999 ) then
       pft_biome1(iq,1) = 7
     else if (pft_biome1(iq,2)==999 ) then
       pft_biome1(iq,2) = 7
     else if (pft_biome1(iq,3)==999 ) then
       pft_biome1(iq,3) = 7  
     else
       pft_biome1(iq,4) = 7
     end if
   end if
   
   ! sclerophyll/succulent (8)     
   if ( climate%mtemp_min20(iq)>=5._8 .and. alpha_PT_scaled>=0.2_8 .and. &
        pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 8
   end if
   
   if (climate%mtemp_max20(iq)>=22._8 .and.alpha_PT_scaled>=0.1_8 &
        .and. pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 9
   end if
   
   if ( climate%agdd5(iq)>=500. .and.alpha_PT_scaled>=0.33_8 .and. &
        pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 10
   end if
   
   if ( climate%agdd0(iq)>=100. .and. alpha_PT_scaled>=0.33_8 ) then 
     if ( pft_biome1(iq,1)==999 ) then
       pft_biome1(iq,1) = 11
     else if ( pft_biome1(iq,1)==10 ) then
       pft_biome1(iq,2) = 11
     end if
   end if
   
   ! hot desert shrub (12)
   if ( climate%mtemp_max20(iq)>=22._8 .and. pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 12
   end if
   
   if ( climate%agdd0(iq)>=100. .and. pft_biome1(iq,1)==999 ) then
     pft_biome1(iq,1) = 13
   end if
   
end do

! end of evironmental constraints on pft
npft = count( pft_biome1(:,:)/=999, dim=2 )
  
! MAP to Biome1 biome and CABLE pft
  
! (1) Tropical Rainforest
where ( pft_biome1(:,1)==1 .and. npft==1 )
    climate%biome(:) = 1
    climate%iveg(:) = 2
end where

! (2) Tropical Seasonal forest
where ( pft_biome1(:,1)==1 .and. pft_biome1(:,2)==2 .and. npft==2 )
  climate%biome(:) = 2
  climate%iveg(:) = 2
end where

! (3) Tropical dry forest/savanna
where ( pft_biome1(:,1)==2 .and. npft==1 )
  climate%biome(:) = 3
  climate%iveg(:) = 2  ! N.B. need to include c4 grass
end where

! (4) Broad-leaved evergreen/warm mixed-forest
where ( pft_biome1(:,1)==3 .and. npft==1 )
  climate%biome(:) = 4
  climate%iveg(:) = 2
end where

! (5) Temperate deciduous forest
where ( pft_biome1(:,1)==4 .and. pft_biome1(:,2)==5 .and. &
        pft_biome1(:,3)==7 .and. npft==3 )
  climate%biome(:) = 5
  climate%iveg(:) = 4
end where

! (6) Cool mixed forest
where ( pft_biome1(:,1)==4 .and. pft_biome1(:,2)==5 .and. &
        pft_biome1(:,3)==6 .and. pft_biome1(:,4)==7 .and. &
        npft==4 )
  climate%biome(:) = 6
  climate%iveg(:) = 4
end where

! (7) Cool conifer forest
where ( pft_biome1(:,1)==5 .and. pft_biome1(:,2)==6 .and. &
        pft_biome1(:,3)==7 .and. npft==3 )
  climate%biome(:) = 7
  climate%iveg(:) = 1
end where
      
! (8) Taiga
where ( pft_biome1(:,1)==6 .and. pft_biome1(:,2)==7 .and. npft==2 )
  climate%biome(:) = 8
  climate%iveg(:) = 1
end where

! (9) Cold mixed forest
where ( pft_biome1(:,1)==5 .and. pft_biome1(:,2)==7 .and. npft==2 )
  climate%biome(:) = 9
  climate%iveg(:) = 1
end where

! (10) Cold deciduous forest
where ( pft_biome1(:,1)==7 .and. npft==1 )
  climate%biome(:) = 10
  climate%iveg(:) = 3
end where

! (11) Xerophytic woods/scrub
where ( pft_biome1(:,1)==8 .and. npft==1 )
  climate%biome(:) = 11
  climate%iveg(:) = 5
end where

! (12) Warm grass/shrub
where ( pft_biome1(:,1)==9 .and. npft==1 )
  climate%biome(:) = 12
  climate%iveg(:) = 5  ! include C4 grass tile ?
end where

! (13) Cool grass/shrub
where ( pft_biome1(:,1)==10 .and. pft_biome1(:,2)==11 .and. npft==2 )
  climate%biome(:) = 13
  climate%iveg(:) = 5  ! include C3 grass tile ?
end where

! (14) Tundra
where ( pft_biome1(:,1)==11 .and. npft==1 )
  climate%biome(:) = 14
  climate%iveg(:) = 8
end where

! (15) Hot desert
where ( pft_biome1(:,1)==12 .and. npft==1 )
  climate%biome(:) = 15
  climate%iveg(:) = 14
end where

! (16) Semidesert
where ( pft_biome1(:,1)==13 .and. npft==1 )
  climate%biome(:) = 16
  climate%iveg(:) = 5
end where

! (17) Ice/polar desert
where ( climate%biome(:)==999 )
  climate%biome(:) = 17
  climate%iveg(:) = 17
end where

! check for DBL or NEL in SH: set to EBL instead
where ( (climate%iveg(:)==1.or.climate%iveg(:)==3.or.climate%iveg(:)==4) .and. &
        rad%latitude(:)<0._8 )
  climate%iveg(:) = 2
  climate%biome(:) = 4
end where

return
end subroutine biome1_pft

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
real, dimension(ifull) :: svs_sum
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
svs_sum = sum(svs,dim=2)
do n = 1,maxtile
  svs(:,n) = svs(:,n)/svs_sum(:)
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

if ( any(vlinprev>10.) .or. any(vlin>10.) .or. any(vlinnext>10.) .or. any(vlinnext2>10.) ) then
  write(6,*) "ERROR: LAI is out of range"
  write(6,*) "vlinprev,vlin,vlinnext,vlinnext2 ",maxval(vlinprev),maxval(vlin),maxval(vlinnext),maxval(vlinnext2)
  call ccmpi_abort(-1)
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
    !do while ( any(newgrid(:)<minfrac.and.newgrid(:)>0.) )
    !  pos = minloc(newgrid(:), newgrid(:)>0.)
    !  newgrid(pos(1)) = 0.
    !  nsum = sum(newgrid(:))
    !  newgrid(:) = newgrid(:)/nsum
    !end do

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
                   casapoint,greenup,fall,phendoy1,fcasapft)

use carbpools_m
use cc_mpi
use cc_omp, only : ntiles, imax
use const_phys
use infile
use latlong_m
use newmpar_m
use nsibd_m
use parm_m
use pop_constants, only : NPATCH, NLAYER, NCOHORT_MAX, HEIGHT_BINS, NDISTURB, AGEMAX
use pbl_m
use sigs_m
use soil_m
use soilsnow_m
use soilv_m
use vegpar_m
  
implicit none
  
integer, dimension(ifull,maxtile), intent(in) :: ivs
integer, dimension(271,mxvt), intent(in) :: greenup, fall, phendoy1
integer, dimension(:), allocatable, save :: cveg
integer(kind=4), dimension(:), allocatable, save :: Iwood
integer(kind=4), dimension(:,:), allocatable, save :: disturbance_interval
integer i,iq,n,k,ipos,ilat,ivp,is,ie
integer jyear,jmonth,jday,jhour,jmin,mins
integer landcount
integer(kind=4) mp_POP
real ivmax
real, dimension(ifull,maxtile), intent(in) :: svs,vlin,vlinprev,vlinnext,vlinnext2
real, dimension(ifull,5), intent(in) :: casapoint
real, dimension(ifull,2) :: albsoilsn
real, dimension(ifull) :: dummy_pack
real, dimension(ifull) :: albsoil
real, dimension(:), allocatable, save :: dummy_unpack
logical, dimension(:), allocatable, save :: pmap_temp
integer :: tile, popcount
character(len=*), intent(in) :: fcasapft

if ( myid==0 ) write(6,*) "Initialising CABLE"

if ( cbm_ms/=ms ) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  call ccmpi_abort(-1)
end if

POP_NPATCH = NPATCH
POP_NLAYER = NLAYER
POP_NCOHORT = NCOHORT_MAX
POP_HEIGHT_BINS = HEIGHT_BINS
POP_NDISTURB = NDISTURB
POP_AGEMAX = AGEMAX

if ( myid==0 ) write(6,*) "Define CABLE and CASA CNP arrays"

! default values (i.e., no land)  
ivegt = 0
albsoilsn = 0.08  
albsoil   = 0.08
albvisdir = 0.08
albvisdif = 0.08
albnirdir = 0.08
albnirdif = 0.08
zolnd     = 0.
mvtype = mxvt
mstype = mxst

! calculate CABLE vector length
allocate( tdata(ntiles) )
do tile = 1,ntiles
  tdata(tile)%mp = 0
  tdata(tile)%np = 0
  is = 1 + (tile-1)*imax
  ie = tile*imax
  do iq = is,ie
    if ( land(iq) ) then
      landcount = count(svs(iq,:)>0.)
      tdata(tile)%mp = tdata(tile)%mp + landcount
      if ( landcount==0 ) then
        write(6,*) "ERROR: No CABLE tiles assigned to land point: myid,iq,tile",myid,iq,tile
        call ccmpi_abort(-1)
      end if
    end if
  end do
end do
mp_global = tdata(1)%mp
tdata(1)%toffset = 0
tdata(1)%poffset = 0
do tile = 2,ntiles
  mp_global = mp_global + tdata(tile)%mp
  tdata(tile)%toffset = tdata(tile-1)%toffset + tdata(tile-1)%mp
  tdata(tile)%poffset = 0 ! disable for now
end do
mp = 0 ! defined when CABLE model is integrated

ktau_gl = 900
kend_gl = 999

! if CABLE is present on this processor, then start allocating arrays
! Write messages here in case myid==0 has no land-points (mp_global==0)
if ( myid==0 ) then
  write(6,*) "Allocating CABLE and CASA CNP arrays"
  if ( soil_struc==1 ) then 
    write(6,*) "Using SLI soil model"
  end if
  if ( ccycle==0 ) then
    write(6,*) "Using CABLE without carbon cycle"
  else if ( ccycle==1 ) then
    write(6,*) "Using CASA C"
  else if ( ccycle==2 ) then
    write(6,*) "Using CASA CN"
  else
    write(6,*) "Using CASA CNP"
  end if
end if

maxnb = 0

do tile = 1,ntiles
  allocate( tdata(tile)%tind(maxtile,2) )
  tdata(tile)%tind(:,1) = 1
  tdata(tile)%tind(:,2) = 0
  allocate( tdata(tile)%pind(maxtile,2) )
  tdata(tile)%pind(:,1) = 1
  tdata(tile)%pind(:,2) = 0
  tdata(tile)%maxnb = 0
end do

if ( mp_global>0 ) then

   climate%nyear_average = 20
   climate%nday_average = 31
    
  allocate( sv(mp_global) )
  allocate( vl1(mp_global), vl2(mp_global), vl3(mp_global), vl4(mp_global) )
  allocate( cveg(mp_global) )  
  call alloc_cbm_var(air, mp_global)
  call alloc_cbm_var(bgc, mp_global)
  call alloc_cbm_var(canopy, mp_global)
  call alloc_cbm_var(met, mp_global)
  call alloc_cbm_var(bal, mp_global)
  call alloc_cbm_var(rad, mp_global)
  call alloc_cbm_var(rough, mp_global)
  call alloc_cbm_var(soil, mp_global)
  call alloc_cbm_var(ssnow, mp_global)
  call alloc_cbm_var(sum_flux, mp_global)
  call alloc_cbm_var(veg, mp_global)
  if ( cable_climate==1 ) then
    call alloc_cbm_var(climate, mp_global, 24)
    allocate( climate_save%APAR_leaf_sun(mp_global) )
    allocate( climate_save%APAR_leaf_shade(mp_global) )
    allocate( climate_save%Dleaf_sun(mp_global) )
    allocate( climate_save%fwsoil(mp_global) )
    allocate( climate_save%Dleaf_shade(mp_global) )
    allocate( climate_save%Tleaf_sun(mp_global) )
    allocate( climate_save%Tleaf_shade(mp_global) )
    allocate( climate_save%cs_sun(mp_global) )
    allocate( climate_save%cs_shade(mp_global) )
    allocate( climate_save%scalex_sun(mp_global) )
    allocate( climate_save%scalex_shade(mp_global) )
  end if
  allocate( dummy_unpack(mp_global) )
  
  ! Cable configuration
  !cable_user%ssnow_POTEV = "P-M"
  cable_user%ssnow_POTEV = "H-D" ! default Humidity Deficit
  cable_user%MetType = "defa" ! only 4 characters for "default"
  cable_user%diag_soil_resp = "ON"
  cable_user%leaf_respiration = "ON"
  cable_user%run_diag_level = "NONE"
  cable_user%consistency_check = .false.
  cable_user%logworker = .false.
  cable_user%l_new_roughness_soil = .false.
  cable_user%l_rev_corr = .true.
  cable_user%gw_model = cable_gw_model==1
  select case ( cable_climate )
    case(1)
      cable_user%call_climate = .true.
      cable_user%phenology_switch = "climate"
      cable_user%finite_gm = .true.
    case default
      cable_user%call_climate = .false.
      cable_user%phenology_switch = "MODIS"
      cable_user%finite_gm = .false.
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
      cable_user%fwsoil_switch = "Haverd2013"
    case(2)
      cable_user%fwsoil_switch = "Lai and Ktaul 2000"  
    case(1)
      cable_user%fwsoil_switch = "non-linear extrapola" ! only 20 characters for "non-linear extrapolation"
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
  select case ( progvcmax )
    case(2)
      cable_user%vcmax = "Walker2014"        
    case default
      cable_user%vcmax = "standard"    
  end select
  select case ( smrf_switch )
    case(1)
      cable_user%SMRF_NAME = "CASA-CNP"
    case(2)
      cable_user%SMRF_NAME = "SOILN"
    case(3)
      cable_user%SMRF_NAME = "TRIFFID"
    case(4)
      cable_user%SMRF_NAME = "Trudinger2016"
    case(5)
      cable_user%SMRF_NAME = "DAMM"
    case default
      cable_user%SMRF_NAME = "Trudinger2016"
  end select
  select case ( strf_switch )
    case(1)
      cable_user%STRF_NAME = "CASA-CNP"
    case(2)
      cable_user%STRF_NAME = "K1995"
    case(3)
      cable_user%STRF_NAME = "PnET-CN"
    case(4)
      cable_user%STRF_NAME = "LT1994"
    case(5)
      cable_user%STRF_NAME = "DAMM"
    case default
      cable_user%STRF_NAME = "LT1994"
  end select
  kwidth_gl = nint(dt) ! MJT notes - what happens when the timestep is less than a second?
  if ( kwidth_gl==0) then
    write(6,*) "ERROR: Timestep too small for CABLE"
    call ccmpi_abort(-1)
  end if
  
  ! soil parameters
  soil%zse        = real(zse,8) ! soil layer thickness
  soil%zshh(1)    = 0.5_8 * soil%zse(1)
  soil%zshh(ms+1) = 0.5_8 * soil%zse(ms)
  soil%zshh(2:ms) = 0.5_8 * (soil%zse(1:ms-1) + soil%zse(2:ms))
 
  sv = 0.
  vl1 = 0.
  vl2 = 0.
  vl3 = 0.
  vl4 = 0.
  cveg = 0

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  do tile = 1,ntiles
    allocate(tdata(tile)%tmap(imax,maxtile))
    tdata(tile)%tmap = .false.
  end do
  allocate( qmap(ifull,maxtile) )
  qmap = -1 ! missing

  ipos = 0
  do tile = 1,ntiles
    is = 1 + (tile-1)*imax
    ie = tile*imax
    do n = 1,maxtile
      tdata(tile)%tind(n,1) = ipos + 1
      do iq = is,ie
        if ( land(iq) ) then
          if ( svs(iq,n)>0. ) then
            ipos = ipos + 1
            tdata(tile)%tmap(iq-is+1,n) = .true.
            qmap(iq,n) = ipos
          end if
        end if
      end do
      tdata(tile)%tind(n,2) = ipos
    end do
  end do
  
  if ( ipos/=mp_global ) then
    write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
    call ccmpi_abort(-1)
  end if

  do n = 1,maxtile
    call cable_pack(svs(:,n),sv,n)
    call cable_pack(ivs(:,n),cveg,n)
    call cable_pack(isoilm,soil%isoilm,n)
    dummy_pack = max( vlinprev(:,n), 0.01 )
    call cable_pack(dummy_pack,vl1,n)
    dummy_pack = max( vlin(:,n), 0.01 )
    call cable_pack(dummy_pack,vl2,n)
    dummy_pack = max( vlinnext(:,n), 0.01 )
    call cable_pack(dummy_pack,vl3,n)
    dummy_pack = max( vlinnext2(:,n), 0.01 )
    call cable_pack(dummy_pack,vl4,n)
  end do

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

  do tile = 1,ntiles
    tdata(tile)%maxnb = 0
    do n = 1,maxtile
      if ( tdata(tile)%tind(n,2)>=tdata(tile)%tind(n,1) ) then
        tdata(tile)%maxnb = n
      end if
    end do
  end do

  ! calculate actual max tile number
  do tile = 1,ntiles
    maxnb = max( tdata(tile)%maxnb, maxnb )
  end do
  
  call cable_soil_parm(soil)

  ! store bare soil albedo and define snow free albedo
  call cable_pack(albvisnir(:,1),soil%albsoil(:,1))
  call cable_pack(albvisnir(:,2),soil%albsoil(:,2))
  soil%albsoil(:,3) = 0.05_8
    
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

  ! MJT patch
  call cable_pack(albsoil,soil%albsoil(:,1))
  call cable_pack(albsoil,soil%albsoil(:,2))
  call cable_pack(albsoilsn(:,1),ssnow%albsoilsn(:,1)) ! overwritten by CABLE
  call cable_pack(albsoilsn(:,2),ssnow%albsoilsn(:,2)) ! overwritten by CABLE
  call cable_pack(albsoil,rad%albedo_T(:))
  dummy_pack = rlatt*180./pi
  call cable_pack(dummy_pack,rad%latitude(:))
  dummy_pack = rlongg*180./pi
  call cable_pack(dummy_pack,rad%longitude(:))
  
  veg%vcmax_shade = veg%vcmax
  veg%ejmax_shade = veg%ejmax
  veg%vcmax_sun   = veg%vcmax
  veg%ejmax_sun   = veg%ejmax
 
  ssnow%albsoilsn(:,3)=0.05_8    
  ssnow%t_snwlr=0.05_8
  ssnow%pudsmx=0._8
  ssnow%zdelta = 0._8
  ssnow%le = 0._8
  
  canopy%oldcansto=0._8  
  canopy%ghflux=0._8
  canopy%sghflux=0._8
  canopy%ga=0._8
  canopy%dgdtg=0._8
  canopy%fhs_cor=0._8
  canopy%fes_cor=0._8
  canopy%fns_cor=0._8
  canopy%ga_cor=0._8
  canopy%us=0.01_8
  ssnow%wb_lake=0._8 ! not used when mlo.f90 is active
  ssnow%fland=1._8
  ssnow%ifland=soil%isoilm
    
  ! Initialise sum flux variables
  sum_flux%sumpn=0._8
  sum_flux%sumrp=0._8
  sum_flux%sumrs=0._8
  sum_flux%sumrd=0._8
  sum_flux%sumrpw=0._8
  sum_flux%sumrpr=0._8
  sum_flux%dsumpn=0._8
  sum_flux%dsumrp=0._8
  sum_flux%dsumrs=0._8
  sum_flux%dsumrd=0._8
  
  bal%evap_tot=0._8
  bal%precip_tot=0._8
  bal%ebal_tot=0._8
  bal%rnoff_tot=0._8
  
  if ( cable_climate==1 ) then
    climate%doy = 0
    climate%nyears = 0
    climate%biome = 999
    climate%iveg = 0
    climate%chilldays = 0 ! used by phenology
    climate%gdd5 = 0._8   ! used by phenology
    climate%gdd0 = 0._8
    climate%agdd5 = 0._8
    climate%agdd0 = 0._8
    climate%dmoist = 0._8
    climate%dtemp = 0._8
    climate%mtemp = 0._8
    climate%mmoist = 0._8
    climate%qtemp = 0._8
    climate%evap_PT = 0._8
    climate%aevap = 0._8
    climate%mtemp_max = 0._8
    climate%mtemp_min = 0._8
    climate%qtemp_max = 0._8
    climate%qtemp_max_last_year = 0._8 ! used by cable_canopy and CASA-CNP
    climate%dmoist_min = 9999._8
    climate%dmoist_max = -9999._8
    climate%dtemp_91 = 0._8
    climate%dtemp_31 = 0._8
    climate%dmoist_31 = 0._8
    climate%dmoist_min_20 = 0._8
    climate%dmoist_max_20 = 0._8
    climate%dmoist_min20 = 0._8
    climate%dmoist_max20 = 0._8
    climate%atemp_mean = 0._8
    climate%mtemp_min_20 = 0._8 ! missing flag
    climate%mtemp_min20 = 0._8
    climate%mtemp_max_20 = 0._8 ! missing flag
    climate%mtemp_max20 = 0._8
    climate%alpha_PT_20 = 0._8  ! missing flag
    climate%alpha_PT20 = 0._8
    climate%alpha_PT = 0._8
    climate%gdd0_rec = 0._8
    climate%dtemp_min = 0._8
    climate%fdorm = 1._8
    climate%fapar_ann_max = 0._8
    climate%fapar_ann_max_last_year = 0._8
    climate%frec = 1._8
    climate%gmd = 0
    climate%APAR_leaf_sun = 0._8
    climate%APAR_leaf_shade = 0._8
    climate%Dleaf_sun = 0._8
    climate%fwsoil = 0._8
    climate%Dleaf_shade = 0._8
    climate%Tleaf_sun = 0._8
    climate%Tleaf_shade = 0._8
    climate%cs_sun = 0._8
    climate%cs_shade = 0._8
    climate%scalex_sun = 0._8
    climate%scalex_shade = 0._8
    climate_save%APAR_leaf_sun = 0._8
    climate_save%APAR_leaf_shade = 0._8
    climate_save%Dleaf_sun = 0._8
    climate_save%fwsoil = 0._8
    climate_save%Dleaf_shade = 0._8
    climate_save%Tleaf_sun = 0._8
    climate_save%Tleaf_shade = 0._8
    climate_save%cs_sun = 0._8
    climate_save%cs_shade = 0._8
    climate_save%scalex_sun = 0._8
    climate_save%scalex_shade = 0._8
  end if
  
  if ( ccycle==0 ) then
    ! Initialise CABLE carbon pools
    if ( cable_pop==1 ) then
      write(6,*) "ERROR: cable_pop=1 requires ccycle>0"
      call ccmpi_abort(-1)
    end if    
  else if ( ccycle>=1 .and. ccycle<=3 ) then
    ! CASA CNP
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp_global)
    call alloc_phenvariable(phen,mp_global)
    
    casamet%lat = rad%latitude
    
    call cable_pack(casapoint(:,1),casamet%isorder)
    dummy_pack = casapoint(:,2)/365.*1.E-3
    call cable_pack(dummy_pack,casaflux%Nmindep)
    dummy_pack = casapoint(:,3)/365.
    call cable_pack(dummy_pack,casaflux%Nminfix)
    dummy_pack = casapoint(:,4)/365.
    call cable_pack(dummy_pack,casaflux%Pdep)
    dummy_pack = casapoint(:,5)/365.
    call cable_pack(dummy_pack,casaflux%Pwea)

    where ( veg%iveg==9 .or. veg%iveg==10 ) ! crops
      ! P fertilizer =13 Mt P globally in 1994
      casaflux%Pdep = casaflux%Pdep + 0.7_8/365._8
      ! N fertilizer =86 Mt N globally in 1994
      casaflux%Nminfix = casaflux%Nminfix + 4.3_8/365._8
    end where

    if ( any( casamet%isorder<1 .or. casamet%isorder>12 ) ) then
      write(6,*) "ERROR: Invalid isorder in CASA-CNP"
      call ccmpi_abort(-1)
    end if

    call casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)
    
    do n = 1,mp_global
      ilat = nint((rad%latitude(n)+55.25)*2.) + 1
      ilat = min( 271, max( 1, ilat ) )
      ivp = veg%iveg(n)
      phen%phen(n)       = 1._8
      phen%aphen(n)      = 0._8
      phen%phase(n)      = phendoy1(ilat,ivp)
      phen%doyphase(n,1) = greenup(ilat,ivp)          ! DOY for greenup
      phen%doyphase(n,2) = phen%doyphase(n,1) + 14    ! DOY for steady LAI
      phen%doyphase(n,3) = fall(ilat,ivp)             ! DOY for leaf senescence
      phen%doyphase(n,4) = phen%doyphase(n,3) + 14    ! DOY for minimal LAI season
      if ( phen%doyphase(n,2) > 365 ) phen%doyphase(n,2) = phen%doyphase(n,2) - 365
      if ( phen%doyphase(n,4) > 365 ) phen%doyphase(n,4) = phen%doyphase(n,4) - 365
    end do
    
    casamet%tairk         = 0._8
    casamet%tsoil         = 0._8
    casamet%moist         = 0._8
    
    casaflux%cgpp         = 0._8
    casaflux%Crsoil       = 0._8
    casaflux%crgplant     = 0._8
    casaflux%crmplant     = 0._8
    casaflux%clabloss     = 0._8
    casaflux%frac_sapwood = 1._8
    casaflux%sapwood_area = 0._8
    casaflux%stemnpp      = 0._8
    casaflux%Cnpp         = 0._8
    !casaflux%fHarvest     = 0._8
    !casaflux%NHarvest     = 0._8
    !casaflux%CHarvest     = 0._8
    !casaflux%fcrop        = 0._8

    canopy%fnee = 0._8
    canopy%fpn = 0._8
    canopy%frday = 0._8
    canopy%frp = 0._8
    canopy%frpw = 0._8
    canopy%frpr = 0._8
    canopy%frs = 0._8
   
    casabal%LAImax = 0._8
    casabal%Cleafmean = 0._8
    casabal%Crootmean = 0._8
    
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
    do k = 1,mplant
      dummy_unpack = sv*real(casapool%cplant(:,k))  
      call cable_unpack(dummy_unpack,cplant(:,k))
      dummy_unpack = sv*real(casapool%nplant(:,k))
      call cable_unpack(dummy_unpack,niplant(:,k))
      dummy_unpack = sv*real(casapool%pplant(:,k))
      call cable_unpack(dummy_unpack,pplant(:,k))
    end do
    do k = 1,mlitter
      dummy_unpack = sv*real(casapool%clitter(:,k))  
      call cable_unpack(dummy_unpack,clitter(:,k))
      dummy_unpack = sv*real(casapool%nlitter(:,k)) 
      call cable_unpack(dummy_unpack,nilitter(:,k))
      dummy_unpack = sv*real(casapool%plitter(:,k))
      call cable_unpack(dummy_unpack,plitter(:,k))
    end do
    do k = 1,msoil
      dummy_unpack = sv*real(casapool%csoil(:,k))   
      call cable_unpack(dummy_unpack,csoil(:,k))
      dummy_unpack = sv*real(casapool%nsoil(:,k))
      call cable_unpack(dummy_unpack,nisoil(:,k))
      dummy_unpack = sv*real(casapool%psoil(:,k))
      call cable_unpack(dummy_unpack,psoil(:,k))
    end do
    !dummy_unpack = sv*real(casamet%glai(:)) 
    !call cable_unpack(dummy_unpack,glai(:))

    
    ! POP
    if ( cable_pop==1 ) then
      mp_POP = count(casamet%iveg2==forest.or.casamet%iveg2==shrub)
      allocate( pmap_temp(mp_global) )      
      allocate( Iwood(mp_POP), disturbance_interval(mp_POP,2) )

      do tile = 1,ntiles
        allocate( tdata(tile)%pmap(imax,maxtile) )
        tdata(tile)%pmap = .false.
      end do
      pmap_temp(:) = .false.
      ipos = 0
      do tile = 1,ntiles
        popcount = 0
        do n = 1,maxtile
          is = tdata(tile)%tind(n,1)
          ie = tdata(tile)%tind(n,2)
          tdata(tile)%pind(n,1) = ipos + 1
          do i = is,ie
            if ( casamet%iveg2(i)==forest .or. casamet%iveg2(i)==shrub ) then
              ipos = ipos + 1
              popcount = popcount + 1
              Iwood(ipos) = i
              pmap_temp(i) = .true.
            end if
          end do    
          tdata(tile)%pind(n,2) = ipos
          tdata(tile)%pmap(:,n) = unpack(pmap_temp(is:ie),tdata(tile)%tmap(:,n),.false.)
        end do  
        tdata(tile)%np = tdata(tile)%np + popcount
      end do  
      tdata(1)%poffset = 0
      do tile=2,ntiles
        tdata(tile)%poffset=tdata(tile-1)%poffset+tdata(tile-1)%np
      end do

      disturbance_interval(:,:) = veg%disturbance_interval(Iwood(1:mp_POP),:)  

      !Iwood only used inside threaded region from now on, remap so it is relative per tile
      ipos = 0
      do tile = 1,ntiles
        do n = 1,maxtile
          is = tdata(tile)%tind(n,1)
          ie = tdata(tile)%tind(n,2)
          do i = is,ie
            if ( casamet%iveg2(i)==forest .or. casamet%iveg2(i)==shrub ) then
              ipos = ipos + 1
              Iwood(ipos) = i - tdata(tile)%toffset
            end if
          end do    
        end do  
      end do  

      call POP_init(POP, disturbance_interval, mp_POP, Iwood(1:mp_POP)) 
      deallocate( pmap_temp )
      deallocate( Iwood, disturbance_interval )

    end if
    
  else  
    write(6,*) "ERROR: Unknown option ccycle ",ccycle
    call ccmpi_abort(-1)
  end if ! ccycle>0

  
  ! Calculate LAI and veg fraction diagnostics
  ! (needs to occur after CASA-CNP in case prognostic LAI is required)
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jmonth,jday,jhour,jmin,mp_global,sv,vl1,vl2,vl3,vl4,casamet,veg,ifull)
  vlai(:) = 0.
  dummy_unpack(1:mp_global) = sv(1:mp_global)*real(veg%vlai(1:mp_global))
  call cable_unpack(dummy_unpack,vlai)

  ! MJT suggestion to get an approx inital albedo (before cable is called)
  where ( land(1:ifull) )
    albvisnir(:,1) = albsoilsn(:,1)*(1.-sigmf) + 0.03*sigmf
    albvisnir(:,2) = albsoilsn(:,2)*(1.-sigmf) + 0.40*sigmf
  end where
  albvisdir(:) = albvisnir(:,1) ! To be updated by CABLE
  albvisdif(:) = albvisnir(:,1) ! To be updated by CABLE
  albnirdir(:) = albvisnir(:,2) ! To be updated by CABLE
  albnirdif(:) = albvisnir(:,2) ! To be updated by CABLE

  deallocate( dummy_unpack )
  
else

  allocate( cveg(0) )
  call cable_biophysic_parm(cveg)
  deallocate( cveg )
  call cable_soil_parm(soil)
  if ( ccycle>=1 .and. ccycle<=3 ) then
    call casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)
  end if
  
end if
  
if (myid==0) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine cbmparm

subroutine casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)

use cc_mpi     ! CC MPI routines
use newmpar_m

implicit none
type(veg_parameter_type), intent(in) :: veg
type(casa_biome),         intent(inout) :: casabiome
type(casa_pool),          intent(inout) :: casapool
type(casa_flux),          intent(inout) :: casaflux
type(casa_met),           intent(inout) :: casamet
type(phen_variable),      intent(inout) :: phen
character(len=*), intent(in) :: fcasapft

real(kind=r_2), dimension(mxvt,mplant) :: ratiocnplant
real(kind=r_2), dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real(kind=r_2), dimension(mso,msoil) :: rationpsoil

real(kind=r_2), dimension(mxvt) :: leafage,woodage,frootage,metage
real(kind=r_2), dimension(mxvt) :: strage,cwdage,micage,slowage,passage,slax
real(kind=r_2), dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real(kind=r_2), dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real(kind=r_2), dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real(kind=r_2), dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real(kind=r_2), dimension(mxvt) :: cleaf,cwood,cfroot,cmet,cstr,ccwd,cmic,cslow,cpass,nleaf
real(kind=r_2), dimension(mxvt) :: nwood,nfroot,nmet,nstr,ncwd,nmic,nslow,npass,xpleaf,xpwood
real(kind=r_2), dimension(mxvt) :: xpfroot,xpmet,xpstr,xpcwd,xpmic,xpslow,xppass,clabileage
real(kind=r_2), dimension(mxvt) :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase
real(kind=r_2), dimension(mxvt) :: xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
real(kind=r_2), dimension(mxvt) :: xla_to_sa,xvcmax_scalar,xdisturbance_interval
real(kind=r_2), dimension(mxvt) :: xDAMM_EnzPool,xDAMM_KMO2,xDAMM_KMcp,xDAMM_Ea,xDAMM_alpha
real(kind=r_2), dimension(mso) :: xxkplab,xxkpsorb,xxkpocc
real(kind=r_2), dimension(mso) :: xkmlabp,xpsorbmax,xfPleach

integer :: i, iso, nv, ierr
integer :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12,nv13
integer :: fflag=0

integer, dimension(mxvt) :: ivt2
real(kind=r_2), dimension(mxvt) :: kroot
real(kind=r_2), dimension(mxvt) :: rootdepth
real(kind=r_2), dimension(mxvt) :: kuptake
real(kind=r_2), dimension(mxvt) :: krootlen
real(kind=r_2), dimension(mxvt) :: kminn
real(kind=r_2), dimension(mxvt) :: kuplabp
real(kind=r_2), dimension(mxvt,mplant) :: fracnpptop
real(kind=r_2), dimension(mxvt,mplant) :: rmplant
real(kind=r_2), dimension(mxvt,mplant) :: ftransnptol
real(kind=r_2), dimension(mxvt,mplant) :: fracligninplant
real(kind=r_2), dimension(mxvt) :: glaimax
real(kind=r_2), dimension(mxvt) :: glaimin
real(kind=r_2), dimension(mxvt) :: xkleafcoldexp
real(kind=r_2), dimension(mxvt) :: xkleafdryexp
real(kind=r_2), dimension(mxvt,mplant) :: rationcplantmin
real(kind=r_2), dimension(mxvt,mplant) :: rationcplantmax
real(kind=r_2), dimension(mxvt,mplant) :: ftranspptol
real(kind=r_2), dimension(mxvt) :: tkshed

if ( trim(fcasapft) /= '' ) fflag = 1
call ccmpi_bcast(fflag,0,comm_world)

if ( fflag == 1  ) then
  if ( myid == 0 ) then
    write(6,*) "Using user defined CASA PFT parameter tables"
  end if

  if ( myid == 0 ) then
    open(86,file=fcasapft,status='old',action='read',iostat=ierr)
  end if
  call ccmpi_bcast(ierr,0,comm_world)
  if ( ierr/=0 ) then
    if ( myid == 0 ) then
      write(6,*) "ERROR: Cannot open casapftfile=",trim(fcasapft)
    end if
    call ccmpi_abort(-1)
  end if

  if ( myid == 0 ) then
    write(6,*) "Reading ",trim(fcasapft)

    do i=1,3
      read(86,*)
    enddo

    do nv=1,mxvt
      read(86,*,iostat=ierr) nv0,ivt2(nv)
      if ( ierr/=0 .or. nv0/=nv  ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+0*(2+mxvt)+nv,"reading nv0"
        if ( nv0/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv1,kroot(nv),rootdepth(nv),      &
                  kuptake(nv),krootlen(nv),         &
                  kminn(nv), kuplabp(nv),           &
                  xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                  metage(nv),strage(nv),cwdage(nv),  &
                  micage(nv),slowage(nv),passage(nv),clabileage(nv),slax(nv)
      if ( ierr/=0 .or. nv1/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+1*(2+mxvt)+nv,"reading nv1"
        if ( nv1/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv2, &
                  fracnpptop(nv,leaf),fracnpptop(nv,wood), &
                  fracnpptop(nv,xroot),rmplant(nv,leaf),   &
                  rmplant(nv,wood),rmplant(nv,xroot)
      if ( ierr/=0 .or. nv2/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+2*(2+mxvt)+nv,"reading nv2"
        if ( nv2/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv3, ratiocnplant(nv,leaf),ratiocnplant(nv,wood),   &
           ratiocnplant(nv,xroot),                                         &
           ftransnptol(nv,leaf), ftransnptol(nv,wood), &
           ftransnptol(nv,xroot),                                &
           fracligninplant(nv,leaf),                             &
           fracligninplant(nv,wood),                             &
           fracligninplant(nv,xroot),                            &
           ratiocnsoil(nv,mic),ratiocnsoil(nv,slow),ratiocnsoil(nv,pass),  &
           ratiocnsoilmin(nv,mic),ratiocnsoilmin(nv,slow),ratiocnsoilmin(nv,pass),  &
           ratiocnsoilmax(nv,mic),ratiocnsoilmax(nv,slow),ratiocnsoilmax(nv,pass),  &
           glaimax(nv),glaimin(nv)
      if ( ierr/=0 .or. nv3/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+3*(2+mxvt)+nv,"reading nv3"
        if ( nv3/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv4, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
                  cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
      if ( ierr/=0 .or. nv4/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+4*(2+mxvt)+nv,"reading nv4"
        if ( nv4/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv5, &
           tkshed(nv),xxkleafcoldmax(nv),xkleafcoldexp(nv), &
           xxkleafdrymax(nv),xkleafdryexp(nv)
      if ( ierr/=0 .or. nv5/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+5*(2+mxvt)+nv,"reading nv5"
        if ( nv5/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv6, &
        rationcplantmin(nv,leaf),rationcplantmax(nv,leaf), &
        rationcplantmin(nv,wood),rationcplantmax(nv,wood), &
        rationcplantmin(nv,xroot),rationcplantmax(nv,xroot), &
        xfnminloss(nv), xfnminleach(nv),xnfixrate(nv)
      if ( ierr/=0 .or. nv6/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+6*(2+mxvt)+nv,"reading nv6"
        if ( nv6/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                  nmet(nv),nstr(nv), ncwd(nv), &
                  nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
      if ( ierr/=0 .or. nv7/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+7*(2+mxvt)+nv,"reading nv7"
        if ( nv7/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv8,xrationpleafmin(nv),xrationpleafmax(nv),    &
           xrationpwoodmin(nv),xrationpwoodmax(nv),              &
           xrationpfrootmin(nv),xrationpfrootmax(nv),            &
           ftranspptol(nv,leaf), ftranspptol(nv,wood), &
           ftranspptol(nv,xroot)
      if ( ierr/=0 .or. nv8/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+8*(2+mxvt)+nv,"reading nv8"
        if ( nv8/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do iso=1,mso
      read(86,*,iostat=ierr) nv9,xkmlabp(iso),xpsorbmax(iso),xfpleach(iso), &
                  rationpsoil(iso,mic),rationpsoil(iso,slow),rationpsoil(iso,pass), &
                  xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)
      if ( ierr/=0 .or. nv9/=iso ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+9*(2+mxvt)+nv,"reading nv9"
        if ( nv9/=iso ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo
    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv10, &
           xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
           xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
      if ( ierr/=0 .or. nv10/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+10*(2+mxvt)+nv,"reading nv10"
        if ( nv10/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv11, &
           xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
           xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)
      if ( ierr/=0 .or. nv11/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+11*(2+mxvt)+nv,"reading nv11"
        if ( nv11/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv12, &
           xla_to_sa(nv),xdisturbance_interval(nv),xvcmax_scalar(nv)
      if ( ierr/=0 .or. nv12/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+12*(2+mxvt)+nv,"reading nv12"
        if ( nv12/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    read(86,*)
    read(86,*)
    do nv=1,mxvt
      read(86,*,iostat=ierr) nv13, &
           xdamm_enzpool(nv), xdamm_kmo2(nv),xdamm_kmcp(nv), xdamm_ea(nv), xdamm_alpha(nv)
      if ( ierr/=0 .or. nv13/=nv ) then
        write(6,*) "ERROR: Cannot read casapftfile file ",trim(fcasapft)
        write(6,*) "Formatting error in line",3+13*(2+mxvt)+nv,"reading nv13"
        if ( nv13/=nv ) then
          write(6,*) "pft index is not sequential"
        end if
        call ccmpi_abort(-1)
       end if
    enddo

    close(86)
  end if

  call ccmpi_bcast(ivt2,0,comm_world)
  call ccmpi_bcastr8(kroot,0,comm_world)
  call ccmpi_bcastr8(rootdepth,0,comm_world)
  call ccmpi_bcastr8(kuptake,0,comm_world)
  call ccmpi_bcastr8(krootlen,0,comm_world)
  call ccmpi_bcastr8(kminn,0,comm_world)
  call ccmpi_bcastr8(kuplabp,0,comm_world)
  call ccmpi_bcastr8(xfherbivore,0,comm_world)
  call ccmpi_bcastr8(leafage,0,comm_world)
  call ccmpi_bcastr8(woodage,0,comm_world)
  call ccmpi_bcastr8(frootage,0,comm_world)
  call ccmpi_bcastr8(metage,0,comm_world)
  call ccmpi_bcastr8(strage,0,comm_world)
  call ccmpi_bcastr8(cwdage,0,comm_world)
  call ccmpi_bcastr8(micage,0,comm_world)
  call ccmpi_bcastr8(slowage,0,comm_world)
  call ccmpi_bcastr8(passage,0,comm_world)
  call ccmpi_bcastr8(clabileage,0,comm_world)
  call ccmpi_bcastr8(slax,0,comm_world)

  call ccmpi_bcastr8(fracnpptop,0,comm_world)
  call ccmpi_bcastr8(rmplant,0,comm_world)

  call ccmpi_bcastr8(ratiocnplant,0,comm_world)
  call ccmpi_bcastr8(ftransnptol,0,comm_world)
  call ccmpi_bcastr8(fracligninplant,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoil,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoilmin,0,comm_world)
  call ccmpi_bcastr8(ratiocnsoilmax,0,comm_world)
  call ccmpi_bcastr8(glaimax,0,comm_world)
  call ccmpi_bcastr8(glaimin,0,comm_world)

  call ccmpi_bcastr8(cleaf,0,comm_world)
  call ccmpi_bcastr8(cwood,0,comm_world)
  call ccmpi_bcastr8(cfroot,0,comm_world)
  call ccmpi_bcastr8(cmet,0,comm_world)
  call ccmpi_bcastr8(cstr,0,comm_world)
  call ccmpi_bcastr8(ccwd,0,comm_world)
  call ccmpi_bcastr8(cmic,0,comm_world)
  call ccmpi_bcastr8(cslow,0,comm_world)
  call ccmpi_bcastr8(cpass,0,comm_world)

  call ccmpi_bcastr8(tkshed,0,comm_world)
  call ccmpi_bcastr8(xxkleafcoldmax,0,comm_world)
  call ccmpi_bcastr8(xkleafcoldexp,0,comm_world)
  call ccmpi_bcastr8(xxkleafdrymax,0,comm_world)
  call ccmpi_bcastr8(xkleafdryexp,0,comm_world)

  call ccmpi_bcastr8(rationcplantmin,0,comm_world)
  call ccmpi_bcastr8(rationcplantmax,0,comm_world)
  call ccmpi_bcastr8(xfnminloss,0,comm_world)
  call ccmpi_bcastr8(xfnminleach,0,comm_world)
  call ccmpi_bcastr8(xnfixrate,0,comm_world)

  call ccmpi_bcastr8(nleaf,0,comm_world)
  call ccmpi_bcastr8(nwood,0,comm_world)
  call ccmpi_bcastr8(nfroot,0,comm_world)
  call ccmpi_bcastr8(nmet,0,comm_world)
  call ccmpi_bcastr8(nstr,0,comm_world)
  call ccmpi_bcastr8(ncwd,0,comm_world)
  call ccmpi_bcastr8(nmic,0,comm_world)
  call ccmpi_bcastr8(nslow,0,comm_world)
  call ccmpi_bcastr8(npass,0,comm_world)
  call ccmpi_bcastr8(xnsoilmin,0,comm_world)

  call ccmpi_bcastr8(xrationpleafmin,0,comm_world)
  call ccmpi_bcastr8(xrationpleafmax,0,comm_world)
  call ccmpi_bcastr8(xrationpwoodmin,0,comm_world)
  call ccmpi_bcastr8(xrationpwoodmax,0,comm_world)
  call ccmpi_bcastr8(xrationpfrootmin,0,comm_world)
  call ccmpi_bcastr8(xrationpfrootmax,0,comm_world)
  call ccmpi_bcastr8(ftranspptol,0,comm_world)

  call ccmpi_bcastr8(xkmlabp,0,comm_world)
  call ccmpi_bcastr8(xpsorbmax,0,comm_world)
  call ccmpi_bcastr8(xfpleach,0,comm_world)
  call ccmpi_bcastr8(rationpsoil,0,comm_world)
  call ccmpi_bcastr8(xxkplab,0,comm_world)
  call ccmpi_bcastr8(xxkpsorb,0,comm_world)
  call ccmpi_bcastr8(xxkpocc,0,comm_world)

  call ccmpi_bcastr8(xpleaf,0,comm_world)
  call ccmpi_bcastr8(xpwood,0,comm_world)
  call ccmpi_bcastr8(xpfroot,0,comm_world)
  call ccmpi_bcastr8(xpmet,0,comm_world)
  call ccmpi_bcastr8(xpstr,0,comm_world)
  call ccmpi_bcastr8(xpcwd,0,comm_world)
  call ccmpi_bcastr8(xpmic,0,comm_world)
  call ccmpi_bcastr8(xpslow,0,comm_world)
  call ccmpi_bcastr8(xppass,0,comm_world)
  call ccmpi_bcastr8(xplab,0,comm_world)
  call ccmpi_bcastr8(xpsorb,0,comm_world)
  call ccmpi_bcastr8(xpocc,0,comm_world)

  call ccmpi_bcastr8(xxnpmax,0,comm_world)
  call ccmpi_bcastr8(xq10soil,0,comm_world)
  call ccmpi_bcastr8(xxkoptlitter,0,comm_world)
  call ccmpi_bcastr8(xxkoptsoil,0,comm_world)
  call ccmpi_bcastr8(xprodptase,0,comm_world)
  call ccmpi_bcastr8(xcostnpup,0,comm_world)
  call ccmpi_bcastr8(xmaxfinelitter,0,comm_world)
  call ccmpi_bcastr8(xmaxcwd,0,comm_world)
  call ccmpi_bcastr8(xnintercept,0,comm_world)
  call ccmpi_bcastr8(xnslope,0,comm_world)

  call ccmpi_bcastr8(xla_to_sa,0,comm_world)
  call ccmpi_bcastr8(xdisturbance_interval,0,comm_world)
  call ccmpi_bcastr8(xvcmax_scalar,0,comm_world)

  call ccmpi_bcastr8(xdamm_enzpool,0,comm_world)
  call ccmpi_bcastr8(xdamm_kmo2,0,comm_world)
  call ccmpi_bcastr8(xdamm_kmcp,0,comm_world)
  call ccmpi_bcastr8(xdamm_ea,0,comm_world)
  call ccmpi_bcastr8(xdamm_alpha,0,comm_world)
  if ( mp_global>0 ) then
    casabiome%ivt2 = ivt2
    casabiome%kroot = kroot
    casabiome%rootdepth = rootdepth
    casabiome%kuptake = kuptake
    casabiome%krootlen = krootlen
    casabiome%kminn = kminn
    casabiome%kuplabp = kuplabp
    casabiome%fracnpptop = fracnpptop
    casabiome%rmplant = rmplant
    casabiome%ftransnptol = ftransnptol
    casabiome%fracligninplant = fracligninplant
    casabiome%glaimax = glaimax
    casabiome%glaimin = glaimin
    phen%tkshed = tkshed
    casabiome%xkleafcoldexp = xkleafcoldexp
    casabiome%xkleafdryexp = xkleafdryexp
    casabiome%rationcplantmin = rationcplantmin
    casabiome%rationcplantmax = rationcplantmax
    casabiome%ftranspptol = ftranspptol
  end if
else
  if ( myid == 0 ) then
    write(6,*) "Using default CASA PFT parameter tables"
  end if
  if ( mp_global>0 ) then

    leafage =(/ 2.0_r_2, 1.5_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2, 0.8_r_2, 0.8_r_2, 1.0_r_2,      0.8_r_2,      0.8_r_2, 1.0_r_2, &
                1.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2 /)
    woodage =(/ 70._r_2, 60._r_2, 80._r_2, 40._r_2, 40._r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2,      1.0_r_2,      1.0_r_2, 1.0_r_2, &
                1.0_r_2, 1.0_r_2, 5.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2 /)
    frootage=(/ 18._r_2, 10._r_2, 10._r_2, 10._r_2, 5.0_r_2, 3.0_r_2, 3.0_r_2, 3.0_r_2, 0.884227_r_2, 0.884227_r_2, 1.0_r_2, &
                1.0_r_2, 1.0_r_2, 4.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2 /)
    metage=0.04_r_2
    strage=0.23_r_2
    cwdage=0.824_r_2
    micage=0.137_r_2
    slowage=5._r_2
    passage=222.22_r_2
    clabileage=0.2_r_2
    slax = (/ 0.007178742_r_2, 0.015319264_r_2, 0.023103346_r_2, 0.026464564_r_2, 0.009919764_r_2, 0.029590494_r_2, &
              0.022417511_r_2, 0.026704118_r_2, 0.029590494_r_2, 0.022417511_r_2,        0.02_r_2,        0.02_r_2, &
              0.02_r_2,        0.024470894_r_2,        0.02_r_2,        0.02_r_2,        0.02_r_2 /)

    xfherbivore   =(/ 0.068_r_2, 0.406_r_2, 0.068_r_2, 0.134_r_2, 0.022_r_2, 0.109_r_2, 0.109_r_2, 0.109_r_2, 0.140_r_2, &
                      0.140_r_2, 0.000_r_2, 0.000_r_2, 0.000_r_2, 0.010_r_2, 0.000_r_2, 0.000_r_2, 0.000_r_2 /)
    xxkleafcoldmax=(/   0.2_r_2,  0.1_r_2,  0.1_r_2,  0.6_r_2,       1._r_2,  0.2_r_2,    0.2_r_2,   0.2_r_2,   0.3_r_2, &
                        0.3_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,      0.1_r_2,  0.1_r_2,    0.1_r_2,   0.1_r_2 /)
    xxkleafdrymax =(/   0.1_r_2,  0.1_r_2,  0.1_r_2,   1._r_2,      0.1_r_2,  0.1_r_2,    0.1_r_2,   0.1_r_2,   0.1_r_2, &
                        0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,      0.1_r_2,  0.1_r_2,    0.1_r_2,   0.1_r_2 /)
    xratioNPleafmin =(/ 10.92308_r_2, 15.95339_r_2, 9.254839_r_2, 12.73848_r_2, 12.07217_r_2, 13.51473_r_2,    14.05_r_2, &
                        12.57800_r_2, 15.12262_r_2,      10._r_2,      13._r_2,      10._r_2,      10._r_2,  16.2336_r_2, &
                             10._r_2,      10._r_2,      10._r_2 /)
    xratioNPleafmax =(/ 12.07288_r_2,  17.6327_r_2, 10.22903_r_2, 14.07938_r_2, 13.34292_r_2, 14.93733_r_2, 15.52895_r_2, &
                          13.902_r_2, 16.71447_r_2,      10._r_2,      13._r_2,      10._r_2,      10._r_2,  17.9424_r_2, &
                             10._r_2,      10._r_2,      10._r_2 /)
    xratioNPwoodmin =(/ 20.30167_r_2, 15.89425_r_2, 17.48344_r_2, 19.08018_r_2, 22.46035_r_2,      15._r_2,      15._r_2, &
                           15.96_r_2,    20.52_r_2,      15._r_2,      15._r_2,      15._r_2,      15._r_2,  17.5275_r_2, &
                             15._r_2,      15._r_2,      15._r_2 /)
    xratioNPwoodmax =(/ 22.43869_r_2, 17.56733_r_2,  19.3238_r_2, 21.08862_r_2,  24.8246_r_2,      15._r_2,      15._r_2, &
                           17.64_r_2,    20.52_r_2,      15._r_2,      15._r_2,      15._r_2,      15._r_2,  19.3725_r_2, &
                             15._r_2,      15._r_2,      15._r_2 /)
    xratioNPfrootmin=(/ 20.29341_r_2, 15.87155_r_2, 17.39767_r_2,  19.0601_r_2, 22.49363_r_2, 15.63498_r_2, 16.08255_r_2, &
                        14.49241_r_2, 22.69109_r_2,      15._r_2,      15._r_2,      15._r_2,      15._r_2, 22.13268_r_2, &
                             15._r_2,      15._r_2,      15._r_2 /)
    xratioNPfrootmax=(/ 22.42955_r_2, 17.54224_r_2,   19.229_r_2, 21.06643_r_2, 24.86138_r_2, 17.28077_r_2, 17.77545_r_2, &
                        16.01793_r_2, 25.07962_r_2,      15._r_2,      15._r_2,      15._r_2,      15._r_2, 24.46244_r_2, &
                             15._r_2,      15._r_2,      15._r_2 /)
    xfNminloss=0.05_r_2
    xfNminleach=0.05_r_2
    xnfixrate=(/ 0.08_r_2, 2.6_r_2, 0.21_r_2, 1.64_r_2, 0.37_r_2, 0.95_r_2, 0.95_r_2, 0.95_r_2, 4._r_2, 4._r_2, 0._r_2, &
                   0._r_2,  0._r_2, 0.35_r_2,   0._r_2,   0._r_2,   0._r_2 /)
    xnsoilmin=1000._r_2
  
    ratiocnplant(:,leaf)=(/  49.8_r_2,  23.1_r_2,  59.3_r_2,  31.4_r_2,  37.6_r_2, 34.8_r_2,  44._r_2,  49.2_r_2,  21.6_r_2, &
                              25._r_2,   30._r_2,   30._r_2,   30._r_2,   50._r_2,  40._r_2,  40._r_2,   40._r_2 /)
    ratiocnplant(:,wood)=(/ 238.1_r_2, 134.9_r_2, 243.8_r_2, 156.2_r_2, 142.1_r_2, 150._r_2, 150._r_2, 147.3_r_2,  150._r_2, &
                             125._r_2,  150._r_2,  150._r_2,  150._r_2,  150._r_2, 150._r_2, 135._r_2,  150._r_2 /)
    ratiocnplant(:,xroot)=(/ 73.7_r_2,  61.2_r_2,   75._r_2,  63.2_r_2,  67.1_r_2, 64.5_r_2, 62.7_r_2,   69._r_2,  60.7_r_2, &
                              71._r_2,   71._r_2,   71._r_2,   71._r_2,   71._r_2,  71._r_2,  71._r_2,   71._r_2 /)
    ratiocnsoil(:,mic)=8._r_2
    ratiocnsoil(:,slow)=(/ 16.1_r_2, 12.8_r_2, 24.8_r_2,  30._r_2, 19.3_r_2, 13.1_r_2, 13.1_r_2, 13.1_r_2, 13.2_r_2, 13.2_r_2, &
                           13.1_r_2, 13.1_r_2, 13.1_r_2, 26.8_r_2,  20._r_2,  20._r_2,  20._r_2 /)
    ratiocnsoil(:,pass)=(/ 16.1_r_2, 12.8_r_2, 24.8_r_2,  30._r_2, 19.3_r_2, 13.1_r_2, 13.1_r_2, 13.1_r_2, 13.2_r_2, 13.2_r_2, &
                           13.1_r_2, 13.1_r_2, 13.1_r_2, 26.8_r_2,  20._r_2,  20._r_2,  20._r_2 /)
    ratiocnsoilmin(:,mic)=3._r_2
    ratiocnsoilmin(:,slow)=12._r_2
    ratiocnsoilmin(:,pass)=7._r_2
    ratiocnsoilmax(:,mic)=15._r_2
    ratiocnsoilmax(:,slow)=30._r_2
    ratiocnsoilmax(:,pass)=15._r_2
   
    ! Initial values for CNP pools over 3*plant, 3*litter and 3*soil (=27 pools in total)
    cleaf  =(/ 384.6037_r_2,     273._r_2, 96.59814_r_2, 150.2638_r_2,      88._r_2, 137.1714_r_2, 137.1714_r_2, 137.1714_r_2, &
                   160._r_2,     160._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cwood  =(/ 7865.396_r_2,   11451._r_2, 5683.402_r_2, 10833.74_r_2,     372._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cfroot =(/     250._r_2,    2586._r_2,     220._r_2,     220._r_2,     140._r_2,     263._r_2,     263._r_2,     263._r_2, &
                   240._r_2,     240._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cmet   =(/ 6.577021_r_2, 44.63457_r_2, 7.127119_r_2, 10.97797_r_2, 3.229374_r_2, 28.57245_r_2, 28.57245_r_2, 28.57245_r_2, &
               28.57245_r_2, 28.57245_r_2,       0._r_2,       0._r_2,       0._r_2, 1.457746_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cstr   =(/ 209.1728_r_2, 433.7626_r_2, 277.7733_r_2, 312.5492_r_2, 39.44449_r_2, 50.91091_r_2, 50.91091_r_2, 50.91091_r_2, &
               50.91091_r_2, 50.91091_r_2,       0._r_2,       0._r_2,       0._r_2, 4.956338_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    ccwd   =(/ 606.0255_r_2, 1150.765_r_2, 776.7331_r_2, 888.5864_r_2, 111.5864_r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, 28.44085_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /) 
    cmic   =(/  528.664_r_2, 11.37765_r_2, 597.0785_r_2, 405.5554_r_2, 168.0451_r_2, 425.6431_r_2, 425.6431_r_2, 425.6431_r_2, &
               512.4247_r_2, 512.4247_r_2,       0._r_2,       0._r_2,       0._r_2, 57.77585_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cslow  =(/ 13795.94_r_2, 311.8092_r_2, 16121.12_r_2, 11153.25_r_2, 4465.478_r_2, 5694.437_r_2, 5694.437_r_2, 5694.437_r_2, &
               6855.438_r_2, 6855.438_r_2,       0._r_2,       0._r_2,       0._r_2, 1325.052_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    cpass  =(/ 4425.396_r_2 ,13201.81_r_2, 5081.802_r_2, 5041.192_r_2, 1386.477_r_2,  4179.92_r_2,  4179.92_r_2,  4179.92_r_2, &
               5032.137_r_2, 5032.137_r_2,       0._r_2,       0._r_2,       0._r_2, 517.1719_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nleaf  =(/ 7.541249_r_2,      9.9_r_2, 1.609969_r_2, 3.756594_r_2, 2.933333_r_2, 4.572381_r_2, 4.572381_r_2, 4.572381_r_2, &
               5.333333_r_2, 5.333333_r_2       ,0._r_2,       0._r_2,       0._r_2,      0.5_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nwood  =(/ 31.46159_r_2,    102._r_2,  22.73361_r_2, 80.24989_r_2, 2.755555_r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,      0._r_2,        0._r_2,       0._r_2,       0._r_2, 0.125926_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nfroot =(/ 6.097561_r_2,      38._r_2, 5.365854_r_2, 5.365854_r_2, 3.414634_r_2, 6.414634_r_2, 6.414634_r_2, 6.414634_r_2, &
               5.853659_r_2, 5.853659_r_2,       0._r_2,       0._r_2,       0._r_2, 1.536585_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nmet   =(/ 0.064481_r_2,  0.74391_r_2, 0.059393_r_2, 0.137225_r_2, 0.053823_r_2, 0.476208_r_2, 0.476208_r_2, 0.476208_r_2, &
               0.476208_r_2, 0.476208_r_2,       0._r_2,       0._r_2,       0._r_2, 0.018222_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nstr   =(/ 1.394485_r_2, 2.891751_r_2, 1.851822_r_2, 2.083661_r_2, 0.262963_r_2, 0.339406_r_2, 0.339406_r_2, 0.339406_r_2, &
               0.339406_r_2, 0.339406_r_2,       0._r_2,       0._r_2,       0._r_2, 0.033042_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    ncwd   =(/ 2.424102_r_2, 8.524183_r_2, 3.106932_r_2, 6.581996_r_2, 0.826566_r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, 0.210673_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    nmic   =(/  52.8664_r_2, 1.137765_r_2, 59.70785_r_2, 40.55554_r_2, 16.80451_r_2, 42.56431_r_2, 42.56431_r_2, 42.56431_r_2, &
               51.24247_r_2, 51.24247_r_2,       0._r_2,       0._r_2,       0._r_2, 5.777585_r_2       ,0._r_2,       0._r_2, &
                     0._r_2 /)
    nslow  =(/ 919.7293_r_2, 20.78728_r_2, 1074.741_r_2, 743.5501_r_2, 297.6985_r_2, 379.6291_r_2, 379.6291_r_2, 379.6291_r_2, &
               457.0292_r_2, 457.0292_r_2,       0._r_2,       0._r_2,       0._r_2, 88.33682_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    npass  =(/ 295.0264_r_2, 880.1209_r_2, 338.7868_r_2, 336.0795_r_2,  92.4318_r_2, 278.6613_r_2, 278.6613_r_2, 278.6613_r_2, &
               335.4758_r_2, 335.4758_r_2,       0._r_2,       0._r_2,       0._r_2, 34.47813_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpleaf =(/ 0.191648_r_2,   0.415_r_2,  0.115988_r_2, 0.135453_r_2, 0.022821_r_2,  0.15125_r_2,  0.15125_r_2,  0.15125_r_2, &
                0.15125_r_2, 0.15125_r_2,        0._r_2,       0._r_2,       0._r_2,    0.007_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpwood =(/ 0.953979_r_2,    5.88_r_2, 0.64438_r_2,   2.424778_r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,      0._r_2,      0._r_2,         0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpfroot=(/ 0.076659_r_2,    1.95_r_2, 0.080548_r_2,  0.141097_r_2, 0.037083_r_2,  0.15125_r_2,  0.15125_r_2,  0.15125_r_2, &
                0.15125_r_2, 0.15125_r_2,       0._r_2,        0._r_2,       0._r_2,  0.00875_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpmet  =(/ 0.004385_r_2, 0.029756_r_2, 0.004751_r_2, 0.007319_r_2, 0.002153_r_2, 0.019048_r_2, 0.019048_r_2, 0.019048_r_2, &
               0.019048_r_2, 0.019048_r_2,       0._r_2,       0._r_2,       0._r_2, 0.000972_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpstr  =(/ 0.069724_r_2, 0.144588_r_2, 0.092591_r_2, 0.104183_r_2, 0.013148_r_2,  0.01697_r_2,  0.01697_r_2,  0.01697_r_2, &
                0.01697_r_2,  0.01697_r_2,       0._r_2,       0._r_2,       0._r_2, 0.001652_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpcwd  =(/ 0.101004_r_2, 0.191794_r_2, 0.129456_r_2, 0.148095_r_2, 0.018598_r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpmic  =(/ 6.872632_r_2 , 0.14791_r_2, 7.762021_r_2,  5.27222_r_2, 2.184586_r_2, 5.533361_r_2, 5.533361_r_2, 5.533361_r_2, &
               6.661522_r_2, 6.661522_r_2,       0._r_2,       0._r_2,       0._r_2, 0.751086_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xpslow =(/ 119.5648_r_2, 2.702347_r_2, 139.7164_r_2, 96.66152_r_2, 38.70081_r_2, 49.35178_r_2, 49.35178_r_2, 49.35178_r_2, &
                59.4138_r_2,  59.4138_r_2,       0._r_2,       0._r_2,       0._r_2, 11.48379_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xppass =(/ 38.35343_r_2, 114.4157_r_2, 44.04228_r_2, 43.69033_r_2, 12.01613_r_2, 36.22598_r_2, 36.22598_r_2, 36.22598_r_2, &
               43.61185_r_2, 43.61185_r_2,       0._r_2,       0._r_2,       0._r_2, 4.482157_r_2,       0._r_2,       0._r_2, &
                     0._r_2 /)
    xplab  =(/   26.737_r_2,   19.947_r_2,   29.107_r_2,   30.509_r_2,   23.206_r_2,   25.538_r_2,   25.538_r_2,   25.538_r_2, &
                 27.729_r_2,   27.729_r_2,       0._r_2,       0._r_2,       0._r_2,   21.038_r_2,       0._r_2,       0._r_2, &
                  0.103_r_2 /)
    xpsorb =(/   126.73_r_2,   92.263_r_2,  134.639_r_2,  132.012_r_2,   173.47_r_2,  186.207_r_2,  186.207_r_2,  186.207_r_2, &
                155.518_r_2,  155.518_r_2,       0._r_2,       0._r_2,       0._r_2,   255.79_r_2,       0._r_2,       0._r_2, &
                  1.176_r_2 /)
    xpocc  =(/  138.571_r_2,  120.374_r_2,   138.22_r_2,  148.083_r_2,  114.496_r_2,  145.163_r_2,  145.163_r_2,  145.163_r_2, &
                158.884_r_2,  158.884_r_2,       0._r_2,       0._r_2,       0._r_2,  108.897_r_2,       0._r_2,       0._r_2, &
                  0.688_r_2 /)

    xkmlabp  =(/ 74.5408_r_2,  68.1584_r_2,   77.952_r_2, 64.41918_r_2, 64.41918_r_2, 70.5856_r_2,  64.5888_r_2, 54.1692_r_2, &
                  9.7704_r_2,    28.29_r_2,   63.963_r_2,   32.402_r_2 /)
    xpsorbmax=(/ 745.408_r_2, 788.0815_r_2, 1110.816_r_2,  744.847_r_2,  744.847_r_2, 816.146_r_2, 746.8081_r_2, 722.256_r_2, &
                 293.112_r_2,   311.19_r_2, 373.1175_r_2, 615.6381_r_2 /)
    xfPleach =0.0005_r_2
    ratioNPsoil(:,mic)=4._r_2
    ratioNPsoil(:,slow)=(/ 5._r_2, 5._r_2, 5._r_2, 15._r_2, 5._r_2, 5._r_2, 5._r_2, 5._r_2, 7._r_2, 7._r_2, 7._r_2, 7._r_2 /)
    ratioNPsoil(:,pass)=(/ 5._r_2, 5._r_2, 5._r_2, 15._r_2, 5._r_2, 5._r_2, 5._r_2, 5._r_2, 7._r_2, 7._r_2, 7._r_2, 7._r_2 /)
  
    xxnpmax = (/ 1.510856726_r_2, 1.27916225_r_2, 1.591076159_r_2, 1.186066584_r_2, 1.358075681_r_2,  1.45621905_r_2, &
                  1.45621905_r_2, 1.45621905_r_2, 1.210382326_r_2, 1.210382326_r_2,  1.45621905_r_2, 1.365993164_r_2, &
                 1.210382326_r_2,         1._r_2, 1.399652677_r_2,          1._r_2,          1._r_2 /)
    xq10soil = 1.72_r_2
    xxkoptlitter = 0.4_r_2
    xxkoptsoil = (/ 0.33_r_2, 0.6_r_2, 0.15_r_2, 0.6_r_2, 0.16_r_2, 0.4_r_2, 0.3_r_2, 0.2_r_2, 0.2_r_2, 0.25_r_2,  1._r_2, &
                    0.65_r_2, 0.5_r_2,   2._r_2, 0.5_r_2,   1._r_2,  1._r_2 /)
    xprodptase = (/  0.5_r_2, 0.2_r_2,  0.5_r_2, 0.5_r_2,  0.5_r_2, 0.5_r_2, 0.5_r_2, 0.5_r_2, 0.5_r_2,  0.5_r_2, 0.5_r_2, &
                      4._r_2, 0.5_r_2,  0.5_r_2, 0.5_r_2,  0.5_r_2, 0.5_r_2 /)
    xcostnpup = (/   40._r_2, 25._r_2,  40._r_2, 40._r_2,  40._r_2, 40._r_2, 40._r_2, 40._r_2, 40._r_2,  40._r_2, 40._r_2, &
                     40._r_2, 40._r_2,  40._r_2, 40._r_2,  40._r_2, 40._r_2 /)
    xmaxfinelitter = (/ 1524._r_2, 384._r_2, 1527._r_2, 887._r_2, 157._r_2, 361._r_2, 225._r_2, 913._r_2, 660._r_2, 100._r_2, &
                         100._r_2, 100._r_2,  100._r_2,  83._r_2, 100._r_2, 100._r_2, 100._r_2 /)
    xmaxcwd = (/ 1795._r_2, 613._r_2, 1918._r_2, 1164._r_2, 107._r_2, 420._r_2, 228._r_2, 573._r_2, 811._r_2, 100._r_2, &
                  100._r_2, 100._r_2,  100._r_2,   23._r_2, 100._r_2, 100._r_2, 100._r_2 /)
    xnintercept = (/ 6.32_r_2, 4.19_r_2,  6.32_r_2,  5.73_r_2, 14.71_r_2,  6.42_r_2,    2._r_2, 14.71_r_2, 4.71_r_2, 14.71_r_2, &
                    14.71_r_2,   7._r_2, 14.71_r_2, 14.71_r_2, 14.71_r_2, 14.71_r_2, 14.71_r_2 /)
    xnslope = (/ 18.15_r_2, 26.19_r_2, 18.15_r_2, 29.81_r_2, 23.15_r_2, 40.96_r_2, 8._r_2, 23.15_r_2, 59.23_r_2, 23.15_r_2, &
                 23.15_r_2,   10._r_2, 23.15_r_2, 23.15_r_2, 23.15_r_2, 23.15_r_2, 23.15_r_2 /)
    xla_to_sa = (/ 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, &
                   5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2, 5000._r_2 /)
    xvcmax_scalar = (/ 0.92_r_2, 1.10_r_2, 0.92_r_2, 0.92_r_2, 1.25_r_2, 1.25_r_2, 1.25_r_2, 1.25_r_2, 1.0_r_2, 1.0_r_2, &
                        1.0_r_2,  1.0_r_2,  1.0_r_2,  1.0_r_2,  1.0_r_2,  1.0_r_2,  1.0_r_2 /)
    xdisturbance_interval = (/ 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, &
                               100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2, 100._r_2 /)
    xDAMM_EnzPool = (/ 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, &
                       10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2, 10.0_r_2 /)
    xDAMM_KMO2 = (/ 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, &
                    0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2, 0.01_r_2 /)
    xDAMM_KMcp = (/  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2, &
                     0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2,  0.1_r_2 /)
    xDAMM_Ea = (/ 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2, &
                  62._r_2, 62._r_2, 62._r_2, 62._r_2, 62._r_2 /)
    xDAMM_alpha = (/ 10.6_r_2, 10.4_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, &
                     10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2, 10.6_r_2 /)
  
    xxkplab = 0.001369863_r_2
    xxkpsorb = (/ 1.8356191E-05_r_2, 2.0547975E-05_r_2, 1.3698650E-05_r_2, 1.4794542E-05_r_2, 2.1369894E-05_r_2, &
                  2.3835651E-05_r_2, 1.9452083E-05_r_2, 2.1095921E-05_r_2, 2.7123327E-05_r_2, 2.1095921E-05_r_2, &
                  2.7123327E-05_r_2, 2.1095921E-05_r_2 /)
    xxkpocc = 2.73973E-05_r_2

    casabiome%ivt2     =(/        3,        3,        3,        3,        2,        1,   1,   2, &
                                  1,        1,        0,        0,        0,        1,   0,   0, &
                                  0 /)
    casabiome%kroot    =(/      5.5_r_2,      3.9_r_2,      5.5_r_2,      3.9_r_2,      2.0_r_2,      5.5_r_2, 5.5_r_2, 5.5_r_2, &
                                5.5_r_2,      5.5_r_2,      5.5_r_2,      5.5_r_2,      5.5_r_2,      2.0_r_2, 2.0_r_2, 5.5_r_2, &
                                5.5_r_2 /)
    casabiome%rootdepth=(/      1.5_r_2,      1.5_r_2,      1.5_r_2,      1.5_r_2,      0.5_r_2,      0.5_r_2, 0.5_r_2, 0.5_r_2, &
                                0.5_r_2,      0.5_r_2,      0.5_r_2,      0.5_r_2,      0.5_r_2,      0.5_r_2, 0.5_r_2, 1.5_r_2, &
                                0.5_r_2 /)
    casabiome%kuptake  =(/      2.0_r_2,      1.9_r_2,      2.0_r_2,      2.0_r_2,      1.8_r_2,      2.0_r_2, 2.0_r_2, 2.0_r_2, &
                                1.6_r_2,      1.6_r_2,      1.6_r_2,      1.8_r_2,      1.8_r_2,      1.8_r_2, 1.8_r_2, 1.8_r_2, &
                                1.8_r_2 /)
    casabiome%krootlen =(/ 14.87805_r_2, 14.38596_r_2, 14.02597_r_2, 18.94737_r_2, 32.30769_r_2,      84._r_2, 84._r_2, 84._r_2, &
                              120.5_r_2,    120.5_r_2,       0._r_2,       0._r_2,       0._r_2, 30.76923_r_2,  0._r_2,  0._r_2, &
                                 0._r_2 /)
    casabiome%kminN=2.0_r_2
    casabiome%kuplabP=0.5_r_2
    casabiome%fracnpptoP(:,leaf) =(/ 0.25_r_2, 0.20_r_2, 0.40_r_2, 0.35_r_2, 0.35_r_2, 0.35_r_2, 0.35_r_2, 0.50_r_2, 0.50_r_2, &
                                     0.50_r_2, 0.50_r_2, 0.50_r_2, 0.50_r_2, 0.25_r_2, 0.50_r_2, 0.60_r_2, 0.50_r_2 /)
    casabiome%fracnpptoP(:,wood) =(/ 0.40_r_2, 0.35_r_2, 0.30_r_2, 0.25_r_2, 0.25_r_2, 0.00_r_2, 0.00_r_2, 0.10_r_2, 0.00_r_2, &
                                     0.00_r_2, 0.00_r_2, 0.00_r_2, 0.00_r_2, 0.25_r_2, 0.00_r_2, 0.40_r_2, 0.00_r_2 /)
    casabiome%fracnpptoP(:,xroot)=(/ 0.35_r_2, 0.45_r_2, 0.30_r_2, 0.40_r_2, 0.40_r_2, 0.65_r_2, 0.65_r_2, 0.40_r_2, 0.50_r_2, &
                                     0.50_r_2, 0.50_r_2, 0.50_r_2, 0.50_r_2, 0.50_r_2, 0.50_r_2, 0.00_r_2, 0.50_r_2 /)
    casabiome%rmplant(:,leaf)    =0.1_r_2
    casabiome%rmplant(:,wood)    =(/ 2.0_r_2, 1.0_r_2, 1.5_r_2, 0.8_r_2, 0.5_r_2, 0.5_r_2, 0.4_r_2, 1.8_r_2, 2.0_r_2, 1.0_r_2, &
                                     1.0_r_2, 1.0_r_2, 1.0_r_2, 2.0_r_2, 1.0_r_2, 1.0_r_2, 1.0_r_2 /)
    casabiome%rmplant(:,xroot)   =(/ 10._r_2, 2.0_r_2, 7.5_r_2, 2.5_r_2, 4.5_r_2, 4.5_r_2, 4.0_r_2, 15._r_2, 25._r_2, 10._r_2, &
                                     10._r_2, 10._r_2, 10._r_2, 10._r_2, 10._r_2, 10._r_2, 10._r_2 /)
    casabiome%ftransNPtoL(:,leaf) =0.5_r_2
    casabiome%ftransNPtoL(:,wood) =0.95_r_2
    casabiome%ftransNPtoL(:,xroot)=0.9_r_2
    casabiome%fracligninplant(:,leaf) =(/ 0.25_r_2, 0.20_r_2, 0.20_r_2, 0.20_r_2, 0.20_r_2, 0.10_r_2, 0.10_r_2, 0.10_r_2, &
                                          0.10_r_2, 0.10_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.25_r_2, &
                                          0.10_r_2 /)
    casabiome%fracligninplant(:,wood) =0.4_r_2
    casabiome%fracligninplant(:,xroot)=(/ 0.25_r_2, 0.20_r_2, 0.20_r_2, 0.20_r_2, 0.20_r_2, 0.10_r_2, 0.10_r_2, 0.10_r_2, &
                                          0.10_r_2, 0.10_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.15_r_2, 0.25_r_2, &
                                          0.10_r_2 /)
    !casabiome%glaimax=(/ 7._r_2, 7._r_2, 7._r_2, 7._r_2, 3._r_2, 3._r_2, 3._r_2, 3._r_2, 6._r_2, 6._r_2,  5._r_2,  5._r_2, &
    !                     5._r_2, 1._r_2, 6._r_2, 1._r_2, 0._r_2 /)
    casabiome%glaimax=(/ 10._r_2, 10._r_2, 10._r_2, 10._r_2, 10._r_2, 3._r_2, 3._r_2, 3._r_2, 6._r_2, 6._r_2,  5._r_2, 5._r_2, &
                          5._r_2, 1._r_2,  6._r_2,   1._r_2,  0._r_2 /)
    casabiome%glaimin=(/ 1._r_2,  1._r_2, .5_r_2,  .5_r_2, .1_r_2, .1_r_2, .1_r_2, .1_r_2, .1_r_2, .1_r_2, .05_r_2, .05_r_2, &
                        .05_r_2, .05_r_2, 0._r_2, .05_r_2, 0._r_2 /)
    phen%TKshed=(/ 268._r_2,   260._r_2, 263.15_r_2, 268.15_r_2, 277.15_r_2, 275.15_r_2, 275.15_r_2, 275.15_r_2, 278.15_r_2, &
                 278.15_r_2, 277.15_r_2, 277.15_r_2, 277.15_r_2, 277.15_r_2, 277.15_r_2, 277.15_r_2, 283.15_r_2 /)
    casabiome%xkleafcoldexp=3._r_2
    casabiome%xkleafdryexp=3._r_2
    casabiome%ratioNCplantmin(:,leaf) =(/     0.02_r_2,     0.04_r_2, 0.016667_r_2, 0.028571_r_2,    0.025_r_2,  0.02631_r_2, &
                                              0.02_r_2,     0.02_r_2,     0.04_r_2,     0.04_r_2, 0.033333_r_2,    0.025_r_2, &
                                             0.025_r_2, 0.018182_r_2,    0.025_r_2,    0.025_r_2,    0.025_r_2 /)
    casabiome%ratioNCplantmax(:,leaf) =(/    0.024_r_2,    0.048_r_2,     0.02_r_2, 0.034286_r_2,     0.03_r_2, 0.031572_r_2, &
                                             0.024_r_2,    0.024_r_2,    0.048_r_2,    0.048_r_2,     0.04_r_2,     0.03_r_2, &
                                              0.03_r_2, 0.022222_r_2,     0.03_r_2,     0.03_r_2,     0.03_r_2 /)
    casabiome%ratioNCplantmin(:,wood) =(/    0.004_r_2, 0.006667_r_2,    0.004_r_2, 0.005714_r_2, 0.006667_r_2, 0.006667_r_2, &
                                          0.006667_r_2, 0.006667_r_2,    0.008_r_2,    0.008_r_2, 0.006667_r_2, 0.006667_r_2, &
                                          0.006667_r_2, 0.006667_r_2, 0.006667_r_2, 0.007307_r_2, 0.006667_r_2 /)
    casabiome%ratioNCplantmax(:,wood) =(/   0.0048_r_2,    0.008_r_2,   0.0048_r_2, 0.006857_r_2,    0.008_r_2,    0.008_r_2, &
                                             0.008_r_2,    0.008_r_2,   0.0096_r_2,   0.0096_r_2,    0.008_r_2,    0.008_r_2, &
                                             0.008_r_2,    0.008_r_2,    0.008_r_2, 0.008889_r_2,    0.008_r_2 /)
    casabiome%ratioNCplantmin(:,xroot)=(/ 0.012821_r_2, 0.014706_r_2, 0.012821_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2, &
                                          0.014085_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2, &
                                          0.014085_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2, 0.014085_r_2 /)
    casabiome%ratioNCplantmax(:,xroot)=(/ 0.015385_r_2, 0.017647_r_2, 0.015385_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2, &
                                          0.016901_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2, &
                                          0.016901_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2, 0.016901_r_2 /)
    casabiome%ftransPPtoL(:,leaf)=0.5_r_2
    casabiome%ftransPPtoL(:,wood)=0.95_r_2
    casabiome%ftransPPtoL(:,xroot)=0.9_r_2

  end if

end if
    
if ( mp_global>0 ) then
  casabiome%ratioPcplantmin(:,leaf)  = 1._r_2/(xratioNPleafmax*ratioCNplant(:,leaf))
  casabiome%ratioPcplantmax(:,leaf)  = 1._r_2/(xratioNPleafmin*ratioCNplant(:,leaf))
  casabiome%ratioPcplantmin(:,wood)  = 1._r_2/(xratioNPwoodmax*ratioCNplant(:,wood))
  casabiome%ratioPcplantmax(:,wood)  = 1._r_2/(xratioNPwoodmin*ratioCNplant(:,wood))
  casabiome%ratioPcplantmin(:,xroot) = 1._r_2/(xratioNPfrootmax*ratioCNplant(:,xroot))
  casabiome%ratioPcplantmax(:,xroot) = 1._r_2/(xratioNPfrootmin*ratioCNplant(:,xroot))

  casabiome%ratioNPplantmin(:,leaf)  = xratioNPleafmin
  casabiome%ratioNPplantmax(:,leaf)  = xratioNPleafmax
  casabiome%ratioNPplantmin(:,wood)  = xratioNPwoodmin
  casabiome%ratioNPplantmax(:,wood)  = xratioNPwoodmax
  casabiome%ratioNPplantmin(:,xroot) = xratioNPfrootmin
  casabiome%ratioNPplantmax(:,xroot) = xratioNPfrootmax    

  casabiome%sla                = slax
  casabiome%fraclabile(:,leaf) = deltcasa*0.6_r_2    !1/day
  casabiome%fraclabile(:,xroot)= deltcasa*0.4_r_2    !1/day
  casabiome%fraclabile(:,wood) = 0._r_2
  casabiome%plantrate(:,leaf)  = deltcasa/(leafage*(1._r_2-xfherbivore))
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

  casabiome%xnpmax(:)          = xxnpmax(:)
  casabiome%q10soil(:)         = xq10soil(:)
  casabiome%xkoptlitter(:)     = xxkoptlitter(:)
  casabiome%xkoptsoil(:)       = xxkoptsoil(:)
  casabiome%prodptase(:)       = xprodptase(:)/365._r_2   ! convert from yearly to daily
  casabiome%costnpup(:)        = xcostnpup(:)
  casabiome%maxfinelitter(:)   = xmaxfinelitter(:)
  casabiome%maxcwd(:)          = xmaxcwd(:)
  casabiome%nintercept(:)      = xnintercept(:)
  casabiome%nslope(:)          = xnslope(:)    

  !casabiome%la_to_sa(:)             = xla_to_sa(:)
  !casabiome%vcmax_scalar(:)         = xvcmax_scalar(:)
  !casabiome%disturbance_interval(:) = xdisturbance_interval(:)
  !casabiome%DAMM_EnzPool(:)         = xDAMM_EnzPool(:)p
  !casabiome%DAMM_KMO2(:)            = xDAMM_KMO2(:)
  !casabiome%DAMM_KMcp(:)            = xDAMM_KMcp(:)
  !casabiome%DAMM_Ea(:)              = xDAMM_Ea(:)
  !casabiome%DAMM_alpha(:)           = xDAMM_alpha(:)

  casabiome%xkplab = xxkplab
  casabiome%xkpsorb = xxkpsorb
  casabiome%xkpocc = xxkpocc

  casamet%iveg2 = casabiome%ivt2(veg%iveg)
  where (casamet%iveg2==forest.or.casamet%iveg2==shrub)
    casamet%lnonwood = 0
    casapool%cplant(:,wood)  = cwood(veg%iveg) 
    casapool%clitter(:,cwd)  = ccwd(veg%iveg)
    casapool%nplant(:,wood)  = nwood(veg%iveg) 
    casapool%nlitter(:,cwd)  = ncwd(veg%iveg)
    casapool%pplant(:,wood)  = xpwood(veg%iveg)
    casapool%plitter(:,cwd)  = xpcwd(veg%iveg)
  elsewhere
    casamet%lnonwood = 1
    casapool%cplant(:,wood)  = 0._r_2
    casapool%clitter(:,cwd)  = 0._r_2
    casapool%nplant(:,wood)  = 0._r_2
    casapool%nlitter(:,cwd)  = 0._r_2
    casapool%pplant(:,wood)  = 0._r_2
    casapool%plitter(:,cwd)  = 0._r_2
  end where
  if ( cable_pop==1 ) then
   where (casamet%iveg2==forest.or.casamet%iveg2==shrub)
      casapool%cplant(:,wood)  = 0.01_r_2
      casapool%nplant(:,wood)  = casabiome%ratioNCplantmin(veg%iveg,wood)*casapool%cplant(:,wood)
      casapool%pplant(:,wood)  = casabiome%ratioPCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
    end where
  end if
  casapool%cplant(:,leaf)     = cleaf(veg%iveg)
  casapool%cplant(:,xroot)    = cfroot(veg%iveg)
  casapool%clabile            = 0._r_2
  casapool%clitter(:,metb)    = cmet(veg%iveg)
  casapool%clitter(:,str)     = cstr(veg%iveg)
  casapool%csoil(:,mic)       = cmic(veg%iveg)
  casapool%csoil(:,slow)      = cslow(veg%iveg)
  casapool%csoil(:,pass)      = cpass(veg%iveg)
  if ( ccycle==1 ) then
    casapool%ratioNCplant     = 1._r_2/ratioCNplant(veg%iveg,:)  
  end if
  casapool%dclabiledt         = 0._r_2

  ! initializing glai in case not reading pool file (eg. during spin)
  casamet%glai = max(casabiome%glaimin(veg%iveg), casabiome%sla(veg%iveg)*casapool%cplant(:,leaf))
  casaflux%fNminloss   = xfNminloss(veg%iveg)
  casaflux%fNminleach  = 10._r_2*xfNminleach(veg%iveg)*deltcasa
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
  casaflux%fpleach        = xfPleach(casamet%isorder)/365._r_2

  casapool%ratioNCplant   = 1._r_2/ratioCNplant(veg%iveg,:)
  casapool%ratioNPplant   = casabiome%ratioNPplantmin(veg%iveg,:)
  casapool%ratioNClitter  = casapool%nlitter/(casapool%clitter+1.0e-10_r_2)
  casapool%ratioNPlitter  = casapool%nlitter/(casapool%plitter+1.0e-10_r_2)
  casapool%ratioNCsoil    = 1._r_2/ratioCNsoil(veg%iveg,:)
  casapool%ratioNPsoil    = ratioNPsoil(casamet%isorder,:)
  casapool%ratioNCsoilmin = 1._r_2/ratioCNsoilmax(veg%iveg,:)
  casapool%ratioNCsoilmax = 1._r_2/ratioCNsoilmin(veg%iveg,:)
  casapool%ratioNCsoilnew = casapool%ratioNCsoilmax

  casapool%ratioPCplant   = casabiome%ratioPcplantmax(veg%iveg,:)
  casapool%ratioPClitter  = casapool%plitter/(casapool%clitter(:,:)+1.0e-10_r_2)
  casapool%ratioPCsoil    = 1._r_2/(ratioCNsoil(veg%iveg,:)*ratioNPsoil(casamet%isorder,:))

  if ( ccycle<2 ) then
    casapool%Nplant         = casapool%Cplant*casapool%ratioNCplant
    casapool%Nsoil          = casapool%ratioNCsoil*casapool%Csoil
  end if
  if ( ccycle<3 ) then
    casapool%Psoil          = casapool%Nsoil/casapool%ratioNPsoil
    casapool%psoilsorb      = casaflux%psorbmax*casapool%psoillab &
                            /(casaflux%kmlabp+casapool%psoillab)
  end if
end if

end subroutine casa_readbiome

subroutine cable_biophysic_parm(cveg)

use cc_mpi     ! CC MPI routines
use darcdf_m   ! Netcdf data
use infile     ! Input file routines

implicit none

integer k
integer, dimension(mp_global), intent(in) :: cveg
integer, dimension(1) :: nstart, ncount
integer, dimension(:), allocatable, save :: csiropft
real totdepth
real, dimension(:), allocatable, save :: hc, xfang, leaf_w, leaf_l, canst1
real, dimension(:), allocatable, save :: shelrb, extkn, vcmax, rpcoef
real, dimension(:), allocatable, save :: rootbeta, c4frac, vbeta
real, dimension(:), allocatable, save :: a1gs, d0gs, alpha, convex, cfrd
real, dimension(:), allocatable, save :: gswmin, conkc0, conko0, ekc, eko, g0, g1
real, dimension(:), allocatable, save :: zr, clitt
real, dimension(:,:), allocatable, save :: refl, taul, froot2

if ( lncveg_numpft<1 ) then
    
  ! default biophysical parameter tables
  if ( myid==0 ) then
    write(6,*) "Using default CABLE biophysical parameter tables"
  end if
  lncveg_numpft = 18
  allocate( csiropft(lncveg_numpft), hc(lncveg_numpft), xfang(lncveg_numpft), leaf_w(lncveg_numpft), leaf_l(lncveg_numpft) )
  allocate( canst1(lncveg_numpft), shelrb(lncveg_numpft), extkn(lncveg_numpft), refl(lncveg_numpft,2), taul(lncveg_numpft,2) )
  allocate( vcmax(lncveg_numpft), rpcoef(lncveg_numpft), rootbeta(lncveg_numpft), c4frac(lncveg_numpft) )
  allocate( froot2(lncveg_numpft,ms), vbeta(lncveg_numpft) )
  allocate( a1gs(lncveg_numpft), d0gs(lncveg_numpft), alpha(lncveg_numpft), convex(lncveg_numpft), cfrd(lncveg_numpft) )
  allocate( gswmin(lncveg_numpft), conkc0(lncveg_numpft), conko0(lncveg_numpft), ekc(lncveg_numpft), eko(lncveg_numpft) )
  allocate( g0(lncveg_numpft), g1(lncveg_numpft), zr(lncveg_numpft), clitt(lncveg_numpft) )
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
  allocate( csiropft(lncveg_numpft), hc(lncveg_numpft), xfang(lncveg_numpft), leaf_w(lncveg_numpft), leaf_l(lncveg_numpft) )
  allocate( canst1(lncveg_numpft), shelrb(lncveg_numpft), extkn(lncveg_numpft), refl(lncveg_numpft,2), taul(lncveg_numpft,2) )
  allocate( vcmax(lncveg_numpft), rpcoef(lncveg_numpft), rootbeta(lncveg_numpft), c4frac(lncveg_numpft) )
  allocate( froot2(lncveg_numpft,ms), vbeta(lncveg_numpft) )
  allocate( a1gs(lncveg_numpft), d0gs(lncveg_numpft), alpha(lncveg_numpft), convex(lncveg_numpft), cfrd(lncveg_numpft) )
  allocate( gswmin(lncveg_numpft), conkc0(lncveg_numpft), conko0(lncveg_numpft), ekc(lncveg_numpft), eko(lncveg_numpft) )
  allocate( g0(lncveg_numpft), g1(lncveg_numpft), zr(lncveg_numpft), clitt(lncveg_numpft) )

  if ( myid==0 ) then
    nstart(1) = 1
    ncount(1) = lncveg_numpft
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

if ( mp_global>0 ) then

  ! froot is now calculated from soil depth and the new parameter rootbeta 
  ! according to Jackson et al. 1996, Oceologica, 108:389-411
  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + real(soil%zse(k))*100.
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

  if ( maxval(cveg)>lncveg_numpft .or. minval(cveg)<1 ) then
    write(6,*) "ERROR: Invalid range of vegetation classes for CABLE"
    write(6,*) "cveg min,max           = ",minval(cveg),maxval(cveg)
    write(6,*) "Expected range min,max = ",1,lncveg_numpft
    call ccmpi_abort(-1)
  end if

  veg%meth      = 1
  veg%iveg      = csiropft(cveg)
  veg%hc        = real(hc(cveg),8)
  veg%xfang     = real(xfang(cveg),8)  
  veg%dleaf     = real(sqrt(max(leaf_w(cveg)*leaf_l(cveg),1.e-20)),8)
  veg%canst1    = real(canst1(cveg),8)
  veg%shelrb    = real(shelrb(cveg),8)
  veg%extkn     = real(extkn(cveg),8)
  veg%refl(:,1) = real(refl(cveg,1),8)
  veg%refl(:,2) = real(refl(cveg,2),8)  
  veg%taul(:,1) = real(taul(cveg,1),8)
  veg%taul(:,2) = real(taul(cveg,2),8)  
  veg%vcmax     = real(vcmax(cveg),8)
  veg%ejmax     = real(2.*veg%vcmax,8)
  veg%rpcoef    = real(rpcoef(cveg),8)
  do k = 1,ms
    veg%froot(:,k)=real(froot2(cveg,k),8)
  end do
  veg%frac4     = real(c4frac(cveg),8)
  veg%xalbnir   = 1._8 ! not used
  veg%vbeta     = real(vbeta(cveg),8)
  veg%a1gs      = real(a1gs(cveg),8)   
  veg%d0gs      = real(d0gs(cveg),8)
  veg%g0        = real(g0(cveg),8)
  veg%g1        = real(g1(cveg),8)
  veg%alpha     = real(alpha(cveg),8)
  veg%convex    = real(convex(cveg),8) 
  veg%cfrd      = real(cfrd(cveg),8)
  veg%gswmin    = real(gswmin(cveg),8)
  veg%conkc0    = real(conkc0(cveg),8)
  veg%conko0    = real(conko0(cveg),8)
  veg%ekc       = real(ekc(cveg),8)
  veg%eko       = real(eko(cveg),8)
  veg%zr        = real(zr(cveg),8)
  veg%clitt     = real(clitt(cveg),8)
  
  veg%gamma     = 3.e-2_8
  veg%F10       = 0.85_8
  !veg%ZR        = 5._8
  veg%disturbance_interval = 100
  veg%disturbance_intensity = 0._8
  
  ! depeciated
  !veg%tminvj    = real(tminvj(veg%iveg),8)
  !veg%tmaxvj    = real(tmaxvj(veg%iveg),8)
  !veg%rp20      = real(rp20(veg%iveg),8)
  !veg%rs20      = real(rs20(veg%iveg),8)
  !veg%vegcf     = real(vegcf(veg%iveg),8)
  
  ! patch
  if ( gs_switch==1 ) then
    if ( any( veg%g0<1.e-8 ) ) then
      if ( myid==0 ) then
        write(6,*) "WARN: Replacing g0=0. with g0=0.01 for gs_switch=1"
      end if  
      where ( veg%g0<1.e-8 ) 
        veg%g0 = 0.01
      end where
    end if    
  end if    

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
                   
subroutine cable_soil_parm(soil)

use cc_mpi     ! CC MPI routines
use darcdf_m   ! Netcdf data
use infile     ! Input file routines
use newmpar_m
use soilv_m

implicit none

type(soil_parameter_type), intent(inout) :: soil
integer isoil, k
integer, dimension(1) :: nstart, ncount

if ( lncveg_numsoil<1 ) then
    
  ! default soil parameter tables
  if ( myid==0 ) then
    write(6,*) "Using default CABLE soil parameter tables"
  end if

  ! redefine rhos
  rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

  if ( myid==0 ) then
    do isoil = 1,mxst
      write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
    end do
  end if

else

  ! user defined soil parameter tables  
  if ( myid==0 ) then
    write(6,*) "Using user defined CABLE soil parameter tables"
  end if
  if ( lncveg_numsoil > mxst ) then
    write(6,*) "ERROR: Number of soil types larger than maximum, mxst,numsoil:",mxst,lncveg_numsoil
    call ccmpi_abort(-1)
  end if

  silt = 0.
  clay = 0.
  sand = 0.
  swilt = 0.
  sfc = 0.
  ssat = 0.
  bch = 0.
  hyds = 0.
  sucs = 0.
  rhos = 0.
  css = 0.
  swilt(0) = 0.
  sfc(0) = 1.
  ssat(0) = 2.

  if ( myid==0 ) then
    nstart(1) = 1
    ncount(1) = lncveg_numsoil
    call ccnf_get_vara(ncidveg,'silt',nstart,ncount,silt)
    call ccnf_get_vara(ncidveg,'clay',nstart,ncount,clay)
    call ccnf_get_vara(ncidveg,'sand',nstart,ncount,sand)
    call ccnf_get_vara(ncidveg,'swilt',nstart,ncount,swilt(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sfc',nstart,ncount,sfc(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'ssat',nstart,ncount,ssat(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'bch',nstart,ncount,bch)
    call ccnf_get_vara(ncidveg,'hyds',nstart,ncount,hyds)
    call ccnf_get_vara(ncidveg,'sucs',nstart,ncount,sucs)
    call ccnf_get_vara(ncidveg,'rhosoil',nstart,ncount,rhos)
    call ccnf_get_vara(ncidveg,'css',nstart,ncount,css)
  end if
  call ccmpi_bcast(silt,0,comm_world)
  call ccmpi_bcast(clay,0,comm_world)
  call ccmpi_bcast(sand,0,comm_world)
  call ccmpi_bcast(swilt,0,comm_world)
  call ccmpi_bcast(sfc,0,comm_world)
  call ccmpi_bcast(ssat,0,comm_world)
  call ccmpi_bcast(bch,0,comm_world)
  call ccmpi_bcast(hyds,0,comm_world)
  call ccmpi_bcast(sucs,0,comm_world)
  call ccmpi_bcast(rhos,0,comm_world)
  call ccmpi_bcast(css,0,comm_world)
  
  !redo from insoil
  do isoil = 1,mxst
    cnsd(isoil)  = sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
    hsbh(isoil)  = hyds(isoil)*abs(sucs(isoil))*bch(isoil) !difsat*etasat
    ibp2(isoil)  = nint(bch(isoil))+2
    i2bp3(isoil) = 2*nint(bch(isoil))+3
    if ( myid==0 ) then
      write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
    end if
  end do
  cnsd(9) = 2.51

end if

if ( mp_global>0 ) then
  ! Load CABLE soil data
  soil%bch       = real(bch(soil%isoilm),8)
  soil%css       = real(css(soil%isoilm),8)
  soil%rhosoil   = real(rhos(soil%isoilm),8)
  soil%cnsd      = real(cnsd(soil%isoilm),8)
  soil%hyds      = real(hyds(soil%isoilm),8)
  soil%sucs      = real(sucs(soil%isoilm),8)
  soil%hsbh      = real(hsbh(soil%isoilm),8)
  soil%sfc       = real(sfc(soil%isoilm),8)
  soil%ssat      = real(ssat(soil%isoilm),8)
  soil%swilt     = real(swilt(soil%isoilm),8)
  soil%ibp2      = real(ibp2(soil%isoilm),8)
  soil%i2bp3     = real(i2bp3(soil%isoilm),8)
  soil%pwb_min   = (soil%swilt/soil%ssat)**soil%ibp2
  soil%clay      = real(clay(soil%isoilm),8)
  soil%sand      = real(sand(soil%isoilm),8)
  soil%silt      = real(silt(soil%isoilm),8)
  soil%zeta      = 0._8
  soil%fsatmax   = 0._8
  soil%nhorizons = 1
  soil%ishorizon = 1
  do k = 1,ms
    soil%swilt_vec(:,k)   = soil%swilt
    soil%ssat_vec(:,k)    = soil%ssat
    soil%sfc_vec(:,k)     = soil%sfc
    soil%rhosoil_vec(:,k) = soil%rhosoil
    soil%sucs_vec(:,k)    = soil%sucs
    soil%bch_vec(:,k)     = soil%bch
    soil%hyds_vec(:,k)    = soil%hyds
    soil%watr(:,k)        = 0.01_8
  end do
  soil%GWsucs_vec    = soil%sucs*1000._8
  soil%GWhyds_vec    = soil%hyds*1000._8
  soil%GWbch_vec     = soil%bch
  soil%GWrhosoil_vec = soil%rhosoil
  soil%GWssat_vec    = soil%ssat
  soil%GWwatr        = 0.01_8
  soil%GWdz          = 20._8
  soil%GWdz = max( 1._8, min( 20._8, soil%GWdz - sum(soil%zse,dim=1) ) )
  
  ! depeciated
  !bgc%ratecp(:) = real(ratecp(:),8)/
  !bgc%ratecs(:) = real(ratecs(:),8)
end if

return
end subroutine cable_soil_parm
                   
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
    write(6,*) "Please upgrade igbpveg to fix this error"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_attg(ncidveg,'cableformat',cableformat,ierr=iernc)
  if ( iernc/=0 ) then
    cableformat=0.
  end if
  do n = 1,maxtile
    vling(:,n) = 0.
    svsg(:,n) = 0.
    ivsg(:,n) = 0
    write(vname,"(A,I1.1)") "lai",n
    call ccnf_inq_varid(ncidveg,vname,varid,tst)
    if ( .not.tst ) then
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
    end if  
  end do
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
      vling(:,n) = 0.  
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid,tst)
      if ( .not.tst ) then
        call ccnf_inq_varndims(ncidveg,varid,ndims)
        call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n))
      end if  
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
    do n = 1,maxtile
      vling(:,n) = 0.  
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid,tst)
      if ( .not.tst ) then
        call ccnf_inq_varndims(ncidveg,varid,ndims)
        call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n))
      end if  
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
      vling(:,n) = 0.  
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid,tst)
      if ( .not.tst ) then
        call ccnf_inq_varndims(ncidveg,varid,ndims)
        call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n))
      end if  
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
  ivsg = 0
  svsg = 0.
  vling = 0.
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
    ivsg = 0
    svsg = 0.
    vling = 0.
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
    ivsg = 0
    svsg = 0.
    vling = 0.
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
    ivsg = 0
    svsg = 0.
    vling = 0.
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinnext2,vling)
  else
    vlinnext2 = -1.
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
  
integer k
real, dimension(ifull) :: dummy_pack

if ( mp_global>0 ) then
  do k = 1,ms
    call cable_pack(tgg(:,k),ssnow%tgg(:,k))
    call cable_pack(wb(:,k),ssnow%wb(:,k))
    call cable_pack(wbice(:,k),ssnow%wbice(:,k))
  end do
  do k = 1,3
    call cable_pack(tggsn(:,k),ssnow%tggsn(:,k))
    call cable_pack(smass(:,k),ssnow%smass(:,k))
    call cable_pack(ssdn(:,k),ssnow%ssdn(:,k))
    dummy_pack = smass(:,k)/ssdn(:,k)
    call cable_pack(dummy_pack,ssnow%sdepth(:,k))
    ssnow%sconds(:,k) = 0.3_8
  end do      
  call cable_pack(tss,rad%trad(:))
  call cable_pack(ssdnn,ssnow%ssdnn)
  call cable_pack(isflag,ssnow%isflag)
  call cable_pack(snowd,ssnow%snowd)
  call cable_pack(snage,ssnow%snage)
  ssnow%Qrecharge = 0._8
  canopy%sublayer_dz = 0._8
  ssnow%rtevap_sat = 0._8
  ssnow%rtevap_unsat = 0._8
  ssnow%satfrac = 0.5_8
  ssnow%wbliq = ssnow%wb - ssnow%wbice
  ssnow%GWwb = 0.95_8*soil%ssat 
  if ( soil_struc==0 ) then
    do k = 1,ms
      ssnow%wb(:,k) = max(ssnow%wb(:,k), 0.5_8*soil%swilt)
      ssnow%wb(:,k) = min(ssnow%wb(:,k), soil%sfc)
    end do    
  end if
  dummy_pack = real(1-isflag)*tgg(:,1) + real(isflag)*tggsn(:,1) - 273.15
  call cable_pack(dummy_pack,ssnow%tsurface)
  ssnow%rtsoil = 50._8
  canopy%cansto = 0._8
  canopy%ga = 0._8
  canopy%us = 0.01_8
  ssnow%pudsto = 0._8
  ssnow%wetfac = 0._8
  ssnow%osnowd = ssnow%snowd
  canopy%fhs_cor = 0._8
  canopy%fes_cor = 0._8
  canopy%fns_cor = 0._8
  canopy%ga_cor = 0._8
  
  ! default value for fwsoil.  Recaculated by cable_canopy or by SLI
  canopy%fwsoil = max( 1.e-9_8, sum( veg%froot*max(1.e-9_8,min(1._8,ssnow%wb-spread(soil%swilt,2,ms))),2) &
      / ( soil%sfc-soil%swilt ) )
  
  call defaulttile_sli
  call defaulttile_casa
  
  call fixtile
  
end if

return
end subroutine defaulttile

subroutine defaulttile_sli

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
  
integer k

if ( mp_global>0 ) then
    
  if ( soil_struc==1 ) then  
    
    ssnow%h0 = 0._8      ! pond height
    ssnow%snowliq = 0._8 ! liquid snow
    ssnow%Tsoil = ssnow%tgg - 273.15_8
    ssnow%thetai = ssnow%wbice
    do k = 1,ms
      ssnow%S(:,k) = ssnow%wb(:,k)/soil%ssat
    end do
    where ( ssnow%snowd>0. )
      ssnow%nsnow = 1
      ssnow%sdepth(:,1) = ssnow%snowd
      ssnow%smass(:,1) = ssnow%snowd*ssnow%ssdn(:,1)
    elsewhere
      ssnow%nsnow = 0
      ssnow%sdepth(:,1) = 0._8
      ssnow%smass(:,1) = 0._8
    end where
    !ssnow%nsnow = 0
    !ssnow%snowd = 0._8
  
  end if  
    
  call fixtile
  
end if

return
end subroutine defaulttile_sli

subroutine defaulttile_casa

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

if ( mp_global>0 ) then

! MJT notes - may be better to use the default values for
! carbon pools than interpolate low resolution data
    
!  if ( ccycle>=1 .and. ccycle<=3 ) then
!    do k = 1,mplant
!      call cable_pack(cplant(:,k),casapool%cplant(:,k))
!      call cable_pack(niplant(:,k),casapool%nplant(:,k))
!      call cable_pack(pplant(:,k),casapool%pplant(:,k))
!    end do
!    do k = 1,mlitter
!      call cable_pack(clitter(:,k),casapool%clitter(:,k))
!      call cable_pack(nilitter(:,k),casapool%nlitter(:,k))
!      call cable_pack(plitter(:,k),casapool%plitter(:,k))
!    end do
!    do k = 1,msoil
!      call cable_pack(csoil(:,k),casapool%csoil(:,k))
!      call cable_pack(nisoil(:,k),casapool%nsoil(:,k))
!      call cable_pack(psoil(:,k),casapool%psoil(:,k))
!    end do
!  else
!    write(6,*) "ERROR: Unknown ccycle option ",ccycle
!    call ccmpi_abort(-1)
!  end if
 
  call fixtile
  
end if

return
end subroutine defaulttile_casa

subroutine newcbmwb

use soilsnow_m

implicit none

integer k

if ( mp_global>0 ) then
  do k = 1,ms
    call cable_pack(wb(:,k),ssnow%wb(:,k))
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
integer k, n, ierr, idv, ierr_casa, ierr_sli, ierr_pop, ierr_svs
integer jyear,jmonth,jday,jhour,jmin,mins, ll, cc, hh, dd
integer, dimension(5) :: ierr_check
integer, dimension(ifull) :: dati
real, dimension(mp_global) :: dummy_unpack
real, dimension(mp_global) :: old_sv
real, dimension(ifull) :: datr
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(ifull,ms) :: datms
real(kind=8), dimension(ifull,3) :: dat3
real(kind=8), dimension(ifull,mplant) :: datmplant
real(kind=8), dimension(ifull,mlitter) :: datmlitter
real(kind=8), dimension(ifull,msoil) :: datmsoil
real(kind=8), dimension(:,:), allocatable, save :: datpatch
real(kind=8), dimension(:,:), allocatable, save :: datage
real(kind=8), dimension(:,:,:), allocatable, save :: datpc
real(kind=8), dimension(:,:), allocatable, save :: dat91days
real(kind=8), dimension(:,:), allocatable, save :: dat31days
real(kind=8), dimension(:,:), allocatable, save :: dat20years
real(kind=8), dimension(:,:), allocatable, save :: dat5days
logical tst
logical defaultmode
character(len=80) vname
character(len=21) testname

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
ierr_casa = 1
ierr_sli = 1
ierr_pop = 1
ierr_svs = 1
! io_in==1 ensures no interpolation is required
if ( io_in==1 .and. .not.defaultmode ) then
  if ( myid==0 .or. pfall ) then
    write(testname,'("t",I1.1,"_tgg1")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr = 0
    end if
    write(testname,'("t",I1.1,"_cplant1")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_casa = 0
    end if
    write(testname,'("t",I1.1,"_hzero")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_sli = 0
    end if
    write(testname,'("t",I1.1,"_pop_grid_cmass_sum")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_pop = 0
    end if
    write(testname,'("t",I1.1,"_svs")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_svs = 0
    end if
  end if
  if ( .not.pfall ) then
    ierr_check(1) = ierr
    ierr_check(2) = ierr_casa
    ierr_check(3) = ierr_sli
    ierr_check(4) = ierr_pop
    ierr_check(5) = ierr_svs
    call ccmpi_bcast(ierr_check(1:5),0,comm_world)
    ierr      = ierr_check(1)
    ierr_casa = ierr_check(2)
    ierr_sli  = ierr_check(3)
    ierr_pop  = ierr_check(4)
    ierr_svs  = ierr_check(5)
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
    if ( ierr_svs /= 0 ) then
      old_sv = sv    
    else    
      write(vname,'("t",I1.1,"_svs")') n
      call histrd(iarchi-1,ierr,vname,datr,ifull)
      if ( n<=maxnb ) then
        call cable_pack(datr,old_sv,n)
      end if
    end if  
    write(vname,'("t",I1.1,"_tgg")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%tgg(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_wb")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%wb(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_wbice")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    if ( n<=maxnb ) then
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%wbice(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_tggsn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        call cable_pack(dat3(:,k),ssnow%tggsn(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_smass")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        call cable_pack(dat3(:,k),ssnow%smass(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_ssdn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        call cable_pack(dat3(:,k),ssnow%ssdn(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        call cable_pack(dat3(:,k),ssnow%sdepth(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_sconds")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    if ( n<=maxnb ) then
      do k = 1,3
        call cable_pack(dat3(:,k),ssnow%sconds(:,k),n)
      end do
    end if
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%ssdnn(:),n)
    write(vname,'("t",I1.1,"_sflag")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)
    if ( n<=maxnb ) call cable_pack(dati,ssnow%isflag(:),n)
    write(vname,'("t",I1.1,"_snd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%snowd(:),n)
    write(vname,'("t",I1.1,"_osnd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%osnowd(:),n)
    write(vname,'("t",I1.1,"_snage")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%snage(:),n)
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%rtsoil(:),n)
    write(vname,'("t",I1.1,"_GWwb")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%GWwb(:),n)
    write(vname,'("t",I1.1,"_cansto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,canopy%cansto(:),n)
    write(vname,'("t",I1.1,"_us")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,canopy%us(:),n)
    write(vname,'("t",I1.1,"_pudsto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%pudsto(:),n)
    write(vname,'("t",I1.1,"_wetfac")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,ssnow%wetfac(:),n)
    write(vname,'("t",I1.1,"_ga")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    if ( n<=maxnb ) call cable_pack(dat,canopy%ga(:),n)
  end do
  if ( soil_struc==1 ) then
    if ( ierr_sli/=0 ) then
      if ( myid==0 ) write(6,*) "Use gridbox averaged data to initialise SLI"
      call defaulttile_sli
    else    
      do n = 1,maxtile
        write(vname,'("t",I1.1,"_hzero")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,ssnow%h0(:),n)
        write(vname,'("t",I1.1,"_s")') n
        call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
        if ( n<=maxnb ) then
          do k = 1,ms
            call cable_pack(datms(:,k),ssnow%S(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_tsoil")') n
        call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
        if ( n<=maxnb ) then
          do k = 1,ms
            call cable_pack(datms(:,k),ssnow%tsoil(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_thetai")') n
        call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
        if ( n<=maxnb ) then
          do k = 1,ms
            call cable_pack(datms(:,k),ssnow%thetai(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,ssnow%snowliq(:,1),n) ! currently nsnow_max=1
        write(vname,'("t",I1.1,"_tsurface")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,ssnow%tsurface(:),n)
        write(vname,'("t",I1.1,"_nsnow")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)
        if ( n<=maxnb ) call cable_pack(dati,ssnow%nsnow(:),n)
        write(vname,'("t",I1.1,"_fwsoil")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,canopy%fwsoil(:),n)
      end do  
    end if
  end if
  if ( ccycle/=0 ) then
    if ( ierr_casa/=0 ) then
      if ( myid==0 ) write(6,*) "Use gridbox averaged data to initialise CASA-CNP"
      call defaulttile_casa
    else
      if ( myid==0 ) write(6,*) "Use tiled data to initialise CASA-CNP"  
      do n = 1,maxtile
        write(vname,'("t",I1.1,"_cplant")') n
        call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
        if ( n<=maxnb ) then
          do k = 1,mplant
            call cable_pack(datmplant(:,k),casapool%cplant(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_nplant")') n
        call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
        if ( n<=maxnb ) then
          do k = 1,mplant
            call cable_pack(datmplant(:,k),casapool%nplant(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_pplant")') n
        call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
        if ( n<=maxnb ) then
          do k = 1,mplant
            call cable_pack(datmplant(:,k),casapool%pplant(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_clitter")') n
        call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
        if ( n<=maxnb ) then
          do k = 1,mlitter
            call cable_pack(datmlitter(:,k),casapool%clitter(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_nlitter")') n
        call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
        if ( n<=maxnb ) then
          do k = 1,mlitter
            call cable_pack(datmlitter(:,k),casapool%nlitter(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_plitter")') n
        call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
        if ( n<=maxnb ) then
          do k = 1,mlitter
            call cable_pack(datmlitter(:,k),casapool%plitter(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_csoil")') n
        call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
        if ( n<=maxnb ) then
          do k = 1,msoil
            call cable_pack(datmsoil(:,k),casapool%csoil(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_nsoil")') n
        call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
        if ( n<=maxnb ) then
          do k = 1,msoil
            call cable_pack(datmsoil(:,k),casapool%nsoil(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_psoil")') n
        call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
        if ( n<=maxnb ) then
          do k = 1,msoil
            call cable_pack(datmsoil(:,k),casapool%psoil(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_glai")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casamet%glai,n)
        write(vname,'("t",I1.1,"_phen")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,phen%phen,n)
        write(vname,'("t",I1.1,"_aphen")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,phen%aphen,n)
        write(vname,'("t",I1.1,"_phenphase")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)
        if ( n<=maxnb ) call cable_pack(dati,phen%phase,n)
        write(vname,'("t",I1.1,"_doyphase3")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)
        if ( n<=maxnb ) call cable_pack(dati,phen%doyphase(:,3),n)
        write(vname,'("t",I1.1,"_clabile")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casapool%clabile,n)
        write(vname,'("t",I1.1,"_nsoilmin")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casapool%nsoilmin,n)
        write(vname,'("t",I1.1,"_psoillab")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casapool%psoillab,n)
        write(vname,'("t",I1.1,"_psoilsorb")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casapool%psoilsorb,n)
        write(vname,'("t",I1.1,"_psoilocc")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casapool%psoilocc,n)
        write(vname,'("t",I1.1,"_crmplant")') n
        call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
        if ( n<=maxnb ) then
          do k = 1,mplant
            call cable_pack(datmplant(:,k),casaflux%crmplant(:,k),n)
          end do
        end if
        write(vname,'("t",I1.1,"_fracsapwood")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%frac_sapwood,n)
        write(vname,'("t",I1.1,"_sapwoodarea")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%sapwood_area,n)
        write(vname,'("t",I1.1,"_crsoil")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%crsoil,n)
        write(vname,'("t",I1.1,"_cnpp")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%cnpp,n)
        write(vname,'("t",I1.1,"_clabloss")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%clabloss,n)
        write(vname,'("t",I1.1,"_crgplant")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%crgplant,n)
        write(vname,'("t",I1.1,"_stemnpp")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casaflux%stemnpp,n)
        write(vname,'("t",I1.1,"_LAImax")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casabal%laimax,n)
        write(vname,'("t",I1.1,"_Cleafmean")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casabal%cleafmean,n)
        write(vname,'("t",I1.1,"_Crootmean")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,casabal%crootmean,n)
        write(vname,'("t",I1.1,"_fpn")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,canopy%fpn,n)
        write(vname,'("t",I1.1,"_frday")') n
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call cable_pack(dat,canopy%frday,n)
      end do  
    end if
  end if
  if ( cable_climate==1 ) then
    allocate( dat91days(ifull,91) )  
    allocate( dat31days(ifull,31) )
    allocate( dat20years(ifull,20) )
    allocate( dat5days(ifull,120) )
    dat91days = 0._8
    dat31days = 0._8
    dat20years = 0._8
    dat5days = 0._8
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_climateiveg")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      if ( n<=maxnb ) call cable_pack(dati,climate%iveg,n)
      write(vname,'("t",I1.1,"_climatebiome")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      if ( n<=maxnb ) call cable_pack(dati,climate%biome,n)
      write(vname,'("t",I1.1,"_climateevapPT")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%evap_PT,n)
      write(vname,'("t",I1.1,"_climateaevap")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%aevap,n)
      write(vname,'("t",I1.1,"_climatemtempmax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%mtemp_max,n)
      write(vname,'("t",I1.1,"_climatemtempmin")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%mtemp_min,n)
      write(vname,'("t",I1.1,"_climatedmoistmax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%dmoist_max,n)
      write(vname,'("t",I1.1,"_climatedmoistmin")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%dmoist_min,n)
      write(vname,'("t",I1.1,"_climateqtempmax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%qtemp_max,n)
      write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%qtemp_max_last_year,n)
      write(vname,'("t",I1.1,"_climategdd0")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%gdd0,n)
      write(vname,'("t",I1.1,"_climategdd5")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%gdd5,n)
      write(vname,'("t",I1.1,"_climateagdd0")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%agdd0,n)
      write(vname,'("t",I1.1,"_climateagdd5")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%agdd5,n)
      write(vname,'("t",I1.1,"_climatechilldays")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)      
      if ( n<=maxnb ) call cable_pack(dati,climate%chilldays,n)
      write(vname,'("t",I1.1,"_climategdd0_rec")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%gdd0_rec,n)
      write(vname,'("t",I1.1,"_climatemtempmin20")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%mtemp_min20,n)
      write(vname,'("t",I1.1,"_climatemtempmax20")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%mtemp_max20,n)
      write(vname,'("t",I1.1,"_climatealphaPT20")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%alpha_PT20,n)
      write(vname,'("t",I1.1,"_climatedmoistmin20")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%dmoist_min20,n)
      write(vname,'("t",I1.1,"_climatedmoistmax20")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%dmoist_max20,n)
      write(vname,'("t",I1.1,"_climatefrec")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%frec,n)
      write(vname,'("t",I1.1,"_climatefdorm")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,climate%fdorm,n)
      write(vname,'("t",I1.1,"_climategmd")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)      
      if ( n<=maxnb ) call cable_pack(dati,climate%gmd,n)
      write(vname,'("t",I1.1,"_climatedtemp")') n
      call histrd(iarchi-1,ierr,vname,dat91days,ifull)
      if ( n<=maxnb ) then
        do k = 1,91  
          call cable_pack(dat91days(:,k),climate%dtemp_91(:,k),n)
          if ( k>60 ) then
            call cable_pack(dat91days(:,k),climate%dtemp_31(:,k-60),n)
          end if
        end do
      end if
      write(vname,'("t",I1.1,"_climatedmoist")') n
      call histrd(iarchi-1,ierr,vname,dat31days,ifull)
      if ( n<=maxnb ) then
        do k = 1,31  
          call cable_pack(dat31days(:,k),climate%dmoist_31(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatemtemp_min_20")') n
      call histrd(iarchi-1,ierr,vname,dat20years,ifull)
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_pack(dat20years(:,k),climate%mtemp_min_20(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatemtemp_max_20")') n
      call histrd(iarchi-1,ierr,vname,dat20years,ifull)
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_pack(dat20years(:,k),climate%mtemp_max_20(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatealpha_PT_20")') n
      call histrd(iarchi-1,ierr,vname,dat20years,ifull)
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_pack(dat20years(:,k),climate%alpha_PT_20(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatedmoist_min_20")') n
      call histrd(iarchi-1,ierr,vname,dat20years,ifull)
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_pack(dat20years(:,k),climate%dmoist_min_20(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatedmoist_max_20")') n
      call histrd(iarchi-1,ierr,vname,dat20years,ifull)
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_pack(dat20years(:,k),climate%dmoist_max_20(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climateapar_leaf_sun")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%apar_leaf_sun(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climateapar_leaf_shade")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%apar_leaf_shade(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatedleaf_sun")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%dleaf_sun(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatedleaf_shade")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%dleaf_shade(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatetleaf_sun")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%tleaf_sun(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatetleaf_shade")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%tleaf_shade(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatecs_sun")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%cs_sun(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatecs_shade")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%cs_shade(:,k),n)
        end do  
      end if
          write(vname,'("t",I1.1,"_climatescalex_sun")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%scalex_sun(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_climatescalex_shade")') n
      call histrd(iarchi-1,ierr,vname,dat5days,ifull)
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_pack(dat5days(:,k),climate%scalex_shade(:,k),n)
        end do  
      end if
      write(vname,'("t",I1.1,"_vegvcmax_sun")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,veg%vcmax_sun,n)
      write(vname,'("t",I1.1,"_vegvcmax_shade")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,veg%vcmax_shade,n)
      write(vname,'("t",I1.1,"_vegejmax_sun")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,veg%ejmax_sun,n)
      write(vname,'("t",I1.1,"_vegejmax_shade")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      if ( n<=maxnb ) call cable_pack(dat,veg%ejmax_shade,n)
    end do
    deallocate( dat91days )
    deallocate( dat31days )
    deallocate( dat20years )
    deallocate( dat5days )
  end if
  if ( cable_pop==1 ) then
    if ( ierr_pop/=0 ) then
      if ( myid==0 ) write(6,*) "Use default data to initialise POP"
    else
      if ( myid==0 ) write(6,*) "Use tiled data to initialise POP"    
      allocate( datpatch(ifull,POP_NPATCH) )  
      allocate( datage(ifull,POP_AGEMAX) )  
      allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
      datpatch = 0._8
      datage = 0._8
      datpc = 0._8
      do n = 1,maxtile  
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%cmass_sum,n)
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call  pop_pack(dat,pop%pop_grid(:)%cmass_sum_old,n)
        write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%cheartwood_sum,n)
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%csapwood_sum,n)
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%csapwood_sum_old,n)
        write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%densindiv,n)
        write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%height_mean,n)
        write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%height_max,n)        
        write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%basal_area,n)
        write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%sapwood_loss,n)
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%sapwood_area_loss,n)
        write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%stress_mortality,n)
        write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%crowding_mortality,n)
        write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%fire_mortality,n)
        write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%cat_mortality,n)
        write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%res_mortality,n)
        write(vname,'("t",I1.1,"_pop_grid_growth")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%growth,n)        
        write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%area_growth,n)
        write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%crown_cover,n)
        write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%crown_area,n)
        write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%crown_volume,n)
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%sapwood_area,n)        
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%sapwood_area_old,n)
        write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%KClump,n)    
        write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
        call histrd(iarchi-1,ierr,vname,datage,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_AGEMAX
            call pop_pack(datage(:,k),pop%pop_grid(:)%freq_age(k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
        call histrd(iarchi-1,ierr,vname,datage,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_AGEMAX
            call pop_pack(datage(:,k),pop%pop_grid(:)%biomass_age(k),n)
          end do
        end if  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%biomass(ll),n)
        end do
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%density(ll),n)
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%hmean(ll),n)
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%hmax(ll),n)
        end do
        do hh = 1,POP_HEIGHT_BINS
          write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%cmass_stem_bin(hh),n) 
        end do  
        do hh = 1,POP_HEIGHT_BINS
          write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%densindiv_bin(hh),n)
        end do
        do hh = 1,POP_HEIGHT_BINS
          write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh  
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%height_bin(hh),n)        
        end do
        do hh = 1,POP_HEIGHT_BINS
          write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          if ( n<=maxnb ) call pop_pack(dat,pop%pop_grid(:)%diameter_bin(hh),n)
        end do
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd  
          call histrd(iarchi-1,ierr,vname,dat,ifull)
          dati = nint(dat)
          if ( n<=maxnb ) call pop_pack(dati,pop%pop_grid(:)%n_age(dd),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb )  then
          do k = 1,POP_NPATCH
            dati = nint(datpatch(:,k))  
            call pop_pack(dati,pop%pop_grid(:)%patch(k)%id,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_freq")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%freq(k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%freq_old(k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%factor_recruit,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%pgap,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%lai,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%biomass,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%biomass_old,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%heartwood,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood_old,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood_area,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood_area_old,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%stress_mortality,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%fire_mortality,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%cat_mortality,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%crowding_mortality,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%cpc,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%mortality,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood_loss,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%sapwood_area_loss,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%growth,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%area_growth,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%frac_NPP,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%frac_respiration,n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        if ( n<=maxnb ) then
          do k = 1,POP_NPATCH
            call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%frac_light_uptake,n)
          end do
        end if  
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH
              dati = nint(datpatch(:,k))  
              call pop_pack(dati,pop%pop_grid(:)%patch(k)%disturbance_interval(dd),n)
            end do
          end if  
        end do  
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH
              dati = nint(datpatch(:,k))   
              call pop_pack(dati,pop%pop_grid(:)%patch(k)%first_disturbance_year(dd),n)
            end do
          end if  
        end do  
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH
              dati = nint(datpatch(:,k)) 
              call pop_pack(dati,pop%pop_grid(:)%patch(k)%age(dd),n)
            end do
          end if  
        end do  
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH
              dati = nint(datpatch(:,k))  
              call pop_pack(dati,pop%pop_grid(:)%ranked_age_unique(k,dd),n)
            end do
          end if  
        end do  
        do dd = 1,POP_NDISTURB
          write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH
              call pop_pack(datpatch(:,k),pop%pop_grid(:)%freq_ranked_age_unique(k,dd),n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH  
              dati = nint(datpatch(:,k))  
              call pop_pack(dati,pop%pop_grid(:)%patch(k)%layer(ll)%ncohort,n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH  
              call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%layer(ll)%biomass,n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH  
              call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%layer(ll)%density,n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH  
              call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%layer(ll)%hmean,n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          call histrd(iarchi-1,ierr,vname,datpatch,ifull)
          if ( n<=maxnb ) then
            do k = 1,POP_NPATCH  
              call pop_pack(datpatch(:,k),pop%pop_grid(:)%patch(k)%layer(ll)%hmax,n)
            end do
          end if  
        end do  
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT                  
              do k = 1,POP_NPATCH
                dati = nint(datpc(:,k,cc))  
                call pop_pack(dati,pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT    
              do k = 1,POP_NPATCH
                dati = nint(datpc(:,k,cc))  
                call pop_pack(dati,pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT    
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT    
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT    
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT    
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake,n)
              end do  
            end do
          end if  
        end do        
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar,n)
              end do  
            end do
          end if  
        end do         
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height,n)
              end do  
            end do
          end if  
        end do         
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area,n)
              end do  
            end do
          end if  
        end do         
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf,n)
              end do  
            end do
          end if  
        end do 
        do ll = 1,POP_NLAYER
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          call histrd(iarchi-1,ierr,vname,datpc,ifull)
          if ( n<=maxnb ) then
            do cc = 1,POP_NCOHORT
              do k = 1,POP_NPATCH  
                call pop_pack(datpc(:,k,cc),pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot,n)
              end do  
            end do
          end if  
        end do 
      end do
      deallocate( datpatch )
      deallocate( datage )
      deallocate( datpc )
    end if      
  end if    
  ! CABLE correction terms
  !do n = 1,maxtile
  !  write(vname,'("t",I1.1,"_fhscor")') n
  !  call histrd(iarchi-1,ierr,vname,dat,ifull)
  !  if ( n<=maxnb ) call cable_pack(dat,canopy%fhs_cor,n)
  !  write(vname,'("t",I1.1,"_fescor")') n
  !  call histrd(iarchi-1,ierr,vname,dat,ifull)
  !  if ( n<=maxnb ) call cable_pack(dat,canopy%fes_cor,n)
  !  write(vname,'("t",I1.1,"_fnscor")') n
  !  call histrd(iarchi-1,ierr,vname,dat,ifull)
  !  if ( n<=maxnb ) call cable_pack(dat,canopy%fns_cor,n)
  !  write(vname,'("t",I1.1,"_gacor")') n
  !  call histrd(iarchi-1,ierr,vname,dat,ifull)
  !  if ( n<=maxnb ) call cable_pack(dat,canopy%ga_cor,n)
  !end do
  ! albvisdir, albvisdif, albnirdir, albnirdif are used when nrad=5
  vname = 'albvisdir'
  call histrd(iarchi-1,ierr,vname,albvisdir,ifull)
  vname = 'albvisdif'
  call histrd(iarchi-1,ierr,vname,albvisdif,ifull)
  vname = 'albnirdir'
  call histrd(iarchi-1,ierr,vname,albnirdir,ifull)
  vname = 'albnirdif'
  call histrd(iarchi-1,ierr,vname,albnirdif,ifull)
  ! albvissav and albnirsav are used when nrad=4
  vname = 'albvis'
  call histrd(iarchi-1,ierr,vname,albvissav,ifull)
  vname = 'albnir'
  call histrd(iarchi-1,ierr,vname,albnirsav,ifull)
  
  call redistribute_tile(old_sv)
  call fixtile
  
end if

! Calculate LAI and veg fraction diagnostics
vlai(:) = 0.
sigmf(:) = 0.
if ( mp_global>0 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jmonth,jday,jhour,jmin,mp_global,sv,vl1,vl2,vl3,vl4,casamet,veg,ifull)
  dummy_unpack = sv*real(veg%vlai)
  call cable_unpack(dummy_unpack,vlai)
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
if ( mp_global>0 ) then

  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + real(soil%zse(k))*100.
  enddo
  
  ssnow%tgg = max(ssnow%tgg, 200._8)
  ssnow%tggsn = max(ssnow%tggsn, 200._8)

  ssnow%wb = max(ssnow%wb,0._8)
  ssnow%wbice = min( max(ssnow%wbice, 0._8), ssnow%wb )
  ssnow%smass = max(ssnow%smass,0._8)
  ssnow%rtsoil = max(ssnow%rtsoil,0._8)
  ssnow%snowd = max(ssnow%snowd,0._8)
  ssnow%osnowd = max(ssnow%osnowd,0._8)
  ssnow%wetfac = min(max(ssnow%wetfac,0._8),1._8)
  canopy%cansto = max(canopy%cansto,0._8)
  ssnow%wbliq = ssnow%wb - ssnow%wbice

  ssnow%wbtot = 0._8
  ssnow%wbtot1 = 0._8
  ssnow%wbtot2 = 0._8
  ssnow%tggav = 0._8
  do k = 1,ms
    ssnow%wbtot = ssnow%wbtot+ssnow%wb(:,k)*1000._8*soil%zse(k)
    ssnow%tggav = ssnow%tggav+soil%zse(k)*ssnow%tgg(:,k)/real(totdepth/100.,8)
    ssnow%gammzz(:,k) = max((1._8-soil%ssat)*soil%css* soil%rhosoil                         &
        + real(ssnow%wb(:,k)-ssnow%wbice(:,k))*4.218e3_8* 1000._8                           &
        + real(ssnow%wbice(:,k))*2.100e3_8*1000._8*0.9_8,soil%css*soil%rhosoil)*soil%zse(k) &
        + (1.-ssnow%isflag)*2090._8*ssnow%snowd
  end do

  if ( ccycle == 0 ) then
    !bgc%cplant = max(bgc%cplant,0.)
    !bgc%csoil = max(bgc%csoil,0.)
  else
    casapool%cplant     = max(0._8,casapool%cplant)
    casapool%clitter    = max(0._8,casapool%clitter)
    casapool%csoil      = max(0._8,casapool%csoil)
    casabal%cplantlast  = casapool%cplant
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal      = 0._8
    casabal%FCgppyear    = 0._8
    casabal%FCrpyear     = 0
    casabal%FCrmleafyear = 0._8
    casabal%FCrmwoodyear = 0._8
    casabal%FCrmrootyear = 0._8
    casabal%FCrgrowyear  = 0._8
    casabal%FCnppyear    = 0._8
    casabal%FCrsyear     = 0._8
    casabal%FCneeyear    = 0._8
    casabal%dCdtyear     = 0._8
    casapool%nplant      = max(1.e-6_8,casapool%nplant)
    casapool%nlitter     = max(1.e-6_8,casapool%nlitter)
    casapool%nsoil       = max(1.e-6_8,casapool%nsoil)
    casapool%nsoilmin    = max(1.e-6_8,casapool%nsoilmin)
    casabal%nplantlast   = casapool%nplant
    casabal%nlitterlast  = casapool%nlitter
    casabal%nsoillast    = casapool%nsoil
    casabal%nsoilminlast = casapool%nsoilmin
    casabal%sumnbal      = 0._8
    casabal%FNdepyear    = 0._8
    casabal%FNfixyear    = 0._8
    casabal%FNsnetyear   = 0._8
    casabal%FNupyear     = 0._8
    casabal%FNleachyear  = 0._8
    casabal%FNlossyear   = 0._8
    casapool%pplant      = max(1.0e-7_8,casapool%pplant)
    casapool%plitter     = max(1.0e-7_8,casapool%plitter)
    casapool%psoil       = max(1.0e-7_8,casapool%psoil)
    casapool%Psoillab    = max(1.0e-7_8,casapool%psoillab)
    casapool%psoilsorb   = max(1.0e-7_8,casapool%psoilsorb)
    casapool%psoilocc    = max(1.0e-7_8,casapool%psoilocc)
    casabal%pplantlast   = casapool%pplant
    casabal%plitterlast  = casapool%plitter
    casabal%psoillast    = casapool%psoil
    casabal%psoillablast = casapool%psoillab
    casabal%psoilsorblast= casapool%psoilsorb
    casabal%psoilocclast = casapool%psoilocc
    casabal%sumpbal      = 0._8
    casabal%FPweayear    = 0._8
    casabal%FPdustyear   = 0._8
    casabal%FPsnetyear   = 0._8
    casabal%FPupyear     = 0._8
    casabal%FPleachyear  = 0._8
    casabal%FPlossyear   = 0._8
    !casamet%glai         = max(min( casamet%glai, casabiome%glaimax(veg%iveg)), casabiome%glaimin(veg%iveg))
  end if
end if
  
return
end subroutine fixtile

subroutine redistribute_tile(old_sv)

implicit none

integer k
real, dimension(mp_global), intent(in) :: old_sv

if ( mp_global>0 ) then

  if ( any( abs(sv-old_sv)>1.e-8 ) ) then
    
    ! assume common soil texture and soil heat capacity
    do k = 1,ms
      call redistribute_work(old_sv,ssnow%tgg(:,k))
      call redistribute_work(old_sv,ssnow%wb(:,k))
      call redistribute_work(old_sv,ssnow%wbice(:,k))
    end do
    call redistribute_work(old_sv,ssnow%GWwb)
    if ( soil_struc==1 ) then
      do k = 1,ms
        call redistribute_work(old_sv,ssnow%tsoil(:,k))
      end do
    end if
    !do k = 1,3
    !  call redistribute_work(old_sv,ssnow%tggsn(:,k))
    !  call redistribute_work(old_sv,ssnow%smass(:,k))
    !  call redistribute_work(old_sv,ssnow%ssdn(:,k))
    !  call redistribute_work(old_sv,ssnow%sdepth(:,k))
    !end do
    !call redistribute_work(old_sv,ssnow%ssdn(:))
    !call redistribute_work(old_sv,ssnow%snowd(:))
    !call redistribute_work(old_sv,ssnow%osnowd(:))
  
  end if
  
end if
  
return
end subroutine redistribute_tile

subroutine redistribute_work(old_sv,vdata)

use cc_omp, only : ntiles, imax
use newmpar_m
use soil_m

implicit none

integer tile, nb, iq, is, ie
real, dimension(mp_global), intent(in) :: old_sv
real, dimension(ifull,maxnb) :: up_new_svs, up_old_svs
real, dimension(ifull) :: svs_sum
real, dimension(maxtile) :: adj_pos_frac, adj_neg_frac, new_svs, old_svs
real(kind=8) :: ave_neg_vdata
real(kind=8), dimension(mp_global), intent(inout) :: vdata
real(kind=8), dimension(ifull,maxnb) :: up_vdata
real(kind=8), dimension(maxtile) :: old_vdata

up_vdata = 0._8
up_new_svs = 0.
up_old_svs = 0.

do nb = 1,maxnb
  call cable_unpack(sv,up_new_svs(:,nb),nb)
  call cable_unpack(old_sv,up_old_svs(:,nb),nb)
  call cable_unpack(vdata,up_vdata(:,nb),nb)
end do  
svs_sum = sum(up_new_svs,dim=2)
do nb = 1,maxnb
  up_new_svs(:,nb) = up_new_svs(:,nb)/svs_sum
end do
svs_sum = sum(up_old_svs,dim=2)
do nb = 1,maxnb
  up_old_svs(:,nb) = up_old_svs(:,nb)/svs_sum
end do

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  do iq = is,ie
    new_svs(:) = up_new_svs(iq,:)
    old_svs(:) = up_old_svs(iq,:)      
    if ( land(iq) .and. any( abs(new_svs-old_svs)>1.e-8 ) ) then 
      adj_pos_frac(:) = max( new_svs(:)-old_svs(:), 0. ) ! tiles with increasing fraction
      adj_neg_frac(:) = max( old_svs(:)-new_svs(:), 0. ) ! tiles with decreasing fraction
      old_vdata(:) = up_vdata(iq,:)
      ! combine tiles with decreasing area fraction into one tile
      ave_neg_vdata = sum( adj_neg_frac(:)*old_vdata )/sum( adj_neg_frac )
      ! Only change tiles that are increasing in area fraction
      do nb = 1,tdata(tile)%maxnb
        if ( adj_pos_frac(nb)>0. ) then  
          up_vdata(iq,nb) = (old_vdata(nb)*old_svs(nb) + ave_neg_vdata*adj_pos_frac(nb)) &
                        /(old_svs(nb) + adj_pos_frac(nb))
        end if  
      end do
    end if    
  end do  
end do

do nb = 1,maxnb
  call cable_pack(up_vdata(:,nb),vdata,nb)
end do  

return
end subroutine redistribute_work    
    
! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetiledef(idnc,local,jdim,jsize,c1dim,c2dim,c3dim,c4dim,c5dim,c6dim,c7dim,c1size,c2size,itype)

use carbpools_m
use cc_mpi, only : myid
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
  
implicit none
  
integer, intent(in) :: idnc, jsize, c1size, c2size
integer k,n
integer ll,dd,hh
integer, dimension(jsize), intent(in) :: jdim  
integer, dimension(c1size), intent(in) :: c1dim, c3dim, c4dim, c5dim, c6dim, c7dim
integer, dimension(c2size), intent(in) :: c2dim
character(len=80) vname
character(len=80) lname
logical, intent(in) :: local
integer, intent(in) :: itype
  
if (myid==0.or.local) then
  if (myid==0) then
    write(6,*) "Defining CABLE tile data"
  end if
  if ( itype==-1 ) then !just for restart file
    do n = 1,maxtile
      write(lname,'("Veg fraction tile ",I1.1)') n
      write(vname,'("t",I1.1,"_svs")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.,0,2) ! kind=8        
      do k=1,ms
        write(lname,'("Soil temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,0,2) ! kind=8
      end do
      do k=1,ms
        write(lname,'("Soil moisture tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wb",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,0,2) ! kind=8
      end do
      do k=1,ms
        write(lname,'("Soil ice tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wbice",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,0,2) ! kind=8
      end do
      do k=1,3
        write(lname,'("Snow temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,0,2) ! kind=8
      end do
      do k=1,3
        write(lname,'("Snow mass tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_smass",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'K',0.,650.,0,2) ! kind=8
      end do
      do k=1,3
        write(lname,'("Snow density tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,0,2) ! kind=8
      end do
      do k=1,3
        write(lname,'("Snow depth tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,0,2) ! kind=8
      end do
      do k=1,3
        write(lname,'("Snow sconds tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sconds",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,0,2) ! kind=8
      end do
      write(lname,'("Snow ssdnn tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ssdnn")') n
      call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,0,2) ! kind=8
      write(lname,'("Snow flag tile ",I1.1)') n
      write(vname,'("t",I1.1,"_sflag")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6.5,0,2) ! kind=8
      write(lname,'("Snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,0,2) ! kind=8
      write(lname,'("Old snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_osnd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,0,2) ! kind=8
      write(lname,'("Snow age tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snage")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,26.,0,2) ! kind=8
      write(lname,'("Soil turbulent resistance tile ",I1.1)') n
      write(vname,'("t",I1.1,"_rtsoil")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3e5,0,2) ! kind=8
      write(lname,'("Aquifer moisture content ",I1.1)') n
      write(vname,'("t",I1.1,"_GWwb")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3e5,0,2) ! kind=8
      write(lname,'("cansto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_cansto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,0,2) ! kind=8
      write(lname,'("us tile ",I1.1)') n
      write(vname,'("t",I1.1,"_us")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m/s',0.,13.,0,2) ! kind=8   
      write(lname,'("pudsto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_pudsto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,0,2) ! kind=8
      write(lname,'("wetfac tile ",I1.1)') n
      write(vname,'("t",I1.1,"_wetfac")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,0,2) ! kind=8
      write(lname,'("ga tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ga")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',-6500.,6500.,0,2) ! kind=8
    end do
    if ( soil_struc==1 ) then
      do n=1,maxtile  
        write(lname,'("hzero tile ",I1.1)') n
        write(vname,'("t",I1.1,"_hzero")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        do k=1,ms
          write(lname,'("S tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_s",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k=1,ms
          write(lname,'("tsoil tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k=1,ms
          write(lname,'("thetai tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        write(lname,'("snowliq tile ",I1.1," lev ",I1.1)') n,1
        write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("tsurface tile ",I1.1)') n
        write(vname,'("t",I1.1,"_tsurface")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("nsnow tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsnow")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("fwsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fwsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
      end do
    end if
    if ( ccycle==0 ) then
      !do n=1,maxtile  
        !write(lname,'("Carbon leaf pool tile ",I1.1)') n
        !write(vname,'("t",I1.1,"_cplant1")') n    
        !call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        !write(lname,'("Carbon wood pool tile ",I1.1)') n
        !write(vname,'("t",I1.1,"_cplant2")') n
        !call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        !write(lname,'("Carbon root pool tile ",I1.1)') n
        !write(vname,'("t",I1.1,"_cplant3")') n
        !call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        !write(lname,'("Carbon soil fast pool tile ",I1.1)') n
        !write(vname,'("t",I1.1,"_csoil1")') n
        !call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        !write(lname,'("Carbon soil slow pool tile ",I1.1)') n
        !write(vname,'("t",I1.1,"_csoil2")') n
        !call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
      !end do  
    else
      do n=1,maxtile  
        do k = 1,mplant
          write(lname,'("C leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,mplant
          write(lname,'("N leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,mplant
          write(lname,'("P leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,mlitter
          write(lname,'("C litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,mlitter
          write(lname,'("N litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,mlitter
          write(lname,'("P litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,msoil
          write(lname,'("C soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,msoil
          write(lname,'("N soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        do k = 1,msoil
          write(lname,'("P soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        write(lname,'("Prognostic LAI tile ",I1.1)') n
        write(vname,'("t",I1.1,"_glai")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Leaf phenology phen tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Leaf phenology rainfall ",I1.1)') n
        write(vname,'("t",I1.1,"_aphen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Leaf phenology phase tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phenphase")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Leaf phenology doyphase3 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_doyphase3")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("C labile tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabile")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("N soilmin tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsoilmin")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("P soillab tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoillab")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("P soilsorb tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilsorb")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("P soilocc tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilocc")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        do k = 1,mplant
          write(lname,'("crmplant tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        end do
        write(lname,'("frac_sapwood tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fracsapwood")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("sapwoodarea tile ",I1.1)') n
        write(vname,'("t",I1.1,"_sapwoodarea")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("crsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("cnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_cnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("clabloss tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabloss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("crgplant tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crgplant")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("stemnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_stemnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("LAImax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_LAImax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Cleafmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Cleafmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("Crootmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Crootmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("fpn ",I1.1)') n
        write(vname,'("t",I1.1,"_fpn")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("frday ",I1.1)') n
        write(vname,'("t",I1.1,"_frday")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
      end do
    end if
    if ( cable_climate==1 ) then
      do n=1,maxtile  
        write(lname,'("climate iveg tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateiveg")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate biome tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatebiome")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate evapPT tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateevapPT")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate aevap tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateaevap")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtempmax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtempmax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtempmin tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtempmin")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoistmax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoistmax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoistmin tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoistmin")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate qtempmax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateqtempmax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate qtempmaxlastyear tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate gdd0 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climategdd0")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate gdd5 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climategdd5")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate agdd0 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateagdd0")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate agdd5 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateagdd5")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate chilldays tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatechilldays")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate gdd0_rec tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climategdd0_rec")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtemp_min20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtempmin20")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtemp_max20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtempmax20")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate alpha_PT20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatealphaPT20")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoist_min20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoistmin20")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoist_max20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoistmax20")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate frec tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatefrec")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate fdorm tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatefdorm")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate gmd tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climategmd")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dtemp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedtemp")') n
        call attrib(idnc,c3dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoist tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoist")') n
        call attrib(idnc,c4dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtemp_min_20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtemp_min_20")') n
        call attrib(idnc,c5dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate mtemp_max_20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatemtemp_max_20")') n
        call attrib(idnc,c5dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate alpha_PT_20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatealpha_PT_20")') n
        call attrib(idnc,c5dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoist_min_20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoist_min_20")') n
        call attrib(idnc,c5dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate dmoist_max_20 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedmoist_max_20")') n
        call attrib(idnc,c5dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate APAR_leaf_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateapar_leaf_sun")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate APAR_leaf_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climateapar_leaf_shade")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate Dleaf_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedleaf_sun")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate fwsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatefwsoil")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate Dleaf_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatedleaf_shade")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate Tleaf_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatetleaf_sun")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate Tleaf_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatetleaf_shade")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate cs_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatecs_sun")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate cs_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatecs_shade")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate scalex_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatescalex_sun")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("climate scalex_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_climatescalex_shade")') n
        call attrib(idnc,c6dim,c1size,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("veg vcmax_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_vegvcmax_sun")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("veg vcmax_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_vegvcmax_shade")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("veg ejmax_sun tile ",I1.1)') n
        write(vname,'("t",I1.1,"_vegejmax_sun")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
        write(lname,'("veg ejmax_shade tile ",I1.1)') n
        write(vname,'("t",I1.1,"_vegejmax_shade")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,0,2) ! kind=8
      end do  
    end if
  end if
  if ( cable_pop==1 ) then
    do n = 1,maxtile  
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( diaglevel_pop > 0 ) then
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
      end if
      if ( itype==-1 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_densindiv")') n
        write(vname,'("t",I1.1,"_pop_grid_densindiv")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_height_mean")') n
        write(vname,'("t",I1.1,"_pop_grid_height_mean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_height_max")') n
        write(vname,'("t",I1.1,"_pop_grid_height_max")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_basal_area")') n
        write(vname,'("t",I1.1,"_pop_grid_basal_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_res_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_area_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_crown_cover")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_crown_area")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_crown_volume")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_KClump")') n
        write(vname,'("t",I1.1,"_pop_grid_KClump")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_freq_age")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
        call attrib(idnc,c7dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_biomass_age")') n
        write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
        call attrib(idnc,c7dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do      
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do         
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do         
        do dd = 1,POP_NDISTURB
          write(lname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do         
      end if
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_patch_id")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_freq")') n
        write(vname,'("t",I1.1,"_pop_grid_freq")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_freq_old")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_lai")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          call attrib(idnc,c1dim,c1size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          call attrib(idnc,c2dim,c2size,vname,lname,'none',0.,6500.,0,2) ! kind=8
        end do
      end if
    end do  
  end if   
  if ( itype==-1 ) then !just for restart file
    !do n=1,maxtile
      !write(lname,'("Sensible correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fhscor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,0,2) ! kind=8
      !write(lname,'("Latent correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fescor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,0,2) ! kind=8
    !end do
    lname='DIR VIS albedo'
    vname='albvisdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
    lname='DIF VIS albedo'
    vname='albvisdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
    lname='DIR NIR albedo'
    vname='albnirdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
    lname='DIF NIR albedo'
    vname='albnirdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
    lname='VIS albedo'
    vname='albvis'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
    lname='NIR albedo'
    vname='albnir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,0,2) ! kind=8
  end if
end if
  
return
end subroutine savetiledef

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,iarch,itype)

use carbpools_m
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
integer, intent(in) :: idnc,iarch,itype
integer k,n
integer cc,dd,hh,ll
real, dimension(ifull) :: datr
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(mp_global) :: dummy_unpack
real(kind=8), dimension(:,:), allocatable, save :: datpatch
real(kind=8), dimension(:,:), allocatable, save :: datage
real(kind=8), dimension(:,:,:), allocatable, save :: datpc
real(kind=8), dimension(:,:), allocatable, save :: dat91days
real(kind=8), dimension(:,:), allocatable, save :: dat31days
real(kind=8), dimension(:,:), allocatable, save :: dat20years
real(kind=8), dimension(:,:), allocatable, save :: dat5days
character(len=80) vname
logical, intent(in) :: local
  
if ( itype==-1 ) then !just for restart file
  do n = 1,maxtile  ! tile
    datr = 0.
    if ( n<=maxnb ) call cable_unpack(sv,datr,n)
    write(vname,'("t",I1.1,"_svs")') n
    call histwrt(datr,vname,idnc,iarch,local,.true.)      
    do k = 1,ms     ! soil layer
      dat = real(tgg(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%tgg(:,k),dat,n)
      write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,ms
      dat = real(wb(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%wb(:,k),dat,n)
      write(vname,'("t",I1.1,"_wb",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,ms
      dat = real(wbice(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%wbice(:,k),dat,n)
      write(vname,'("t",I1.1,"_wbice",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3 ! snow layer
      dat = real(tggsn(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%tggsn(:,k),dat,n)
      write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(smass(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%smass(:,k),dat,n)
      write(vname,'("t",I1.1,"_smass",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(ssdn(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%ssdn(:,k),dat,n)
      write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do  
    do k = 1,3
      dat = real(snowd/3.,8)
      if ( n<=maxnb ) call cable_unpack(ssnow%sdepth(:,k),dat,n)
      write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = 0.2_8
      if ( n<=maxnb ) call cable_unpack(ssnow%sconds(:,k),dat,n)
      write(vname,'("t",I1.1,"_sconds",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    dat = real(ssdnn,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%ssdnn,dat,n)
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(isflag,8)
    if ( n<=maxnb ) then
      dummy_unpack = real(ssnow%isflag,8)    
      call cable_unpack(dummy_unpack,dat,n)
    end if    
    write(vname,'("t",I1.1,"_sflag")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%snowd,dat,n)
    write(vname,'("t",I1.1,"_snd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%osnowd,dat,n)
    write(vname,'("t",I1.1,"_osnd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snage,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%snage,dat,n)
    write(vname,'("t",I1.1,"_snage")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 100._8
    if ( n<=maxnb ) call cable_unpack(ssnow%rtsoil,dat,n)
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 0._8
    if ( n<=maxnb ) call cable_unpack(ssnow%GWwb,dat,n)
    write(vname,'("t",I1.1,"_GWwb")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(canopy%cansto,dat,n)
    write(vname,'("t",I1.1,"_cansto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0.01_8 ! ustar
    if (n<=maxnb) call cable_unpack(canopy%us,dat,n)
    write(vname,'("t",I1.1,"_us")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)  
    dat=0._8
    if (n<=maxnb) call cable_unpack(ssnow%pudsto,dat,n)
    write(vname,'("t",I1.1,"_pudsto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(ssnow%wetfac,dat,n)
    write(vname,'("t",I1.1,"_wetfac")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(canopy%ga,dat,n)
    write(vname,'("t",I1.1,"_ga")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do
  if ( soil_struc==1 ) then
    do n = 1,maxtile  ! tile
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%h0,dat,n)
      write(vname,'("t",I1.1,"_hzero")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      do k = 1,ms     ! soil layer
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%S(:,k),dat,n)
        write(vname,'("t",I1.1,"_s",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,ms
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%tsoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,ms
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%thetai(:,k),dat,n)
        write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%snowliq(:,1),dat,n)
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%tsurface,dat,n)
      write(vname,'("t",I1.1,"_tsurface")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(ssnow%nsnow,8)  
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_nsnow")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%fwsoil,dat,n)
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.) 
    end do
  end if
  if ( ccycle/=0 ) then
    do n = 1,maxtile  ! tile
      do k = 1,mplant     
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%cplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%pplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%clitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nlitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do  
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%plitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%csoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nsoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%psoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(casamet%glai(:),dat,n)
      write(vname,'("t",I1.1,"_glai")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(phen%phen(:),dat,n)
      write(vname,'("t",I1.1,"_phen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(phen%aphen(:),dat,n)
      write(vname,'("t",I1.1,"_aphen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(phen%phase(:),8)   
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_phenphase")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(phen%doyphase(:,3),8)    
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%clabile(:),dat,n)
      write(vname,'("t",I1.1,"_clabile")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%nsoilmin(:),dat,n)
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoillab(:),dat,n)
      write(vname,'("t",I1.1,"_psoillab")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoilsorb(:),dat,n)
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoilocc(:),dat,n)
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do k = 1,mplant     
        dat=0._8
        if (n<=maxnb) call cable_unpack(casaflux%crmplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%frac_sapwood(:),dat,n)
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%sapwood_area(:),dat,n)
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%Crsoil(:),dat,n)
      write(vname,'("t",I1.1,"_crsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%cnpp(:),dat,n)
      write(vname,'("t",I1.1,"_cnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%clabloss(:),dat,n)
      write(vname,'("t",I1.1,"_clabloss")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%crgplant(:),dat,n)
      write(vname,'("t",I1.1,"_crgplant")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%stemnpp(:),dat,n)
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%laimax(:),dat,n)
      write(vname,'("t",I1.1,"_LAImax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%cleafmean(:),dat,n)
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%crootmean(:),dat,n)
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%fpn(:),dat,n)
      write(vname,'("t",I1.1,"_fpn")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%frday(:),dat,n)
      write(vname,'("t",I1.1,"_frday")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
  end if
  if ( cable_climate==1 ) then
    allocate( dat91days(ifull,91) )
    allocate( dat31days(ifull,31) )
    allocate( dat20years(ifull,20) )
    allocate( dat5days(ifull,120) )
    dat91days = 0._8
    dat31days = 0._8
    dat20years = 0._8
    dat5days = 0._8
    do n = 1,maxtile  ! tile
      dat=0._8
      if (n<=maxnb) dummy_unpack = real(climate%iveg,8)
      if (n<=maxnb) call cable_unpack(dummy_unpack,dat,n)
      write(vname,'("t",I1.1,"_climateiveg")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) dummy_unpack = real(climate%biome,8)
      if (n<=maxnb) call cable_unpack(dummy_unpack,dat,n)
      write(vname,'("t",I1.1,"_climatebiome")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%evap_PT,dat,n)
      write(vname,'("t",I1.1,"_climateevapPT")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%aevap,dat,n)
      write(vname,'("t",I1.1,"_climateaevap")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%mtemp_max,dat,n)
      write(vname,'("t",I1.1,"_climatemtempmax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%mtemp_min,dat,n)
      write(vname,'("t",I1.1,"_climatemtempmin")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%dmoist_max,dat,n)
      write(vname,'("t",I1.1,"_climatedmoistmax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%dmoist_min,dat,n)
      write(vname,'("t",I1.1,"_climatedmoistmin")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%qtemp_max,dat,n)
      write(vname,'("t",I1.1,"_climateqtempmax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%qtemp_max_last_year,dat,n)
      write(vname,'("t",I1.1,"_climateqtempmaxlastyear")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%gdd0,dat,n)
      write(vname,'("t",I1.1,"_climategdd0")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%gdd5,dat,n)
      write(vname,'("t",I1.1,"_climategdd5")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%agdd0,dat,n)
      write(vname,'("t",I1.1,"_climateagdd0")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%agdd5,dat,n)
      write(vname,'("t",I1.1,"_climateagdd5")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(climate%chilldays,8)  
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_climatechilldays")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%gdd0_rec,dat,n)
      write(vname,'("t",I1.1,"_climategdd0_rec")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%mtemp_min20,dat,n)
      write(vname,'("t",I1.1,"_climatemtempmin20")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%mtemp_max20,dat,n)
      write(vname,'("t",I1.1,"_climatemtempmax20")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%alpha_PT20,dat,n)
      write(vname,'("t",I1.1,"_climatealphaPT20")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%dmoist_min20,dat,n)
      write(vname,'("t",I1.1,"_climatedmoistmin20")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%dmoist_max20,dat,n)
      write(vname,'("t",I1.1,"_climatedmoistmax20")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%frec,dat,n)
      write(vname,'("t",I1.1,"_climatefrec")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(climate%fdorm,dat,n)
      write(vname,'("t",I1.1,"_climatefdorm")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(climate%gmd,8)  
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_climategmd")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat91days = 0._8
      if ( n<=maxnb ) then
        do k = 1,91  
          call cable_unpack(climate%dtemp_91(:,k),dat91days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedtemp")') n
      call histwrt(dat91days,vname,idnc,iarch,local,.true.)
      dat31days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,31  
          call cable_unpack(climate%dmoist_31(:,k),dat31days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedmoist")') n
      call histwrt(dat31days,vname,idnc,iarch,local,.true.)
      dat20years = 0._8  
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_unpack(climate%mtemp_min_20(:,k),dat20years(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatemtemp_min_20")') n
      call histwrt(dat20years,vname,idnc,iarch,local,.true.)
      dat20years = 0._8  
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_unpack(climate%mtemp_max_20(:,k),dat20years(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatemtemp_max_20")') n
      call histwrt(dat20years,vname,idnc,iarch,local,.true.)
      dat20years = 0._8  
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_unpack(climate%alpha_PT_20(:,k),dat20years(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatealpha_PT_20")') n
      call histwrt(dat20years,vname,idnc,iarch,local,.true.)
      dat20years = 0._8  
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_unpack(climate%dmoist_min_20(:,k),dat20years(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedmoist_min_20")') n
      call histwrt(dat20years,vname,idnc,iarch,local,.true.)
      dat20years = 0._8  
      if ( n<=maxnb ) then
        do k = 1,20  
          call cable_unpack(climate%dmoist_max_20(:,k),dat20years(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedmoist_max_20")') n
      call histwrt(dat20years,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%APAR_leaf_sun(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climateapar_leaf_sun")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%APAR_leaf_shade(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climateapar_leaf_shade")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%Dleaf_sun(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedleaf_sun")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%fwsoil(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatefwsoil")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%Dleaf_shade(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatedleaf_shade")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%Tleaf_sun(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatetleaf_sun")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%Tleaf_shade(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatetleaf_shade")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%cs_sun(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatecs_sun")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%cs_shade(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatecs_shade")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%scalex_sun(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatescalex_sun")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat5days = 0._8  
      if ( n<=maxnb ) then
        do k = 1,120  
          call cable_unpack(climate%scalex_shade(:,k),dat5days(:,k),n)
        end do  
      end if  
      write(vname,'("t",I1.1,"_climatescalex_shade")') n
      call histwrt(dat5days,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(veg%vcmax_sun,dat,n)
      write(vname,'("t",I1.1,"_vegvcmax_sun")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(veg%vcmax_shade,dat,n)
      write(vname,'("t",I1.1,"_vegvcmax_shade")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(veg%ejmax_sun,dat,n)
      write(vname,'("t",I1.1,"_vegejmax_sun")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(veg%ejmax_shade,dat,n)
      write(vname,'("t",I1.1,"_vegejmax_shade")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do  
    deallocate( dat91days )
    deallocate( dat31days )
    deallocate( dat20years )
    deallocate( dat5days )
  end if
end if
if ( cable_pop==1 ) then
  allocate( datpatch(ifull,POP_NPATCH) )  
  allocate( datage(ifull,POP_AGEMAX) )  
  allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
  do n = 1,maxtile
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    dat = 0._8
    if ( diaglevel_pop > 0 ) then
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%cmass_sum,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end if
    if ( itype==-1 ) then !just for restart file
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%cmass_sum_old,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%cheartwood_sum,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%csapwood_sum,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%csapwood_sum_old,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%densindiv,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%height_mean,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%height_max,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%basal_area,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%sapwood_loss,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%sapwood_area_loss,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%stress_mortality,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%crowding_mortality,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%fire_mortality,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%cat_mortality,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%res_mortality,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%growth,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%area_growth,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%crown_cover,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%crown_area,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%crown_volume,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%sapwood_area,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%sapwood_area_old,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%KClump,dat,n)
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do k = 1,POP_AGEMAX
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%freq_age(k),datage(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do k = 1,POP_AGEMAX
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%biomass_age(k),datage(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%biomass(ll),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%density(ll),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%hmean(ll),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%hmax(ll),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%cmass_stem_bin(hh),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%densindiv_bin(hh),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%height_bin(hh),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%diameter_bin(hh),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%n_age(dd),dat,n)
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
    end if
    if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%id,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%freq(k),datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%freq_old(k),datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%factor_recruit,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%pgap,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%lai,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%biomass,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%biomass_old,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)    
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%heartwood,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood_old,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood_area,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood_area_old,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%stress_mortality,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%fire_mortality,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%cat_mortality,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%crowding_mortality,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%cpc,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%mortality,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood_loss,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%sapwood_area_loss,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%growth,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%area_growth,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%frac_NPP,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%frac_respiration,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%frac_light_uptake,datpatch(:,k),n)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%disturbance_interval(dd),datpatch(:,k),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%first_disturbance_year(dd),datpatch(:,k),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%age(dd),datpatch(:,k),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%ranked_age_unique(k,dd),datpatch(:,k),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%freq_ranked_age_unique(k,dd),datpatch(:,k),n)
        end do  
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%ncohort,datpatch(:,k),n)
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%biomass,datpatch(:,k),n)
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%density,datpatch(:,k),n)
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%hmean,datpatch(:,k),n)
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%hmax,datpatch(:,k),n)
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT  
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT  
          do k = 1,POP_NPATCH  
           if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP,datpc(:,k,cc),n)
          end do
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT  
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter,datpc(:,k,cc),n)
          end do 
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            if ( n<=maxnb ) call pop_unpack(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot,datpc(:,k,cc),n)
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
    end if
  end do
  deallocate( datpatch )
  deallocate( datage )
  deallocate( datpc )
end if    
if ( itype==-1 ) then !just for restart file
  !do n = 1,maxtile  ! tile
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fhs_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fhscor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fes_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fescor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fns_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fnscor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%ga_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_gacor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
  !end do
  vname='albvisdir'
  call histwrt(albvisdir,vname,idnc,iarch,local,.true.)
  vname='albvisdif'
  call histwrt(albvisdif,vname,idnc,iarch,local,.true.)
  vname='albnirdir'
  call histwrt(albnirdir,vname,idnc,iarch,local,.true.)
  vname='albnirdif'
  call histwrt(albnirdif,vname,idnc,iarch,local,.true.)
  vname='albvis'
  call histwrt(albvissav,vname,idnc,iarch,local,.true.)
  vname='albnir'
  call histwrt(albnirsav,vname,idnc,iarch,local,.true.)
end if
  
return
end subroutine savetile 

! *************************************************************************************
! Water inflow from river routing
subroutine cableinflow(inflow)

use newmpar_m
use soil_m

implicit none

integer k
real, dimension(ifull), intent(inout) :: inflow
real, dimension(ifull) :: delflow
real, dimension(mp_global) :: xx, ll, delxx
real, dimension(mp_global) :: dummy_unpack

if ( mp_global<=0 ) return

call cable_pack( inflow(1:ifull), xx )
delxx(1:mp_global) = 0.
do k = 1,cbm_ms
  ll(1:mp_global) = real( max( soil%sfc(1:mp_global)-ssnow%wb(1:mp_global,k), 0._8 )*1000._8*soil%zse(k) )
  ll(1:mp_global) = min( xx(1:mp_global), ll(1:mp_global) )
  ssnow%wb(1:mp_global,k) = ssnow%wb(1:mp_global,k) + real(ll(1:mp_global),8)/(1000._8*soil%zse(k))
  delxx(1:mp_global) = delxx(1:mp_global) - ll(1:mp_global)
end do
delflow(1:ifull) = 0.
dummy_unpack = sv(1:mp_global)*delxx(1:mp_global)
call cable_unpack(dummy_unpack,delflow(1:ifull))
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
real, dimension(:,:), allocatable, save :: dumg
real, dimension(ifull,5), intent(out) :: casapoint
character(len=*), intent(in) :: casafile
logical tst

if ( myid==0 ) then
  allocate( dumg(ifull_g,5) )
  dumg = 0.
  if ( casafile==" " ) then
    write(6,*) "ERROR: casafile is not specified"
    call ccmpi_abort(-1)
  end if    
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
use infile
use newmpar_m

implicit none

integer, parameter :: nphen = 8
integer ilat, ierr
integer ncid
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
integer, dimension(nphen) :: greenupx, fallx, xphendoy1
integer, dimension(nphen) :: ivtx
real :: xlat
logical :: ncfile
character(len=*), intent(in) :: fphen

! initilize for evergreen PFTs
greenup(:,:)  = -50
fall(:,:)     = 367
phendoy1(:,:) = 2

! MJT notes - if the CSIRO PFTs are redefined, then this phenology data will be mismatched

if ( myid==0 ) then
  write(6,*) "Reading CASA leaf phenology data"
  call ccnf_open(fphen,ncid,ierr)
  if ( ierr==0 ) then
    ncfile = .true.
    write(6,*) "Found netcdf file ",trim(fphen)
  else
    call ccnf_open(fphen//'.nc',ncid,ierr)
    if ( ierr==0 ) then
      ncfile = .true.
      write(6,*) "Found netcdf file ",trim(fphen)
    else  
      open(87,file=fphen,status='old',iostat=ierr)
      if ( ierr/=0 ) then
        write(6,*) "ERROR: Cannot open phenfile=",trim(fphen)
        call ccmpi_abort(-1)
      end if
      read(87,*)
      read(87,*) ivtx
    end if
  end if
  
  if ( ncfile ) then
    call ccnf_get_vara(ncid,"greenup",(/1,1/),(/271,mxvt/),greenup)
    call ccnf_get_vara(ncid,"fall",(/1,1/),(/271,mxvt/),fall)
    call ccnf_get_vara(ncid,"phendoy1",(/1,1/),(/271,mxvt/),phendoy1)
    call ccnf_close(ncid)
  else    
    do ilat = 271,1,-1
      read(87,*) xlat,greenupx,fallx,xphendoy1 
      greenup(ilat,ivtx(:))  = greenupx(:)
      fall(ilat,ivtx(:))     = fallx(:)
      phendoy1(ilat,ivtx(:)) = xphendoy1(:)
    end do
    close(87)
  end if  
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

ssnow%tss = real(tset(:),8)
do k = 1,3
  ssnow%tggsn(:,k) = real(tset(:),8)
end do
do k = 1,ms
  ssnow%tgg(:,k) = real(tset(:),8)
end do

return
end subroutine cablesettemp

end module cable_ccam

