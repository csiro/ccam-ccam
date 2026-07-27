! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! Subsequently modified by MJT for tile mosaic and SEAESF radiation scheme
    
! - Currently all tiles within a gridbox have the same soil texture, but independent soil
!   temperatures, moisture, etc.
! - LAI can be prescribed or using CASA-CNP can predict LAI and vcmax, but requires
!   considerable time to spin-up.
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.
! - Options exist for using SLI soil model (soil_struc=1), the POP model (cable_pop=1)
!   and the Ground Water model (cable_gw_model=1).  POPLUC has not yet been implemented.

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

use cable_ccam_common, cbm_mp => mp
use cable_ccam_file
use cable_ccam_setup

implicit none

private
public sib4
public newcbmwb, cbmemiss, cable_casatile

! Imported from cable_ccam_common
public proglai, progvcmax, maxtile, soil_struc, cable_pop, ccycle, cable_potev
public fwsoil_switch, cable_litter, gs_switch, cable_enablefao
public smrf_switch, strf_switch, cable_gw_model, cable_roughness
public cable_runoff, cable_soilevap
public POP_NPATCH, POP_NCOHORT, POP_AGEMAX
public mplant, mlitter, msoil

! Imported from cable_ccam_file
public loadtile, defaulttile, savetiledef, savetile

! Imported from cable_ccam_setup
public loadcbmparm, cbmparm

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
use raddiag_m
use radisw_m
use riverarrays_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use tracers_m
use vegpar_m
use work2_m, only : qsttg,zo,zoh,zoq,theta,vmod,wetfac
use work3_m, only : ga
use zenith_m

integer ico2, igas
integer iyr, imo, iday, mdays, k
real x
integer tile, is, ie
real, dimension(imax,mlitter) :: lclitter, lnilitter, lplitter
real, dimension(imax,mplant) :: lcplant, lniplant, lpplant
real, dimension(imax,msoil) :: lcsoil, lnisoil, lpsoil
real, dimension(imax,ms) :: lwb_clim, ltgg, lwb, lwbice
real, dimension(imax,3) :: lsmass, lssdn, ltggsn
real, dimension(imax) :: lcnbp, lcnpp, lfnee, lfpn, lfrd, lfrp, lfrpr, lfrpw, lfrs
real, dimension(imax) :: latmco2
real, dimension(imax) :: lfevc, lplant_turnover, lplant_turnover_wood

!$omp do schedule(static) private(is,ie,ico2,igas,x,k,iyr,imo,iday)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  if ( tdata(tile)%mp>0 ) then
      
    cbm_mp = tdata(tile)%mp
  
    ltgg(:,:) = tgg(is:ie,:)
    ltggsn(:,:) = tggsn(is:ie,:)
    lwb(:,:) = wb(is:ie,:)
    lwbice(:,:) = wbice(is:ie,:)
    lssdn(:,:) = ssdn(is:ie,:)
    lsmass(:,:) = smass(is:ie,:)
    if ( ccycle/=0 ) then
      lcnbp(:) = cnbp(is:ie)
      lcnpp(:) = cnpp(is:ie)
      lcplant(:,:) = cplant(is:ie,:)
      lclitter(:,:) = clitter(is:ie,:)
      lcsoil(:,:) = csoil(is:ie,:)
      lniplant(:,:) = niplant(is:ie,:)
      lnilitter(:,:) = nilitter(is:ie,:)
      lnisoil(:,:) = nisoil(is:ie,:)
      lpplant(:,:) = pplant(is:ie,:)
      lplitter(:,:) = plitter(is:ie,:)
      lpsoil(:,:) = psoil(is:ie,:)
      lfnee(:) = fnee(is:ie)
      lfpn(:) = fpn(is:ie)
      lfrd(:) = frd(is:ie)
      lfrp(:) = frp(is:ie)
      lfrpr(:) = frpr(is:ie)
      lfrpw(:) = frpw(is:ie)
      lfrs(:) = frs(is:ie)
      if ( diaglevel_carbon>0 ) then
        lfevc(:) = fevc(is:ie)
        lplant_turnover(:) = plant_turnover(is:ie)
        lplant_turnover_wood(:) = plant_turnover_wood(is:ie)
      end if
    end if 
    if ( nspecial == 51 ) then
      call time_of_month(iyr,imo,iday,mdays,kdate,mtimer,leap)
      x = (real(iday)-0.5)/real(mdays)
      do k = 1,ms
        call time_interpolate(lwb_clim(:,k),wb_clim(is:ie,k,1),wb_clim(is:ie,k,2),      &
                              wb_clim(is:ie,k,3),wb_clim(is:ie,k,4),wb_clim(is:ie,k,5), &
                              x,24)
      end do  
    end if
    
    ! set co2 forcing for cable
    ! host: atmospheric co2 follows that from CCAM radiation scheme
    ! interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)
    latmco2(:) = 1.E6*rrvco2          ! from radiative CO2 forcings
    if ( ngas>0 ) then
      ico2 = 0
      do igas = 1,ngas
        if ( trim(tractype(igas))=='online' .and. trim(tracname(igas))=='cbmnep' ) then
          ico2 = igas
          exit
        end if
      end do
      if ( ico2>0 ) then
        latmco2(:) = tr(is:ie,1,ico2) ! use interactive tracers
      end if
    end if

    call sib4_work(albnirdif(is:ie),albnirdir(is:ie),albnirsav(is:ie),albvisdif(is:ie),albvisdir(is:ie),                        &
                   albvissav(is:ie),cansto(is:ie),cdtq(is:ie),cduv(is:ie),lclitter,lcnbp,                                       &
                   lcnpp,condg(is:ie),conds(is:ie),condx(is:ie),lcplant,lcsoil,eg(is:ie),epot(is:ie),fbeamnir(is:ie),           &
                   fbeamvis(is:ie),fg(is:ie),lfnee,lfpn,lfrd,lfrp,lfrpr,lfrpw,lfrs,fwet(is:ie),ga(is:ie),                       &
                   isflag(is:ie),land(is:ie),lnilitter,lniplant,lnisoil,lplitter,lpplant,                                       &
                   ps(is:ie),lpsoil,qg(is:ie,1),qsttg(is:ie),rgsave(is:ie),rlatt(is:ie),rlongg(is:ie),rnet(is:ie),rsmin(is:ie), &
                   runoff(is:ie),runoff_surface(is:ie),sigmf(is:ie),lsmass,snage(is:ie),snowd(is:ie),snowmelt(is:ie),           &
                   lssdn,ssdnn(is:ie),sgdn(is:ie),swrsave(is:ie),t(is:ie,1),ltgg,ltggsn,theta(is:ie),latmco2,                   & 
                   tss(is:ie),ustar(is:ie),vlai(is:ie),vmod(is:ie),lwb,lwbice,wetfac(is:ie),zo(is:ie),zoh(is:ie),zoq(is:ie),    &
                   evspsbl(is:ie),sbl(is:ie),wtd(is:ie),air(tile),bal(tile),canopy(tile),casabal(tile),casabiome(tile),         &
                   casaflux(tile),casamet(tile),casapool(tile),climate(tile),met(tile),phen(tile),pop(tile),rad(tile),          &
                   rough(tile),soil(tile),ssnow(tile),sum_flux(tile),veg(tile),lfevc,lplant_turnover,lplant_turnover_wood,      &
                   lwb_clim,tile)
 
    tgg(is:ie,:) = ltgg(:,:)
    tggsn(is:ie,:) = ltggsn(:,:)
    wb(is:ie,:) = lwb(:,:)
    wbice(is:ie,:) = lwbice(:,:)
    ssdn(is:ie,:) = lssdn(:,:)
    smass(is:ie,:) = lsmass(:,:)
    if ( ccycle/=0 ) then
      cnbp(is:ie) = lcnbp(:)
      cnpp(is:ie) = lcnpp(:)
      cplant(is:ie,:) = lcplant(:,:)
      clitter(is:ie,:) = lclitter(:,:)
      csoil(is:ie,:) = lcsoil(:,:)
      niplant(is:ie,:) = lniplant(:,:)
      nilitter(is:ie,:) = lnilitter(:,:)
      nisoil(is:ie,:) = lnisoil(:,:)
      pplant(is:ie,:) = lpplant(:,:)
      plitter(is:ie,:) = lplitter(:,:)
      psoil(is:ie,:) = lpsoil(:,:)
      fnee(is:ie) = lfnee(:)
      fpn(is:ie) = lfpn(:)
      frd(is:ie) = lfrd(:)
      frp(is:ie) = lfrp(:)
      frpr(is:ie) = lfrpr(:)
      frpw(is:ie) = lfrpw(:)
      frs(is:ie) = lfrs(:)
      if ( diaglevel_carbon > 0 ) then
        fevc(is:ie) = lfevc(:)
        plant_turnover(is:ie) = lplant_turnover(:)
        plant_turnover_wood(is:ie) = lplant_turnover_wood(:)
      end if
    end if
    
  end if  

end do
!$omp end do nowait

return
end subroutine sib4

! CABLE-CCAM interface
subroutine sib4_work(albnirdif,albnirdir,albnirsav,albvisdif,albvisdir,albvissav,cansto,cdtq,cduv,clitter,cnbp, &
                     cnpp,condg,conds,condx,cplant,csoil,eg,epot,fbeamnir,fbeamvis,fg,fnee,fpn,frd,frp,frpr,    &
                     frpw,frs,fwet,ga,isflag,land,nilitter,niplant,nisoil,plitter,pplant,ps,                    &
                     psoil,qg,qsttg,rgsave,rlatt,rlongg,rnet,rsmin,runoff,runoff_surface,sigmf,smass,           &
                     snage,snowd,snowmelt,ssdn,ssdnn,sgdn,swrsave,t,tgg,tggsn,theta,atmco2,tss,                 &
                     ustar,vlai,vmod,wb,wbice,wetfac,zo,zoh,zoq,evspsbl,sbl,wtd,air,bal,canopy,casabal,         &
                     casabiome,casaflux,casamet,casapool,climate,met,phen,pop,rad,rough,soil,ssnow,sum_flux,    &
                     veg,fevc,plant_turnover,plant_turnover_wood,wb_clim,tile)

use const_phys
use dates_m
use estab, only : qsat
use infile, only : getzinp
use newmpar_m, only : ms, imax, ifull, ntiles
use parm_m
use sigs_m
use soil_m, only : zmin
use zenith_m, only : solargh, zenith
  
integer, intent(in) :: tile
integer jyear, jmonth, jday, jhour, jmin, idoy
integer k, mins, j, nb
integer casaperiod, npercasa, lalloc, mp_POP
integer, dimension(imax), intent(inout) :: isflag
integer, dimension(cbm_mp) :: i_doy
real fjd, r1, dlt, slag, dhr, alp
real, dimension(imax), intent(in) :: atmco2
real, dimension(imax), intent(in) :: qg, t
real, dimension(imax,mplant), intent(inout) :: cplant, niplant, pplant
real, dimension(imax,mlitter), intent(inout) :: clitter, nilitter, plitter
real, dimension(imax,msoil), intent(inout) :: csoil, nisoil, psoil
real, dimension(imax,ms), intent(inout) :: tgg, wb, wbice
real, dimension(imax,ms), intent(in) :: wb_clim
real, dimension(imax,3), intent(inout) :: smass, ssdn, tggsn
real, dimension(imax), intent(inout) :: albnirdif, albnirdir, albnirsav, albvisdif, albvisdir, albvissav
real, dimension(imax), intent(inout) :: cansto, cdtq, cduv, cnbp, cnpp, eg, epot, fg, fnee, fpn, frd
real, dimension(imax), intent(inout) :: frp, frpr, frpw, frs, fwet, ga, qsttg, rnet, rsmin
real, dimension(imax), intent(inout) :: sigmf, snage, snowd, ssdnn, tss, ustar
real, dimension(imax), intent(inout) :: vlai, wetfac, zo, zoh, zoq
real, dimension(imax), intent(inout) :: fevc, plant_turnover, plant_turnover_wood
real, dimension(imax), intent(inout) :: wtd
real, dimension(imax), intent(in) :: condg, conds, condx, fbeamnir, fbeamvis, ps, rgsave, rlatt, rlongg
real, dimension(imax), intent(in) :: sgdn, swrsave, theta, vmod
real, dimension(imax), intent(inout) :: runoff, runoff_surface, snowmelt
real, dimension(imax), intent(inout) :: evspsbl, sbl
real, dimension(imax) :: coszro2, taudar2, tmps, qsttg_land, dummy_pack
real, dimension(imax) :: evspsbl_l, sbl_l
real, dimension(cbm_mp) :: dumtr4, dumpr4, qsatfvar
real, dimension(cbm_mp) :: wbclim_pack
real(kind=8) dtr8
real(kind=8), dimension(cbm_mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
real(kind=8), dimension(cbm_mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
real(kind=8), dimension(cbm_mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
real(kind=8), dimension(cbm_mp) :: xKNlimiting, xkleafcold, xkleafdry
real(kind=8), dimension(cbm_mp) :: xkleaf, xnplimit, xNPuptake, xklitter
real(kind=8), dimension(cbm_mp) :: xksoil
real(kind=8), dimension(cbm_mp) :: tmp, veg_wt, veg_trad, soil_wt, soil_trad, trad_corr
real(kind=8), dimension(cbm_mp) :: dumt, cc1, cc2, ssnowpotev
real(kind=8), dimension(cbm_mp,cbm_nrb) :: c1, rhoch, xk
logical, dimension(imax), intent(in) :: land
logical, dimension(cbm_mp) :: veg_mask, sunlit_mask, sunlit_veg_mask
logical, parameter :: cbl_standalone = .false.
logical, parameter :: jls_standalone = .false.
logical, parameter :: jls_radiation  = .false.
character(len=4), parameter :: subr_name = "ccam"
type(air_type), intent(inout) :: air
type(balances_type), intent(inout) :: bal
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

grid_mp = cbm_mp
dtr8 = real(dt, 8)

! run CASA CNP once per day
casaperiod = nint(86400._8*deltpool)
npercasa = max( nint(real(casaperiod,8)/dtr8), 1 )

! calculate time from beginning of the simulation
! jyear, jmonth, jday, jhour, jmin are for the start date of the simulation
! mins is the time from the start of the year
! mtimer is the time elapsed from the start of the simulation
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
! calculate zenith angle
fjd = real(mins)/1440.
call solargh(fjd,bpyear,r1,dlt,alp,slag)
dhr = dt/3600.
call zenith(fjd,r1,dlt,slag,rlatt(:),rlongg(:),dhr,imax,coszro2(:),taudar2(:))

!! store soil temperature and moisture for budget calculation
!bwb(:,:) = ssnow%wb(:,:)
!btgg(:,:) = ssnow%tgg(:,:)

! checks
if ( any(atmco2(:)<1.) ) then
  write(6,*) "ERROR: Invalid CO2 mixing ratio in cable_ccam2 ",minval(atmco2(:))
  stop
end if

! set meteorological forcing
albvissav(:) = fbeamvis(:)*albvisdir(:) + (1.-fbeamvis(:))*albvisdif(:) ! for nrad=4
albnirsav(:) = fbeamnir(:)*albnirdir(:) + (1.-fbeamnir(:))*albnirdif(:) ! for nrad=4
call cable_pack(theta(:),met%tk,tile)
call cable_pack(vmod(:),met%ua,tile)
dummy_pack(:) = atmco2(:)*1.e-6
call cable_pack(dummy_pack(:),met%ca,tile)
call cable_pack(coszro2(:),met%coszen,tile)
call cable_pack(qg(:),met%qv,tile)
dummy_pack(:) = ps(:)*0.01
call cable_pack(dummy_pack(:),met%pmb,tile)
call cable_pack(condx(:),met%precip,tile)
dummy_pack(:) = conds(:) + condg(:)
call cable_pack(dummy_pack(:),met%precip_sn,tile)
! swrsave indicates the fraction of net VIS radiation (compared to NIR)
! fbeamvis indicates the beam fraction of downwelling direct radiation (compared to diffuse) for VIS
! fbeamnir indicates the beam fraction of downwelling direct radiation (compared to diffuse) for NIR
dummy_pack(:) = swrsave(:)*sgdn(:)
call cable_pack(dummy_pack(:),met%fsd(:,1),tile)
dummy_pack(:) = (1.-swrsave(:))*sgdn(:)
call cable_pack(dummy_pack(:),met%fsd(:,2),tile)
dummy_pack(:) = -rgsave(:)
call cable_pack(dummy_pack(:),met%fld,tile)
call cable_pack(fbeamvis(:),rad%fbeam(:,1),tile)
call cable_pack(fbeamnir(:),rad%fbeam(:,2),tile)
rad%fbeam(:,3) = 0._8 ! dummy for now
dummy_pack(:) = bet(1)*t(:)
call cable_pack(dummy_pack(:),rough%za_tq,tile)
rough%za_tq(:) = max( rough%za_tq(:), veg%hc(:)+1._8 )
met%doy     = real(fjd, 8)
i_doy(:)    = int(met%doy)
met%hod     = (met%doy-int(met%doy))*24._8 + rad%longitude/15._8
met%hod     = mod(met%hod, 24._8)
met%tvair   = met%tk
met%tvrad   = met%tk
met%qvair   = met%qv
met%ua      = max(met%ua, umin)
met%coszen  = max(met%coszen, 1.e-8_8) 
rough%za_uv = rough%za_tq

if ( cable_gw_model==1 ) then
  call cable_pack(wtd(:),ssnow%wtd,tile)
  ssnow%wtd = ssnow%wtd*1000._8
end if

! Interpolate LAI.  Also need sigmf for LDR prognostic aerosols.
call setlai(tdata(tile)%sv,tdata(tile)%vl2,casamet,veg,cbm_mp)

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

! CCAM does not use cable lakes
!WHERE( veg%iveg == lakes_cable .AND. ssnow%wb(:,1) < soil%sfc ) 
!  ssnow%wbtot1(:)  = REAL( ssnow%wb(:,1) ) * density_liq * soil%zse(1)
!  ssnow%wb(:,1) = soil%sfc
!  ssnow%wbtot2  = REAL( ssnow%wb(:,1) ) * density_liq * soil%zse(1)
!ENDWHERE
!ssnow%wb_lake = ssnow%wb_lake + MAX( ssnow%wbtot2 - ssnow%wbtot1, 0.)

call ruff_resist(veg,rough,ssnow,canopy,veg%vlai,veg%hc,canopy%vlaiw)
call define_air(met,air)
call fveg_mask(veg_mask,cbm_mp,lai_thresh,canopy%vlaiw)
tmp(:) = met%fsd(:,1) + met%fsd(:,2)
call fsunlit_mask(sunlit_mask,cbm_mp,rad_thresh,tmp)
call fsunlit_veg_mask(sunlit_veg_mask,veg_mask,sunlit_mask,cbm_mp)
call init_radiation(rad%extkb,rad%extkd,rad%extkbm,rad%extkdm,rad%Fbeam,c1,rhoch,xk, &
                    cbm_mp,cbm_nrb,lai_thresh,coszen_tols,Gauss_w,cbm_pi,cbm_pi180,  &
                    cbl_standalone,jls_standalone,jls_radiation,subr_name,veg_mask,  &
                    veg%Xfang,veg%taul,veg%refl,met%coszen,i_doy,met%fsd,            &
                    canopy%vlaiw)
call albedo(ssnow%AlbSoilsn,soil%AlbSoil,cbm_mp,cbm_nrb,ICE_SoilType,lakes_cable,    &
            jls_radiation,veg_mask,coszen_tols,GAUSS_W,veg%iveg,soil%isoilm,         &
            veg%refl,veg%taul,met%coszen,canopy%vlaiw,ssnow%snowd,ssnow%ssdnn,       &
            ssnow%tgg(:,1),ssnow%snage,xk,c1,rhoch,rad%fbeam,rad%albedo,rad%extkb,   &
            rad%extkd,rad%extkbm,rad%extkdm,rad%rhocbm,rad%rhocdf,rad%cexpkbm,       &
            rad%cexpkdm,rad%reffbm,rad%reffdf)
!rad%albedo_T = (rad%albedo(:, 1) + rad%albedo(:, 2)) * 0.5
call define_canopy(bal,rad,rough,air,met,dtr8,ssnow,soil,veg,canopy,climate,         &
                   sunlit_veg_mask,canopy%vlaiw)
ssnow%otss_0 = ssnow%otss
ssnow%otss = ssnow%tss
ssnow%owetfac = ssnow%wetfac
select case ( soil_struc )
  case(0)  
!    if ( cable_gw_model==1 ) then
!      call soil_snow_gw(dtr8,soil,ssnow,canopy,met,bal,veg)
!    else
       call soil_snow(dtr8,soil,ssnow,canopy,met,bal,veg)
       call snow_aging(ssnow%snage,cbm_mp,dtr8,ssnow%snowd,ssnow%osnowd, &
               ssnow%tggsn(:,1),ssnow%tgg(:,1),ssnow%isflag,             &
               veg%iveg,soil%isoilm) 
!    end if  
  case(1)
     call sli_main(999,dtr8,veg,soil,ssnow,met,canopy,air,rad,0)
  case default
    write(6,*) "ERROR: Unknown option soil_struc ",soil_struc
    stop
end select
! adjust for new soil temperature
ssnow%deltss = ssnow%tss - ssnow%otss
if ( soil_struc==0 ) then
  canopy%fhs = canopy%fhs + ssnow%deltss*ssnow%dfh_dtg
  canopy%fes = canopy%fes + ssnow%deltss*ssnow%dfe_dtg
  canopy%fess = canopy%fess + ssnow%deltss*ssnow%dfe_dtg
  canopy%fns = canopy%fns + ssnow%deltss*ssnow%dfn_dtg
  canopy%ga  = canopy%ga  + ssnow%deltss*canopy%dgdtg
  !! MJT fix
  !ssnow%wb(:,cbm_ms) = ssnow%wb(:,cbm_ms) - (ssnow%deltss*ssnow%dfe_dtg)*dtr8 &
  !                                 /(ssnow%cls*air%rlam*soil%zse(cbm_ms)*1000._8) 
end if
canopy%fh   = canopy%fhv + canopy%fhs
canopy%fev  = canopy%fevc + canopy%fevw
canopy%fe   = canopy%fev + canopy%fes
canopy%rnet = canopy%fns + canopy%fnv

!! Calculate radiative/skin temperature:
!!Jan 2018: UM assumes a single emissivity for the surface in the radiation scheme
!!To accommodate this a single value of is 1. is assumed in ACCESS
!! any leaf/soil emissivity /=1 must be incorporated into rad%trad.  
!! check that emissivities (pft and nvg) set = 1 within the UM i/o configuration
!! CM2 - further adapted to pass the correction term onto %trad correctly
!!rad%trad = ( ( 1.-rad%transd ) * Cemleaf * canopy%tv**4 +                      &
!!       rad%transd * Cemsoil * ssnow%otss**4 + canopy%fns_cor/CSBOLTZ )**0.25
!
!veg_wt(:)    = 1._8 - rad%transd
!veg_trad(:)  = emleaf*canopy%tv**4
!soil_wt(:)   = rad%transd 
!soil_trad(:) = emsoil*ssnow%otss**4
!trad_corr(:) = canopy%fns_cor/SBOLTZ
!
!rad%trad = veg_wt(:)*veg_trad(:) + soil_wt(:)*soil_trad(:) + trad_corr(:)  
!rad%trad = rad%trad**0.25_8

rad%trad = ( (1._8-rad%transd)*canopy%tv**4 + rad%transd*ssnow%tss**4 )**0.25_8

! note that conservation is still preserved at this point
! canopy%ga    = canopy%rnet - canopy%fh - canopy%fe
! canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq*ssnow%ddq_dtg

! MJT suggestion
canopy%cdtq =  max( 0._8, canopy%cdtq )
ssnow%wbice = max( ssnow%wbice, 0._8 )

! change in soil temperature
!deltgg(:) = 0._8
!do k = 1,cbm_ms
!  deltgg(:) = deltgg(:) + (ssnow%tgg(:,k)-btgg(:,k))*ssnow%gammzz(:,k)/dtr8
!end do  
!! energy balance
!tbal = canopy%rnet - canopy%fh - canopy%fe - deltgg
!write(6,*) "energy ",minval(tbal),maxval(tbal)
!! change in soil moisture
!delwb(:) = 0._8
!do k = 1,cbm_ms
!  delwb(:) = delwb(:) + (ssnow%wb(:,k)-bwb(:,k))*soil%zse(k)*1000.
!end do  
!! net water into soil
!bal%wbal = met%precip - canopy%delwc - ssnow%snowd + ssnow%osnowd - ssnow%runoff    &
!          - (canopy%fevw + canopy%fevc + canopy%fes/ssnow%cls)*dtr8/air%rlam - delwb
!write(6,*) "bal%wbal ",minval(bal%wbal),maxval(bal%wbal)

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
if ( any( ssnow%tgg>425. ) ) then
  write(6,*) "WARN: tgg>425. after CABLE"
end if
if ( any( ssnow%tggsn>425. ) ) then
  write(6,*) "WARN: tggsn>425. after CABLE"
end if
if ( any( canopy%tv>425. ) ) then
  write(6,*) "WARN: tv>425. after CABLE" 
end if
if ( any( canopy%tv<100. ) ) then
  write(6,*) "WARN: tv<100. after CABLE" 
end if

!--------------------------------------------------------------
! CASA CNP
select case (ccycle)
  case(0) ! off
    !call plantcarb(veg,bgc,met,canopy)
    !call soilcarb(soil,ssnow,veg,bgc,met,canopy)
    !call carbon_pl(dtr8,soil,ssnow,veg,canopy,bgc)
    !canopy%fnpp = -1.0 * canopy%fpn - canopy%frp
    !canopy%fgpp = -1.0 * canopy%fpn + canopy%frday
    !canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    !canopy%fra = canopy%frp + canopy%frday
  case(1,2,3) ! C, C+N, C+N+P
    ! update casamet
    if ( cable_pop==1 ) then
      lalloc = 3
    else
      lalloc = 0  
    end if
    casamet%tairk = casamet%tairk + met%tk
    casamet%tsoil = casamet%tsoil + ssnow%tgg
    casamet%moist = casamet%moist + ssnow%wb
    casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dtr8
    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dtr8
    if ( mod(ktau,npercasa)==0 ) then
      casamet%tairk=casamet%tairk/real(npercasa,8)
      casamet%tsoil=casamet%tsoil/real(npercasa,8)
      casamet%moist=casamet%moist/real(npercasa,8)
      xKNlimiting(:) = 1._8
      call phenology(idoy,veg,phen)
      call avgsoil(veg,soil,casamet)
      call casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
      if ( cable_pop/=1 ) then
        call casa_allocation(veg,soil,casabiome,casaflux,casapool,       &
                             casamet,phen,lalloc)
      end if
      call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casamet,phen)
      call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool,casaflux,casamet,phen)
      call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)
      if ( cable_pop==1 ) then
        call casa_pop_firstcall(lalloc,casabiome,casaflux,casamet,casapool,phen,pop,soil,veg)
      end if  
      call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
      call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
      if ( ccycle>1 ) then
        call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
        do j = 1,mlitter
          casaflux%klitter(:,j) = casaflux%klitter(:,j)*xkNlimiting(:)
        end do
        call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
        if ( ccycle>2 ) then
          call casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
        end if  
      end if  
      call casa_delplant(veg,casabiome,casapool,casaflux,casamet,      &
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,        &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,        &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
      casaflux%Cplant_turnover_disturbance = 0
      casaflux%Cplant_turnover_crowding = 0
      casaflux%Cplant_turnover_resource_limitation = 0
      if ( cable_pop==1 ) then
        mp_POP = size(POP%pop_grid)  
        call casa_pop_secondcall(mp_POP,casaflux,pop)  
      end if
      call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
      call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet,casabal,lalloc)
      if ( ccycle<2 ) call casa_ndummy(casamet,casabal,casapool)
      if ( ccycle<3 ) call casa_pdummy(casamet,casabal,casaflux,casapool)
      call casa_cnpbal(veg,casamet,casapool,casaflux,casabal)
      casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp*deltpool
      casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp*deltpool
      casabal%FCrmleafyear = casabal%FCrmleafyear + casaflux%Crmplant(:,leaf)*deltpool
      casabal%FCrmwoodyear = casabal%FCrmwoodyear + casaflux%Crmplant(:,wood)*deltpool
      casabal%FCrmrootyear = casabal%FCrmrootyear + casaflux%Crmplant(:,xroot)*deltpool
      casabal%FCrgrowyear  = casabal%FCrgrowyear  + casaflux%Crgplant*deltpool
      casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp*deltpool
      casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil*deltpool
      casabal%FCneeyear = casabal%FCneeyear + (casaflux%Cnpp-casaflux%Crsoil)*deltpool
      casabal%dCdtyear  = casabal%dCdtyear  + (casapool%Ctot-casapool%Ctot_0)*deltpool
      casabal%FNdepyear = casabal%FNdepyear + casaflux%Nmindep*deltpool
      casabal%FNfixyear = casabal%FNfixyear + casaflux%Nminfix*deltpool
      casabal%FNsnetyear = casabal%FNsnetyear + casaflux%Nsnet*deltpool
      casabal%FNupyear = casabal%FNupyear + casaflux%Nminuptake*deltpool
      casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach*deltpool
      casabal%FNlossyear = casabal%FNlossyear + casaflux%Nminloss*deltpool
      casabal%FPweayear = casabal%FPweayear + casaflux%Pwea*deltpool
      casabal%FPdustyear = casabal%FPdustyear + casaflux%Pdep*deltpool
      casabal%FPsnetyear = casabal%FPsnetyear + casaflux%Psnet*deltpool
      casabal%FPupyear = casabal%FPupyear + casaflux%Plabuptake*deltpool
      casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach*deltpool  
      casabal%FPlossyear = casabal%FPlossyear + casaflux%Ploss*deltpool 
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
    canopy%frp  = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,xroot)+casaflux%crgplant(:)) &
                       /real(casaperiod,8)
    canopy%frs  = casaflux%Crsoil(:)/real(casaperiod,8)
    canopy%frpw = casaflux%crmplant(:,wood)/real(casaperiod,8)
    canopy%frpr = casaflux%crmplant(:,xroot)/real(casaperiod,8)
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnpp = -canopy%fpn - canopy%frp
    if ( ccycle<=1 ) then
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    else
      if ( progvcmax>0 ) then
        canopy%fnee = canopy%fpn + canopy%frs + canopy%frp &
                          + casaflux%clabloss/real(casaperiod,8)
      else
        canopy%fnee = (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss) &
                           /real(casaperiod,8)
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
! Special
if ( nspecial==51 ) then
  ! reset soil moisture to prescribed climatology  
  do k = 1,cbm_ms
    call cable_pack(wb_clim(:,k),wbclim_pack(:),tile)  
    ! 0 >= rad%longitude >= 360 and -90 >= rad%latitude >= 90
    where( rad%longitude>=wbclim_lonn .and. rad%longitude<=wbclim_lonx .and. &
           rad%latitude>=wbclim_latn .and. rad%latitude<=wbclim_latx )
      ssnow%wb(:,k) = real(wbclim_pack(:),8)
    end where
  end do
end if

!--------------------------------------------------------------
 
! Unpack tiles into grid point averages.
! Note that albsav and albnirsav are the VIS and NIR albedo output from CABLE to
! be used by the radiadiation scheme at the next time step.
do k = 1,cbm_ms
  where ( land(:) )
    tgg(:,k) = 0.
    wb(:,k) = 0.
    wbice(:,k) = 0.
  end where
end do
do k = 1,3
  where ( land(:) )
    tggsn(:,k) = 0.
    smass(:,k) = 0.
    ssdn(:,k) = 0.
  end where
end do
where ( land(:) )
  sigmf(:) = 0.  
  albvisdir(:) = 0.
  albvisdif(:) = 0.
  albnirdir(:) = 0.
  albnirdif(:) = 0.
  fg(:) = 0.
  eg(:) = 0.
  ga(:) = 0.
  rnet(:) = 0.
  epot(:) = 0.
  tss(:) = 0.
  zo(:) = 0.
  zoh(:) = 0.
  zoq(:) = 0.
  cduv(:) = 0.
  cdtq(:) = 0.
  ustar(:) = 0.
  wetfac(:) = 0.
  rsmin(:) = 0.
  ssdnn(:) = 0.
  snowd(:) = 0.
  snage(:) = 0.
end where
evspsbl_l(:) = 0.
sbl_l(:) = 0.
tmps(:) = 0. ! average isflag

! veg fraction
dumt(:) = 1._8-exp(-0.4_8*veg%vlai)
call cable_update(sigmf(:),dumt(:),tile)
! albedo
call cable_update(albvisdir(:),rad%reffbm(:,1),tile)
call cable_update(albnirdir(:),rad%reffbm(:,2),tile)
call cable_update(albvisdif(:),rad%reffdf(:,1),tile)
call cable_update(albnirdif(:),rad%reffdf(:,2),tile)
! fluxes
call cable_update(fg(:),canopy%fh,tile)
call cable_update(eg(:),canopy%fe,tile)
call cable_update(ga(:),canopy%ga,tile)
call cable_update(rnet(:),canopy%rnet,tile)
dumt(:) = rad%trad**4
call cable_update(tss(:),dumt(:),tile)
! drag and mixing
dumt(:) = 1._8/(log(real(zmin,8)/rough%z0m)**2)
call cable_update(zo(:),dumt(:),tile)
call cable_update(cduv(:),canopy%cduv,tile)
call cable_update(cdtq(:),canopy%cdtq,tile)
! soil
do k = 1,cbm_ms
  call cable_update(tgg(:,k),ssnow%tgg(:,k),tile)
  call cable_update(wb(:,k),ssnow%wb(:,k),tile)
  call cable_update(wbice(:,k),ssnow%wbice(:,k),tile)
end do
! hydrology
dumt(:) = ssnow%runoff(:)*dtr8 ! convert mm/s to mm
call cable_update(runoff(:),dumt(:),tile)
dumt(:) = ssnow%rnof1(:)*dtr8  ! convert mm/s to mm
call cable_update(runoff_surface(:),dumt(:),tile)
call cable_update(fwet(:),canopy%fwet,tile)
call cable_update(wetfac(:),ssnow%wetfac,tile)
call cable_update(cansto(:),canopy%cansto,tile)
! snow
dumt(:) = real(ssnow%isflag,8)
call cable_update(tmps(:),dumt(:),tile) ! used in radiation (for nsib==3)
do k = 1,3
  call cable_update(tggsn(:,k),ssnow%tggsn(:,k),tile)
  call cable_update(smass(:,k),ssnow%smass(:,k),tile)
  call cable_update(ssdn(:,k),ssnow%ssdn(:,k),tile)
end do
call cable_update(ssdnn(:),ssnow%ssdnn,tile)
call cable_update(snage(:),ssnow%snage,tile)
call cable_update(snowd(:),ssnow%snowd,tile)
call cable_update(snowmelt(:),ssnow%smelt,tile)
dumt(:) = (canopy%fev+canopy%fesp+canopy%fess/ssnow%cls)/real(hl,8)
call cable_update(evspsbl_l(:),dumt(:),tile)
dumt(:) = ssnow%evapsn/dtr8
call cable_update(sbl_l(:),dumt(:),tile)
! Replace potev with Penman_Monteith
if ( cable_potev==1 .and. cable_enablefao==1 ) then
  dumtr4(:) = real( met%tvair(:) )
  dumpr4(:) = real( met%pmb(:) )*100. ! convert from mb to Pa
  qsatfvar(:) = qsat(dumpr4(:),dumtr4(:))
  cc1(:) = air%dsatdk/(air%dsatdk+air%psyc)
  cc2(:) = air%psyc/(air%dsatdk+air%psyc)
  ssnowpotev(:) = cc1(:)*(canopy%fns - canopy%ga) +                  &
                      cc2(:)*(air%rho*air%rlam)*(real(qsatfvar(:),8) &
                    - met%qvair)/ssnow%rtsoil
  call cable_update(epot(:),ssnowpotev(:),tile)
else
  call cable_update(epot(:),ssnow%potev,tile)
end if
call cable_update(vlai(:),veg%vlai,tile)
call cable_update(rsmin(:),canopy%gswx_T,tile)

if ( cable_gw_model==1 ) then
  wtd(:) = 0.
  dumt(:) = ssnow%wtd/1000._8
  call cable_update(wtd(:),dumt(:),tile)
end if ! wt_transport==1

if ( ccycle/=0 ) then
  cplant(:,:) = 0.
  niplant(:,:) = 0.
  pplant(:,:) = 0.
  clitter(:,:) = 0.
  nilitter(:,:) = 0.
  plitter(:,:) = 0.
  csoil(:,:) = 0.
  nisoil(:,:) = 0.
  psoil(:,:) = 0.
  fnee(:) = 0.
  fpn(:) = 0.
  frd(:) = 0.
  frp(:) = 0.
  frpw(:) = 0.
  frpr(:) = 0.
  frs(:) = 0.
  cnpp(:) = 0.
  cnbp(:) = 0.
  if ( diaglevel_carbon > 0 ) then
    fevc(:) = 0.
    plant_turnover(:) = 0.
    plant_turnover_wood(:) = 0.
  end if
  do k = 1,mplant
    call cable_update(cplant(:,k),casapool%cplant(:,k),tile)  
    call cable_update(niplant(:,k),casapool%nplant(:,k),tile)  
    call cable_update(pplant(:,k),casapool%pplant(:,k),tile)  
  end do
  do k = 1,mlitter
    call cable_update(clitter(:,k),casapool%clitter(:,k),tile)  
    call cable_update(nilitter(:,k),casapool%nlitter(:,k),tile)  
    call cable_update(plitter(:,k),casapool%plitter(:,k),tile)  
  end do
  do k = 1,msoil
    call cable_update(csoil(:,k),casapool%csoil(:,k),tile)  
    call cable_update(nisoil(:,k),casapool%nsoil(:,k),tile)  
    call cable_update(psoil(:,k),casapool%psoil(:,k),tile)  
  end do
  !call cable_update(glai(:),casamet%glai,tile)
  ! carbon cycle
  call cable_update(fnee(:),canopy%fnee,tile)
  call cable_update(fpn(:),canopy%fpn,tile)
  call cable_update(frd(:),canopy%frday,tile)
  call cable_update(frp(:),canopy%frp,tile)
  call cable_update(frpw(:),canopy%frpw,tile)
  call cable_update(frpr(:),canopy%frpr,tile)
  call cable_update(frs(:),canopy%frs,tile)
  dumt(:) = casaflux%cnpp/real(casaperiod,8)
  call cable_update(cnpp(:),dumt(:),tile)
  dumt(:) = (casaflux%Crsoil-casaflux%cnpp-casapool%dClabiledt) &
     /real(casaperiod,8)
  call cable_update(cnbp(:),dumt(:),tile)
  if ( diaglevel_carbon > 0 ) then
    call cable_update(fevc(:),canopy%fevc,tile)  
    dumt(:) = sum(casaflux%Cplant_turnover,2)/real(casaperiod,8)
    call cable_update(plant_turnover(:),dumt(:),tile)
    dumt(:) = casaflux%Cplant_turnover(:,2)/real(casaperiod,8)
    call cable_update(plant_turnover_wood(:),dumt(:),tile)
  end if
end if


! MJT notes - ustar, cduv, fg and eg are passed to the boundary layer turbulence scheme
! zoh, zoq and zo are passed to the scrnout diagnostics routines
! rsmin is typically used by CTM

where ( land(:) )
  zo(:)      = zmin*exp(-1./sqrt(zo(:)))
  zoh(:)     = zo(:)/7.4
  zoq(:)     = zoh(:)
  ustar(:)   = sqrt(cduv(:))*vmod(:)  
  cduv(:)    = cduv(:)*vmod(:)           ! cduv is Cd*vmod in CCAM
  cdtq(:)    = cdtq(:)*vmod(:)
  tss(:)     = tss(:)**0.25
  rsmin(:)   = 1./rsmin(:)
  evspsbl(:) = evspsbl(:) + dt*evspsbl_l(:)
  sbl(:)     = sbl(:) + dt*sbl_l(:)
  ! update albedo and tss before calculating net radiation
  albvissav(:) = fbeamvis(:)*albvisdir(:) + (1.-fbeamvis(:))*albvisdif(:) ! for nrad=4
  albnirsav(:) = fbeamnir(:)*albnirdir(:) + (1.-fbeamnir(:))*albnirdif(:) ! for nrad=4 
  isflag(:) = nint(tmps(:)) ! tmps is average isflag
end where
qsttg_land(:) = qsat(ps(:),tss(:)) ! must wait for tss to be updated first
where ( land(:) )
  qsttg(:)  = qsttg_land(:)
end where

return
end subroutine sib4_work

! *************************************************************************************
subroutine cbmemiss(trsrc,mvegt,mode,tile,imax)
  
use parm_m
  
integer, intent(in) :: tile,imax
integer, intent(in) :: mvegt,mode
real, dimension(imax), intent(out) :: trsrc

call cbmemiss_work(trsrc,mvegt,mode,imax,tdata(tile)%tind,tdata(tile)%tmap,tdata(tile)%maxnb,tile)
  
return
end subroutine cbmemiss

subroutine cbmemiss_work(trsrc,mvegt,mode,imax,tind,tmap,maxnb,tile)
  
use parm_m
  
integer, intent(in) :: imax, tile
integer, intent(in) :: mvegt,mode
integer nb
real, dimension(imax), intent(out) :: trsrc
integer, dimension(maxtile,2), intent(in) :: tind
logical, dimension(imax,maxtile), intent(in) :: tmap
integer, intent(in) :: maxnb
real, dimension(imax) :: fpn,frd,frp,frs

if ( nsib/=6 .and. nsib/=7 ) then
  write(6,*) "ERROR: Attempted to read CABLE emissions with CABLE disabled"
  stop
end if

fpn=0.
frd=0.
frp=0.
frs=0.
do nb = 1,maxnb
  where ( veg(tile)%iveg(:)==mvegt )
    fpn(:) = fpn(:) + unpack(tdata(tile)%sv(:)*real(canopy(tile)%fpn(:)),  tmap(:,nb),0.)
    frd(:) = frd(:) + unpack(tdata(tile)%sv(:)*real(canopy(tile)%frday(:)),tmap(:,nb),0.)
    frp(:) = frp(:) + unpack(tdata(tile)%sv(:)*real(canopy(tile)%frp(:)),  tmap(:,nb),0.)
    frs(:) = frs(:) + unpack(tdata(tile)%sv(:)*real(canopy(tile)%frs(:)),  tmap(:,nb),0.)
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
! Calculate prognostic vcmax from CASA-CNP
subroutine vcmax_feedback(casabiome,casamet,casapool,veg,climate,ktau)

use newmpar_m, only : mxvt
use parm_m, only : nperday

type(casa_biome), intent(in) :: casabiome
type(casa_met), intent(in) :: casamet
type(casa_pool), intent(in) :: casapool
type(veg_parameter_type), intent(inout) :: veg
type(climate_type), intent(in) :: climate
integer, intent(in) :: ktau
integer np, ivt
real(kind=8) bjvref
real(kind=8), dimension(cbm_mp) :: ncleafx, npleafx, pleafx, nleafx
real(kind=8), dimension(mxvt), parameter :: xnslope = (/ 0.8,1.,2.,1.,1.,1.,0.5,1.,0.34,1.,1.,1.,1.,1.,1.,1.,1. /)

if ( progvcmax>0 .and. ccycle>=2 ) then

  ! initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx(:) = casabiome%ratioNPplantmin(veg%iveg(:),leaf)

  do np = 1,cbm_mp
    ivt = veg%iveg(np)
    if ( casamet%iveg2(np)/=icewater .and. casamet%glai(np)>casabiome%glaimin(ivt) .and. &
         casapool%cplant(np,leaf)>0._8 ) then

      if (ccycle>1 .and. casapool%cplant(np,leaf)>0._8) then
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
      do np = 1,cbm_mp
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
  
!    case(2)
!      !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax – 
!      !to leaf nitrogen, leaf phosphorus, and specific leaf area: 
!      !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
!      bjvref = 1.7_8
!      do np = 1,mp
!        ivt = veg%iveg(np)  
!        nleafx(np) = ncleafx(np)/casabiome%sla(ivt) ! leaf N in g N m-2 leaf
!        pleafx(np) = nleafx(np)/npleafx(np)         ! leaf P in g P m-2 leaf
!
!        if (ivt==7 .or. ivt==9  ) then
!          veg%vcmax(np) = 1.0e-5_8 ! special for C4 grass: set here to value from  parameter file
!          veg%ejmax(np) = 2._8*veg%vcmax(np)
!        elseif (ivt==1) then
!          ! account here for spring recovery
!          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))*climate%frec(np)  
!          veg%ejmax(np) =bjvref*veg%vcmax(np)
!        else
!          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))
!          veg%ejmax(np) =bjvref*veg%vcmax(np)
!        endif
!      end do
!
!      if ( cable_user%finite_gm ) then
!        ! vcmax and jmax modifications according to Sun et al. 2014 Table S3  
!        where( veg%iveg==1 )
!          veg%vcmax = veg%vcmax*2.2_8
!          veg%ejmax = veg%vcmax*1.1_8
!        elsewhere ( veg%iveg==2 )
!          veg%vcmax = veg%vcmax*1.9_8
!          veg%ejmax = veg%vcmax*1.2_8
!        elsewhere ( veg%iveg==3 )
!          veg%vcmax = veg%vcmax*1.4_8
!          veg%ejmax = veg%vcmax*1.5_8
!        elsewhere ( veg%iveg==4 )
!          veg%vcmax = veg%vcmax*1.45_8
!          veg%ejmax = veg%vcmax*1.3_8
!        elsewhere ( veg%iveg==5 )
!          veg%vcmax = veg%vcmax*1.7_8
!          veg%ejmax = veg%vcmax*1.2_8
!        elsewhere ( veg%iveg==6 .or. veg%iveg==8 .or. veg%iveg==9 )
!          veg%vcmax = veg%vcmax*1.6_8
!          veg%ejmax = veg%vcmax*1.2_8
!        end where
!      end if    
      
    case default
      write(6,*) "ERROR: Invalid progvcmax ",progvcmax
      stop
  end select
   
end if
  
return
end subroutine vcmax_feedback

! *************************************************************************************
! POP subroutines
subroutine casa_pop_firstcall(lalloc,casabiome,casaflux,casamet,casapool,phen,pop,soil,veg)

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
  
type(casa_balance), intent(in) :: casabal
type(casa_flux), intent(in) :: casaflux
type(pop_type), intent(inout) :: pop
type(veg_parameter_type), intent(in) :: veg
integer, dimension(POP%np) :: Iw
integer(kind=i4b), dimension(POP%np,2) :: disturbance_interval_tmp
real(kind=dp), dimension(POP%np,2) :: disturbance_intensity_tmp, stemNPP_tmp
real(kind=dp), dimension(POP%np) :: LAImax_tmp, Cleafmean_tmp, Crootmean_tmp, NPPtoGPP_tmp
real(kind=dp), dimension(cbm_mp) :: NPPtoGPP

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

subroutine newcbmwb

use soilsnow_m
use newmpar_m

integer k, tile, js, je, nb

do k = 1,cbm_ms
  do tile = 1,ntiles
    js = (tile-1)*imax + 1
    je = tile*imax
    wb(js:je,k) = 0.
    call cable_update(wb(js:je,k),ssnow(tile)%wb(:,k),tile)
  end do
end do

return
end subroutine newcbmwb

subroutine cable_casatile(dat_out,vname,k,n)

use newmpar_m

implicit none

integer, intent(in) :: k, n
integer tile, is, ie
real, dimension(ifull), intent(out) :: dat_out
character(len=*), intent(in) :: vname

dat_out = 0.

select case(vname)
  case("cplant")
    select case(k)
      case(1)
        dat_out(:) = real(cleaf(n))
      case(2)
        dat_out(:) = real(cwood(n))
      case(3)
        dat_out(:) = real(cfroot(n))
    end select  
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%cplant(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,             &
                               tdata(tile)%maxnb)
    end do  
  case("nplant")
    select case(k)
      case(1)
        dat_out(:) = real(nleaf(n))
      case(2)
        dat_out(:) = real(nwood(n))
      case(3)
        dat_out(:) = real(nfroot(n))
    end select
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%nplant(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,             &
                               tdata(tile)%maxnb)
    end do  
  case("pplant")
    select case(k)
      case(1)
        dat_out(:) = real(xpleaf(n))
      case(2)
        dat_out(:) = real(xpwood(n))
      case(3)
        dat_out(:) = real(xpfroot(n))  
      end select      
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%pplant(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,             &
                               tdata(tile)%maxnb)
    end do  
  case("clitter")
    select case(k)
      case(1)
        dat_out(:) = real(cmet(n))
      case(2)
        dat_out(:) = real(cstr(n))
      case(3)
        dat_out(:) = real(ccwd(n))
    end select
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%clitter(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,              &
                               tdata(tile)%maxnb)
    end do  
  case("nlitter")
    select case(k)
      case(1)
        dat_out(:) = real(nmet(n))
      case(2)
        dat_out(:) = real(nstr(n))
      case(3)
        dat_out(:) = real(ncwd(n))
    end select      
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%nlitter(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,              &
                               tdata(tile)%maxnb)
    end do  
  case("plitter")
    select case(k)
      case(1)
        dat_out(:) = real(xpmet(n))
      case(2)
        dat_out(:) = real(xpstr(n)) 
      case(3)
        dat_out(:) = real(xpcwd(n))
    end select      
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%plitter(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,              &
                               tdata(tile)%maxnb)
    end do  
  case("csoil")
    select case(k)
      case(1)
        dat_out(:) = real(cmic(n))
      case(2)
        dat_out(:) = real(cslow(n))
      case(3)
        dat_out(:) = real(cpass(n))  
    end select
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%csoil(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,            &
                               tdata(tile)%maxnb)
    end do  
  case("nsoil")
    select case(k)
      case(1)
        dat_out(:) = real(nmic(n))
      case(2)
        dat_out(:) = real(nslow(n))
      case(3)
        dat_out(:) = real(npass(n))
    end select      
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%nsoil(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,            &
                               tdata(tile)%maxnb)
    end do  
  case("psoil")
    select case(k)
      case(1)
        dat_out(:) = real(xpmic(n))
      case(2)
        dat_out(:) = real(xpslow(n))
      case(3)
        dat_out(:) = real(xppass(n))
    end select      
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      call cable_casatile_work(n,casapool(tile)%psoil(:,k),dat_out(is:ie), &
                               veg(tile)%iveg,tdata(tile)%tmap,            &
                               tdata(tile)%maxnb)
    end do  
  case default
    write(6,*) "ERROR: Unknown option for cable_tile ",trim(vname)
    stop -1
end select

end subroutine cable_casatile

subroutine cable_casatile_work(n,dat_in,dat_out,iveg,tmap,maxnb)

use newmpar_m

implicit none

integer, intent(in) :: n, maxnb
integer, dimension(:), intent(in) :: iveg
integer nb
real, dimension(:), intent(inout) :: dat_out
real(kind=8), dimension(:), intent(in) :: dat_in
integer, dimension(size(dat_out)) :: iveg_local
real(kind=8), dimension(size(dat_out)) :: dat_local
logical, dimension(:,:), intent(in) :: tmap

iveg_local(:) = 0.
dat_local(:) = 0.
do nb = 1,maxnb
  iveg_local(:) = unpack(iveg,tmap(:,nb),0)
  dat_local(:) = unpack(dat_in,tmap(:,nb),0._8)
  where ( iveg_local(:)==n .and. dat_local(:)>1.e-8_8 )
     dat_out(:) = real(dat_local(:))
  end where
end do  

return
end subroutine cable_casatile_work

end module cable_ccam

