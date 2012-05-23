module cable_ccam

! CABLE interface originally developed by the CABLE group
! Subsequently modified by MJT for 5 tile mosaic and SEAESF radiation scheme
  
! - Currently all tiles have the same soil texture, but independent soil temperatures,
!   moisture, etc.
! - LAI can be interpolated between timesteps using a PWCB fit to the LAI integral
!   or LAI can be taken as constant for the month
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.  However,
!   since the input files can be modified at runtime (not all months are loaded
!   at once), then we can support evolving/dynamic vegetation, etc.

! ivegt   type
! 1       Evergreen Needleleaf Forest
! 2       Evergreen Broadleaf Forest
! 3       Deciduous Needleaf Forest
! 4       Deciduous Broadleaf Forest
! 5       Mixed Forest
! 6       Closed Shrublands
! 7       Open Shrublands
! 8       Woody Savannas
! 9       Savannas
! 10      Grasslands
! 11      Permanent Wetlands
! 12      Croplands
! 13      Urban and Built-up
! 14      Cropland/Natural Vegetation Mosaic
! 15      Snow and Ice
! 16      Barren or Sparsely Vegetated
! 17      Water Bodies
  
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

use air_module
use albedo_module
use canopy_module
use carbon_module
use casa_cnp
use casadimension
use casaparm, xroot => froot
use casavariable
use cable_common_module
use define_dimensions, cbm_ms => ms
use define_types
use physical_constants
use radiation_module
use roughness_module
use soil_snow_module

implicit none

private
public sib4,loadcbmparm,loadtile,savetile,cableinflow,cbmemiss

integer, parameter :: hruffmethod    = 1 ! Method for max hruff
integer, parameter :: CO2forcingtype = 2 ! CO2 input source (1 constant, 2 use radiation CO2 forcing, 3 interactive tracer)
integer, parameter :: proglai        = 0 ! 0 prescribed LAI, 1 prognostic LAI 
integer, dimension(:), allocatable, save :: cmap
integer, dimension(5,2), save :: pind  
real, dimension(:), allocatable, save :: sv,vl1,vl2,vl3
type (air_type), save :: air
type (bgc_pool_type), save :: bgc
type (met_type), save :: met
type (balances_type), save :: bal
type (radiation_type), save :: rad
type (roughness_type), save :: rough
type (soil_parameter_type), save :: soil
type (soil_snow_type), save :: ssoil
type (sum_flux_type), save :: sum_flux
type (veg_parameter_type), save :: veg
type (canopy_type), save :: canopy
type (casa_balance), save :: casabal
type (casa_biome), save :: casabiome
type (casa_flux), save :: casaflux
type (casa_met), save :: casamet
type (casa_pool), save :: casapool
type (phen_variable), save :: phen

contains
! ****************************************************************************

! CABLE-CCAM interface
subroutine sib4

use arrays_m
use carbpools_m
use extraout_m
use infile
use latlong_m
use morepbl_m
use nharrs_m
use nsibd_m
use pbl_m
use permsurf_m
use prec_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use vegpar_m
use work2_m, only : qsttg,zo,theta,vmod
use work3_m, only : ga
use zenith_m
  
implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'
include 'establ.h'
include 'parm.h'

! for calculation of zenith angle
real fjd,r1,dlt,slag,dhr,alp,x,esatf
real, dimension(ifull) :: coszro2,taudar2,tmps,hruff_grmx,atmco2
real, dimension(ifull) :: wetfac
real, dimension(mp) :: swdwn,deltat
real(r_2), dimension(mp) :: xKNlimiting,xkleafcold,xkleafdry
real(r_2), dimension(mp) :: xkleaf,xnplimit,xNPuptake,xklitter
real(r_2), dimension(mp) :: xksoil
integer jyear,jmonth,jday,jhour,jmin
integer k,mins,nb,iq,j
integer(i_d) :: idoy

! abort calculation if no land points on this processor  
if (mp.le.0) return

! set meteorological forcing
dhr = dt/3600.
call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
fjd = float(mod(mins,525600))/1440.
call solargh(fjd,bpyear,r1,dlt,alp,slag)
call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro2,taudar2)
met%doy=fjd
met%tk=theta(cmap)
met%tvair=met%tk
met%tvrad=met%tk
met%ua=vmod(cmap)
met%ua=max(met%ua,umin)
call setco2for(atmco2)
met%ca=1.e-6*atmco2(cmap)
met%coszen=max(1.e-8,coszro2(cmap)) ! use instantaneous value
met%qv=qg(cmap,1)         ! specific humidity in kg/kg
met%pmb=0.01*ps(cmap)     ! pressure in mb at ref height
met%precip=condx(cmap)    ! in mm not mm/sec
met%precip_sn=conds(cmap) ! in mm not mm/sec
met%hod=real(mtimer+jhour*60+jmin)/60.+rlongg(cmap)*12./pi
met%hod=mod(met%hod,24.)
rough%za_tq=(bet(1)*t(cmap,1)+phi_nh(cmap,1))/grav ! reference height
rough%za_uv=rough%za_tq
! swrsave indicates the fraction of net VIS radiation (compared to NIR)
! fbeamvis indicates the beam fraction of downwelling direct radiation (compared to diffuse) for VIS
! fbeamnir indicates the beam fraction of downwelling direct radiation (compared to diffuse) for NIR
swdwn=sgsave(cmap)/(1.-swrsave(cmap)*albvisnir(cmap,1)-(1.-swrsave(cmap))*albvisnir(cmap,2)) ! short wave down (positive) W/m^2
met%fsd(:,1)=swrsave(cmap)*swdwn
met%fsd(:,2)=(1.-swrsave(cmap))*swdwn
rad%fbeam(:,1)=fbeamvis(cmap)
rad%fbeam(:,2)=fbeamnir(cmap)
rad%fbeam(:,3)=0. ! dummy for now
met%fld=-rgsave(cmap) ! long wave down (positive) W/m^2
! Interpolate lai.  Also need sigmf for LDR prognostic aerosols.
call setlai(sigmf,jyear,jmonth,jday,jhour,jmin)

rough%hruff=max(0.01,veg%hc-1.2*ssoil%snowd/max(ssoil%ssdnn,100.))
select case(hruffmethod)
  case(0) ! hruff is mixed in a tile (find max hruff for tile)
    hruff_grmx=0.01
    do nb=1,5
      if (pind(nb,1).le.mp) then
        hruff_grmx(cmap(pind(nb,1):pind(nb,2)))=max( &
          hruff_grmx(cmap(pind(nb,1):pind(nb,2))),rough%hruff(pind(nb,1):pind(nb,2)))
      end if
    end do
    rough%hruff_grmx=hruff_grmx(cmap)
  case(1) ! hruff is seperate in a tile (no max hruff for tile)
    rough%hruff_grmx=rough%hruff
  case DEFAULT
    write(6,*) "ERROR: Unsupported hruffmethod ",hruffmethod
    stop
end select
  
!--------------------------------------------------------------
! CABLE
ktau_gl=900
kend_gl=999
ssoil%owetfac = ssoil%wetfac
canopy%oldcansto=canopy%cansto
call ruff_resist(veg,rough,ssoil,soil,met,canopy)
call define_air(met,air)
call init_radiation(met,rad,veg,canopy)
call surface_albedo(ssoil,veg,met,rad,soil,canopy)
call define_canopy(bal,rad,rough,air,met,dt,ssoil,soil,veg,canopy)
ssoil%otss = ssoil%tss
call soil_snow(dt,soil,ssoil,canopy,met,bal,veg)
! adjust for new soil temperature
deltat = ssoil%tss - ssoil%otss
canopy%fhs     = canopy%fhs + deltat*ssoil%dfh_dtg
canopy%fes     = canopy%fes + deltat*ssoil%cls*ssoil%dfe_ddq*ssoil%ddq_dtg
canopy%fhs_cor = canopy%fhs_cor + deltat*ssoil%dfh_dtg
canopy%fes_cor = canopy%fes_cor + deltat*ssoil%cls*ssoil%dfe_ddq*ssoil%ddq_dtg
canopy%fh      = canopy%fhv + canopy%fhs
canopy%fev     = canopy%fevc + canopy%fevw
canopy%fe      = canopy%fev + canopy%fes
canopy%rnet    = canopy%fns + canopy%fnv
rad%trad       = ( (1.-rad%transd)*canopy%tv**4 + rad%transd*ssoil%tss**4 )**0.25

!--------------------------------------------------------------
! CASA CNP
select case (icycle)
  case(0) ! off
    call plantcarb(veg,bgc,met,canopy)
    call soilcarb(soil,ssoil,veg,bgc,met,canopy)
    call carbon_pl(dt,soil,ssoil,veg,canopy,bgc)
    canopy%fnpp = -canopy%fpn - canopy%frp
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
  case(3) ! C+N+P
    ! update casamet
    casamet%tairk = casamet%tairk + met%tk
    casamet%tsoil = casamet%tsoil + ssoil%tgg
    casamet%moist = casamet%moist + ssoil%wb
    casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dt
    casaflux%crmplant(:,leaf) =casaflux%crmplant(:,leaf) + canopy%frday*dt
    ! run CASA CNP once per day
    if (mod(ktau,nperday).eq.0) then
      casamet%tairk=casamet%tairk/real(nperday)
      casamet%tsoil=casamet%tsoil/real(nperday)
      casamet%moist=casamet%moist/real(nperday)
      xKNlimiting = 1.
      idoy=nint(fjd)
      call phenology(idoy,veg,phen)
      call avgsoil(veg,soil,casamet)
      call casa_rplant(veg,casabiome,casapool,casaflux,casamet)
      call casa_allocation(veg,soil,casabiome,casaflux,casamet,phen)
      call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                           casamet,phen)
      call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                           casaflux,casamet)
      call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)
      call casa_xratesoil(xklitter,xksoil,veg,soil,casamet)
      call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
      call casa_xkN(xkNlimiting,casapool,casaflux,casamet,veg)
      do j=1,mlitter
        casaflux%klitter(:,j) = casaflux%klitter(:,j)*xkNlimiting
      end do
      call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      call casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      call casa_delplant(veg,casabiome,casapool,casaflux,casamet)
      call casa_delsoil(veg,casapool,casaflux,casamet)
      call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)
      call casa_cnpbal(casapool,casaflux,casabal)
      casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
      casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
      casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp   * deltpool
      casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
      casabal%FCneeyear = casabal%FCneeyear + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
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
    canopy%frp  = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,xroot)+casaflux%crgplant(:))/86400.
    canopy%frs  = casaflux%Crsoil(:)/86400.
    canopy%frpw = casaflux%crmplant(:,wood)/86400.
    canopy%frpr = casaflux%crmplant(:,xroot)/86400.
    ! Set net ecosystem exchange after adjustments to frs:
    !if (l_vcmaxFeedbk) then
    !  canopy%fnee = canopy%fpn+canopy%frs+canopy%frp+casaflux%clabloss/86400.
    !else
      canopy%fnee = (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/86400.
    !end if
  case default
    write(6,*) "ERROR: Unknown icycle ",icycle
    stop
end select  
  
sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dt
sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dt
sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dt
sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dt
sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dt
sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dt
sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dt
sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dt
sum_flux%sumrs  = sum_flux%sumrs  + canopy%frs*dt
!--------------------------------------------------------------
      
! Unpack tiles into grid point averages.
! Note that albsav and albnirsav are the VIS and NIR albedo output from CABLE to
! be used by the radiadiation scheme at the next time step.  albvisnir(:,1) and
! albvisnir(:,2) are the VIS and NIR albedo used by the radiation scheme for the
! current time step.
tgg(iperm(1:ipland),:)=0.
wb(iperm(1:ipland),:)=0.
wbice(iperm(1:ipland),:)=0.
cplant(iperm(1:ipland),:)=0.
csoil(iperm(1:ipland),:)=0.
albsav(iperm(1:ipland))=0.
albnirsav(iperm(1:ipland))=0.
albvisdir(iperm(1:ipland))=0.
albvisdif(iperm(1:ipland))=0.
albnirdir(iperm(1:ipland))=0.
albnirdif(iperm(1:ipland))=0.
! From 11/8/98 runoff() is accumulated & zeroed with precip
rnet(iperm(1:ipland))=0.
fg(iperm(1:ipland))=0.
eg(iperm(1:ipland))=0.
ga(iperm(1:ipland))=0.
epot(iperm(1:ipland))=0.
tss(iperm(1:ipland))=0.
cansto=0.
fwet=0.
fnee=0.
fpn=0.
frd=0.
frp=0.
frpw=0.
frpr=0.
frs=0.
zo(iperm(1:ipland))=0.
cduv(iperm(1:ipland))=0.
cdtq(iperm(1:ipland))=0.
ustar(iperm(1:ipland))=0.
wetfac(iperm(1:ipland))=0.
tmps=0. ! average isflag
vlai=0.
      
! screen and 10m diagnostics - rhscrn calculated in sflux.f
tscrn(iperm(1:ipland))=0.
uscrn(iperm(1:ipland))=0.
qgscrn(iperm(1:ipland))=0.
u10(iperm(1:ipland))=0.

if (icycle==0) then
  cplant=0.
  csoil=0.
else
  cplant=0.
  niplant=0.
  pplant=0.
  clitter=0.
  nilitter=0.
  plitter=0.
  csoil=0.
  nisoil=0.
  psoil=0.
  glai=0.
end if
 
do nb=1,5
  if (pind(nb,1).le.mp) then
    do k=1,ms
      tgg(cmap(pind(nb,1):pind(nb,2)),k)=tgg(cmap(pind(nb,1):pind(nb,2)),k) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),k)
      wb(cmap(pind(nb,1):pind(nb,2)),k)=wb(cmap(pind(nb,1):pind(nb,2)),k) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%wb(pind(nb,1):pind(nb,2),k)
      wbice(cmap(pind(nb,1):pind(nb,2)),k)=wbice(cmap(pind(nb,1):pind(nb,2)),k) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%wbice(pind(nb,1):pind(nb,2),k)
    end do
    if (icycle==0) then
      do k=1,ncp
        cplant(cmap(pind(nb,1):pind(nb,2)),k)=cplant(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*bgc%cplant(pind(nb,1):pind(nb,2),k)
      enddo
      do k=1,ncs
        csoil(cmap(pind(nb,1):pind(nb,2)),k)=csoil(cmap(pind(nb,1):pind(nb,2)),k) &
                                      +sv(pind(nb,1):pind(nb,2))*bgc%csoil(pind(nb,1):pind(nb,2),k)
      enddo
    else
      do k=1,mplant
        cplant(cmap(pind(nb,1):pind(nb,2)),k)=cplant(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%cplant(pind(nb,1):pind(nb,2),k)
        niplant(cmap(pind(nb,1):pind(nb,2)),k)=niplant(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%nplant(pind(nb,1):pind(nb,2),k)
        pplant(cmap(pind(nb,1):pind(nb,2)),k)=pplant(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%pplant(pind(nb,1):pind(nb,2),k)
      end do
      do k=1,mlitter
        clitter(cmap(pind(nb,1):pind(nb,2)),k)=clitter(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%clitter(pind(nb,1):pind(nb,2),k)
        nilitter(cmap(pind(nb,1):pind(nb,2)),k)=nilitter(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%nlitter(pind(nb,1):pind(nb,2),k)
        plitter(cmap(pind(nb,1):pind(nb,2)),k)=plitter(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%plitter(pind(nb,1):pind(nb,2),k)
      end do
      do k=1,msoil
        csoil(cmap(pind(nb,1):pind(nb,2)),k)=csoil(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%csoil(pind(nb,1):pind(nb,2),k)
        nisoil(cmap(pind(nb,1):pind(nb,2)),k)=nisoil(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%nsoil(pind(nb,1):pind(nb,2),k)
        psoil(cmap(pind(nb,1):pind(nb,2)),k)=psoil(cmap(pind(nb,1):pind(nb,2)),k) &
                                        +sv(pind(nb,1):pind(nb,2))*casapool%psoil(pind(nb,1):pind(nb,2),k)
      end do
      glai(cmap(pind(nb,1):pind(nb,2)))=glai(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*casamet%glai(pind(nb,1):pind(nb,2))
    end if
    albsav(cmap(pind(nb,1):pind(nb,2)))=albsav(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%albedo(pind(nb,1):pind(nb,2),1)
    albnirsav(cmap(pind(nb,1):pind(nb,2)))=albnirsav(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%albedo(pind(nb,1):pind(nb,2),2)
    albvisdir(cmap(pind(nb,1):pind(nb,2)))=albvisdir(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%reffbm(pind(nb,1):pind(nb,2),1)
    albnirdir(cmap(pind(nb,1):pind(nb,2)))=albnirdir(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%reffbm(pind(nb,1):pind(nb,2),2)
    albvisdif(cmap(pind(nb,1):pind(nb,2)))=albvisdif(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%reffdf(pind(nb,1):pind(nb,2),1)
    albnirdif(cmap(pind(nb,1):pind(nb,2)))=albnirdif(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%reffdf(pind(nb,1):pind(nb,2),2)
    runoff(cmap(pind(nb,1):pind(nb,2)))=runoff(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%runoff(pind(nb,1):pind(nb,2))*dt ! convert mm/s to mm
    rnet(cmap(pind(nb,1):pind(nb,2)))=rnet(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%rnet(pind(nb,1):pind(nb,2))
    fg(cmap(pind(nb,1):pind(nb,2)))=fg(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fh(pind(nb,1):pind(nb,2))
    eg(cmap(pind(nb,1):pind(nb,2)))=eg(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fe(pind(nb,1):pind(nb,2))
    ga(cmap(pind(nb,1):pind(nb,2)))=ga(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%ga(pind(nb,1):pind(nb,2))
    epot(cmap(pind(nb,1):pind(nb,2)))=epot(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%potev(pind(nb,1):pind(nb,2))
    tss(cmap(pind(nb,1):pind(nb,2)))=tss(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*rad%trad(pind(nb,1):pind(nb,2))**4 ! ave longwave radiation
    cansto(cmap(pind(nb,1):pind(nb,2)))=cansto(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%cansto(pind(nb,1):pind(nb,2))
    fwet(cmap(pind(nb,1):pind(nb,2)))=fwet(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fwet(pind(nb,1):pind(nb,2))
    fnee(cmap(pind(nb,1):pind(nb,2)))=fnee(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fnee(pind(nb,1):pind(nb,2))
    fpn(cmap(pind(nb,1):pind(nb,2)))=fpn(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fpn(pind(nb,1):pind(nb,2))
    frd(cmap(pind(nb,1):pind(nb,2)))=frd(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frday(pind(nb,1):pind(nb,2))
    frp(cmap(pind(nb,1):pind(nb,2)))=frp(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frp(pind(nb,1):pind(nb,2))
    frpw(cmap(pind(nb,1):pind(nb,2)))=frpw(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frpw(pind(nb,1):pind(nb,2))
    frs(cmap(pind(nb,1):pind(nb,2)))=frs(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frs(pind(nb,1):pind(nb,2))
    zo(cmap(pind(nb,1):pind(nb,2)))=zo(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))/log(zmin/max(rough%z0m(pind(nb,1):pind(nb,2)),zobgin))**2
    cduv(cmap(pind(nb,1):pind(nb,2)))=cduv(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%cduv(pind(nb,1):pind(nb,2))
    cdtq(cmap(pind(nb,1):pind(nb,2)))=cdtq(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%cdtq(pind(nb,1):pind(nb,2))
    ustar(cmap(pind(nb,1):pind(nb,2)))=ustar(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%us(pind(nb,1):pind(nb,2))
    wetfac(cmap(pind(nb,1):pind(nb,2)))=wetfac(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*ssoil%wetfac(pind(nb,1):pind(nb,2))
    tmps(cmap(pind(nb,1):pind(nb,2)))=tmps(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*real(ssoil%isflag(pind(nb,1):pind(nb,2)))
    vlai(cmap(pind(nb,1):pind(nb,2)))=vlai(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*veg%vlai(pind(nb,1):pind(nb,2))

    tscrn(cmap(pind(nb,1):pind(nb,2)))=tscrn(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%tscrn(pind(nb,1):pind(nb,2))
    uscrn(cmap(pind(nb,1):pind(nb,2)))=uscrn(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%uscrn(pind(nb,1):pind(nb,2))
    qgscrn(cmap(pind(nb,1):pind(nb,2)))=qgscrn(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%qscrn(pind(nb,1):pind(nb,2))
    !u10(cmap(pind(nb,1):pind(nb,2)))=u10(cmap(pind(nb,1):pind(nb,2))) &
    !                                  +sv(pind(nb,1):pind(nb,2))*canopy%ua_10m(pind(nb,1):pind(nb,2))
  end if
end do
where (land)
  zo=max(zmin*exp(-sqrt(1./zo)),zobgin)
  cduv=cduv*vmod     ! cduv is Cd*vmod in CCAM
  tscrn=tscrn+273.16 ! convert from degC to degK
  tss=tss**0.25
end where
do iq=1,ifull
  if (land(iq)) then
    esatf=establ(tss(iq))
    qsttg(iq)=.622*esatf/(ps(iq)-esatf)
  end if
end do
      
! The following lines unpack snow.  This is more complicated as we need to decide
! how to unpack tiles with 1 layer or 3 layers of snow in the same grid point.
! Here we estimate whether the majority of snow points is 1 layer or 3 layers and then
! convert each snow tile to that number of layers.  Note these calculations are purely
! diagnoistic.  They are not fed back into the CCAM simulation.
tggsn(iperm(1:ipland),:)=0.
smass(iperm(1:ipland),:)=0.
ssdn(iperm(1:ipland),:)=0.
ssdnn(iperm(1:ipland))=0.
snowd(iperm(1:ipland))=0.
snage(iperm(1:ipland))=0.
where (land.and.tmps.ge.0.5) ! tmps is average isflag
  isflag=1
elsewhere
  isflag=0
endwhere
do nb=1,5 ! update snow (diagnostic only)
  if (pind(nb,1).le.mp) then      
    do k=1,3
      where (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.1)          ! pack 1-layer into 3-layer
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)       ! pack 1-layer into 3-layer
        smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1) ! pack 1-layer into 3-layer
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)      ! pack 1-layer into 3-layer
      elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.2)      ! pack 1-layer into 3-layer
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)       ! pack 1-layer into 3-layer
        smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*(ssoil%snowd(pind(nb,1):pind(nb,2)) &    ! pack 1-layer into 3-layer
                                         -0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1))*0.4                      ! pack 1-layer into 3-layer
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)      ! pack 1-layer into 3-layer
      elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.3)      ! pack 1-layer into 3-layer
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)       ! pack 1-layer into 3-layer
        smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*(ssoil%snowd(pind(nb,1):pind(nb,2)) &    ! pack 1-layer into 3-layer
                                         -0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1))*0.6                      ! pack 1-layer into 3-layer
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! pack 1-layer into 3-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)      ! pack 1-layer into 3-layer
      elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).gt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.1)      ! pack 3-layer into 1-layer
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 3-layer into 1-layer
                                         +sv(pind(nb,1):pind(nb,2))*273.16                                   ! pack 3-layer into 1-layer
        smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 3-layer into 1-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%snowd(pind(nb,1):pind(nb,2))       ! pack 3-layer into 1-layer
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! pack 3-layer into 1-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdnn(pind(nb,1):pind(nb,2))       ! pack 3-layer into 1-layer
      elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).gt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.ge.2)      ! pack 3-layer into 1-layer
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! pack 3-layer into 1-layer
                                         +sv(pind(nb,1):pind(nb,2))*273.16                                   ! pack 3-layer into 1-layer
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! pack 3-layer into 1-layer
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),k)      ! pack 3-layer into 1-layer
      elsewhere                                                                                              ! no change in layers
        tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &                          ! no change in layers
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%tggsn(pind(nb,1):pind(nb,2),k)     ! no change in layers
        smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &                          ! no change in layers
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%smass(pind(nb,1):pind(nb,2),k)     ! no change in layers
        ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &                            ! no change in layers
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),k)      ! no change in layers
      end where
    end do
    ssdnn(cmap(pind(nb,1):pind(nb,2)))=ssdnn(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdnn(pind(nb,1):pind(nb,2))
    snage(cmap(pind(nb,1):pind(nb,2)))=snage(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%snage(pind(nb,1):pind(nb,2))
    snowd(cmap(pind(nb,1):pind(nb,2)))=snowd(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%snowd(pind(nb,1):pind(nb,2))
  end if
end do

return
end subroutine sib4

! *************************************************************************************
subroutine setco2for(atmco2)
! set co2 forcing for cable
! constant: atmospheric co2 = 360 ppm 
! host: atmospheric co2 follows that from CCAM radiation scheme
! interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)

use radisw_m, only : rrco2,ssolar,rrvco2
use tracermodule, only : tractype,tracname
use tracers_m

implicit none

include 'newmpar.h'

integer ico2,igas
real, dimension(ifull), intent(out) :: atmco2

! replace these options when 3d CO2 is fed to the radiation scheme
select case (CO2forcingtype)
  case (1)    ! constant
    atmco2 = 360.
  case (2)    ! from radiative CO2 forcings
    atmco2 = 1.E6*rrvco2
  case (3)    ! use interactive tracers
    ico2=0
    do igas=1,ngas
      if (trim(tractype(igas)).eq.'online') then
        if (trim(tracname(igas)).eq.'cbmnep') then
          ico2=igas
          exit
        end if
      end if
    end do
    if (ico2.le.0) then
      write(6,*) "ERROR: Cannot locate co2 in tracers"
      stop
    end if
    atmco2 = tr(1:ifull,1,ico2)
end select

return
end subroutine setco2for

! *************************************************************************************
subroutine cbmemiss(trsrc,mvegt,mode)
  
implicit none
  
include 'newmpar.h'
include 'parm.h'
  
integer, intent(in) :: mvegt,mode
integer nb
real, dimension(ifull), intent(out) :: trsrc
real, dimension(ifull) :: fpn,frd,frp,frs
  
if (nsib.ne.6.and.nsib.ne.7) then
  write(6,*) "ERROR: Attempted to read CABLE emissions with CABLE disabled"
  stop
end if
  
fpn=0.
frd=0.
frp=0.
frs=0.
  
do nb=1,5
  if (pind(nb,1).le.mp) then
    where (veg%iveg(pind(nb,1):pind(nb,2)).eq.mvegt)
      fpn(cmap(pind(nb,1):pind(nb,2)))=fpn(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%fpn(pind(nb,1):pind(nb,2))
      frd(cmap(pind(nb,1):pind(nb,2)))=frd(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frday(pind(nb,1):pind(nb,2))
      frp(cmap(pind(nb,1):pind(nb,2)))=frp(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frp(pind(nb,1):pind(nb,2))
      frs(cmap(pind(nb,1):pind(nb,2)))=frs(cmap(pind(nb,1):pind(nb,2))) &
                                      +sv(pind(nb,1):pind(nb,2))*canopy%frs(pind(nb,1):pind(nb,2))
    end where
  end if
end do
  
select case(mode)
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
subroutine setlai(sigmf,jyear,jmonth,jday,jhour,jmin)
  
implicit none

include 'newmpar.h'
include 'dates.h'
  
integer, intent(in) :: jyear,jmonth,jday,jhour,jmin
integer monthstart,nb,leap
integer, dimension(12) :: imonth
real, parameter :: vextkn = 0.4
real, dimension(ifull), intent(out) :: sigmf
real x
common/leap_yr/leap  ! 1 to allow leap years

select case(proglai)
  case(0)
    imonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    if (leap.eq.1) then
      if (mod(jyear,4)  .eq.0) imonth(2)=29
      if (mod(jyear,100).eq.0) imonth(2)=28
      if (mod(jyear,400).eq.0) imonth(2)=29
    end if

    monthstart=1440*(jday-1) + 60*jhour + jmin ! mins from start month
    x=min(max(real(mtimer+monthstart)/real(1440.*imonth(jmonth)),0.),1.)
    veg%vlai=vl1+vl2*x+vl3*x*x ! LAI as a function of time
    veg%vlai=max(veg%vlai,0.1)
    where (veg%iveg.ge.15.and.veg%iveg.le.17)
      veg%vlai=0.
    end where
  case(1)
    write(6,*) "ERROR: CASA CNP LAI not operational"
    stop
    !veg%vlai=max(casamet%glai,0.1)
  case default
    write(6,*) "ERROR: Unknown proglai option ",proglai
    stop
end select

sigmf=0.
do nb=1,5
  if (pind(nb,1).le.mp) then
    sigmf(cmap(pind(nb,1):pind(nb,2)))=sigmf(cmap(pind(nb,1):pind(nb,2))) &
      +sv(pind(nb,1):pind(nb,2))*(1.-exp(-vextkn*veg%vlai(pind(nb,1):pind(nb,2))))
  end if
end do
  
return
end subroutine setlai
  
  

! *************************************************************************************
subroutine loadcbmparm(fveg,fvegprev,fvegnext,fphen,casafile)

use carbpools_m
use cc_mpi
use infile
use latlong_m
use nsibd_m
use pbl_m
use sigs_m
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none
  
include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer(i_d), dimension(ifull,5) :: ivs
integer(i_d) iq,n,k,ipos,isoil
integer jyear,jmonth,jday,jhour,jmin,mins
real(r_1) :: totdepth
real(r_1), dimension(mxvt,ms) :: froot2
real(r_1), dimension(ifull,5) :: svs,vlin,vlinprev,vlinnext
real(r_1), dimension(ncp) :: ratecp
real(r_1), dimension(ncs) :: ratecs
real(r_1), dimension(mxvt,ncp) :: tcplant
real(r_1), dimension(mxvt,ncs) :: tcsoil
real(r_1), dimension(mxvt)   :: canst1,dleaf,ejmax,frac4,hc,rp20
real(r_1), dimension(mxvt)   :: rpcoef,shelrb,vcmax,xfang
real(r_1), dimension(mxvt)   :: tminvj,tmaxvj,vbeta
real(r_1), dimension(mxvt)   :: extkn,rootbeta,vegcf
real(r_1), dimension(mxvt,2) :: taul,refl  
real(r_1), dimension(mxvt)   :: leafage,woodage,frootage,metage
real(r_1), dimension(mxvt)   :: strage,cwdage,micage,slowage,passage
real(r_1), dimension(mxvt)   :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real(r_1), dimension(mxvt)   :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real(r_1), dimension(mxvt)   :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real(r_1), dimension(mxvt)   :: xnsoilmin,xplab,xpsorb,xpocc
real(r_1), dimension(mxvt)   :: cleaf,cwood,cfroot,cmet,cstr,ccwd,cmic,cslow,cpass,nleaf
real(r_1), dimension(mxvt)   :: nwood,nfroot,nmet,nstr,ncwd,nmic,nslow,npass,xpleaf,xpwood
real(r_1), dimension(mxvt)   :: xpfroot,xpmet,xpstr,xpcwd,xpmic,xpslow,xppass,clabileage
real(r_1), dimension(mxvt,mplant) :: ratiocnplant
real(r_1), dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real(r_1), dimension(12) :: xkmlabp,xpsorbmax,xfPleach
real(r_1), dimension(12,msoil) :: rationpsoil
real(r_1), dimension(ifull) :: dumr
real, dimension(ifull) :: albsoil,c4frac
real, dimension(ifull,2) :: albsoilsn
real fjd  
character(len=*), intent(in) :: fveg,fvegprev,fvegnext,fphen,casafile

if (myid==0) write(6,*) "Initialising CABLE"

if (cbm_ms.ne.ms) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  stop
end if

! redefine rhos
rhos=(/ 1600.,1595.,1381.,1373.,1476.,1521.,1373.,1537.,1455.,2600.,2600.,2600.,2600. /)

hc=(/ 17.,35.,15.5,20.,19.25,0.6,0.6,7.0426,8.,0.567,0.5,0.55,6.017,0.55,0.2,0.2,0.2 /)
xfang=(/ 0.01,0.1,0.01,0.25,0.125,0.,0.,-0.14,-0.01,-0.3,0.,0.,-0.17,0.,0.01,0.01,0. /)
dleaf=(/ 0.055,0.1,0.04,0.15,0.1,0.1,0.1,0.232,0.129,0.3,0.3,0.3,0.242,0.3,0.03,0.03,0.03 /)
canst1=0.1
shelrb=2.
extkn=0.001 ! new definition for nitrogen (CABLE v1.9b)
refl(:,1)=(/ 0.062,0.076,0.062,0.092,0.069,0.100,0.100,0.091,0.075,0.107,0.107,0.101,0.097,0.101,0.159,0.159,0.159 /)
refl(:,2)=(/ 0.302,0.350,0.302,0.380,0.336,0.400,0.400,0.414,0.347,0.469,0.469,0.399,0.396,0.399,0.305,0.305,0.305 /)
taul(:,1)=(/ 0.050,0.050,0.050,0.050,0.050,0.054,0.054,0.062,0.053,0.070,0.070,0.067,0.062,0.067,0.026,0.026,0.026 /)
taul(:,2)=(/ 0.100,0.250,0.100,0.250,0.158,0.240,0.240,0.221,0.166,0.250,0.250,0.225,0.232,0.250,0.126,0.126,0.126 /)
vegcf=(/ 0.91, 1.95, 0.73, 1.50, 1.55, 0.60, 2.05, 2.80, 2.75, 2.75, 0.00, 2.80, 0.00, 2.80, 0.00, 0.40, 0.40 /)
vcmax=(/ 65.2E-6,65.E-6,70.E-6,85.0E-6,80.E-6,20.E-6,20.E-6,10.E-6,20.E-6,10.E-6,50.E-6,80.E-6,1.E-6,80.E-6, &
         17.E-6,17.E-6,17.E-6 /)
ejmax=2.*vcmax
rp20=(/ 3.3039,1.1342,3.2879,1.4733,2.3704,5.,5.,1.05,2.,0.8037,1.,3.,1.,3.,0.1,0.1,0.1 /)
rpcoef=0.0832
rs20=(/ 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0. /)
tminvj=(/ -15.,-15.,5.,5.,5.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15. /)
tmaxvj=(/ -10.,-10.,10.,15.,10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10. /)
vbeta=(/ 2.,2.,2.,2.,2.,4.,4.,4.,4.,4.,2.,2.,2.,2.,4.,4.,4. /)
rootbeta=(/ 0.943,0.962,0.966,0.961,0.966,0.914,0.964,0.972,0.943,0.943,0.961,0.961,0.961,0.961,0.961,0.975,0.975 /)
tcplant(:,1)=(/ 200.  ,300.  ,200.  ,300.  ,200.  ,150. ,150. ,250. ,250. ,250. ,250.,150.,0.1,150.,0.,1.,0. /)
tcplant(:,2)=(/ 10217.,16833.,5967. ,12000.,10217.,5000.,5000.,5247.,5247.,0.   ,0.  ,0.  ,0. ,0.  ,0.,0.,0. /)
tcplant(:,3)=(/ 876.  ,1443. ,511.  ,1029. ,876.  ,500. ,500. ,1124.,1124.,500. ,500.,607.,0.1,607.,0.,1.,0. /)
tcsoil(:,1)=(/ 184.   ,303.  ,107.  ,216.  ,184.  ,100. ,100. ,275. ,275. ,275. ,275.,149.,0.1,149.,0.,1.,0. /)
tcsoil(:,2)=(/ 367.   ,606.  ,214.  ,432.  ,367.  ,250. ,250. ,314. ,314. ,314. ,314.,300.,0.1,300.,0.,1.,0. /)
ratecp(1)=1.
ratecp(2)=0.03
ratecp(3)=0.14
ratecs(1)=2.
ratecs(2)=0.5

! read CABLE biome and LAI data
if (myid==0) then
  write(6,*) "Reading tiled surface data for CABLE"
  call vegta(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
else
  call vegtb(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
end if
do n=1,5
  svs(:,n)=svs(:,n)/sum(svs,2)
end do

icycle=ccycle
cable_user%fwsoil_switch="standard"

if (myid==0) write(6,*) "Define CABLE and CASA CNP arrays"

! default values (i.e., no land)  
ivegt=0
albsoilsn=0.08  
albsoil=0.08
albvisdir=0.08
albvisdif=0.08
albnirdir=0.08
albnirdif=0.08
zolnd=0.
vlai=0.
cplant=0.
csoil=0.
pind=ifull+1
mvtype=mxvt
mstype=mxst

! calculate length of CABLE vectors
mp=0
mland=0
do iq=1,ifull
  if (land(iq)) then
    mp=mp+count(svs(iq,:).gt.0.)
    mland=mland+1
  end if
end do
  
! if CABLE is present on this processor, then start allocating arrays
if (mp.gt.0) then
  
  allocate(sv(mp),cmap(mp))
  allocate(vl1(mp),vl2(mp),vl3(mp))
  call alloc_cbm_var(air, mp)
  call alloc_cbm_var(bgc, mp)
  call alloc_cbm_var(canopy, mp)
  call alloc_cbm_var(met, mp)
  call alloc_cbm_var(bal, mp)
  call alloc_cbm_var(rad, mp)
  call alloc_cbm_var(rough, mp)
  call alloc_cbm_var(soil, mp)
  call alloc_cbm_var(ssoil, mp)
  call alloc_cbm_var(sum_flux, mp)
  call alloc_cbm_var(veg, mp)

  ! soil parameters
  soil%zse = zse ! soil layer thickness
  soil%zshh(1)    = 0.5 * soil%zse(1)
  soil%zshh(ms+1) = 0.5 * soil%zse(ms)
  soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
  
  ! froot is now calculated from soil depth and the new parameter rootbeta 
  ! according to Jackson et al. 1996, Oceologica, 108:389-411
  !totdepth = 0.
  !do k=1,ms
  !  totdepth = totdepth + soil%zse(k)*100.
  !  froot2(:,k) = min(1.,1.-rootbeta(:)**totdepth)
  !enddo
  !do k = ms-1, 2, -1
  !  froot2(:,k) = froot2(:,k) - froot2(:,k-1)
  !enddo
  !froot2(:,ms)=1.-sum(froot(:,1:ms-1),2)
  
  froot2(:,1)=0.05
  froot2(:,2)=0.20
  froot2(:,3)=0.20
  froot2(:,4)=0.20
  froot2(:,5)=0.20
  froot2(:,6)=0.15
 
  sv=0.
  vl1=0.
  vl2=0.
  vl3=0.

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  ipos=0
  do n=1,5
    dumr=rlatt(:)*180./pi
    call getc4(ifull,ivs(:,n),dumr,c4frac)
    pind(n,1)=ipos+1
    do iq=1,ifull
      if (land(iq).and.(svs(iq,n).gt.0.)) then
        ipos=ipos+1
        if (ivs(iq,n).lt.1) then
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq=",myid,iq
          stop
        end if
        cmap(ipos)=iq
        sv(ipos)=svs(iq,n)
        veg%iveg(ipos)=ivs(iq,n)
        soil%isoilm(ipos)=isoilm(iq)
        veg%frac4(ipos)=c4frac(iq)
        if (fvegprev.ne.' '.and.fvegnext.ne.' ') then
          vlin(iq,n)=vlin(iq,n)+vlinprev(iq,n)
          vlinnext(iq,n)=vlinnext(iq,n)+vlin(iq,n)
          vl1(ipos)=0.5*vlin(iq,n)
          vl2(ipos)=4.*vlin(iq,n)-5.*vlinprev(iq,n)-vlinnext(iq,n)
          vl3(ipos)=1.5*(vlinnext(iq,n)+3.*vlinprev(iq,n)-3.*vlin(iq,n))
        else
          vl1(ipos)=vlin(iq,n)
          vl2(ipos)=0.
          vl3(ipos)=0.
        end if
        if (veg%iveg(ipos).eq.15.or.veg%iveg(ipos).eq.16) then
          vl1(ipos)=0.001
          vl2(ipos)=0.
          vl3(ipos)=0.
        end if
      end if
    end do
    pind(n,2)=ipos
  end do
  
  if (ipos.ne.mp) then
    write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
    stop
  end if

  ! Load CABLE arrays
  ivegt=ivs(:,1) ! diagnostic
  veg%meth      = 1
  veg%hc        = hc(veg%iveg)
  veg%canst1    = canst1(veg%iveg)
  veg%ejmax     = ejmax(veg%iveg)
  veg%tminvj    = tminvj(veg%iveg)
  veg%tmaxvj    = tmaxvj(veg%iveg)
  veg%vbeta     = vbeta(veg%iveg)
  veg%rp20      = rp20(veg%iveg)
  veg%rpcoef    = rpcoef(veg%iveg)
  veg%shelrb    = shelrb(veg%iveg)
  veg%vcmax     = vcmax(veg%iveg)
  veg%xfang     = xfang(veg%iveg)
  veg%dleaf     = dleaf(veg%iveg)
  veg%xalbnir   = 1. ! not used
  veg%taul(:,1) = taul(veg%iveg,1)
  veg%taul(:,2) = taul(veg%iveg,2)  
  veg%refl(:,1) = refl(veg%iveg,1)
  veg%refl(:,2) = refl(veg%iveg,2)  
  veg%extkn     = extkn(veg%iveg)
  veg%rs20      = rs20(veg%iveg)
  veg%vegcf     = vegcf(veg%iveg)
  !veg%wai       = wai(veg%iveg)
  do k=1,ms
    veg%froot(:,k)=froot2(veg%iveg,k)
  end do

  ! Calculate LAI and veg fraction diagnostics
  call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jyear,jmonth,jday,jhour,jmin)
  do n=1,5
    if (pind(n,1).le.mp) then
      vlai(cmap(pind(n,1):pind(n,2)))=vlai(cmap(pind(n,1):pind(n,2))) &
                                      +sv(pind(n,1):pind(n,2))*veg%vlai(pind(n,1):pind(n,2))
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
  bgc%ratecp(:) = ratecp(:)
  bgc%ratecs(:) = ratecs(:)

  ! store bare soil albedo and define snow free albedo
  soil%albsoil(:,1)=albvisnir(cmap,1)
  soil%albsoil(:,2)=albvisnir(cmap,2)
  soil%albsoil(:,3)=0.05
    
  where (land)
    albsoil(:)=0.5*sum(albvisnir,2)
  end where
  where ((albsoil.le.0.14).and.land)
    !sfact=0.5 for alb <= 0.14
    albsoilsn(:,1)=(1.00/1.50)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.50)*albsoil(:)
  elsewhere ((albsoil(:).le.0.2).and.land)
    !sfact=0.62 for 0.14 < alb <= 0.20
    albsoilsn(:,1)=(1.24/1.62)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.62)*albsoil(:)
  elsewhere (land)
    !sfact=0.68 for 0.2 < alb
    albsoilsn(:,1)=(1.36/1.68)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.68)*albsoil(:)
  end where
  ! MJT suggestion to get an approx inital albedo (before cable is called)
  where (land)
    albvisnir(:,1)=albsoilsn(:,1)*(1.-sigmf)+0.03*sigmf
    albvisnir(:,2)=albsoilsn(:,2)*(1.-sigmf)+0.20*sigmf
  end where
  albvisdir=albvisnir(:,1) ! To be updated by CABLE
  albvisdif=albvisnir(:,1) ! To be updated by CABLE
  albnirdir=albvisnir(:,2) ! To be updated by CABLE
  albnirdif=albvisnir(:,2) ! To be updated by CABLE
    
  ! MJT patch
  soil%albsoil(:,1)=albsoil(cmap)
  soil%albsoil(:,2)=albsoil(cmap)

  ssoil%albsoilsn(:,1)=albsoilsn(cmap,1) ! overwritten by CABLE
  ssoil%albsoilsn(:,2)=albsoilsn(cmap,2) ! overwritten by CABLE
  ssoil%albsoilsn(:,3)=0.05
  
  ssoil%t_snwlr=0.05
  ssoil%pudsmx=0.
  
  rad%albedo_T=albsoil(cmap)
  rad%trad=tss(cmap)
  rad%latitude=rlatt(cmap)*180./pi
  rad%longitude=rlongg(cmap)*180./pi

  canopy%oldcansto=0.  
  canopy%ghflux=0.
  canopy%sghflux=0.
  canopy%ga=0.
  canopy%dgdtg=0.
  canopy%fhs_cor=0.
  canopy%fes_cor=0.
  ssoil%wb_lake=0. ! not used when mlo.f90 is active
  canopy%ga=0.
  canopy%dgdtg=0.
  ssoil%fland=1.
  ssoil%ifland=soil%isoilm
  
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
  
  
  if (icycle==0) then
    ! Initialise CABLE carbon pools
    do k=1,ncp
      where (land)
        cplant(:,k)=tcplant(ivegt,k)
      end where
    end do
    do k=1,ncs
      where (land)        
        csoil(:,k)=tcsoil(ivegt,k)
      end where
    end do
  else
    ! CASA CNP
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    call alloc_phenvariable(phen,mp)
    
    casamet%lat=rad%latitude
    
    call casa_readpoint(casafile) ! read point sources

    leafage =(/ 2.0,1.5,1.0,1.0,1.4,1.4,1.0,1.0,1.0,0.8,1.0,     0.8,1.0,1.0,1.0,1.0,1.0 /)
    woodage =(/ 70.,60.,80.,40.,63.,63.,40.,40.,40.,1.0,1.0,     1.0,1.0,1.0,1.0,1.0,1.0 /)
    frootage=(/ 18.,10.,10.,10.,12.,12.,5.0,5.0,5.0,3.0,1.0,0.884227,1.0,1.0,1.0,1.0,1.0 /)
    metage=0.04
    strage=0.23
    cwdage=0.824
    micage=0.137
    slowage=5.
    passage=222.22
    clabileage=0.2

    xfherbivore   =(/ 0.068,0.406,0.068,0.134,0.169,0.169,0.022,0.022,0.022,0.109,0.000,0.140,0.000,0.000,0.000,0.000,0.000 /)
    xxkleafcoldmax=(/   0.2,  0.1,  0.1,  0.6, 0.25, 0.25,   1.,   1.,   1.,  0.2,  0.1,  0.3,  0.1,  0.1,  0.1,  0.1,  0.1 /)
    xxkleafdrymax =(/   0.1,  0.1,  0.1,   1.,0.325,0.325,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1 /)
    xratioNPleafmin =(/ 10.92308,15.95339,9.254839,12.73848,12.21745,12.21745,12.07217,12.07217,12.07217,13.51374,13.,15.12262,10.,10.,10.,10.,10. /)
    xratioNPleafmax =(/ 12.07288, 17.6327,10.22903,14.07938,13.50350,13.50350,12.07217,12.07217,12.07217,14.93733,13.,16.71447,10.,10.,10.,10.,10. /)
    xratioNPwoodmin =(/ 20.30167,15.89425,17.48344,19.08018,18.18989,18.18989,22.46035,22.46035,22.46035,     15.,15.,   20.52,15.,15.,15.,15.,15. /)
    xratioNPwoodmax =(/ 22.43869,17.56733, 19.3238,21.08862,20.10461,20.10461, 24.8246, 24.8246, 24.8246,     15.,15.,   20.52,15.,15.,15.,15.,15. /)
    xratioNPfrootmin=(/ 20.29341,15.87155,17.39767, 19.0601,18.15568,18.15568,22.39464,22.39464,22.39464,15.63498,15.,22.69109,15.,15.,15.,15.,15. /)
    xratioNPfrootmax=(/ 22.42955,17.54224,  19.229,21.06643,20.06681,20.06681,24.86138,24.86138,24.86138,17.28077,15.,25.07962,15.,15.,15.,15.,15. /)
    xfNminloss=0.05
    xfNminleach=0.05
    xnfixrate=(/ 0.08,2.6,0.21,1.64,1.13,1.13,0.37,0.37,0.37,0.95,0.,4.,0.,0.,0.,0.,0. /)
    xnsoilmin=1000.
    
    
    ratiocnplant(:,1)=(/  49.8, 23.1, 59.3, 31.4, 40.9, 40.9, 37.6, 37.6, 37.6,34.8, 30.,21.6, 40., 30., 40., 30., 40. /)
    ratiocnplant(:,2)=(/ 238.1,134.9,243.8,156.2,193.3,193.3,142.1,142.1,142.1,150.,150.,150.,150.,150.,150.,150.,135. /)
    ratiocnplant(:,3)=(/  73.7, 61.2,  75., 63.2, 68.3, 68.3, 67.1, 67.1, 67.1,64.5, 71.,60.7, 71., 71., 71., 71., 71. /)
    ratiocnsoil(:,1)=8.
    ratiocnsoil(:,2)=(/ 16.1,12.8,24.8, 30.,20.7,20.7,19.3,19.3,19.3,13.1,13.1,13.2, 20.,13.1, 20.,31.1, 20. /)
    ratiocnsoil(:,3)=(/ 16.1,12.8,24.8, 30.,20.7,20.7,19.3,19.3,19.3,13.1,13.1,13.2, 20.,13.1, 20.,31.1, 20. /)
    ratiocnsoilmin(:,1)=3.
    ratiocnsoilmin(:,2)=12.
    ratiocnsoilmin(:,3)=7.
    ratiocnsoilmax(:,1)=15.
    ratiocnsoilmax(:,2)=30.
    ratiocnsoilmax(:,3)=15.
     
    ! Initial values for CNP pools over 3*plant, 3*litter and 3*soil (=27 pools in total)
    cleaf =(/ 384.6037,    273.,96.59814,150.2638,226.1164,226.1164,     88.,     88.,     88.,137.1714,0.,    160.,0.,0.,0.,0.,0. /)
    cwood =(/ 7865.396,  11451.,5683.402,10833.74,8958.385,8958.385,    372.,    372.,    372.,      0.,0.,      0.,0.,0.,0.,0.,0. /)
    cfroot=(/     250.,   2586.,    220.,    220.,    819.,    819.,    140.,    140.,    140.,    263.,0.,    240.,0.,0.,0.,0.,0. /)
    cmet  =(/ 6.577021,44.63457,7.127119,10.97797,17.32917,17.32917,3.229374,3.229374,3.229374,28.57245,0.,28.57245,0.,0.,0.,0.,0. /)
    cstr  =(/ 209.1728,433.7626,277.7733,312.5492,308.3145,308.3145,39.44449,39.44449,39.44449,50.91091,0.,50.91091,0.,0.,0.,0.,0. /)
    ccwd  =(/ 606.0255,1150.765,776.7331,888.5864,855.5233,855.5233,111.5864,111.5864,111.5864,      0.,0.,      0.,0.,0.,0.,0.,0. /) 
    cmic  =(/  528.664,11.37765,597.0785,405.5554,385.6689,385.6689,168.0451,168.0451,168.0451,425.6431,0.,512.3247,0.,0.,0.,0.,0. /)
    cslow =(/ 13795.94,311.8092,16121.12,11153.25,10345.53,10345.53,4465.478,4465.478,4465.478,5694.437,0.,6855.438,0.,0.,0.,0.,0. /)
    cpass =(/ 4425.396,13201.81,5081.802,5041.192, 5937.55, 5937.55,1386.477,1386.477,1386.477, 4179.92,0.,5032.137,0.,0.,0.,0.,0. /)
    nleaf =(/ 7.541249,     9.9,1.609969,3.756594,5.701953,5.701953,2.933333,2.933333,2.933333,4.572381,0.,5.333333,0.,0.,0.,0.,0. /)
    nwood =(/ 31.46159,    102.,22.73361,80.24989,59.11127,59.11127,2.755555,2.755555,2.755555,      0.,0.,      0.,0.,0.,0.,0.,0. /)
    nfroot=(/ 6.097561,     38.,5.365854,5.365854,13.70732,13.70732,3.414634,3.414634,3.414634,6.414634,0.,5.853659,0.,0.,0.,0.,0. /)
    nmet  =(/ 0.064481, 0.74391,0.059393,0.137225,0.251252,0.251252,0.053823,0.053823,0.053823,0.476208,0.,0.476208,0.,0.,0.,0.,0. /)
    nstr  =(/ 1.394485,2.891751,1.851822,2.083661, 2.05543, 2.05543,0.262963,0.262963,0.262963,0.339406,0.,0.339406,0.,0.,0.,0.,0. /)
    ncwd  =(/ 2.424102,8.524183,3.106932,6.581996,5.159303,5.159303,0.826566,0.826566,0.826566,      0.,0.,      0.,0.,0.,0.,0.,0. /)
    nmic  =(/  52.8664,1.137765,59.70785,40.55554,38.56689,38.56689,16.80451,16.80451,16.80451,42.56431,0.,51.23247,0.,0.,0.,0.,0. /)
    nslow =(/ 919.7293,20.78728,1074.741,743.5501,689.7019,689.7019,297.6985,297.6985,297.6985,379.6291,0.,457.0292,0.,0.,0.,0.,0. /)
    npass =(/ 295.0264,880.1209,338.7868,336.0795,462.5034,462.5034, 92.4318, 92.4318, 92.4318,278.6613,0.,335.4758,0.,0.,0.,0.,0. /)
    xpleaf =(/ 0.191648,   0.415,0.115988,0.135453,0.214522,0.214522,0.022821,0.022821,0.022821, 0.15125,0., 0.15125,0.,0.,0.,0.,0. /)
    xpwood =(/ 0.953979,    5.88, 0.64438,2.424778,2.476034,2.476034,      0.,      0.,      0.,      0.,0.,      0.,0.,0.,0.,0.,0. /)
    xpfroot=(/ 0.076659,    1.95,0.080548,0.141097,0.562076,0.562076,0.037083,0.037083,0.037083, 0.15125,0., 0.15125,0.,0.,0.,0.,0. /)
    xpmet  =(/ 0.004385,0.029756,0.004751,0.007319,0.011553,0.011553,0.002153,0.002153,0.002153,0.019048,0.,0.019048,0.,0.,0.,0.,0. /)
    xpstr  =(/ 0.069724,0.144588,0.092591,0.104183,0.102772,0.102772,0.013148,0.013148,0.013148, 0.01697,0., 0.01697,0.,0.,0.,0.,0. /)
    xpcwd  =(/ 0.101004,0.191794,0.129456,0.148095,0.142587,0.142587,0.018598,0.018598,0.018598,      0.,0.,      0.,0.,0.,0.,0.,0. /)
    xpmic  =(/ 6.872632, 0.14791,7.762021, 5.27222,5.013696,5.013696,2.184586,2.184586,2.184586,5.533361,0.,6.661522,0.,0.,0.,0.,0. /)
    xpslow =(/ 119.5648,2.702347,139.7164,96.66152,89.66127,89.66127,38.70081,38.70081,38.70081,49.35178,0., 59.4138,0.,0.,0.,0.,0. /)
    xppass =(/ 38.35343,114.4157,44.04228,43.69033,60.12544,60.12544,12.01613,12.01613,12.01613,36.22598,0.,43.61185,0.,0.,0.,0.,0. /)
    xplab =(/  26.737, 19.947, 29.107, 30.509, 26.575, 26.575, 23.206, 23.206, 23.206, 25.538, 0., 27.729, 0., 0., 0.103, 0., 0. /)
    xpsorb=(/  126.73, 92.263,134.639,132.012,121.411,121.411, 173.47, 173.47, 173.47,186.207, 0.,155.518, 0., 0., 1.176, 0., 0. /)
    xpocc =(/ 138.571,120.374, 138.22,148.083,136.312,136.312,114.496,114.496,114.496,145.163, 0.,158.884, 0., 0., 0.688, 0., 0. /)
 
    xkmlabp  =(/ 74.5408, 68.1584,  77.952,64.41918,64.41918,70.5856, 64.5888,54.1692, 9.7704, 28.29,  63.963,  32.402 /)
    xpsorbmax=(/ 745.408,788.0815,1110.816, 744.847, 744.847,816.146,746.8081,722.256,293.112,311.19,373.1175,615.6381 /)
    xfPleach =0.0005
    rationpsoil(:,1)=4.
    rationpsoil(:,2)=(/ 5.,5.,5.,15.,5.,5.,5.,5.,7.,7.,7.,7. /)
    rationpsoil(:,3)=(/ 5.,5.,5.,15.,5.,5.,5.,5.,7.,7.,7.,7. /)
 
    casabiome%ivt2     =(/   3,  3,  3,  3,  3,  3,  2,  2,  2,  1,  0,  1,  0,  1,  0,  1,  0 /)
    casabiome%kroot    =(/ 5.5,3.9,5.5,3.9,4.7,4.7,2.0,2.0,2.0,5.5,5.5,5.5,2.0,5.5,5.5,5.5,5.5 /)
    casabiome%rootdepth=(/ 1.5,1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.5 /)
    casabiome%kuptake  =(/ 2.0,1.9,2.0,2.0,2.0,2.0,1.8,1.8,1.8,2.0,1.6,1.6,1.8,1.8,1.8,1.8,1.8 /)
    casabiome%krootlen =(/ 14.87805,14.38596,14.02597,18.94737,15.55933,15.55933,32.30769,32.30769,32.30769,84., &
                           0.,120.5,0.,0.,0.,0.,0. /)
    casabiome%kminN=2
    casabiome%kuplabP=0.5
    casabiome%fracnpptoP(:,leaf) =(/ 0.25,0.20,0.40,0.35,0.30,0.30,0.35,0.35,0.35,0.35,0.50,0.50,0.50,0.50,0.50,0.50,0.60 /)
    casabiome%fracnpptoP(:,wood) =(/ 0.40,0.35,0.30,0.25,0.33,0.33,0.25,0.25,0.25,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.40 /)
    casabiome%fracnpptoP(:,xroot)=(/ 0.35,0.45,0.30,0.40,0.38,0.38,0.40,0.40,0.40,0.65,0.50,0.50,0.50,0.50,0.50,0.50,0.00 /) 
    casabiome%rmplant(:,leaf)    =0.1
    casabiome%rmplant(:,wood)    =(/ 2.0,1.0,1.5,0.8,1.3,1.3,0.5,0.5,0.5,0.5,1.0,2.0,1.0,1.0,1.0,1.0,1.0 /)
    casabiome%rmplant(:,xroot)   =(/ 10.,2.0,7.5,2.5,5.5,5.5,4.5,4.5,4.5,4.5,10.,25.,10.,10.,10.,10.,10. /)
    casabiome%ftransNPtoL(:,leaf) =0.5
    casabiome%ftransNPtoL(:,wood) =0.95
    casabiome%ftransNPtoL(:,xroot)=0.9
    casabiome%fracligninplant(:,leaf) =(/ 0.25,0.20,0.20,0.20,0.21,0.21,0.20,0.20,0.20,0.10,0.15,0.10,0.15,0.15,0.10,0.15,0.25 /)
    casabiome%fracligninplant(:,wood) =0.4
    casabiome%fracligninplant(:,xroot)=(/ 0.25,0.20,0.20,0.20,0.21,0.21,0.20,0.20,0.20,0.10,0.15,0.10,0.15,0.15,0.10,0.15,0.25 /)
    casabiome%glaimax=(/ 7.,7.,7.,7.,7.,7.,3.,3.,3.,3., 5.,6., 6., 5.,0., 5., 1. /)
    casabiome%glaimin=(/ 1.,1.,.5,.5,.8,.8,.1,.1,.1,.1,.05,.1,.05,.05,0.,.05,.05 /)
    phen%TKshed=(/ 268.,260.,263.15,268.15,264.83,264.83,277.15,277.15,277.15,275.15,277.15,278.15,277.15,277.15,283.15,277.15,277.15 /)
    casabiome%xkleafcoldexp=3.
    casabiome%xkleafdryexp=3.
    casabiome%ratioNCplantmin(:,leaf) =(/     0.02,    0.04,0.016667,0.028571, 0.02631, 0.02631,   0.025,   0.025,   0.025,0.026316,0.033333,    0.04,   0.025,   0.025,   0.025,   0.025,   0.025 /)
    casabiome%ratioNCplantmax(:,leaf) =(/    0.024,   0.048,    0.02,0.034286,0.031572,0.031572,    0.03,    0.03,    0.03,0.031579,    0.04,   0.048,    0.03,    0.03,    0.03,    0.03,    0.03 /)
    casabiome%ratioNCplantmin(:,wood) =(/    0.004,0.006667,   0.004,0.005714,0.005095,0.005095,0.006667,0.006667,0.006667,0.006667,0.006667,   0.008,0.006667,0.006667,0.006667,0.006667,0.007407 /)
    casabiome%ratioNCplantmax(:,wood) =(/   0.0048,   0.008,  0.0048,0.006857,0.006114,0.006114,   0.008,   0.008,   0.008,   0.008,   0.008,  0.0096,   0.008,   0.008,   0.008,   0.008,0.008889 /)
    casabiome%ratioNCplantmin(:,xroot)=(/ 0.012821,0.014706,0.012821,0.014085,0.013608,0.013608,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085,0.014085 /)
    casabiome%ratioNCplantmax(:,xroot)=(/ 0.015385,0.017647,0.015385,0.016901, 0.01633, 0.01633,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901,0.016901 /)
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
    casamet%lnonwood = 1

    casapool%cplant(:,wood)  = 0.
    casapool%clitter(:,cwd)  = 0.
    casapool%nplant(:,wood)  = 0.
    casapool%nlitter(:,cwd)  = 0.
    casapool%pplant(:,wood)  = 0.
    casapool%plitter(:,cwd)  = 0.
    where (casamet%iveg2==3.or.casamet%iveg2==2)
      casamet%lnonwood = 0
      casapool%cplant(:,wood)  = cwood(veg%iveg) 
      casapool%clitter(:,cwd)  = ccwd(veg%iveg)
      casapool%nplant(:,wood)  = nwood(veg%iveg) 
      casapool%nlitter(:,cwd)  = ncwd(veg%iveg)
      casapool%pplant(:,wood)  = xpwood(veg%iveg)
      casapool%plitter(:,cwd)  = xpcwd(veg%iveg)
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
    casamet%glai = max(casabiome%glaimin(veg%iveg),casabiome%sla(veg%iveg)*casapool%cplant(:,leaf))
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
    casapool%plitter(:,str) = cstr(veg%iveg)*ratioPCstrfix
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
    call casa_readphen(fphen) ! read MODIS leaf phenology
    casamet%tairk     = 0.
    casamet%tsoil     = 0.
    casamet%moist     = 0.
    casaflux%cgpp     = 0.
    casaflux%Crsoil   = 0.
    casaflux%crgplant = 0.
    casaflux%crmplant = 0.
    casaflux%clabloss = 0.
  end if ! icycle>0

end if
  
if (myid==0) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine loadcbmparm

! *************************************************************************************
! define C4 crops.  Since we do not have a high-resolution dataset, instead we
! diagnose C4 crops from biome and location.
subroutine getc4(ifull,ivegt,rlatt,c4frac)
  
implicit none
  
integer, intent(in) :: ifull
integer, dimension(ifull), intent(in) :: ivegt
real, dimension(ifull), intent(in) :: rlatt
real, dimension(ifull), intent(out) :: c4frac

where ((ivegt.eq.7).and.(rlatt.ge.-30.).and.(rlatt.le.30.))
  c4frac=0.95
elsewhere ((ivegt.eq.8).and.(rlatt.ge.-30.).and.(rlatt.le.0.))
  c4frac=0.5
elsewhere ((ivegt.eq.8).and.(rlatt.ge.0.).and.(rlatt.le.20.))
  c4frac=0.8
elsewhere ((ivegt.eq.8).and.(rlatt.ge.20.).and.(rlatt.le.30.))
  c4frac=0.5
elsewhere (ivegt.eq.9)
  c4frac=0.75
elsewhere ((ivegt.eq.10).and.(rlatt.ge.-30.).and.(rlatt.le.-20.))
  c4frac=0.5
elsewhere ((ivegt.eq.10).and.(rlatt.ge.-20.).and.(rlatt.le.20.))
  c4frac=0.95
elsewhere ((ivegt.eq.10).and.(rlatt.ge.20.).and.(rlatt.le.35.))
  c4frac=0.5
elsewhere ((ivegt.eq.12).and.(rlatt.ge.0.).and.(rlatt.le.40.))
  c4frac=0.3
elsewhere
  c4frac=0.
end where
  
return
end subroutine getc4

! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi

implicit none
  
include 'newmpar.h'
include 'parmgeom.h'  ! rlong0,rlat0,schmidt  
  
character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer(i_d), dimension(ifull,5), intent(out) :: ivs
integer(i_d), dimension(ifull_g,5) :: ivsg  
integer n,iq,ilx,jlx,iad  
real(r_1), dimension(ifull,5), intent(out) :: svs,vlinprev,vlin,vlinnext
real(r_1), dimension(ifull_g,5) :: svsg,vling
real rlong0x,rlat0x,schmidtx,dsx,ra,rb
character(len=47) header  

write(6,*) "Reading land-use data for CABLE"
if (fvegprev.ne.' '.and.fvegnext.ne.' ') then
  open(87,file=fvegprev,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt) stop 'wrong data file supplied'
  do iq=1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  do n=1,5
    call ccmpi_distribute(vlinprev(:,n),vling(:,n))
  end do
  open(87,file=fvegnext,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt) stop 'wrong data file supplied'
  do iq=1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  do n=1,5
    call ccmpi_distribute(vlinnext(:,n),vling(:,n))
  end do      
else
  vlinprev=-1.
  vlinnext=-1.    
end if
open(87,file=fveg,status='old')
read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt) stop 'wrong data file supplied'
do iq=1,ifull_g
  read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
             ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
end do
close(87)
do n=1,5
  call ccmpi_distribute(ivs(:,n),ivsg(:,n))
  call ccmpi_distribute(svs(:,n),svsg(:,n))
  call ccmpi_distribute(vlin(:,n),vling(:,n))
end do
  
return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi
  
implicit none

include 'newmpar.h'

character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer(i_d), dimension(ifull,5), intent(out) :: ivs
integer n,iq
real(r_1), dimension(ifull,5), intent(out) :: svs,vlinprev,vlin,vlinnext

if (fvegprev.ne.' '.and.fvegnext.ne.' ') then
  do n=1,5
    call ccmpi_distribute(vlinprev(:,n))
  end do
  do n=1,5
    call ccmpi_distribute(vlinnext(:,n))
  end do
else
  vlinprev=-1.
  vlinnext=-1.
end if    
do n=1,5
  call ccmpi_distribute(ivs(:,n))
  call ccmpi_distribute(svs(:,n))
  call ccmpi_distribute(vlin(:,n))
end do
  
return
end subroutine vegtb

! *************************************************************************************  
! This subroutine loads CABLE tile data
subroutine loadtile

use carbpools_m
use cc_mpi, only : myid
use infile
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none

include 'newmpar.h'
include 'darcdf.h'
include 'mpif.h'
include 'netcdf.inc'
include 'parm.h'  
  
integer k,n,ierr,ierr2,idv
real, dimension(ifull) :: dat
real totdepth
character(len=9) vname
  
if (io_in.eq.1) then
  if (myid.eq.0) idv = ncvid(ncid,"tgg1_1",ierr)
  call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr2)    
else
  ierr=1
end if
  
! Cannot locate tile data, use diagnostic data instead
if (ierr.ne.0) then
  if (myid==0) write(6,*) "Use gridbox averaged data to initialise CABLE"
  if (mp.gt.0) then
    do k = 1,ms
      ssoil%tgg(:,k) = tgg(cmap,k)
      ssoil%wb(:,k) = wb(cmap,k)
      ssoil%wbice(:,k) = wbice(cmap,k)
    enddo
    do k = 1,3
      ssoil%tggsn(:,k) = tggsn(cmap,k)
      ssoil%smass(:,k) = smass(cmap,k)
      ssoil%ssdn(:,k) = ssdn(cmap,k)
    enddo      
    ssoil%isflag=isflag(cmap)
    ssoil%snowd=snowd(cmap)
    ssoil%snage=snage(cmap)
    ssoil%rtsoil=100.
    canopy%cansto=0.
    ssoil%pudsto=0.
    ssoil%wetfac=0.
    if (icycle==0) then
      do k=1,ncp
        bgc%cplant(:,k) = cplant(cmap,k)
      enddo
      do k=1,ncs
        bgc%csoil(:,k) = csoil(cmap,k)
      enddo
    else
      do k=1,mplant
        casapool%cplant(:,k)=cplant(cmap,k)
        casapool%nplant(:,k)=niplant(cmap,k)
        casapool%pplant(:,k)=pplant(cmap,k)
      end do
      do k=1,mlitter
        casapool%clitter(:,k)=clitter(cmap,k)
        casapool%nlitter(:,k)=nilitter(cmap,k)
        casapool%plitter(:,k)=plitter(cmap,k)
      end do
      do k=1,msoil
        casapool%csoil(:,k)=csoil(cmap,k)
        casapool%nsoil(:,k)=nisoil(cmap,k)
        casapool%psoil(:,k)=psoil(cmap,k)
      end do
      casamet%glai=glai(cmap)
      ! may need the following fields if allowing prognostic LAI to feedback into CABLE
      !phen%phase=
      !casapool%clabile=
      !casapool%nsoilmin=
      !casapool%psoillab=
      !casapool%psoilsorb=
      !casapool%psoilocc=
    end if
  end if
else
  ! Located CABLE tile data
  if (myid==0) write(6,*) "Use tiled data to initialise CABLE"
  do n=1,5
    do k=1,ms
      write(vname,'("tgg",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%tgg(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("wb",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%wb(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("wbice",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%wbice(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
    end do
    do k=1,3
      write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%tggsn(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("smass",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%smass(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) ssoil%ssdn(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
    end do
    write(vname,'("sflag_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%isflag(pind(n,1):pind(n,2))=nint(dat(cmap(pind(n,1):pind(n,2))))
    write(vname,'("snd_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%snowd(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    write(vname,'("snage_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%snage(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    write(vname,'("rtsoil_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%rtsoil(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    write(vname,'("cansto_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) canopy%cansto(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    write(vname,'("pudsto_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%pudsto(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    write(vname,'("wetfac_",I1.1)') n
    call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
    if (pind(n,1).le.mp) ssoil%wetfac(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    if (icycle==0) then
      do k=1,ncp
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) bgc%cplant(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2))) 
      enddo
      do k=1,ncs
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) bgc%csoil(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
      enddo
    else
      do k=1,mplant
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%cplant(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("nplant",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%nplant(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("pplant",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%pplant(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
      end do
      do k=1,mlitter
        write(vname,'("clitter",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%clitter(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("nlitter",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%nlitter(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("plitter",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%plitter(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
      end do
      do k=1,msoil
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%csoil(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("nsoil",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%nsoil(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("psoil",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
        if (pind(n,1).le.mp) casapool%psoil(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
      end do
      write(vname,'("glai_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casamet%glai(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("phenphase_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) phen%phase(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("clabile_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casapool%clabile(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("nsoilmin_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casapool%nsoilmin(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("psoillab_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casapool%psoillab(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("psoilsorb_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casapool%psoilsorb(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("psoilocc_",I1.1)') n
      call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,dat,ifull)
      if (pind(n,1).le.mp) casapool%psoilocc(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
    end if
  end do
  ! albvisdir, albvisdif, albnirdir, albnirdif are used when nrad=5
  vname='albvisdir'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albvisdir,ifull)
  vname='albvisdif'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albvisdif,ifull)
  vname='albnirdir'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albnirdir,ifull)
  vname='albnirdif'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albnirdif,ifull)
  ! albvis and albnir are used when nrad=4
  vname='albvis'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albvisnir(:,1),ifull)
  vname='albnir'
  call histrd1(ncid,iarchi-1,ierr,vname,il_g,jl_g,albvisnir(:,2),ifull)        
end if
  
! Some fixes for rounding errors
if (mp.gt.0) then

  totdepth = 0.
  do k=1,ms
    totdepth = totdepth + soil%zse(k)*100.
  enddo

  ssoil%wb=max(ssoil%wb,0.)
  ssoil%wbice=max(ssoil%wbice,0.)
  ssoil%smass=max(ssoil%smass,0.)
  ssoil%rtsoil=max(ssoil%rtsoil,0.)
  ssoil%snowd=max(ssoil%snowd,0.)
  ssoil%wetfac=min(max(ssoil%wetfac,0.),1.)
  canopy%cansto=max(canopy%cansto,0.)

  ! overwritten by CABLE
  ssoil%osnowd=ssoil%snowd                                ! overwritten by CABLE
  bal%osnowd0=ssoil%snowd                                 ! overwritten by CABLE
  where (ssoil%isflag.gt.0)
    ssoil%sdepth(:,1)=ssoil%smass(:,1)/ssoil%ssdn(:,1)    ! overwritten by CABLE
    ssoil%ssdnn=(ssoil%ssdn(:,1)*ssoil%smass(:,1)+ssoil%ssdn(:,2) &
         & *ssoil%smass(:,2)+ssoil%ssdn(:,3)*ssoil%smass(:,3))    &
         & /ssoil%snowd
  elsewhere
    ssoil%sdepth(:,1)=ssoil%snowd/ssoil%ssdn(:,1)         ! overwritten by CABLE
    ssoil%ssdnn=max(120.,ssoil%ssdn(:,1))                 ! overwritten by CABLE
  end where
  do k=2,3
    where (ssoil%isflag.gt.0)
      ssoil%sdepth(:,k)=ssoil%smass(:,k)/ssoil%ssdn(:,k)  ! overwritten by CABLE
    elsewhere
      ssoil%sdepth(:,k)=ssoil%sconds(:,1)                 ! overwritten by CABLE
    end where
  end do  
  ssoil%owetfac=ssoil%wetfac
  ssoil%wbtot=0.
  ssoil%wbtot1=0.
  ssoil%wbtot2=0.
  ssoil%tggav=0.
  do k = 1,ms
    ssoil%wbtot=ssoil%wbtot+ssoil%wb(:,k)*1000.0*soil%zse(k)
    ssoil%tggav=ssoil%tggav+soil%zse(k)*ssoil%tgg(:,k)/(totdepth/100.)
    ssoil%gammzz(:,k)=max((1.-soil%ssat)*soil%css* soil%rhosoil &
       & +(ssoil%wb(:,k)-ssoil%wbice(:,k))*4.218e3* 1000.       &
       & +ssoil%wbice(:,k)*2.100e3*1000.*0.9,soil%css*soil%rhosoil)*soil%zse(k)
  end do
  bal%wbtot0=ssoil%wbtot

  if (icycle==0) then
    bgc%cplant=max(bgc%cplant,0.)
    bgc%csoil=max(bgc%csoil,0.)
  else
    casapool%cplant     = max(0.,casapool%cplant)
    casapool%clitter    = max(0.,casapool%clitter)
    casapool%csoil      = max(0.,casapool%csoil)
    casabal%cplantlast  = casapool%cplant
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal     = 0.
    casabal%FCgppyear   = 0.
    casabal%FCrpyear    = 0.
    casabal%FCnppyear   = 0.
    casabal%FCrsyear    = 0.
    casabal%FCneeyear   = 0.
    casapool%nplant     = max(1.e-6,casapool%nplant)
    casapool%nlitter    = max(1.e-6,casapool%nlitter)
    casapool%nsoil      = max(1.e-6,casapool%nsoil)
    casapool%nsoilmin   = max(1.e-6,casapool%nsoilmin)
    casabal%nplantlast  = casapool%nplant
    casabal%nlitterlast = casapool%nlitter
    casabal%nsoillast   = casapool%nsoil       
    casabal%nsoilminlast= casapool%nsoilmin
    casabal%sumnbal     = 0.
    casabal%FNdepyear   = 0.
    casabal%FNfixyear   = 0.
    casabal%FNsnetyear  = 0.
    casabal%FNupyear    = 0.
    casabal%FNleachyear = 0.
    casabal%FNlossyear  = 0.
    casapool%pplant     = max(1.0e-7,casapool%pplant)
    casapool%plitter    = max(1.0e-7,casapool%plitter)
    casapool%psoil      = max(1.0e-7,casapool%psoil)
    casapool%Psoillab   = max(1.0e-7,casapool%psoillab)
    casapool%psoilsorb  = max(1.0e-7,casapool%psoilsorb)
    casapool%psoilocc   = max(1.0e-7,casapool%psoilocc)
    casabal%pplantlast  = casapool%pplant
    casabal%plitterlast = casapool%plitter
    casabal%psoillast   = casapool%psoil       
    casabal%psoillablast= casapool%psoillab
    casabal%psoilsorblast=casapool%psoilsorb
    casabal%psoilocclast= casapool%psoilocc
    casabal%sumpbal     = 0.
    casabal%FPweayear   = 0.
    casabal%FPdustyear  = 0.
    casabal%FPsnetyear  = 0.
    casabal%FPupyear    = 0.
    casabal%FPleachyear = 0.
    casabal%FPlossyear  = 0.
  end if
end if
  
return
end subroutine loadtile

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,idim,iarch)

use carbpools_m
use cc_mpi, only : myid
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none

include 'newmpar.h'
  
integer, intent(in) :: idnc,iarch
integer k,n,ierr
integer, dimension(3), intent(in) :: idim  
real, dimension(ifull) :: dat
character(len=9) vname
character(len=40) lname
logical, intent(in) :: local
  
if (myid.eq.0) then
  write(6,*) "Storing CABLE tile data"
  call ncredf(idnc,ierr)
  do n=1,5
    do k=1,ms
      write(lname,'("Soil temperature lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("tgg",I1.1,"_",I1.1)') k,n
      call attrib(idnc,idim,3,vname,lname,'K',100.,400.,0,-1)
      write(lname,'("Soil moisture lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("wb",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'m3/m3',0.,2.6,0,-1)
      write(lname,'("Soil ice lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("wbice",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'m3/m3',0.,2.6,0,-1)
    end do
    do k=1,3
      write(lname,'("Snow temperature lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
      call attrib(idnc,idim,3,vname,lname,'K',100.,400.,0,-1)
      write(lname,'("Snow mass lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("smass",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'K',0.,650.,0,-1)
      write(lname,'("Snow density lev ",I1.1," tile ",I1.1)') k,n
      write(vname,'("ssdn",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'kg/m3',0.,650.,0,-1)
    end do
    write(lname,'("Snow flag tile ",I1.1)') n
    write(vname,'("sflag_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'mm',0.,6.5,0,-1)
    write(lname,'("Snow depth tile ",I1.1)') n
    write(vname,'("snd_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'mm',0.,6500.,0,-1)  ! -1=long
    write(lname,'("Snow age tile ",I1.1)') n
    write(vname,'("snage_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'none',0.,26.,0,-1)
    write(lname,'("Soil turbulent resistance tile ",I1.1)') n
    write(vname,'("rtsoil_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'none',0.,1.3e5,0,-1)
    write(lname,'("cansto tile ",I1.1)') n
    write(vname,'("cansto_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'none',0.,13.,0,-1)
    write(lname,'("pudsto tile ",I1.1)') n
    write(vname,'("pudsto_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'none',0.,13.,0,-1)
    write(lname,'("wetfac tile ",I1.1)') n
    write(vname,'("wetfac_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'none',0.,6.5,0,-1)
    if (icycle==0) then
      write(lname,'("Carbon leaf pool tile ",I1.1)') n
      write(vname,'("cplant1_",I1.1)') n    
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Carbon wood pool tile ",I1.1)') n
      write(vname,'("cplant2_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Carbon root pool tile ",I1.1)') n
      write(vname,'("cplant3_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Carbon soil fast pool tile ",I1.1)') n
      write(vname,'("csoil1_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Carbon soil slow pool tile ",I1.1)') n
      write(vname,'("csoil2_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
    else
      do k=1,mplant
        write(lname,'("C leaf pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N leaf pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("nplant",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P leaf pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("pplant",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      end do
      do k=1,mlitter
        write(lname,'("C litter pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("clitter",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N litter pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("nlitter",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P litter pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("plitter",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      end do
      do k=1,msoil
        write(lname,'("C soil pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("N soil pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("nsoil",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
        write(lname,'("P soil pool ",I1.1," tile ",I1.1)') k,n
        write(vname,'("psoil",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      end do
      write(lname,'("Prognostic LAI tile ",I1.1)') n
      write(vname,'("glai_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("Leaf phenology phase tile ",I1.1)') n
      write(vname,'("phenphase_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("C labile tile ",I1.1)') n
      write(vname,'("clabile_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("N soilmin tile ",I1.1)') n
      write(vname,'("nsoilmin_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soillab tile ",I1.1)') n
      write(vname,'("psoillab_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soilsorb tile ",I1.1)') n
      write(vname,'("psoilsorb_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
      write(lname,'("P soilocc tile ",I1.1)') n
      write(vname,'("psoilocc_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,65000.,0,-1)
    end if
  end do
  lname='DIR VIS albedo'
  vname='albvisdir'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  lname='DIF VIS albedo'
  vname='albvisdif'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  lname='DIR NIR albedo'
  vname='albnirdir'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  lname='DIF NIR albedo'
  vname='albnirdif'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  lname='VIS albedo'
  vname='albvis'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  lname='NIR albedo'
  vname='albnir'
  call attrib(idnc,idim,3,vname,lname,'none',0.,1.3,0,-1)
  call ncendf(idnc,ierr)
end if
do n=1,5
  do k=1,ms
    dat=tgg(:,k)
    if (pind(n,1).le.mp) then      
      dat(cmap(pind(n,1):pind(n,2)))=ssoil%tgg(pind(n,1):pind(n,2),k)
    end if
    write(vname,'("tgg",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wb(:,k)
    if (pind(n,1).le.mp) then  
      dat(cmap(pind(n,1):pind(n,2)))=ssoil%wb(pind(n,1):pind(n,2),k)
    end if
    write(vname,'("wb",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wbice(:,k)
    if (pind(n,1).le.mp) then  
      dat(cmap(pind(n,1):pind(n,2)))=ssoil%wbice(pind(n,1):pind(n,2),k)
    end if
    write(vname,'("wbice",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  do k=1,3
    dat=tggsn(:,k)
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%tggsn(pind(n,1):pind(n,2),k)
    write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=smass(:,k)
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%smass(pind(n,1):pind(n,2),k)
    write(vname,'("smass",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=ssdn(:,k)
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%ssdn(pind(n,1):pind(n,2),k)
    write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  dat=real(isflag)
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=real(ssoil%isflag(pind(n,1):pind(n,2)))
  write(vname,'("sflag_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snowd
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%snowd(pind(n,1):pind(n,2))
  write(vname,'("snd_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)  ! long write    
  dat=snage
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%snage(pind(n,1):pind(n,2))
  write(vname,'("snage_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=100.
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%rtsoil(pind(n,1):pind(n,2))
  write(vname,'("rtsoil_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)   
  dat=cansto
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=canopy%cansto(pind(n,1):pind(n,2))
  write(vname,'("cansto_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%pudsto(pind(n,1):pind(n,2))
  write(vname,'("pudsto_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=ssoil%wetfac(pind(n,1):pind(n,2))
  write(vname,'("wetfac_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  if (icycle==0) then
    do k=1,ncp
      dat=cplant(:,k)
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=bgc%cplant(pind(n,1):pind(n,2),k)
      write(vname,'("cplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)    
    end do
    do k=1,ncs
      dat=csoil(:,k)
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=bgc%csoil(pind(n,1):pind(n,2),k)
      write(vname,'("csoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
  else
    do k=1,mplant     
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%cplant(pind(n,1):pind(n,2),k)
      write(vname,'("cplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%nplant(pind(n,1):pind(n,2),k)
      write(vname,'("nplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%pplant(pind(n,1):pind(n,2),k)
      write(vname,'("pplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k=1,mlitter
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%clitter(pind(n,1):pind(n,2),k)
      write(vname,'("clitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%nlitter(pind(n,1):pind(n,2),k)
      write(vname,'("nlitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%plitter(pind(n,1):pind(n,2),k)
      write(vname,'("plitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k=1,msoil
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%csoil(pind(n,1):pind(n,2),k)
      write(vname,'("csoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%nsoil(pind(n,1):pind(n,2),k)
      write(vname,'("nsoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%psoil(pind(n,1):pind(n,2),k)
      write(vname,'("psoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casamet%glai(pind(n,1):pind(n,2))
    write(vname,'("glai_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=phen%phase(pind(n,1):pind(n,2))
    write(vname,'("phenphase_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%clabile(pind(n,1):pind(n,2))
    write(vname,'("clabile_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%nsoilmin(pind(n,1):pind(n,2))
    write(vname,'("nsoilmin_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%psoillab(pind(n,1):pind(n,2))
    write(vname,'("psoillab_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%psoilsorb(pind(n,1):pind(n,2))
    write(vname,'("psoilsorb_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (pind(n,1).le.mp) dat(cmap(pind(n,1):pind(n,2)))=casapool%psoilocc(pind(n,1):pind(n,2))
    write(vname,'("psoilocc_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end if
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

subroutine cableinflow(iqin,inflow)
  
implicit none
  
integer, intent(in) :: iqin
integer n,iq,i
real, intent(inout) :: inflow
real, dimension(5) :: xx
real yy
  
xx(:)=inflow
inflow=0.
do n=1,5
  iq=-1
  do i=pind(n,1),pind(n,2)
    if (cmap(i).eq.iqin) then
      iq=i
      exit
    elseif (cmap(i).gt.iqin) then
      exit
    end if
  end do
  if (iq.gt.0) then
    yy=min(xx(n),(soil%ssat(iq)-ssoil%wb(iq,cbm_ms))*1000.*soil%zse(cbm_ms))
    ssoil%wb(iq,cbm_ms)=ssoil%wb(iq,cbm_ms)+yy/(1000.*soil%zse(cbm_ms))
    xx(n)=max(xx(n)-yy,0.)
    inflow=inflow+sv(iq)*xx(n)
  else
    exit
  end if
end do
  
return
end subroutine cableinflow

! *************************************************************************************
! Transfer grid information from CABLE internally, read N&P input from
! integral NETCDF file
subroutine casa_readpoint(casafile)

use cc_mpi
use infile

implicit none
  
include 'newmpar.h'
include 'mpif.h'
include 'netcdf.inc'
include 'parmgeom.h'
  
integer ncstatus,ncid,varid,tilg
integer, dimension(2) :: spos,npos
real tlat,tlon,tschmidt
real, dimension(:), allocatable :: dumg
real, dimension(ifull) :: duma
character(len=*), intent(in) :: casafile

if (myid==0) then
  allocate(dumg(ifull_g))
  write(6,*) "Reading ",trim(casafile)
  ncstatus=nf_open(casafile,nf_nowrite,ncid)
  call ncmsg('CASA_readpoint',ncstatus)
  ! check dimensions and location
  ncstatus=nf_get_att_real(ncid,nf_global,'lat0',tlat)
  call ncmsg('lat0',ncstatus)
  ncstatus=nf_get_att_real(ncid,nf_global,'lon0',tlon)
  call ncmsg('lon0',ncstatus)
  ncstatus=nf_get_att_real(ncid,nf_global,'schmidt0',tschmidt)
  call ncmsg('schmidt0',ncstatus)
  if (rlong0.ne.tlon.or.rlat0.ne.tlat.or.schmidt.ne.tschmidt) then
    write(6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    stop
  end if
  ncstatus = nf_inq_dimid(ncid,'longitude',varid)
  call ncmsg('longitude',ncstatus)
  ncstatus = nf_inq_dimlen(ncid,varid,tilg)
  call ncmsg('longitude',ncstatus)
  if (tilg.ne.il_g) then
    write (6,*) "ERROR: Grid mismatch for ",trim(casafile)
    write (6,*) "il_g,tilg ",il_g,tilg
    stop
  end if
  ! load casa fields
  spos=1
  npos(1)=il_g
  npos(2)=il_g*6
  write(6,*) "Loading soil order"
  ncstatus = nf_inq_varid(ncid,'sorder',varid)
  call ncmsg('sorder',ncstatus)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ncmsg('sorder',ncstatus)
  call ccmpi_distribute(duma,dumg)
  casamet%isorder=duma(cmap)
  write(6,*) "Loading N deposition rate"
  ncstatus = nf_inq_varid(ncid,'ndep',varid)
  call ncmsg('ndep',ncstatus)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ncmsg('ndep',ncstatus)
  call ccmpi_distribute(duma,dumg)
  casaflux%Nmindep=duma(cmap)/365.*1.E-3
  write(6,*) "Loading N fixation rate"
  ncstatus = nf_inq_varid(ncid,'nfix',varid)
  call ncmsg('nfix',ncstatus)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ncmsg('nfix',ncstatus)
  call ccmpi_distribute(duma,dumg)
  casaflux%Nminfix=duma(cmap)/365.
  write(6,*) "Loading P dust deposition"
  ncstatus = nf_inq_varid(ncid,'pdust',varid)
  call ncmsg('pdust',ncstatus)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ncmsg('pdust',ncstatus)
  call ccmpi_distribute(duma,dumg)
  casaflux%Pdep=duma(cmap)/365.
  write(6,*) "Loading P weathering rate"
  ncstatus = nf_inq_varid(ncid,'pweather',varid)
  call ncmsg('pweather',ncstatus)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ncmsg('pweather',ncstatus)
  call ccmpi_distribute(duma,dumg)
  casaflux%Pwea=duma(cmap)/365.
  deallocate(dumg)
  ncstatus=nf_close(ncid)
  call ncmsg('CASA_readpoint',ncstatus)
else
  call ccmpi_distribute(duma)
  casamet%isorder=duma(cmap)
  call ccmpi_distribute(duma)
  casaflux%Nmindep=duma(cmap)/365.*1.E-3  
  call ccmpi_distribute(duma)
  casaflux%Nminfix=duma(cmap)/365.
  call ccmpi_distribute(duma)
  casaflux%Pdep=duma(cmap)/365.
  call ccmpi_distribute(duma)
  casaflux%Pwea=duma(cmap)/365.
end if

where (veg%iveg==12)
  ! P fertilizer =13 Mt P globally in 1994
  casaflux%Pdep=casaflux%Pdep+0.7/365.
  ! N fertilizer =86 Mt N globally in 1994
  casaflux%Nminfix=casaflux%Nminfix+4.3/365.
end where

if (any(casamet%isorder.le.0.or.casamet%isorder.gt.12)) then
  write(6,*) "ERROR: Invalid isorder in ",trim(casafile)
  stop
end if

end subroutine casa_readpoint

! *************************************************************************************  
! This subroutine reads the MODIS derived leaf phenology data
subroutine casa_readphen(fphen)

use cc_mpi

implicit none

include 'newmpar.h'
include 'mpif.h'

integer, parameter :: csirovt=17 ! number of CSIRO functional plant types
integer, parameter :: nphen=8 ! was 10(IGBP). changed by Q.Zhang @01/12/2011
integer np,nx,ilat,ierr,ivp
integer, dimension(271,csirovt) :: greenup,fall,phendoy1
integer, dimension(nphen) :: greenupx,fallx,xphendoy1
integer, dimension(nphen) :: ivtx
integer, dimension(mxvt) :: igbp2csiromap
real, dimension(271) :: xlat
character(len=*), intent(in) :: fphen

! CSIRO types (only 3 to 10 have a non-trivial leaf phenology, others are set to steady LAI)
!1  evergreen_needleleaf
!2  evergreen_broadleaf  (same leaf phenology as 1)
!3  deciduous_needleleaf
!4  deciduous_broadleaf
!5  shrub
!6  C3 grassland
!7  C4 grassland (same leaf phenology as 6)
!8  Tundra       (same leaf phenology as 6)
!9  C3 cropland
!10 C4 cropland  (same leaf phenology as 9)
!11 wetland
!12 empty
!13 empty
!14 barren
!15 urban
!16 lakes
!17 ice

! note - igbp 5 and 6 should be a mixture of csiro 1,2,3,4
!        Here we use a evergreen proxy for igbp 5 and 6
igbp2csiromap=(/ 1,2,3,4,2,2,5,5,5,8,11,9,15,13,17,12,16 /)


! initilize for evergreen PFTs
greenup = -50
fall    = 367
phendoy1= 2

if (myid==0) then
  write(6,*) "Reading CASA leaf phenology data"
  open(87,file=fphen,status='old')
  read(87,*)
  read(87,*) ivtx
  do ilat=271,1,-1
    read(87,*) xlat(ilat),greenup,fallx,xphendoy1 
    do nx=1,nphen
      greenup(ilat,ivtx) = greenupx
      fall(ilat,ivtx)    = fallx
      phendoy1(ilat,ivtx)= xphendoy1
    end do
  end do
  close(87)
end if
call MPI_Bcast(greenup,271*csirovt,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(fall,271*csirovt,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(phendoy1,271*csirovt,MPI_REAL,0,MPI_COMM_WORLD,ierr)

do np=1,mp
  ilat=(rad%latitude(np)+55.25)/0.5+1
  ilat=min(271,max(1,ilat))
  ivp=igbp2csiromap(veg%iveg(np))
  phen%phase(np)      = phendoy1(ilat,ivp)
  phen%doyphase(np,1) = greenup(ilat,ivp)          ! DOY for greenup
  phen%doyphase(np,2) = phen%doyphase(np,1)+14     ! DOY for steady LAI
  phen%doyphase(np,3) = fall(ilat,ivp)             ! DOY for leaf senescence
  phen%doyphase(np,4) = phen%doyphase(np,3)+14     ! DOY for minimal LAI season
  if (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
  if (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365
end do

return
end subroutine casa_readphen

end module cable_ccam

