! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! Subsequently modified by MJT for 5 tile mosaic and SEAESF radiation scheme
  
! - Currently all tiles have the same soil texture, but independent soil temperatures,
!   moisture, etc.
! - LAI can be interpolated between timesteps using a PWCB fit to the LAI integral
!   or LAI can be taken as constant for the month
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.  However,
!   since the input files can be modified at runtime (not all months are loaded
!   at once), then we can support off-line evolving/dynamic vegetation, etc.

! The following mappings between IGBP and CSIRO PFT were recommended by RL
    
! ivegt   IGBP type                             CSIRO PFT
! 1       Evergreen Needleleaf Forest           1.  Evergreen Needleleaf
! 2       Evergreen Broadleaf Forest            1.  Evergreen Broadleaf
! 3       Deciduous Needleaf Forest             1.  Deciduous Needleleaf
! 4       Deciduous Broadleaf Forest            1.  Deciduous Broadleaf
! 5       Mixed Forest                          1.  Deciduous Broadlead                              when -25<lat<25
!                                               0.5 Evergreen Needleleaf    0.5 Deciduous Broadlead  when lat<-25 or lat>25
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

implicit none

private
public sib4,loadcbmparm,loadtile,savetiledef,savetile,cableinflow,cbmemiss
public proglai

! The following options will eventually be moved to the globpe.f namelist
integer, save :: proglai             = 0 ! 0 prescribed LAI, 1 prognostic LAI
integer, parameter :: tracerco2      = 0 ! 0 use radiation CO2, 1 use tracer CO2 
real, parameter :: minfrac = 0.01        ! minimum non-zero tile fraction (improves load balancing)

integer, dimension(5,2), save :: pind  
integer, save :: maxnb
real, dimension(:), allocatable, save :: sv,vl1,vl2,vl3
logical, dimension(:,:), allocatable, save :: tmap
type (air_type), save            :: air
type (bgc_pool_type), save       :: bgc
type (met_type), save            :: met
type (balances_type), save       :: bal
type (radiation_type), save      :: rad
type (roughness_type), save      :: rough
type (soil_parameter_type), save :: soil
type (soil_snow_type), save      :: ssnow
type (sum_flux_type), save       :: sum_flux
type (veg_parameter_type), save  :: veg
type (canopy_type), save         :: canopy
type (casa_balance), save        :: casabal
type (casa_biome), save          :: casabiome
type (casa_flux), save           :: casaflux
type (casa_met), save            :: casamet
type (casa_pool), save           :: casapool
type (phen_variable), save       :: phen
type (physical_constants), save  :: c

contains
! ****************************************************************************

! CABLE-CCAM interface
subroutine sib4

use arrays_m
use carbpools_m
use cc_mpi
use estab
use extraout_m
use infile
use latlong_m
use liqwpar_m
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
use work2_m, only : qsttg,zo,zoh,zoq,theta,vmod,wetfac
use work3_m, only : ga
use zenith_m
  
implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'
include 'parm.h'

real fjd,r1,dlt,slag,dhr,alp,esatf
real, dimension(ifull) :: coszro2,taudar2,tmps,atmco2
real, dimension(ifull) :: tv,swdwn,alb
real(r_2), dimension(mp) :: xKNlimiting,xkleafcold,xkleafdry
real(r_2), dimension(mp) :: xkleaf,xnplimit,xNPuptake,xklitter
real(r_2), dimension(mp) :: xksoil
integer jyear,jmonth,jday,jhour,jmin
integer k,mins,nb,iq,j
integer idoy,is,ie

cansto=0.
fwet=0.
fnee=0.
fpn=0.
frd=0.
frp=0.
frpw=0.
frpr=0.
frs=0.
vlai=0.

! abort calculation if no land points on this processor  
if (mp<=0) return

! calculate zenith angle
dhr = dt/3600.
call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
fjd = float(mod(mins,525600))/1440.
call solargh(fjd,bpyear,r1,dlt,alp,slag)
call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro2,taudar2)

! calculate CO2 concentration
call setco2for(atmco2)

! set meteorological forcing
tv(:) = t(1:ifull,1)*(1.+0.61*qg(1:ifull,1)-qlg(1:ifull,1)-qfg(1:ifull,1) &
                     -qrg(1:ifull,1)-qsng(1:ifull,1)-qgrg(1:ifull,1))
! swdwn is downwelling shortwave (positive) W/m^2
albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif
albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif
alb   = swrsave*albvissav + (1.-swrsave)*albnirsav
swdwn = sgsave/(1.-alb)
do nb=1,maxnb
  is = pind(nb,1)
  ie = pind(nb,2)
  met%tk(is:ie)          =pack(theta,  tmap(:,nb))
  met%ua(is:ie)          =pack(vmod,   tmap(:,nb))
  met%ca(is:ie)          =pack(atmco2, tmap(:,nb))*1.e-6
  met%coszen(is:ie)      =pack(coszro2,tmap(:,nb))             ! use instantaneous value
  met%qv(is:ie)          =pack(qg(1:ifull,1),tmap(:,nb))       ! specific humidity in kg/kg
  met%pmb(is:ie)         =pack(ps(1:ifull),  tmap(:,nb))*0.01  ! pressure in mb at ref height
  met%precip(is:ie)      =pack(condx,  tmap(:,nb))             ! in mm not mm/sec
  met%precip_sn(is:ie)   =pack(conds+condg,  tmap(:,nb))       ! in mm not mm/sec
  met%hod(is:ie)         =pack(rlongg, tmap(:,nb))*12./pi+real(mtimer+jhour*60+jmin)/60.
  ! swrsave indicates the fraction of net VIS radiation (compared to NIR)
  ! fbeamvis indicates the beam fraction of downwelling direct radiation (compared to diffuse) for VIS
  ! fbeamnir indicates the beam fraction of downwelling direct radiation (compared to diffuse) for NIR
  met%fsd(is:ie,1)       =pack(swrsave*swdwn,        tmap(:,nb))
  met%fsd(is:ie,2)       =pack((1.-swrsave)*swdwn,   tmap(:,nb))
  rad%fbeam(is:ie,1)     =pack(fbeamvis,             tmap(:,nb))
  rad%fbeam(is:ie,2)     =pack(fbeamnir,             tmap(:,nb))
  met%fld(is:ie)         =pack(-rgsave,              tmap(:,nb))      ! long wave down (positive) W/m^2
  rough%za_tq(is:ie)     =pack(bet(1)*tv+phi_nh(:,1),tmap(:,nb))/grav ! reference height
end do
met%doy         =fjd
met%tvair       =met%tk
met%tvrad       =met%tk
met%ua          =max(met%ua,c%umin)
met%coszen      =max(met%coszen,1.e-8) 
met%hod         =mod(met%hod,24.)
rough%za_uv     =rough%za_tq
rad%fbeam(:,3)  =0.            ! dummy for now
!rough%hruff     =max(1.e-6,veg%hc-1.2*ssnow%snowd/max(ssnow%ssdnn,100.))

! Interpolate LAI.  Also need sigmf for LDR prognostic aerosols.
call setlai(sigmf,jyear,jmonth,jday,jhour,jmin)

!--------------------------------------------------------------
! CABLE
ktau_gl          = 900
kend_gl          = 999
ssnow%owetfac    = ssnow%wetfac
canopy%oldcansto = canopy%cansto
!call point2constants(C)
call ruff_resist(veg,rough,ssnow,canopy)
!met%tk=met%tk+C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
call define_air(met,air)
call init_radiation(met,rad,veg,canopy)
call surface_albedo(ssnow,veg,met,rad,soil,canopy)
call define_canopy(bal,rad,rough,air,met,dt,ssnow,soil,veg,canopy)
ssnow%otss_0     = ssnow%otss
ssnow%otss       = ssnow%tss
ssnow%owetfac    = ssnow%wetfac
call soil_snow(dt,soil,ssnow,canopy,met,bal,veg)
! adjust for new soil temperature
ssnow%deltss     = ssnow%tss - ssnow%otss
canopy%fhs       = canopy%fhs + ssnow%deltss*ssnow%dfh_dtg
!canopy%fhs_cor   = canopy%fhs_cor + ssnow%deltss*ssnow%dfh_dtg
!canopy%fes_cor   = canopy%fes_cor + ssnow%deltss*ssnow%cls*ssnow%dfe_ddq*ssnow%ddq_dtg
canopy%fh        = canopy%fhv + canopy%fhs
canopy%fev       = real(canopy%fevc + canopy%fevw)
canopy%fe        = real(canopy%fev + canopy%fes)
canopy%rnet      = canopy%fns + canopy%fnv
rad%trad         = ( (1.-rad%transd)*canopy%tv**4 + rad%transd*ssnow%tss**4 )**0.25

! EK suggestion
!canopy%cdtq =  max( 0.1*canopy%cduv, canopy%cdtq )
! MJT suggestion
canopy%cdtq =  max( 0., canopy%cdtq )

!--------------------------------------------------------------
! CASA CNP
select case (icycle)
  case(0) ! off
    call plantcarb(veg,bgc,met,canopy)
    call soilcarb(soil,ssnow,veg,bgc,met,canopy)
    call carbon_pl(dt,soil,ssnow,veg,canopy,bgc)
    canopy%fnpp = -canopy%fpn - canopy%frp
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
  case(3) ! C+N+P
    ! update casamet
    casamet%tairk = casamet%tairk + met%tk
    casamet%tsoil = casamet%tsoil + ssnow%tgg
    casamet%moist = casamet%moist + ssnow%wb
    casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dt
    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dt
    ! run CASA CNP once per day
    if (mod(ktau,nperday)==0) then
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
    canopy%frp  = real((casaflux%crmplant(:,wood)+casaflux%crmplant(:,xroot)+casaflux%crgplant(:))/86400._8)
    canopy%frs  = real(casaflux%Crsoil(:)/86400._8)
    canopy%frpw = real(casaflux%crmplant(:,wood)/86400._8)
    canopy%frpr = real(casaflux%crmplant(:,xroot)/86400._8)
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnee = real((casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/86400._8)
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
do k=1,ms
  where ( land )
    tgg(:,k)=0.
    wb(:,k)=0.
    wbice(:,k)=0.
  end where
end do
do k=1,3
  where ( land )
    tggsn(:,k)=0.
    smass(:,k)=0.
    ssdn(:,k)=0.
  end where
end do
where ( land )
  albvisdir=0.
  albvisdif=0.
  albnirdir=0.
  albnirdif=0.
  fg=0.
  eg=0.
  ga=0.
  epot=0.
  tss=0.
  zo=0.
  zoh=0.
  cduv=0.
  cdtq=0.
  ustar=0.
  wetfac=0.
  rsmin=0.
  ssdnn=0.
  snowd=0.
  snage=0.
  ! screen and 10m diagnostics - rhscrn calculated in sflux.f
  !tscrn=0.
  !uscrn=0.
  !qgscrn=0.
  !u10=0.
end where
tmps=0. ! average isflag

do nb=1,maxnb
  is = pind(nb,1)
  ie = pind(nb,2)
  ! radiation
  albvisdir=albvisdir+unpack(sv(is:ie)*rad%reffbm(is:ie,1),tmap(:,nb),0.)
  albnirdir=albnirdir+unpack(sv(is:ie)*rad%reffbm(is:ie,2),tmap(:,nb),0.)
  albvisdif=albvisdif+unpack(sv(is:ie)*rad%reffdf(is:ie,1),tmap(:,nb),0.)
  albnirdif=albnirdif+unpack(sv(is:ie)*rad%reffdf(is:ie,2),tmap(:,nb),0.)
  ! fluxes
  fg=fg+unpack(sv(is:ie)*canopy%fh(is:ie),tmap(:,nb),0.)
  eg=eg+unpack(sv(is:ie)*canopy%fe(is:ie),tmap(:,nb),0.)
  ga=ga+unpack(sv(is:ie)*canopy%ga(is:ie),tmap(:,nb),0.)
  tss=tss+unpack(sv(is:ie)*rad%trad(is:ie)**4,tmap(:,nb),0.) ! ave longwave radiation
  ! drag and mixing
  zo  =zo  +unpack(sv(is:ie)/log(zmin/rough%z0m(is:ie))**2,tmap(:,nb),0.)
  cduv=cduv+unpack(sv(is:ie)*canopy%cduv(is:ie),tmap(:,nb),0.)
  cdtq=cdtq+unpack(sv(is:ie)*canopy%cdtq(is:ie),tmap(:,nb),0.)
  ! soil
  do k=1,ms
    tgg(:,k)  =tgg(:,k)  +unpack(sv(is:ie)*ssnow%tgg(is:ie,k),        tmap(:,nb),0.)
    wb(:,k)   =wb(:,k)   +unpack(sv(is:ie)*real(ssnow%wb(is:ie,k)),   tmap(:,nb),0.)
    wbice(:,k)=wbice(:,k)+unpack(sv(is:ie)*real(ssnow%wbice(is:ie,k)),tmap(:,nb),0.)
  end do
  ! hydrology
  runoff=runoff+unpack(sv(is:ie)*ssnow%runoff(is:ie)*dt,tmap(:,nb),0.) ! convert mm/s to mm
  fwet=fwet+unpack(sv(is:ie)*canopy%fwet(is:ie),tmap(:,nb),0.)         ! used for aerosols
  wetfac=wetfac+unpack(sv(is:ie)*ssnow%wetfac(is:ie),tmap(:,nb),0.)    ! used for aerosols
  cansto=cansto+unpack(sv(is:ie)*canopy%cansto(is:ie),tmap(:,nb),0.)   ! not used
  ! diagnostic
  epot=epot+unpack(sv(is:ie)*ssnow%potev(is:ie),tmap(:,nb),0.)         ! diagnostic in history file
  vlai=vlai+unpack(sv(is:ie)*veg%vlai(is:ie),tmap(:,nb),0.)
  rsmin=rsmin+unpack(sv(is:ie)*canopy%gswx_T(is:ie),tmap(:,nb),0.)     ! diagnostic in history file
  ! carbon cycle
  fnee=fnee+unpack(sv(is:ie)*canopy%fnee(is:ie), tmap(:,nb),0.)
  fpn =fpn +unpack(sv(is:ie)*canopy%fpn(is:ie),  tmap(:,nb),0.)
  frd =frd +unpack(sv(is:ie)*canopy%frday(is:ie),tmap(:,nb),0.)
  frp =frp +unpack(sv(is:ie)*canopy%frp(is:ie),  tmap(:,nb),0.)
  frpw=frpw+unpack(sv(is:ie)*canopy%frpw(is:ie), tmap(:,nb),0.)
  frs =frs +unpack(sv(is:ie)*canopy%frs(is:ie),  tmap(:,nb),0.)
  ! snow
  tmps=tmps+unpack(sv(is:ie)*real(ssnow%isflag(is:ie)),tmap(:,nb),0.)  ! used in radiation (for nsib==3)
  do k=1,3
    tggsn(:,k)=tggsn(:,k)+unpack(sv(is:ie)*ssnow%tgg(is:ie,k),tmap(:,nb),0.)   ! for restart file
    smass(:,k)=smass(:,k)+unpack(sv(is:ie)*ssnow%smass(is:ie,k),tmap(:,nb),0.) ! for restart file
    ssdn(:,k) =ssdn(:,k) +unpack(sv(is:ie)*ssnow%ssdn(is:ie,k),tmap(:,nb),0.)  ! for restart file
  end do
  ssdnn=ssdnn+unpack(sv(is:ie)*ssnow%ssdnn(is:ie),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snage=snage+unpack(sv(is:ie)*ssnow%snage(is:ie),tmap(:,nb),0.)      ! used in radiation (for nsib==3)
  snowd=snowd+unpack(sv(is:ie)*ssnow%snowd(is:ie),tmap(:,nb),0.)
  
  !tscrn=tscrn+unpack(sv(pind(nb,1):pind(nb,2))*canopy%tscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
  !uscrn=uscrn+unpack(sv(pind(nb,1):pind(nb,2))*canopy%uscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
  !qgscrn=qgscrn+unpack(sv(pind(nb,1):pind(nb,2))*canopy%qscrn(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
end do

if (icycle==0) then
  cplant=0.
  csoil=0.
  do nb=1,maxnb
    is = pind(nb,1)
    ie = pind(nb,2)
    do k=1,ncp
      cplant(:,k)=cplant(:,k)+unpack(sv(is:ie)*bgc%cplant(is:ie,k),tmap(:,nb),0.)
    end do
    do k=1,ncs
      csoil(:,k)=csoil(:,k)+unpack(sv(is:ie)*bgc%csoil(is:ie,k),   tmap(:,nb),0.)
    end do
  end do
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
  do nb=1,maxnb
    is = pind(nb,1)
    ie = pind(nb,2)
    do k=1,mplant
      cplant(:,k) =cplant(:,k) +unpack(sv(is:ie)*real(casapool%cplant(is:ie,k)),tmap(:,nb),0.)
      niplant(:,k)=niplant(:,k)+unpack(sv(is:ie)*real(casapool%nplant(is:ie,k)),tmap(:,nb),0.)
      pplant(:,k) =pplant(:,k) +unpack(sv(is:ie)*real(casapool%pplant(is:ie,k)),tmap(:,nb),0.)
    end do
    do k=1,mlitter
      clitter(:,k) =clitter(:,k) +unpack(sv(is:ie)*real(casapool%clitter(is:ie,k)),tmap(:,nb),0.)
      nilitter(:,k)=nilitter(:,k)+unpack(sv(is:ie)*real(casapool%nlitter(is:ie,k)),tmap(:,nb),0.)
      plitter(:,k) =plitter(:,k) +unpack(sv(is:ie)*real(casapool%plitter(is:ie,k)),tmap(:,nb),0.)
    end do
    do k=1,msoil
      csoil(:,k) =csoil(:,k) +unpack(sv(is:ie)*real(casapool%csoil(is:ie,k)),tmap(:,nb),0.)
      nisoil(:,k)=nisoil(:,k)+unpack(sv(is:ie)*real(casapool%nsoil(is:ie,k)),tmap(:,nb),0.)
      psoil(:,k) =psoil(:,k) +unpack(sv(is:ie)*real(casapool%psoil(is:ie,k)),tmap(:,nb),0.)
    end do
    glai=glai+unpack(sv(is:ie)*real(casamet%glai(is:ie)),tmap(:,nb),0.)
  end do
end if

! MJT notes - ustar, cduv, fg and eg are passed to the boundary layer turbulence scheme
! zoh, zoq and zo are passed to the scrnout diagnostics routines
! rsmin is typically used by CTM

where ( land )
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
  rnet      = sgsave-rgsave-stefbo*tss**4
  !tscrn    = tscrn+273.16       ! convert from degC to degK
end where
where ( land .and. tmps>=0.5 ) ! tmps is average isflag
  isflag = 1
elsewhere
  isflag = 0
endwhere
do iq=1,ifull
  if ( land(iq) ) then
    esatf = establ(tss(iq))
    qsttg(iq) = 0.622*esatf/(ps(iq)-esatf)
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

use cc_mpi, only : myid
use radisw_m, only : rrvco2
use tracermodule, only : tractype,tracname
use tracers_m, only : tr,ngas

implicit none

include 'newmpar.h'
include 'parm.h'

integer ico2,igas
real, dimension(ifull), intent(out) :: atmco2

ico2=0
if ( tracerco2==1 ) then
  do igas=1,ngas
    if ( trim(tractype(igas))=='online' .and. trim(tracname(igas))=='cbmnep' ) then
      ico2=igas
      exit
    end if
  end do
  if ( ico2>0 ) then
    atmco2 = tr(1:ifull,1,ico2) ! use interactive tracers
  else
    atmco2 = 1.E6*rrvco2        ! from radiative CO2 forcings
  end if
else
  atmco2 = 1.E6*rrvco2          ! from radiative CO2 forcings
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
  
implicit none
  
include 'newmpar.h'
include 'parm.h'
  
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
subroutine setlai(sigmf,jyear,jmonth,jday,jhour,jmin)

use cc_mpi
  
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

select case( proglai )
  case(0) ! PWCB interpolated LAI
    imonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    if ( leap==1 ) then
      if ( mod(jyear,4)  ==0 ) imonth(2)=29
      if ( mod(jyear,100)==0 ) imonth(2)=28
      if ( mod(jyear,400)==0 ) imonth(2)=29
    end if
    monthstart = 1440*(jday-1) + 60*jhour + jmin ! mins from start month
    x = min(max(real(mtimer+monthstart)/real(1440*imonth(jmonth)),0.),1.)
    veg%vlai = vl1+vl2*x+vl3*x*x     ! LAI as a function of time
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.01 )
    elsewhere
      veg%vlai = 1.E-8
    end where

  case(1) ! prognostic LAI
    if (icycle==0) then
      write(6,*) "ERROR: CASA CNP LAI is not operational"
      call ccmpi_abort(-1)
    end if
    veg%vlai(:)=real(casamet%glai(:))

  case default
    write(6,*) "ERROR: Unknown proglai option ",proglai
    call ccmpi_abort(-1)
end select

sigmf(:)=0.
do nb=1,maxnb
  sigmf=sigmf+unpack(sv(pind(nb,1):pind(nb,2))*(1.-exp(-vextkn*veg%vlai(pind(nb,1):pind(nb,2)))),tmap(:,nb),0.)
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

integer, dimension(ifull,5) :: ivs
integer iq,n,k,ipos,iv,ncount
integer, dimension(1) :: pos
integer jyear,jmonth,jday,jhour,jmin,mins
integer, dimension(1) :: lndtst,lndtst_g
real fc3,fc4,ftu,fg3,fg4,clat,nsum
real fjd,xp
real, dimension(ifull,mxvt,0:2) :: newlai
real, dimension(mxvt,ms) :: froot2
real, dimension(mxvt,ncp) :: tcplant
real, dimension(mxvt,ncs) :: tcsoil
real, dimension(mxvt,mplant) :: ratiocnplant
real, dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real, dimension(mxvt,2) :: taul,refl  
real, dimension(ifull,mxvt) :: newgrid
real, dimension(ifull,5) :: svs,vlin,vlinprev,vlinnext
real, dimension(ifull,2) :: albsoilsn
real, dimension(12,msoil) :: rationpsoil
real, dimension(ncp) :: ratecp
real, dimension(ncs) :: ratecs
real, dimension(mxvt) :: canst1,leaf_w,leaf_l,ejmax,hc,rp20
real, dimension(mxvt) :: rpcoef,shelrb,vcmax,xfang
real, dimension(mxvt) :: tminvj,tmaxvj,vbeta
real, dimension(mxvt) :: extkn,rootbeta,vegcf,c4frac
real, dimension(mxvt) :: leafage,woodage,frootage,metage
real, dimension(mxvt) :: strage,cwdage,micage,slowage,passage
real, dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real, dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real, dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real, dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real, dimension(mxvt) :: cleaf,cwood,cfroot,cmet,cstr,ccwd,cmic,cslow,cpass,nleaf
real, dimension(mxvt) :: nwood,nfroot,nmet,nstr,ncwd,nmic,nslow,npass,xpleaf,xpwood
real, dimension(mxvt) :: xpfroot,xpmet,xpstr,xpcwd,xpmic,xpslow,xppass,clabileage
real, dimension(ifull) :: albsoil, savannafrac
real, dimension(12) :: xkmlabp,xpsorbmax,xfPleach
character(len=*), intent(in) :: fveg,fvegprev,fvegnext,fphen,casafile

if ( myid==0 ) write(6,*) "Initialising CABLE"

if ( cbm_ms/=ms ) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  call ccmpi_abort(-1)
end if

! redefine rhos
rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

! biophysical parameter tables
hc    =(/   17.,  35.,  15.5,  20.,   0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567,  0.2, 6.017,  0.2,  0.2,  0.2,  0.2 /)
xfang =(/  0.01,  0.1,  0.01, 0.25,  0.01,  -0.3,  -0.3,  -0.3, -0.3, -0.3,  -0.3,  0.1,    0.,   0.,   0.,   0.,   0. /)
leaf_w=(/ 0.001, 0.05, 0.001, 0.08, 0.005,  0.01,  0.01,  0.01, 0.01, 0.01,  0.01, 0.03, 0.015, 0.00,   0.,   0.,   0. /)
leaf_l=(/ 0.055, 0.10, 0.040, 0.15, 0.100,  0.30,  0.30,  0.30, 0.30, 0.30,  0.30, 0.30, 0.242, 0.03, 0.03, 0.03, 0.03 /)
canst1=0.1
shelrb=2.
extkn=0.001 ! new definition for nitrogen (since CABLE v1.9b)
!refl(:,1)=(/ 0.043,0.053,0.050,0.064,0.120,0.099,0.080,0.058,0.080,0.050,0.060,0.031,0.051,0.175,0.079,0.079,0.159 /)
!refl(:,2)=(/ 0.211,0.245,0.242,0.266,0.480,0.423,0.320,0.171,0.320,0.200,0.191,0.105,0.172,0.336,0.153,0.153,0.305 /)
!taul(:,1)=(/ 0.035,0.035,0.040,0.035,0.060,0.063,0.080,0.040,0.080,0.050,0.041,0.013,0.033,0.029,0.013,0.013,0.026 /)
!taul(:,2)=(/ 0.120,0.250,0.120,0.250,0.192,0.200,0.165,0.124,0.165,0.125,0.081,0.110,0.181,0.126,0.063,0.063,0.063 /)
refl(:,1)=(/ 0.062,0.076,0.056,0.092,0.100,0.110,0.100,0.117,0.100,0.090,0.108,0.055,0.091,0.238,0.143,0.143,0.159 /)
refl(:,2)=(/ 0.302,0.350,0.275,0.380,0.400,0.470,0.400,0.343,0.400,0.360,0.343,0.190,0.310,0.457,0.275,0.275,0.305 /)
taul(:,1)=(/ 0.050,0.050,0.045,0.050,0.050,0.070,0.100,0.080,0.100,0.090,0.075,0.023,0.059,0.039,0.023,0.023,0.026 /)
taul(:,2)=(/ 0.100,0.250,0.144,0.250,0.240,0.250,0.150,0.124,0.150,0.225,0.146,0.198,0.163,0.189,0.113,0.113,0.113 /)
vegcf    =(/    9.,  14.,   9.,   8.,   5.,   7.,   7.,   5.,   7.,   1.,   7.,   1.,   1.,   1.,   1.,   1.,   1. /)
vcmax=(/ 40.E-6,55.E-6,40.E-6,60.E-6,40.E-6,60.E-6,10.E-6,40.E-6,80.E-6,80.E-6,60.E-6,17.E-6,1.E-6,17.E-6,17.E-6,17.E-6,17.E-6 /)
ejmax=2.*vcmax
rp20=(/ 3., 0.6, 3., 2.2, 1., 1.5, 2.8, 2.5, 1.5, 1., 1.5, 1., 1., 1., 1., 1., 1. /)
rpcoef=0.0832
rs20  =(/   1.,   1.,  1.,  1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   0.,   0.,   0. /)
tminvj=(/ -15., -15.,  5.,  5., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15. /)
tmaxvj=(/ -10., -10., 10., 15., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10. /)
vbeta=1.
rootbeta=(/ 0.943,0.962,0.966,0.961,0.964,0.943,0.943,0.943,0.961,0.961,0.943,0.975,0.961,0.961,0.961,0.961,0.961 /)
tcplant(:,1)=(/ 200.  , 300.  , 200. , 300.  , 159. , 250., 250., 250., 150., 150., 250., 1., 0.1, 0., 1., 1., 0. /)
tcplant(:,2)=(/ 10217., 16833., 5967., 12000., 5000., 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0., 0.,  0., 0., 0., 0. /)
tcplant(:,3)=(/ 876.  , 1443. , 511. , 1029. , 500. , 500., 500., 500., 607., 607., 500., 1., 0.1, 0., 1., 1., 0. /)
tcsoil(:,1) =(/ 184.  , 303.  , 107. , 216.  , 100. , 275., 275., 275., 149., 149., 275., 1., 0.1, 1., 1., 1., 1. /)
tcsoil(:,2) =(/ 367.  , 606.  , 214. , 432.  , 250. , 314., 314., 314., 300., 300., 314., 1., 0.1, 1., 1., 1., 1. /)
ratecp(1:3)=(/ 1., 0.03, 0.14 /)
ratecs(1:2)=(/ 2., 0.5 /)
c4frac=(/ 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0. /)

! read CABLE biome and LAI data
if ( myid==0 ) then
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
cplant=0.
csoil=0.
pind=ifull+1
mvtype=mxvt
mstype=mxst

if ( myid==0 ) write(6,*) "Mapping IGBP classes to CSIRO PFTs"
mp=0
newgrid=0.
newlai=0.
savannafrac=0.
do iq=1,ifull
  if ( land(iq) ) then
    clat=rlatt(iq)*180./pi
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
    do n=1,5
      select case (ivs(iq,n))
        case (1,2,3,4,11)
          newgrid(iq,ivs(iq,n))=newgrid(iq,ivs(iq,n))+svs(iq,n)
          newlai(iq,ivs(iq,n),0)=newlai(iq,ivs(iq,n),0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,ivs(iq,n),1)=newlai(iq,ivs(iq,n),1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,ivs(iq,n),2)=newlai(iq,ivs(iq,n),2)+svs(iq,n)*vlinnext(iq,n)
        case (5)
          if (abs(clat)>25.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.5
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.5*vlinnext(iq,n)
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)*0.5
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*0.5*vlinnext(iq,n)
          else if (abs(clat)>24.5) then
            xp=abs(clat)-24.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.5*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.5*vlinprev(iq,n)*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.5*vlin(iq,n)*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.5*vlinnext(iq,n)*xp
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)*(1.-0.5*xp)
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*vlinprev(iq,n)*(1.-0.5*xp)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*vlin(iq,n)*(1.-0.5*xp)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*vlinnext(iq,n)*(1.-0.5*xp)
          else
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*vlinprev(iq,n)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*vlin(iq,n)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*vlinnext(iq,n)
          end if
        case (6)
          newgrid(iq,5)=newgrid(iq,5)+svs(iq,n)*0.8
          newlai(iq,5,0)=newlai(iq,5,0)+svs(iq,n)*0.8*vlinprev(iq,n)
          newlai(iq,5,1)=newlai(iq,5,1)+svs(iq,n)*0.8*vlin(iq,n)
          newlai(iq,5,2)=newlai(iq,5,2)+svs(iq,n)*0.8*vlinnext(iq,n)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.2*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.2*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.2*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.2*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.2*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.2*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.2*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.2*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.2*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.2*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.2*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.2*ftu*vlinnext(iq,n)
        case (7)
          newgrid(iq,5)=newgrid(iq,5)+svs(iq,n)*0.2
          newlai(iq,5,0)=newlai(iq,5,0)+svs(iq,n)*0.2*vlinprev(iq,n)
          newlai(iq,5,1)=newlai(iq,5,1)+svs(iq,n)*0.2*vlin(iq,n)
          newlai(iq,5,2)=newlai(iq,5,2)+svs(iq,n)*0.2*vlinnext(iq,n)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.8*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.8*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.8*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.8*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.8*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.8*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.8*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.8*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.8*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.8*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.8*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.8*ftu*vlinnext(iq,n)
        case (8)
          if (abs(clat)>40.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.4
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.4*vlinnext(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.4*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*vlinprev(iq,n)*0.4*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*vlin(iq,n)*0.4*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*vlinnext(iq,n)*0.4*xp   
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.4*(1.-xp)
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.4*(1.-xp)
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*vlinprev(iq,n)*0.4*(1.-xp)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*vlin(iq,n)*0.4*(1.-xp)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*vlinnext(iq,n)*0.4*(1.-xp)
          else
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.4
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.4
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*0.4*vlinnext(iq,n)
          end if
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.6*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.6*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.6*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.6*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.6*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.6*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.6*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.6*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.6*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.6*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.6*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.6*ftu*vlinnext(iq,n)
        case (9)
          if (abs(clat)>40.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.1
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.1*vlinnext(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.1*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*vlinprev(iq,n)*0.1*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*vlin(iq,n)*0.1*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*vlinnext(iq,n)*0.1*xp
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.1*(1.-xp)
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.1*(1.-xp)
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*vlinprev(iq,n)*0.1*(1.-xp)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*vlin(iq,n)*0.1*(1.-xp)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*vlinnext(iq,n)*0.1*(1.-xp)
          else
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.1
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.1
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*0.1*vlinnext(iq,n)
          end if
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.9*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.9*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.9*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.9*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.9*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.9*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.9*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.9*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.9*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.9*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.9*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.9*ftu*vlinnext(iq,n)
        case (10)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*ftu*vlinnext(iq,n)
        case (12,14)
          newgrid(iq,9)=newgrid(iq,9)+svs(iq,n)*fc3
          newlai(iq,9,0)=newlai(iq,9,0)+svs(iq,n)*fc3*vlinprev(iq,n)
          newlai(iq,9,1)=newlai(iq,9,1)+svs(iq,n)*fc3*vlin(iq,n)
          newlai(iq,9,2)=newlai(iq,9,2)+svs(iq,n)*fc3*vlinnext(iq,n)
          newgrid(iq,10)=newgrid(iq,10)+svs(iq,n)*fc4
          newlai(iq,10,0)=newlai(iq,10,0)+svs(iq,n)*fc4*vlinprev(iq,n)
          newlai(iq,10,1)=newlai(iq,10,1)+svs(iq,n)*fc4*vlin(iq,n)
          newlai(iq,10,2)=newlai(iq,10,2)+svs(iq,n)*fc4*vlinnext(iq,n)
        case (13)
          newgrid(iq,15)=newgrid(iq,15)+svs(iq,n)
          newlai(iq,15,0)=newlai(iq,15,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,15,1)=newlai(iq,15,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,15,2)=newlai(iq,15,2)+svs(iq,n)*vlinnext(iq,n)
        case (15)
          newgrid(iq,17)=newgrid(iq,17)+svs(iq,n)
          newlai(iq,17,0)=newlai(iq,17,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,17,1)=newlai(iq,17,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,17,2)=newlai(iq,17,2)+svs(iq,n)*vlinnext(iq,n)
        case (16)
          newgrid(iq,14)=newgrid(iq,14)+svs(iq,n)
          newlai(iq,14,0)=newlai(iq,14,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,14,1)=newlai(iq,14,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,14,2)=newlai(iq,14,2)+svs(iq,n)*vlinnext(iq,n)
        case (17)
          newgrid(iq,16)=newgrid(iq,16)+svs(iq,n)
          newlai(iq,16,0)=newlai(iq,16,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,16,1)=newlai(iq,16,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,16,2)=newlai(iq,16,2)+svs(iq,n)*vlinnext(iq,n)
        case DEFAULT
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq,ivs,land=",myid,iq,ivs(iq,n),land(iq)
          call ccmpi_abort(-1)
      end select
    end do
    if (newgrid(iq,2)>0.) then
      savannafrac(iq)=savannafrac(iq)/newgrid(iq,2)
    end if
    where (newgrid(iq,:)>0.)
      newlai(iq,:,0)=newlai(iq,:,0)/newgrid(iq,:)
      newlai(iq,:,1)=newlai(iq,:,1)/newgrid(iq,:)
      newlai(iq,:,2)=newlai(iq,:,2)/newgrid(iq,:)
    end where
    ipos=count(newgrid(iq,:)>0.)
    do while (ipos>5)
      pos=minloc(newgrid(iq,:),newgrid(iq,:)>0.)
      newgrid(iq,pos(1))=0.
      nsum=sum(newgrid(iq,:))
      newgrid(iq,:)=newgrid(iq,:)/nsum
      ipos=count(newgrid(iq,:)>0.)
    end do    
    do while (any(newgrid(iq,:)<minfrac.and.newgrid(iq,:)>0.))
      pos=minloc(newgrid(iq,:),newgrid(iq,:)>0.)
      newgrid(iq,pos(1))=0.
      nsum=sum(newgrid(iq,:))
      newgrid(iq,:)=newgrid(iq,:)/nsum
    end do
    ipos=count(newgrid(iq,:)>0.)
    !if (ipos>5) then
    !  write(6,*) "ERROR: Too many CABLE tiles"
    !  call ccmpi_abort(-1)
    !end if
    !if (ipos<1) then
    !  write(6,*) "ERROR: Missing CABLE tiles"
    !  call ccmpi_abort(-1)
    !end if    
    mp=mp+ipos
  end if
end do

if (nmaxpr==1) then
  write(6,*) "myid,landtile ",myid,mp

  lndtst(1)=0
  if (mp>0) lndtst(1)=1
  call ccmpi_reduce(lndtst(1:1),lndtst_g(1:1),"sum",0,comm_world)
  if (myid==0) then
    write(6,*) "Processors with land ",lndtst_g(1),nproc
  end if
end if

! if CABLE is present on this processor, then start allocating arrays
! Write messages here in case myid==0 has no land-points (mp==0)
if (myid==0) then
  write(6,*) "Allocating CABLE and CASA CNP arrays"
  if (icycle==0) then
    write(6,*) "Using CABLE carbon cycle"
  else
    write(6,*) "Using CASA CNP"
  end if
end if

if (mp>0) then
  
  allocate(sv(mp))
  allocate(vl1(mp),vl2(mp),vl3(mp))
  allocate(tmap(ifull,5))
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

  ! Cable configuration
  cable_user%ssnow_POTEV = ""
  knode_gl = myid
  kwidth_gl = nint(dt) ! MJT notes - what happens when the timestep is less than a second?
  if (kwidth_gl == 0) then
    write(6,*) "ERROR: Timestep too small for CABLE"
    call ccmpi_abort(-1)
  end if
  
  ! soil parameters
  soil%zse        = zse ! soil layer thickness
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
  
  ! Eva's method for ACCESS1.3
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
  tmap=.false.

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  ! now up to 5 PFT tiles from 5 IGBP classes (need correct order for vectorisation)
  ipos=0
  do n=1,5
    pind(n,1)=ipos+1
    do iq=1,ifull
      if (land(iq)) then
        ncount=0
        do iv=1,mxvt
          if (newgrid(iq,iv)>0.) then
            ncount=ncount+1
            if (ncount==n) exit
          end if
        end do
        if (ncount==n) then
          ipos=ipos+1
          tmap(iq,n)=.true.
          sv(ipos)=newgrid(iq,iv)
          veg%iveg(ipos)=iv
          soil%isoilm(ipos)=isoilm(iq)
          newlai(iq,iv,:)=max(newlai(iq,iv,:),0.01)
          if (fvegprev/=' '.and.fvegnext/=' ') then
            newlai(iq,iv,1)=newlai(iq,iv,1)+newlai(iq,iv,0)
            newlai(iq,iv,2)=newlai(iq,iv,2)+newlai(iq,iv,1)
            vl1(ipos)=0.5*newlai(iq,iv,1)
            vl2(ipos)=4.*newlai(iq,iv,1)-5.*newlai(iq,iv,0)-newlai(iq,iv,2)
            vl3(ipos)=1.5*(newlai(iq,iv,2)+3.*newlai(iq,iv,0)-3.*newlai(iq,iv,1))
          else
            vl1(ipos)=newlai(iq,iv,1)
            vl2(ipos)=0.
            vl3(ipos)=0.
          end if
          if (veg%iveg(ipos)>=14.and.veg%iveg(ipos)<=17) then
            vl1(ipos)=1.E-8
            vl2(ipos)=0.
            vl3(ipos)=0.
          end if
        end if
      end if
    end do
    pind(n,2)=ipos
  end do
  
  if (ipos/=mp) then
    write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
    call ccmpi_abort(-1)
  end if

  ! Load CABLE arrays
  ivegt=ivs(:,1) ! diagnostic (usually IGBP, not CSIRO pft)
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
  veg%dleaf     = sqrt(max(leaf_w(veg%iveg)*leaf_l(veg%iveg),1.e-20))
  veg%xalbnir   = 1. ! not used
  veg%taul(:,1) = taul(veg%iveg,1)
  veg%taul(:,2) = taul(veg%iveg,2)  
  veg%refl(:,1) = refl(veg%iveg,1)
  veg%refl(:,2) = refl(veg%iveg,2)  
  veg%extkn     = extkn(veg%iveg)
  veg%rs20      = rs20(veg%iveg)
  veg%vegcf     = vegcf(veg%iveg)
  veg%frac4     = c4frac(veg%iveg)
  do k=1,ms
    veg%froot(:,k)=froot2(veg%iveg,k)
  end do

  ! calculate max tile number
  do n=1,5
    if (pind(n,1)<=mp) then
      maxnb=n
    end if
  end do
  
  ! MJT special case for (woody) savannas
  do n=1,maxnb
    where (veg%iveg(pind(n,1):pind(n,2))==2)
      veg%hc(pind(n,1):pind(n,2))=veg%hc(pind(n,1):pind(n,2))+pack((17.-hc(2))*savannafrac,tmap(:,n))
    end where
  end do
  
  ! Calculate LAI and veg fraction diagnostics
  call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jyear,jmonth,jday,jhour,jmin)
  vlai=0.
  do n=1,maxnb
    vlai=vlai+unpack(sv(pind(n,1):pind(n,2))*veg%vlai(pind(n,1):pind(n,2)),tmap(:,n),0.)
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
  do n=1,maxnb
    soil%albsoil(pind(n,1):pind(n,2),1)=pack(albvisnir(:,1),tmap(:,n))
    soil%albsoil(pind(n,1):pind(n,2),2)=pack(albvisnir(:,2),tmap(:,n))
  end do
  soil%albsoil(:,3)=0.05
    
  where (land)
    albsoil(:)=0.5*sum(albvisnir,2)
  end where
  where (albsoil<=0.14.and.land)
    !sfact=0.5 for alb <= 0.14
    albsoilsn(:,1)=(1.00/1.50)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.50)*albsoil(:)
  elsewhere ((albsoil(:)<=0.2).and.land)
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

  do n=1,maxnb
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
  canopy%ga=0.
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
  
  
  if (icycle==0) then
    ! Initialise CABLE carbon pools
    cplant=0.
    csoil=0.
    do n=1,maxnb
      do k=1,ncp
        cplant(:,k)=cplant(:,k)+unpack(sv(pind(n,1):pind(n,2))*tcplant(veg%iveg(pind(n,1):pind(n,2)),k),tmap(:,n),0.)
      end do
      do k=1,ncs
        csoil(:,k) =csoil(:,k) +unpack(sv(pind(n,1):pind(n,2))*tcsoil(veg%iveg(pind(n,1):pind(n,2)),k), tmap(:,n),0.)
      end do
    end do
  else
    ! CASA CNP
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    call alloc_phenvariable(phen,mp)
    
    casamet%lat=rad%latitude
    
    call casa_readpoint(casafile) ! read point sources

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

    cplant=0.
    clitter=0.
    csoil=0.
    niplant=0.
    nilitter=0.
    nisoil=0.
    pplant=0.
    plitter=0.
    psoil=0.
    glai=0.
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
      glai(:)=glai(:)+unpack(sv(pind(n,1):pind(n,2))*real(casamet%glai(pind(n,1):pind(n,2))),tmap(:,n),0.)
    end do

  end if ! icycle>0

end if
  
if (myid==0) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine loadcbmparm

! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi
use infile

implicit none
  
include 'newmpar.h'
include 'darcdf.h'
include 'parmgeom.h'  ! rlong0,rlat0,schmidt  
  
character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer, dimension(ifull,5), intent(out) :: ivs
integer, dimension(ifull_g,5) :: ivsg  
integer, dimension(3) :: spos,npos
integer n,iq,ilx,jlx,iad 
integer ncidx,iernc,varid,ndims
real, dimension(ifull,5), intent(out) :: svs,vlinprev,vlin,vlinnext
real, dimension(ifull_g,5) :: svsg,vling
real rlong0x,rlat0x,schmidtx,dsx,ra,rb,cablever
character(len=47) header  
character(len=6) vname
real, parameter :: cableversion = 223. ! version id for input data

write(6,*) "Reading land-use parameters for CABLE"
if (lncveg == 1) then
  ! assume this file grid has been tested when opened
  spos(1:3)=1
  npos(1)=il_g
  npos(2)=6*il_g
  npos(3)=1
  call ccnf_inq_dimlen(ncidveg,'longitude',ilx)
  call ccnf_inq_dimlen(ncidveg,'latitude',jlx)
  call ccnf_get_attg(ncidveg,'lon0',rlong0x)
  call ccnf_get_attg(ncidveg,'lat0',rlat0x)
  call ccnf_get_attg(ncidveg,'schmidt',schmidtx)
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  call ccnf_get_attg(ncidveg,'cableversion',cablever,ierr=iernc)
  if (iernc /= 0) then
    write(6,*) "Missing version of CABLE data"
    write(6,*) "Regenerate land-use data with up-to-date version of igbpveg"
    call ccmpi_abort(-1)
  end if
  if (cablever /= cableversion) then
    write(6,*) "Wrong version of CABLE data"
    write(6,*) "Expecting ",cableversion
    write(6,*) "Found     ",cablever
    call ccmpi_abort(-1)
  end if
  do n=1,5
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
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if (fvegprev/=' '.and.fvegnext/=' ') then
    call ccnf_open(fvegprev,ncidx,iernc)
    if (iernc/=0) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do n=1,5
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
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    do n=1,5
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinnext,vling)
  else
    vlinprev=-1.
    vlinnext=-1.    
  end if

else
  open(87,file=fveg,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  do iq=1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if (fvegprev/=' '.and.fvegnext/=' ') then
    open(87,file=fvegprev,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
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
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
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
end if
return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi
  
implicit none

include 'newmpar.h'

character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer, dimension(ifull,5), intent(out) :: ivs
real, dimension(ifull,5), intent(out) :: svs, vlinprev, vlin, vlinnext

call ccmpi_distribute(ivs)
call ccmpi_distribute(svs)
call ccmpi_distribute(vlin)
if (fvegprev/=' '.and.fvegnext/=' ') then
  call ccmpi_distribute(vlinprev)
  call ccmpi_distribute(vlinnext)
else
  vlinprev=-1.
  vlinnext=-1.
end if    
  
return
end subroutine vegtb

! *************************************************************************************  
! This subroutine loads CABLE tile data
subroutine loadtile

use carbpools_m
use cc_mpi
use infile
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none

include 'newmpar.h'
include 'darcdf.h'
include 'parm.h'  
  
integer k, n, ierr, idv
integer, dimension(1) :: dum
real, dimension(ifull) :: dat
real totdepth
logical tst
character(len=11) vname

! check that CABLE data exists in restart file
! and communicate the result to all processors
! as not all processors are assigned an input file
ierr = 1
if ( io_in == 1 ) then
  if ( myid==0 .or. pfall ) then
    call ccnf_inq_varid(ncid,"tgg1_5",idv,tst)
    if ( tst ) then
      ierr = 1
    else
      ierr = 0
    end if
  end if
  if ( .not.pfall ) then
    dum(1) = ierr
    call ccmpi_bcast(dum(1:1),0,comm_world)
    ierr = dum(1)
  end if
end if
  
! Cannot locate tile data, use diagnostic data instead
if ( ierr /= 0 ) then
  if ( myid == 0 ) write(6,*) "Use gridbox averaged data to initialise CABLE"
  if ( mp > 0 ) then
    do n = 1,maxnb
      do k = 1,ms
        ssnow%tgg(pind(n,1):pind(n,2),k)   = pack(tgg(:,k),  tmap(:,n))
        ssnow%wb(pind(n,1):pind(n,2),k)    = pack(wb(:,k),   tmap(:,n))
        ssnow%wbice(pind(n,1):pind(n,2),k) = pack(wbice(:,k),tmap(:,n))
      end do
      do k = 1,3
        ssnow%tggsn(pind(n,1):pind(n,2),k)  = pack(tggsn(:,k),tmap(:,n))
        ssnow%smass(pind(n,1):pind(n,2),k)  = pack(smass(:,k),tmap(:,n))
        ssnow%ssdn(pind(n,1):pind(n,2),k)   = pack(ssdn(:,k), tmap(:,n))
        ssnow%sdepth(pind(n,1):pind(n,2),k) = pack(snowd/3.,  tmap(:,n))
        ssnow%sconds(pind(n,1):pind(n,2),k) = 0.2
      end do      
      ssnow%ssdnn(pind(n,1):pind(n,2))  = pack(ssdnn, tmap(:,n))
      ssnow%isflag(pind(n,1):pind(n,2)) = pack(isflag,tmap(:,n))
      ssnow%snowd(pind(n,1):pind(n,2))  = pack(snowd, tmap(:,n))
      ssnow%snage(pind(n,1):pind(n,2))  = pack(snage, tmap(:,n))
    end do
    ssnow%rtsoil=50.
    canopy%cansto=0.
    canopy%us=0.01
    ssnow%pudsto=0.
    ssnow%wetfac=0.
    ssnow%osnowd=ssnow%snowd
    if ( icycle == 0 ) then
      do n = 1,maxnb
        do k = 1,ncp
          bgc%cplant(pind(n,1):pind(n,2),k) = pack(cplant(:,k),tmap(:,n))
        end do
        do k = 1,ncs
          bgc%csoil(pind(n,1):pind(n,2),k) = pack(csoil(:,k),tmap(:,n))
        end do
      end do
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
        casamet%glai(pind(n,1):pind(n,2)) = pack(glai,tmap(:,n))
      end do
    end if
  end if
else
  ! Located CABLE tile data
  if ( myid == 0 ) write(6,*) "Use tiled data to initialise CABLE"
  do n = 1,5
    do k = 1,ms
      write(vname,'("tgg",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%tgg(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("wb",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%wb(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("wbice",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%wbice(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
    end do
    do k = 1,3
      write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%tggsn(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("smass",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%smass(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%ssdn(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("sdepth",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%sdepth(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      write(vname,'("sconds",I1.1,"_",I1.1)') k,n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) ssnow%sconds(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
    end do
    write(vname,'("ssdnn_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%ssdnn(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("sflag_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%isflag(pind(n,1):pind(n,2))=nint(pack(dat,tmap(:,n)))
    write(vname,'("snd_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%snowd(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("osnd_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%osnowd(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("snage_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%snage(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("rtsoil_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%rtsoil(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("cansto_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) canopy%cansto(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("us_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) canopy%us(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("pudsto_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%pudsto(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    write(vname,'("wetfac_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) ssnow%wetfac(pind(n,1):pind(n,2))=pack(dat,tmap(:,n))
    if ( icycle == 0 ) then
      do k = 1,ncp
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) bgc%cplant(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      enddo
      do k = 1,ncs
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) bgc%csoil(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      enddo
    else
      do k = 1,mplant
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%cplant(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("nplant",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%nplant(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("pplant",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%pplant(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      end do
      do k = 1,mlitter
        write(vname,'("clitter",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%clitter(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("nlitter",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%nlitter(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("plitter",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%plitter(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      end do
      do k = 1,msoil
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%csoil(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("nsoil",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%nsoil(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
        write(vname,'("psoil",I1.1,"_",I1.1)') k,n
        call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
        if ( n <= maxnb ) casapool%psoil(pind(n,1):pind(n,2),k) = pack(dat,tmap(:,n))
      end do
      write(vname,'("glai_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casamet%glai(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("phenphase_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) phen%phase(pind(n,1):pind(n,2)) = nint(pack(dat,tmap(:,n)))
      write(vname,'("clabile_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casapool%clabile(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("nsoilmin_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casapool%nsoilmin(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("psoillab_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casapool%psoillab(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("psoilsorb_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casapool%psoilsorb(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
      write(vname,'("psoilocc_",I1.1)') n
      call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
      if ( n <= maxnb ) casapool%psoilocc(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
    end if
    ! CABLE correction terms
    write(vname,'("fhscor_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) canopy%fhs_cor(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
    write(vname,'("fescor_",I1.1)') n
    call histrd1(iarchi-1,ierr,vname,il_g,dat,ifull)
    if ( n <= maxnb ) canopy%fes_cor(pind(n,1):pind(n,2)) = pack(dat,tmap(:,n))
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
end if
  
! Some fixes for rounding errors
if ( mp > 0 ) then

  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + soil%zse(k)*100.
  enddo

  ssnow%wb = max(ssnow%wb,0._r_2)
  ssnow%wbice = max(ssnow%wbice,0._r_2)
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
    bgc%cplant = max(bgc%cplant,0.)
    bgc%csoil = max(bgc%csoil,0.)
  else
    casapool%cplant     = max(0._r_2,casapool%cplant)
    casapool%clitter    = max(0._r_2,casapool%clitter)
    casapool%csoil      = max(0._r_2,casapool%csoil)
    casabal%cplantlast(1:mp,1:mplant)   = casapool%cplant(1:mp,1:mplant)
    casabal%clitterlast(1:mp,1:mlitter) = casapool%clitter(1:mp,1:mlitter)
    casabal%csoillast(1:mp,1:msoil)     = casapool%csoil(1:mp,1:msoil)
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal     = 0.
    casabal%FCgppyear   = 0.
    casabal%FCrpyear    = 0.
    casabal%FCnppyear   = 0.
    casabal%FCrsyear    = 0.
    casabal%FCneeyear   = 0.
    casapool%nplant     = max(1.e-6_r_2,casapool%nplant)
    casapool%nlitter    = max(1.e-6_r_2,casapool%nlitter)
    casapool%nsoil      = max(1.e-6_r_2,casapool%nsoil)
    casapool%nsoilmin   = max(1.e-6_r_2,casapool%nsoilmin)
    casabal%nplantlast(1:mp,1:mplant)   = casapool%nplant(1:mp,1:mplant)
    casabal%nlitterlast(1:mp,1:mlitter) = casapool%nlitter(1:mp,1:mlitter)
    casabal%nsoillast(1:mp,1:msoil)     = casapool%nsoil(1:mp,1:msoil)      
    casabal%nsoilminlast= casapool%nsoilmin
    casabal%sumnbal     = 0.
    casabal%FNdepyear   = 0.
    casabal%FNfixyear   = 0.
    casabal%FNsnetyear  = 0.
    casabal%FNupyear    = 0.
    casabal%FNleachyear = 0.
    casabal%FNlossyear  = 0.
    casapool%pplant     = max(1.0e-7_r_2,casapool%pplant)
    casapool%plitter    = max(1.0e-7_r_2,casapool%plitter)
    casapool%psoil      = max(1.0e-7_r_2,casapool%psoil)
    casapool%Psoillab   = max(1.0e-7_r_2,casapool%psoillab)
    casapool%psoilsorb  = max(1.0e-7_r_2,casapool%psoilsorb)
    casapool%psoilocc   = max(1.0e-7_r_2,casapool%psoilocc)
    casabal%pplantlast(1:mp,1:mplant)   = casapool%pplant(1:mp,1:mplant)
    casabal%plitterlast(1:mp,1:mlitter) = casapool%plitter(1:mp,1:mlitter)
    casabal%psoillast(1:mp,1:msoil)     = casapool%psoil(1:mp,1:msoil)       
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
subroutine savetiledef(idnc,local,idim)

use carbpools_m
use cc_mpi, only : myid
use infile
  
implicit none

include 'newmpar.h'
  
integer, intent(in) :: idnc
integer k,n
integer, dimension(3), intent(in) :: idim  
character(len=11) vname
character(len=40) lname
logical, intent(in) :: local
  
if (myid==0.or.local) then
  if (myid==0) then
    write(6,*) "Defining CABLE tile data"
  end if
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
      write(lname,'("Snow depth ",I1.1," tile ",I1.1)') k,n
      write(vname,'("sdepth",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'mm',0.,6500.,0,-1)
      write(lname,'("Snow sconds ",I1.1," tile ",I1.1)') k,n
      write(vname,'("sconds",I1.1,"_",I1.1)') k,n 
      call attrib(idnc,idim,3,vname,lname,'none',0.,6.5,0,-1)
    end do
    write(lname,'("Snow ssdnn tile ",I1.1)') n
    write(vname,'("ssdnn_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'kg/m3',0.,650.,0,-1)
    write(lname,'("Snow flag tile ",I1.1)') n
    write(vname,'("sflag_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'mm',0.,6.5,0,-1)
    write(lname,'("Snow depth tile ",I1.1)') n
    write(vname,'("snd_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'mm',0.,6500.,0,-1)  ! -1=long
    write(lname,'("Old snow depth tile ",I1.1)') n
    write(vname,'("osnd_",I1.1)') n
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
    write(lname,'("us tile ",I1.1)') n
    write(vname,'("us_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'m/s',0.,13.,0,-1)    
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
    write(lname,'("Sensible correction term ",I1.1)') n
    write(vname,'("fhscor_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'W/m2',-3000.,3000.,0,-1)
    write(lname,'("Latent correction term ",I1.1)') n
    write(vname,'("fescor_",I1.1)') n
    call attrib(idnc,idim,3,vname,lname,'W/m2',-3000.,3000.,0,-1)
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
end if
  
return
end subroutine savetiledef

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,iarch)

use carbpools_m
use infile
use soil_m
use soilsnow_m
use vegpar_m
  
implicit none

include 'newmpar.h'
  
integer, intent(in) :: idnc,iarch
integer k,n
real, dimension(ifull) :: dat
character(len=11) vname
logical, intent(in) :: local
  
do n=1,5
  do k=1,ms
    dat=tgg(:,k)
    if (n<=maxnb) dat=unpack(ssnow%tgg(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("tgg",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wb(:,k)
    if (n<=maxnb) dat=unpack(real(ssnow%wb(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
    write(vname,'("wb",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=wbice(:,k)
    if (n<=maxnb) dat=unpack(real(ssnow%wbice(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
    write(vname,'("wbice",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  do k=1,3
    dat=tggsn(:,k)
    if (n<=maxnb) dat=unpack(ssnow%tggsn(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=smass(:,k)
    if (n<=maxnb) dat=unpack(ssnow%smass(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("smass",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=ssdn(:,k)
    if (n<=maxnb) dat=unpack(ssnow%ssdn(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=snowd/3.
    if (n<=maxnb) dat=unpack(ssnow%sdepth(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("sdepth",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.2
    if (n<=maxnb) dat=unpack(ssnow%sconds(pind(n,1):pind(n,2),k),tmap(:,n),dat)
    write(vname,'("sconds",I1.1,"_",I1.1)') k,n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end do
  dat=ssdnn
  if (n<=maxnb) dat=unpack(ssnow%ssdnn(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("ssdnn_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=real(isflag)
  if (n<=maxnb) dat=unpack(real(ssnow%isflag(pind(n,1):pind(n,2))),tmap(:,n),dat)
  write(vname,'("sflag_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snowd
  if (n<=maxnb) dat=unpack(ssnow%snowd(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("snd_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snowd
  if (n<=maxnb) dat=unpack(ssnow%osnowd(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("osnd_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=snage
  if (n<=maxnb) dat=unpack(ssnow%snage(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("snage_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=100.
  if (n<=maxnb) dat=unpack(ssnow%rtsoil(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("rtsoil_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)   
  dat=0.
  if (n<=maxnb) dat=unpack(canopy%cansto(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("cansto_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.01 ! ustar
  if (n<=maxnb) dat=unpack(canopy%us(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("us_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)  
  dat=0.
  if (n<=maxnb) dat=unpack(ssnow%pudsto(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("pudsto_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (n<=maxnb) dat=unpack(ssnow%wetfac(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("wetfac_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  if (icycle==0) then
    do k=1,ncp
      dat=cplant(:,k)
      if (n<=maxnb) dat=unpack(bgc%cplant(pind(n,1):pind(n,2),k),tmap(:,n),dat)
      write(vname,'("cplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)    
    end do
    do k=1,ncs
      dat=csoil(:,k)
      if (n<=maxnb) dat=unpack(bgc%csoil(pind(n,1):pind(n,2),k),tmap(:,n),dat)
      write(vname,'("csoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
  else
    do k=1,mplant     
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%cplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("cplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("nplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%pplant(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("pplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k=1,mlitter
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%clitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("clitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nlitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("nlitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%plitter(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("plitter",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    do k=1,msoil
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%csoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("csoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%nsoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("nsoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
      dat=0.
      if (n<=maxnb) dat=unpack(real(casapool%psoil(pind(n,1):pind(n,2),k)),tmap(:,n),dat)
      write(vname,'("psoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local,.true.)
    end do
    dat=0.
    if (n<=maxnb) dat=unpack(real(casamet%glai(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("glai_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(phen%phase(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("phenphase_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%clabile(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("clabile_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%nsoilmin(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("nsoilmin_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoillab(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("psoillab_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoilsorb(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("psoilsorb_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
    dat=0.
    if (n<=maxnb) dat=unpack(real(casapool%psoilocc(pind(n,1):pind(n,2))),tmap(:,n),dat)
    write(vname,'("psoilocc_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local,.true.)
  end if
  dat=0.
  if (n<=maxnb) dat=unpack(canopy%fhs_cor(pind(n,1):pind(n,2)),tmap(:,n),dat)
  write(vname,'("fhscor_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
  dat=0.
  if (n<=maxnb) dat=unpack(real(canopy%fes_cor(pind(n,1):pind(n,2))),tmap(:,n),dat)
  write(vname,'("fescor_",I1.1)') n
  call histwrt3(dat,vname,idnc,iarch,local,.true.)
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
subroutine cableinflow(inflow,rate)

use soil_m

implicit none

include 'newmpar.h'

integer nb, k
real, dimension(ifull), intent(in) :: rate
real, dimension(ifull), intent(inout) :: inflow
real, dimension(ifull) :: delflow
real, dimension(mp) :: xx, ll, delxx, ratepack

if ( mp <= 0 ) return

do nb = 1,maxnb
  xx(pind(nb,1):pind(nb,2)) = pack( inflow(1:ifull), tmap(:,nb) )
  ratepack(pind(nb,1):pind(nb,2)) = pack( rate(1:ifull), tmap(:,nb) )
end do
delxx(1:mp) = 0.
do k = 1,cbm_ms
  ll(1:mp) = max( soil%sfc(1:mp)-real(ssnow%wb(1:mp,k)), 0. )*1000.*soil%zse(k)
  ll(1:mp) = ll(1:mp)*ratepack(1:mp)
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
subroutine casa_readpoint(casafile)

use cc_mpi
use infile

implicit none

include 'newmpar.h'
include 'parmgeom.h'

integer ncstatus, ncid, varid, tilg
integer n
integer, dimension(2) :: spos, npos
real tlat, tlon, tschmidt
real, dimension(:,:), allocatable :: dumg
real, dimension(ifull,5) :: duma
character(len=*), intent(in) :: casafile
logical tst

if ( myid == 0 ) then
  allocate( dumg(ifull_g,5) )
  write(6,*) "Reading ",trim(casafile)
  call ccnf_open(casafile,ncid,ncstatus)
  call ncmsg('CASA_readpoint',ncstatus)
  ! check dimensions and location
  call ccnf_get_attg(ncid,'lat0',tlat)
  call ccnf_get_attg(ncid,'lon0',tlon)
  call ccnf_get_attg(ncid,'schmidt0',tschmidt)
  if ( rlong0/=tlon .or. rlat0/=tlat .or. schmidt/=tschmidt ) then
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
  call ccmpi_distribute(duma,dumg)
  deallocate(dumg)
else
  call ccmpi_distribute(duma)
end if
do n = 1,maxnb
  casamet%isorder(pind(n,1):pind(n,2))  = nint(pack(duma(:,1),tmap(:,n)))
  casaflux%Nmindep(pind(n,1):pind(n,2)) = pack(duma(:,2),tmap(:,n))/365.*1.E-3
  casaflux%Nminfix(pind(n,1):pind(n,2)) = pack(duma(:,3),tmap(:,n))/365.
  casaflux%Pdep(pind(n,1):pind(n,2))    = pack(duma(:,4),tmap(:,n))/365.
  casaflux%Pwea(pind(n,1):pind(n,2))    = pack(duma(:,5),tmap(:,n))/365.
end do

where ( veg%iveg==9 .or. veg%iveg==10 ) ! crops
  ! P fertilizer =13 Mt P globally in 1994
  casaflux%Pdep = casaflux%Pdep + 0.7/365.
  ! N fertilizer =86 Mt N globally in 1994
  casaflux%Nminfix = casaflux%Nminfix + 4.3/365.
end where

if ( any(casamet%isorder<1.or.casamet%isorder>12) ) then
  write(6,*) "ERROR: Invalid isorder in ",trim(casafile)
  call ccmpi_abort(-1)
end if

end subroutine casa_readpoint

! *************************************************************************************  
! This subroutine reads the MODIS derived leaf phenology data
subroutine casa_readphen(fphen)

use cc_mpi

implicit none

include 'newmpar.h'

integer, parameter :: nphen = 8 ! was 10(IGBP). changed by Q.Zhang @01/12/2011
integer np, ilat, ivp
integer, dimension(271,mxvt) :: greenup, fall, phendoy1
integer, dimension(nphen) :: greenupx, fallx, xphendoy1
integer, dimension(nphen) :: ivtx
real :: xlat
character(len=*), intent(in) :: fphen

! initilize for evergreen PFTs
greenup  = -50
fall     = 367
phendoy1 = 2

if ( myid == 0 ) then
  write(6,*) "Reading CASA leaf phenology data"
  open(87,file=fphen,status='old')
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
call ccmpi_bcast(greenup,0,comm_world)
call ccmpi_bcast(fall,0,comm_world)
call ccmpi_bcast(phendoy1,0,comm_world)

do np = 1,mp
  ilat = nint((rad%latitude(np)+55.25)*2.) + 1
  ilat = min( 271, max( 1, ilat ) )
  ivp = veg%iveg(np)
  phen%phase(np)      = phendoy1(ilat,ivp)
  phen%doyphase(np,1) = greenup(ilat,ivp)          ! DOY for greenup
  phen%doyphase(np,2) = phen%doyphase(np,1) + 14   ! DOY for steady LAI
  phen%doyphase(np,3) = fall(ilat,ivp)             ! DOY for leaf senescence
  phen%doyphase(np,4) = phen%doyphase(np,3) + 14   ! DOY for minimal LAI season
  if ( phen%doyphase(np,2) > 365 ) phen%doyphase(np,2) = phen%doyphase(np,2) - 365
  if ( phen%doyphase(np,4) > 365 ) phen%doyphase(np,4) = phen%doyphase(np,4) - 365
end do

return
end subroutine casa_readphen

end module cable_ccam

