! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! Thanks to Paul Ryan for OMP routines
    
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

use cable_air_module
use cable_albedo_module
use cable_canopy_module
use cable_ccam3
use cable_ccam4
use cable_common_module
use cable_data_module
use cable_def_types_mod, cbm_ms => ms
use cable_gw_hydro_module
use cable_optimise_JV_module, only : optimise_JV
use cable_radiation_module
use cable_roughness_module
use cable_soil_snow_module
use casa_cnp_module
use casadimension
use casaparm, xroot => froot
use casavariable
use newmpar_m, only : mxvt
use phenvariable
use popmodule, only : pop_init, popstep
use pop_types
use sli_main_mod

implicit none

private
public sib4, cable_version
public loadcbmparm, cbmparm, loadtile, defaulttile, savetiledef, savetile, newcbmwb
public cablesettemp, cableinflow, cbmemiss
public proglai, progvcmax, maxtile, soil_struc, cable_pop, ccycle, cable_potev
public fwsoil_switch, cable_litter, gs_switch, cable_enablefao
public smrf_switch, strf_switch, cable_gw_model, cable_roughness
public POP_NPATCH, POP_NCOHORT, POP_AGEMAX
public calc_wt_ave, calc_wt_flux
public mplant,mlitter,msoil
public cable_casatile

! CABLE biophysical options
real, save :: cable_version = 6608. ! expected version id for input data
                                    ! 6608 includes update to ice albedo
                                    ! 3939 had fixes for soil albedo

integer, save :: soil_struc      = 0          ! 0 default, 1 SLI soil model
integer, save :: fwsoil_switch   = 0          ! 0 default, 1 non-linear, 2 Lai and Ktaul, 3 Haverd2013
integer, save :: cable_litter    = 0          ! 0 off, 1 on
integer, save :: gs_switch       = 0          ! 0 leuning, 1 medlyn
integer, save :: smrf_switch     = 4          ! 1 CASA-CNP, 2 SOLIN, 3 TRIFFID, 4 Trudinger2016(default),
                                              ! 5 DAMM (Soil Moist Respiration Function)
integer, save :: strf_switch     = 4          ! 1 CASA-CNP, 2 K1995, 3 PnET-CN, 4 LT1994(default),
                                              ! 5 DAMM (Soil Temp Respiration Function)
integer, save :: cable_gw_model  = 0          ! 0 off, 1 GW_Hydro
integer, save :: cable_roughness = 0          ! 0 default, 1 new
integer, save :: cable_potev     = 1          ! 0 Penman Monteith, 1 Humidity Deficit
integer, save :: cable_enablefao = 1          ! 0 off, 1 on when calculating potential evaporation
! CABLE biochemical options
integer, save :: ccycle          = 0          ! 0 off, 1 (C), 2 (CN), 3 (CNP)
integer, save :: proglai         = 0          ! 0 prescribed, 1 prognostic LAI
integer, save :: progvcmax       = 0          ! 0 prescribed, 1 prognostic vcmax (standard)
! CABLE POP options
integer, save :: cable_pop       = 0          ! 0 off, 1 on
! CABLE POP parameters
integer, parameter :: maxtile    = 9          ! maximum possible number of tiles in a grid box
                                              ! (1-5=natural/secondary, 6-7=pasture/rangeland, 8-9=crops)
integer, parameter :: COLDEST_DAY_NHEMISPHERE = 355
integer, parameter :: COLDEST_DAY_SHEMISPHERE = 172
integer, save :: POP_NPATCH      = -1
integer, save :: POP_NLAYER      = -1
integer, save :: POP_NCOHORT     = -1
integer, save :: POP_HEIGHT_BINS = -1
integer, save :: POP_NDISTURB    = -1
integer, save :: POP_AGEMAX      = -1

integer, save :: maxnb                        ! maximum number of tiles within a gridbox
integer, save :: mp_global                    ! maximum number of points on this process (sum of all land tiles)

integer, dimension(:), allocatable, target, save :: cveg ! CABLE vegetation index for each point
real, dimension(:), allocatable, target, save :: sv      ! area fraction for each point
real, dimension(:), allocatable, target, save :: vl2     ! prescribed leaf area index (LAI) for point

! Initial values for CNP pools over 3*plant, 3*litter and 3*soil (=27 pools in total)
real(kind=8), dimension(mxvt), save :: cleaf  =(/ 384.6037_8,     273._8, 96.59814_8, 150.2638_8,      88._8, 137.1714_8, &
                                                  137.1714_8, 137.1714_8,     160._8,     160._8,       0._8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cwood  =(/ 7865.396_8,   11451._8, 5683.402_8, 10833.74_8,     372._8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cfroot =(/     250._8,    2586._8,     220._8,     220._8,     140._8,     263._8, &
                                                      263._8,     263._8,     240._8,     240._8,       0._8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cmet   =(/ 6.577021_8, 44.63457_8, 7.127119_8, 10.97797_8, 3.229374_8, 28.57245_8, &
                                                  28.57245_8, 28.57245_8, 28.57245_8, 28.57245_8,       0._8,       0._8, &
                                                        0._8, 1.457746_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cstr   =(/ 209.1728_8, 433.7626_8, 277.7733_8, 312.5492_8, 39.44449_8, 50.91091_8, &
                                                  50.91091_8, 50.91091_8, 50.91091_8, 50.91091_8,       0._8,       0._8, &
                                                        0._8, 4.956338_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: ccwd   =(/ 606.0255_8, 1150.765_8, 776.7331_8, 888.5864_8, 111.5864_8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8,       0._8, &
                                                        0._8, 28.44085_8,       0._8,       0._8,       0._8 /) 
real(kind=8), dimension(mxvt), save :: cmic   =(/  528.664_8, 11.37765_8, 597.0785_8, 405.5554_8, 168.0451_8, 425.6431_8, &
                                                  425.6431_8, 425.6431_8, 512.4247_8, 512.4247_8,       0._8,       0._8, &
                                                        0._8, 57.77585_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cslow  =(/ 13795.94_8, 311.8092_8, 16121.12_8, 11153.25_8, 4465.478_8, 5694.437_8, &
                                                  5694.437_8, 5694.437_8, 6855.438_8, 6855.438_8,       0._8,       0._8, &
                                                        0._8, 1325.052_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: cpass  =(/ 4425.396_8 ,13201.81_8, 5081.802_8, 5041.192_8, 1386.477_8,  4179.92_8, &
                                                   4179.92_8,  4179.92_8, 5032.137_8, 5032.137_8,       0._8,       0._8, &
                                                        0._8, 517.1719_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nleaf  =(/ 7.541249_8,      9.9_8, 1.609969_8, 3.756594_8, 2.933333_8, 4.572381_8, &
                                                  4.572381_8, 4.572381_8, 5.333333_8, 5.333333_8,       0._8,       0._8, &
                                                        0._8,      0.5_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nwood  =(/ 31.46159_8,    102._8,  22.73361_8, 80.24989_8, 2.755555_8,       0._8, &
                                                        0._8,      0._8,        0._8,       0._8,       0._8,       0._8, &
                                                        0._8, 0.125926_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nfroot =(/ 6.097561_8,      38._8, 5.365854_8, 5.365854_8, 3.414634_8, 6.414634_8, &
                                                  6.414634_8, 6.414634_8, 5.853659_8, 5.853659_8,       0._8,       0._8, &
                                                        0._8, 1.536585_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nmet   =(/ 0.064481_8,  0.74391_8, 0.059393_8, 0.137225_8, 0.053823_8, 0.476208_8, &
                                                  0.476208_8, 0.476208_8, 0.476208_8, 0.476208_8,       0._8,       0._8, &
                                                        0._8, 0.018222_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nstr   =(/ 1.394485_8, 2.891751_8, 1.851822_8, 2.083661_8, 0.262963_8, 0.339406_8, &
                                                  0.339406_8, 0.339406_8, 0.339406_8, 0.339406_8,       0._8,       0._8, &
                                                        0._8, 0.033042_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: ncwd   =(/ 2.424102_8, 8.524183_8, 3.106932_8, 6.581996_8, 0.826566_8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8,       0._8, &
                                                        0._8, 0.210673_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nmic   =(/  52.8664_8, 1.137765_8, 59.70785_8, 40.55554_8, 16.80451_8, 42.56431_8, &
                                                  42.56431_8, 42.56431_8, 51.24247_8, 51.24247_8,       0._8,       0._8, &
                                                        0._8, 5.777585_8       ,0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: nslow  =(/ 919.7293_8, 20.78728_8, 1074.741_8, 743.5501_8, 297.6985_8, 379.6291_8, &
                                                  379.6291_8, 379.6291_8, 457.0292_8, 457.0292_8,       0._8,       0._8, &
                                                        0._8, 88.33682_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: npass  =(/ 295.0264_8, 880.1209_8, 338.7868_8, 336.0795_8,  92.4318_8, 278.6613_8, &
                                                  278.6613_8, 278.6613_8, 335.4758_8, 335.4758_8,       0._8,       0._8, &
                                                        0._8, 34.47813_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpleaf =(/ 0.191648_8,   0.415_8,  0.115988_8, 0.135453_8, 0.022821_8,  0.15125_8, &
                                                   0.15125_8,  0.15125_8,  0.15125_8, 0.15125_8,        0._8,       0._8, &
                                                        0._8,    0.007_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpwood =(/ 0.953979_8,    5.88_8, 0.64438_8,   2.424778_8,       0._8,       0._8, &
                                                        0._8,      0._8,      0._8,         0._8,       0._8,       0._8, &
                                                        0._8,      0._8,      0._8,         0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpfroot=(/ 0.076659_8,    1.95_8, 0.080548_8,  0.141097_8, 0.037083_8,  0.15125_8, &
                                                   0.15125_8,  0.15125_8, 0.15125_8,   0.15125_8,       0._8,       0._8, &
                                                        0._8,  0.00875_8,      0._8,        0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpmet  =(/ 0.004385_8, 0.029756_8, 0.004751_8, 0.007319_8, 0.002153_8, 0.019048_8, &
                                                  0.019048_8, 0.019048_8, 0.019048_8, 0.019048_8,       0._8,       0._8, &
                                                        0._8, 0.000972_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpstr  =(/ 0.069724_8, 0.144588_8, 0.092591_8, 0.104183_8, 0.013148_8,  0.01697_8, &
                                                   0.01697_8,  0.01697_8,  0.01697_8,  0.01697_8,       0._8,       0._8, &
                                                        0._8, 0.001652_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpcwd  =(/ 0.101004_8, 0.191794_8, 0.129456_8, 0.148095_8, 0.018598_8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8,       0._8, &
                                                        0._8,       0._8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpmic  =(/ 6.872632_8 , 0.14791_8, 7.762021_8,  5.27222_8, 2.184586_8, 5.533361_8, &
                                                  5.533361_8, 5.533361_8, 6.661522_8, 6.661522_8,       0._8,       0._8, &
                                                        0._8, 0.751086_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xpslow =(/ 119.5648_8, 2.702347_8, 139.7164_8, 96.66152_8, 38.70081_8, 49.35178_8, &
                                                  49.35178_8, 49.35178_8,  59.4138_8,  59.4138_8,       0._8,       0._8, &
                                                        0._8, 11.48379_8,       0._8,       0._8,       0._8 /)
real(kind=8), dimension(mxvt), save :: xppass =(/ 38.35343_8, 114.4157_8, 44.04228_8, 43.69033_8, 12.01613_8, 36.22598_8, &
                                                  36.22598_8, 36.22598_8, 43.61185_8, 43.61185_8,       0._8,       0._8, &
                                                        0._8, 4.482157_8,       0._8,       0._8,       0._8 /)

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
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use tracers_m
use vegpar_m
use work2_m, only : qsttg,zo,zoh,zoq,theta,vmod,wetfac
use work3_m, only : ga
use zenith_m

integer, dimension(maxtile,2) :: ltind
integer :: lmaxnb, tile, is, ie, js, je
integer :: ico2, igas
integer :: iyr, imo, iday, mdays, k
real, dimension(imax,mlitter) :: lclitter, lnilitter, lplitter
real, dimension(imax,mplant) :: lcplant, lniplant, lpplant
real, dimension(imax,msoil) :: lcsoil, lnisoil, lpsoil
real, dimension(imax,ms) :: ltgg, lwb, lwbice
real, dimension(imax,ms) :: lwb_clim
real, dimension(imax,3) :: lsmass, lssdn, ltggsn
real, dimension(imax) :: lcnbp, lcnpp, lfnee, lfpn, lfrd, lfrp, lfrpr, lfrpw, lfrs
real, dimension(imax) :: latmco2
real, dimension(imax) :: lfevc,lplant_turnover,lplant_turnover_wood
real :: x
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
      lcnbp = cnbp(is:ie)
      lcnpp = cnpp(is:ie)
      lcplant = cplant(is:ie,:)
      lclitter = clitter(is:ie,:)
      lcsoil = csoil(is:ie,:)
      lniplant = niplant(is:ie,:)
      lnilitter = nilitter(is:ie,:)
      lnisoil = nisoil(is:ie,:)
      lpplant = pplant(is:ie,:)
      lplitter = plitter(is:ie,:)
      lpsoil = psoil(is:ie,:)
      lfnee = fnee(is:ie)
      lfpn = fpn(is:ie)
      lfrd = frd(is:ie)
      lfrp = frp(is:ie)
      lfrpr = frpr(is:ie)
      lfrpw = frpw(is:ie)
      lfrs = frs(is:ie)
      if ( diaglevel_carbon>0 ) then
        lfevc = fevc(is:ie)
        lplant_turnover = plant_turnover(is:ie)
        lplant_turnover_wood = plant_turnover_wood(is:ie)
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
    latmco2 = 1.E6*rrvco2          ! from radiative CO2 forcings
    ico2 = 0
    if ( ngas>0 ) then
      do igas = 1,ngas
        if ( trim(tractype(igas))=='online' .and. trim(tracname(igas))=='cbmnep' ) then
          ico2 = igas
          exit
        end if
      end do
      if ( ico2>0 ) then
        latmco2 = tr(1:imax,1,ico2) ! use interactive tracers
      end if
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
    !call setp(climate,lclimate,tile)  ! disable climate

    call sib4_work(albnirdif(is:ie),albnirdir(is:ie),albnirsav(is:ie),albvisdif(is:ie),albvisdir(is:ie),                 &
                   albvissav(is:ie),cansto(is:ie),cdtq(is:ie),cduv(is:ie),lclitter,lcnbp,                                &
                   lcnpp,condg(is:ie),conds(is:ie),condx(is:ie),lcplant,lcsoil,eg(is:ie),epot(is:ie),fbeamnir(is:ie),    &
                   fbeamvis(is:ie),fg(is:ie),lfnee,lfpn,lfrd,lfrp,                                                       &
                   lfrpr,lfrpw,lfrs,fwet(is:ie),ga(is:ie),isflag(is:ie),land(is:ie),lmaxnb,mp,lnilitter,lniplant,        &
                   lnisoil,lplitter,lpplant,ps(is:ie),lpsoil,qg(is:ie,1),qsttg(is:ie),rgsave(is:ie),                     &
                   rlatt(is:ie),rlongg(is:ie),rnet(is:ie),rsmin(is:ie),runoff(is:ie),runoff_surface(is:ie),              &
                   sigmf(is:ie),lsmass,snage(is:ie),snowd(is:ie),snowmelt(is:ie),lssdn,ssdnn(is:ie),                     &
                   sv(js:je),sgdn(is:ie),swrsave(is:ie),t(is:ie,1),ltgg,ltggsn,theta(is:ie),ltind,ltmap,                 &
                   latmco2,tss(is:ie),ustar(is:ie),vlai(is:ie),vl2(js:je),vmod(is:ie),                                   &
                   lwb,lwbice,wetfac(is:ie),zo(is:ie),zoh(is:ie),zoq(is:ie),evspsbl(is:ie),sbl(is:ie),lair,lbal,c,       &
                   lcanopy,lcasabal,casabiome,lcasaflux,lcasamet,lcasapool,lclimate,lmet,lphen,lpop,lrad,lrough,lsoil,   &
                   lssnow,lsum_flux,lveg,lfevc,lplant_turnover,lplant_turnover_wood,lwb_clim,wtd(is:ie),imax)

    smass(is:ie,:) = lsmass
    ssdn(is:ie,:) = lssdn
    tgg(is:ie,:) = ltgg
    tggsn(is:ie,:) = ltggsn
    wb(is:ie,:) = lwb
    wbice(is:ie,:) = lwbice
    if ( ccycle/=0 ) then
      cnbp(is:ie) = lcnbp
      cnpp(is:ie) = lcnpp
      cplant(is:ie,:) = lcplant
      clitter(is:ie,:) = lclitter
      csoil(is:ie,:) = lcsoil
      niplant(is:ie,:) = lniplant
      nilitter(is:ie,:) = lnilitter
      nisoil(is:ie,:) = lnisoil
      pplant(is:ie,:) = lpplant
      plitter(is:ie,:) = lplitter
      psoil(is:ie,:) = lpsoil
      fnee(is:ie) = lfnee
      fpn(is:ie) = lfpn
      frd(is:ie) = lfrd
      frp(is:ie) = lfrp
      frpr(is:ie) = lfrpr
      frpw(is:ie) = lfrpw
      frs(is:ie) = lfrs
      if ( diaglevel_carbon > 0 ) then
        fevc(is:ie) = lfevc
        plant_turnover(is:ie) = lplant_turnover
        plant_turnover_wood(is:ie) = lplant_turnover_wood
      end if
    end if  
  end if ! mp>0

end do

return
end subroutine sib4

! CABLE-CCAM interface
subroutine sib4_work(albnirdif,albnirdir,albnirsav,albvisdif,albvisdir,albvissav,cansto,cdtq,cduv,clitter,cnbp, &
                     cnpp,condg,conds,condx,cplant,csoil,eg,epot,fbeamnir,fbeamvis,fg,fnee,fpn,frd,frp,frpr,    &
                     frpw,frs,fwet,ga,isflag,land,maxnb,mp,nilitter,niplant,nisoil,plitter,pplant,ps,           &
                     psoil,qg,qsttg,rgsave,rlatt,rlongg,rnet,rsmin,runoff,runoff_surface,sigmf,smass,           &
                     snage,snowd,snowmelt,ssdn,ssdnn,sv,sgdn,swrsave,t,tgg,tggsn,theta,tind,tmap,atmco2,tss,    &
                     ustar,vlai,vl2,vmod,wb,wbice,wetfac,zo,zoh,zoq,evspsbl,sbl,air,bal,c,canopy,               &
                     casabal,casabiome,casaflux,casamet,casapool,climate,met,phen,pop,rad,rough,soil,ssnow,     &
                     sum_flux,veg,fevc,plant_turnover,plant_turnover_wood,wb_clim,wtd,imax)

use const_phys
use dates_m
use estab, only : qsat
use infile, only : getzinp
use parm_m
use sigs_m
use soil_m, only : zmin
use zenith_m, only : solargh, zenith
  
integer, intent(in) :: imax, maxnb, mp
integer jyear, jmonth, jday, jhour, jmin, idoy
integer k, mins, nb, j
integer is, ie, casaperiod, npercasa
integer lalloc
integer mp_POP
integer, dimension(maxtile,2), intent(in) :: tind
integer, dimension(imax), intent(inout) :: isflag
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
real, dimension(imax), intent(out) :: wtd
real, dimension(imax), intent(in) :: condg, conds, condx, fbeamnir, fbeamvis, ps, rgsave, rlatt, rlongg
real, dimension(imax), intent(in) :: sgdn, swrsave, theta, vmod
real, dimension(imax) :: coszro2, taudar2, tmps, qsttg_land
real, dimension(imax) :: evspsbl_l, sbl_l
real, dimension(mp), intent(in) :: sv, vl2
real, dimension(mp) :: cc1, cc2, dumt, dump, qsatfvar, ssnowpotev
real, dimension(mp) :: wbclim_pack
real, dimension(imax), intent(inout) :: runoff, runoff_surface, snowmelt
real, dimension(imax), intent(inout) :: evspsbl, sbl
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
fjd = real(mins)/1440.
if ( leap==0 ) then ! 365 day calendar
  idoy = int(fjd)
else if ( leap==1 ) then ! 365/366 day calendar
  idoy = mod(int(fjd),365)
else if ( leap==2 ) then ! 360 day calendar
  idoy = int(fjd*365./360.) ! CABLE expects standard calendar  
else    
  write(6,*) "ERROR: Unknown option for leap = ",leap
  stop
end if
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
  met%qv(is:ie)        = real(pack(qg(1:imax), tmap(:,nb)), 8)        ! specific humidity in kg/kg
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
  rough%za_tq(is:ie)   = max( rough%za_tq(is:ie), veg%hc(is:ie)+1._8 )
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
call setlai(sigmf,jmonth,jday,jhour,jmin,mp,sv,vl2,casamet,veg,imax,tind,tmap,maxnb)

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
    if ( cable_gw_model==1 ) then
      call soil_snow_gw(dtr8,soil,ssnow,canopy,met,bal,veg)        
    else    
      call soil_snow(dtr8,soil,ssnow,canopy,met,bal,veg)
    end if  
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

! MJT suggestion
canopy%cdtq =  max( 0._8, canopy%cdtq )
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

#ifdef debug
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
if ( any( ssnow%tgg(:,1)>425. ) ) then
  write(6,*) "WARN: tgg1>425. after CABLE"
end if
if ( any( ssnow%tggsn(:,1)>425. ) ) then
  write(6,*) "WARN: tggsn1>425. after CABLE"
end if
if ( any( canopy%tv(:)>425. ) ) then
  write(6,*) "WARN: tv>425. after CABLE" 
end if
#endif


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
      call phenology(idoy,veg,phen)
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
! Special
if ( nspecial == 51 ) then
  ! reset soil moisture to prescribed climatology  
  do k = 1,ms
    do nb = 1,maxnb  
      is = tind(nb,1)
      ie = tind(nb,2)        
      wbclim_pack(is:ie) = pack(wb_clim(:,k),tmap(:,nb))
    end do  
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
wtd = 0.

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
  ! Replace potev with Penman_Monteith
  if ( cable_potev==1 .and. cable_enablefao==1 ) then
    dumt(is:ie) = real( met%tvair(is:ie) )
    dump(is:ie) = real( met%pmb(is:ie) )*100. ! convert from mb to Pa
    qsatfvar(is:ie) = qsat(dump(is:ie),dumt(is:ie))
    cc1(is:ie) = real( air%dsatdk(is:ie)/(air%dsatdk(is:ie)+air%psyc(is:ie) ) )
    cc2(is:ie) = real( air%psyc(is:ie)/(air%dsatdk(is:ie)+air%psyc(is:ie) ) )
    ssnowpotev(is:ie) = cc1(is:ie) * real(canopy%fns(is:ie) - canopy%ga(is:ie)) +            &
                        cc2(is:ie) * real(air%rho(is:ie) * air%rlam(is:ie))*(qsatfvar(is:ie) &
                        - real(met%qvair(is:ie)))/real(ssnow%rtsoil(is:ie))
    epot = epot + unpack(sv(is:ie)*ssnowpotev(is:ie),tmap(:,nb),0.)             ! diagnostic
  else
    epot = epot + unpack(sv(is:ie)*real(ssnow%potev(is:ie)),tmap(:,nb),0.)      ! diagnostic in history file
  end if
  vlai = vlai + unpack(sv(is:ie)*real(veg%vlai(is:ie)),tmap(:,nb),0.)
  rsmin = rsmin + unpack(sv(is:ie)*real(canopy%gswx_T(is:ie)),tmap(:,nb),0.)    ! diagnostic in history file
  wtd = wtd + unpack(sv(is:ie)*real(ssnow%wtd(is:ie))/1000.,tmap(:,nb),0.)  
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


! MJT notes - ustar, cduv, fg and eg are passed to the boundary layer turbulence scheme
! zoh, zoq and zo are passed to the scrnout diagnostics routines
! rsmin is typically used by CTM

where ( land(1:imax) )
  zo      = zmin*exp(-1./sqrt(zo))
  zoh     = zo/7.4
  zoq     = zoh
  ustar   = sqrt(cduv)*vmod  
  cduv    = cduv*vmod           ! cduv is Cd*vmod in CCAM
  cdtq    = cdtq*vmod
  tss     = tss**0.25
  rsmin   = 1./rsmin
  evspsbl = evspsbl + dt*evspsbl_l
  sbl     = sbl + dt*sbl_l
  ! update albedo and tss before calculating net radiation
  albvissav = fbeamvis*albvisdir + (1.-fbeamvis)*albvisdif ! for nrad=4
  albnirsav = fbeamnir*albnirdir + (1.-fbeamnir)*albnirdif ! for nrad=4 
  isflag(:) = nint(tmps(:)) ! tmps is average isflag
end where
qsttg_land(:) = qsat(ps,tss) ! must wait for tss to be updated first
where ( land(1:imax) )
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

call cbmemiss_work(trsrc,mvegt,mode,imax,tdata(tile)%tind,tdata(tile)%tmap,tdata(tile)%maxnb)
  
return
end subroutine cbmemiss

subroutine cbmemiss_work(trsrc,mvegt,mode,imax,tind,tmap,maxnb)
  
use parm_m
  
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
subroutine setlai(sigmf,jmonth,jday,jhour,jmin,mp,sv,vl2,casamet,veg,imax,tind,tmap,maxnb)

use cc_mpi
use dates_m
use parm_m
  
integer, intent(in) :: jmonth,jday,jhour,jmin,mp
integer, intent(in) :: imax
integer, optional :: maxnb
integer monthstart, nb, is, ie
integer, dimension(12), parameter :: imonth = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
integer, dimension(maxtile,2), intent(in), optional :: tind
real, dimension(mp), intent(in) :: sv, vl2
real, parameter :: vextkn = 0.4
real, dimension(imax), intent(out) :: sigmf
real, dimension(mp) :: a0, a1, a2, aa, bb, cc, mp1, mp2, c2, c3, c4
real, dimension(mp) :: dummy_unpack
real x
logical, dimension(imax,maxtile), intent(in), optional :: tmap
type(casa_met), intent(in) :: casamet
type(veg_parameter_type), intent(inout) :: veg

select case( proglai )
  case(-1,0) ! Prescribed LAI
    veg%vlai = vl2(:)  
    where ( veg%iveg<14 )
      veg%vlai = max( veg%vlai, 0.02_8 )
    elsewhere
      veg%vlai = 1.E-8_8
    end where

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

  do np = 1,mp
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
subroutine loadcbmparm(fveg,fphen,casafile, &
                       ivs,svs,vlin,casapoint,greenup,fall,phendoy1)

use cc_mpi
use newmpar_m
use parm_m

integer n
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs,vlin
real, dimension(ifull,maxtile), intent(out) :: casapoint
real, dimension(ifull) :: svs_sum
real cableformat
integer, dimension(271,mxvt), intent(out) :: greenup, fall, phendoy1
character(len=*), intent(in) :: fveg,fphen,casafile

! read CABLE biome and LAI data
if ( myid==0 ) then
  write(6,*) "Reading tiled surface data for CABLE"
  call vegta(ivs,svs,vlin,fveg,cableformat)
else
  call vegtb(ivs,svs,vlin,cableformat)
end if
svs_sum = sum(svs,dim=2)
do n = 1,maxtile
  svs(:,n) = svs(:,n)/svs_sum(:)
end do

if ( abs(cableformat-1.)<1.e-20 ) then
  if ( myid==0 ) write(6,*) "Procesing CSIRO PFTs"    
else
  if ( myid==0 ) write(6,*) "Processing IGBP and converting to CSIRO PFTs"
  call convertigbp(ivs,svs,vlin)
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

if ( any(vlin>10.) ) then
  write(6,*) "ERROR: LAI is out of range"
  write(6,*) "vlin ",maxval(vlin)
  call ccmpi_abort(-1)
end if

return
end subroutine loadcbmparm

! legacy code for IGBP vegetation classes    
subroutine convertigbp(ivs,svs,vlin)

use cc_mpi
use const_phys
use latlong_m
use newmpar_m
use soil_m

integer, dimension(ifull,maxtile), intent(inout) :: ivs
integer, dimension(1) :: pos
integer iq, n, ipos, iv
real, dimension(ifull,maxtile), intent(inout) :: svs,vlin
real, dimension(18) :: newlai
real, dimension(18) :: newgrid
real fc3, fc4, ftu, fg3, fg4, clat, nsum
real xp

do iq = 1,ifull
  if ( land(iq) ) then
    newgrid(:)  = 0.
    newlai(:) = 0.      
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
          newlai(ivs(iq,n))=newlai(ivs(iq,n))+svs(iq,n)*vlin(iq,n)
        case (5)
          if (abs(clat)>25.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5
            newlai(1)=newlai(1)+svs(iq,n)*0.5*vlin(iq,n)
            newgrid(4)=newgrid(4)+svs(iq,n)*0.5
            newlai(4)=newlai(4)+svs(iq,n)*0.5*vlin(iq,n)
          else if (abs(clat)>24.5) then
            xp=abs(clat)-24.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.5*xp
            newlai(1)=newlai(1)+svs(iq,n)*0.5*vlin(iq,n)*xp
            newgrid(4)=newgrid(4)+svs(iq,n)*(1.-0.5*xp)
            newlai(4)=newlai(4)+svs(iq,n)*vlin(iq,n)*(1.-0.5*xp)
          else
            newgrid(4)=newgrid(4)+svs(iq,n)
            newlai(4)=newlai(4)+svs(iq,n)*vlin(iq,n)
          end if
        case (6)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.8
          newlai(5)=newlai(5)+svs(iq,n)*0.8*vlin(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.2*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.2*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.2*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.2*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.2*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.2*ftu*vlin(iq,n)
        case (7)
          newgrid(5)=newgrid(5)+svs(iq,n)*0.2
          newlai(5)=newlai(5)+svs(iq,n)*0.2*vlin(iq,n)
          newgrid(6)=newgrid(6)+svs(iq,n)*0.8*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.8*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.8*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.8*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.8*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.8*ftu*vlin(iq,n)
        case (8)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4
            newlai(1)=newlai(1)+svs(iq,n)*0.4*vlin(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.4*xp
            newlai(1)=newlai(1)+svs(iq,n)*vlin(iq,n)*0.4*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4*(1.-xp)
            newlai(18)=newlai(18)+svs(iq,n)*vlin(iq,n)*0.4*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.4
            newlai(18)=newlai(18)+svs(iq,n)*0.4*vlin(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.6*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.6*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.6*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.6*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.6*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.6*ftu*vlin(iq,n)
        case (9)
          if (abs(clat)>40.5) then
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1
            newlai(1)=newlai(1)+svs(iq,n)*0.1*vlin(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(1)=newgrid(1)+svs(iq,n)*0.1*xp
            newlai(1)=newlai(1)+svs(iq,n)*vlin(iq,n)*0.1*xp
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1*(1.-xp)
            newlai(18)=newlai(18)+svs(iq,n)*vlin(iq,n)*0.1*(1.-xp)
          else
            newgrid(18)=newgrid(18)+svs(iq,n)*0.1
            newlai(18)=newlai(18)+svs(iq,n)*0.1*vlin(iq,n)
          end if
          newgrid(6)=newgrid(6)+svs(iq,n)*0.9*fg3
          newlai(6)=newlai(6)+svs(iq,n)*0.9*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*0.9*fg4
          newlai(7)=newlai(7)+svs(iq,n)*0.9*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*0.9*ftu
          newlai(8)=newlai(8)+svs(iq,n)*0.9*ftu*vlin(iq,n)
        case (10)
          newgrid(6)=newgrid(6)+svs(iq,n)*fg3
          newlai(6)=newlai(6)+svs(iq,n)*fg3*vlin(iq,n)
          newgrid(7)=newgrid(7)+svs(iq,n)*fg4
          newlai(7)=newlai(7)+svs(iq,n)*fg4*vlin(iq,n)
          newgrid(8)=newgrid(8)+svs(iq,n)*ftu
          newlai(8)=newlai(8)+svs(iq,n)*ftu*vlin(iq,n)
        case (12,14)
          newgrid(9)=newgrid(9)+svs(iq,n)*fc3
          newlai(9)=newlai(9)+svs(iq,n)*fc3*vlin(iq,n)
          newgrid(10)=newgrid(10)+svs(iq,n)*fc4
          newlai(10)=newlai(10)+svs(iq,n)*fc4*vlin(iq,n)
        case (13)
          newgrid(15)=newgrid(15)+svs(iq,n)
          newlai(15)=newlai(15)+svs(iq,n)*vlin(iq,n)
        case (15)
          newgrid(17)=newgrid(17)+svs(iq,n)
          newlai(17)=newlai(17)+svs(iq,n)*vlin(iq,n)
        case (16)
          newgrid(14)=newgrid(14)+svs(iq,n)
          newlai(14)=newlai(14)+svs(iq,n)*vlin(iq,n)
        case (17)
          newgrid(16)=newgrid(16)+svs(iq,n)
          newlai(16)=newlai(16)+svs(iq,n)*vlin(iq,n)
        case DEFAULT
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq,ivs,land=",myid,iq,ivs(iq,n),land(iq)
          call ccmpi_abort(-1)
      end select
    end do
    where ( newgrid(:)>0. )
      newlai(:) = newlai(:)/newgrid(:)
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
    vlin(iq,:)     = 0.
    do iv = 1,18
      if ( newgrid(iv)>0. ) then
        n = n + 1
        ivs(iq,n)      = iv
        svs(iq,n)      = newgrid(iv)
        vlin(iq,n)     = newlai(iv)
      end if
    end do

  end if
end do

return
end subroutine convertigbp

subroutine cbmparm(ivs,svs,vlin,casapoint,greenup,fall,phendoy1,fcasapft)

use carbpools_m
use cc_mpi
use const_phys
use infile
use latlong_m
use newmpar_m
use nsibd_m
use parm_m
use pop_constants, only : NPATCH, NLAYER, NCOHORT_MAX, HEIGHT_BINS, NDISTURB, AGEMAX
use pbl_m
use riverarrays_m
use sigs_m
use soil_m
use soilsnow_m
use soilv_m
use vegpar_m
  
integer, dimension(ifull,maxtile), intent(in) :: ivs
integer, dimension(271,mxvt), intent(in) :: greenup, fall, phendoy1
integer(kind=4), dimension(:), allocatable, save :: Iwood
integer(kind=4), dimension(:,:), allocatable, save :: disturbance_interval
integer i,iq,n,k,ipos,ilat,ivp,is,ie
integer jyear,jmonth,jday,jhour,jmin,mins
integer landcount
integer(kind=4) mp_POP
real ivmax, landsum
real, dimension(ifull,maxtile), intent(in) :: svs,vlin
real, dimension(ifull,5), intent(in) :: casapoint
real, dimension(ifull,2) :: albsoilsn
real, dimension(ifull) :: dummy_pack
real, dimension(ifull) :: albsoil
real, dimension(0:maxtile) :: stat_count, global_stat_count
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

if ( myid==0 .and. nmaxpr==1 ) write(6,*) "-> Define CABLE and CASA CNP arrays"

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
      landsum = sum(svs(iq,:))
      if ( landsum<0.99 .or. landsum>1.01 ) then
        write(6,*) "ERROR: CABLE tiles do not sum to 1."
        write(6,*) "myid,iq,tile,landsum ",myid,iq,tile,landsum
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
  if ( nmaxpr==1 ) write(6,*) "-> Allocating CABLE and CASA CNP arrays"
  if ( soil_struc==1 ) write(6,*) "-> Using SLI soil model"
  if ( ccycle==0 ) then
    write(6,*) "-> Using CABLE without carbon cycle"
  else if ( ccycle==1 ) then
    write(6,*) "-> Using CASA C"
  else if ( ccycle==2 ) then
    write(6,*) "-> Using CASA CN"
  else
    write(6,*) "-> Using CASA CNP"
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

icycle = ccycle

if ( mp_global>0 ) then

  climate%nyear_average = 20
  climate%nday_average = 31
    
  allocate( sv(mp_global) )
  allocate( vl2(mp_global) )
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
  allocate( dummy_unpack(mp_global) )
  
  ! Cable configuration
  select case( cable_potev )
    case(0)
      cable_user%ssnow_POTEV = "P-M" ! Penman Monteith
    case default  
      cable_user%ssnow_POTEV = "HDM" ! default Humidity Deficit
  end select    
  cable_user%MetType = "defa"    ! only 4 characters for "default"
  cable_user%diag_soil_resp = "ON"
  cable_user%leaf_respiration = "ON"
  cable_user%run_diag_level = "NONE"
  cable_user%consistency_check = .false.
  cable_user%logworker = .false.
  select case ( cable_roughness )
    case(1)
      cable_user%l_new_roughness_soil = .true.
    case default  
      cable_user%l_new_roughness_soil = .false.
  end select
  cable_user%l_rev_corr = .true.
  cable_user%gw_model = cable_gw_model==1
  cable_user%soil_thermal_fix = .true.
  cable_user%call_climate = .false.
  cable_user%phenology_switch = "MODIS"
  cable_user%finite_gm = .false.
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
  vl2 = 0.
  cveg = 0

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  do tile = 1,ntiles
    ! tile is the spatial decomposition and maxtile is the mosaic of vegetation PFTs  
    allocate(tdata(tile)%tmap(imax,maxtile))
    tdata(tile)%tmap = .false.
  end do
  
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
    call cable_pack(slope_ave,soil%slope,n)
    call cable_pack(slope_std,soil%slope_std,n)
    dummy_pack = max( vlin(:,n), 0.01 )
    call cable_pack(dummy_pack,vl2,n)
  end do
  
  soil%slope=min(max(soil%slope,1.e-8_8),0.95_8)
  soil%slope_std=min(max(soil%slope_std,1.e-8_8),0.95_8)

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
    vl2(:) = 1.e-8
  end where
  
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
  call setlai(sigmf,jmonth,jday,jhour,jmin,mp_global,sv,vl2,casamet,veg,ifull)
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
  call alloc_cbm_var(air, 0)
  call alloc_cbm_var(bgc, 0)
  call alloc_cbm_var(canopy, 0)
  call alloc_cbm_var(met, 0)
  call alloc_cbm_var(bal, 0)
  call alloc_cbm_var(rad, 0)
  call alloc_cbm_var(rough, 0)
  call alloc_cbm_var(soil, 0)
  call alloc_cbm_var(ssnow, 0)
  call alloc_cbm_var(sum_flux, 0)
  call alloc_cbm_var(veg, 0)  
  call cable_biophysic_parm(cveg)
  call cable_soil_parm(soil)
  if ( ccycle>=1 .and. ccycle<=3 ) then
    call casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)
  end if
  
end if
  
! statistics
stat_count(:) = 0.
global_stat_count(:) = 0.
do iq = 1,ifull
  if ( land(iq) ) then  
    landcount = count( svs(iq,:)>0. )
    stat_count(landcount) = stat_count(landcount) + 1.
    stat_count(0) = stat_count(0) + 1.
  end if  
end do  
call ccmpi_reduce(stat_count,global_stat_count,"sum",0,comm_world)
if ( myid==0 ) then
  write(6,*) "CABLE statistics:"
  do n = 1,maxtile
    write(6,'(A,I1.1,A,F5.1,A)') "   Percentage of gridpoints with ",n," tile(s) is ", &
        100.*global_stat_count(n)/global_stat_count(0),"%"
  end do  
end if
  
if ( myid==0 .and. nmaxpr==1 ) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine cbmparm

subroutine casa_readbiome(veg,casabiome,casapool,casaflux,casamet,phen,fcasapft)

use cc_mpi     ! CC MPI routines
use newmpar_m

type(veg_parameter_type), intent(in) :: veg
type(casa_biome),         intent(inout) :: casabiome
type(casa_pool),          intent(inout) :: casapool
type(casa_flux),          intent(inout) :: casaflux
type(casa_met),           intent(inout) :: casamet
type(phen_variable),      intent(inout) :: phen
character(len=*), intent(in) :: fcasapft

real(kind=8), dimension(mxvt,mplant) :: ratiocnplant
real(kind=8), dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real(kind=8), dimension(mso,msoil) :: rationpsoil

real(kind=8), dimension(mxvt) :: leafage,woodage,frootage,metage
real(kind=8), dimension(mxvt) :: strage,cwdage,micage,slowage,passage,slax
real(kind=8), dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real(kind=8), dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real(kind=8), dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real(kind=8), dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real(kind=8), dimension(mxvt) :: clabileage
real(kind=8), dimension(mxvt) :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase
real(kind=8), dimension(mxvt) :: xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
real(kind=8), dimension(mxvt) :: xla_to_sa,xvcmax_scalar,xdisturbance_interval
real(kind=8), dimension(mxvt) :: xDAMM_EnzPool,xDAMM_KMO2,xDAMM_KMcp,xDAMM_Ea,xDAMM_alpha
real(kind=8), dimension(mso) :: xxkplab,xxkpsorb,xxkpocc
real(kind=8), dimension(mso) :: xkmlabp,xpsorbmax,xfPleach

integer :: i, iso, nv, ierr
integer :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12,nv13
integer :: fflag=0

integer, dimension(mxvt) :: ivt2
real(kind=8), dimension(mxvt) :: kroot
real(kind=8), dimension(mxvt) :: rootdepth
real(kind=8), dimension(mxvt) :: kuptake
real(kind=8), dimension(mxvt) :: krootlen
real(kind=8), dimension(mxvt) :: kminn
real(kind=8), dimension(mxvt) :: kuplabp
real(kind=8), dimension(mxvt,mplant) :: fracnpptop
real(kind=8), dimension(mxvt,mplant) :: rmplant
real(kind=8), dimension(mxvt,mplant) :: ftransnptol
real(kind=8), dimension(mxvt,mplant) :: fracligninplant
real(kind=8), dimension(mxvt) :: glaimax
real(kind=8), dimension(mxvt) :: glaimin
real(kind=8), dimension(mxvt) :: xkleafcoldexp
real(kind=8), dimension(mxvt) :: xkleafdryexp
real(kind=8), dimension(mxvt,mplant) :: rationcplantmin
real(kind=8), dimension(mxvt,mplant) :: rationcplantmax
real(kind=8), dimension(mxvt,mplant) :: ftranspptol
real(kind=8), dimension(mxvt) :: tkshed

if ( trim(fcasapft) /= '' ) fflag = 1
call ccmpi_bcast(fflag,0,comm_world)

if ( fflag == 1  ) then
  if ( myid == 0 ) then
    write(6,*) "-> Using user defined CASA PFT parameter tables"
    open(86,file=fcasapft,status='old',action='read',iostat=ierr)
    call ccmpi_abort(-1)

    write(6,*) "-> Reading ",trim(fcasapft)

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
    
    if ( cable_pop==1 ) then
      if ( cwood(1)>0.02_8 ) then
        write(6,*) "ERROR: casapftfile is not configured for POP"
        call ccmpi_abort(-1)
      end if
    else
      if ( cwood(1)<0.02_8 ) then
        write(6,*) "ERROR: casapftfile is configured for POP"
        call ccmpi_abort(-1)          
      end if
    end if

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

    leafage =(/ 2.0_8, 1.5_8, 1.0_8, 1.0_8, 1.0_8, 0.8_8, 0.8_8, 1.0_8,      0.8_8,      0.8_8, 1.0_8, &
                1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    woodage =(/ 70._8, 60._8, 80._8, 40._8, 40._8, 1.0_8, 1.0_8, 1.0_8,      1.0_8,      1.0_8, 1.0_8, &
                1.0_8, 1.0_8, 5.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    frootage=(/ 18._8, 10._8, 10._8, 10._8, 5.0_8, 3.0_8, 3.0_8, 3.0_8, 0.884227_8, 0.884227_8, 1.0_8, &
                1.0_8, 1.0_8, 4.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    metage=0.04_8
    strage=0.23_8
    cwdage=0.824_8
    micage=0.137_8
    slowage=5._8
    passage=222.22_8
    clabileage=0.2_8
    slax = (/ 0.007178742_8, 0.015319264_8, 0.023103346_8, 0.026464564_8, 0.009919764_8, 0.029590494_8, &
              0.022417511_8, 0.026704118_8, 0.029590494_8, 0.022417511_8,        0.02_8,        0.02_8, &
              0.02_8,        0.024470894_8,        0.02_8,        0.02_8,        0.02_8 /)

    xfherbivore   =(/ 0.068_8, 0.406_8, 0.068_8, 0.134_8, 0.022_8, 0.109_8, 0.109_8, 0.109_8, 0.140_8, &
                      0.140_8, 0.000_8, 0.000_8, 0.000_8, 0.010_8, 0.000_8, 0.000_8, 0.000_8 /)
    xxkleafcoldmax=(/   0.2_8,  0.1_8,  0.1_8,  0.6_8,       1._8,  0.2_8,    0.2_8,   0.2_8,   0.3_8, &
                        0.3_8,  0.1_8,  0.1_8,  0.1_8,      0.1_8,  0.1_8,    0.1_8,   0.1_8 /)
    xxkleafdrymax =(/   0.1_8,  0.1_8,  0.1_8,   1._8,      0.1_8,  0.1_8,    0.1_8,   0.1_8,   0.1_8, &
                        0.1_8,  0.1_8,  0.1_8,  0.1_8,      0.1_8,  0.1_8,    0.1_8,   0.1_8 /)
    xratioNPleafmin =(/ 10.92308_8, 15.95339_8, 9.254839_8, 12.73848_8, 12.07217_8, 13.51473_8,    14.05_8, &
                        12.57800_8, 15.12262_8,      10._8,      13._8,      10._8,      10._8,  16.2336_8, &
                             10._8,      10._8,      10._8 /)
    xratioNPleafmax =(/ 12.07288_8,  17.6327_8, 10.22903_8, 14.07938_8, 13.34292_8, 14.93733_8, 15.52895_8, &
                          13.902_8, 16.71447_8,      10._8,      13._8,      10._8,      10._8,  17.9424_8, &
                             10._8,      10._8,      10._8 /)
    xratioNPwoodmin =(/ 20.30167_8, 15.89425_8, 17.48344_8, 19.08018_8, 22.46035_8,      15._8,      15._8, &
                           15.96_8,    20.52_8,      15._8,      15._8,      15._8,      15._8,  17.5275_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPwoodmax =(/ 22.43869_8, 17.56733_8,  19.3238_8, 21.08862_8,  24.8246_8,      15._8,      15._8, &
                           17.64_8,    20.52_8,      15._8,      15._8,      15._8,      15._8,  19.3725_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPfrootmin=(/ 20.29341_8, 15.87155_8, 17.39767_8,  19.0601_8, 22.49363_8, 15.63498_8, 16.08255_8, &
                        14.49241_8, 22.69109_8,      15._8,      15._8,      15._8,      15._8, 22.13268_8, &
                             15._8,      15._8,      15._8 /)
    xratioNPfrootmax=(/ 22.42955_8, 17.54224_8,   19.229_8, 21.06643_8, 24.86138_8, 17.28077_8, 17.77545_8, &
                        16.01793_8, 25.07962_8,      15._8,      15._8,      15._8,      15._8, 24.46244_8, &
                             15._8,      15._8,      15._8 /)
    xfNminloss=0.05_8
    xfNminleach=0.05_8
    xnfixrate=(/ 0.08_8, 2.6_8, 0.21_8, 1.64_8, 0.37_8, 0.95_8, 0.95_8, 0.95_8, 4._8, 4._8, 0._8, &
                   0._8,  0._8, 0.35_8,   0._8,   0._8,   0._8 /)
    xnsoilmin=1000._8
  
    ratiocnplant(:,leaf)=(/  49.8_8,  23.1_8,  59.3_8,  31.4_8,  37.6_8, 34.8_8,  44._8,  49.2_8,  21.6_8, &
                              25._8,   30._8,   30._8,   30._8,   50._8,  40._8,  40._8,   40._8 /)
    ratiocnplant(:,wood)=(/ 238.1_8, 134.9_8, 243.8_8, 156.2_8, 142.1_8, 150._8, 150._8, 147.3_8,  150._8, &
                             125._8,  150._8,  150._8,  150._8,  150._8, 150._8, 135._8,  150._8 /)
    ratiocnplant(:,xroot)=(/ 73.7_8,  61.2_8,   75._8,  63.2_8,  67.1_8, 64.5_8, 62.7_8,   69._8,  60.7_8, &
                              71._8,   71._8,   71._8,   71._8,   71._8,  71._8,  71._8,   71._8 /)
    ratiocnsoil(:,mic)=8._8
    ratiocnsoil(:,slow)=(/ 16.1_8, 12.8_8, 24.8_8,  30._8, 19.3_8, 13.1_8, 13.1_8, 13.1_8, 13.2_8, 13.2_8, &
                           13.1_8, 13.1_8, 13.1_8, 26.8_8,  20._8,  20._8,  20._8 /)
    ratiocnsoil(:,pass)=(/ 16.1_8, 12.8_8, 24.8_8,  30._8, 19.3_8, 13.1_8, 13.1_8, 13.1_8, 13.2_8, 13.2_8, &
                           13.1_8, 13.1_8, 13.1_8, 26.8_8,  20._8,  20._8,  20._8 /)
    ratiocnsoilmin(:,mic)=3._8
    ratiocnsoilmin(:,slow)=12._8
    ratiocnsoilmin(:,pass)=7._8
    ratiocnsoilmax(:,mic)=15._8
    ratiocnsoilmax(:,slow)=30._8
    ratiocnsoilmax(:,pass)=15._8
   
    xplab  =(/   26.737_8,   19.947_8,   29.107_8,   30.509_8,   23.206_8,   25.538_8,   25.538_8,   25.538_8, &
                 27.729_8,   27.729_8,       0._8,       0._8,       0._8,   21.038_8,       0._8,       0._8, &
                  0.103_8 /)
    xpsorb =(/   126.73_8,   92.263_8,  134.639_8,  132.012_8,   173.47_8,  186.207_8,  186.207_8,  186.207_8, &
                155.518_8,  155.518_8,       0._8,       0._8,       0._8,   255.79_8,       0._8,       0._8, &
                  1.176_8 /)
    xpocc  =(/  138.571_8,  120.374_8,   138.22_8,  148.083_8,  114.496_8,  145.163_8,  145.163_8,  145.163_8, &
                158.884_8,  158.884_8,       0._8,       0._8,       0._8,  108.897_8,       0._8,       0._8, &
                  0.688_8 /)

    xkmlabp  =(/ 74.5408_8,  68.1584_8,   77.952_8, 64.41918_8, 64.41918_8, 70.5856_8,  64.5888_8, 54.1692_8, &
                  9.7704_8,    28.29_8,   63.963_8,   32.402_8 /)
    xpsorbmax=(/ 745.408_8, 788.0815_8, 1110.816_8,  744.847_8,  744.847_8, 816.146_8, 746.8081_8, 722.256_8, &
                 293.112_8,   311.19_8, 373.1175_8, 615.6381_8 /)
    xfPleach =0.0005_8
    ratioNPsoil(:,mic)=4._8
    ratioNPsoil(:,slow)=(/ 5._8, 5._8, 5._8, 15._8, 5._8, 5._8, 5._8, 5._8, 7._8, 7._8, 7._8, 7._8 /)
    ratioNPsoil(:,pass)=(/ 5._8, 5._8, 5._8, 15._8, 5._8, 5._8, 5._8, 5._8, 7._8, 7._8, 7._8, 7._8 /)
  
    xxnpmax = (/ 1.510856726_8, 1.27916225_8, 1.591076159_8, 1.186066584_8, 1.358075681_8,  1.45621905_8, &
                  1.45621905_8, 1.45621905_8, 1.210382326_8, 1.210382326_8,  1.45621905_8, 1.365993164_8, &
                 1.210382326_8,         1._8, 1.399652677_8,          1._8,          1._8 /)
    xq10soil = 1.72_8
    xxkoptlitter = 0.4_8
    xxkoptsoil = (/ 0.33_8, 0.6_8, 0.15_8, 0.6_8, 0.16_8, 0.4_8, 0.3_8, 0.2_8, 0.2_8, 0.25_8,  1._8, &
                    0.65_8, 0.5_8,   2._8, 0.5_8,   1._8,  1._8 /)
    xprodptase = (/  0.5_8, 0.2_8,  0.5_8, 0.5_8,  0.5_8, 0.5_8, 0.5_8, 0.5_8, 0.5_8,  0.5_8, 0.5_8, &
                      4._8, 0.5_8,  0.5_8, 0.5_8,  0.5_8, 0.5_8 /)
    xcostnpup = (/   40._8, 25._8,  40._8, 40._8,  40._8, 40._8, 40._8, 40._8, 40._8,  40._8, 40._8, &
                     40._8, 40._8,  40._8, 40._8,  40._8, 40._8 /)
    xmaxfinelitter = (/ 1524._8, 384._8, 1527._8, 887._8, 157._8, 361._8, 225._8, 913._8, 660._8, 100._8, &
                         100._8, 100._8,  100._8,  83._8, 100._8, 100._8, 100._8 /)
    xmaxcwd = (/ 1795._8, 613._8, 1918._8, 1164._8, 107._8, 420._8, 228._8, 573._8, 811._8, 100._8, &
                  100._8, 100._8,  100._8,   23._8, 100._8, 100._8, 100._8 /)
    xnintercept = (/ 6.32_8, 4.19_8,  6.32_8,  5.73_8, 14.71_8,  6.42_8,    2._8, 14.71_8, 4.71_8, 14.71_8, &
                    14.71_8,   7._8, 14.71_8, 14.71_8, 14.71_8, 14.71_8, 14.71_8 /)
    xnslope = (/ 18.15_8, 26.19_8, 18.15_8, 29.81_8, 23.15_8, 40.96_8, 8._8, 23.15_8, 59.23_8, 23.15_8, &
                 23.15_8,   10._8, 23.15_8, 23.15_8, 23.15_8, 23.15_8, 23.15_8 /)
    xla_to_sa = (/ 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, &
                   5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8, 5000._8 /)
    xvcmax_scalar = (/ 0.92_8, 1.10_8, 0.92_8, 0.92_8, 1.25_8, 1.25_8, 1.25_8, 1.25_8, 1.0_8, 1.0_8, &
                        1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8 /)
    xdisturbance_interval = (/ 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, &
                               100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8, 100._8 /)
    xDAMM_EnzPool = (/ 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, &
                       10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8, 10.0_8 /)
    xDAMM_KMO2 = (/ 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, &
                    0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8, 0.01_8 /)
    xDAMM_KMcp = (/  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8, &
                     0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8,  0.1_8 /)
    xDAMM_Ea = (/ 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, 62._8, &
                  62._8, 62._8, 62._8, 62._8, 62._8 /)
    xDAMM_alpha = (/ 10.6_8, 10.4_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, &
                     10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8, 10.6_8 /)
  
    xxkplab = 0.001369863_8
    xxkpsorb = (/ 1.8356191E-05_8, 2.0547975E-05_8, 1.3698650E-05_8, 1.4794542E-05_8, 2.1369894E-05_8, &
                  2.3835651E-05_8, 1.9452083E-05_8, 2.1095921E-05_8, 2.7123327E-05_8, 2.1095921E-05_8, &
                  2.7123327E-05_8, 2.1095921E-05_8 /)
    xxkpocc = 2.73973E-05_8

    casabiome%ivt2     =(/        3,        3,        3,        3,        2,        1,   1,   2, &
                                  1,        1,        0,        0,        0,        1,   0,   0, &
                                  0 /)
    casabiome%kroot    =(/      5.5_8,      3.9_8,      5.5_8,      3.9_8,      2.0_8,      5.5_8, 5.5_8, 5.5_8, &
                                5.5_8,      5.5_8,      5.5_8,      5.5_8,      5.5_8,      2.0_8, 2.0_8, 5.5_8, &
                                5.5_8 /)
    casabiome%rootdepth=(/      1.5_8,      1.5_8,      1.5_8,      1.5_8,      0.5_8,      0.5_8, 0.5_8, 0.5_8, &
                                0.5_8,      0.5_8,      0.5_8,      0.5_8,      0.5_8,      0.5_8, 0.5_8, 1.5_8, &
                                0.5_8 /)
    casabiome%kuptake  =(/      2.0_8,      1.9_8,      2.0_8,      2.0_8,      1.8_8,      2.0_8, 2.0_8, 2.0_8, &
                                1.6_8,      1.6_8,      1.6_8,      1.8_8,      1.8_8,      1.8_8, 1.8_8, 1.8_8, &
                                1.8_8 /)
    casabiome%krootlen =(/ 14.87805_8, 14.38596_8, 14.02597_8, 18.94737_8, 32.30769_8,      84._8, 84._8, 84._8, &
                              120.5_8,    120.5_8,       0._8,       0._8,       0._8, 30.76923_8,  0._8,  0._8, &
                                 0._8 /)
    casabiome%kminN=2.0_8
    casabiome%kuplabP=0.5_8
    casabiome%fracnpptoP(:,leaf) =(/ 0.25_8, 0.20_8, 0.40_8, 0.35_8, 0.35_8, 0.35_8, 0.35_8, 0.50_8, 0.50_8, &
                                     0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.25_8, 0.50_8, 0.60_8, 0.50_8 /)
    casabiome%fracnpptoP(:,wood) =(/ 0.40_8, 0.35_8, 0.30_8, 0.25_8, 0.25_8, 0.00_8, 0.00_8, 0.10_8, 0.00_8, &
                                     0.00_8, 0.00_8, 0.00_8, 0.00_8, 0.25_8, 0.00_8, 0.40_8, 0.00_8 /)
    casabiome%fracnpptoP(:,xroot)=(/ 0.35_8, 0.45_8, 0.30_8, 0.40_8, 0.40_8, 0.65_8, 0.65_8, 0.40_8, 0.50_8, &
                                     0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.50_8, 0.00_8, 0.50_8 /)
    casabiome%rmplant(:,leaf)    =0.1_8
    casabiome%rmplant(:,wood)    =(/ 2.0_8, 1.0_8, 1.5_8, 0.8_8, 0.5_8, 0.5_8, 0.4_8, 1.8_8, 2.0_8, 1.0_8, &
                                     1.0_8, 1.0_8, 1.0_8, 2.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    casabiome%rmplant(:,xroot)   =(/ 10._8, 2.0_8, 7.5_8, 2.5_8, 4.5_8, 4.5_8, 4.0_8, 15._8, 25._8, 10._8, &
                                     10._8, 10._8, 10._8, 10._8, 10._8, 10._8, 10._8 /)
    casabiome%ftransNPtoL(:,leaf) =0.5_8
    casabiome%ftransNPtoL(:,wood) =0.95_8
    casabiome%ftransNPtoL(:,xroot)=0.9_8
    casabiome%fracligninplant(:,leaf) =(/ 0.25_8, 0.20_8, 0.20_8, 0.20_8, 0.20_8, 0.10_8, 0.10_8, 0.10_8, &
                                          0.10_8, 0.10_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.25_8, &
                                          0.10_8 /)
    casabiome%fracligninplant(:,wood) =0.4_8
    casabiome%fracligninplant(:,xroot)=(/ 0.25_8, 0.20_8, 0.20_8, 0.20_8, 0.20_8, 0.10_8, 0.10_8, 0.10_8, &
                                          0.10_8, 0.10_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.15_8, 0.25_8, &
                                          0.10_8 /)
    casabiome%glaimax=(/ 10._8, 10._8, 10._8, 10._8, 10._8, 3._8, 3._8, 3._8, 6._8, 6._8,  5._8, 5._8, &
                          5._8, 1._8,  6._8,   1._8,  0._8 /)
    casabiome%glaimin=(/ 1._8,  1._8, .5_8,  .5_8, .1_8, .1_8, .1_8, .1_8, .1_8, .1_8, .05_8, .05_8, &
                        .05_8, .05_8, 0._8, .05_8, 0._8 /)
    phen%TKshed=(/ 268._8,   260._8, 263.15_8, 268.15_8, 277.15_8, 275.15_8, 275.15_8, 275.15_8, 278.15_8, &
                 278.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 277.15_8, 283.15_8 /)
    casabiome%xkleafcoldexp=3._8
    casabiome%xkleafdryexp=3._8
    casabiome%ratioNCplantmin(:,leaf) =(/     0.02_8,     0.04_8, 0.016667_8, 0.028571_8,    0.025_8,  0.02631_8, &
                                              0.02_8,     0.02_8,     0.04_8,     0.04_8, 0.033333_8,    0.025_8, &
                                             0.025_8, 0.018182_8,    0.025_8,    0.025_8,    0.025_8 /)
    casabiome%ratioNCplantmax(:,leaf) =(/    0.024_8,    0.048_8,     0.02_8, 0.034286_8,     0.03_8, 0.031572_8, &
                                             0.024_8,    0.024_8,    0.048_8,    0.048_8,     0.04_8,     0.03_8, &
                                              0.03_8, 0.022222_8,     0.03_8,     0.03_8,     0.03_8 /)
    casabiome%ratioNCplantmin(:,wood) =(/    0.004_8, 0.006667_8,    0.004_8, 0.005714_8, 0.006667_8, 0.006667_8, &
                                          0.006667_8, 0.006667_8,    0.008_8,    0.008_8, 0.006667_8, 0.006667_8, &
                                          0.006667_8, 0.006667_8, 0.006667_8, 0.007307_8, 0.006667_8 /)
    casabiome%ratioNCplantmax(:,wood) =(/   0.0048_8,    0.008_8,   0.0048_8, 0.006857_8,    0.008_8,    0.008_8, &
                                             0.008_8,    0.008_8,   0.0096_8,   0.0096_8,    0.008_8,    0.008_8, &
                                             0.008_8,    0.008_8,    0.008_8, 0.008889_8,    0.008_8 /)
    casabiome%ratioNCplantmin(:,xroot)=(/ 0.012821_8, 0.014706_8, 0.012821_8, 0.014085_8, 0.014085_8, 0.014085_8, &
                                          0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, &
                                          0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8, 0.014085_8 /)
    casabiome%ratioNCplantmax(:,xroot)=(/ 0.015385_8, 0.017647_8, 0.015385_8, 0.016901_8, 0.016901_8, 0.016901_8, &
                                          0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, &
                                          0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8, 0.016901_8 /)
    casabiome%ftransPPtoL(:,leaf)=0.5_8
    casabiome%ftransPPtoL(:,wood)=0.95_8
    casabiome%ftransPPtoL(:,xroot)=0.9_8

  end if

end if
    
if ( mp_global>0 ) then
  casabiome%ratioPcplantmin(:,leaf)  = 1._8/(xratioNPleafmax*ratioCNplant(:,leaf))
  casabiome%ratioPcplantmax(:,leaf)  = 1._8/(xratioNPleafmin*ratioCNplant(:,leaf))
  casabiome%ratioPcplantmin(:,wood)  = 1._8/(xratioNPwoodmax*ratioCNplant(:,wood))
  casabiome%ratioPcplantmax(:,wood)  = 1._8/(xratioNPwoodmin*ratioCNplant(:,wood))
  casabiome%ratioPcplantmin(:,xroot) = 1._8/(xratioNPfrootmax*ratioCNplant(:,xroot))
  casabiome%ratioPcplantmax(:,xroot) = 1._8/(xratioNPfrootmin*ratioCNplant(:,xroot))

  casabiome%ratioNPplantmin(:,leaf)  = xratioNPleafmin
  casabiome%ratioNPplantmax(:,leaf)  = xratioNPleafmax
  casabiome%ratioNPplantmin(:,wood)  = xratioNPwoodmin
  casabiome%ratioNPplantmax(:,wood)  = xratioNPwoodmax
  casabiome%ratioNPplantmin(:,xroot) = xratioNPfrootmin
  casabiome%ratioNPplantmax(:,xroot) = xratioNPfrootmax    

  casabiome%sla                = slax
  casabiome%fraclabile(:,leaf) = deltcasa*0.6_8    !1/day
  casabiome%fraclabile(:,xroot)= deltcasa*0.4_8    !1/day
  casabiome%fraclabile(:,wood) = 0._8
  casabiome%plantrate(:,leaf)  = deltcasa/(leafage*(1._8-xfherbivore))
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
  casabiome%rmplant            = deltcasa*casabiome%rmplant
  casabiome%kclabrate          = deltcasa/clabileage

  casabiome%xnpmax(:)          = xxnpmax(:)
  casabiome%q10soil(:)         = xq10soil(:)
  casabiome%xkoptlitter(:)     = xxkoptlitter(:)
  casabiome%xkoptsoil(:)       = xxkoptsoil(:)
  casabiome%prodptase(:)       = xprodptase(:)/365._8   ! convert from yearly to daily
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
    casapool%cplant(:,wood)  = 0._8
    casapool%clitter(:,cwd)  = 0._8
    casapool%nplant(:,wood)  = 0._8
    casapool%nlitter(:,cwd)  = 0._8
    casapool%pplant(:,wood)  = 0._8
    casapool%plitter(:,cwd)  = 0._8
  end where
  if ( cable_pop==1 ) then
   where (casamet%iveg2==forest.or.casamet%iveg2==shrub)
      casapool%cplant(:,wood)  = 0.01_8
      casapool%nplant(:,wood)  = casabiome%ratioNCplantmin(veg%iveg,wood)*casapool%cplant(:,wood)
      casapool%pplant(:,wood)  = casabiome%ratioPCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
    end where
  end if
  casapool%cplant(:,leaf)     = cleaf(veg%iveg)
  casapool%cplant(:,xroot)    = cfroot(veg%iveg)
  casapool%clabile            = 0._8
  casapool%clitter(:,metb)    = cmet(veg%iveg)
  casapool%clitter(:,str)     = cstr(veg%iveg)
  casapool%csoil(:,mic)       = cmic(veg%iveg)
  casapool%csoil(:,slow)      = cslow(veg%iveg)
  casapool%csoil(:,pass)      = cpass(veg%iveg)
  if ( ccycle==1 ) then
    casapool%ratioNCplant     = 1._8/ratioCNplant(veg%iveg,:)  
  end if
  casapool%dclabiledt         = 0._8

  ! initializing glai in case not reading pool file (eg. during spin)
  casamet%glai = max(casabiome%glaimin(veg%iveg), casabiome%sla(veg%iveg)*casapool%cplant(:,leaf))
  casaflux%fNminloss   = xfNminloss(veg%iveg)
  casaflux%fNminleach  = 10._8*xfNminleach(veg%iveg)*deltcasa
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
  casaflux%fpleach        = xfPleach(casamet%isorder)/365._8

  casapool%ratioNCplant   = 1._8/ratioCNplant(veg%iveg,:)
  casapool%ratioNPplant   = casabiome%ratioNPplantmin(veg%iveg,:)
  casapool%ratioNClitter  = casapool%nlitter/(casapool%clitter+1.0e-10_8)
  casapool%ratioNPlitter  = casapool%nlitter/(casapool%plitter+1.0e-10_8)
  casapool%ratioNCsoil    = 1._8/ratioCNsoil(veg%iveg,:)
  casapool%ratioNPsoil    = ratioNPsoil(casamet%isorder,:)
  casapool%ratioNCsoilmin = 1._8/ratioCNsoilmax(veg%iveg,:)
  casapool%ratioNCsoilmax = 1._8/ratioCNsoilmin(veg%iveg,:)
  casapool%ratioNCsoilnew = casapool%ratioNCsoilmax

  casapool%ratioPCplant   = casabiome%ratioPcplantmax(veg%iveg,:)
  casapool%ratioPClitter  = casapool%plitter/(casapool%clitter(:,:)+1.0e-10_8)
  casapool%ratioPCsoil    = 1._8/(ratioCNsoil(veg%iveg,:)*ratioNPsoil(casamet%isorder,:))

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
    write(6,*) "-> Using default CABLE biophysical parameter tables"
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
    write(6,*) "-> Using user defined CABLE biophysical parameter tables"
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
        write(6,*) "-> WARN: Replacing g0=0. with g0=0.01 for gs_switch=1"
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
use parm_m, only : nmaxpr
use soilv_m

type(soil_parameter_type), intent(inout) :: soil
!real, dimension(mp_global) :: ssat_bounded, rho_soil_bulk
!real, parameter :: ssat_hi = 0.65
!real, parameter :: ssat_lo = 0.15
!real, parameter :: rhob_hi = 2300.
!real, parameter :: rhob_lo = 810.
integer isoil, k
integer, dimension(1) :: nstart, ncount

if ( lncveg_numsoil<1 ) then
    
  ! default soil parameter tables
  if ( myid==0 ) then
    write(6,*) "-> Using default CABLE soil parameter tables"
  end if

  ! redefine rhos
  rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

  if ( myid==0 .and. nmaxpr==1 ) then
    do isoil = 1,mxst
      write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
    end do
  end if

else

  ! user defined soil parameter tables  
  if ( myid==0 ) then
    write(6,*) "-> Using user defined CABLE soil parameter tables"
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
    call ccnf_get_vara(ncidveg,'silt',nstart,ncount,silt(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'clay',nstart,ncount,clay(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sand',nstart,ncount,sand(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'swilt',nstart,ncount,swilt(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sfc',nstart,ncount,sfc(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'ssat',nstart,ncount,ssat(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'bch',nstart,ncount,bch(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'hyds',nstart,ncount,hyds(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'sucs',nstart,ncount,sucs(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'rhosoil',nstart,ncount,rhos(1:lncveg_numsoil))
    call ccnf_get_vara(ncidveg,'css',nstart,ncount,css(1:lncveg_numsoil))
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
    if ( myid==0 .and. nmaxpr==1 ) then
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
    soil%sucs_vec(:,k)    = 1000._8*abs(soil%sucs)
    soil%bch_vec(:,k)     = soil%bch
    soil%hyds_vec(:,k)    = 1000._8*soil%hyds
    soil%watr(:,k)        = 0.05_8
    soil%cnsd_vec(:,k)    = soil%cnsd
    soil%clay_vec(:,k)    = soil%clay
    soil%sand_vec(:,k)    = soil%sand
    soil%silt_vec(:,k)    = soil%silt
    soil%zse_vec(:,k)     = soil%zse(k)
    soil%css_vec(:,k)     = soil%css
  end do
  soil%GWhyds_vec    = soil%hyds*1000._8  
  soil%GWsucs_vec    = abs(soil%sucs)*1000._8
  soil%GWbch_vec     = soil%bch
  soil%GWrhosoil_vec = soil%rhosoil
  soil%GWssat_vec    = soil%ssat
  soil%GWwatr        = soil%watr(:,ms) !residual water content of the aquifer [mm3/mm3]
  soil%GWdz          = 20._8           !thickness of the aquifer   [m]
  soil%GWdz          = max( 1._8, min( 20._8, soil%GWdz - sum(soil%zse,dim=1) ) )
  soil%drain_dens    = 0.008_8         !  drainage density ( mean dist to rivers/streams )
  
  soil%heat_cap_lower_limit = 0.01 ! recalculated in cable_soilsnow.F90
  
  !if ( cable_gw_model==1 ) then
  !
  !  ! note 17 hard coded vegetation PFTs  
  !  psi_o(1:3)  = -66000._r_2
  !  psi_o(4)    = -35000._r_2
  !  psi_o(5)    = -83000._r_2
  !  psi_o(6:17) = -74000._r_2
  !  psi_c(1:3)  = -2550000._r_2
  !  psi_c(4)    = -2240000._r_2
  !  psi_c(5)    = -4280000._r_2
  !  psi_c(6:17) = -2750000._r_2      
  !    
  !  soil_depth(1) = REAL(soil%zse(1),r_2)
  !  DO k=2,ms
  !     soil_depth(k) = soil_depth(k-1) + REAL(soil%zse(k),r_2)
  !  END DO    
  !  
  !  do k=1,ms
  !    soil%hyds_vec(:,k) = 0.0070556_r_2*10.0_r_2**(-0.884_r_2 + 0.0153_r_2*soil%sand_Vec(:,k)*100.0_r_2)* &
  !      EXP(-gw_params%hkrz*(MAX(0._r_2,soil_depth(k)-gw_params%zdepth)))
  !    soil%sucs_vec(:,k) = 10.0_r_2 * 10.0_r_2**(1.88_r_2 -0.0131_r_2*soil%Sand_Vec(:,k)*100.0_r_2)
  !    soil%bch_vec(:,k) = 2.91_r_2 + 0.159_r_2*soil%Clay_Vec(:,k)*100.0_r_2
  !    soil%ssat_vec(:,k) = 0.489_r_2 - 0.00126_r_2*soil%Sand_Vec(:,k)*100.0_r_2
  !    soil%watr(:,k) = 0.02_r_2 + 0.00018_r_2*soil%Clay_Vec(:,k)*100.0_r_2
  !  end do
  !  !aquifer share non-organic with last layer if not found in param file
  !  soil%GWhyds_vec(:) = soil%hyds_vec(:,ms)
  !  soil%GWsucs_vec(:) = soil%sucs_vec(:,ms)
  !  soil%GWbch_vec(:)  = soil%bch_vec(:,ms)
  !  soil%GWssat_vec(:) = soil%ssat_vec(:,ms)
  !  soil%GWwatr(:)     = soil%watr(:,ms)
  !  !include organin impact.  fraction of grid cell where percolation through
  !  !organic macropores dominates
  !  !soil%Org_Vec = MAX(0._r_2,soil%Org_Vec)
  !  !soil%Org_Vec = MIN(1._r_2,soil%Org_Vec)
  !  !do k=1,3  !0-23.3 cm, data really is to 30cm
  !  !  soil%hyds_vec(:,k)  = (1.-soil%Org_Vec(:,k))*soil%hyds_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%hyds_vec_organic
  !  !  soil%sucs_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%sucs_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%sucs_vec_organic
  !  !  soil%bch_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%bch_vec(:,k) +&
  !  !    soil%Org_Vec(:,k)*gw_params%org%clappb_organic
  !  !  soil%ssat_vec(:,k) = (1.-soil%Org_Vec(:,k))*soil%ssat_vec(:,k) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%ssat_vec_organic
  !  !  soil%watr(:,k)   = (1.-soil%Org_Vec(:,klev))*soil%watr(:,klev) + &
  !  !    soil%Org_Vec(:,k)*gw_params%org%watr_organic
  !  !end do
  !
  !  !vegetation dependent field capacity (point plants get stressed) and
  !  !wilting point
  !  do k = 1,ms
  !    psi_tmp(:,k) = -psi_c(veg%iveg(:))
  ! end do
  ! soil%sfc_vec = (soil%ssat_vec-soil%watr) * (ABS(psi_tmp(:,:)) &
  !    /(ABS(soil%sucs_vec)))**(-1.0/soil%bch_vec)+soil%watr
  !  do k = 1,ms
  !    psi_tmp(:,k) = -psi_c(veg%iveg(:))
  !  end do
  !  soil%swilt_vec = (soil%ssat_vec-soil%watr) * (ABS(psi_tmp(:,:)) &
  !    /(ABS(soil%sucs_vec)))**(-1.0/soil%bch_vec)+soil%watr
  !
  !  !set the non-vectored values to srf value
  !  soil%sfc(:) = REAL(soil%sfc_vec(:,1))
  !  soil%swilt(:) = REAL(soil%swilt_vec(:,1))
  !
  !  !convert the units back to what default uses and GW only uses the
  !  !vectored versions
  !  soil%hyds = REAL(soil%hyds_vec(:,1))/1000.0
  !  soil%sucs = REAL(soil%sucs_vec(:,1))/1000.0
  !  soil%ssat = REAL(soil%ssat_vec(:,1))
  !  soil%bch  = REAL(soil%bch_vec(:,1))  
  !  
  !end if 
  
  !if ( cable_user%gw_model ) then
  !  do k = 1,ms
  !    soil%hyds_vec(:,k) = soil%hyds_vec(:,k) * exp(-gw_parms%hkrz*(znode(k)-gw_params%zdepth) )
  !  end do
  !  soil%hyds(:) = soil%hyds_vec(:,1)
  !end if
  
  !if ( cable_user%soil_thermal_fix ) then
  !  do k = 1,ms
  !    ssat_bounded = min( ssat_hi, max( ssat_lo, soil%ssat_vec(:,k) ) )  
  !    rho_soil_bulk = min( rhob_hi, max( rhob_lo, (2700.*(1.-ssat_bounded)) ) )
  !    where ( soil%isoilm(:) /= 9 )
  !     soil%rhosoil_vec(:,k) = 2700.
  !      soil%cnsd_vec(:,k) = ( 0.135*(1.-ssat_bounded(:)) + 64.7/rho_soil_bulk(:) ) &
  !                         / ( 1.-0.947*(1.-ssat_bounded(:)) )
  !    end where    
  !  end do
  !  where ( soil%isoilm(:) /= 9 )
  !    soil%rhosoil(:) = soil%rhosoil_vec(:,1)
  !    soil%cnsd(:) = soil%cnsd_vec(:,1)
  !  end where  
  !end if
  
end if

return
end subroutine cable_soil_parm
                   
! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlin,fveg,cableformat)
  
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parmgeom_m

character(len=*), intent(in) :: fveg
integer, dimension(ifull,maxtile), intent(out) :: ivs
integer, dimension(ifull_g,maxtile) :: ivsg  
integer, dimension(3) :: spos,npos
integer n,iq,ilx,jlx,iad 
integer ncidx,iernc,varid,ndims
real, dimension(ifull,maxtile), intent(out) :: svs, vlin
real, dimension(ifull_g,maxtile) :: svsg, vling
real, dimension(ifull_g) :: savannafrac_g
real rlong0x,rlat0x,schmidtx,dsx,ra,rb
real cablever
real, intent(out) :: cableformat
character(len=47) header  
character(len=7) vname
logical tst

cableformat=0.

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
  if (abs(cablever-cable_version)>1.e-20 .and. abs(cablever-6608.)>1.e-20 .and. &
      abs(cablever-3939.)>1.e-20 ) then
    write(6,*) "Wrong version of CABLE data"
    write(6,*) "Expecting 3939. or 6608. or ",cable_version
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
end if

call ccmpi_bcast(cableformat,0,comm_world)

return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlin,cableformat)
  
use cc_mpi
use newmpar_m
  
integer, dimension(ifull,maxtile), intent(out) :: ivs
real, dimension(ifull,maxtile), intent(out) :: svs, vlin
real, intent(out) :: cableformat

cableformat = 0.

call ccmpi_distribute(ivs)
call ccmpi_distribute(svs)
call ccmpi_distribute(vlin)
  
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
  ssnow%GWwb = 0.5_8*soil%ssat 
  ssnow%wtd = 20000._8
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

subroutine newcbmwb

use soilsnow_m

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
use nsibd_m, only : sigmf, carb_plant, carb_litter, carb_soil
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
logical, intent(in), optional :: usedefault
integer k, n, ierr, idv, ierr_casa, ierr_sli, ierr_pop, ierr_svs, ierr_cvc
integer jyear,jmonth,jday,jhour,jmin,mins, ll, cc, hh, dd
integer np_pop, iq, m
integer, dimension(6) :: ierr_check
integer, dimension(ifull) :: dati
integer, dimension(mp_global) :: old_cv
integer, dimension(ifull,maxtile) :: nmp
integer, dimension(:), allocatable :: dati_out
real, dimension(ifull) :: datr
real, dimension(mp_global) :: dummy_unpack, old_sv
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(ifull,ms) :: datms
real(kind=8), dimension(ifull,3) :: dat3
real(kind=8), dimension(ifull,mplant) :: datmplant
real(kind=8), dimension(ifull,mlitter) :: datmlitter
real(kind=8), dimension(ifull,msoil) :: datmsoil
real(kind=8), dimension(:), allocatable :: dat_out
real(kind=8), dimension(:,:), allocatable :: datpatch
real(kind=8), dimension(:,:), allocatable :: datage
real(kind=8), dimension(:,:,:), allocatable :: datpc
logical tst
logical defaultmode
character(len=80) vname
character(len=21) testname

if ( myid==0 ) write(6,*) 'Read CABLE and CASA initial conditions'

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
ierr_cvc = 1

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
    write(testname,'("t",I1.1,"_cvc")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_cvc = 0
    end if
  end if
end if

! Communicate with processes if not all processes are reading the input file.
if ( .not.pfall ) then
  ierr_check(1) = ierr
  ierr_check(2) = ierr_casa
  ierr_check(3) = ierr_sli
  ierr_check(4) = ierr_pop
  ierr_check(5) = ierr_svs
  ierr_check(6) = ierr_cvc
  call ccmpi_bcast(ierr_check(1:6),0,comm_world)
  ierr       = ierr_check(1)
  ierr_casa  = ierr_check(2)
  ierr_sli   = ierr_check(3)
  ierr_pop   = ierr_check(4)
  ierr_svs   = ierr_check(5)
  ierr_cvc   = ierr_check(6)
end if

if ( myid==0 ) then
  write(6,*) "-> Found ierr,ierr_casa_ierr_sli ",ierr,ierr_casa,ierr_sli
  write(6,*) "->    ierr_pop,ierr_svs,ierr_cvc ",ierr_pop,ierr_svs,ierr_cvc
end if
  
call defaulttile ! initially use default values before overwriting

! default
old_sv = sv
old_cv = cveg

! check for changes
if ( ierr_cvc==0 ) then
  do n = 1,maxtile      
    write(vname,'("t",I1.1,"_cvc")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)  
    call cable_pack(dati,old_cv,n)
  end do
end if
call create_new_tile_map(old_cv,nmp)

if ( ierr/=0 ) then
  ! Cannot locate tile data, use diagnostic data instead    
  if ( myid==0 ) write(6,*) "-> Use gridbox averaged data to initialise CABLE"
else
  ! read tile data
  if ( myid==0 ) write(6,*) "-> Use tiled data to initialise CABLE"  
  do n = 1,maxtile
    if ( ierr_svs == 0 ) then
      write(vname,'("t",I1.1,"_svs")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      datr = real( dat )  
      call cable_pack(datr,old_sv,nmp(:,n))
    end if        
    write(vname,'("t",I1.1,"_tgg")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      do iq = 1,ifull
        if ( land(iq) .and. (datms(iq,k)<100._8.or.datms(iq,k)>400._8) ) then
          ! change in land-sea mask?
          datms(iq,k) = tss(iq) ! use surface temperature
        end if
      end do  
      call cable_pack(datms(:,k),ssnow%tgg(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_wb")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      call cable_pack(datms(:,k),ssnow%wb(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_wbice")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      call cable_pack(datms(:,k),ssnow%wbice(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_tggsn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%tggsn(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_smass")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%smass(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_ssdn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%ssdn(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%sdepth(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_sconds")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%sconds(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%ssdnn(:),nmp(:,n))
    write(vname,'("t",I1.1,"_sflag")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)
    call cable_pack(dati,ssnow%isflag(:),nmp(:,n))
    write(vname,'("t",I1.1,"_snd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%snowd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_osnd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%osnowd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_snage")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%snage(:),nmp(:,n))
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%rtsoil(:),nmp(:,n))
    write(vname,'("t",I1.1,"_GWwb")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%GWwb(:),nmp(:,n))
    write(vname,'("t",I1.1,"_wtd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%wtd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_cansto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%cansto(:),nmp(:,n))
    write(vname,'("t",I1.1,"_us")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%us(:),nmp(:,n))
    write(vname,'("t",I1.1,"_pudsto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%pudsto(:),nmp(:,n))
    write(vname,'("t",I1.1,"_wetfac")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%wetfac(:),nmp(:,n))
    write(vname,'("t",I1.1,"_ga")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%ga(:),nmp(:,n))
  end do
  
  ! soil temperature check
  if ( mp_global>0 ) then
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid CABLE temperature when reading tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      stop -1
    end if
  end if 

end if ! ierr/=0 ..else..
  
if ( soil_struc==1 ) then
  if ( ierr_sli/=0 ) then
    if ( myid==0 ) write(6,*) "-> Use gridbox averaged data to initialise SLI"
  else 
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise SLI"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_hzero")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%h0(:),nmp(:,n))
      write(vname,'("t",I1.1,"_s")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%S(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_tsoil")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%tsoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_thetai")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%thetai(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%snowliq(:,1),nmp(:,n)) ! currently nsnow_max=1
      write(vname,'("t",I1.1,"_tsurface")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%tsurface(:),nmp(:,n))
      write(vname,'("t",I1.1,"_nsnow")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,ssnow%nsnow(:),nmp(:,n))
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%fwsoil(:),nmp(:,n))
    end do  
  end if ! ierr_sli/=0 ..else..
end if   ! soil_struc==1
  
if ( ccycle/=0 ) then
  if ( ierr_casa/=0 ) then
    if ( .not.allocated( carb_plant ) ) then
      if ( myid==0 ) write(6,*) "-> Use default data to initialise CASA-CNP"  
    else  
      if ( myid==0 ) write(6,*) "-> Use interpolated tiled data to initialise CASA-CNP"
      do m = 1,10
        call loadtile_carbonpools(carb_plant(:,:,1,m),casapool%cplant(:,:),m)
        call loadtile_carbonpools(carb_plant(:,:,2,m),casapool%nplant(:,:),m)
        call loadtile_carbonpools(carb_plant(:,:,3,m),casapool%pplant(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,1,m),casapool%clitter(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,2,m),casapool%nlitter(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,3,m),casapool%plitter(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,1,m),casapool%csoil(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,2,m),casapool%nsoil(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,3,m),casapool%psoil(:,:),m)
      end do   ! mm = 1,10
      m = 14 ! use index 11 to store pft 14
      call loadtile_carbonpools(carb_plant(:,:,1,11),casapool%cplant(:,:),m)
      call loadtile_carbonpools(carb_plant(:,:,2,11),casapool%nplant(:,:),m)
      call loadtile_carbonpools(carb_plant(:,:,3,11),casapool%pplant(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,1,11),casapool%clitter(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,2,11),casapool%nlitter(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,3,11),casapool%plitter(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,1,11),casapool%csoil(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,2,11),casapool%nsoil(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,3,11),casapool%psoil(:,:),m)
      deallocate( carb_plant, carb_litter, carb_soil )
    end if  
  else
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise CASA-CNP"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_cplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%cplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%nplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_pplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%pplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_clitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%clitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nlitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%nlitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_plitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%plitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_csoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%csoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nsoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%nsoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_psoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%psoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_glai")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casamet%glai,nmp(:,n))
      write(vname,'("t",I1.1,"_phen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,phen%phen,nmp(:,n))
      write(vname,'("t",I1.1,"_aphen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,phen%aphen,nmp(:,n))
      write(vname,'("t",I1.1,"_phenphase")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,phen%phase,nmp(:,n))
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,phen%doyphase(:,3),nmp(:,n))
      write(vname,'("t",I1.1,"_clabile")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%clabile,nmp(:,n))
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%nsoilmin,nmp(:,n))
      write(vname,'("t",I1.1,"_psoillab")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoillab,nmp(:,n))
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoilsorb,nmp(:,n))
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoilocc,nmp(:,n))
      write(vname,'("t",I1.1,"_crmplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casaflux%crmplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%frac_sapwood,nmp(:,n))
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%sapwood_area,nmp(:,n))
      write(vname,'("t",I1.1,"_crsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%crsoil,nmp(:,n))
      write(vname,'("t",I1.1,"_cnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%cnpp,nmp(:,n))
      write(vname,'("t",I1.1,"_clabloss")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%clabloss,nmp(:,n))
      write(vname,'("t",I1.1,"_crgplant")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%crgplant,nmp(:,n))
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%stemnpp,nmp(:,n))
      write(vname,'("t",I1.1,"_LAImax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%laimax,nmp(:,n))
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%cleafmean,nmp(:,n))
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%crootmean,nmp(:,n))
      write(vname,'("t",I1.1,"_fpn")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%fpn,nmp(:,n))
      write(vname,'("t",I1.1,"_frday")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%frday,nmp(:,n))
    end do  
  end if ! ierr_casa/=0 ..else..
end if   ! ccycle/=0

if ( cable_pop==1 ) then
  if ( ierr_pop/=0 ) then
    if ( myid==0 ) write(6,*) "-> Use default data to initialise POP"
  else
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise POP"    
    allocate( datpatch(ifull,POP_NPATCH) )  
    allocate( datage(ifull,POP_AGEMAX) )  
    allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    np_pop = size(pop%pop_grid)
    allocate( dat_out(np_pop), dati_out(np_pop) )
    dat_out = 0._8
    dati_out = 0
    do n = 1,maxtile  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cmass_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cmass_sum_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cheartwood_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%csapwood_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%csapwood_sum_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%densindiv = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%height_mean = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%height_max = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%basal_area = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_loss = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area_loss = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%stress_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crowding_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%fire_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cat_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%res_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%growth = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%area_growth = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_cover = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_area = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_volume = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%KClump = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        call pop_pack(datage(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq_age(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        call pop_pack(datage(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%biomass_age(k) = dat_out(1:np_pop)
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%biomass(ll) = dat_out(1:np_pop)
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%density(ll) = dat_out(1:np_pop)
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%hmean(ll) = dat_out(1:np_pop)
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%hmax(ll) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n)) 
        pop%pop_grid(:)%cmass_stem_bin(hh) = dat_out(1:np_pop)
      end do  
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%densindiv_bin(hh) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%height_bin(hh) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%diameter_bin(hh) = dat_out(1:np_pop)
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)  
        call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%n_age(dd) = dati_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        dati = nint(datpatch(:,k))  
        call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%id = dati_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq_old(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%factor_recruit = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%pgap = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%lai = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%biomass = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%biomass_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%heartwood = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%stress_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%fire_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%cat_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%crowding_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%cpc = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_loss = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area_loss = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%growth = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%area_growth = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_NPP = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_respiration = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_light_uptake = dat_out(1:np_pop)
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%disturbance_interval(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))   
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%first_disturbance_year(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k)) 
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%age(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%ranked_age_unique(k,dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%freq_ranked_age_unique(k,dd) = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%ncohort = dati_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%biomass = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%density = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%hmean = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%hmax = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT                  
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age = dati_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id = dati_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake = dat_out(1:np_pop)
          end do  
        end do
      end do        
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar = dat_out(1:np_pop)
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height = dat_out(1:np_pop)
          end do  
        end do
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area = dat_out(1:np_pop)
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot = dat_out(1:np_pop)
          end do  
        end do
      end do 
    end do
    deallocate( datpatch )
    deallocate( datage )
    deallocate( datpc )
    deallocate( dat_out, dati_out )
  end if ! ierr_pop/=0
end if   ! cable_pop==1  
  
if ( ierr==0 ) then
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
end if ! ierr==0  
  
call redistribute_tile(old_sv)
call fixtile

! Calculate LAI and veg fraction diagnostics
vlai(:) = 0.
sigmf(:) = 0.
if ( mp_global>0 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jmonth,jday,jhour,jmin,mp_global,sv,vl2,casamet,veg,ifull)
  dummy_unpack = sv*real(veg%vlai)
  call cable_unpack(dummy_unpack,vlai)
end if

return
end subroutine loadtile

subroutine loadtile_carbonpools(dat_in,dat_out,mm)

use newmpar_m

implicit none

integer, intent(in) :: mm
integer nb, k, msize
real, dimension(:,:), intent(in) :: dat_in
real(kind=8), dimension(:,:), intent(inout) :: dat_out
real(kind=8), dimension(size(dat_out,1)) :: dummy_unpack

msize = size(dat_in,2)

do k = 1,msize
  do nb = 1,maxnb
    dummy_unpack(:) = 0._8  
    call cable_pack(dat_in(:,k),dummy_unpack(:),nb)
  end do
  where ( veg%iveg(:)==mm .and. dummy_unpack(:)>1.e-8_8 )
    dat_out(:,k) = dummy_unpack(:)
  end where
end do

return
end subroutine loadtile_carbonpools

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

  if ( ccycle > 0 ) then
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
    
    where ( .not.( casamet%iveg2==forest.or.casamet%iveg2==shrub ) )
      casapool%cplant(:,wood)  = 0._8
      casapool%clitter(:,cwd)  = 0._8
      casapool%nplant(:,wood)  = 0._8
      casapool%nlitter(:,cwd)  = 0._8
      casapool%pplant(:,wood)  = 0._8
      casapool%plitter(:,cwd)  = 0._8
    end where

    ! initializing glai in case not reading pool file (eg. during spin)
    casapool%ratioNClitter = casapool%nlitter/(casapool%clitter+1.0e-10_8)
    casapool%ratioNPlitter = casapool%nlitter/(casapool%plitter+1.0e-10_8)
    casapool%ratioPClitter = casapool%plitter/(casapool%clitter+1.0e-10_8)

    if ( ccycle<2 ) then
      casapool%Nplant = casapool%Cplant*casapool%ratioNCplant
      casapool%Nsoil  = casapool%ratioNCsoil*casapool%Csoil
    else if ( ccycle<3 ) then
      casapool%Psoil  = casapool%Nsoil/casapool%ratioNPsoil
    end if
    
    where ( veg%iveg==6 .or. veg%iveg==7 .or. veg%iveg==8 .or. veg%iveg==9 .or. veg%iveg==10 )
      casapool%cplant(:,wood) = 0._8
      casapool%clitter(:,cwd) = 0._8
      casapool%nplant(:,wood) = 0._8
      casapool%nlitter(:,cwd) = 0._8  
      casapool%pplant(:,wood) = 0._8
      casapool%plitter(:,cwd) = 0._8   
    end where
    where ( veg%iveg==11 .or. veg%iveg==12 .or. veg%iveg==13 .or. veg%iveg==15 .or. veg%iveg==16 .or. &
            veg%iveg==17 )
      casapool%cplant(:,leaf) = 0._8
      casapool%cplant(:,wood) = 0._8
      casapool%cplant(:,froot) = 0._8
      casapool%clitter(:,metb) = 0._8
      casapool%clitter(:,str) = 0._8
      casapool%clitter(:,cwd) = 0._8
      casapool%csoil(:,mic) = 0._8
      casapool%csoil(:,slow) = 0._8
      casapool%csoil(:,pass) = 0._8
      casapool%nplant(:,leaf) = 0._8
      casapool%nplant(:,wood) = 0._8
      casapool%nplant(:,froot) = 0._8
      casapool%nlitter(:,metb) = 0._8
      casapool%nlitter(:,str) = 0._8
      casapool%nlitter(:,cwd) = 0._8
      casapool%nsoil(:,mic) = 0._8
      casapool%nsoil(:,slow) = 0._8
      casapool%nsoil(:,pass) = 0._8
      casapool%pplant(:,leaf) = 0._8
      casapool%pplant(:,wood) = 0._8
      casapool%pplant(:,froot) = 0._8
      casapool%plitter(:,metb) = 0._8
      casapool%plitter(:,str) = 0._8
      casapool%plitter(:,cwd) = 0._8
      casapool%psoil(:,mic) = 0._8
      casapool%psoil(:,slow) = 0._8
      casapool%psoil(:,pass) = 0._8
    end where
    where ( veg%iveg==14 )        
      casapool%cplant(:,leaf) = 0._8
      casapool%cplant(:,wood) = 0._8
      casapool%cplant(:,froot) = 0._8
      casapool%pplant(:,leaf) = 0._8
      casapool%pplant(:,wood) = 0._8
      casapool%pplant(:,froot) = 0._8
    end where
    
  end if ! ccycle>0
end if   ! mp_global>0
  
return
end subroutine fixtile

! if the vegetation class has changed for a tile, then attempt to find the tile associated
! with the equivilent vegetation class in the restart data
subroutine create_new_tile_map(old_cv,nmp)

use cc_mpi
use newmpar_m, only : ifull

integer n, m, iq, store_n
integer, dimension(mp_global), intent(in) :: old_cv
integer, dimension(ifull,maxtile), intent(inout) :: nmp
integer, dimension(ifull,maxtile) :: oldv_up, newv_up
integer, dimension(maxtile) :: oldv_v, newv_v
logical :: cv_test

if ( myid==0 ) write(6,*) "-> Create map to vegetation types"
if ( mp_global>0 ) then
  cv_test = any( cveg/=old_cv )
else
  cv_test = .false.
end if

oldv_up = 0
newv_up = 0

! default map  
do iq = 1,ifull
  do n = 1,maxtile  
    nmp(iq,n) = n  
  end do  
end do  

if ( mp_global>0 ) then
  if ( cv_test ) then
    
    do n = 1,maxtile
      call cable_unpack(cveg,newv_up(:,n),n)
      call cable_unpack(old_cv,oldv_up(:,n),n)
    end do     
    
    do iq = 1,ifull
      newv_v = newv_up(iq,:) ! current vegetation class
      oldv_v = oldv_up(iq,:) ! vegetation class in restart file
      ! check each tile for mismatch
      do n = 1,maxtile
        if ( oldv_v(n)/=newv_v(n) ) then
          ! search over all tiles for replacement
          do m = 1,maxtile
            if ( oldv_v(m)==newv_v(n) ) then
              ! shuffle tiles to keep area fraction sum equal to 1.  
              store_n   = nmp(iq,m)
              nmp(iq,m) = nmp(iq,n) ! found required tile
              nmp(iq,n) = store_n
              store_n   = oldv_v(m)
              oldv_v(m) = oldv_v(n)
              oldv_v(n) = store_n !=newv_v(n)   
              exit
            end if  
          end do
        end if  
      end do
    end do
    
  end if
end if

return
end subroutine create_new_tile_map


! redistribute temperature and water with a gridbox if tile area fraction changes
subroutine redistribute_tile(old_sv)

use cc_mpi

integer k
real, dimension(mp_global), intent(in) :: old_sv
logical :: sv_test

! check if any point have land-cover change
if ( mp_global>0 ) then
  sv_test = any( abs(sv-old_sv)>1.e-8 )
else
  sv_test = .false. 
end if

if ( mp_global>0 ) then
  if ( sv_test ) then

    ! check for errors prior to redistribution
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid input temperature for CABLE redistribute_tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      call ccmpi_abort(-1)
    end if     
    
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
    
    ! Do we need a special treatment for snow?
    do k = 1,3
      call redistribute_work(old_sv,ssnow%tggsn(:,k))
      call redistribute_work(old_sv,ssnow%smass(:,k))
      call redistribute_work(old_sv,ssnow%ssdn(:,k))
      call redistribute_work(old_sv,ssnow%sdepth(:,k))
      call redistribute_work(old_sv,ssnow%sconds(:,k))
    end do
    call redistribute_work(old_sv,ssnow%ssdnn(:))
    call redistribute_work(old_sv,ssnow%snowd(:))
    call redistribute_work(old_sv,ssnow%osnowd(:))
    call redistribute_work(old_sv,ssnow%snage(:))
    call redistribute_work(old_sv,ssnow%rtsoil(:))
    call redistribute_work(old_sv,ssnow%GWwb(:))
    call redistribute_work(old_sv,ssnow%wtd(:))
    call redistribute_work(old_sv,canopy%cansto(:))
    call redistribute_work(old_sv,ssnow%pudsto(:))
    call redistribute_work(old_sv,ssnow%wetfac(:))
    
    ! also treatment for carbon
    if ( ccycle/=0 ) then
      do k = 1,3    
        call redistribute_work(old_sv,casapool%cplant(:,k))  
        call redistribute_work(old_sv,casapool%clitter(:,k))
        call redistribute_work(old_sv,casapool%csoil(:,k))
        call redistribute_work(old_sv,casapool%nplant(:,k))
        call redistribute_work(old_sv,casapool%nlitter(:,k))
        call redistribute_work(old_sv,casapool%nsoil(:,k))
        call redistribute_work(old_sv,casapool%pplant(:,k))
        call redistribute_work(old_sv,casapool%plitter(:,k))
        call redistribute_work(old_sv,casapool%psoil(:,k))
      end do    
      !call redistribute_work(old_sv,casamet%glai)
      !call redistribute_work(old_sv,phen%phen)
      !call redistribute_work(old_sv,phen%aphen)
    end if

    ! check for errors after redistribution
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid output temperature for CABLE redistribute_tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      call ccmpi_abort(-1)
    end if  
    
  end if
end if
  
return
end subroutine redistribute_tile

subroutine redistribute_work(old_sv,vdata)

use newmpar_m
use soil_m

integer tile, nb, iq, is, ie
real, dimension(mp_global), intent(in) :: old_sv
real, dimension(ifull,maxnb) :: up_new_svs, up_old_svs
real, dimension(ifull) :: svs_sum
real, dimension(maxnb) :: adj_pos_frac, adj_neg_frac, new_svs, old_svs
real(kind=8) :: ave_neg_vdata, adj_neg_frac_sum
real(kind=8), dimension(mp_global), intent(inout) :: vdata
real(kind=8), dimension(ifull,maxnb) :: up_vdata
real(kind=8), dimension(maxnb) :: old_vdata

! update data
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
  where ( land(1:ifull) )  
    up_new_svs(:,nb) = up_new_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do
svs_sum = sum(up_old_svs,dim=2)
do nb = 1,maxnb
  where ( land(1:ifull) )  
    up_old_svs(:,nb) = up_old_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  do iq = is,ie
    new_svs(:) = up_new_svs(iq,:)
    old_svs(:) = up_old_svs(iq,:)      
    adj_pos_frac(:) = max( new_svs(:)-old_svs(:), 0. ) ! tiles with increasing fraction
    adj_neg_frac(:) = max( old_svs(:)-new_svs(:), 0. ) ! tiles with decreasing fraction
    adj_neg_frac_sum = sum( adj_neg_frac )   
    ! check adj_neg_frac_sum>0 which can occur if water has changed to land
    if ( land(iq) .and. any( abs(new_svs-old_svs)>1.e-8 ) .and. adj_neg_frac_sum>1.e-8 ) then 
      old_vdata(:) = up_vdata(iq,:)
      ! summarise tiles with decreasing area fraction with an average value
      ave_neg_vdata = sum( adj_neg_frac(:)*old_vdata(:) )/adj_neg_frac_sum
      ! Only change tiles that are increasing in area fraction
      do nb = 1,maxnb
        if ( adj_pos_frac(nb)>1.e-8 ) then  
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
subroutine savetiledef(idnc,local,jdim,jsize,cdim,csize,itype)

use carbpools_m
use cc_mpi, only : myid
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
  
integer, intent(in) :: idnc, jsize
integer k,n
integer ll,dd,hh
integer, dimension(2), intent(in) :: csize
integer, dimension(jsize), intent(in) :: jdim  
integer, dimension(6,7), intent(in) :: cdim
character(len=80) vname
character(len=80) lname
logical, intent(in) :: local
integer, intent(in) :: itype
  
if (myid==0.or.local) then
  if (myid==0) write(6,*) "-> define CABLE tile data"
  if ( itype==-1 ) then !just for restart file
    do n = 1,maxtile  
      write(lname,'("Veg type tile ",I1.1)') n
      write(vname,'("t",I1.1,"_cvc")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.,any_m,point_m,double_m)
    end do      
    do n = 1,maxtile
      write(lname,'("Veg fraction tile ",I1.1)') n
      write(vname,'("t",I1.1,"_svs")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.,any_m,point_m,double_m)
      do k = 1,ms
        write(lname,'("Soil temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,any_m,point_m,double_m)
      end do
      do k = 1,ms
        write(lname,'("Soil moisture tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wb",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,any_m,point_m,double_m)
      end do
      do k = 1,ms
        write(lname,'("Soil ice tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wbice",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,any_m,point_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,any_m,point_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow mass tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_smass",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'K',0.,650.,any_m,point_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow density tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,any_m,point_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow depth tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow sconds tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sconds",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,any_m,point_m,double_m)
      end do
      write(lname,'("Snow ssdnn tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ssdnn")') n
      call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,any_m,point_m,double_m)
      write(lname,'("Snow flag tile ",I1.1)') n
      write(vname,'("t",I1.1,"_sflag")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6.5,any_m,point_m,double_m)
      write(lname,'("Snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,double_m)
      write(lname,'("Old snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_osnd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,double_m)
      write(lname,'("Snow age tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snage")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,26.,any_m,point_m,double_m)
      write(lname,'("Soil turbulent resistance tile ",I1.1)') n
      write(vname,'("t",I1.1,"_rtsoil")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3e5,any_m,point_m,double_m)
      write(lname,'("Aquifer moisture content ",I1.1)') n
      write(vname,'("t",I1.1,"_GWwb")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,1.3e5,any_m,point_m,double_m)
      write(lname,'("Water table depth ",I1.1)') n
      write(vname,'("t",I1.1,"_wtd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m',0.,6.5e4,any_m,point_m,double_m)
      write(lname,'("cansto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_cansto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,any_m,point_m,double_m)
      write(lname,'("us tile ",I1.1)') n
      write(vname,'("t",I1.1,"_us")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m/s',0.,13.,any_m,point_m,double_m) 
      write(lname,'("pudsto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_pudsto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,any_m,point_m,double_m)
      write(lname,'("wetfac tile ",I1.1)') n
      write(vname,'("t",I1.1,"_wetfac")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,any_m,point_m,double_m)
      write(lname,'("ga tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ga")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',-6500.,6500.,any_m,point_m,double_m)
    end do
    if ( soil_struc==1 ) then
      do n = 1,maxtile  
        write(lname,'("hzero tile ",I1.1)') n
        write(vname,'("t",I1.1,"_hzero")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        do k = 1,ms
          write(lname,'("S tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_s",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,ms
          write(lname,'("tsoil tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,ms
          write(lname,'("thetai tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do  
        write(lname,'("snowliq tile ",I1.1," lev ",I1.1)') n,1
        write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("tsurface tile ",I1.1)') n
        write(vname,'("t",I1.1,"_tsurface")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("nsnow tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsnow")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("fwsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fwsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
      end do
    end if
    if ( ccycle/=0 ) then
      do n = 1,maxtile  
        do k = 1,mplant
          write(lname,'("C leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,mplant
          write(lname,'("N leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,mplant
          write(lname,'("P leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("C litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("N litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("P litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("C soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("N soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("P soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        write(lname,'("Prognostic LAI tile ",I1.1)') n
        write(vname,'("t",I1.1,"_glai")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Leaf phenology phen tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Leaf phenology rainfall ",I1.1)') n
        write(vname,'("t",I1.1,"_aphen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Leaf phenology phase tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phenphase")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Leaf phenology doyphase3 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_doyphase3")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("C labile tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabile")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("N soilmin tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsoilmin")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("P soillab tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoillab")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("P soilsorb tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilsorb")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("P soilocc tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilocc")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        do k = 1,mplant
          write(lname,'("crmplant tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        end do
        write(lname,'("frac_sapwood tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fracsapwood")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("sapwoodarea tile ",I1.1)') n
        write(vname,'("t",I1.1,"_sapwoodarea")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("crsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("cnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_cnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("clabloss tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabloss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("crgplant tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crgplant")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("stemnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_stemnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("LAImax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_LAImax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Cleafmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Cleafmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("Crootmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Crootmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("fpn ",I1.1)') n
        write(vname,'("t",I1.1,"_fpn")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
        write(lname,'("frday ",I1.1)') n
        write(vname,'("t",I1.1,"_frday")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,double_m)
      end do
    end if
  end if ! itype==-1
  if ( cable_pop==1 ) then
    do n = 1,maxtile  
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( diaglevel_pop > 0 ) then
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
      end if
      if ( itype==-1 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_densindiv")') n
        write(vname,'("t",I1.1,"_pop_grid_densindiv")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_height_mean")') n
        write(vname,'("t",I1.1,"_pop_grid_height_mean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_height_max")') n
        write(vname,'("t",I1.1,"_pop_grid_height_max")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_basal_area")') n
        write(vname,'("t",I1.1,"_pop_grid_basal_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_res_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_area_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_cover")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_area")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_volume")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_KClump")') n
        write(vname,'("t",I1.1,"_pop_grid_KClump")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq_age")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
        call attrib(idnc,cdim(:,7),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_biomass_age")') n
        write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
        call attrib(idnc,cdim(:,7),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do      
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do         
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do         
        do dd = 1,POP_NDISTURB
          write(lname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do         
      end if
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_patch_id")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq")') n
        write(vname,'("t",I1.1,"_pop_grid_freq")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq_old")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_lai")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,double_m)
        end do
      end if
    end do  
  end if   
  if ( itype==-1 ) then !just for restart file
    !do n=1,maxtile
      !write(lname,'("Sensible correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fhscor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,any_m,point_m,double_m)
      !write(lname,'("Latent correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fescor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,any_m,point_m,double_m)
    !end do
    lname='DIR VIS albedo'
    vname='albvisdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
    lname='DIF VIS albedo'
    vname='albvisdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
    lname='DIR NIR albedo'
    vname='albnirdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
    lname='DIF NIR albedo'
    vname='albnirdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
    lname='VIS albedo'
    vname='albvis'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
    lname='NIR albedo'
    vname='albnir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,double_m)
  end if
end if
  
return
end subroutine savetiledef

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,iarch,itype)

use carbpools_m
use cc_mpi, only : ccmpi_abort
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
use soil_m
use soilsnow_m
use vegpar_m
  
integer, intent(in) :: idnc,iarch,itype
integer k,n,np_pop,iq
integer cc,dd,hh,ll
integer, dimension(ifull) :: dati
real, dimension(ifull) :: datr
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(mp_global) :: dummy_unpack
real(kind=8), dimension(:), allocatable :: dat_in
real(kind=8), dimension(:,:), allocatable :: datpatch
real(kind=8), dimension(:,:), allocatable :: datage
real(kind=8), dimension(:,:,:), allocatable :: datpc
character(len=80) vname
logical, intent(in) :: local
  
if ( itype==-1 ) then !just for restart file
  do n = 1,maxtile
    dat = 0
    if ( n<=maxnb ) then
      dummy_unpack = real(cveg,8)    
      call cable_unpack(dummy_unpack,dat,n)
    end if
    write(vname,'("t",I1.1,"_cvc")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do      
  
  ! soil temperature check
  if ( mp_global>0 ) then
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid CABLE temperature when writing tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      stop -1
    end if
  end if
  
  do n = 1,maxtile  ! tile
    datr = 0.
    if ( n<=maxnb ) call cable_unpack(sv,datr,n)
    dat = real( datr, 8 )
    write(vname,'("t",I1.1,"_svs")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    do k = 1,ms     ! soil layer
      dat = real(tgg(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%tgg(:,k),dat,n)
      if ( any( dat(1:ifull)>400._8 .and. land(1:ifull) ) ) then
        write(6,*) "ERROR: Invalid CABLE soil temperature when writing tile"
        write(6,*) "n,k = ",n,k
        write(6,*) "tgg ",maxval(dat(1:ifull),mask=land(1:ifull))
        stop -1
      end if
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
    dat = 0._8
    if ( n<=maxnb ) call cable_unpack(ssnow%wtd,dat,n)
    write(vname,'("t",I1.1,"_wtd")') n
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
end if ! itype==-1
if ( cable_pop==1 ) then
  allocate( datpatch(ifull,POP_NPATCH) )  
  allocate( datage(ifull,POP_AGEMAX) )  
  allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
  np_pop = size( pop%pop_grid(:) )
  allocate( dat_in(np_pop) )
  do n = 1,maxtile
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    dat = 0._8
    dat_in = 0._8
    if ( diaglevel_pop > 0 ) then
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cmass_sum 
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end if
    if ( itype==-1 ) then !just for restart file
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cmass_sum_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cheartwood_sum 
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%csapwood_sum
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%csapwood_sum_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%densindiv  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%height_mean  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%height_max  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%basal_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_loss  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area_loss  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%stress_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crowding_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%fire_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cat_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%res_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%growth  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%area_growth  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_cover  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_volume  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%KClump  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then      
        do k = 1,POP_AGEMAX
          dat_in(1:np_pop) = pop%pop_grid(:)%freq_age(k)  
          call pop_unpack(dat_in(1:np_pop),datage(:,k),n)
        end do
      end if          
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then      
        do k = 1,POP_AGEMAX
          dat_in(1:np_pop) = pop%pop_grid(:)%biomass_age(k)  
          call pop_unpack(dat_in(1:np_pop),datage(:,k),n)
        end do
      end if          
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%biomass(ll)
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%density(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%hmean(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%hmax(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%cmass_stem_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%densindiv_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%height_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%diameter_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = real(pop%pop_grid(:)%n_age(dd),8)
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
    end if
    if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
      if ( n<=maxnb ) then  
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%id,8)
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%freq(k)  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%freq_old(k)   
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%factor_recruit  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%pgap  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%lai  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%biomass  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%biomass_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)    
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%heartwood  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%stress_mortality   
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%fire_mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%cat_mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%crowding_mortality 
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%cpc  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_loss  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area_loss  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%growth  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%area_growth  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_NPP  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_respiration  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_light_uptake  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do dd = 1,POP_NDISTURB  
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%disturbance_interval(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%first_disturbance_year(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%age(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%ranked_age_unique(k,dd),8)  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%freq_ranked_age_unique(k,dd)  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%ncohort,8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%biomass  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%density  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%hmean 
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%hmax  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age,8)
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT 
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id,8)
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT  
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration 
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT  
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH 
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do 
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then 
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH 
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
    end if
  end do
  deallocate( datpatch )
  deallocate( datage )
  deallocate( datpc )
  deallocate( dat_in )
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

subroutine cable_casatile(dat_out,vname,k,n)

use newmpar_m

implicit none

integer, intent(in) :: k, n
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
    call cable_casatile_work(n,casapool%cplant(:,k),dat_out)
  case("nplant")
    select case(k)
      case(1)
        dat_out(:) = real(nleaf(n))
      case(2)
        dat_out(:) = real(nwood(n))
      case(3)
        dat_out(:) = real(nfroot(n))
    end select
    call cable_casatile_work(n,casapool%nplant(:,k),dat_out)
  case("pplant")
    select case(k)
      case(1)
        dat_out(:) = real(xpleaf(n))
      case(2)
        dat_out(:) = real(xpwood(n))
      case(3)
        dat_out(:) = real(xpfroot(n))  
    end select      
    call cable_casatile_work(n,casapool%pplant(:,k),dat_out)
  case("clitter")
    select case(k)
      case(1)
        dat_out(:) = real(cmet(n))
      case(2)
        dat_out(:) = real(cstr(n))
      case(3)
        dat_out(:) = real(ccwd(n))
    end select
    call cable_casatile_work(n,casapool%clitter(:,k),dat_out)
  case("nlitter")
    select case(k)
      case(1)
        dat_out(:) = real(nmet(n))
      case(2)
        dat_out(:) = real(nstr(n))
      case(3)
        dat_out(:) = real(ncwd(n))
    end select      
    call cable_casatile_work(n,casapool%nlitter(:,k),dat_out)
  case("plitter")
    select case(k)
      case(1)
        dat_out(:) = real(xpmet(n))
      case(2)
        dat_out(:) = real(xpstr(n)) 
      case(3)
        dat_out(:) = real(xpcwd(n))
    end select      
    call cable_casatile_work(n,casapool%plitter(:,k),dat_out)
  case("csoil")
    select case(k)
      case(1)
        dat_out(:) = real(cmic(n))
      case(2)
        dat_out(:) = real(cslow(n))
      case(3)
        dat_out(:) = real(cpass(n))  
    end select
    call cable_casatile_work(n,casapool%csoil(:,k),dat_out)
  case("nsoil")
    select case(k)
      case(1)
        dat_out(:) = real(nmic(n))
      case(2)
        dat_out(:) = real(nslow(n))
      case(3)
        dat_out(:) = real(npass(n))
    end select      
    call cable_casatile_work(n,casapool%nsoil(:,k),dat_out)
  case("psoil")
    select case(k)
      case(1)
        dat_out(:) = real(xpmic(n))
      case(2)
        dat_out(:) = real(xpslow(n))
      case(3)
        dat_out(:) = real(xppass(n))
    end select      
    call cable_casatile_work(n,casapool%psoil(:,k),dat_out)
  case default
    write(6,*) "ERROR: Unknown option for cable_tile ",trim(vname)
    stop -1
end select

end subroutine cable_casatile

subroutine cable_casatile_work(n,dat_in,dat_out)

use newmpar_m

implicit none

integer, intent(in) :: n
integer nb
real, dimension(:), intent(inout) :: dat_out
real(kind=8), dimension(:), intent(in) :: dat_in
integer, dimension(size(dat_out)) :: iveg_local
real, dimension(size(dat_out)) :: dat_local

do nb = 1,maxnb
  iveg_local(:) = 0  
  call cable_unpack(veg%iveg,iveg_local,nb)
  dat_local(:) = 0.
  call cable_unpack(dat_in,dat_local,nb)
  where ( iveg_local(:)==n .and. dat_local(:)>1.e-8 )
    dat_out(:) = dat_local(:)
  end where
end do

return
end subroutine cable_casatile_work

! *************************************************************************************
! Water inflow from river routing
subroutine cableinflow(inflow)

use newmpar_m
use soil_m

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
use parm_m

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
  write(6,*) "-> Reading CASA leaf phenology data ",trim(fphen)
  call ccnf_open(fphen,ncid,ierr)
  if ( ierr==0 ) then
    ncfile = .true.
    if ( nmaxpr==1 ) write(6,*) "-> Found netcdf file ",trim(fphen)
  else
    call ccnf_open(fphen//'.nc',ncid,ierr)
    if ( ierr==0 ) then
      ncfile = .true.
      if ( nmaxpr==1 ) write(6,*) "-> Found netcdf file ",trim(fphen)
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

! Subroutines for ground water transport - average water table height

subroutine calc_wt_ave(wth_ave,gwwb_min,gwdz)

use arrays_m
use const_phys
use map_m
use newmpar_m
use parm_m
use soil_m

integer nb, js, je, is, ie, ks, ke, tile, ls, le
integer lmp, lmaxnb, k
integer, dimension(maxtile,2) :: ltind
real, dimension(:), intent(out) :: wth_ave, gwwb_min, gwdz
real, dimension(size(wth_ave)) :: wtd_ave
logical, dimension(imax,maxtile) :: ltmap
type(soil_parameter_type) :: lsoil
type(soil_snow_type) :: lssnow

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  lmp = tdata(tile)%mp

  wth_ave(is:ie) = 0.
  wtd_ave(is:ie) = 0.
  gwwb_min(is:ie) = 9.e9
  gwdz(is:ie) = 0.

  if ( lmp>0 ) then

    lmaxnb = tdata(tile)%maxnb
    ltmap = tdata(tile)%tmap
    js = tdata(tile)%toffset + 1
    je = tdata(tile)%toffset + tdata(tile)%mp
    ltind = tdata(tile)%tind - tdata(tile)%toffset
    call setp(soil,lsoil,tile)
    call setp(ssnow,lssnow,tile)

    ! soil layer thickness
    do k = 1,ms
      where ( land(is:ie) )
        gwdz(is:ie) = gwdz(is:ie) + lsoil%zse(k)
      end where
    end do

    ! calculate water table depth and other information from tiles
    do nb = 1,maxnb
      ks = ltind(nb,1)
      ke = ltind(nb,2)
      ls = js + ks - 1
      le = js + ke - 1
      wtd_ave(is:ie) = wtd_ave(is:ie) + unpack( sv(ls:le)*real(lssnow%wtd(ks:ke))/1000., ltmap(:,nb), 0. )
      gwwb_min(is:ie) = min( gwwb_min(is:ie), unpack( real(lssnow%GWwb(ks:ke)), ltmap(:,nb), 9.e9 ) )
      gwdz(is:ie) = gwdz(is:ie) + unpack( sv(ls:le)*real(lsoil%GWdz(ks:ke)), ltmap(:,nb), 0. )
    end do

  end if ! mp>0

  ! adjust water table depth for orography
  wth_ave(is:ie) = zs(is:ie)/grav - wtd_ave(is:ie)

  ! set minimum ground water to zero where undefined (e.g., ocean)
  where ( gwwb_min(is:ie)>=1.e9 )
    gwwb_min(is:ie) = 0.
  elsewhere
    gwwb_min(is:ie) = (ds/em(is:ie))**2*gwwb_min(is:ie)
  end where

  ! redistribute ground water across tiles within the gridbox?

end do ! tile

return
end subroutine calc_wt_ave

subroutine calc_wt_flux(flux,flux_m,flux_c,ddt)

use arrays_m
use const_phys
use map_m
use newmpar_m
use parm_m

integer nb, is, ie, ks, ke, tile, ls, le, js, je
integer lmp, lmaxnb
integer, dimension(maxtile,2) :: ltind
real, intent(in) :: ddt
real, dimension(ifull), intent(in) :: flux, flux_m, flux_c
real(kind=r_2), dimension(ifull) :: tot_flux, flux_corr
real(kind=r_2), dimension(ifull) :: wth, mm, cc, flux_adj, indxsq ! ifull>=mp
real(kind=r_2), dimension(ifull,maxtile) :: new_flux              ! ifull>=mp
logical, dimension(imax,maxtile) :: ltmap
type(soil_snow_type) :: lssnow

tot_flux(:) = 0._r_2

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  lmp = tdata(tile)%mp

  if ( lmp>0 ) then

    lmaxnb = tdata(tile)%maxnb
    ltmap = tdata(tile)%tmap
    js = tdata(tile)%toffset + 1
    je = tdata(tile)%toffset + tdata(tile)%mp
    ltind = tdata(tile)%tind(:,:) - tdata(tile)%toffset
    call setp(ssnow,lssnow,tile)

    ! remove remaing flux from ground water

    do nb = 1,lmaxnb
      ks = ltind(nb,1)
      ke = ltind(nb,2)
      ls = js + ks - 1
      le = js + ke - 1
      wth(ks:ke) = pack( real(zs(is:ie)/grav,r_2), ltmap(:,nb) ) - lssnow%wtd(ks:ke)/1000._r_2
      mm(ks:ke) = pack( real(flux_m(is:ie),r_2), ltmap(:,nb) )
      cc(ks:ke) = pack( real(flux_c(is:ie),r_2), ltmap(:,nb) )
      new_flux(ks:ke,nb) = mm(ks:ke)*wth(ks:ke) + cc(ks:ke)
      tot_flux(is:ie) = tot_flux(is:ie) + &
        unpack( real(sv(ls:le),r_2)*new_flux(ks:ke,nb), ltmap(:,nb), 0._r_2 )
    end do

  end if
end do

where( abs(tot_flux)>0._r_2 )
  flux_corr = real(abs(flux),r_2)/abs(tot_flux)
elsewhere
  flux_corr = 0._r_2
end where

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  lmp = tdata(tile)%mp

  if ( lmp>0 ) then

    lmaxnb = tdata(tile)%maxnb
    ltmap = tdata(tile)%tmap
    ltind = tdata(tile)%tind - tdata(tile)%toffset
    call setp(ssnow,lssnow,tile)

    do nb = 1,lmaxnb
      ks = ltind(nb,1)
      ke = ltind(nb,2)
      flux_adj(ks:ke) = pack( real(flux_corr(is:ie),r_2 ), ltmap(:,nb) )
      indxsq(ks:ke) = pack( real((em(is:ie)/ds)**2,r_2), ltmap(:,nb) )
      new_flux(ks:ke,nb) = new_flux(ks:ke,nb)*flux_adj(ks:ke)
      lssnow%GWwb(ks:ke) = lssnow%GWwb(ks:ke) - real(ddt,r_2)*indxsq(ks:ke)*new_flux(ks:ke,nb)
    end do

  end if ! lmp>0.

end do ! tile

return
end subroutine calc_wt_flux

end module cable_ccam

