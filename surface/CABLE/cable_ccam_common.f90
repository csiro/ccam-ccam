! Conformal Cubic Atmospheric Model
    
! Copyright 2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Common variables and subroutines for CABLE-CCAM interface
    
module cable_ccam_common

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
use pop_constants, only : NPATCH, NLAYER, NCOHORT_MAX, HEIGHT_BINS, NDISTURB, AGEMAX
use popmodule, only : pop_init, popstep
use pop_types
use sli_main_mod

implicit none

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
    
end module cable_ccam_common

