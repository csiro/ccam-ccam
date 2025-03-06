! UCLEM urban canopy model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the UCLEM urban canopy model
!
! UCLEM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UCLEM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UCLEM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! definition of types used by UCLEM urban climate and energy model
    
module uclem_types

implicit none

private
public facetdata, facetparams, hydrodata, vegdata, intldata, fparmdata, pdiagdata

type facetdata
  real(kind=8), dimension(:,:), allocatable :: nodetemp ! Temperature of node (prognostic)       [K]
  real(kind=8), dimension(:,:), allocatable :: storage  ! Facet energy storage (diagnostic)
end type facetdata

type facetparams
  real, dimension(:,:), allocatable :: depth         ! Layer depth                              [m]
  real, dimension(:,:), allocatable :: volcp         ! Layer volumetric heat capacity           [J m^-3 K-1]
  real, dimension(:,:), allocatable :: lambda        ! Layer conductivity                       [W m^-1 K^-1]
  real, dimension(:),   allocatable :: alpha         ! Facet albedo (internal & external)
  real, dimension(:),   allocatable :: emiss         ! Facet emissivity (internal & external)
end type facetparams

type hydrodata
  real, dimension(:), allocatable   :: surfwater
  real, dimension(:), allocatable   :: leafwater
  real, dimension(:), allocatable   :: soilwater
  real, dimension(:), allocatable   :: snow
  real, dimension(:), allocatable   :: den
  real, dimension(:), allocatable   :: snowalpha
end type hydrodata

type vegdata
  real, dimension(:), allocatable :: temp          ! Temperature of veg (prognostic)  [K]
  real, dimension(:), allocatable :: sigma         ! Fraction of veg on roof/canyon
  real, dimension(:), allocatable :: alpha         ! Albedo of veg
  real, dimension(:), allocatable :: emiss         ! Emissivity of veg
  real, dimension(:), allocatable :: lai           ! Leaf area index of veg
  real, dimension(:), allocatable :: zo            ! Roughness of veg
  real, dimension(:), allocatable :: rsmin         ! Minimum stomatal resistance of veg
end type vegdata

type intldata
  real, dimension(:,:,:), allocatable :: psi   ! internal radiation
  real, dimension(:,:,:), allocatable :: viewf ! internal radiation
end type intldata

type fparmdata
  real, dimension(:), allocatable :: hwratio,effhwratio,sigmabld
  real, dimension(:), allocatable :: industryfg,intgains_flr,trafficfg,bldheight,bldwidth
  real, dimension(:), allocatable :: ctime,vangle,hangle,fbeam,weekdayload
  real, dimension(:), allocatable :: bldairtemp
  real, dimension(:), allocatable :: swilt,sfc,ssat,rfvegdepth
  real, dimension(:), allocatable :: infilach,ventilach,heatprop,coolprop
  real, dimension(:), allocatable :: sigmau
  integer, dimension(:), allocatable :: intmassn
end type fparmdata

type pdiagdata
  real, dimension(:), allocatable :: lzom, lzoh, cndzmin, cduv, cdtq
  real, dimension(:), allocatable :: tscrn, qscrn, uscrn, u10, emiss, snowmelt
  real, dimension(:), allocatable :: bldheat, bldcool, traf, intgains_full
  real, dimension(:), allocatable :: irrig,acond_vegw
  real, dimension(:), allocatable :: surfrunoff,soilwetness,soilwater
  real, dimension(:), allocatable :: transveg,soilmoist,delsoilmoist
  real, dimension(:), allocatable :: rootmoistc
  real, dimension(:), allocatable :: frac_sigma
  real, dimension(:), allocatable :: delintercept, snowt, vegt, swe, surfstor
  real, dimension(:), allocatable :: snowfrac, salbedo, calbedo, taircanyon, ulai
  real, dimension(:), allocatable :: delswe
  real, dimension(:), allocatable :: roof_water_runoff, roof_snow_runoff, roof_soil_runoff
  real, dimension(:), allocatable :: road_water_runoff, road_snow_runoff, road_soil_runoff
  real, dimension(:), allocatable :: ac_heat_on, ac_cool_on
  real(kind=8), dimension(:), allocatable :: surferr, atmoserr, surferr_bias, atmoserr_bias
  real(kind=8), dimension(:), allocatable :: storage_flux
  ! real(kind=8), dimension(:,:), allocatable :: storagetot_road, storagetot_walle, storagetot_wallw, storagetot_roof
  logical :: first_call
end type pdiagdata

end module uclem_types
