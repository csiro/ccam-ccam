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
    
module cable_ccam3
 
!use cable_data_module
use cable_def_types_mod
use casavariable
use phenvariable
use pop_types

implicit none

private
public air, bgc, met, bal, rad, rough, ssnow
public sum_flux, climate, veg, soil, canopy
public casabal, casabiome, casaflux, casamet
public casapool, phen, pop
   
type (air_type), dimension(:), allocatable, save            :: air
type (bgc_pool_type), dimension(:), allocatable, save       :: bgc
type (met_type), dimension(:), allocatable, save            :: met
type (balances_type), dimension(:), allocatable, save       :: bal
type (radiation_type), dimension(:), allocatable, save      :: rad
type (roughness_type), dimension(:), allocatable, save      :: rough
type (soil_snow_type), dimension(:), allocatable, save      :: ssnow
type (sum_flux_type), dimension(:), allocatable, save       :: sum_flux
type (climate_type), dimension(:), allocatable, save        :: climate
type (veg_parameter_type), dimension(:), allocatable, save  :: veg
type (soil_parameter_type), dimension(:), allocatable, save :: soil
type (canopy_type), dimension(:), allocatable, save         :: canopy
type (casa_balance), dimension(:), allocatable, save        :: casabal
type (casa_biome), dimension(:), allocatable, save          :: casabiome
type (casa_flux), dimension(:), allocatable, save           :: casaflux
type (casa_met), dimension(:), allocatable, save            :: casamet
type (casa_pool), dimension(:), allocatable, save           :: casapool
type (phen_variable), dimension(:), allocatable, save       :: phen
type (pop_type), dimension(:), allocatable, save            :: pop

end module cable_ccam3