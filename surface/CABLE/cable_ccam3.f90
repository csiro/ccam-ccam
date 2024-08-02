! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
 
use cable_data_module
use cable_def_types_mod
use casavariable
use phenvariable
use pop_types

implicit none

private
public air, bgc, met, bal, rad, rough, ssnow
public sum_flux, climate, veg, soil, canopy
public casabal, casabiome, casaflux, casamet
public casapool, phen, pop, c
public climate_save, climate_save_type

type climate_save_type
  real(kind=r_2), dimension(:), pointer :: APAR_leaf_sun
  real(kind=r_2), dimension(:), pointer :: APAR_leaf_shade
  real(kind=r_2), dimension(:), pointer :: Dleaf_sun
  real(kind=r_2), dimension(:), pointer :: fwsoil
  real(kind=r_2), dimension(:), pointer :: Dleaf_shade
  real(kind=r_2), dimension(:), pointer :: Tleaf_sun 
  real(kind=r_2), dimension(:), pointer :: Tleaf_shade
  real(kind=r_2), dimension(:), pointer :: cs_sun
  real(kind=r_2), dimension(:), pointer :: cs_shade
  real(kind=r_2), dimension(:), pointer :: scalex_sun
  real(kind=r_2), dimension(:), pointer :: scalex_shade
end type climate_save_type
    
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
type (climate_save_type), save   :: climate_save

end module cable_ccam3