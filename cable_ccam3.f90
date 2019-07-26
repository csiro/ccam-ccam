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
  real, dimension(:), pointer :: APAR_leaf_sun
  real, dimension(:), pointer :: APAR_leaf_shade
  real, dimension(:), pointer :: Dleaf_sun
  real, dimension(:), pointer :: fwsoil
  real, dimension(:), pointer :: Dleaf_shade
  real, dimension(:), pointer :: Tleaf_sun 
  real, dimension(:), pointer :: Tleaf_shade
  real, dimension(:), pointer :: cs_sun
  real, dimension(:), pointer :: cs_shade
  real, dimension(:), pointer :: scalex_sun
  real, dimension(:), pointer :: scalex_shade
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