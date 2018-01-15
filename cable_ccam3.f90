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

end module cable_ccam3