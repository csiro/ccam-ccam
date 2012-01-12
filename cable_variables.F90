      module cable_variables

      use define_types
      implicit none

      private

      TYPE (air_type)             :: air
      TYPE (bgc_pool_type)        :: bgc
      TYPE (met_type)             :: met
      TYPE (balances_type)        :: bal
      TYPE (radiation_type)       :: rad
      TYPE (roughness_type)       :: rough
      TYPE (soil_parameter_type)  :: soil       ! soil parameters
      TYPE (soil_snow_type)       :: ssoil
      TYPE (sum_flux_type)        :: sum_flux
      TYPE (veg_parameter_type)   :: veg        ! vegetation parameters
      TYPE (canopy_type)          :: canopy
         
      ! Save these so only have to do the allocation once.
      save air, bgc, met, bal, rad, rough, soil, ssoil, &
           sum_flux, veg, canopy
      public air, bgc, met, bal, rad, rough, soil, ssoil, &
           sum_flux, veg, canopy

      end module cable_variables
   
