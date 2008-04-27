! eak 16/03/06
      real condxpr,fev, fes, ga, dgdtg
      common /soilsn/ condxpr(ifull) ! precip. reaching ground 
      common /soilsn/ fev(ifull)     ! canopy transpiration 
      common /soilsn/ fes(ifull)     ! soil evaporation 
!      common /soilsn/ fwtop(ifull)   ! net water flux at the surface 
      common /soilsn/ ga(ifull)      ! net heat flux at the surface 
      common /soilsn/ dgdtg(ifull)   ! d(ga)/ d(tgg) 

