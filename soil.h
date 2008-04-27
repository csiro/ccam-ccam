      real zoland,zmin,zolnd,zolog,albsav,so4t ! MJT CHANGE albedo delete alb, add albnirsav
	  real albnirsav                           ! alb is now in soilsnow.h for cable
      logical land
      common/soil1/zoland,zmin
      common/soil2/zolnd(ifull),zolog(ifull),land(ifull)
      !common/soil4/alb(ifull),albsav(ifull),so4t(ifull)
	  common/soil4/albsav(ifull),so4t(ifull),albnirsav(ifull) 

