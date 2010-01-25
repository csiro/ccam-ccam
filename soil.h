      real zoland,zmin,zolnd,zolog,albsav,so4t,albnirsav ! MJT CHANGE albedo delete alb, add albnirsav
      real albvisdif,albnirdif,albvisdir,albnirdir ! MJT radiation
      logical land
      common/soil1/zoland,zmin
      common/soil2/zolnd(ifull),zolog(ifull),land(ifull)
      common/soil4/albsav(ifull),so4t(ifull),albnirsav(ifull) ! MJT CHANGE albedo delete alb, add albnirsav
      common/soil4/albvisdif(ifull),albnirdif(ifull)  ! MJT radiation
      common/soil4/albvisdir(ifull),albnirdir(ifull)  ! MJT radiation
