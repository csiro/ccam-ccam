      real zoland,zmin,zolnd,zolog,alb,albsav,so4t
      logical land
      common/soil1/zoland,zmin
! rml & after column 72 so that it will compile with both .f and .f90 files
      common/soil2/zolnd(ifull),zolog(ifull),land(ifull)                
      common/soil4/alb(ifull),albsav(ifull),so4t(ifull)
