!     this one forces tggsn to directly precede tgg
      common /soilcom/
     .     tggsn(ifull,3) ! snow temperatures
     .   , tgg(ifull,ms)  ! soil temperature (in K)
     .   , wb(ifull,ms)   ! soil moisture (volumetric)
     .   , wbice(ifull,ms)! soil ice
     .   , smass(ifull,3) ! snow masses
     .   , ssdn(ifull,3)  ! snow densities
     .   , ssdnn(ifull)   ! average snow density
     .   , snowd(ifull)   ! snow depth (liquid water)    Fri  09-10-1999
     .   , osnowd(ifull)  ! snow depth from previous time step
     .   , snage(ifull)   ! snow age
     .   , sno(ifull)     ! accum. snow in mm since last write (like precip)
     .   , isflag(ifull)
     .   , gflux(ifull)   !
     .   , sgflux(ifull)  !
     .   , otgsoil(ifull) ! soil or snow surface temper. at the previous dt
!        N.B. otgsoil was missing from common block before 19/8/99
