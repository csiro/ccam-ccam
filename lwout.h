c     common block lwout contains the quantities outputted by the 
c     longwave radiation code to the external module: 

      real heatra(imax,l) ! Heating rate at data levels (K/day) 
      real grnflx(imax)   ! Net longwave flux at the ground (CGS units) 
      real topflx(imax)   ! Net longwave flux at the top    (CGS units) 
      real grnflxclr(imax)
      real exctsclr(imax,l)
      real ctso3clr(imax,l)
c 
      common / lwout / heatra, grnflx, topflx,
     &                 grnflxclr, exctsclr, ctso3clr
