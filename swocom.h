c     common block swocom contains quantities outputted by the 
c     shortwave radiation code to the external module:
      real fsw(imax,lp1)  ! Net radiation (up-down) in CGS units at all
                          ! pressure levels
      real dfsw(imax,lp1) ! Downward radiation at all pressure levels
      real ufsw(imax,lp1) ! Upward radiation at all pressure levels
      real hsw(imax,l)    ! Shortwave heating rates in K/day for pressure
                          ! layers. 
      common / swocom / fsw, dfsw, ufsw, hsw
