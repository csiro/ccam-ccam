c   The following common blocks contain pretabulated Co2 transmission 
c   functions, evaluated using the methods of Fel and
c   Schwarzkopf (1981) and Schwarzkopf and Fels (1985), 
c***Common co2bd3 contains Co2 transmission functions and temperature 
c   and pressure derivatives for the 560-800 cm-1 band. Also included 
c   are the standard temperatures and the weighting function. These 
c   data are in block data bd3: 

      real co251(lp1,lp1) ! Transmission fctns for t0 (std. profile)
                          ! with p(sfc)=1013.25 mb
      real co258(lp1,lp1) ! Transmission fctns. for t0 (std. profile) 
                          ! with p(sfc)= 810 mb
      real cdt51(lp1,lp1) ! First temperature derivative of co251 
      real cdt58(lp1,lp1) ! First temperature derivative of co258 
      real c2d51(lp1,lp1) ! Second temperature derivative of co251
      real c2d58(lp1,lp1) ! Second temperature derivative of co251
      real co2m51(l)      ! Transmission fctns for t0 for adjacent pressure 
                          ! levels, with no pressure quadrature. Used for
                          ! nearby layer computations. p(sfc)=1013.25 mb 
      real co2m58(l)      ! Same as co2m51,with p(sfc)= 810 mb 
      real cdtm51(l)      ! First temperature derivative of co2m51
      real cdtm58(l)      ! First temperature derivative of co2m58
      real c2dm51(l)      ! Second temperature derivative of co2m51 
      real c2dm58(l)      ! Second temperature derivative of co2m58 
      real stemp(lp1)     ! Standard temperatures for model pressure level
                          ! structure with p(sfc)=1013.25 mb 
      real gtemp(lp1)     ! Weighting function for model pressure level 
                          ! structure with p(sfc)=1013.25 mb.
      real b0             ! Temp. coefficient used for CO2 trans. fctn. 
                          ! correction for t(k). (see ref. 4 and bd3)
      real b1             ! Temp. coefficient, used along with b0 
      real b2             ! Temp. coefficient, used along with b0 
      real b3             ! Temp. coefficient, used along with b0 
 
      common /co2bd3/ co251, co258, cdt51, cdt58, c2d51, c2d58, co2m51,
     &                co2m58, cdtm51, cdtm58, c2dm51, c2dm58, stemp,
     &                gtemp, b0, b1, b2, b3
c 
c***common co2bd2 contains CO2 transmission functions and temperature 
c   and pressure derivatives for the 560-670 cm-1 part of the 15 um 
c   CO2 band.  

      real co231(lp1)     ! Transmission fctns for t0 (std. profile)
                          ! with p(sfc)=1013.25 mb
      real co238(lp1)     ! Transmission fctns. for t0 (std. profile) 
                          ! with p(sfc)= 810 mb
      real cdt31(lp1)     ! First temperature derivative of co231 
      real cdt38(lp1)     ! First temperature derivative of co238 
      real c2d31(lp1)     ! Second temperature derivative of co231
      real c2d38(lp1)     ! Second temperature derivative of co231
 
      common / co2bd2 / co231, co238, cdt31, cdt38, c2d31, c2d38
c 
c***common co2bd4 contains CO2 transmission functions and temperature 
c   and pressure derivatives for the 670-800 cm-1 part of the 15 um 
c   CO2 band.  

      real co271(lp1)     ! Transmission fctns for t0 (std. profile)
                          ! with p(sfc)=1013.25 mb
      real co278(lp1)     ! Transmission fctns. for t0 (std. profile) 
                          ! with p(sfc)= 810 mb
      real cdt71(lp1)     ! First temperature derivative of co271 
      real cdt78(lp1)     ! First temperature derivative of co278 
      real c2d71(lp1)     ! Second temperature derivative of co271
      real c2d78(lp1)     ! Second temperature derivative of co271
 
      common / co2bd4 / co271, co278, cdt71, cdt78, c2d71, c2d78
c 
c***common co2bd5 contains co2 transmission functions for the 2270- 
c   2380 part of the 4.3 um co2 band. 

      real co211(lp1)     ! Transmission fctns for t0 (std. profile)
                          !with p(sfc)=1013.25 mb
      real co218(lp1)     ! Transmission fctns. for t0 (std. profile) 
                          !with p(sfc)= ^810 mb
 
      common / co2bd5 / co211, co218

