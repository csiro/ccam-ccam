c      common block tfcom contains transmission functions used for
c      radiative computations, and output heating rates and fluxes, 
c      except those needed out of the radiative module:
 
      real to3(imax,lp1,lp1)    ! Transmission fctn for the 990-1070 cm-1 band
                                ! O3(9.6 um) + H2O continuum (no lines) 
      real co21(imax,lp1,lp1)   ! Transmission fctn for the 560-800 cm-1 band 
                                ! (as 1 band). Includes CO2 (in lwr88) and
                                ! h2o(l+c) after multiplication with "over"
                                ! in fst88       
      real emiss(imax,lp1,lp1)  ! e2 emissivityy fctn for H2O lines
                                ! (0-560,1200-2200 cm-1). Obtained in e1e288.
      real emiss2(imax,lp1,lp1) ! Transmission fctn for H2O continuum in the
                                ! 800-990 and 1070-1200 cm-1 region, taken 
                                ! as 1 band 
      real avephi(imax,lp1,lp1) ! H2O optical paths bet. flux pressures: 
                                ! input to emissivity calculations. 
      real cts(imax,l)          ! Approx CTS heating rates for 160-560 and
                                ! 800-990, 1070-1200 cm-1 ranges
      real ctso3(imax,l)        ! Approx CTS heating rates for 560-800, 
                                ! 990-1070 cm-1 ranges 
      real excts(imax,l)        ! Exact CTS heating rates for 160-1200 cm-1
                                ! range 
      real exctsn(imax,l,nbly)  ! Exact CTS heating rates, by bands 
      real e1flx(imax,lp1)      ! e1 emissivity fctn for H2O lines
                                ! (0-560,1200-cm-1)
      real co2nbl(imax,l)       ! CO2 trans. fctns. (not pressure-integrated)
                                ! for adjacent levels,over the 560-800 cm-1 
                                ! range. 
      real co2sp1(imax,lp1)     ! CO2 trans. fctns. (not pressure-integrated)
                                ! bet. a flux level and space, for the 560-670
                                ! cm-1 range. Used for exact CTS calcs. 
      real co2sp2(imax,lp1)     ! Same as co2sp1, but for the 670-800 cm-1
                                ! range. 
      real co2sp(imax,lp1)      ! Same as co2sp1, but for the 560-800 cm-1
                                ! band. Used for approx CTS calcs. 
      real to3spc(imax,l)       ! O3 optical depths bet. a level and space. 
                                ! Used for exact CTS calcs. 
      real totvo2(imax,lp1)     ! H2O continuum optical paths bet. space and a
                                ! level, using the cnt. coefficient for the
                                ! 1-band 800-990,1070-1200 cm-1 band. Used for 
                                ! CTS calcs. 
 
      common /tfcom/ to3, co21, emiss, emiss2, avephi, cts, ctso3,
     &               excts, exctsn, e1flx, co2nbl, co2sp1, co2sp2,
     &               co2sp, to3spc, totvo2
