c     common block kdacom contains basic quantities computed in 
c     subroutine lwr88 and used in remaining longwave routines: 

      real qh2o(imax,lp1)   ! H2O mass mixing ratio,multiplied by the 
                            ! diffusivity factor (diffctr)
      real p(imax,lp1)      ! Pressure at flux levels of model
      real delp2(imax,l)    ! Pressure difference between flux levels 
      real delp(imax,l)     ! Inverse of delp2
      real t(imax,lp1)      ! Temperature assigned to model flux levels 
      real var1(imax,l)     ! H2O optical path in model layers (between 
                            ! flux levels)
      real var2(imax,l)     ! Pressure-weighted H2O optical path in model layers
      real var3(imax,l)     ! O3 optical path in model layers 
      real var4(imax,l)     ! Pressure-weighted O3 optical path in model layers 
      real cntval(imax,lp1) ! H2O continuum path in model layers for the
                            ! 800-990 and 1070-1200 cm-1 combined band
 
      common /kdacom/ qh2o, p, delp2, delp, t, var1, var2,
     &                var3, var4, cntval
c 
