      subroutine optm(dt,phibar,omega)
      include 'newmpar.h'

! Calculate the optimum acceleration factor for SOR on the cubic-conformal
! grid. This uses a function derived from fitting a set of numerical results.

! Arguments
!     integer il  ! Number of grid points on a panel
      real dt     ! Time step
      real phibar ! Equivalent depth
      real omega
      parameter(radsq=6.371e6**2)

      real lambda
      lambda = 4 * radsq / ( phibar * dt**2 )

      if(npanels.eq.5)then
        if ( il.eq.20 .or. il.eq.21 ) then
           omega = 1 + 0.8715 * (1 + 0.0291874*lambda) /
     &               (1 + 0.0540404*lambda + 153.109e-6*lambda**2)
        else if ( il.eq.30 ) then
           omega = 1 + 0.9083 * (1 + 0.019229*lambda) /
     &               (1 + 0.0340753*lambda + 41.772e-6*lambda**2)
        else if ( il.eq.37 ) then
           omega = 1 + 0.9250 * (1 + 0.0152456*lambda) /
     &                  (1 + 0.0261757*lambda + 19.6856e-6*lambda**2)
        else if ( il.eq.48 ) then
           omega = 1 + 0.94239 * (1 + 0.01129 *lambda) /
     &                  (1 + 0.01853*lambda + 7.6e-6*lambda**2)
        else if ( il.eq.63 ) then  ! jlm by extrap of above coeffs
           omega = 1 + 0.95505 * (1 + 0.010215 *lambda) /
     &                  (1 + 0.015466*lambda + 4.33e-6*lambda**2)
        else if ( il.eq.157 ) then  ! still to do this
           omega = 1 + 0.97966 * (1 + 0.005798 *lambda) /
     &                  (1 + 0.00735*lambda + 5.200e-7*lambda**2)
        else
           print*, ' Error: optm not set for il = ', il
           stop
        endif     !  ( il.eq.20 .or. il.eq.21 )
      elseif(npanels.eq.13)then
        if ( il.eq.13 ) then   ! 6*13=78 approx equals 80=4*20
           omega = 1 + 0.8715 * (1 + 0.0291874*lambda) /
     &               (1 + 0.0540404*lambda + 0.000153109*lambda**2)
        else if ( il.eq.20 ) then  ! 6*20=120 equals 120=4*30
           omega = 1 + 0.9083 * (1 + 0.019229*lambda) /
     &               (1 + 0.0340753*lambda + 4.1772e-05*lambda**2)
        else if ( il.eq.25 ) then  ! 6*25=150 approx equals 148=4*37
           omega = 1 + 0.9250 * (1 + 0.0152456*lambda) /
     &                  (1 + 0.0261757*lambda + 1.96856e-05*lambda**2)
        else if ( il.eq.32 ) then
!  John's interpolated version
!          omega = 1 + 0.934  * (1 + 0.012    *lambda) /
!    &                  (1 + 0.021    *lambda + 1.7    e-05*lambda**2)
           omega = 1 + 0.94239 * (1 + 0.01129 *lambda) /
     &                  (1 + 0.01853*lambda + 7.6e-06*lambda**2)
        else
           print*, ' Error: optm not set for il = ', il
           stop
        endif     !  ( il.eq.13 )
      endif       !  (npanels.eq.13)

      return
      end
