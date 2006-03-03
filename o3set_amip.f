      subroutine o3set_amip ( alat, npts, mins, sigh, ps, qo3 )

!     Interpolate the AMIP2 ozone data in time, latitude and pressure to 
!     produce ozone mixing ratio for model radiation code. 

!     In CCAM version latitude may be different for every point.

!     Fixed format f90

      implicit none
      include 'newmpar.h'
      integer, parameter :: klp = kl+1
      include 'o3amip.h'
      include 'const_phys.h'
      integer, intent(in) :: npts
      real, dimension(npts), intent(in) :: alat     ! Latitude
      integer, intent(in) :: mins  ! Time
      real, intent(in),  dimension(klp) :: sigh ! Half level sigma
      real, intent(in),  dimension(npts) :: ps  ! Surface pressure
      real, intent(out), dimension(npts,kl) :: qo3  ! Ozone mixing ratio
!     Day no of middle of month
      real, parameter, dimension(12) :: monmid =  
     &  (/ 15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 
     &     288.5, 319.0, 349.5 /)
      real, dimension(kg) :: oz    ! Column for this date and lat.
      real, dimension(lg) :: ozcol ! Integrated column
      real, dimension(klp) :: qo3p
      integer :: k1
      integer :: j, m1, m2, j1, j2, k, kk, k1min
      real, dimension(klp) :: prh ! Half level pressure
      real :: fac1, fac2, tfac1, tfac2, date, theta

!     Time interpolation factors (assume year of exactly 365 days)
      date = real(modulo(mins,525600)) / 1440.0
      if ( date <= monmid(1) ) then ! First half of Jan
         m1 = 12
         m2 = 1
         tfac1 = (monmid(1) - date)/31.0
         tfac2 = 1.0 - tfac1
      else if ( date >= monmid(12) ) then
         m1 = 12
         m2 = 1
         tfac2 = (date - monmid(12))/31.0
         tfac1 = 1 - tfac2
      else
!        Search for bracketing dates, i such that monmid(i) >= date
         do m2=2,12
            if ( monmid(m2) >= date ) exit
         end do
         m1 = m2 - 1
         tfac1 = ( monmid(m2) - date ) / ( monmid(m2) - monmid(m1) )
         tfac2 = 1.0 - tfac1
      end if


      do j=1,npts
!        Factors for interpolation in latitude, 
!        Note that ozone grid latitudes run S to N.
         theta = alat(j)*180.0/pi
         if ( theta <= glat(1) ) then
            j1 = 1
            j2 = 1
            fac1 = 1.0
            fac2 = 0.0
         else if ( theta >= glat(jg) ) then
            j1 = jg
            j2 = jg
            fac1 = 1.0
            fac2 = 0.0
         else
!           Input data isn't exactly equally spaced, so search to find the spanning
!           latitudes. Want to find first j such that glat(j) >= theta
            j2 = 1 + count ( glat < theta )
!            do j2=2,jg
!               if ( glat(j2) >= theta ) exit
!            end do
            j1 = j2 - 1
            fac1 = ( glat(j2) - theta ) / ( glat(j2) - glat(j1) )
            fac2 = 1.0 - fac1
         end if

!     Now interpolate in latitude and time.
         oz = tfac1 * ( fac1*gdat(j1,:,m1) + fac2*gdat(j2,:,m1 ) ) + 
     &        tfac2 * ( fac1*gdat(j1,:,m2) + fac2*gdat(j2,:,m2 ) ) 

!     Each longitude has same vertical profile as a function of pressure,
!     but model levels depend on surface pressure so each must be 
!     interpolated separately. To ensure the column amount is maintained
!     correctly, first calculate ozone column between every half level and
!     TOA.
         ozcol(1) = 0.0
         do k=1,kg
            ozcol(k+1) = ozcol(k) + oz(k)*dp(k)
         end do

!     To calculate model amounts, find the data half level immediately above
!     the given model half level. Note that model levels start at surface while
!     data levels start at the top.
!
         prh = sigh * ps(j)*0.01 ! Input PS in Pa

         qo3p = 0.
         qo3p(klp) = 0.0      ! TOA
         do k=kl,1,-1           ! Start at the top model level

!        Find largest k1 such that gpri(k1) <= prh(k)
!        Restrict to range 1:kg so as not to go out of bounds of oz
            k1 = max(1, count( gpri(1:kg) <= prh(k) ) )
            qo3p(k) = ozcol(k1) + (prh(k)-gpri(k1))*oz(k1)
         end do

!     Finally get the model ozone mixing ratios by differencing the half level
!     path amounts. The vertical order of the qo3 array is reversed to 
!     match the radiation code.
!         if ( j==1 ) print*, 'QO3P', qo3p
         do k=1,kl
            qo3(j,kl+1-k) = ( qo3p(k+1) - qo3p(k) ) / 
     &                      ( prh(k+1) - prh(k) )
         end do
      end do

      end subroutine o3set_amip
