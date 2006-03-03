!     Include file for AMIP2 ozone

!     Resolution of dataset
      integer, parameter :: mo=12    ! Months
      integer, parameter :: jg=64    ! Latitudes
      integer, parameter :: kg=59    ! Levels
      integer, parameter :: lg=kg+1  ! Layer Interfaces
      
      real, dimension(jg)       :: glat
      real, dimension(kg)       :: dp
      real, dimension(lg)       :: gpri
      real, dimension(jg,kg,mo) :: gdat

      common / o3amip2 / gdat, glat, dp, gpri
