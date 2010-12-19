!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
      integer :: nproc ! Number of processors to use

      integer :: kl ! Vertical levels
      integer, parameter :: ms=6  ! Levels in surface scheme

      integer, parameter :: npanels = 5
      integer :: il_g
      integer :: nrows_rad    ! usually 8, but 6 for C63/3
      integer :: jl_g
      integer :: ifull_g, ijk_g
      ! Note that iquad is only used globally
      integer :: iquad
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

! eak addition 16/03/06  ! MJT cable
      integer, parameter :: mxst=13       ! max_no_of_soil_types
      integer, parameter :: mxvt=17       ! max_no_of_vegetation_types

      integer, parameter :: nprocmax = 384
!     This array defines the split up of processors. Zero values are 
!     for values of nproc that won't work.
      integer, parameter, dimension(nprocmax) :: nxp = (/                &
     &         1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,                       & 
     &         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 4,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,                       &
     &         0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,                       &
     &         0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 6,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,                       &
     &         0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 5,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8 /)    
      integer :: nyp

      integer :: il, jl

      integer :: npan
      integer :: ifull, ijk

!     The perimeter of the processor region has length 2*(il+jl).
!     The first row has 8 possible corner points per panel and the 
!     second has 16. In practice these are not all distinct so there could
!     be some optimisation.
      integer :: iextra

      common/newmpar/nproc,kl,il_g,jl_g,ifull_g,nrows_rad,ijk_g,         &
     &               iquad,nyp,il,jl,ifull,npan,ijk,iextra