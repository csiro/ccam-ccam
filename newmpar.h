!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
      integer :: nproc            ! Number of processors to use
      integer :: kl               ! Vertical levels
      integer :: ol               ! Ocean levels
      integer, parameter :: ms=6  ! Levels in surface scheme
      integer, parameter :: npanels = 5
      integer :: il_g
      integer :: nrows_rad        ! usually 8, but 6 for C63/3
      integer :: jl_g
      integer :: ifull_g, ijk_g
      ! Note that iquad is only used globally
      integer :: iquad
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

      integer :: nxp,nyp
      integer :: il, jl
      integer :: npan
      integer :: ifull, ijk
      integer :: iextra

! eak addition 16/03/06  ! MJT cable
      integer, parameter :: mxst=13       ! max_no_of_soil_types
      integer, parameter :: mxvt=17       ! max_no_of_vegetation_types

      common/newmpar/nproc,kl,ol,il_g,jl_g,ifull_g,nrows_rad,ijk_g,       &
     &               iquad,nxp,nyp,il,jl,ifull,npan,ijk,iextra