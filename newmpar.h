!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
      integer, parameter :: nproc = 1  ! Number of processors to use

      integer, parameter :: kl=18  ! Vertical levels
      integer, parameter :: ms=6  ! Levels in surface scheme

      integer, parameter :: npanels = 5
      integer, parameter :: il_g = 48
      integer, parameter :: nrows_rad = 8    ! usually 8, but 6 for C63/3
      integer, parameter :: jl_g = il_g + npanels*il_g
      integer, parameter :: ifull_g = il_g*jl_g, ijk_g = il_g*jl_g*kl
      ! Note that iquad is only used globally
      integer, parameter :: iquad=1+il_g*((8*npanels)/(npanels+4))
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

! eak addition 16/03/06
      integer, parameter :: mxst=13       ! max_no_of_soil_types
      integer, parameter :: mxvt=13       ! max_no_of_vegetation_types
      integer, parameter :: ncp=3         ! no of vegetation carbon pools
      integer, parameter :: ncs=3         ! no of soil carbon pools

      integer, parameter :: nprocmax = 96
!     This array defines the split up of processors. Zero values are 
!     for values of nproc that won't work.
      integer, parameter, dimension(nprocmax) :: nxp = (/                &
     &         1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,                       & 
     &         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                       &
     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4    /)    
      integer, parameter :: nyp = nproc/nxp(nproc)

      integer, parameter :: il=il_g/nxp(nproc), jl=jl_g/nyp

      integer, parameter :: npan=max(1,(npanels+1)/nproc)
      integer, parameter :: ifull = il*jl, ijk = il*jl*kl

!     The perimeter of the processor region has length 2*(il+jl).
!     The first row has 8 possible corner points per panel and the 
!     second has 16. In practice these are not all distinct so there could
!     be some optimisation.
      integer, parameter :: iextra = 4*(il+jl)+24*npan
