!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
      integer, parameter :: nproc = 2  ! Number of processors to use

      integer, parameter :: kl=18  ! Vertical levels
      integer, parameter :: ms=6  ! Levels in surface scheme

      integer, parameter :: npanels = 5
      integer, parameter :: il_g = 48
      integer, parameter :: jl_g = il_g + npanels*il_g
      integer, parameter :: ifull_g = il_g*jl_g, ijk_g = il_g*jl_g*kl
      ! Note that iquad is only used globally
      integer, parameter :: iquad=1+il_g*((8*npanels)/(npanels+4))
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

      ! For 1 - 6 processors
      integer, parameter :: il = il_g, jl = jl_g/nproc
      integer, parameter :: ifull = il*jl, ijk = il*jl*kl
      ! For 1-6
      integer, parameter :: iextra = (4*(npanels+1)*(il_g+il_g+8))/nproc
      integer, parameter :: npan=(npanels+1)/nproc
      !  iextra = 4*npan*(ipan+jpan+8)

