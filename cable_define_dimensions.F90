
module define_dimensions
   implicit none
   
   !---CABLE default KINDs for representing INTEGER/REAL values   
   integer, parameter :: i_d = KIND(9)
   integer, parameter :: r_1  = KIND(1.0) 
   !---at least 10-digit precision
   integer, parameter :: r_2  = SELECTED_REAL_KIND(12, 50)

   integer, parameter :: n_tiles = 17        ! # possible no of different tiles/patches
   integer :: mp                             ! # total no of patches/tiles 
   integer, parameter :: ncp = 3             ! # vegetation carbon stores
   integer, parameter :: ncs = 2             ! # soil carbon stores
   integer, parameter :: mf = 2              ! # leaves (sunlit, shaded)
   integer, parameter :: nrb = 3             ! # radiation bands
   integer, parameter :: msn = 3              ! max # snow layers
   integer, parameter :: swb = 2             ! # shortwave bands 
   INTEGER, PARAMETER :: ms = 6  ! # soil layers

   integer            :: mvtype               ! total # vegetation types,   from input
!  integer            :: mvtype               ! total # vegetation types,   from input
   integer            :: mstype               ! total # soil types,         from input

   integer :: mland                           ! # land grid cells

end module define_dimensions
