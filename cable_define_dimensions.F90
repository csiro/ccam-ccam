!!$ cable_variables.f90
!!$
!!$ This file declares all non-local variables for CABLE, 
!!$ CSIRO land surface model
!!$
!!$ Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
!!$ Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!!$
!!$ Fortran-95 coding by Harvey Davies, Gab Abramowitz and Martin Dix
!!$ bugs to gabsun@gmail.com.

!=========================================================================
MODULE define_dimensions
!les  INTEGER 	        :: mg      ! # no of grid land points
  INTEGER 	        :: mp      ! # total no of patches/tiles 
  INTEGER, PARAMETER	:: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER	:: nrb = 3 ! # radiation bands
!  INTEGER, PARAMETER	:: ms = 4  ! # soil layers
  INTEGER, PARAMETER	:: ms = 6  ! # soil layers
  INTEGER, PARAMETER	:: ncp = 3 ! # vegetation carbon stores
  INTEGER, PARAMETER	:: ncs = 2 ! # soil carbon stores
  INTEGER, PARAMETER	:: niter = 4 ! # soil carbon stores
  ! i_d is default kind for representing integer values.
  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)
END MODULE define_dimensions
