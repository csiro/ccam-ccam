! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------


!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
      integer :: nproc                  ! Number of processors to use
      integer :: kl                     ! Atmosphere vertical levels
      integer :: ol                     ! Ocean vertical levels
      integer, parameter :: ms = 6      ! Soil levels in surface scheme

      integer :: il_g                   ! Global grid size in X-dir
      integer :: jl_g                   ! Global grid size in Y-dir
      integer :: ifull_g                ! Number of global grid points
      integer :: nrows_rad              ! Subset of grid for radiation
      integer, parameter :: npanels = 5 ! Cubic panels (0-5)
      integer :: iquad                  ! iquad is only used globally
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

      integer :: nxp, nyp               ! Number of processors for decompostion
      integer :: il, jl                 ! Local processor grid size
      integer :: npan                   ! Number of panels for processor
      integer :: ifull                  ! Number of grid points for processor
      integer :: iextra                 ! Size of halo for processor

      integer, parameter :: mxst = 13   ! max_no_of_soil_types
      integer, parameter :: mxvt = 17   ! max_no_of_vegetation_types

      common/newmpar/nproc,kl,ol,il_g,jl_g,ifull_g,nrows_rad,iquad,nxp,   &
     &               nyp,il,jl,ifull,npan,iextra