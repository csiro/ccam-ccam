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
      integer :: nproc            ! Number of processors to use
      integer :: kl               ! Vertical levels
      integer :: ol               ! Ocean levels
      integer, parameter :: ms=6  ! Levels in surface scheme
      integer, parameter :: npanels = 5
      integer :: il_g
      integer :: nrows_rad        ! usually 8, but 6 for C63/3
      integer :: jl_g
      integer :: ifull_g
      ! Note that iquad is only used globally
      integer :: iquad
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

      integer :: nxp,nyp
      integer :: il, jl
      integer :: npan
      integer :: ifull
      integer :: iextra

      integer, parameter :: mxst=13       ! max_no_of_soil_types
      integer, parameter :: mxvt=17       ! max_no_of_vegetation_types

      common/newmpar/nproc,kl,ol,il_g,jl_g,ifull_g,nrows_rad,iquad,nxp,   &
     &               nyp,il,jl,ifull,npan,iextra