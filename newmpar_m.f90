! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module newmpar_m

implicit none

private
public nproc, kl, ol, ms
public il_g, jl_g, ifull_g, nrows_rad, npanels, iquad
public nxp, nyp, il, jl, npan, ifull, iextra
public mxst, mxvt

!     This version is for the MPI code. Variables with suffix _g
!     are truly global, others refer to a processor's own region.
integer, save :: nproc                  ! Number of processors to use
integer, save :: kl                     ! Atmosphere vertical levels
integer, save :: ol                     ! Ocean vertical levels
integer, parameter :: ms = 6      ! Soil levels in surface scheme
!$acc declare create(kl)

integer, save :: il_g                   ! Global grid size in X-dir
integer, save :: jl_g                   ! Global grid size in Y-dir
integer, save :: ifull_g                ! Number of global grid points
integer, save :: nrows_rad              ! Subset of grid for radiation
integer, parameter :: npanels = 5 ! Cubic panels (0-5)
integer, save :: iquad                  ! iquad is only used globally
!$acc declare create(ifull_g)
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1

integer, save :: nxp, nyp               ! Number of processors for decompostion
integer, save :: il, jl                 ! Local processor grid size
integer, save :: npan                   ! Number of panels for processor
integer, save :: ifull                  ! Number of grid points for processor
integer, save :: iextra                 ! Size of halo for processor
!$acc declare create(ifull,iextra)

integer, parameter :: mxst = 13   ! max_no_of_soil_types
integer, parameter :: mxvt = 17   ! max_no_of_vegetation_types

end module newmpar_m