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
    
module utilities
   implicit none
contains
   pure function calc_rotpole(rlong0,rlat0) result (rotpole)
!     Coefficients for rotation of pole
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
      use const_phys  ! Physical constants
      real, intent(in) :: rlong0, rlat0
      real, dimension(3,3) :: rotpole
      real ::  coslong, sinlong, coslat, sinlat
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      rotpole(1,1)=coslong*sinlat
      rotpole(1,2)=-sinlong
      rotpole(1,3)=coslong*coslat
      rotpole(2,1)=sinlong*sinlat
      rotpole(2,2)=coslong
      rotpole(2,3)=sinlong*coslat
      rotpole(3,1)=-coslat
      rotpole(3,2)=0.
      rotpole(3,3)=sinlat
   end function calc_rotpole
end module utilities

