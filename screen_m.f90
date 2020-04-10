! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module screen_m

implicit none

private
public tscrn,qgscrn,uscrn,rhscrn,u10
public screen_init,screen_end

real, dimension(:), allocatable, save :: uscrn, rhscrn
real, dimension(:), allocatable, save :: tscrn, qgscrn, u10

contains

subroutine screen_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(tscrn(ifull),qgscrn(ifull),uscrn(ifull),rhscrn(ifull),u10(ifull))

tscrn = 0.
qgscrn = 0.
uscrn = 0.
rhscrn = 0.
u10 = 0.

return
end subroutine screen_init

subroutine screen_end

implicit none

deallocate(tscrn,qgscrn,uscrn,rhscrn,u10)

return
end subroutine screen_end

end module screen_m