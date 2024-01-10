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
    
module permsurf_m

implicit none

private
public ipsice,ipsea,ipland,iperm
public permsurf_init,permsurf_end

integer, save :: ipsice,ipsea,ipland
integer, dimension(:), allocatable, save :: iperm

contains

subroutine permsurf_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(iperm(ifull))

return
end subroutine permsurf_init

subroutine permsurf_end

implicit none

deallocate(iperm)

return
end subroutine permsurf_end

end module permsurf_m