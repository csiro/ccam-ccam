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
    
module bigxy4_m

implicit none

private
public xx4,yy4
public bigxy4_init,bigxy4_end

real(kind=8), dimension(:,:), allocatable, save :: xx4,yy4

contains

subroutine bigxy4_init(iquad)

implicit none

integer, intent(in) :: iquad

allocate(xx4(iquad,iquad),yy4(iquad,iquad))

return
end subroutine bigxy4_init

subroutine bigxy4_end

implicit none

deallocate(xx4,yy4)

return
end subroutine bigxy4_end

end module bigxy4_m