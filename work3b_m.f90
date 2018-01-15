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
    
module work3b_m

implicit none

private
public wblf,wbfice,sdepth
public work3b_init,work3b_end

real, dimension(:,:), allocatable, save :: wblf,wbfice,sdepth

contains

subroutine work3b_init(ifull,ms)

implicit none

integer, intent(in) :: ifull,ms

allocate(wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3))

return
end subroutine work3b_init

subroutine work3b_end

implicit none

deallocate(wblf,wbfice,sdepth)

return
end subroutine work3b_end

end module work3b_m