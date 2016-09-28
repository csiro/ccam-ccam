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
    
module vvel_m

implicit none

private
public sdot,dpsldt
public vvel_init,vvel_end

real, dimension(:,:), allocatable, save :: sdot,dpsldt

contains

subroutine vvel_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(sdot(ifull,kl+1),dpsldt(ifull,kl))
sdot=0.
dpsldt=0.

return
end subroutine vvel_init

subroutine vvel_end

implicit none

deallocate(sdot,dpsldt)

return
end subroutine vvel_end

end module vvel_m