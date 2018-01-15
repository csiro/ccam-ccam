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
    
module unn_m

implicit none

private
public unn,vnn
public unn_init,unn_end

real, dimension(:,:), allocatable, save :: unn,vnn

contains

subroutine unn_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(unn(ifull,kl),vnn(ifull,kl))

return
end subroutine unn_init

subroutine unn_end

implicit none

deallocate(unn,vnn)

return
end subroutine unn_end

end module unn_m