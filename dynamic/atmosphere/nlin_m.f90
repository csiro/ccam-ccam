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
    
module nlin_m

implicit none

private
public tn,un,vn
public nlin_init,nlin_end

real, dimension(:,:), allocatable, save :: tn,un,vn

contains

subroutine nlin_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(tn(ifull,kl),un(ifull,kl),vn(ifull,kl))

return
end subroutine nlin_init

subroutine nlin_end

implicit none

deallocate(tn,un,vn)

return
end subroutine nlin_end

end module nlin_m