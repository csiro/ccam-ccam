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
    
module work3f_m

implicit none

private
public qccon,qlrad,qfrad
public nface,xg,yg
public work3f_init,work3f_end

real, dimension(:,:), allocatable, save :: qccon,qlrad,qfrad
real, dimension(:,:), allocatable, save :: xg,yg
integer, dimension(:,:), allocatable, save :: nface

contains

subroutine work3f_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl))
allocate(nface(ifull,kl),xg(ifull,kl),yg(ifull,kl))

qlrad = 0.
qfrad = 0.

return
end subroutine work3f_init

subroutine work3f_end

implicit none

deallocate(qccon,qlrad,qfrad)
deallocate(nface,xg,yg)

return
end subroutine work3f_end

end module work3f_m