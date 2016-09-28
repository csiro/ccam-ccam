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
    
module vegpar_m

implicit none

private
public cansto,vlai,fwet
public vegpar_init,vegpar_end

real, dimension(:), allocatable, save :: cansto,vlai,fwet

contains

subroutine vegpar_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(cansto(ifull),vlai(ifull),fwet(ifull))
cansto=0.
vlai=0.
fwet=0.

return
end subroutine vegpar_init

subroutine vegpar_end

implicit none

deallocate(cansto,vlai,fwet)

return
end subroutine vegpar_end

end module vegpar_m