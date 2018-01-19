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
    
module pbl_m

implicit none

private
public cduv,cdtq,tss,slwa
public pbl_init,pbl_end

real, dimension(:), allocatable, save :: slwa
real, dimension(:), allocatable, save :: tss, cduv, cdtq

contains

subroutine pbl_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(cduv(ifull),cdtq(ifull),tss(ifull),slwa(ifull))
cduv = 0.
cdtq = 0.
tss = 300.
slwa = 0.

return
end subroutine pbl_init

subroutine pbl_end

implicit none

deallocate(cduv,cdtq,tss,slwa)

return
end subroutine pbl_end

end module pbl_m