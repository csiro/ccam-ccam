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
    
module kuocomb_m

implicit none

private
public kbsav,ktsav
public convpsav,fluxtot
public kuocomb_init,kuocomb_end

integer, dimension(:), allocatable, target, save :: kbsav, ktsav
real, dimension(:,:), allocatable, save :: fluxtot
real, dimension(:), allocatable, save :: convpsav

contains

subroutine kuocomb_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(kbsav(ifull),ktsav(ifull))
allocate(convpsav(ifull),fluxtot(ifull,kl))
kbsav=kl-1
ktsav=kl-1
convpsav=0.
fluxtot=0.

return
end subroutine kuocomb_init

subroutine kuocomb_end

implicit none

deallocate(kbsav,ktsav)
deallocate(convpsav,fluxtot)

return
end subroutine kuocomb_end

end module kuocomb_m