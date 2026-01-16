! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public sdot, dpsldt, wvel
public updraft_helicity
public vvel_init, vvel_end

real, dimension(:,:), allocatable, save :: sdot,dpsldt
real, dimension(:,:), allocatable, save :: wvel
real, dimension(:), allocatable, save :: updraft_helicity

contains

subroutine vvel_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate( sdot(ifull,kl+1), dpsldt(ifull,kl) )
allocate( wvel(ifull+iextra,kl) )
allocate( updraft_helicity(ifull) )
sdot = 0.
dpsldt = 0.
wvel = 0.
updraft_helicity = 0.

return
end subroutine vvel_init

subroutine vvel_end

implicit none

deallocate( sdot, dpsldt )
deallocate( wvel )
deallocate( updraft_helicity )

return
end subroutine vvel_end

end module vvel_m