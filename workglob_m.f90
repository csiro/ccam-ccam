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
    
module workglob_m

implicit none

private
public rlong4,rlat4
public rlong4_l,rlat4_l
public workglob_init,workglob_end

real, dimension(:,:), allocatable, save :: rlong4,rlat4
real, dimension(:,:), allocatable, save :: rlong4_l,rlat4_l

contains

subroutine workglob_init(ifull_g,ifull,myid)

implicit none

integer, intent(in) :: ifull_g,ifull,myid

if ( myid==0 ) then
  if (.not.allocated(rlong4)) then
    allocate(rlong4(ifull_g,4),rlat4(ifull_g,4))
  end if
end if
if (.not.allocated(rlong4_l)) then
  allocate(rlong4_l(ifull,4),rlat4_l(ifull,4))
end if

return
end subroutine workglob_init

subroutine workglob_end

implicit none

if ( allocated(rlong4) ) then
  deallocate(rlong4)
end if
if ( allocated(rlat4) ) then
  deallocate(rlat4)
end if
deallocate(rlong4_l,rlat4_l)

return
end subroutine workglob_end

end module workglob_m