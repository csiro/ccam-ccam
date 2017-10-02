! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module latlong_m

implicit none

private
public rlatt_g,rlongg_g
public rlatt, rlongg
public latlong_init,latlong_end

real, dimension(:), allocatable, save :: rlatt_g,rlongg_g
real, dimension(:), allocatable, save :: rlatt, rlongg

contains

subroutine latlong_init(ifull_g,ifull,myid)

implicit none

integer, intent(in) :: ifull_g,ifull,myid

if (myid==0) then
  allocate(rlatt_g(ifull_g),rlongg_g(ifull_g))
end if
allocate(rlatt(ifull),rlongg(ifull))

return
end subroutine latlong_init

subroutine latlong_end

implicit none

if (allocated(rlatt_g)) deallocate(rlatt_g)
if (allocated(rlongg_g)) deallocate(rlongg_g)
deallocate(rlatt,rlongg)

return
end subroutine latlong_end

end module latlong_m