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
    
module vecsuv_m

implicit none

private
public ax_g,bx_g,ay_g,by_g,az_g,bz_g
public ax,bx,ay,by,az,bz
public vecsuv_init,vecsuv_end

real, dimension(:), allocatable, save :: ax_g,bx_g,ay_g,by_g,az_g,bz_g
real, dimension(:), allocatable, save :: ax,bx,ay,by,az,bz

contains

subroutine vecsuv_init(ifull_g,ifull,iextra,myid)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,myid

if (myid==0) then
  allocate(ax_g(ifull_g),bx_g(ifull_g))
  allocate(ay_g(ifull_g),by_g(ifull_g))
  allocate(az_g(ifull_g),bz_g(ifull_g))
end if
allocate(ax(ifull+iextra),bx(ifull+iextra))
allocate(ay(ifull+iextra),by(ifull+iextra))
allocate(az(ifull+iextra),bz(ifull+iextra))

return
end subroutine vecsuv_init

subroutine vecsuv_end

implicit none

if (allocated(ax_g)) deallocate(ax_g)
if (allocated(bx_g)) deallocate(bx_g)
if (allocated(ay_g)) deallocate(ay_g)
if (allocated(by_g)) deallocate(by_g)
if (allocated(az_g)) deallocate(az_g)
if (allocated(bz_g)) deallocate(bz_g)
deallocate(ax,bx,ay,by,az,bz)

return
end subroutine vecsuv_end

end module vecsuv_m