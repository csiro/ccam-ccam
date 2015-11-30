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
    
module xyzinfo_m

implicit none

private
public x_g, y_g, z_g, wts_g
public x, y, z, wts
#ifdef usempi3
public x_g_win, y_g_win, z_g_win
#endif
public xyzinfo_init, xyzinfo_end


real, dimension(:), allocatable, save :: wts
real(kind=8), dimension(:), allocatable, save :: x, y, z
real, dimension(:), allocatable, save :: wts_g
#ifdef usempi3
real(kind=8), dimension(:), pointer, save :: x_g, y_g, z_g
integer, save :: x_g_win, y_g_win, z_g_win
#else
real(kind=8), dimension(:), allocatable, save :: x_g, y_g, z_g
#endif

contains

subroutine xyzinfo_init(ifull_g,ifull,iextra,myid)

implicit none

integer, intent(in) :: ifull_g, ifull, iextra, myid

if ( myid==0 ) then
  allocate( wts_g(ifull_g) )  
end if
allocate( x(ifull), y(ifull), z(ifull) )
allocate( wts(ifull))

return
end subroutine xyzinfo_init

subroutine xyzinfo_end

implicit none

if ( allocated(wts_g) ) then
  deallocate( wts_g )
end if
deallocate( x, y, z )
deallocate( wts )

return
end subroutine xyzinfo_end

end module xyzinfo_m
