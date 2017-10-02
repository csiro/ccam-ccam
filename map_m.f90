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
    
module map_m

implicit none

private
public em_g, emu_g, emv_g, f_g, fu_g, fv_g
public dmdx_g, dmdy_g
public em, emu, emv, f, fu, fv
public em_g_win
public em_g_dummy
public map_init, map_end

real, dimension(:), allocatable, save :: emu_g, emv_g, f_g, fu_g, fv_g
real, dimension(:), allocatable, save :: dmdx_g, dmdy_g
real, dimension(:), allocatable, save :: emu, emv, fu, fv
real, dimension(:), allocatable, save :: em, f
real, dimension(:), pointer, save :: em_g
real, dimension(:), allocatable, target, save :: em_g_dummy
integer, save :: em_g_win


contains

subroutine map_init(ifull_g,ifull,iextra,myid)

implicit none

integer, intent(in) :: ifull_g, ifull, iextra, myid

if ( myid==0 ) then
  allocate( emu_g(ifull_g), emv_g(ifull_g) )
  allocate( f_g(ifull_g), fu_g(ifull_g), fv_g(ifull_g) )
  allocate( dmdx_g(ifull_g), dmdy_g(ifull_g) )
end if
allocate( em(ifull+iextra), emu(ifull+iextra), emv(ifull+iextra) )
allocate( f(ifull+iextra), fu(ifull+iextra), fv(ifull+iextra) )

return
end subroutine map_init

subroutine map_end

implicit none

if ( allocated(emu_g) ) then
  deallocate( emu_g, emv_g )
  deallocate( f_g, fu_g, fv_g )
  deallocate( dmdx_g, dmdy_g )
end if
deallocate( em, emu, emv )
deallocate( f, fu, fv )

return
end subroutine map_end

end module map_m