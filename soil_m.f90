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
    
module soil_m

implicit none

private
public zoland,zmin,zolnd,zolog,albvissav,so4t,albnirsav
public albvisdif,albnirdif,albvisdir,albnirdir
public land
public soil_init,soil_end

real, dimension(:), allocatable, save :: zolnd,zolog,albvissav,albnirsav
real, dimension(:), allocatable, save :: albvisdif,albnirdif,albvisdir,albnirdir
real, dimension(:), allocatable, save :: so4t
real, save :: zoland,zmin
logical, dimension(:), allocatable, save :: land

contains

subroutine soil_init(ifull,iaero,nsib)

implicit none

integer, intent(in) :: ifull,iaero,nsib

allocate(zolnd(ifull),albvissav(ifull),albnirsav(ifull))
allocate(albvisdif(ifull),albnirdif(ifull),albvisdir(ifull),albnirdir(ifull))
allocate(land(ifull))
if (iaero/=0) then
  allocate(so4t(ifull))
end if
if (nsib==3.or.nsib==5) then
  allocate(zolog(ifull))
end if

zoland=0.16
albvissav = 0.
albnirsav = 0.
albvisdif = 0.
albnirdif = 0.
albvisdir = 0.
albnirdir = 0.

return
end subroutine soil_init

subroutine soil_end

implicit none

deallocate(zolnd,zolog,albvissav,albnirsav)
deallocate(albvisdif,albnirdif,albvisdir,albnirdir)
deallocate(land)
if (allocated(so4t)) deallocate(so4t)
if (allocated(zolog)) deallocate(zolog)

return
end subroutine soil_end

end module soil_m
