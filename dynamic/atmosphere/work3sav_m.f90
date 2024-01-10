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
    
module work3sav_m

implicit none

private
public tsav, qgsav, qfgsav, qlgsav, trsav
public work3sav_init, work3sav_end

real, dimension(:,:), allocatable, save :: tsav
real, dimension(:,:), allocatable, save :: qgsav, qfgsav, qlgsav
real, dimension(:,:,:), allocatable, save :: trsav

contains

subroutine work3sav_init(ifull,kl,ngas)

implicit none

integer, intent(in) :: ifull,kl,ngas

allocate( tsav(ifull,kl) )
allocate( qgsav(ifull,kl), qfgsav(ifull,kl), qlgsav(ifull,kl) )
if ( ngas>0 ) then
  allocate( trsav(ifull,kl,ngas) )
end if

return
end subroutine work3sav_init

subroutine work3sav_end

implicit none

deallocate( tsav )
deallocate( qgsav, qfgsav, qlgsav )
if ( allocated(trsav) ) then
  deallocate( trsav )
end if

return
end subroutine work3sav_end

end module work3sav_m