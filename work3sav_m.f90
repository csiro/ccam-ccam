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
public qgsav,qfgsav,qlgsav,trsav
public qrgsav,qsngsav,qgrgsav
public work3sav_init,work3sav_end

real, dimension(:,:), allocatable, save :: qgsav,qfgsav,qlgsav
real, dimension(:,:), allocatable, save :: qrgsav,qsngsav,qgrgsav
real, dimension(:,:,:), allocatable, save :: trsav

contains

subroutine work3sav_init(ifull,iextra,kl,ilt,jlt,klt,ngasmax)

implicit none

integer, intent(in) :: ifull,iextra,kl,ilt,jlt,klt,ngasmax

allocate(qgsav(ifull,kl),qfgsav(ifull,kl),qlgsav(ifull,kl))
allocate(qrgsav(ifull,kl))
allocate(qsngsav(ifull,kl),qgrgsav(ifull,kl))
if (ilt>0) allocate(trsav(ilt*jlt,klt,ngasmax))

return
end subroutine work3sav_init

subroutine work3sav_end

implicit none

deallocate(qgsav,qfgsav,qlgsav)
deallocate(qrgsav)
deallocate(qsngsav,qgrgsav)
if (allocated(trsav)) deallocate(trsav)

return
end subroutine work3sav_end

end module work3sav_m