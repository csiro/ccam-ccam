! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module mlodynamicsarrays_m

implicit none

private
public ee, eeu, eev
public dd, ddu, ddv
public stwgtu, stwgtv
public gosig, gosigh, godsig
public godsigu, godsigv, gosighu, gosighv
public oldu1, oldu2, oldv1, oldv2
public ipice
public mlodynamicsarrays_init, mlodynamicsarrays_end

real, dimension(:,:), allocatable, save :: ee, eeu, eev
real, dimension(:), allocatable, save :: dd, ddu, ddv
real, dimension(:,:), allocatable, save :: stwgtu, stwgtv
real, dimension(:,:), allocatable, save :: gosig, gosigh, godsig
real, dimension(:,:), allocatable, save :: godsigu, godsigv, gosighu, gosighv
real, dimension(:,:), allocatable, save :: oldu1, oldu2, oldv1, oldv2
real, dimension(:), allocatable, save :: ipice

contains

subroutine mlodynamicsarrays_init(ifull,iextra,wlev)

implicit none

integer, intent(in) :: ifull, iextra, wlev

allocate(oldu1(ifull,wlev),oldv1(ifull,wlev))
allocate(oldu2(ifull,wlev),oldv2(ifull,wlev))
allocate(ipice(ifull+iextra))
oldu1(:,:)=0.
oldv1(:,:)=0.
oldu2(:,:)=0.
oldv2(:,:)=0.
ipice(:)=0.
  
return
end subroutine mlodynamicsarrays_init

subroutine mlodynamicsarrays_end

implicit none

if ( allocated( oldu1 ) ) then
  deallocate( oldu1, oldu2, oldv1, oldv2 )
  deallocate( ipice )
end if

return
end subroutine mlodynamicsarrays_end

end module mlodynamicsarrays_m