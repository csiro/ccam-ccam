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
    
module srccom_m

implicit none

private
public sorc,csour1,csour2,osour,csour,ss1
public srccom_init,srccom_end

real, dimension(:,:), allocatable, save :: csour1,csour2,osour,csour,ss1
real, dimension(:,:,:), allocatable, save :: sorc

contains

subroutine srccom_init(ifull,iextra,kl,imax,nbly)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax,nbly

allocate(sorc(imax,kl+1,nbly),csour1(imax,kl+1),csour2(imax,kl+1))
allocate(osour(imax,kl+1),csour(imax,kl+1),ss1(imax,kl+1))

return
end subroutine srccom_init

subroutine srccom_end

implicit none

deallocate(sorc,csour1,csour2,osour,csour,ss1)

return
end subroutine srccom_end

end module srccom_m