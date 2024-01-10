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
    
module swocom_m

implicit none

private
public fsw,dfsw,ufsw,hsw
public swocom_init,swocom_end

real, dimension(:,:), allocatable, save :: fsw,dfsw,ufsw,hsw

contains

subroutine swocom_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(fsw(imax,kl+1),dfsw(imax,kl+1),ufsw(imax,kl+1),hsw(imax,kl))

return
end subroutine swocom_init

subroutine swocom_end

implicit none

deallocate(fsw,dfsw,ufsw,hsw)

return
end subroutine swocom_end

end module swocom_m