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
    
module xarrs_m

implicit none

private
public ux,vx,tx,pslx
public xarrs_init,xarrs_end

real, dimension(:,:), allocatable, save :: ux,vx,tx,pslx

contains

subroutine xarrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(ux(ifull,kl),vx(ifull,kl),tx(ifull+iextra,kl),pslx(ifull+iextra,kl))

return
end subroutine xarrs_init

subroutine xarrs_end

implicit none

deallocate(ux,vx,tx,pslx)

return
end subroutine xarrs_end

end module xarrs_m