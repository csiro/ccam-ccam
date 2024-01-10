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
    
module lwout_m

implicit none

private
public heatra,exctsclr,ctso3clr
public grnflx,topflx,grnflxclr
public lwout_init,lwout_end

real, dimension(:,:), allocatable, save :: heatra,exctsclr,ctso3clr
real, dimension(:), allocatable, save :: grnflx,topflx,grnflxclr

contains

subroutine lwout_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(heatra(imax,kl),exctsclr(imax,kl),ctso3clr(imax,kl))
allocate(grnflx(imax),topflx(imax),grnflxclr(imax))

return
end subroutine lwout_init

subroutine lwout_end

implicit none

deallocate(heatra,exctsclr,ctso3clr)
deallocate(grnflx,topflx,grnflxclr)

return
end subroutine lwout_end

end module lwout_m