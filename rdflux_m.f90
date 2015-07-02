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
    
module rdflux_m

implicit none

private
public flx1e1,gxcts,fctsg,flx1e1clr,gxctsclr
public rdflux_init,rdflux_end

real, dimension(:), allocatable, save :: flx1e1,gxcts,flx1e1clr,gxctsclr
real, dimension(:,:), allocatable, save :: fctsg

contains

subroutine rdflux_init(ifull,iextra,kl,imax,nbly)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax,nbly

allocate(flx1e1(imax),gxcts(imax),fctsg(imax,nbly),flx1e1clr(imax),gxctsclr(imax))

return
end subroutine rdflux_init

subroutine rdflux_end

implicit none

deallocate(flx1e1,gxcts,fctsg,flx1e1clr,gxctsclr)

return
end subroutine rdflux_end

end module rdflux_m