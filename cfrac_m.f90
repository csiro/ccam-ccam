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
    
module cfrac_m

! cfrac is total cloud fraction (including convection)
! rfrac is the total rain fraction (not including convection)
! sfrac is the total snow fraction
! gfrac is the total graupel fraction
! see cloudmod for large scale cloud fraction with ncloud>3
    
implicit none

private
public cfrac,rfrac
public sfrac,gfrac
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac,rfrac
real, dimension(:,:), allocatable, save :: sfrac,gfrac

contains

subroutine cfrac_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cfrac(ifull,kl),rfrac(ifull+iextra,kl))
allocate(sfrac(ifull+iextra,kl),gfrac(ifull+iextra,kl))
cfrac=0.
rfrac=0.
sfrac=0.
gfrac=0.

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac,rfrac)
deallocate(sfrac,gfrac)

return
end subroutine cfrac_end

end module cfrac_m