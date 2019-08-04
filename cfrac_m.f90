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

! cfrac is total cloud fraction (including convection) for the radiation code
! stratcloud is the stratiform cloud fraction for the condensate
! rfrac is the total rain fraction (not including convection)
! sfrac is the total snow fraction
! gfrac is the total graupel fraction
! nettend is the temperature tendency for prognostic cloud (ncloud>3)
    
implicit none

private
public cfrac, rfrac
public sfrac, gfrac
public stratcloud
public nettend
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac,rfrac
real, dimension(:,:), allocatable, save :: sfrac,gfrac
real, dimension(:,:), allocatable, save :: stratcloud
real, dimension(:,:), allocatable, save :: nettend

contains

subroutine cfrac_init(ifull,iextra,kl,ncloud)

implicit none

integer, intent(in) :: ifull, iextra, kl, ncloud

allocate(cfrac(ifull,kl),rfrac(ifull,kl))
allocate(sfrac(ifull,kl),gfrac(ifull,kl))
allocate(stratcloud(ifull+iextra,kl))
cfrac = 0.
rfrac = 0.
sfrac = 0.
gfrac = 0.
stratcloud = 0.

if ( ncloud>=4 ) then
  allocate(nettend(ifull,kl))
  nettend = 0.
end if

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac,rfrac)
deallocate(sfrac,gfrac)
deallocate(stratcloud)

if ( allocated(nettend) ) then
  deallocate(nettend)
end if

return
end subroutine cfrac_end

end module cfrac_m