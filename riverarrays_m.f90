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
    
module riverarrays_m

implicit none

private
public watbdy, outflowmask
public riverarrays_init, riverarrays_end

real, dimension(:), allocatable, save :: watbdy
logical, dimension(:), allocatable, save :: outflowmask

contains

subroutine riverarrays_init(ifull,iextra,nriver)

implicit none

integer, intent(in) :: ifull, iextra, nriver

if ( nriver>0 ) then
  allocate( watbdy(ifull+iextra) )
  allocate( outflowmask(ifull) )
  watbdy(1:ifull+iextra) = 0.
  outflowmask(1:ifull) = .false.
end if
  
return
end subroutine riverarrays_init

subroutine riverarrays_end

implicit none

if ( allocated( watbdy ) ) then
  deallocate( watbdy )
  deallocate( outflowmask )
end if

return
end subroutine riverarrays_end

end module riverarrays_m