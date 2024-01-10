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
    
module parmhdff_m

implicit none

private
public nhor,nhorps,khor,khdif,nhorjlm
public hdiff
public parmhdff_init,parmhdff_end

integer, save :: nhor,nhorps,khor,khdif,nhorjlm
real, dimension(:), allocatable, save :: hdiff

contains

subroutine parmhdff_init(kl)

implicit none

integer, intent(in) :: kl

allocate(hdiff(kl))
hdiff=0.1*khdif

return
end subroutine parmhdff_init

subroutine parmhdff_end

implicit none

deallocate(hdiff)

return
end subroutine parmhdff_end

end module parmhdff_m