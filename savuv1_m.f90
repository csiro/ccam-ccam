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
    
module savuv1_m

implicit none

private
public savs1,savu1,savv1,savu2,savv2
public savuv1_init,savuv1_end

real, dimension(:,:), allocatable, save :: savs1,savu1,savv1,savu2,savv2

contains

subroutine savuv1_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl))
allocate(savu2(ifull,kl),savv2(ifull,kl))

return
end subroutine savuv1_init

subroutine savuv1_end

implicit none

deallocate(savs1,savu1,savv1)
deallocate(savu2,savv2)

return
end subroutine savuv1_end

end module savuv1_m