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

module savuvt_m

implicit none

private
public savt, savpsl
public savs, savu, savv
public savuvt_init, savuvt_end

real, dimension(:), allocatable, save :: savpsl
real, dimension(:,:), allocatable, save :: savt, savs, savu, savv

contains

subroutine savuvt_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate( savt(ifull,kl), savpsl(ifull) )
allocate( savs(ifull,2:kl), savu(ifull,kl), savv(ifull,kl) )

return
end subroutine savuvt_init

subroutine savuvt_end

implicit none

deallocate( savt, savpsl )
deallocate( savs, savu, savv )

return
end subroutine savuvt_end

end module savuvt_m