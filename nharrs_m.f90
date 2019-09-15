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
    
module nharrs_m

implicit none

private
public phi,phi_nh,h_nh
public lrestart, lrestart_radiation, always_mspeca
public nharrs_init,nharrs_end

real, dimension(:,:), allocatable, save :: phi,phi_nh,h_nh
logical, save :: lrestart, lrestart_radiation, always_mspeca

contains

subroutine nharrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(phi(ifull,kl),phi_nh(ifull,kl),h_nh(ifull+iextra,kl))
phi=-999.
phi_nh=0.
lrestart=.false.
lrestart_radiation=.false.
always_mspeca = .false.

return
end subroutine nharrs_init

subroutine nharrs_end

implicit none

deallocate(phi,phi_nh,h_nh)

return
end subroutine nharrs_end

end module nharrs_m