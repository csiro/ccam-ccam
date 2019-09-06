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
    
module sigs_m

implicit none

private
public sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,dtmax,betm,ratha,rathb
public sigs_init,sigs_end

real, dimension(:), allocatable, save :: sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,betm,ratha,rathb
real, save :: dtmax
!$acc declare create(sig,dsig,sigmh,bet,betm,ratha,rathb)

contains

subroutine  sigs_init(kl)

implicit none

integer, intent(in) :: kl

allocate(sig(kl),sigmh(kl),dsig(kl),rata(kl),ratb(kl))
allocate(bet(kl),tbar(kl),tbardsig(kl),betm(kl),ratha(kl),rathb(kl))

return
end subroutine sigs_init

subroutine sigs_end

implicit none

deallocate(sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,betm,ratha,rathb)

return
end subroutine sigs_end

end module sigs_m