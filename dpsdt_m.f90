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
    
module dpsdt_m

implicit none

private
public dpsdt,dpsdtb,dpsdtbb
public dpsdt_init,dpsdt_end

real, dimension(:), allocatable, save :: dpsdt,dpsdtb,dpsdtbb

contains

subroutine dpsdt_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(dpsdt(ifull),dpsdtb(ifull),dpsdtbb(ifull))
dpsdt=0.
dpsdtb=0.
dpsdtbb=0.

return
end subroutine dpsdt_init

subroutine dpsdt_end

implicit none

deallocate(dpsdt,dpsdtb,dpsdtbb)

return
end subroutine dpsdt_end

end module dpsdt_m