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
    
module kdacom_m

implicit none

private
public qh2o,p,delp2,delp,t,var1,var2
public var3,var4,cntval
public kdacom_init,kdacom_end

real, dimension(:,:), allocatable, save :: qh2o,p,delp2,delp,t,var1,var2
real, dimension(:,:), allocatable, save :: var3,var4,cntval

contains

subroutine kdacom_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(qh2o(imax,kl+1),p(imax,kl+1),delp2(imax,kl),delp(imax,kl),t(imax,kl+1),var1(imax,kl),var2(imax,kl))
allocate(var3(imax,kl),var4(imax,kl),cntval(imax,kl+1))

return
end subroutine kdacom_init

subroutine kdacom_end

implicit none

deallocate(qh2o,p,delp2,delp,t,var1,var2)
deallocate(var3,var4,cntval)

return
end subroutine kdacom_end

end module kdacom_m