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
    
module work3lwr_m

implicit none

private
public co2r,dift
public co2r1,dco2d1
public d2cd21,d2cd22
public co2r2,dco2d2
public co2mr,co2md,co2m2d
public tdav,tstdav
public vv,vsum3,vsum1,vsum2
public a1,a2
public work3lwr_init,work3lwr_end

real, dimension(:,:,:), allocatable, save :: co2r,dift
real, dimension(:,:), allocatable, save :: co2r1,dco2d1
real, dimension(:,:), allocatable, save :: d2cd21,d2cd22
real, dimension(:,:), allocatable, save :: co2r2,dco2d2
real, dimension(:,:), allocatable, save :: co2mr,co2md,co2m2d
real, dimension(:,:), allocatable, save :: tdav,tstdav
real, dimension(:,:), allocatable, save :: vv,vsum3
real, dimension(:), allocatable, save :: vsum1,vsum2,a1,a2

contains

subroutine work3lwr_init(kl,imax)

implicit none

integer, intent(in) :: kl,imax

allocate(co2r(imax,kl+1,kl+1),dift(imax,kl+1,kl+1))
allocate(co2r1(imax,kl+1),dco2d1(imax,kl+1))
allocate(d2cd21(imax,kl+1),d2cd22(imax,kl+1))
allocate(co2r2(imax,kl+1),dco2d2(imax,kl+1))
allocate(co2mr(imax,kl),co2md(imax,kl),co2m2d(imax,kl))
allocate(tdav(imax,kl+1),tstdav(imax,kl+1))
allocate(vv(imax,kl),vsum3(imax,kl+1),vsum1(imax),vsum2(imax))
allocate(a1(imax),a2(imax))

return
end subroutine work3lwr_init

subroutine work3lwr_end

implicit none

deallocate(co2r,dift)
deallocate(co2r1,dco2d1)
deallocate(d2cd21,d2cd22)
deallocate(co2r2,dco2d2)
deallocate(co2mr,co2md,co2m2d)
deallocate(tdav,tstdav)
deallocate(vv,vsum3,vsum1,vsum2)
deallocate(a1,a2)

return
end subroutine work3lwr_end

end module work3lwr_m