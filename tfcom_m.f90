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
    
module tfcom_m

implicit none

private
public to3,co21,emiss,emiss2,avephi,cts,ctso3
public excts,exctsn,e1flx,co2nbl,co2sp1,co2sp2
public co2sp,to3spc,totvo2
public tfcom_init,tfcom_end

real, dimension(:,:,:), allocatable, save :: to3,co21,emiss,emiss2,avephi,exctsn
real, dimension(:,:), allocatable, save :: cts,ctso3,excts
real, dimension(:,:), allocatable, save :: e1flx,co2nbl,co2sp1,co2sp2
real, dimension(:,:), allocatable, save :: co2sp,to3spc,totvo2

contains

subroutine tfcom_init(kl,imax,nbly)

implicit none

integer, intent(in) :: kl,imax,nbly

allocate(to3(imax,kl+1,kl+1),co21(imax,kl+1,kl+1),emiss(imax,kl+1,kl+1),emiss2(imax,kl+1,kl+1),avephi(imax,kl+1,kl+1))
allocate(cts(imax,kl),ctso3(imax,kl),excts(imax,kl),exctsn(imax,kl,nbly))
allocate(e1flx(imax,kl+1),co2nbl(imax,kl),co2sp1(imax,kl+1),co2sp2(imax,kl+1))
allocate(co2sp(imax,kl+1),to3spc(imax,kl),totvo2(imax,kl+1))

return
end subroutine tfcom_init

subroutine tfcom_end

implicit none

deallocate(to3,co21,emiss,emiss2,avephi,cts,ctso3)
deallocate(excts,exctsn,e1flx,co2nbl,co2sp1,co2sp2)
deallocate(co2sp,to3spc,totvo2)

return
end subroutine tfcom_end

end module tfcom_m