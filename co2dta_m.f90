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
    
module co2dta_m

implicit none

private
public co251,co258,cdt51,cdt58,c2d51,c2d58,co2m51
public co2m58,cdtm51,cdtm58,c2dm51,c2dm58,stemp
public gtemp,b0,b1,b2,b3
public co231,co238,cdt31,cdt38,c2d31,c2d38
public co271,co278,cdt71,cdt78,c2d71,c2d78
public co211,co218
public co2dta_init,co2dta_end

real, save :: b0,b1,b2,b3
real, dimension(:), allocatable, save :: co2m51,co2m58,cdtm51,cdtm58,c2dm51,c2dm58,stemp,gtemp
real, dimension(:), allocatable, save :: co231,co238,cdt31,cdt38,c2d31,c2d38
real, dimension(:), allocatable, save :: co271,co278,cdt71,cdt78,c2d71,c2d78
real, dimension(:), allocatable, save :: co211,co218
real, dimension(:,:), allocatable, save :: co251,co258,cdt51,cdt58,c2d51,c2d58

contains

subroutine co2dta_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(co251(kl+1,kl+1),co258(kl+1,kl+1),cdt51(kl+1,kl+1))
allocate(cdt58(kl+1,kl+1),c2d51(kl+1,kl+1),c2d58(kl+1,kl+1))
allocate(co2m51(kl),co2m58(kl),cdtm51(kl),cdtm58(kl))
allocate(c2dm51(kl),c2dm58(kl),stemp(kl+1),gtemp(kl+1))
allocate(co231(kl+1),co238(kl+1),cdt31(kl+1))
allocate(cdt38(kl+1),c2d31(kl+1),c2d38(kl+1))
allocate(co271(kl+1),co278(kl+1),cdt71(kl+1))
allocate(cdt78(kl+1),c2d71(kl+1),c2d78(kl+1))
allocate(co211(kl+1),co218(kl+1))

b0=-.51926410e-4
b1=-.18113332e-3
b2=-.10680132e-5
b3=-.67303519e-7

return
end subroutine co2dta_init

subroutine co2dta_end

implicit none

deallocate(co2m51,co2m58,cdtm51,cdtm58,c2dm51,c2dm58,stemp,gtemp)
deallocate(co251,co258,cdt51,cdt58,c2d51,c2d58)
deallocate(co231,co238,cdt31,cdt38,c2d31,c2d38)
deallocate(co271,co278,cdt71,cdt78,c2d71,c2d78)
deallocate(co211,co218)

return
end subroutine co2dta_end

end module co2dta_m