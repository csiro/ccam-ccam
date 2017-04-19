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
    
module radisw_m

implicit none

private
public ktop,kbtm,nclds,ktopsw,kbtmsw,emcld
public temp,temp2,press,press2,rh2o,qo3
public camt,cuvrf,cirrf,cirab,coszro,taudar
public rrco2,ssolar,rrvco2,rrvch4,rrvn2o
public rrvf11,rrvf12,rrvf113,rrvf22
public radisw_init,radisw_end

integer, dimension(:), allocatable, save :: nclds
integer, dimension(:,:), allocatable, save :: ktop,kbtm,ktopsw,kbtmsw
real, dimension(:), allocatable, save :: coszro,taudar
real, dimension(:,:), allocatable, save :: emcld,camt,cuvrf,cirrf,cirab
real, dimension(:,:), allocatable, save :: temp,temp2,press,press2,rh2o,qo3
real, save :: rrco2=0.,ssolar=0.,rrvco2=0.,rrvch4=0.,rrvn2o=0.
real, save :: rrvf11=0.,rrvf12=0.,rrvf113=0.,rrvf22=0.

contains

subroutine radisw_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(ktop(imax,kl+1),kbtm(imax,kl+1),nclds(imax),ktopsw(imax,kl+1),kbtmsw(imax,kl+1),emcld(imax,kl+1))
allocate(temp(imax,kl+1),temp2(imax,kl+1),press(imax,kl+1),press2(imax,kl+1),rh2o(imax,kl),qo3(imax,kl))
allocate(camt(imax,kl+1),cuvrf(imax,kl+1),cirrf(imax,kl+1),cirab(imax,kl+1),coszro(imax),taudar(imax))

return
end subroutine radisw_init

subroutine radisw_end

implicit none

deallocate(ktop,kbtm,nclds,ktopsw,kbtmsw,emcld)
deallocate(temp,temp2,press,press2,rh2o,qo3)
deallocate(camt,cuvrf,cirrf,cirab,coszro,taudar)

return
end subroutine radisw_end

end module radisw_m