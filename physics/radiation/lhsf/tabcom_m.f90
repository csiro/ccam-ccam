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

module tabcom_m

implicit none

private
public em1,em1wde,table1,table2,table3,em3
public source,dsrce,ind,indx2,kmaxv,kmaxvm
public tabcom_init,tabcom_end

real, dimension(:,:), allocatable, save :: em1,em1wde,table1,table2,table3,em3
real, dimension(:,:), allocatable, save :: source,dsrce
integer, dimension(:), allocatable, save :: ind,indx2,kmaxv
integer, save :: kmaxvm

contains

subroutine tabcom_init(kl,imax,nbly)

implicit none

integer :: lp1v

integer, intent(in) :: kl,imax,nbly

lp1v=(kl+1)*(1+2*kl/2)
allocate(em1(28,180),em1wde(28,180),table1(28,180),table2(28,180),table3(28,180),em3(28,180))
allocate(source(28,nbly),dsrce(28,nbly),ind(imax),indx2(lp1v),kmaxv(kl+1))

return
end subroutine tabcom_init

subroutine tabcom_end

implicit none

deallocate(em1,em1wde,table1,table2,table3,em3)
deallocate(source,dsrce,ind,indx2,kmaxv)

return
end subroutine tabcom_end

end module tabcom_m