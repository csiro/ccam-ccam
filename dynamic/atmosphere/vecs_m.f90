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
    
module vecs_m

implicit none

private
public emat,einv,bam,bami,qvec,tmat
public vecs_init,vecs_end

real, dimension(:), allocatable, save :: bam,bami,qvec
real, dimension(:,:), allocatable, save :: emat,einv,tmat

contains

subroutine vecs_init(kl)

implicit none

integer, intent(in) :: kl

allocate(bam(kl),bami(kl),qvec(kl))
allocate(emat(kl,kl),einv(kl,kl),tmat(kl,kl))

bam(:)=0.
bami(:)=0.
qvec(:)=0.
emat(:,:)=0.
einv(:,:)=0.
tmat(:,:)=0.

return
end subroutine vecs_init

subroutine vecs_end

implicit none

deallocate(bam,bami,qvec)
deallocate(emat,einv,tmat)

return
end subroutine vecs_end

end module vecs_m
