! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module arrays_m

implicit none

private
public t,u,v,qg
public psl,ps,zs
public arrays_init,arrays_end

real, dimension(:,:), allocatable, save :: t,u,v,qg
real, dimension(:), allocatable, save :: psl
real, dimension(:), allocatable, save :: ps, zs

contains

subroutine arrays_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(t(ifull+iextra,kl),u(ifull+iextra,kl),v(ifull+iextra,kl),qg(ifull+iextra,kl))
allocate(psl(ifull+iextra),ps(ifull+iextra),zs(ifull+iextra))

t(:,:)=9.e9
u(:,:)=9.e9
v(:,:)=9.e9
qg(:,:)=9.e9
psl(:)=9.e9
ps(:)=9.e9
zs(:)=9.e9

return
end subroutine arrays_init

subroutine arrays_end

implicit none

deallocate(t,u,v,qg)
deallocate(psl,ps,zs)

return
end subroutine arrays_end

end module arrays_m