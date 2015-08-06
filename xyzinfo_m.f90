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
    
module xyzinfo_m

implicit none

private
public xyz_g
public x_g,y_g,z_g,wts_g
public x,y,z,wts
public xyzinfo_init,xyzinfo_end

real(kind=8), dimension(:), allocatable, save :: x_g,y_g,z_g
real, dimension(:), allocatable, save :: wts_g
real(kind=8), dimension(:), allocatable, save :: x,y,z
real, dimension(:), allocatable, save :: wts

contains

subroutine xyzinfo_init(ifull_g,ifull,iextra,myid,mbd,nud_uv)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,myid,mbd,nud_uv

if (myid==0.or.(mbd/=0.and.nud_uv==9)) then
  allocate(x_g(ifull_g),y_g(ifull_g),z_g(ifull_g))
end if
if (myid==0) then
  allocate(wts_g(ifull_g))  
end if
allocate(x(ifull),y(ifull),z(ifull))
allocate(wts(ifull))

return
end subroutine xyzinfo_init

subroutine xyzinfo_end

implicit none

if (allocated(x_g)) deallocate(x_g,y_g,z_g)
if (allocated(wts_g)) deallocate(wts_g)
deallocate(x,y,z)
deallocate(wts)

return
end subroutine xyzinfo_end

function xyz_g(iqg) result(ans)
use bigxy4_m
implicit none
include 'newmpar.h'
include 'parmgeom.h'
integer, intent(in) :: iqg
integer i,j,n,kx,ky
real(kind=8) den, alf
real(kind=8), dimension(3) :: ans
real(kind=8), dimension(3) :: ain=0.
n = (iqg - 1) / (il_g*il_g)
j = 1 + (iqg - n*il_g*il_g - 1)/il_g
i = iqg - (j - 1)*il_g - n*il_g*il_g
kx = 4*i-1
ky = 4*j-1
select case(n)
  case(0)
    ain(1) = 1._8
    ain(2) = xx4(kx,ky)
    ain(3) = yy4(kx,ky)
  case(1)
    ain(1) = -yy4(kx,ky)
    ain(2) = xx4(kx,ky)
    ain(3) = 1._8
  case(2)
    ain(1) = -yy4(kx,ky)
    ain(2) = 1._8
    ain(3) = -xx4(kx,ky)
  case(3)
    ain(1) = -1._8
    ain(2) = -yy4(kx,ky)
    ain(3) = -xx4(kx,ky)
  case(4)
    ain(1) = xx4(kx,ky)
    ain(2) = -yy4(kx,ky)
    ain(3) = -1._8
  case(5)
    ain(1) = xx4(kx,ky)
    ain(2) = -1._8
    ain(3) = yy4(kx,ky)
end select
den = sqrt(ain(1)*ain(1)+ain(2)*ain(2)+ain(3)*ain(3))
ain(:) = ain(:)/den

alf=(1._8-schmidt**2)/(1._8+schmidt**2)
ans(1)=ain(1)*schmidt*(1._8+alf)/(1._8+alf*ain(3))
ans(2)=ain(2)*schmidt*(1._8+alf)/(1._8+alf*ain(3))
ans(3)=(alf+ain(3))/(1._8+alf*ain(3))

end function xyz_g

end module xyzinfo_m