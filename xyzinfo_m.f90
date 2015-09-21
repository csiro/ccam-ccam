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

interface xyz_g
  module procedure xyz_g_s, xyz_g_v
end interface xyz_g

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

function xyz_g_s(iqg) result(ans)
use bigxy4_m
implicit none
include 'newmpar.h'
include 'parmgeom.h'
integer, intent(in) :: iqg
integer i, j, n, kx, ky, il2
real(kind=8) den, alf
real(kind=8), dimension(3) :: ans
real(kind=8), dimension(3) :: ain=0.

il2 = il_g*il_g
n = (iqg - 1) / il2
j = 1 + (iqg - n*il2 - 1)/il_g
i = iqg - (j - 1)*il_g - n*il2
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

alf = (1._8-schmidt*schmidt)/(1._8+schmidt*schmidt)
ans(1) = ain(1)*schmidt*(1._8+alf)/(1._8+alf*ain(3))
ans(2) = ain(2)*schmidt*(1._8+alf)/(1._8+alf*ain(3))
ans(3) = (alf+ain(3))/(1._8+alf*ain(3))

end function xyz_g_s

function xyz_g_v(ibeg,iend,jin,a,b,c) result(ans)
use bigxy4_m
implicit none
include 'newmpar.h'
include 'parmgeom.h'
integer, intent(in) :: ibeg, iend, a, b, c, jin
integer i, j, n, kx, ky, il2, iqg, m, fn
real(kind=8) alf
real(kind=8), dimension(iend-ibeg+1) :: den
real(kind=8), dimension(iend-ibeg+1,3) :: ans
real(kind=8), dimension(iend-ibeg+1,3) :: ain

il2 = il_g*il_g
fn = iend - ibeg + 1
ain(:,:) = 0._8

do m = 1,fn
  iqg = (m-1+ibeg)*a + b*jin + c
  n = (iqg - 1) / il2
  j = 1 + (iqg - n*il2 - 1)/il_g
  i = iqg - (j - 1)*il_g - n*il2
  kx = 4*i-1
  ky = 4*j-1
  select case(n)
    case(0)
      ain(m,1) = 1._8
      ain(m,2) = xx4(kx,ky)
      ain(m,3) = yy4(kx,ky)
    case(1)
      ain(m,1) = -yy4(kx,ky)
      ain(m,2) = xx4(kx,ky)
      ain(m,3) = 1._8
    case(2)
      ain(m,1) = -yy4(kx,ky)
      ain(m,2) = 1._8
      ain(m,3) = -xx4(kx,ky)
    case(3)
      ain(m,1) = -1._8
      ain(m,2) = -yy4(kx,ky)
      ain(m,3) = -xx4(kx,ky)
    case(4)
      ain(m,1) = xx4(kx,ky)
      ain(m,2) = -yy4(kx,ky)
      ain(m,3) = -1._8
    case(5)
      ain(m,1) = xx4(kx,ky)
      ain(m,2) = -1._8
      ain(m,3) = yy4(kx,ky)
  end select
end do

den(1:fn) = sqrt(ain(1:fn,1)*ain(1:fn,1)+ain(1:fn,2)*ain(1:fn,2)+ain(1:fn,3)*ain(1:fn,3))
ain(1:fn,1) = ain(1:fn,1)/den(1:fn)
ain(1:fn,2) = ain(1:fn,2)/den(1:fn)
ain(1:fn,3) = ain(1:fn,3)/den(1:fn)

alf = (1._8-schmidt*schmidt)/(1._8+schmidt*schmidt)
ans(1:fn,1) = ain(1:fn,1)*schmidt*(1._8+alf)/(1._8+alf*ain(1:fn,3))
ans(1:fn,2) = ain(1:fn,2)*schmidt*(1._8+alf)/(1._8+alf*ain(1:fn,3))
ans(1:fn,3) = (alf+ain(1:fn,3))/(1._8+alf*ain(1:fn,3))

end function xyz_g_v

end module xyzinfo_m
