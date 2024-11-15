! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! Tridiagonal solver for turbulent mixing
    
module trimmix_m

private
public trimmix

interface trimmix
  module procedure trimmix2, trimmix3
end interface trimmix

contains
    
pure subroutine trimmix2(a,c,rhs,imax,kl)

implicit none

integer, intent(in) :: imax, kl
integer k, iq
real, dimension(imax,kl), intent(in) :: a, c
real, dimension(imax,kl), intent(inout) :: rhs
real, dimension(imax,kl) :: e, g
real temp, b

! this routine solves the system
!   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!   and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

! the Thomas algorithm is used
do iq = 1,imax
  b=1./(1.-a(iq,1)-c(iq,1))
  e(iq,1)=c(iq,1)*b
  g(iq,1)=rhs(iq,1)*b
end do
do k = 2,kl-1
  do iq = 1,imax
    b=1.-a(iq,k)-c(iq,k)
    temp= 1./(b-a(iq,k)*e(iq,k-1))
    e(iq,k)=c(iq,k)*temp
    g(iq,k)=(rhs(iq,k)-a(iq,k)*g(iq,k-1))*temp
  end do
end do

! do back substitution to give answer now
do iq = 1,imax
  b=1.-a(iq,kl)-c(iq,kl)
  temp = 1./(b-a(iq,kl)*e(iq,kl-1))
  rhs(iq,kl)=(rhs(iq,kl)-a(iq,kl)*g(iq,kl-1))*temp
end do
do k = kl-1,1,-1
  do iq = 1,imax
    rhs(iq,k) = g(iq,k)-e(iq,k)*rhs(iq,k+1)
  end do
end do

return
end subroutine trimmix2

pure subroutine trimmix3(a,c,rhs,imax,kl,ndim)

implicit none

integer, intent(in) :: imax, kl, ndim
integer k, iq, n
real, dimension(imax,kl), intent(in) :: a, c
real, dimension(imax,kl,ndim), intent(inout) :: rhs
real, dimension(imax,ndim,kl) :: e, g
real temp, b

! this routine solves the system
!   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
!   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!   and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

! the Thomas algorithm is used
do n = 1,ndim
  do iq = 1,imax
    b = 1./(1.-a(iq,1)-c(iq,1))
    e(iq,n,1) = c(iq,1)*b
    g(iq,n,1) = rhs(iq,1,n)*b
  end do
end do
do k = 2,kl-1
  do n = 1,ndim
    do iq = 1,imax
      b = 1.-a(iq,k)-c(iq,k)
      temp = 1./(b-a(iq,k)*e(iq,n,k-1))
      e(iq,n,k) = c(iq,k)*temp
      g(iq,n,k) = (rhs(iq,k,n)-a(iq,k)*g(iq,n,k-1))*temp
    end do  
  end do
end do

! do back substitution to give answer
do n = 1,ndim
  do iq = 1,imax
    b = 1.-a(iq,kl)-c(iq,kl)
    temp = 1./(b-a(iq,kl)*e(iq,n,kl-1))
    rhs(iq,kl,n) = (rhs(iq,kl,n)-a(iq,kl)*g(iq,n,kl-1))*temp
  end do  
end do
do k = kl-1,1,-1
  do n = 1,ndim
    do iq = 1,imax
      rhs(iq,k,n) = g(iq,n,k)-e(iq,n,k)*rhs(iq,k+1,n)
    end do  
  end do
end do

return
end subroutine trimmix3

end module trimmix_m