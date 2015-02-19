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
integer, intent(in) :: iqg
integer i,j,n,kx,ky
real(kind=8) den
real(kind=8), dimension(3) :: ans
n = (iqg - 1) / (il_g*il_g)
j = 1 + (iqg - n*il_g*il_g - 1)/il_g
i = iqg - (j - 1)*il_g - n*il_g*il_g
kx = 4*i-1
ky = 4*j-1
select case(n)
  case(0)
    ans(1) = 1._8
    ans(2) = xx4(kx,ky)
    ans(3) = yy4(kx,ky)
  case(1)
    ans(1) = -yy4(kx,ky)
    ans(2) = xx4(kx,ky)
    ans(3) = 1._8
  case(2)
    ans(1) = -yy4(kx,ky)
    ans(2) = 1._8
    ans(3) = -xx4(kx,ky)
  case(3)
    ans(1) = -1._8
    ans(2) = -yy4(kx,ky)
    ans(3) = -xx4(kx,ky)
  case(4)
    ans(1) = xx4(kx,ky)
    ans(2) = -yy4(kx,ky)
    ans(3) = -1._8
  case(5)
    ans(1) = xx4(kx,ky)
    ans(2) = -1._8
    ans(3) = yy4(kx,ky)
end select
den = sqrt(ans(1)*ans(1)+ans(2)*ans(2)+ans(3)*ans(3))
ans(:) = ans(:)/den
end function xyz_g

end module xyzinfo_m