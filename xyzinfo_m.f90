module xyzinfo_m

implicit none

private
public x_g,y_g,z_g,wts_g
public x,y,z,wts
public xyzinfo_init,xyzinfo_end

real(kind=8), dimension(:), allocatable, save :: x_g,y_g,z_g
real, dimension(:), allocatable, save :: wts_g
real(kind=8), dimension(:), allocatable, save :: x,y,z
real, dimension(:), allocatable, save :: wts

contains

subroutine xyzinfo_init(ifull_g,ifull,iextra,myid,mbd)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,myid,mbd

if (myid==0.or.mbd/=0) then
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

end module xyzinfo_m