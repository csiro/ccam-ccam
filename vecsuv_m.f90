module vecsuv_m

implicit none

private
public ax_g,bx_g,ay_g,by_g,az_g,bz_g
public ax,bx,ay,by,az,bz
public vecsuv_init,vecsuv_end

real, dimension(:), allocatable, save :: ax_g,bx_g,ay_g,by_g,az_g,bz_g
real, dimension(:), allocatable, save :: ax,bx,ay,by,az,bz

contains

subroutine vecsuv_init(ifull_g,ifull,iextra,myid,mbd,nud_uv)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,myid,mbd,nud_uv

if (myid==0.or.(mbd/=0.and.nud_uv==9)) then
  allocate(ax_g(ifull_g),bx_g(ifull_g))
  allocate(ay_g(ifull_g),by_g(ifull_g))
  allocate(az_g(ifull_g),bz_g(ifull_g))
end if
allocate(ax(ifull+iextra),bx(ifull+iextra))
allocate(ay(ifull+iextra),by(ifull+iextra))
allocate(az(ifull+iextra),bz(ifull+iextra))

return
end subroutine vecsuv_init

subroutine vecsuv_end

implicit none

if (allocated(ax_g)) deallocate(ax_g)
if (allocated(bx_g)) deallocate(bx_g)
if (allocated(ay_g)) deallocate(ay_g)
if (allocated(by_g)) deallocate(by_g)
if (allocated(az_g)) deallocate(az_g)
if (allocated(bz_g)) deallocate(bz_g)
deallocate(ax,bx,ay,by,az,bz)

return
end subroutine vecsuv_end

end module vecsuv_m