module map_m

implicit none

private
public em_g, emu_g, emv_g, f_g, fu_g, fv_g
public dmdx_g, dmdy_g
public em, emu, emv, f, fu, fv
public dmdx, dmdy
public map_init,map_end

real, dimension(:), allocatable, save :: em_g, emu_g, emv_g, f_g, fu_g, fv_g
real, dimension(:), allocatable, save :: dmdx_g, dmdy_g
real, dimension(:), allocatable, save :: em, emu, emv, f, fu, fv
real, dimension(:), allocatable, save :: dmdx, dmdy

contains

subroutine map_init(ifull_g,ifull,iextra)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra

allocate(em_g(ifull_g),emu_g(ifull_g),emv_g(ifull_g))
allocate(f_g(ifull_g),fu_g(ifull_g),fv_g(ifull_g))
allocate(dmdx_g(ifull_g),dmdy_g(ifull_g))
allocate(em(ifull+iextra),emu(ifull+iextra),emv(ifull+iextra))
allocate(f(ifull+iextra),fu(ifull+iextra),fv(ifull+iextra))
allocate(dmdx(ifull),dmdy(ifull))

return
end subroutine map_init

subroutine map_end

implicit none

deallocate(em_g)
if (allocated(emu_g)) deallocate(emu_g)
if (allocated(emv_g)) deallocate(emv_g)
if (allocated(f_g)) deallocate(f_g)
if (allocated(fu_g)) deallocate(fu_g)
if (allocated(fv_g)) deallocate(fv_g)
deallocate(em,emu,emv,f,fu,fv)
deallocate(dmdx,dmdy)

return
end subroutine map_end

end module map_m