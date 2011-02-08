module rdflux_m

implicit none

private
public flx1e1,gxcts,fctsg,flx1e1clr,gxctsclr
public rdflux_init,rdflux_end

real, dimension(:), allocatable, save :: flx1e1,gxcts,flx1e1clr,gxctsclr
real, dimension(:,:), allocatable, save :: fctsg

contains

subroutine rdflux_init(ifull,iextra,kl,imax,nbly)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax,nbly

allocate(flx1e1(imax),gxcts(imax),fctsg(imax,nbly),flx1e1clr(imax),gxctsclr(imax))

return
end subroutine rdflux_init

subroutine rdflux_end

implicit none

deallocate(flx1e1,gxcts,fctsg,flx1e1clr,gxctsclr)

return
end subroutine rdflux_end

end module rdflux_m