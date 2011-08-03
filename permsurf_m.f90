module permsurf_m

implicit none

private
public ipsice,ipsea,ipland,iperm
public permsurf_init,permsurf_end

integer, save :: ipsice,ipsea,ipland
integer, dimension(:), allocatable, save :: iperm

contains

subroutine permsurf_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(iperm(ifull))

return
end subroutine permsurf_init

subroutine permsurf_end

implicit none

deallocate(iperm)

return
end subroutine permsurf_end

end module permsurf_m