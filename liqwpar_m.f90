module liqwpar_m

implicit none

private
public ifullw,qlg,qfg
public liqwpar_init,liqwpar_end

integer, save :: ifullw
real, dimension(:,:), allocatable, save :: qlg,qfg

contains

subroutine liqwpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(qlg(ifull+iextra,kl),qfg(ifull+iextra,kl))
ifullw=ifull

return
end subroutine liqwpar_init

subroutine liqwpar_end

implicit none

deallocate(qlg,qfg)

return
end subroutine liqwpar_end

end module liqwpar_m