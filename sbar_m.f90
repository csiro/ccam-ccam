module sbar_m

implicit none

private
public sbar
public sbar_init,sbar_end

real, dimension(:,:), allocatable, save :: sbar

contains

subroutine sbar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(sbar(ifull,2:kl))

return
end subroutine sbar_init

subroutine sbar_end

implicit none

deallocate(sbar)

return
end subroutine sbar_end

end module sbar_m