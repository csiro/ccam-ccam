module uvbar_m

implicit none

private
public ubar,vbar
public uvbar_init,uvbar_end

real, dimension(:,:), allocatable, save :: ubar,vbar

contains

subroutine uvbar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(ubar(ifull,kl),vbar(ifull,kl))

return
end subroutine uvbar_init

subroutine uvbar_end

implicit none

deallocate(ubar,vbar)

return
end subroutine uvbar_end

end module uvbar_m