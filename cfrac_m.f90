module cfrac_m

implicit none

private
public cfrac
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac

contains

subroutine cfrac_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cfrac(ifull,kl))

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac)

return
end subroutine cfrac_end

end module cfrac_m