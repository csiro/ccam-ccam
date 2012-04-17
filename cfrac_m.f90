module cfrac_m

implicit none

private
public cfrac,cffall
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac,cffall

contains

subroutine cfrac_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cfrac(ifull,kl),cffall(ifull+iextra,kl))
cfrac=0.
cffall=0.

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac,cffall)

return
end subroutine cfrac_end

end module cfrac_m
