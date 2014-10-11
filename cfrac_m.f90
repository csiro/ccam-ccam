module cfrac_m

! cfrac is total cloud fraction (including convection)
! cffall is the total rain fraction (not including convection)
! see cloudmod for large scale cloud fraction with ncloud>2
    
implicit none

private
public cfrac,rfrac
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac,rfrac

contains

subroutine cfrac_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cfrac(ifull,kl),rfrac(ifull+iextra,kl))
cfrac=0.
rfrac=0.

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac,rfrac)

return
end subroutine cfrac_end

end module cfrac_m