module cfrac_m

! cfrac is total cloud fraction (including convection)
! rfrac is the total rain fraction (not including convection)
! sfrac is the total snow fraction
! gfrac is the total grauple fraction
! see cloudmod for large scale cloud fraction with ncloud>3
    
implicit none

private
public cfrac,rfrac
!public sfrac,gfrac
public cfrac_init,cfrac_end

real, dimension(:,:), allocatable, save :: cfrac,rfrac
!real, dimension(:,:), allocatable, save :: sfrac,gfrac

contains

subroutine cfrac_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cfrac(ifull,kl),rfrac(ifull+iextra,kl))
!allocate(sfrac(ifull,kl),gfrac(ifull+iextra,kl))
cfrac=0.
rfrac=0.
!sfrac=0.
!gfrac=0.

return
end subroutine cfrac_init

subroutine cfrac_end

implicit none

deallocate(cfrac,rfrac)
!deallocate(sfrac,gfrac)

return
end subroutine cfrac_end

end module cfrac_m