module nlin_m

implicit none

private
public tn,un,vn
public nlin_init,nlin_end

real, dimension(:,:), allocatable, save :: tn,un,vn

contains

subroutine nlin_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tn(ifull,kl),un(ifull,kl),vn(ifull,kl))

return
end subroutine nlin_init

subroutine nlin_end

implicit none

deallocate(tn,un,vn)

return
end subroutine nlin_end

end module nlin_m