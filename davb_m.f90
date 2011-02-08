module davb_m

implicit none

private
public psls,qgg,tt,uu,vv
public davb_init,davb_end

real, dimension(:), allocatable, save :: psls
real, dimension(:,:), allocatable, save :: qgg,tt,uu,vv

contains

subroutine davb_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(psls(ifull),qgg(ifull,kl),tt(ifull,kl),uu(ifull,kl),vv(ifull,kl))

return
end subroutine davb_init

subroutine davb_end

implicit none

deallocate(psls,qgg,tt,uu,vv)

return
end subroutine davb_end

end module davb_m