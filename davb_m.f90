module davb_m

implicit none

private
public psls,qgg,tt,uu,vv
public xtgdav
public davb_init,davb_end

real, dimension(:), allocatable, save :: psls
real, dimension(:,:), allocatable, save :: qgg,tt,uu,vv
real, dimension(:,:,:), allocatable, save :: xtgdav

contains

subroutine davb_init(ifull,iextra,kl,naero)

implicit none

integer, intent(in) :: ifull,iextra,kl,naero

allocate(psls(ifull),qgg(ifull,kl),tt(ifull,kl),uu(ifull,kl),vv(ifull,kl))
if (naero>0) then
  allocate(xtgdav(ifull,kl,naero))
end if

return
end subroutine davb_init

subroutine davb_end

implicit none

deallocate(psls,qgg,tt,uu,vv)
if (allocated(xtgdav)) then
  deallocate(xtgdav)
end if

return
end subroutine davb_end

end module davb_m