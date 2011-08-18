module tracers_m

implicit none

private
public ngas,ntrac,ntracmax,nllp
public ilt,jlt,klt,ngasmax
public tr,traver,trback_g,acloss_g !,gasmin
public trpm,npm
public tracers_init,tracers_end

! These parameters should be controlled from namelist
integer, save :: ngas=0
integer, parameter :: nllp=0
integer, save :: ntrac
integer, save :: ntracmax
integer, save :: ngasmax
integer, save :: ilt
integer, save :: jlt
integer, save :: klt
real, dimension(:,:,:), allocatable, save :: tr,traver
real, dimension(:,:,:), allocatable, save :: trpm
real, dimension(:), allocatable, save :: trback_g,acloss_g !,gasmin
integer, dimension(:), allocatable, save :: npm

contains

subroutine tracers_init(il,jl,kl,iextra)

implicit none

integer, intent(in) :: il,jl,kl,iextra

ntrac=ngas+nllp
ntracmax=max(ntrac,1) ! ntracmax >= 1
ngasmax=max(ngas,1)   ! ngasmax >= 1

! old trick for common blocks
! now just set ilt to il or 1, etc
!ilt=il**npwr
!jlt=jl**npwr
!klt=kl**npwr

ilt=il
jlt=jl
klt=kl

allocate(tr(ilt*jlt+iextra,klt,ntracmax),traver(ilt*jlt,klt,ntrac))
allocate(trback_g(ntrac))
!allocate(gasmin(ngasmax))

trback_g=0.
!gasmin=-1000.

return
end subroutine tracers_init

subroutine tracers_end

implicit none

deallocate(tr,traver)
deallocate(trback_g)
!deallocate(gasmin)

return
end subroutine tracers_end

end module tracers_m
