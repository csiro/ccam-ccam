module tracers_m

implicit none

private
public ngas,ntrac,ntracmax,nllp
public ilt,jlt,klt,ngasmax
public tr,traver,acloss_g !,gasmin
public trpm,npm
public tracers_init,tracers_end

! parameters should be controlled from namelist
integer, parameter :: nllp=0
integer, save :: ngas=0
integer, save :: ntrac=0
integer, save :: ntracmax=1
integer, save :: ngasmax=1
integer, save :: ilt=1
integer, save :: jlt=1
integer, save :: klt=1
real, dimension(:,:,:), allocatable, save :: tr,traver
real, dimension(:,:,:), allocatable, save :: trpm
real, dimension(:), allocatable, save :: acloss_g
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

return
end subroutine tracers_init

subroutine tracers_end

implicit none

deallocate(tr,traver)

return
end subroutine tracers_end

end module tracers_m
