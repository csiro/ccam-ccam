module tracers_m

implicit none

private
public iradon,ico2,ngas,nllp,ntrac,ntracmax
public ilt,jlt,klt,ngasmax
public tr,traver,radonem,trback_g,acloss_g,gasmin
public tracers_init,tracers_end

! Remove this line when all common blocks are gone
include 'newmpar.h'

! These parameters should be controlled from namelist
integer, parameter :: iradon=0
integer, parameter :: ico2=0
integer, parameter :: ngas=0
integer, parameter :: nllp=0
integer, parameter :: ntrac=ngas+nllp
integer, parameter :: ntracmax=max(ntrac,1) ! ntracmax >= 1
integer, parameter :: ngasmax=max(ngas,1)   ! ngasmax >= 1
integer, parameter :: npwr=min(ntrac,1)     ! i.e. 0 or 1 for no_gas/gas
integer, save :: ilt
integer, save :: jlt
integer, save :: klt
real, dimension(:,:,:), allocatable, save :: tr,traver
real, dimension(:), allocatable, save :: radonem,trback_g,acloss_g,gasmin

contains

subroutine tracers_init(il,jl,kl)

implicit none

integer, intent(in) :: il,jl,kl

ilt=il**npwr
jlt=jl**npwr
klt=kl**npwr

allocate(tr(ilt*jlt+iextra,klt,ntracmax),traver(ilt*jlt,klt,ntrac))
allocate(radonem(ilt*jlt),trback_g(ntrac),acloss_g(ntrac),gasmin(ngasmax))

gasmin=-1000.

return
end subroutine tracers_init

subroutine tracers_end

implicit none

deallocate(tr,traver)
deallocate(radonem,trback_g,acloss_g,gasmin)

return
end subroutine tracers_end

end module tracers_m