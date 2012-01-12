module vecs_m

implicit none

private
public emat,einv,bam,bami,qvec,tmat
public vecs_init,vecs_end

real, dimension(:), allocatable, save :: bam,bami,qvec
real, dimension(:,:), allocatable, save :: emat,einv,tmat

contains

subroutine vecs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(bam(kl),bami(kl),qvec(kl))
allocate(emat(kl,kl),einv(kl,kl),tmat(kl,kl))

return
end subroutine vecs_init

subroutine vecs_end

implicit none

deallocate(bam,bami,qvec)
deallocate(emat,einv,tmat)

return
end subroutine vecs_end

end module vecs_m