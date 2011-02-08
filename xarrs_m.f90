module xarrs_m

implicit none

private
public ux,vx,tx,pslx
public xarrs_init,xarrs_end

real, dimension(:,:), allocatable, save :: ux,vx,tx,pslx

contains

subroutine xarrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(ux(ifull,kl),vx(ifull,kl),tx(ifull+iextra,kl),pslx(ifull+iextra,kl))

return
end subroutine xarrs_init

subroutine xarrs_end

implicit none

deallocate(ux,vx,tx,pslx)

return
end subroutine xarrs_end

end module xarrs_m