module pbl_m

implicit none

private
public cduv,cdtq,tss,slwa
public pbl_init,pbl_end

real, dimension(:), allocatable, save :: cduv,cdtq,tss,slwa

contains

subroutine pbl_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cduv(ifull),cdtq(ifull),tss(ifull),slwa(ifull))

return
end subroutine pbl_init

subroutine pbl_end

implicit none

deallocate(cduv,cdtq,tss,slwa)

return
end subroutine pbl_end

end module pbl_m