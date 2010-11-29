module sigs_m

implicit none

private
public sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,dtmax,betm,ratha,rathb
public sigs_init,sigs_end

real, dimension(:), allocatable, save :: sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,betm,ratha,rathb
real, save :: dtmax

contains

subroutine  sigs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(sig(kl),sigmh(kl),dsig(kl),rata(kl),ratb(kl))
allocate(bet(kl),tbar(kl),tbardsig(kl),betm(kl),ratha(kl),rathb(kl))

return
end subroutine sigs_init

subroutine sigs_end

implicit none

deallocate(sig,sigmh,dsig,rata,ratb,bet,tbar,tbardsig,betm,ratha,rathb)

return
end subroutine sigs_end

end module sigs_m