module liqwpar_m

implicit none

private
public ifullw,qlg,qfg ! liquid water, ice water
public qrg !,qsg,qgrau ! rain, snow, graupel
public liqwpar_init,liqwpar_end

integer, save :: ifullw
real, dimension(:,:), allocatable, save :: qlg,qfg
real, dimension(:,:), allocatable, save :: qrg !,qsg
!real, dimension(:,:), allocatable, save :: qgrau

contains

subroutine liqwpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(qlg(ifull+iextra,kl),qfg(ifull+iextra,kl))
allocate(qrg(ifull,kl)) !,qsg(ifull,kl))
!allocate(qgrau(ifull,kl))
ifullw=ifull
qlg=0.
qfg=0.
qrg=0.
!qsg=0.
!qgrau=0.

return
end subroutine liqwpar_init

subroutine liqwpar_end

implicit none

deallocate(qlg,qfg)
deallocate(qrg) !,qsg)
!deallocate(qgrau)

return
end subroutine liqwpar_end

end module liqwpar_m