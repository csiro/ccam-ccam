module liqwpar_m

implicit none

private
public ifullw,qlg,qfg ! liquid water, ice water
!public qrg,qsg,qgg ! rain, snow, graupel
public liqwpar_init,liqwpar_end

integer, save :: ifullw
real, dimension(:,:), allocatable, save :: qlg,qfg

contains

subroutine liqwpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(qlg(ifull+iextra,kl),qfg(ifull+iextra,kl))
!allocate(qrg(ifull+iextra,kl),qsg(ifull+iextra,kl))
!allocate(qgg(ifull+iextra,kl))
ifullw=ifull
qlg=0.
qfg=0.

return
end subroutine liqwpar_init

subroutine liqwpar_end

implicit none

deallocate(qlg,qfg)
!deallocate(qrg,qsg)
!deallocate(qgg)

return
end subroutine liqwpar_end

end module liqwpar_m