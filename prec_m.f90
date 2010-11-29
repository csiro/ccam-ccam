module prec_m

implicit none

private
public evap,precip,precc,rnd_3hr,cape
public prec_init,prec_end

real, dimension(:), allocatable, save :: evap,precip,precc,cape
real, dimension(:,:), allocatable, save :: rnd_3hr

contains

subroutine prec_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(evap(ifull),precip(ifull),precc(ifull),rnd_3hr(ifull,8),cape(ifull))

return
end subroutine prec_init

subroutine prec_end

implicit none

deallocate(evap,precip,precc,rnd_3hr,cape)

return
end subroutine prec_end

end module prec_m