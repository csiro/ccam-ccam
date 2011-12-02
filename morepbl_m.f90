module morepbl_m

implicit none

private
public condx,fg,eg,epot,condc,rnet,pblh,epan,tpan
public morepbl_init,morepbl_end

real, dimension(:), allocatable, save :: condx,fg,eg,epot,condc,rnet,pblh,epan,tpan

contains

subroutine morepbl_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(condx(ifull),fg(ifull),eg(ifull),epot(ifull))
allocate(condc(ifull),rnet(ifull),pblh(ifull),epan(ifull))
allocate(tpan(ifull))
fg=0.
eg=0.
pblh=1000.

return
end subroutine morepbl_init

subroutine morepbl_end

implicit none

deallocate(condx,fg,eg,epot,condc,rnet,pblh,epan,tpan)

return
end subroutine morepbl_end

end module morepbl_m
