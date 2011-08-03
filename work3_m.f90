module work3_m

implicit none

private
public egg,evapxf,ewww,fgf
public fgg,ggflux,rdg,rgg,residf
public ga,condxpr,fev,fes
public ism,fwtop,af,extin
public work3_init,work3_end

real, dimension(:), allocatable, save :: egg,evapxf,ewww,fgf
real, dimension(:), allocatable, save :: fgg,ggflux,rdg,rgg,residf
real, dimension(:), allocatable, save :: ga,condxpr,fev,fes
real, dimension(:), allocatable, save :: fwtop,af,extin
integer, dimension(:), allocatable, save :: ism

contains

subroutine work3_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(egg(ifull),evapxf(ifull),ewww(ifull),fgf(ifull))
allocate(fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull))
allocate(ga(ifull),condxpr(ifull),fev(ifull),fes(ifull))
allocate(ism(ifull),fwtop(ifull),af(ifull),extin(ifull))

return
end subroutine work3_init

subroutine work3_end

implicit none

deallocate(egg,evapxf,ewww,fgf)
deallocate(fgg,ggflux,rdg,rgg,residf)
deallocate(ga,condxpr,fev,fes)
deallocate(ism,fwtop,af,extin)

return
end subroutine work3_end

end module work3_m