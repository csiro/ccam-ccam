module vegpar_m

implicit none

private
public cansto,vlai
public sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
public dsumpn,dsumrp,dsumrs,dsumrd
public rlai,c4frac
public vegpar_init,vegpar_end

real, dimension(:), allocatable, save :: cansto,vlai
real, dimension(:), allocatable, save :: sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
real, dimension(:), allocatable, save :: dsumpn,dsumrp,dsumrs,dsumrd
real, dimension(:), allocatable, save :: rlai,c4frac

contains

subroutine vegpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cansto(ifull),vlai(ifull))
allocate(sumpn(ifull),sumrp(ifull),sumrpw(ifull),sumrpr(ifull),sumrs(ifull),sumrd(ifull))
allocate(dsumpn(ifull),dsumrp(ifull),dsumrs(ifull),dsumrd(ifull))
allocate(rlai(ifull),c4frac(ifull))

return
end subroutine vegpar_init

subroutine vegpar_end

implicit none

deallocate(cansto,vlai)
deallocate(sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd)
deallocate(dsumpn,dsumrp,dsumrs,dsumrd)
deallocate(rlai,c4frac)

return
end subroutine vegpar_end

end module vegpar_m