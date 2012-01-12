module dpsdt_m

implicit none

private
public dpsdt,dpsdtb,dpsdtbb
public dpsdt_init,dpsdt_end

real, dimension(:), allocatable, save :: dpsdt,dpsdtb,dpsdtbb

contains

subroutine dpsdt_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(dpsdt(ifull),dpsdtb(ifull),dpsdtbb(ifull))
dpsdt=0.
dpsdtb=0.
dpsdtbb=0.

return
end subroutine dpsdt_init

subroutine dpsdt_end

implicit none

deallocate(dpsdt,dpsdtb,dpsdtbb)

return
end subroutine dpsdt_end

end module dpsdt_m