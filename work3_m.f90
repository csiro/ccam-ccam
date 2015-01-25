module work3_m

implicit none

private
public ga,condxpr,condspr,fes
public fwtop
public work3_init,work3_end

real, dimension(:), allocatable, save :: ga,condxpr,condspr,fes
real, dimension(:), allocatable, save :: fwtop

contains

subroutine work3_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib

allocate(ga(ifull))
if (nsib==3.or.nsib==5) then
  allocate(condxpr(ifull),condspr(ifull),fes(ifull))
  allocate(fwtop(ifull))
end if

return
end subroutine work3_init

subroutine work3_end

implicit none

deallocate(ga)
if (allocated(condxpr)) then
  deallocate(condxpr,condspr,fes)
  deallocate(fwtop)
end if

return
end subroutine work3_end

end module work3_m
