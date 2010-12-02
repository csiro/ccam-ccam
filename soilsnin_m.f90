module soilsnin_m

implicit none

private
public condxpr,fev,fes,ga,dgdtg
public soilsnin_init,soilsnin_end

real, dimension(:), allocatable, save :: condxpr,fev,fes,ga,dgdtg

contains

subroutine soilsnin_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(condxpr(ifull),fev(ifull),fes(ifull),ga(ifull),dgdtg(ifull))

return
end subroutine soilsnin_init

subroutine soilsnin_end

implicit none

deallocate(condxpr,fev,fes,ga,dgdtg)

return
end subroutine soilsnin_end

end module soilsnin_m