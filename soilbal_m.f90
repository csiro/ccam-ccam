module soilbal_m

implicit none

private
public tevap,tprecip,trnoff,totenbal,osnowd0,wbtot0
public soilbal_init,soilbal_end

real, dimension(:), allocatable, save :: tevap,tprecip,trnoff,totenbal,osnowd0,wbtot0

contains

subroutine soilbal_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tevap(ifull),tprecip(ifull),trnoff(ifull),totenbal(ifull),osnowd0(ifull),wbtot0(ifull))

return
end subroutine soilbal_init

subroutine soilbal_end

implicit none

deallocate(tevap,tprecip,trnoff,totenbal,osnowd0,wbtot0)

return
end subroutine soilbal_end

end module soilbal_m