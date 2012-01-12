module screen_m

implicit none

private
public tscrn,qgscrn,uscrn,rhscrn,u10
public screen_init,screen_end

real, dimension(:), allocatable, save :: tscrn,qgscrn,uscrn,rhscrn,u10

contains

subroutine screen_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tscrn(ifull),qgscrn(ifull),uscrn(ifull),rhscrn(ifull),u10(ifull))

return
end subroutine screen_init

subroutine screen_end

implicit none

deallocate(tscrn,qgscrn,uscrn,rhscrn,u10)

return
end subroutine screen_end

end module screen_m