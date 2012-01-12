module vvel_m

implicit none

private
public sdot,dpsldt
public vvel_init,vvel_end

real, dimension(:,:), allocatable, save :: sdot,dpsldt

contains

subroutine vvel_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(sdot(ifull,kl+1),dpsldt(ifull,kl))
sdot=0.
dpsldt=0.

return
end subroutine vvel_init

subroutine vvel_end

implicit none

deallocate(sdot,dpsldt)

return
end subroutine vvel_end

end module vvel_m