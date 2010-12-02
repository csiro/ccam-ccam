module swocom_m

implicit none

private
public fsw,dfsw,ufsw,hsw
public swocom_init,swocom_end

real, dimension(:,:), allocatable, save :: fsw,dfsw,ufsw,hsw

contains

subroutine swocom_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(fsw(imax,kl+1),dfsw(imax,kl+1),ufsw(imax,kl+1),hsw(imax,kl))

return
end subroutine swocom_init

subroutine swocom_end

implicit none

deallocate(fsw,dfsw,ufsw,hsw)

return
end subroutine swocom_end

end module swocom_m