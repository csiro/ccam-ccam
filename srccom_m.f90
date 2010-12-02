module srccom_m

implicit none

private
public sorc,csour1,csour2,osour,csour,ss1
public srccom_init,srccom_end

real, dimension(:,:), allocatable, save :: csour1,csour2,osour,csour,ss1
real, dimension(:,:,:), allocatable, save :: sorc

contains

subroutine srccom_init(ifull,iextra,kl,imax,nbly)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax,nbly

allocate(sorc(imax,kl+1,nbly),csour1(imax,kl+1),csour2(imax,kl+1))
allocate(osour(imax,kl+1),csour(imax,kl+1),ss1(imax,kl+1))

return
end subroutine srccom_init

subroutine srccom_end

implicit none

deallocate(sorc,csour1,csour2,osour,csour,ss1)

return
end subroutine srccom_end

end module srccom_m