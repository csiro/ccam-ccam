module nharrs_m

implicit none

private
public phi,h_nh,lrestart
public nharrs_init,nharrs_end

real, dimension(:,:), allocatable, save :: phi,h_nh
logical, save :: lrestart

contains

subroutine nharrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(phi(ifull,kl),h_nh(ifull+iextra,kl))
lrestart=.false.

return
end subroutine nharrs_init

subroutine nharrs_end

implicit none

deallocate(phi,h_nh)

return
end subroutine nharrs_end

end module nharrs_m