module nharrs_m

implicit none

private
public phi,phi_nh,h_nh,lrestart
public nharrs_init,nharrs_end

real, dimension(:,:), allocatable, save :: phi,phi_nh,h_nh
logical, save :: lrestart

contains

subroutine nharrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(phi(ifull,kl),phi_nh(ifull,kl),h_nh(ifull+iextra,kl))
phi=-999.
phi_nh=0.
lrestart=.false.

return
end subroutine nharrs_init

subroutine nharrs_end

implicit none

deallocate(phi,phi_nh,h_nh)

return
end subroutine nharrs_end

end module nharrs_m