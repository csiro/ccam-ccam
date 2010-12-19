module nharrs_m

implicit none

private
public phi,h_nh
public nharrs_init,nharrs_end

real, dimension(:,:), allocatable, save :: phi,h_nh

contains

subroutine nharrs_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(phi(ifull,kl),h_nh(ifull+iextra,kl))

return
end subroutine nharrs_init

subroutine nharrs_end

implicit none

deallocate(phi,h_nh)

return
end subroutine nharrs_end

end module nharrs_m