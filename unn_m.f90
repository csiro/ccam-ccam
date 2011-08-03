module unn_m

implicit none

private
public unn,vnn
public unn_init,unn_end

real, dimension(:,:), allocatable, save :: unn,vnn

contains

subroutine unn_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(unn(ifull,kl),vnn(ifull,kl))

return
end subroutine unn_init

subroutine unn_end

implicit none

deallocate(unn,vnn)

return
end subroutine unn_end

end module unn_m