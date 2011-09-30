module savuv1_m

implicit none

private
public savs1,savu1,savv1
public savuv1_init,savuv1_end

real, dimension(:,:), allocatable, save :: savs1,savu1,savv1

contains

subroutine savuv1_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl))

return
end subroutine savuv1_init

subroutine savuv1_end

implicit none

deallocate(savs1,savu1,savv1)

return
end subroutine savuv1_end

end module savuv1_m