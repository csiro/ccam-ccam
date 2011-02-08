module savuvt_m

implicit none

private
public savt,savpsl
public savs,savu,savv
public savuvt_init,savuvt_end

real, dimension(:), allocatable, save :: savpsl
real, dimension(:,:), allocatable, save :: savt,savs,savu,savv

contains

subroutine savuvt_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(savt(ifull,kl),savpsl(ifull))
allocate(savs(ifull,2:kl),savu(ifull,kl),savv(ifull,kl))

return
end subroutine savuvt_init

subroutine savuvt_end

implicit none

deallocate(savt,savpsl)
deallocate(savs,savu,savv)

return
end subroutine savuvt_end

end module savuvt_m