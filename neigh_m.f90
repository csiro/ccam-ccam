module neigh_m

implicit none

private
public neigh
public neigh_init,neigh_end

real, dimension(:), allocatable, save :: neigh

contains

subroutine neigh_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(neigh(ifull))

return
end subroutine neigh_init

subroutine neigh_end

implicit none

deallocate(neigh)

return
end subroutine neigh_end

end module neigh_m