module dava_m

implicit none

private
public davt,davu
public dava_init,dava_end

real, dimension(:), allocatable, save :: davt,davu

contains

subroutine dava_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(davt(ifull),davu(ifull))

return
end subroutine dava_init

subroutine dava_end

implicit none

deallocate(davt,davu)

return
end subroutine dava_end

end module dava_m