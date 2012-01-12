module gdrag_m

implicit none

private
public he,helo
public gdrag_init,gdrag_end

real, dimension(:), allocatable, save :: he,helo

contains

subroutine gdrag_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(he(ifull),helo(ifull))

return
end subroutine gdrag_init

subroutine gdrag_end

implicit none

deallocate(he,helo)

return
end subroutine gdrag_end

end module gdrag_m