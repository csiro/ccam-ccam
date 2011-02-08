module tbar2d_m

implicit none

private
public tbar2d
public tbar2d_init,tbar2d_end

real, dimension(:), allocatable, save :: tbar2d

contains

subroutine tbar2d_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tbar2d(ifull))

return
end subroutine tbar2d_init

subroutine tbar2d_end

implicit none

deallocate(tbar2d)

return
end subroutine tbar2d_end

end module tbar2d_m