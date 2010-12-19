module workglob_m

implicit none

private
public rlong4,rlat4
public workglob_init,workglob_end

real, dimension(:,:), allocatable, save :: rlong4,rlat4

contains

subroutine workglob_init(ifull_g)

implicit none

integer, intent(in) :: ifull_g

if (.not.allocated(rlong4)) allocate(rlong4(ifull_g,4))
if (.not.allocated(rlat4)) allocate(rlat4(ifull_g,4))

return
end subroutine workglob_init

subroutine workglob_end

implicit none

deallocate(rlong4,rlat4)

return
end subroutine workglob_end

end module workglob_m