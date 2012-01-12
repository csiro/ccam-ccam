module bigxy4_m

implicit none

private
public xx4,yy4
public bigxy4_init,bigxy4_end

real*8, dimension(:,:), allocatable, save :: xx4,yy4

contains

subroutine bigxy4_init(iquad)

implicit none

integer, intent(in) :: iquad

allocate(xx4(iquad,iquad),yy4(iquad,iquad))

return
end subroutine bigxy4_init

subroutine bigxy4_end

implicit none

deallocate(xx4,yy4)

return
end subroutine bigxy4_end

end module bigxy4_m