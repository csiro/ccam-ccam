module work3b_m

implicit none

private
public wblf,wbfice,sdepth
public work3b_init,work3b_end

real, dimension(:,:), allocatable, save :: wblf,wbfice,sdepth

contains

subroutine work3b_init(ifull,iextra,kl,ms)

implicit none

integer, intent(in) :: ifull,iextra,kl,ms

allocate(wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3))

return
end subroutine work3b_init

subroutine work3b_end

implicit none

deallocate(wblf,wbfice,sdepth)

return
end subroutine work3b_end

end module work3b_m