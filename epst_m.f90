module epst_m

implicit none

private
public epst
public epst_init,epst_end

real, dimension(:), allocatable, save :: epst

contains

subroutine epst_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(epst(ifull))

return
end subroutine epst_init

subroutine epst_end

implicit none

deallocate(epst)

return
end subroutine epst_end

end module epst_m