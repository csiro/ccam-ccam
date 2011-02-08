module parmhdff_m

implicit none

private
public hdiff
public parmhdff_init,parmhdff_end

real, dimension(:), allocatable, save :: hdiff

contains

subroutine parmhdff_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(hdiff(kl))

return
end subroutine parmhdff_init

subroutine parmhdff_end

implicit none

deallocate(hdiff)

return
end subroutine parmhdff_end

end module parmhdff_m