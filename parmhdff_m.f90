module parmhdff_m

implicit none

private
public nhor,nhorps,khor,khdif,nhorjlm
public hdifmax
public hdiff
public parmhdff_init,parmhdff_end

integer, save :: nhor,nhorps,khor,khdif,nhorjlm
real, save :: hdifmax
real, dimension(:), allocatable, save :: hdiff

contains

subroutine parmhdff_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(hdiff(kl))
! Horizontal mixing options
nhor=-157
nhorps=-1
khor=-8
khdif=2
nhorjlm=1
hdifmax=0.
hdiff=0.1*khdif

return
end subroutine parmhdff_init

subroutine parmhdff_end

implicit none

deallocate(hdiff)

return
end subroutine parmhdff_end

end module parmhdff_m