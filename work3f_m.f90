module work3f_m

implicit none

private
public qccon,qlrad,qfrad
public nface,xg,yg
public work3f_init,work3f_end

real, dimension(:,:), allocatable, save :: qccon,qlrad,qfrad
real, dimension(:,:), allocatable, save :: xg,yg
integer, dimension(:,:), allocatable, save :: nface

contains

subroutine work3f_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl))
allocate(nface(ifull,kl),xg(ifull,kl),yg(ifull,kl))

return
end subroutine work3f_init

subroutine work3f_end

implicit none

deallocate(qccon,qlrad,qfrad)
deallocate(nface,xg,yg)

return
end subroutine work3f_end

end module work3f_m