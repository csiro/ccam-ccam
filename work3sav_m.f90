module work3sav_m

implicit none

private
public qgsav,qfgsav,qlgsav,trsav
public work3sav_init,work3sav_end

real, dimension(:,:), allocatable, save :: qgsav,qfgsav,qlgsav
real, dimension(:,:,:), allocatable, save :: trsav

contains

subroutine work3sav_init(ifull,iextra,kl,ilt,jlt,klt,ngasmax)

implicit none

integer, intent(in) :: ifull,iextra,kl,ilt,jlt,klt,ngasmax

allocate(qgsav(ifull,kl),qfgsav(ifull,kl),qlgsav(ifull,kl))
if (ilt.gt.0) allocate(trsav(ilt*jlt,klt,ngasmax))

return
end subroutine work3sav_init

subroutine work3sav_end

implicit none

deallocate(qgsav,qfgsav,qlgsav)
if (allocated(trsav)) deallocate(trsav)

return
end subroutine work3sav_end

end module work3sav_m