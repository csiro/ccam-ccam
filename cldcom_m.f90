module cldcom_m

implicit none

private
public cldfac
public cldcom_init,cldcom_end

real, dimension(:,:,:), allocatable, save :: cldfac

contains

subroutine cldcom_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(cldfac(imax,kl+1,kl+1))

return
end subroutine cldcom_init

subroutine cldcom_end

implicit none

deallocate(cldfac)

return
end subroutine cldcom_end

end module cldcom_m