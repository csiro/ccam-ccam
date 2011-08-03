module kdacom_m

implicit none

private
public qh2o,p,delp2,delp,t,var1,var2
public var3,var4,cntval
public kdacom_init,kdacom_end

real, dimension(:,:), allocatable, save :: qh2o,p,delp2,delp,t,var1,var2
real, dimension(:,:), allocatable, save :: var3,var4,cntval

contains

subroutine kdacom_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(qh2o(imax,kl+1),p(imax,kl+1),delp2(imax,kl),delp(imax,kl),t(imax,kl+1),var1(imax,kl),var2(imax,kl))
allocate(var3(imax,kl),var4(imax,kl),cntval(imax,kl+1))

return
end subroutine kdacom_init

subroutine kdacom_end

implicit none

deallocate(qh2o,p,delp2,delp,t,var1,var2)
deallocate(var3,var4,cntval)

return
end subroutine kdacom_end

end module kdacom_m