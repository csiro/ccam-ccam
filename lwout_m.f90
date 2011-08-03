module lwout_m

implicit none

private
public heatra,exctsclr,ctso3clr
public grnflx,topflx,grnflxclr
public lwout_init,lwout_end

real, dimension(:,:), allocatable, save :: heatra,exctsclr,ctso3clr
real, dimension(:), allocatable, save :: grnflx,topflx,grnflxclr

contains

subroutine lwout_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(heatra(imax,kl),exctsclr(imax,kl),ctso3clr(imax,kl))
allocate(grnflx(imax),topflx(imax),grnflxclr(imax))

return
end subroutine lwout_init

subroutine lwout_end

implicit none

deallocate(heatra,exctsclr,ctso3clr)
deallocate(grnflx,topflx,grnflxclr)

return
end subroutine lwout_end

end module lwout_m