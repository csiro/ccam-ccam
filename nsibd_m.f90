module nsibd_m

implicit none

private
public rsmin,sigmf,tgf,sigmu
public ivegt,isoilm
public nsibd_init,nsibd_end

real, dimension(:), allocatable, save :: rsmin,sigmf,tgf,sigmu
integer, dimension(:), allocatable, save :: ivegt,isoilm

contains

subroutine nsibd_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib


allocate(ivegt(ifull),isoilm(ifull))
allocate(sigmf(ifull),sigmu(ifull))
sigmf=0.
sigmu=0.
if (nsib.eq.3.or.nsib.eq.5) then
  allocate(rsmin(ifull),tgf(ifull))
  rsmin=0.
  tgf=293.
end if

return
end subroutine nsibd_init

subroutine nsibd_end

implicit none

deallocate(ivegt,isoilm)
deallocate(sigmf,sigmu)
if (allocated(rsmin)) then
  deallocate(rsmin,tgf)
end if

return
end subroutine nsibd_end

end module nsibd_m
