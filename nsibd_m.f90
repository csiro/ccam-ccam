module nsibd_m

implicit none

private
public rsmin,sigmf,tgf,res,tsigmf,sigmu
public ivegt,isoilm
public nsibd_init,nsibd_end

real, dimension(:), allocatable, save :: rsmin,sigmf,tgf,res,tsigmf,sigmu
integer, dimension(:), allocatable, save :: ivegt,isoilm

contains

subroutine nsibd_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(sigmf(ifull),tgf(ifull),sigmu(ifull),res(ifull),rsmin(ifull),tsigmf(ifull))
allocate(ivegt(ifull),isoilm(ifull))

return
end subroutine nsibd_init

subroutine nsibd_end

implicit none

deallocate(sigmf,tgf,sigmu,res,rsmin,tsigmf)
deallocate(ivegt,isoilm)

return
end subroutine nsibd_end

end module nsibd_m