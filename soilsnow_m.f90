module soilsnow_m

implicit none

private
public tggsn,tgg,wb,wbice,smass,ssdn,ssdnn,snowd 
public osnowd,snage,sno,gflux,sgflux,snowflx,otgsoil 
public runoff,albvisnir
public fracice,sicedep
public isflag
public soilsnow_init,soilsnow_end

integer, dimension(:), allocatable, save :: isflag
real, dimension(:), allocatable, save :: ssdnn,snowd
real, dimension(:), allocatable, save :: osnowd,snage,sno,gflux,sgflux,snowflx,otgsoil
real, dimension(:), allocatable, save :: runoff
real, dimension(:), allocatable, save :: fracice,sicedep
real, dimension(:,:), allocatable, save :: tggsn,tgg,wb,wbice,smass,ssdn
real, dimension(:,:), allocatable, save :: albvisnir

contains

subroutine soilsnow_init(ifull,iextra,kl,ms,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,ms,nsib

allocate(tggsn(ifull,3),tgg(ifull,ms),wb(ifull,ms),wbice(ifull,ms))
allocate(smass(ifull,3),ssdn(ifull,3),ssdnn(ifull),snowd(ifull))
allocate(snage(ifull),sno(ifull),gflux(ifull))
allocate(runoff(ifull),albvisnir(ifull,2))
allocate(fracice(ifull),sicedep(ifull))
allocate(isflag(ifull))
if (nsib.eq.3.or.nsib.eq.5) then
  allocate(sgflux(ifull))
  allocate(osnowd(ifull),otgsoil(ifull),snowflx(ifull))
end if

return
end subroutine soilsnow_init

subroutine soilsnow_end

implicit none

deallocate(tggsn,tgg,wb,wbice,smass,ssdn,ssdnn,snowd)
deallocate(snage,sno,gflux)
deallocate(runoff,albvisnir)
deallocate(fracice,sicedep)
deallocate(isflag)
if (allocated(sgflux)) then
  deallocate(sgflux)
  deallocate(osnowd,otgsoil,snowflx)
end if

return
end subroutine soilsnow_end

end module soilsnow_m
