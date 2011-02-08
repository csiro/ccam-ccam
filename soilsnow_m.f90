module soilsnow_m

implicit none

private
public tggsn,tgg,wb,wbice,wbtot,smass,ssdn,ssdnn,snowd 
public osnowd,snage,sno,gflux,sgflux,snowflx,otgsoil 
public runoff,rnof1,rnof2,rtsoil,albvisnir,albsoil,albsoilsn
public fracice,sicedep
public isflag
public soilsnow_init,soilsnow_end

integer, dimension(:), allocatable, save :: isflag
real, dimension(:), allocatable, save :: wbtot,ssdnn,snowd
real, dimension(:), allocatable, save :: osnowd,snage,sno,gflux,sgflux,snowflx,otgsoil
real, dimension(:), allocatable, save :: runoff,rnof1,rnof2,rtsoil,albsoil
real, dimension(:), allocatable, save :: fracice,sicedep
real, dimension(:,:), allocatable, save :: tggsn,tgg,wb,wbice,smass,ssdn
real, dimension(:,:), allocatable, save :: albvisnir,albsoilsn

contains

subroutine soilsnow_init(ifull,iextra,kl,ms)

implicit none

integer, intent(in) :: ifull,iextra,kl,ms

allocate(tggsn(ifull,3),tgg(ifull,ms),wb(ifull,ms),wbice(ifull,ms),wbtot(ifull))
allocate(smass(ifull,3),ssdn(ifull,3),ssdnn(ifull),snowd(ifull))
allocate(osnowd(ifull),snage(ifull),sno(ifull),gflux(ifull))
allocate(sgflux(ifull),snowflx(ifull),otgsoil(ifull))
allocate(runoff(ifull),rnof1(ifull),rnof2(ifull),rtsoil(ifull))
allocate(albvisnir(ifull,2),albsoil(ifull),albsoilsn(ifull,2))
allocate(fracice(ifull),sicedep(ifull))
allocate(isflag(ifull))

return
end subroutine soilsnow_init

subroutine soilsnow_end

implicit none

deallocate(tggsn,tgg,wb,wbice,wbtot,smass,ssdn,ssdnn,snowd)
deallocate(osnowd,snage,sno,gflux,sgflux,snowflx,otgsoil)
deallocate(runoff,rnof1,rnof2,rtsoil,albvisnir,albsoil,albsoilsn)
deallocate(fracice,sicedep)
deallocate(isflag)

return
end subroutine soilsnow_end

end module soilsnow_m