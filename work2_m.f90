module work2_m

implicit none

private
public dirad,dfgdt,degdt,wetfac,degdw,cie
public factch,qsttg,rho,zo,aft,fh,ri,theta
public gamm,rg,vmod,dgdtg
public work2_init,work2_end

real, dimension(:), allocatable, save :: dirad,dfgdt,degdt,wetfac,degdw,cie
real, dimension(:), allocatable, save :: factch,qsttg,rho,zo,aft,fh,ri,theta
real, dimension(:), allocatable, save :: gamm,rg,vmod,dgdtg

contains

subroutine work2_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(dirad(ifull),dfgdt(ifull),degdt(ifull),wetfac(ifull),degdw(ifull),cie(ifull))
allocate(factch(ifull),qsttg(ifull),rho(ifull),zo(ifull),aft(ifull),fh(ifull),ri(ifull),theta(ifull))
allocate(gamm(ifull),rg(ifull),vmod(ifull),dgdtg(ifull))
wetfac=0.
rho=1.
zo=0.
gamm=0.
rg=0.
vmod=0.

return
end subroutine work2_init

subroutine work2_end

implicit none

deallocate(dirad,dfgdt,degdt,wetfac,degdw,cie)
deallocate(factch,qsttg,rho,zo,aft,fh,ri,theta)
deallocate(gamm,rg,vmod,dgdtg)

return
end subroutine work2_end

end module work2_m