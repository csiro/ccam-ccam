module work2_m

implicit none

private
public zoh,qsttg,wetfac
public rho,zo,theta
public vmod,dgdtg
public work2_init,work2_end

real, dimension(:), allocatable, save :: zoh,qsttg,wetfac
real, dimension(:), allocatable, save :: rho,zo,theta
real, dimension(:), allocatable, save :: vmod,dgdtg

contains

subroutine work2_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib

allocate(zoh(ifull),qsttg(ifull),wetfac(ifull))
allocate(rho(ifull),zo(ifull),theta(ifull))
allocate(vmod(ifull))
rho=1.
zo=0.
zoh=0.
vmod=0.
qsttg=0.
wetfac=0.
if (nsib.eq.3.or.nsib.eq.5) then
  allocate(dgdtg(ifull))
  dgdtg=0.
end if

return
end subroutine work2_init

subroutine work2_end

implicit none

deallocate(zoh,qsttg,wetfac)
deallocate(rho,zo,theta)
deallocate(vmod)
if (allocated(dgdtg)) then
  deallocate(dgdtg)
end if

return
end subroutine work2_end

end module work2_m