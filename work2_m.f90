module work2_m

implicit none

private
public zoh,zoq,qsttg,wetfac
public zo,theta
public vmod,dgdtg
public work2_init,work2_end

real, dimension(:), allocatable, save :: zoh,zoq,qsttg,wetfac
real, dimension(:), allocatable, save :: zo,theta
real, dimension(:), allocatable, save :: vmod,dgdtg

contains

subroutine work2_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib

allocate(zoh(ifull),zoq(ifull),qsttg(ifull),wetfac(ifull))
allocate(zo(ifull),theta(ifull))
allocate(vmod(ifull))
zo=0.
zoh=0.
zoq=0.
vmod=0.
qsttg=0.
wetfac=1.
if (nsib==3.or.nsib==5) then
  allocate(dgdtg(ifull))
  dgdtg=0.
end if

return
end subroutine work2_init

subroutine work2_end

implicit none

deallocate(zoh,zoq,qsttg,wetfac)
deallocate(zo,theta)
deallocate(vmod)
if (allocated(dgdtg)) then
  deallocate(dgdtg)
end if

return
end subroutine work2_end

end module work2_m