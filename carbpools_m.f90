module carbpools_m

use define_dimensions, only : ncs,ncp   ! CABLE dimensions

implicit none

private
public inyear_carb
public fnee,fpn,frd,frp,frpw,frpr,frs
public cplant,csoil
public carbpools_init,carbpools_end

integer, save :: inyear_carb
real, dimension(:), allocatable, save :: fnee,fpn,frd,frp,frpw,frpr,frs
real, dimension(:,:), allocatable, save :: cplant,csoil

contains

subroutine carbpools_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib

allocate(fnee(ifull),fpn(ifull),frd(ifull),frp(ifull))
allocate(frpw(ifull),frpr(ifull),frs(ifull))
fnee=0.
fpn=0.
frd=0.
frp=0.
frpw=0.
frpr=0.
frs=0.
if (nsib.eq.4.or.nsib.ge.6) then
  allocate(cplant(ifull,ncp),csoil(ifull,ncs))
  cplant=0.
  csoil=0.
end if

return
end subroutine carbpools_init

subroutine carbpools_end

implicit none

deallocate(fnee,fpn,frd,frp,frpw,frpr,frs)
if (allocated(cplant)) then
  deallocate(cplant,csoil)
end if

return
end subroutine carbpools_end

end module carbpools_m