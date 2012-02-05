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

subroutine carbpools_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(fnee(ifull),fpn(ifull),frd(ifull),frp(ifull))
allocate(frpw(ifull),frpr(ifull),frs(ifull))
allocate(cplant(ifull,ncp),csoil(ifull,ncs))
fnee=0.
fpn=0.
frd=0.
frp=0.
frpw=0.
frpr=0.
frs=0.
cplant=0.
csoil=0.

return
end subroutine carbpools_init

subroutine carbpools_end

implicit none

deallocate(fnee,fpn,frd,frp,frpw,frpr,frs)
deallocate(cplant,csoil)

return
end subroutine carbpools_end

end module carbpools_m