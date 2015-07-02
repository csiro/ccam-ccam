! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module carbpools_m

use casadimension, only : mplant,mlitter,msoil ! CASA dimensions
use cable_def_types_mod, only : ncs,ncp        ! CABLE dimensions

implicit none

private
public inyear_carb
public fnee,fpn,frd,frp,frpw,frpr,frs
public cplant,clitter,csoil,niplant,nisoil,nilitter
public pplant,plitter,psoil,glai
public carbpools_init,carbpools_end

integer, save :: inyear_carb
real, dimension(:), allocatable, save :: fnee,fpn,frd,frp,frpw,frpr,frs
real, dimension(:,:), allocatable, save :: cplant,clitter,csoil,niplant,nilitter,nisoil
real, dimension(:,:), allocatable, save :: pplant,plitter,psoil
real, dimension(:), allocatable, save :: glai

contains

subroutine carbpools_init(ifull,iextra,kl,nsib,ccycle)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib,ccycle

allocate(fnee(ifull),fpn(ifull),frd(ifull),frp(ifull))
allocate(frpw(ifull),frpr(ifull),frs(ifull))
fnee=0.
fpn=0.
frd=0.
frp=0.
frpw=0.
frpr=0.
frs=0.
if (nsib==4.or.nsib>=6) then
  if (ccycle==0) then
    allocate(cplant(ifull,ncp),csoil(ifull,ncs))
    cplant=0.
    csoil=0.
  else
    allocate(cplant(ifull,mplant),clitter(ifull,mlitter),csoil(ifull,msoil))
    allocate(niplant(ifull,mplant),nilitter(ifull,mlitter),nisoil(ifull,msoil))
    allocate(pplant(ifull,mplant),plitter(ifull,mlitter),psoil(ifull,msoil))
    allocate(glai(ifull))
    cplant=0.
    clitter=0.
    csoil=0.
    niplant=0.
    nilitter=0.
    nisoil=0.
    pplant=0.
    plitter=0.
    psoil=0.
    glai=0.
  end if
end if

return
end subroutine carbpools_init

subroutine carbpools_end

implicit none

deallocate(fnee,fpn,frd,frp,frpw,frpr,frs)
if (allocated(cplant)) then
  deallocate(cplant,csoil)
end if
if (allocated(clitter)) then
  deallocate(clitter,niplant,nilitter,nisoil)
  deallocate(pplant,plitter,psoil,glai)
end if

return
end subroutine carbpools_end

end module carbpools_m
