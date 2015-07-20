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
    
module morepbl_m

implicit none

private
public condx,fg,eg,epot,condc,rnet,pblh,epan,tpan
public conds,condg
public morepbl_init,morepbl_end

real, dimension(:), allocatable, save :: condx,fg,eg,epot,condc,rnet,pblh,epan,tpan
real, dimension(:), allocatable, save :: conds,condg

contains

subroutine morepbl_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(condx(ifull),fg(ifull),eg(ifull),epot(ifull))
allocate(condc(ifull),rnet(ifull),pblh(ifull),epan(ifull))
allocate(tpan(ifull),conds(ifull),condg(ifull))
fg=0.
eg=0.
epot=0.
epan=0.
tpan=0.
rnet=0.
pblh=1000.
condx=0.
condc=0.
conds=0.
condg=0.

return
end subroutine morepbl_init

subroutine morepbl_end

implicit none

deallocate(condx,fg,eg,epot,condc,rnet,pblh,epan,tpan)
deallocate(conds,condg)

return
end subroutine morepbl_end

end module morepbl_m
