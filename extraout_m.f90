! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module extraout_m

implicit none

private
public cloudlo,cloudmi,cloudhi,cloudtot
public rgsave,sgsave
public taux,tauy,taux_ave,tauy_ave
public ustar, tstar, qstar, thetavstar
public swrsave,fbeamvis,fbeamnir
public u10_3hr,v10_3hr,tscr_3hr,rh1_3hr
public extraout_init,extraout_end

real, dimension(:), allocatable, save :: cloudlo,cloudmi,cloudhi,cloudtot
real, dimension(:), allocatable, save :: rgsave
real, dimension(:), allocatable, save :: taux,tauy,taux_ave,tauy_ave
real, dimension(:), allocatable, save :: tstar, qstar, thetavstar
real, dimension(:), allocatable, save :: swrsave,fbeamvis,fbeamnir
real, dimension(:), allocatable, save :: sgsave, ustar
real, dimension(:,:), allocatable, save :: u10_3hr,v10_3hr,tscr_3hr,rh1_3hr

contains

subroutine extraout_init(ifull,nextout)

implicit none

integer, intent(in) :: ifull,nextout

allocate(cloudlo(ifull),cloudmi(ifull),cloudhi(ifull),cloudtot(ifull))
allocate(rgsave(ifull),sgsave(ifull))
allocate(taux(ifull),tauy(ifull),taux_ave(ifull),tauy_ave(ifull))
allocate( ustar(ifull), tstar(ifull), qstar(ifull), thetavstar(ifull) )
allocate(swrsave(ifull),fbeamvis(ifull),fbeamnir(ifull))
if (nextout>=2) then
  allocate(u10_3hr(ifull,8),v10_3hr(ifull,8),tscr_3hr(ifull,8),rh1_3hr(ifull,8))
end if

sgsave(:) = 0.
rgsave(:) = 0.
fbeamvis(:) = 0.
fbeamnir(:) = 0.
swrsave(:) = 0.
ustar(:) = 0.
tstar(:) = 0.
qstar(:) = 0.
thetavstar(:) = 0.
taux(:) = 0.
tauy(:) = 0.
taux_ave(:) = 0.
tauy_ave(:) = 0.

return
end subroutine extraout_init

subroutine extraout_end

implicit none

deallocate(cloudlo,cloudmi,cloudhi,cloudtot)
deallocate(rgsave,sgsave)
deallocate(taux,tauy,taux_ave,tauy_ave,ustar)
deallocate(swrsave,fbeamvis,fbeamnir)
if (allocated(u10_3hr)) then
  deallocate(u10_3hr,v10_3hr,tscr_3hr,rh1_3hr)
end if

return
end subroutine extraout_end

end module extraout_m