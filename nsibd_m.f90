! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module nsibd_m

implicit none

private
public rsmin,sigmf,tgf,sigmu
public ivegt,isoilm,isoilm_in,iurbant
public climate_ivegt,climate_biome
public climate_min20,climate_max20,climate_alpha20,climate_agdd5
public climate_gmd,climate_dmoist_min20,climate_dmoist_max20
public nsibd_init,nsibd_end

real, dimension(:), allocatable, save :: rsmin,tgf,sigmu
real, dimension(:), allocatable, save :: sigmf
real, dimension(:), allocatable, save :: climate_min20, climate_max20, climate_alpha20
real, dimension(:), allocatable, save :: climate_agdd5
real, dimension(:), allocatable, save :: climate_dmoist_min20, climate_dmoist_max20
integer, dimension(:), allocatable, save :: ivegt,isoilm,isoilm_in,iurbant
integer, dimension(:), allocatable, save :: climate_ivegt,climate_biome
integer, dimension(:), allocatable, save :: climate_gmd

contains

subroutine nsibd_init(ifull,nsib,cable_climate)

implicit none

integer, intent(in) :: ifull,nsib,cable_climate

allocate(ivegt(ifull),isoilm(ifull))
allocate(iurbant(ifull))
allocate(sigmf(ifull),sigmu(ifull))
allocate(rsmin(ifull),isoilm_in(ifull))
ivegt=0
isoilm=0
iurbant=0
isoilm_in=0
sigmf=0.
sigmu=0.
rsmin=995.
if (nsib==3.or.nsib==5) then
  allocate(tgf(ifull))
  tgf=293.
end if
if (cable_climate==1) then
  allocate(climate_ivegt(ifull),climate_biome(ifull))
  allocate(climate_min20(ifull),climate_max20(ifull),climate_alpha20(ifull))
  allocate(climate_agdd5(ifull))
  allocate(climate_gmd(ifull),climate_dmoist_min20(ifull),climate_dmoist_max20(ifull))
  climate_ivegt = 0
  climate_biome = 0
  climate_min20 = 0.
  climate_max20 = 0.
  climate_alpha20 = 0.
  climate_agdd5 = 0.
  climate_gmd = 0
  climate_dmoist_min20 = 0.
  climate_dmoist_max20 = 0.
end if

return
end subroutine nsibd_init

subroutine nsibd_end

implicit none

deallocate(ivegt,isoilm)
deallocate(iurbant)
deallocate(sigmf,sigmu)
deallocate(rsmin,isoilm_in)
if (allocated(tgf)) then
  deallocate(tgf)
end if
if (allocated(climate_ivegt)) then
  deallocate(climate_ivegt,climate_biome)
  deallocate(climate_min20,climate_max20,climate_alpha20)
  deallocate(climate_agdd5)
  deallocate(climate_gmd,climate_dmoist_min20,climate_dmoist_max20)
end if

return
end subroutine nsibd_end

end module nsibd_m