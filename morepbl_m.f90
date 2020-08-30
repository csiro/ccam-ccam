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
    
module morepbl_m

implicit none

private
public condx,fg,eg,epot,condc,rnet,pblh,epan,tpan
public conds,condg
public anthropogenic_flux, urban_tas, urban_ts, urban_wetfac
public urban_storage_flux, urban_elecgas_flux
public urban_heating_flux, urban_cooling_flux
public urban_zom, urban_zoh, urban_zoq, urban_emiss
public ua150, va150, ua250, va250
public morepbl_init, morepbl_end

#ifdef scm
public wth_flux, wq_flux, uw_flux, vw_flux
public tkesave, epssave, mfsave
public rkmsave, rkhsave
public buoyproduction, shearproduction, totaltransport
#endif

real, dimension(:), allocatable, save :: epot,rnet,epan,tpan
real, dimension(:), allocatable, save :: anthropogenic_flux, urban_tas, urban_ts, urban_wetfac
real, dimension(:), allocatable, save :: urban_storage_flux, urban_elecgas_flux
real, dimension(:), allocatable, save :: urban_heating_flux, urban_cooling_flux
real, dimension(:), allocatable, save :: urban_zom, urban_zoh, urban_zoq, urban_emiss
real, dimension(:), allocatable, save :: condc, condx, conds, condg, pblh, fg, eg
real, dimension(:), allocatable, save :: ua150, va150, ua250, va250

#ifdef scm
real, dimension(:,:), allocatable, save :: wth_flux, wq_flux, uw_flux, vw_flux
real, dimension(:,:), allocatable, save :: tkesave, epssave, mfsave
real, dimension(:,:), allocatable, save :: rkmsave, rkhsave
real, dimension(:,:), allocatable, save :: buoyproduction, shearproduction, totaltransport
#endif

contains

#ifdef scm
subroutine morepbl_init(ifull,kl)
#else
subroutine morepbl_init(ifull,diaglevel_pbl)
#endif

implicit none

#ifdef scm
integer, intent(in) :: ifull, kl
#else
integer, intent(in) :: ifull, diaglevel_pbl
#endif

allocate( condx(ifull), fg(ifull), eg(ifull), epot(ifull) )
allocate( condc(ifull), rnet(ifull), pblh(ifull), epan(ifull) )
allocate( tpan(ifull), conds(ifull), condg(ifull) )
allocate( anthropogenic_flux(ifull), urban_tas(ifull), urban_ts(ifull), urban_wetfac(ifull) )
allocate( urban_storage_flux(ifull), urban_elecgas_flux(ifull) )
allocate( urban_heating_flux(ifull), urban_cooling_flux(ifull) )
allocate( urban_zom(ifull), urban_zoh(ifull), urban_zoq(ifull), urban_emiss(ifull) )

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
anthropogenic_flux = 0.
urban_elecgas_flux = 0.
urban_heating_flux = 0.
urban_cooling_flux = 0.
urban_storage_flux = 0.
urban_tas          = 0.
urban_ts           = 0.
urban_wetfac       = 0.
urban_zom          = 0.
urban_zoh          = 0.
urban_zoq          = 0.
urban_emiss        = 0.

#ifdef scm
allocate( wth_flux(ifull,kl), wq_flux(ifull,kl) )
allocate( uw_flux(ifull,kl), vw_flux(ifull,kl) )
allocate( tkesave(ifull,kl), epssave(ifull,kl), mfsave(ifull,kl-1) )
allocate( rkmsave(ifull,kl), rkhsave(ifull,kl) )
allocate( buoyproduction(ifull,kl), shearproduction(ifull,kl) )
allocate( totaltransport(ifull,kl) )
wth_flux=0.
wq_flux=0.
uw_flux=0.
vw_flux=0.
tkesave=0.
epssave=0.
mfsave=0.
rkmsave=0.
rkhsave=0.
buoyproduction=0.
shearproduction=0.
totaltransport=0.
#else
if ( diaglevel_pbl>3 ) then
  allocate( ua150(ifull), va150(ifull), ua250(ifull), va250(ifull) )
  ua150 = 0.
  va150 = 0.
  ua250 = 0.
  va250 = 0.
end if
#endif

return
end subroutine morepbl_init

subroutine morepbl_end

implicit none

deallocate( condx, fg, eg, epot, condc, rnet, pblh, epan, tpan )
deallocate( conds, condg )
deallocate( anthropogenic_flux, urban_tas, urban_ts, urban_wetfac )
deallocate( urban_storage_flux, urban_elecgas_flux )
deallocate( urban_heating_flux, urban_cooling_flux )
deallocate( urban_zom, urban_zoh, urban_zoq, urban_emiss )

#ifdef scm
deallocate( wth_flux, wq_flux, uw_flux, vw_flux )
deallocate( tkesave, epssave, mfsave )
deallocate( rkmsave, rkhsave )
deallocate( buoyproduction, shearproduction )
deallocate( totaltransport )
#else
deallocate( ua150, va150, ua250, va250 )
#endif

return
end subroutine morepbl_end

end module morepbl_m
