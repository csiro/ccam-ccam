! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public anthropogenic_flux
public rkmsave, rkhsave
public morepbl_init,morepbl_end

#ifdef scm
public wth_flux, wq_flux, uw_flux, vw_flux
public tkesave, mfsave
#endif

real, dimension(:), allocatable, save :: condx,fg,eg,epot,condc,rnet,pblh,epan,tpan
real, dimension(:), allocatable, save :: conds,condg
real, dimension(:), allocatable, save :: anthropogenic_flux
real, dimension(:,:), allocatable, save :: rkmsave, rkhsave

#ifdef scm
real, dimension(:,:), allocatable, save :: wth_flux, wq_flux, uw_flux, vw_flux
real, dimension(:,:), allocatable, save :: tkesave, mfsave
#endif

contains

subroutine morepbl_init(ifull,kl)

implicit none

integer, intent(in) :: ifull, kl

allocate( condx(ifull), fg(ifull), eg(ifull), epot(ifull) )
allocate( condc(ifull), rnet(ifull), pblh(ifull), epan(ifull) )
allocate( tpan(ifull), conds(ifull), condg(ifull) )
allocate( anthropogenic_flux(ifull) )
allocate( rkmsave(ifull,kl), rkhsave(ifull,kl) )

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
rkmsave=0.
rkhsave=0.

#ifdef scm
allocate( wth_flux(ifull,kl), wq_flux(ifull,kl) )
allocate( uw_flux(ifull,kl), vw_flux(ifull,kl) )
allocate( tkesave(ifull,kl), mfsave(ifull,kl-1) )
wth_flux=0.
wq_flux=0.
uw_flux=0.
vw_flux=0.
tkesave=0.
mfsave=0.
#endif

return
end subroutine morepbl_init

subroutine morepbl_end

implicit none

deallocate( condx, fg, eg, epot, condc, rnet, pblh, epan, tpan )
deallocate( conds, condg )
deallocate( anthropogenic_flux )
deallocate( rkmsave, rkhsave )

#ifdef scm
deallocate( wth_flux, wq_flux, uw_flux, vw_flux )
deallocate( tkesave, mfsave )
#endif

return
end subroutine morepbl_end

end module morepbl_m
