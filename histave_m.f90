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
    
module histave_m

implicit none

private
public eg_ave,fg_ave,ga_ave,epan_ave,dew_ave
public cbas_ave,ctop_ave,rndmax,qscrn_ave
public tmaxscr,tminscr,tscr_ave
public rhmaxscr,rhminscr
public riwp_ave,rlwp_ave,u10max,v10max,u10mx
public u1max,v1max,u2max,v2max,cape_max,cape_ave,epot_ave
public rnet_ave,mixdep_ave
public wb_ave,wbice_ave,tsu_ave,alb_ave,fbeam_ave,psl_ave,convh_ave
public fpn_ave,frs_ave,frp_ave
!public tgg_ave
public histave_init,histave_end

real, dimension(:), allocatable, save :: eg_ave,fg_ave,ga_ave,epan_ave,dew_ave
real, dimension(:), allocatable, save :: cbas_ave,ctop_ave,rndmax,qscrn_ave
real, dimension(:), allocatable, save :: tmaxscr,tminscr,tscr_ave
real, dimension(:), allocatable, save :: rhmaxscr,rhminscr
real, dimension(:), allocatable, save :: riwp_ave,rlwp_ave,u10max,v10max,u10mx
real, dimension(:), allocatable, save :: u1max,v1max,u2max,v2max,cape_max,cape_ave,epot_ave
real, dimension(:), allocatable, save :: rnet_ave,mixdep_ave
real, dimension(:,:), allocatable, save :: wb_ave,wbice_ave,convh_ave
real, dimension(:), allocatable, save :: tsu_ave,alb_ave,fbeam_ave,psl_ave
real, dimension(:), allocatable, save :: fpn_ave,frs_ave,frp_ave
!real, dimension(:,:), allocatable, save :: tgg_ave

contains

subroutine histave_init(ifull,kl,ms,ccycle)

implicit none

integer, intent(in) :: ifull,kl,ms,ccycle

allocate(eg_ave(ifull),fg_ave(ifull),ga_ave(ifull),epan_ave(ifull),dew_ave(ifull))
allocate(cbas_ave(ifull),ctop_ave(ifull),rndmax(ifull),qscrn_ave(ifull))
allocate(tmaxscr(ifull),tminscr(ifull),tscr_ave(ifull))
allocate(rhmaxscr(ifull),rhminscr(ifull))
allocate(riwp_ave(ifull),rlwp_ave(ifull),u10max(ifull),v10max(ifull),u10mx(ifull))
allocate(u1max(ifull),v1max(ifull),u2max(ifull),v2max(ifull),cape_max(ifull),cape_ave(ifull),epot_ave(ifull))
allocate(rnet_ave(ifull),mixdep_ave(ifull))
allocate(wb_ave(ifull,ms),wbice_ave(ifull,ms),tsu_ave(ifull),alb_ave(ifull),fbeam_ave(ifull),psl_ave(ifull),convh_ave(ifull,kl))
!allocate(tgg_ave(ifull,ms))

! needs to be initialised here for zeroth time-step in outcdf.f90
rndmax(:)      = 0.
tmaxscr(:)     = 0.
tminscr(:)     = 400.
rhmaxscr(:)    = 0.
rhminscr(:)    = 400.
u10max(:)      = 0.
v10max(:)      = 0.
u1max(:)       = 0.
v1max(:)       = 0.
u2max(:)       = 0.
v2max(:)       = 0.
cape_max(:)    = 0.
cape_ave(:)    = 0.
u10mx(:)       = 0.
tscr_ave(:)    = 0.
qscrn_ave(:)   = 0.
dew_ave(:)     = 0.
epan_ave(:)    = 0.
epot_ave(:)    = 0.
eg_ave(:)      = 0.
fg_ave(:)      = 0.
ga_ave(:)      = 0.
rnet_ave(:)    = 0.
riwp_ave(:)    = 0.
rlwp_ave(:)    = 0.
convh_ave(:,:) = 0.
cbas_ave(:)    = 0.
ctop_ave(:)    = 0.
wb_ave(:,:)    = 0.
wbice_ave(:,:) = 0.
tsu_ave(:)     = 0.
alb_ave(:)     = 0.
fbeam_ave(:)   = 0.
psl_ave(:)     = 0.
mixdep_ave(:)  = 0.

if ( ccycle/=0 ) then
  allocate(fpn_ave(ifull),frs_ave(ifull),frp_ave(ifull))
  fpn_ave(:)     = 0.
  frs_ave(:)     = 0.
  frp_ave(:)     = 0.
end if

return
end subroutine histave_init

subroutine histave_end

implicit none

deallocate(eg_ave,fg_ave,ga_ave,epan_ave,dew_ave)
deallocate(cbas_ave,ctop_ave,rndmax,qscrn_ave)
deallocate(tmaxscr,tminscr,tscr_ave)
deallocate(rhmaxscr,rhminscr)
deallocate(riwp_ave,rlwp_ave,u10max,v10max,u10mx)
deallocate(u1max,v1max,u2max,v2max,cape_max,cape_ave,epot_ave)
deallocate(rnet_ave,mixdep_ave)
deallocate(wb_ave,wbice_ave,tsu_ave,alb_ave,fbeam_ave,psl_ave,convh_ave)
!deallocate(tgg_ave)

if ( allocated(fpn_ave) ) then
  deallocate(fpn_ave,frs_ave,frp_ave)
end if

return
end subroutine histave_end

end module histave_m