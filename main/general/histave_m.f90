! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
use parm_m, only : diaglevel_carbon

implicit none

private
public eg_ave,fg_ave,ga_ave,epan_ave,dew_ave
public cbas_ave,ctop_ave,rndmax,prhmax,prhour
public tmaxscr,tminscr,tscr_ave
public rhmaxscr,rhminscr,rhscr_ave
public u10max,v10max
public u1max,v1max,u2max,v2max,cape_max,cape_ave,epot_ave
public rnet_ave
public wb_ave,wbice_ave,convh_ave
public fnee_ave,fpn_ave,frd_ave,frp_ave,frpw_ave,frpr_ave,frs_ave
public cnpp_ave,cnbp_ave
public anthropogenic_ave, tmaxurban, tminurban
public anth_elecgas_ave, anth_heating_ave, anth_cooling_ave
public u_max, v_max, u10m_max, v10m_max ! sub-daily maximums
public fevc_ave,plant_turnover_ave,plant_turnover_wood_ave
public hailradave_ave, hailradmax_max
public histave_init,histave_end


real, dimension(:), allocatable, save :: cbas_ave,ctop_ave,rndmax,prhmax
real, dimension(:), allocatable, save :: tmaxscr,tminscr,tscr_ave
real, dimension(:), allocatable, save :: rhmaxscr,rhminscr,rhscr_ave
!real, dimension(:), allocatable, save :: riwp_ave,rlwp_ave
real, dimension(:), allocatable, save :: u10max,v10max
real, dimension(:), allocatable, save :: u1max,v1max,u2max,v2max,cape_max,cape_ave
real, dimension(:,:), allocatable, save :: wb_ave,wbice_ave,convh_ave
real, dimension(:), allocatable, save :: fnee_ave,fpn_ave,frd_ave,frp_ave,frpw_ave,frpr_ave,frs_ave
real, dimension(:), allocatable, save :: cnpp_ave,cnbp_ave
real, dimension(:), allocatable, save :: tmaxurban, tminurban
real, dimension(:), allocatable, save :: fevc_ave,plant_turnover_ave,plant_turnover_wood_ave
real, dimension(:,:), allocatable, save :: u_max, v_max     ! sub-daily maximums
real, dimension(:), allocatable, save :: u10m_max, v10m_max ! sub-daily maximums
real, dimension(:), allocatable, save :: hailradave_ave, hailradmax_max
real(kind=8), dimension(:), allocatable, save :: eg_ave, fg_ave
real(kind=8), dimension(:), allocatable, save :: epot_ave, prhour
real(kind=8), dimension(:), allocatable, save :: ga_ave, epan_ave, dew_ave
real(kind=8), dimension(:), allocatable, save :: rnet_ave
real(kind=8), dimension(:), allocatable, save :: anthropogenic_ave
real(kind=8), dimension(:), allocatable, save :: anth_elecgas_ave, anth_heating_ave, anth_cooling_ave

contains

subroutine histave_init(ifull,kl,ms,ccycle,output_windmax)

implicit none

integer, intent(in) :: ifull,kl,ms
integer, intent(in) :: output_windmax,ccycle

allocate(eg_ave(ifull),fg_ave(ifull),ga_ave(ifull),epan_ave(ifull),dew_ave(ifull))
allocate(cbas_ave(ifull),ctop_ave(ifull),rndmax(ifull),prhmax(ifull),prhour(ifull))
allocate(tmaxscr(ifull),tminscr(ifull),tscr_ave(ifull))
allocate(rhmaxscr(ifull),rhminscr(ifull),rhscr_ave(ifull))
!allocate(riwp_ave(ifull),rlwp_ave(ifull))
allocate(u10max(ifull),v10max(ifull))
allocate(u1max(ifull),v1max(ifull),u2max(ifull),v2max(ifull),cape_max(ifull),cape_ave(ifull),epot_ave(ifull))
allocate(rnet_ave(ifull))
allocate(wb_ave(ifull,ms),wbice_ave(ifull,ms),convh_ave(ifull,kl))
allocate(anthropogenic_ave(ifull), tmaxurban(ifull), tminurban(ifull))
allocate(anth_elecgas_ave(ifull), anth_heating_ave(ifull), anth_cooling_ave(ifull))
allocate(hailradave_ave(ifull), hailradmax_max(ifull))
!allocate(tgg_ave(ifull,ms))

! needs to be initialised here for zeroth time-step in outcdf.f90
rndmax(:)      = 0.
prhmax(:)      = 0.
prhour(:)      = 0._8
tmaxscr(:)     = 0.
tminscr(:)     = 400.
tscr_ave(:)    = 0.
rhmaxscr(:)    = 0.
rhminscr(:)    = 400.
rhscr_ave(:)   = 0.
u10max(:)      = 0.
v10max(:)      = 0.
u1max(:)       = 0.
v1max(:)       = 0.
u2max(:)       = 0.
v2max(:)       = 0.
cape_max(:)    = 0.
cape_ave(:)    = 0.
dew_ave(:)     = 0._8
epan_ave(:)    = 0._8
epot_ave(:)    = 0._8
eg_ave(:)      = 0._8
fg_ave(:)      = 0._8
ga_ave(:)      = 0._8
rnet_ave(:)    = 0._8
!riwp_ave(:)    = 0.
!rlwp_ave(:)    = 0.
convh_ave(:,:) = 0.
cbas_ave(:)    = 0.
ctop_ave(:)    = 0.
wb_ave(:,:)    = 0.
wbice_ave(:,:) = 0.
anthropogenic_ave(:) = 0._8
anth_elecgas_ave(:)  = 0._8
anth_heating_ave(:)  = 0._8
anth_cooling_ave(:)  = 0._8
tmaxurban(:)   = 0.
tminurban(:)   = 400.
hailradave_ave(:) = 0.
hailradmax_max(:) = 0.

if ( ccycle/=0 ) then
  allocate(fnee_ave(ifull))  
  allocate(fpn_ave(ifull),frd_ave(ifull))
  allocate(frp_ave(ifull),frpw_ave(ifull),frpr_ave(ifull),frs_ave(ifull))
  allocate(cnpp_ave(ifull),cnbp_ave(ifull))
  if ( diaglevel_carbon > 0 ) then
    allocate(fevc_ave(ifull))
    allocate(plant_turnover_ave(ifull))
    allocate(plant_turnover_wood_ave(ifull))
  end if
  fnee_ave(:)    = 0.
  fpn_ave(:)     = 0.
  frd_ave(:)     = 0.
  frp_ave(:)     = 0.
  frpw_ave(:)    = 0.
  frpr_ave(:)    = 0.
  frs_ave(:)     = 0.
  cnpp_ave(:)    = 0.
  cnbp_ave(:)    = 0.
  if ( diaglevel_carbon > 0 ) then
    fevc_ave(:)    = 0.
    plant_turnover_ave = 0.
    plant_turnover_wood_ave = 0.
  end if
end if

if ( output_windmax/=0 ) then
  allocate( u_max(ifull,kl), v_max(ifull,kl) )
  allocate( u10m_max(ifull), v10m_max(ifull) )
  u_max = 0.
  v_max = 0.
  u10m_max = 0.
  v10m_max = 0.
end if

return
end subroutine histave_init

subroutine histave_end

implicit none

deallocate(eg_ave,fg_ave,ga_ave,epan_ave,dew_ave)
deallocate(cbas_ave,ctop_ave,rndmax,prhmax,prhour)
deallocate(tmaxscr,tminscr,tscr_ave)
deallocate(rhmaxscr,rhminscr,rhscr_ave)
!deallocate(riwp_ave,rlwp_ave)
deallocate(u10max,v10max)
deallocate(u1max,v1max,u2max,v2max,cape_max,cape_ave,epot_ave)
deallocate(rnet_ave)
deallocate(wb_ave,wbice_ave,convh_ave)
deallocate(anthropogenic_ave, tmaxurban, tminurban)
deallocate(anth_elecgas_ave, anth_heating_ave, anth_cooling_ave)
deallocate(hailradave_ave, hailradmax_max)
!deallocate(tgg_ave)

if ( allocated(fpn_ave) ) then
  deallocate(fnee_ave)  
  deallocate(fpn_ave,frd_ave)
  deallocate(frp_ave,frpw_ave,frpr_ave,frs_ave)
  deallocate(cnpp_ave,cnbp_ave)
  if ( diaglevel_carbon > 0 ) then
    deallocate(fevc_ave)
  end if
end if

if ( allocated(u_max) ) then
  deallocate(u_max,v_max)
  deallocate(u10m_max,v10m_max)
end if

return
end subroutine histave_end

end module histave_m
