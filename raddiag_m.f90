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
    
module raddiag_m

implicit none

private
public odcalc
public sint_ave,sot_ave,soc_ave,sgn_ave
public sgdn_ave,rgdn_ave,sgdn,rgdn,rgn
public rtu_ave,rtc_ave,rgn_ave,rgc_ave,sgc_ave
public cld_ave,cll_ave,clm_ave,clh_ave,dni_ave
public sunhours,sint,sout,rt,dni
public soutclr,rtclr,rgclr,sgclr
public sgdn_amp, dni_amp, sw_tend_amp, sint_amp, sout_amp
public soutclr_amp, sgclr_amp, sgn_amp
public raddiag_init,raddiag_end
public sw_tend, lw_tend

real, dimension(:), allocatable, save :: sint_ave,sot_ave,soc_ave,sgn_ave
real, dimension(:), allocatable, save :: sgdn_ave,rgdn_ave,sgdn,rgdn,rgn
real, dimension(:), allocatable, save :: rtu_ave,rtc_ave,rgn_ave,rgc_ave,sgc_ave
real, dimension(:), allocatable, save :: cld_ave,cll_ave,clm_ave,clh_ave,dni_ave
real, dimension(:), allocatable, save :: sunhours,sint,sout,rt,dni
real, dimension(:), allocatable, save :: soutclr,rtclr,rgclr,sgclr
real, dimension(:), allocatable, save :: sgdn_amp, dni_amp, sint_amp, sout_amp
real, dimension(:), allocatable, save :: soutclr_amp, sgclr_amp, sgn_amp
real, dimension(:,:), allocatable, save :: sw_tend, lw_tend
real, dimension(:,:), allocatable, save :: sw_tend_amp
logical, save :: odcalc = .false.

contains

subroutine raddiag_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(sint_ave(ifull),sot_ave(ifull),soc_ave(ifull),sgn_ave(ifull))
allocate(sgdn_ave(ifull),rgdn_ave(ifull),sgdn(ifull),rgdn(ifull),rgn(ifull))
allocate(rtu_ave(ifull),rtc_ave(ifull),rgn_ave(ifull),rgc_ave(ifull),sgc_ave(ifull))
allocate(cld_ave(ifull),cll_ave(ifull),clm_ave(ifull),clh_ave(ifull),dni_ave(ifull))
allocate(sunhours(ifull),sint(ifull),sout(ifull),rt(ifull),dni(ifull))
allocate(soutclr(ifull),rtclr(ifull),rgclr(ifull),sgclr(ifull))
allocate(sgdn_amp(ifull),dni_amp(ifull),sint_amp(ifull),sout_amp(ifull))
allocate(soutclr_amp(ifull),sgclr_amp(ifull),sgn_amp(ifull))
allocate(sw_tend(ifull,kl),lw_tend(ifull,kl))
allocate(sw_tend_amp(ifull,kl))

! needs to be initialised here for zeroth time-step in outcdf.f90
sint_ave=0.
sot_ave=0.
soc_ave=0.
sgn_ave=0.
sgdn_ave=0.
rgdn_ave=0.
sgdn=0.
rgdn=0.
rgn=0.
rtu_ave=0.
rtc_ave=0.
rgn_ave=0.
rgc_ave=0.
sgc_ave=0.
cld_ave=0.
cll_ave=0.
clm_ave=0.
clh_ave=0.
dni_ave=0.
sunhours=0.
sint=0.
sout=0.
rt=0.
dni=0.
soutclr=0.
rtclr=0.
rgclr=0.
sgclr=0.
sgdn_amp=0.
dni_amp=0.
sint_amp=0.
sout_amp=0.
soutclr_amp=0.
sgclr_amp=0.
sgn_amp=0.
sw_tend=0.
lw_tend=0.
sw_tend_amp=0.

return
end subroutine raddiag_init

subroutine raddiag_end

implicit none

deallocate(sint_ave,sot_ave,soc_ave,sgn_ave)
deallocate(sgdn_ave,rgdn_ave,sgdn,rgdn,rgn)
deallocate(rtu_ave,rtc_ave,rgn_ave,rgc_ave,sgc_ave)
deallocate(cld_ave,cll_ave,clm_ave,clh_ave,dni_ave)
deallocate(sunhours,sint,sout,rt,dni)
deallocate(soutclr,rtclr,rgclr,sgclr)
deallocate(sgdn_amp,dni_amp,sint_amp,sout_amp)
deallocate(soutclr_amp,sgclr_amp,sgn_amp)
deallocate(sw_tend,lw_tend)
deallocate(sw_tend_amp)

return
end subroutine raddiag_end

end module raddiag_m