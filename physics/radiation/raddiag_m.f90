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
    
module raddiag_m

implicit none

private
public odcalc
public sint_ave,sot_ave,soc_ave,sgn_ave
public sgdn_ave,rgdn_ave,sgdn,rgdn,rgn
public sgdndir_ave,sgdndir
public rtu_ave,rtc_ave,rgn_ave,rgc_ave,sgc_ave
public rgdc_ave,sgdc_ave
public cld_ave,cll_ave,clm_ave,clh_ave,dni_ave
public sunhours,sint,sout,rt,dni
public soutclr,rtclr,rgclr,sgclr
public rgdclr, sgdclr
public sgdn_amp, dni_amp, sw_tend_amp, sint_amp, sout_amp
public soutclr_amp, sgclr_amp, sgn_amp, sgdndir_amp
public sw_tend, lw_tend, sgdclr_amp
public raddiag_init,raddiag_end

real, dimension(:), allocatable, save :: sgdn, sgdndir, sint, sout, rt, dni
real, dimension(:), allocatable, save :: rgdn, rgn
real, dimension(:), allocatable, save :: soutclr, sgclr, sgdclr
real, dimension(:), allocatable, save :: rtclr, rgclr, rgdclr
real, dimension(:), allocatable, save :: sgdn_amp, dni_amp, sint_amp, sout_amp
real, dimension(:), allocatable, save :: soutclr_amp, sgclr_amp, sgn_amp, sgdndir_amp
real, dimension(:), allocatable, save :: sgdclr_amp
real, dimension(:,:), allocatable, save :: sw_tend, lw_tend
real, dimension(:,:), allocatable, save :: sw_tend_amp
real(kind=8), dimension(:), allocatable, save :: sgdn_ave, sint_ave, sot_ave, soc_ave, sgdndir_ave
real(kind=8), dimension(:), allocatable, save :: sgn_ave, sgc_ave, sgdc_ave
real(kind=8), dimension(:), allocatable, save :: rgdn_ave, rtu_ave, rtc_ave, rgn_ave, rgc_ave, rgdc_ave
real(kind=8), dimension(:), allocatable, save :: cld_ave, cll_ave, clm_ave, clh_ave
real(kind=8), dimension(:), allocatable, save :: dni_ave, sunhours
logical, save :: odcalc = .false.

contains

subroutine raddiag_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(sint_ave(ifull),sot_ave(ifull),soc_ave(ifull),sgn_ave(ifull))
allocate(sgdn_ave(ifull),rgdn_ave(ifull),sgdn(ifull),rgdn(ifull),rgn(ifull))
allocate(sgdndir_ave(ifull),sgdndir(ifull))
allocate(rtu_ave(ifull),rtc_ave(ifull),rgn_ave(ifull),rgc_ave(ifull),sgc_ave(ifull))
allocate(rgdc_ave(ifull),sgdc_ave(ifull))
allocate(cld_ave(ifull),cll_ave(ifull),clm_ave(ifull),clh_ave(ifull),dni_ave(ifull))
allocate(sunhours(ifull),sint(ifull),sout(ifull),rt(ifull),dni(ifull))
allocate(soutclr(ifull),rtclr(ifull),rgclr(ifull),sgclr(ifull))
allocate(rgdclr(ifull),sgdclr(ifull))
allocate(sgdn_amp(ifull),dni_amp(ifull),sint_amp(ifull),sout_amp(ifull))
allocate(soutclr_amp(ifull),sgclr_amp(ifull),sgn_amp(ifull),sgdndir_amp(ifull))
allocate(sgdclr_amp(ifull))
allocate(sw_tend(ifull,kl),lw_tend(ifull,kl))
allocate(sw_tend_amp(ifull,kl))

! needs to be initialised here for zeroth time-step in outcdf.f90
sint_ave=0._8
sot_ave=0._8
soc_ave=0._8
sgn_ave=0._8
sgdn_ave=0._8
rgdn_ave=0._8
sgdn=0.
rgdn=0.
rgn=0.
sgdndir_ave=0._8
sgdndir=0.
rtu_ave=0._8
rtc_ave=0._8
rgn_ave=0._8
rgc_ave=0._8
rgdc_ave=0._8
sgc_ave=0._8
sgdc_ave=0._8
cld_ave=0._8
cll_ave=0._8
clm_ave=0._8
clh_ave=0._8
dni_ave=0._8
sunhours=0._8
sint=0.
sout=0.
rt=0.
dni=0.
soutclr=0.
rtclr=0.
rgclr=0.
sgclr=0.
rgdclr=0.
sgdclr=0.
sgdn_amp=0.
dni_amp=0.
sint_amp=0.
sout_amp=0.
soutclr_amp=0.
sgclr_amp=0.
sgn_amp=0.
sgdndir_amp=0.
sgdclr_amp=0.
sw_tend=0.
lw_tend=0.
sw_tend_amp=0.

return
end subroutine raddiag_init

subroutine raddiag_end

implicit none

deallocate(sint_ave,sot_ave,soc_ave,sgn_ave)
deallocate(sgdn_ave,rgdn_ave,sgdn,rgdn,rgn)
deallocate(sgdndir_ave,sgdndir)
deallocate(rtu_ave,rtc_ave,rgn_ave,rgc_ave,sgc_ave)
deallocate(rgdc_ave,sgdc_ave)
deallocate(cld_ave,cll_ave,clm_ave,clh_ave,dni_ave)
deallocate(sunhours,sint,sout,rt,dni)
deallocate(soutclr,rtclr,rgclr,sgclr)
deallocate(rgdclr,sgdclr)
deallocate(sgdn_amp,dni_amp,sint_amp,sout_amp)
deallocate(soutclr_amp,sgclr_amp,sgn_amp,sgdndir_amp)
deallocate(sgdclr_amp)
deallocate(sw_tend,lw_tend)
deallocate(sw_tend_amp)

return
end subroutine raddiag_end

end module raddiag_m