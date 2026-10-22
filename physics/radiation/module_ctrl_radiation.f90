! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module module_ctrl_radiation

use module_aux_rad                         ! Additional cloud and radiation routines
use ozoneread                              ! Ozone input routines
use raddiag_m                              ! Radiation diagnostic
use seaesfrad_m                            ! SEA-ESF radiation

implicit none

private

public ctrl_radiation_init, radiation_arrays_init
public radiation_calc

! module_aux_rad
public liqradmethod, iceradmethod
public cloud3

! ozoneread
public o3_read, o3set, fieldinterpolate, o3regrid
public o3read_amip, o3set_amip
public o3_vert_interpolate

! raddiag
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

! seaesf
public seaesfrad_settime, seaesfrad, sw_resolution, sw_diff_streams
public carbonradmethod, so4radmethod, dustradmethod, seasaltradmethod, lwem_form
public csolar, linecatalog_form, continuum_form, do_co2_10um
public seaice_albvis, seaice_albnir
public do_quench, remain_rayleigh_bug

contains
    
subroutine ctrl_radiation_init(nrad)

implicit none

integer, intent(in) :: nrad

select case(nrad)
  case(5)
    call seaesfrad_init
end select

return
end subroutine ctrl_radiation_init

subroutine radiation_arrays_init(ifull,kl)

implicit none

integer, intent(in) :: ifull, kl

call raddiag_init(ifull,kl)

return
end subroutine radiation_arrays_init

subroutine radiation_arrays_end

implicit none

call raddiag_end

return
end subroutine radiation_arrays_end

subroutine radiation_calc(mtimer_sav,mins_gmt,koundiag)

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use cc_mpi                                 ! CC MPI routines
use cfrac_m                                ! Cloud fraction
use dates_m                                ! Date data
use newmpar_m                              ! Grid parameters
use parm_m                                 ! Model configuration
use pbl_m                                  ! Boundary layer arrays

implicit none

integer, intent(inout) :: mtimer_sav, mins_gmt, koundiag
integer tile, js, je, k

select case ( nrad )
  case(4)
    ! Fels-Schwarzkopf radiation
    !$omp barrier  
    !$omp single          
    if ( nhstest<0 ) then    ! aquaplanet test -1 to -8  
      mtimer_sav = mtimer
      mtimer     = mins_gmt  ! so radn scheme repeatedly works thru same day
      call radrive(il*nrows_rad)
      mtimer = mtimer_sav
    else
      call radrive(il*nrows_rad)  
    end if    ! (nhstest<0)
    !$omp end single
  case(5)
    ! GFDL SEA-EFS radiation
    call seaesfrad(koundiag)
  case DEFAULT
    ! use preset slwa array (use +ve nrad)
    !$omp do schedule(static) private(js,je)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      slwa(js:je) = -real(10*nrad)
    end do
    !$omp end do nowait
end select
!$omp do schedule(static) private(js,je,k)  
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  do k = 1,kl
    t(js:je,k) = t(js:je,k) - dt*(sw_tend(js:je,k)+lw_tend(js:je,k))
  end do
end do
!$omp end do nowait

return
end subroutine radiation_calc

end module module_ctrl_radiation