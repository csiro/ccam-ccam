! Conformal Cubic Atmospheric Model
    
! Copyright 2016-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module dates_m

implicit none

private
public ktime, kdate, mtimer
public timer, timeg
public cal_365, cal_leap, cal_360
public calendar_function, getzinp
public time_of_month, time_interpolate

integer, save :: ktime, kdate
integer, save :: mtimer = 0
real, save :: timer = 0.
real, save :: timeg = 0.

! calendar
integer, parameter :: cal_365 = 0
integer, parameter :: cal_leap = 1
integer, parameter :: cal_360 = 2

contains

!--------------------------------------------------------------------
! Set up number of days per month
subroutine calendar_function(mdays,kdate,leap)

integer, dimension(1:12), intent(out) :: mdays
integer, intent(in) :: kdate, leap
integer iyr, month

iyr = kdate/10000
month = (kdate-10000*iyr)/100
if ( leap==cal_365 ) then ! 365 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
else if ( leap==cal_leap ) then ! 365/366 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  if (mod(iyr,4)==0) mdays(2)=29
  if (mod(iyr,100)==0) mdays(2)=28
  if (mod(iyr,400)==0) mdays(2)=29
else if ( leap==cal_360 ) then ! 360 day calendar
  mdays=(/30,30,30,30,30,30,30,30,30,30,30,30/)
else
  write(6,*) "ERROR: Unknown option for leap = ",leap
  stop -1
end if

return
end subroutine calendar_function

!--------------------------------------------------------------------
! Set up number of minutes from beginning of year
subroutine getzinp(jyear,jmonth,jday,jhour,jmin,mins,allleap)

use parm_m

implicit none

integer, intent(out) :: jyear,jmonth,jday,jhour,jmin ! start date of run
integer, intent(out) :: mins                         ! elapsed time from start of year
integer mstart, n1
integer, dimension(12) :: ndoy
logical, intent(in), optional :: allleap ! force use of leap days even if leap=0
logical lleap

if ( present(allleap) ) then
  lleap = allleap
else
  lleap = .false.    
end if

jyear  = kdate/10000
jmonth = (kdate-jyear*10000)/100
jday   = kdate - jyear*10000 - jmonth*100
jhour  = ktime/100
jmin   = ktime - jhour*100
      
if ( jmonth<1 .or. jmonth>12 ) then
  write(6,*) "ERROR: Invalid month ",jmonth," for kdate ",kdate
  stop
end if 

if ( lleap .or. leap==cal_leap ) then ! 365/366 day calendar
  ndoy(:) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  n1 = 0
  if ( mod(jyear,4)  ==0 ) n1 = 1
  if ( mod(jyear,100)==0 ) n1 = 0
  if ( mod(jyear,400)==0 ) n1 = 1
  ndoy(3:12)=ndoy(3:12)+n1
else if ( leap==cal_365 ) then ! 365 day calendar
  ndoy(:) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)  
else if ( leap==cal_360 ) then ! 360 day calendar  
  ndoy(:) = (/ 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 /)
else
  write(6,*) "ERROR: Unknown option for leap = ",leap
  stop
end if  

mstart = 1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year

! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

return
end subroutine getzinp

!--------------------------------------------------------------------
subroutine time_of_month(iyr,imo,iday,mdays,kdate,mtimer,leap)

implicit none

integer, intent(in) :: kdate, leap, mtimer
integer, intent(out) :: iyr, imo, iday, mdays
integer, dimension(0:13) :: mdays_a

iyr = kdate/10000
imo = (kdate-10000*iyr)/100
iday = kdate - 10000*iyr - 100*imo + mtimer/(60*24)
if ( leap==0 ) then ! 365 day calendar
  mdays_a = (/ 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 /)
else if ( leap==1 ) then ! 365/366 day calendar
  mdays_a = (/ 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 /)
  if ( mod(iyr,4)==0 ) mdays_a(2) = 29
  if ( mod(iyr,100)==0 ) mdays_a(2) = 28
  if ( mod(iyr,400)==0 ) mdays_a(2) = 29
else if ( leap==2 ) then ! 360 day calendar
  mdays_a = (/ 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 /)  
else
  write(6,*) "ERROR: Invalid leap = ",leap
  stop
end if
do while ( iday>mdays_a(imo) )
  iday = iday - mdays_a(imo)
  imo = imo + 1
  if ( imo>12 ) then
    imo = 1
    iyr = iyr + 1
    if ( leap==1 ) then
      if ( mod(iyr,4)==0 ) mdays_a(2) = 29
      if ( mod(iyr,100)==0 ) mdays_a(2) = 28
      if ( mod(iyr,400)==0 ) mdays_a(2) = 29
    end if
  end if
end do
mdays = mdays_a(imo)

return
end subroutine time_of_month

!--------------------------------------------------------------------
subroutine time_interpolate(datout,dat_pm,dat_p,dat_c,dat_n,dat_np,x,method,land)

integer, intent(in) :: method
integer iq
real, intent(in) :: x
real a0, a1, a2, aa, bb, cc, mp1, mp2, c2, c3, c4
real ssta2, ssta3, ssta4
real, dimension(:), intent(out) :: datout
real, dimension(:), intent(in) :: dat_pm, dat_p, dat_c, dat_n, dat_np
logical, dimension(:), intent(in), optional :: land
logical, dimension(size(datout)) :: land_l

! dat_pm - previous month -1 (current month - 2)
! dat_p  - previous month (current month - 1)
! dat_c  - current month
! dat_n  - next month ( current month + 1)
! dat_np - next month +1 ( current month + 2 )

land_l = .false.
if ( present(land) ) then
  land_l = land
end if

! method=3,4,5  Use PWCB interpolation
if ( method==3 .or. method==4 .or. method==5 ) then
  !--------------------------------------------------------------------------------------------------
  ! Piece-wise cubic bessel interpolation
  do iq=1,size(datout)
    if( .not.land_l(iq) )then
      c2=dat_p(iq)
      c3=dat_p(iq)+dat_c(iq)
      c4=c3+dat_n(iq)          
      datout(iq)=.5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x
    endif      ! (.not.land(iq))
  enddo
endif

! method=11,13,14,15  Use JMc interpolation
if ( method==11 .or. method==13 .or. method==14 .or. method==15 ) then
  !--------------------------------------------------------------------------------------------------
  ! John McGregor 5-pt, piece-wise, cubic interpolation
  do iq=1,size(datout)
    if ( .not.land_l(iq) ) then
      ssta2=(24.*dat_p(iq)-dat_pm(iq)-dat_c(iq))/22. 
      ssta3=(24.*dat_c(iq)-dat_p(iq)-dat_n(iq))/22. 
      ssta4=(24.*dat_n(iq)-dat_c(iq)-dat_np(iq))/22. 
      c2=(-dat_pm(iq)+9*ssta2+9*ssta3-dat_n(iq))/16.
      c3=(-dat_p(iq)+9*ssta3+9*ssta4-dat_np(iq))/16.
      datout(iq)=c2+(6.*dat_c(iq)-4.*c2-2.*c3)*x    &
                 +(3.*c2+3.*c3-6.*dat_c(iq))*x*x 
    endif      ! (.not.land(iq))
  enddo
end if

! method=21,24,25  Use approx linear AMIP interpolation
if ( method==21 .or. method==24 .or. method==25 ) then
  !--------------------------------------------------------------------------------------------------
  ! Approximation of piece-wise, linear AMIP interpolation
  if ( x<0.5 ) then
    do iq = 1,size(datout)
      if ( .not.land_l(iq) ) then
        a0 = 0.5*dat_p(iq)
        a1 = -dat_c(iq)
        a2 = 0.5*dat_n(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp1 = 0.25*aa + 0.5*bb + cc ! start of month value
        a0 = 0.5*dat_c(iq)
        a1 = -dat_n(iq)
        a2 = 0.5*dat_np(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp2 = 0.25*aa + 0.5*bb + cc ! end of month value
        c4 = 2.*dat_c(iq) - 0.5*mp1 - 0.5*mp2 ! mid-point value
        c2 = mp1                              ! intercept
        c3 = 2.*(c4-c2)                       ! gradient
        datout(iq) = c3*x + c2
      end if
    end do
  else
    do iq = 1,size(datout)
      if ( .not.land_l(iq) ) then
        a0 = 0.5*dat_p(iq)
        a1 = -dat_c(iq)
        a2 = 0.5*dat_n(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp1 = 0.25*aa + 0.5*bb + cc ! start of month value
        a0 = 0.5*dat_c(iq)
        a1 = -dat_n(iq)
        a2 = 0.5*dat_np(iq)
        aa = a0 + a1 + a2
        bb = -3.*a0 - 2.*a1 - a2
        cc = 2.*a0
        mp2 = 0.25*aa + 0.5*bb + cc ! end of month value
        c4 = 2.*dat_c(iq) - 0.5*mp1 - 0.5*mp2 ! mid-point value
        c3 = 2.*(mp2 - c4)                    ! gradient
        c2 = 2.*c4 - mp2                      ! intercept
        datout(iq) = c3*x + c2
      end if
    end do
  end if
end if

return
end subroutine time_interpolate

end module dates_m