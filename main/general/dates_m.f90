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
public calendar_function

integer, save :: ktime, kdate
integer, save :: mtimer = 0
real, save :: timer = 0.
real, save :: timeg = 0.

contains

subroutine calendar_function(mdays,kdate,leap)

implicit none

integer, dimension(1:12), intent(out) :: mdays
integer, intent(in) :: kdate, leap
integer iyr, month

iyr = kdate/10000
month = (kdate-10000*iyr)/100
if ( leap==0 ) then ! 365 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
else if ( leap==1 ) then ! 365/366 day calendar
  mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  if (mod(iyr,4)==0) mdays(2)=29
  if (mod(iyr,100)==0) mdays(2)=28
  if (mod(iyr,400)==0) mdays(2)=29
else if ( leap==2 ) then ! 360 day calendar
  mdays=(/30,30,30,30,30,30,30,30,30,30,30,30/)
else
  write(6,*) "ERROR: Unknown option for leap = ",leap
  stop -1
end if

return
end subroutine calendar_function

end module dates_m