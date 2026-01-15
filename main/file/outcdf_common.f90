! Conformal Cubic Atmospheric Model
    
! Copyright 2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Common data and tools for outcdf
    
module outcdf_common_m
    
private
public month, cordex_levels, height_levels
public cordex_level_data, height_level_data
public mslp, cordex_name, bisect

character(len=3), dimension(12), parameter :: month = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
integer, parameter :: cordex_levels = 17
integer, dimension(cordex_levels) :: cordex_level_data = &
    (/ 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10 /)
integer, parameter :: height_levels = 6
integer, dimension(height_levels) :: height_level_data = &
    (/ 50, 100, 150, 200, 250, 300 /)

contains

subroutine cordex_name(lname,stringa,press_level,stringb)

use cc_mpi, only : ccmpi_abort

implicit none

integer, intent(in) :: press_level
character(len=*), intent(out) :: lname
character(len=*), intent(in) :: stringa
character(len=*), intent(in), optional :: stringb

if ( present(stringb) ) then
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4,A)') stringa,press_level,stringb
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3,A)') stringa,press_level,stringb
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2,A)') stringa,press_level,stringb
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1,A)') stringa,press_level,stringb
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
else
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4)') stringa,press_level
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3)') stringa,press_level
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2)') stringa,press_level
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1)') stringa,press_level
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
end if

return
end subroutine cordex_name

subroutine mslp(pmsl,psl,zs,t)

use cc_mpi, only : mydiag, myid
use const_phys
use newmpar_m
use parm_m
use sigs_m

implicit none
! this one will ignore negative zs (i.e. over the ocean)

integer, save :: lev = -1
integer, dimension(1) :: pos
!real c, conr, con
real, dimension(ifull), intent(out) :: pmsl
real, dimension(ifull), intent(in) :: psl, zs
real, dimension(ifull) :: phi1, tsurf, tav,  dlnps
real, dimension(:,:), intent(in) :: t
      
!c = grav/stdlapse
!conr = c/rdry
if ( lev<0 ) then
  pos = minloc(abs(sig-0.9),sig>=0.9)
  lev = pos(1)
  if ( myid==0 .and. nmaxpr==1 ) then
    write(6,*) "Reducing ps to MSLP with lev,sig ",lev,sig(lev) 
  end if
end if
!con = sig(lev)**(rdry/c)/c
      
phi1(:) = t(1:ifull,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
tsurf(:) = t(1:ifull,lev)+phi1(:)*stdlapse/grav
tav(:) = tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
dlnps(:) = zs(1:ifull)/(rdry*tav(:))
pmsl(:) = 1.e5*exp(psl(:)+dlnps(:))
      
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'meth,lev,sig(lev) ',1,lev,sig(lev)
  write(6,*) 'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
end if
      
return
end subroutine mslp

! Find pressure level
pure function bisect(press_target, ps, sig) result(ans)
real, intent(in) :: press_target, ps
real, dimension(:), intent(in) :: sig
integer :: ans
integer a, b, i, kx

kx = size(sig)
a = 1
b = kx
do while ( b-a > 1 )
  i = (a+b)/2
  if ( press_target > ps*sig(i) ) then
    b = i
  else
    a = i
  end if
end do
if ( ps*sig(a)>=press_target .and. ps*sig(b)>=press_target ) then
  ans = b
else
  ans = a
end if
ans = min( ans, kx-1 )

end function bisect

end module outcdf_common_m
