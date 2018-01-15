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
    
subroutine gettin(n)

!  gettina reads back preceding u,v,t,psl for n=2 (for nmi, assuming mex=1)
!  saves and re-reads initial arrays of t and psl
!  called only by darlam and vmodes
    
use arrays_m
use newmpar_m
use savuvt_m

implicit none
      
integer, intent(in) :: n
      
if ( n==0 ) then
  savt(:,:)=t(1:ifull,:)
  savpsl(:)=psl(1:ifull)
else if ( n==2 ) then
  t(1:ifull,:)=savt(:,:)   ! for n=1, n=2
  psl(1:ifull)=savpsl(:)   ! for n=1, n=2
  u(1:ifull,:)=savu(:,:)   ! only for n=2 (VMI init)
  v(1:ifull,:)=savv(:,:)
else    ! for n=1
  t(1:ifull,:)=savt(:,:)   ! for n=1, n=2
  psl(1:ifull)=savpsl(:)   ! for n=1, n=2
endif

return
end subroutine gettin
