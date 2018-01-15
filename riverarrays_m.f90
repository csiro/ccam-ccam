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
    
module riverarrays_m

implicit none

private
public watbdy, outflowmask
public river_vel, river_outdir, river_dx
public river_discharge
public riverarrays_init, riverarrays_end, rivervector

integer, dimension(:), allocatable, save :: river_outdir
real, dimension(:), allocatable, save :: watbdy
real, dimension(:), allocatable, save :: river_vel, river_dx
real, dimension(:), allocatable, save :: river_discharge
logical, dimension(:), allocatable, save :: outflowmask

contains

subroutine riverarrays_init(ifull,iextra,nriver)

implicit none

integer, intent(in) :: ifull, iextra, nriver

if ( abs(nriver)>0 ) then
  allocate( watbdy(ifull+iextra) )
  allocate( outflowmask(ifull) )
  allocate( river_vel(ifull), river_dx(ifull) )
  allocate( river_discharge(ifull) )
  allocate( river_outdir(ifull) )
  watbdy(1:ifull+iextra) = 0.
  outflowmask(1:ifull) = .false.
  river_vel(1:ifull) = 0.
  river_dx(1:ifull) = 1.e-9
  river_discharge(1:ifull) = 0.
  river_outdir(1:ifull) = -1
end if
  
return
end subroutine riverarrays_init

subroutine riverarrays_end

implicit none

if ( allocated( watbdy ) ) then
  deallocate( watbdy )
  deallocate( outflowmask )
  deallocate( river_vel, river_dx )
  deallocate( river_discharge )
  deallocate ( river_outdir )
end if

return
end subroutine riverarrays_end

subroutine rivervector(xvec,yvec)

use newmpar_m

implicit none

real, dimension(ifull), intent(out) :: xvec, yvec

xvec(:) = 0.
yvec(:) = 0.

where ( river_outdir(:)==1 )     ! N
  xvec(:) = 0.
  yvec(:) = river_vel(:)
elsewhere ( river_outdir(:)==2 ) ! E
  xvec(:) = river_vel(:)
  yvec(:) = 0.
elsewhere ( river_outdir(:)==3 ) ! S
  xvec(:) = 0.
  yvec(:) = -river_vel(:)
elsewhere ( river_outdir(:)==4 ) ! W
  xvec(:) = -river_vel(:)
  yvec(:) = 0.
elsewhere ( river_outdir(:)==5 ) ! NE
  xvec(:) = sqrt(0.5)*river_vel(:)
  yvec(:) = sqrt(0.5)*river_vel(:)
elsewhere ( river_outdir(:)==6 ) ! SE
  xvec(:) = sqrt(0.5)*river_vel(:)
  yvec(:) = -sqrt(0.5)*river_vel(:)
elsewhere ( river_outdir(:)==7 ) ! SW
  xvec(:) = -sqrt(0.5)*river_vel(:)
  yvec(:) = -sqrt(0.5)*river_vel(:)
elsewhere ( river_outdir(:)==8 ) ! NW
  xvec(:) = -sqrt(0.5)*river_vel(:)
  yvec(:) = sqrt(0.5)*river_vel(:)
end where

return
end subroutine rivervector

end module riverarrays_m