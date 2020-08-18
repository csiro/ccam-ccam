! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

module prec_m

implicit none

private
public precip,precc,rnd_3hr,cape,evspsbl,sbl
public prec_init,prec_end
!public evap

real, dimension(:), allocatable, save :: evspsbl, sbl
real, dimension(:), allocatable, save :: cape, precc, precip
real, dimension(:,:), allocatable, save :: rnd_3hr
!real, dimension(:), allocatable, save :: evap

contains

subroutine prec_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(precip(ifull),precc(ifull),rnd_3hr(ifull,8),cape(ifull))
allocate(evspsbl(ifull),sbl(ifull))
!allocate(evap(ifull))

! needs to be initialised here for zeroth time-step in outcdf.f90
!evap(:)     = 0.
evspsbl(:)   = 0.
sbl(:)       = 0.
precip(:)    = 0.
precc(:)     = 0.
rnd_3hr(:,:) = 0.
cape(:)      = 0.

return
end subroutine prec_init

subroutine prec_end

implicit none

deallocate(precip,precc,rnd_3hr,cape)
deallocate(evspsbl,sbl)
!deallocate(evap)

return
end subroutine prec_end

end module prec_m