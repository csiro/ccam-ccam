! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public precip, precc, rnd_3hr, cape, evspsbl_ave, sbl_ave
public sno, grpl
public cape_d, cin_d, li_d
public prec_init, prec_end

real, dimension(:), allocatable, save :: cape
real, dimension(:), allocatable, save :: cape_d, cin_d, li_d
real(kind=8), dimension(:), allocatable, save :: evspsbl_ave, sbl_ave
real(kind=8), dimension(:), allocatable, save :: precc, precip
real(kind=8), dimension(:), allocatable, save :: sno, grpl
real(kind=8), dimension(:,:), allocatable, save :: rnd_3hr

contains

subroutine prec_init(ifull)

implicit none

integer, intent(in) :: ifull

allocate(precip(ifull),precc(ifull),rnd_3hr(ifull,8),cape(ifull))
allocate(sno(ifull),grpl(ifull))
allocate(evspsbl_ave(ifull),sbl_ave(ifull))
allocate(cape_d(ifull),cin_d(ifull),li_d(ifull))


! needs to be initialised here for zeroth time-step in outcdf.f90
evspsbl_ave(:) = 0._8
sbl_ave(:)   = 0._8
precip(:)    = 0._8
precc(:)     = 0._8
rnd_3hr(:,:) = 0._8
sno(:)       = 0._8
grpl(:)      = 0._8
cape(:)      = 0.
cape_d(:)    = 0.
cin_d(:)     = 0.
li_d(:)      = 0.

return
end subroutine prec_init

subroutine prec_end

implicit none

deallocate(precip,precc,rnd_3hr,cape)
deallocate(sno,grpl)
deallocate(evspsbl_ave,sbl_ave)
deallocate(cape_d,cin_d,li_d)

return
end subroutine prec_end

end module prec_m