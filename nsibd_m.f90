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
    
module nsibd_m

implicit none

private
public rsmin,sigmf,tgf,sigmu
public ivegt,isoilm,isoilm_in
public nsibd_init,nsibd_end

real, dimension(:), allocatable, save :: rsmin,sigmf,tgf,sigmu
integer, dimension(:), allocatable, save :: ivegt,isoilm,isoilm_in

contains

subroutine nsibd_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib


allocate(ivegt(ifull),isoilm(ifull))
allocate(sigmf(ifull),sigmu(ifull))
allocate(rsmin(ifull),isoilm_in(ifull))
sigmf=0.
sigmu=0.
rsmin=995.
if (nsib==3.or.nsib==5) then
  allocate(tgf(ifull))
  tgf=293.
end if

return
end subroutine nsibd_init

subroutine nsibd_end

implicit none

deallocate(ivegt,isoilm)
deallocate(sigmf,sigmu)
deallocate(rsmin,isoilm_in)
if (allocated(tgf)) then
  deallocate(tgf)
end if

return
end subroutine nsibd_end

end module nsibd_m