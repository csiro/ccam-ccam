! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module cloudmod

! prognostic cloud fraction scheme based on Tiedtke from GFDL-CM3.

! This module works with the LDR cloud microphysics scheme with
! prognostic cloud liquid, frozen/snow and rain condensate.
    
implicit none
    
private
public convectivecloudfrac, convectivecloudarea

contains


subroutine convectivecloudfrac(clcon,kbsav,ktsav,condc,cldcon)

use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer k, kl
real, dimension(:,:), intent(out) :: clcon
real, dimension(:), intent(out), optional :: cldcon
real, dimension(size(clcon,1)) :: cldcon_temp
real, dimension(size(clcon,1)) :: n, cldcon_local
integer, dimension(:), intent(in) :: kbsav
integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc

! MJT notes - This is an old parameterisation from NCAR.  acon and
! bcon represent shallow and deep convection, respectively.  It can
! be argued that acon should be zero in CCAM to avoid discontinuous
! evolution of the model.  Furthermore, a much better fit can be
! obtained from the mass flux, rather than rainfall.  It also
! should be noted that acon and bcon are likely to depend on the
! spatial resolution.

kl = size(clcon,2)

cldcon_temp = 0. ! for cray compiler
call convectivecloudarea(cldcon_temp,ktsav,condc)
if ( present(cldcon) ) then
  cldcon = cldcon_temp
end if

! Impose cloud overlap assumption
if ( nmr>=1 ) then
  cldcon_local = cldcon_temp               ! maximum overlap
else
  n = 1./real(max(ktsav-kbsav,1))
  cldcon_local = 1. - (1.-cldcon_temp)**n  ! random overlap
end if

do k = 1,kl
  where( k<kbsav+1 .or. k>ktsav )
    clcon(:,k) = 0.
  elsewhere
    clcon(:,k) = cldcon_local
  end where
end do

return
end subroutine convectivecloudfrac

pure subroutine convectivecloudarea(cldcon,ktsav,condc)

use newmpar_m        ! Grid parameters
use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc
real, dimension(:), intent(out) :: cldcon

where ( ktsav<kl-1 )
  cldcon = min( acon+bcon*log(1.+condc*86400./dt), 0.8 ) !NCAR
elsewhere
  cldcon = 0.
end where

return
end subroutine convectivecloudarea

end module cloudmod
