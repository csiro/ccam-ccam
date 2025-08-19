! Conformal Cubic Atmospheric Model

! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)

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

!
! Calculate integrated vertical advection of moisture
! Formula: -u * (dq/dx) - v * (dq/dy) - q * (du/dx + dv/dy)
!
! Input parameters:
! u(nx,ny)        - u-component of wind (m/s)
! v(nx,ny)        - v-component of wind (m/s)
! q(nx,ny)        - specific humidity (kg/kg)
! dx              - grid spacing in x-direction (m)
! dy              - grid spacing in y-direction (m)
! nx, ny          - grid dimensions
!
! Output:
! moisture_adv(nx,ny) - integrated vertical advection of moisture
!

module calculate_mconv_m

private
public calculate_mconv

contains

subroutine calculate_mconv

use arrays_m
use cc_mpi            ! CC MPI routines
use indices_m         ! Grid index arrays
use kuocom_m          ! JLM convection
use map_m             ! Grid map arrays
use newmpar_m         ! Grid parameters
use parm_m            ! Model configuration

implicit none

integer iq, k
real dq_dx, dq_dy, du_dx, dv_dy
real term1, term2, term3
real dx

! qg, u and v must be ifull+iextra for bounds to work
call bounds(qg(:,:))    ! n, s, e, w
call boundsuv(u(:,:),v(:,:)) ! nv, sv, eu, wu

!call staguv(u,v,u_s,v_s)

! Initialize output array
mconv_save = 0.

! Calculate moisture advection for interior points using centered differences
do k = 1,kl
  do iq = 1,ifull
    
    dx = ds/em(iq) ! =dy

    ! using unstagged coordinates

    ! Calculate moisture gradients using centered differences
    dq_dx = (qg(ie(iq),k) - qg(iw(iq),k))/(2.*dx)
    dq_dy = (qg(in(iq),k) - qg(is(iq),k))/(2.*dx)

    ! Calculate wind divergence using centered differences
    du_dx = (u(ieu(iq),k) - u(iwu(iq),k))/(2.*dx)
    dv_dy = (v(inv(iq),k) - v(isv(iq),k))/(2.*dx)

    ! Calculate the three terms
    term1 = -u(iq,k) * dq_dx
    term2 = -v(iq,k) * dq_dy
    term3 = -qg(iq,k) * (du_dx + dv_dy)

    ! Sum all terms
    mconv_save(iq,k) = term1 + term2 + term3
 
  end do
end do

return
end subroutine calculate_mconv

end module calculate_mconv_m


