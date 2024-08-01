! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module module_ctrl_convection

implicit none

private
public ctrl_convection

contains
    
subroutine ctrl_convection

use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi                        ! CC MPI routines
use convjlm_m                     ! Convection
use convjlm22_m                   ! Convection v2
use kuocom_m                      ! JLM convection
use nlin_m                        ! Atmosphere non-linear dynamics
use soil_m                        ! Soil and surface data

implicit none

select case ( interp_convection(nkuo) )
  case("betts_conv")
    !$omp barrier
    !$omp single
    call betts(t,qg,tn,land,ps) ! not called these days
    !$omp end single
  case("john_conv22")
    call convjlm22              ! split convjlm
  case("john_conv")
    call convjlm                ! split convjlm
  case default
    write(6,*) "ERROR: unknown convection option nkuo=",nkuo
    call ccmpi_abort(-1) 
end select

return
end subroutine ctrl_convection



!====================================================================================================
! SUBROUTINE interp_nconvection
!   
! subroutine to select the cloud convection scheme for CCAM
!====================================================================================================
pure function interp_convection(nconvection) result(cv_physics)

implicit none

integer, intent(in) :: nconvection
character(len=20) :: cv_physics

cv_physics = "ERROR"

select case(nconvection)
  case(5)
    cv_physics = "betts_conv"
  case(21,22)
    cv_physics = "john_conv22"
  case(23,24)
    cv_physics = "john_conv"
  case(31)
    cv_physics = "grell_conv"
end select

return
end function interp_convection  
  
end module module_ctrl_convection
