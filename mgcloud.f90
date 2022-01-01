! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
  
! This module is the Morrison-Gettelman cloud microphysics parameterisation for CCAM
    
module mgcloud_m

private

public mg_progcld, mg_2cond

contains

!       The statistical cloud scheme treatment, which is used as
!       a replacement for the Tiedtke cloud fraction scheme, is based
!       on a number of publications: Tompkins, A., 2002: J. Atmos. 
!       Sci., 59, 1917-1942, Klein et al., 2005: J. Geophys. Res., 
!       110, D15S06, doi:10.1029/2004JD005017. 
    
subroutine mg_progcld

write(6,*) "ERROR: mg_progcld is not supported"
stop

return
end subroutine mg_progcld

subroutine mg_2cond

write(6,*) "ERROR: mg_2cond is not supported"
stop

return
end subroutine mg_2cond

end module mgcloud_m