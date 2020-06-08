! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module darcdf_m

implicit none

private
public idnc, ncid, ncidtopo, ncidveg, ncidbath
public iarchi, lnctopo, lncveg, lncbath
public lncveg_numpft, lncveg_numsoil

integer, save :: idnc, ncid, ncidtopo, ncidveg, ncidbath
integer, save :: iarchi, lnctopo=0, lncveg=0, lncbath=0  
integer, save :: lncveg_numpft=-1, lncveg_numsoil=-1

end module darcdf_m