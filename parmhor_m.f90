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
    
module parmhor_m

implicit none

private
public mh_bs, nt_adv

!     horizontal advection/staggering options

integer, save :: mh_bs=4, nt_adv=7

!     for RMIP1 m_bs was -2; during 2002 it was 2

!     integer, parameter :: m_bs=-2  !  0 for B&S off     usually -2

!                               2 for B&S on (in ints)

!                              -2 on for gases only

!     m_bs is superseded on 23/7/03 by mh_bs

!     integer, parameter :: mh_bs=4  !  5 for B&S off     usually 4

!                               4 for B&S on for qg, gases (in ints)

!                               3 for B&S on for T, qg, gases 

!                               2 for B&S on for u, v, T, qg, gases 

!                               1 for B&S on for psl, u, v, T, qg, gases


end module parmhor_m