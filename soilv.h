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


      real, dimension(0:mxst) :: swilt, ssat, sfc
      real, dimension(mxst)   :: bch, cnsd, css, hsbh, hyds, rhos, sucs,  &
     &                           clay, sand, silt
      real, dimension(mxvt)   :: rs20
      integer, dimension(mxst) :: i2bp3, ibp2
      real, dimension(44) :: rlaim44,rlais44,scveg44,rsunc44,slveg44
      real :: froot(5), zse(ms)
      common/soilpr/swilt,ssat,sfc,bch,cnsd,css,hsbh,hyds,i2bp3,ibp2,     &
     &              rhos,sucs,clay,sand,silt,rlaim44,rlais44,scveg44,     &
     &              rsunc44,slveg44,froot,zse,rs20 ! MJT cable

