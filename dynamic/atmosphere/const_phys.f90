! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

module const_phys
    
private

      real, parameter, public :: stefbo=5.67e-8 !Stefan-Boltzmann constant
      real, parameter, public :: erad=6.37122e6, eradsq=erad*erad !Radius of earth
      real, parameter, public :: cp=1004.64 ! Specific heat of dry air at const P
      real, parameter, public :: rdry=287.04 ! Specific gas const for dry air
      real, parameter, public :: epsil=0.622 ! Ratio molec wt of H2O vapour to dry air
      real, parameter, public :: rvap=461. ! Gas constant for water vapour
      real, parameter, public :: hlf=3.36e5 ! Latent heat of fusion (at 0 deg C)
      real, parameter, public :: hl=2.5104e6 !Latent heat of vaporization (at 0 deg. C)
      real, parameter, public :: hls=hl+hlf !  "      "   " sublimation
      real, parameter, public :: hlcp=hl/cp, hlfcp=hlf/cp, hlscp=hlcp+hlfcp
      real, parameter, public :: grav=9.80616 ! Acceleration of gravity

      real, parameter, public :: cappa=rdry/cp 
      real, parameter, public :: pi=3.141592653589793, tpi=2.0*pi !Good ol' pi

!      real, parameter :: tfrz=273.15
!     CCAM value, Should this be 273.15 or is there a physical difference?
      real, parameter, public :: tfrz=273.1

      real, parameter, public :: cpv=1869.46 ! Specific heat of water vapor at const P
      real, parameter, public :: roncp=cappa ! Just an alias for the CC model
      real, parameter, public :: rearth=erad 
      real, parameter, public :: ars=rvap
      real, parameter, public :: hlars=hl/ars
      real, parameter, public :: stdlapse=6.5e-3

! --- chemical
      real, parameter, public :: fAIR_MolM = 28.965
      real, parameter, public :: fH_MolM   =  1.00794
      real, parameter, public :: fC_MolM   = 12.011
      real, parameter, public :: fO_MolM   = 15.9994
      real, parameter, public :: fS_MolM   = 32.066

      real, parameter, public :: fH2_MolM  = 2*fH_MolM                  ! molecular hydrogen
      real, parameter, public :: fO2_MolM  = 2*fO_MolM                  ! molecular oxygen
      real, parameter, public :: fCO2_MolM = fC_MolM   +   fO2_MolM     ! carbon dioxide
      real, parameter, public :: fSO2_MolM = fS_MolM   +   fO2_MolM     ! sulphur dioxide
      real, parameter, public :: fSO4_MolM = fSO2_MolM +   fO2_MolM     ! sulphate
      real, parameter, public :: fCH4_MolM = fC_MolM   + 2*fH2_MolM     ! methane
      real, parameter, public :: fH2O_MolM = fO_MolM   +   fH2_MolM     ! water

end module const_phys