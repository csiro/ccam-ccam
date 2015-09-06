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

! Imported from CSIRO GCM

      real, parameter :: stefbo=5.67e-8 !Stefan-Boltzmann constant
      real, parameter :: erad=6.37122e6, eradsq=erad*erad !Radius of earth
      real, parameter :: cp=1004.64 ! Specific heat of dry air at const P
      real, parameter :: rdry=287.04 ! Specific gas const for dry air
      real, parameter :: epsil=0.622 ! Ratio molec wt of H2O vapour to dry air
      real, parameter :: rvap=461. ! Gas constant for water vapour
      real, parameter :: hlf=3.36e5 ! Latent heat of fusion (at 0 deg C)
      real, parameter :: hl=2.5104e6 !Latent heat of vaporization (at 0 deg. C)
      real, parameter :: hls=hl+hlf !  "      "   " sublimation
      real, parameter :: hlcp=hl/cp, hlfcp=hlf/cp, hlscp=hlcp+hlfcp
      real, parameter :: grav=9.80616 ! Acceleration of gravity

      real, parameter :: sq2=1.414213562373092 ! Square root of 2
      real, parameter :: cappa=rdry/cp 
      real, parameter :: tomg=2*7.2921233e-5 ! 2*omega
      real, parameter :: pi=3.141592653589793, tpi=2.0*pi !Good ol' pi

      real, parameter :: hcap50=2.095e8 ! Heat capacity of sea water * 50m (J/m**2/K)
      real, parameter :: hdzmlo=100.0   ! Depth of mixed layer ocean points
      real, parameter :: hcap=hcap50/50.0 ! Heat capacity of 1m of sea water (J/m**2/K)

!      real, parameter :: tfrz=273.15
!     CCAM value, Should this be 273.15 or is there a physical difference?
      real, parameter :: tfrz=273.1

      real, parameter :: cpv=1869.46 ! Specific heat of water vapor at const P
      real, parameter :: roncp=cappa ! Just an alias for the CC model
      real, parameter :: rearth=erad 
      real, parameter :: ars=rvap
      real, parameter :: hlars=hl/ars
      real, parameter :: stdlapse=6.5e-3

! --- chemical
      real, parameter :: fAIR_MolM = 28.965
      real, parameter :: fH_MolM   =  1.00794
      real, parameter :: fC_MolM   = 12.011
      real, parameter :: fO_MolM   = 15.9994
      real, parameter :: fS_MolM   = 32.066

      real, parameter :: fH2_MolM  = 2*fH_MolM                  ! molecular hydrogen
      real, parameter :: fO2_MolM  = 2*fO_MolM                  ! molecular oxygen
      real, parameter :: fCO2_MolM = fC_MolM   +   fO2_MolM     ! carbon dioxide
      real, parameter :: fSO2_MolM = fS_MolM   +   fO2_MolM     ! sulphur dioxide
      real, parameter :: fSO4_MolM = fSO2_MolM +   fO2_MolM     ! sulphate
      real, parameter :: fCH4_MolM = fC_MolM   + 2*fH2_MolM     ! methane
      real, parameter :: fH2O_MolM = fO_MolM   +   fH2_MolM     ! water
