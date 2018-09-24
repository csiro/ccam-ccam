!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory cTYPE(casa_flux_type), INTENT(IN) :: casaflux ! casa fluxesontaining CABLE code.
!
! ==============================================================================
! Purpose: Defines input/output related variables for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Development by Gab Abramowitz
!          Additional code to use multiple vegetation types per grid-cell (patches)
!
! ==============================================================================
MODULE cable_IO_vars_module

   IMPLICIT NONE

   PUBLIC

   INTEGER :: wlogn

END MODULE cable_IO_vars_module
