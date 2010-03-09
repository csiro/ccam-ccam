! cable_cbm.f90
!
! Source file containing main routine and canopy code for CABLE, 
! CSIRO land surface model
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 coding by Harvey Davies, Gab Abramowitz and Martin Dix
! bugs to gabsun@gmail.com.
!
! This file contains modules:
! cbm_module, air_module, roughness_module, radiation_module, 
! and canopy_module.
!
! Most user-defined types (e.g. met%tk) are defined in define_types module
! in cable_variables.f90

!=========================================================================
MODULE air_module
  USE physical_constants
  USE define_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_air
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE define_air
!  SUBROUTINE define_air(met,air)
!    TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
!    TYPE (met_type), INTENT(IN)	 :: met ! meteorological variables
    REAL(r_1), DIMENSION(mp)	 :: es	! sat vapour pressure (mb)
!    print *,'air 1',met%tvair,tfrz
    es	 = tetena * EXP(tetenb * (met%tvair-tfrz)/(tetenc + (met%tvair-tfrz)))
    ! Calculate conversion factor from from m/s to mol/m2/s
!    print *,'air 2',es
    air%cmolar = met%pmb * 100.0 / (rgas * (met%tvair))
!    print *,'air 3',met%pmb,rgas,air%cmolar
    ! Calculate dry air density:
    air%rho = min(1.3,rmair * air%cmolar)
!    print *,'air 4',air%rho
    ! molar volume (m^3/mol)
    air%volm = rgas * (met%tvair) / (100.0 * met%pmb)
!    print *,'air 5',air%volm
    ! latent heat for water (j/kg)
    air%rlam = (2501.0 - 2.38 * (met%tvair- tfrz)) * 1000.0
!    print *,'air 6',air%rlam
    air%rlam=2.5104e6
    ! saturation specific humidity
    air%qsat = (rmh2o / rmair) * es / met%pmb
!    print *,'air 7',air%qsat
    ! d(qsat)/dT ((kg/kg)/K)
    air%epsi = (air%rlam / capp) * (rmh2o / rmair) * es * tetenb * tetenc / &
         (tetenc + (met%tvair - tfrz)) ** 2 / met%pmb
!    print *,'air 8',air%epsi
    ! air kinematic viscosity (m^2/s)
    air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - tfrz))
!    print *,'air 9',air%visc
    ! psychrometric constant
    air%psyc = met%pmb * 100.0 * capp * rmair / air%rlam / rmh2o
!    print *,'air 10',air%psyc
    ! d(es)/dT (mb/K)

!sxy
    ! d(es)/dT (mb/K)
    air%dsatdk = 100.0 * (tetena*tetenb*tetenc) / ((met%tvair-tfrz)+tetenc)**2 &
               * EXP(tetenb * (met%tvair-tfrz) / ((met%tvair-tfrz) + tetenc))
!    air%dsatdk = (610.078 * 17.27 * 237.3) / ((met%tvair-tfrz)+237.2)** 2 * &
!         EXP(17.27 * (met%tvair-tfrz) / ((met%tvair-tfrz) + 237.3))
!    print *,'air 11',air%dsatdk
  END SUBROUTINE define_air
END MODULE air_module
