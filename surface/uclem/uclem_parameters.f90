! UCLEM urban canopy model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the UCLEM urban canopy model
!
! UCLEM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UCLEM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UCLEM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! parameters used by UCLEM urban climate and energy model

module uclem_parameters

implicit none

private
public urbtemp,energytol,resmeth,zohmeth,acmeth,nrefl,                                 &
       scrnmeth,wbrelaxc,wbrelaxr,nfgits,tol,infilalpha,                               &
       zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden,refheight,     &
       zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf,        &
       intairtmeth,intmassmeth,statsmeth,lwintmeth,cvcoeffmeth,infilmeth,acfactor,     &
       ac_heatcap,ac_coolcap,ac_deltat,cyc_prop,cyc_base,cyc_traf,cyc_tran
public nl
public waterden, icelambda, aircp, icecp, grav, vkar, lv, lf, ls, pi, rd, rv, sbconst
public icmax, a_1, b_1, c_1, d_1

! model parameters
integer, save      :: resmeth=3            ! Canyon sensible heat transfer (0=Masson, 2=Kusaka, 3=Harman (fixed width))
integer, save      :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save      :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, 3=COP of 1.0)
integer, save      :: intairtmeth=0        ! Internal air temperature (0=prescribed, 1=aggregated varying, 2=fractional varying)
integer, save      :: intmassmeth=2        ! Internal thermal mass (0=none, 1=one floor, 2=dynamic floor number)
integer, save      :: cvcoeffmeth=1        ! Internal surface convection heat transfer coefficient (0=DOE, 1=ISO6946, 2=fixed)
integer, save      :: statsmeth=1          ! Use statistically based diurnal QF amendments (0=off, 1=on) from Thatcher 2007 
integer, save      :: infilmeth=0          ! Method to calculate infiltration rate (0=constant, 1=EnergyPlus/BLAST, 2=ISO)
integer, save      :: nrefl=3              ! Number of canyon reflections for radiation (default=3)
integer, save      :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save      :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save      :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save      :: lwintmeth=2          ! Interior longwave flux method (0=off, 1=explicit, 2=with forward step facet temps)
integer, parameter :: nl=4                 ! Number of layers (default 4, must be a multiple of 4)
integer, save      :: iqt=314              ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, save      :: nfgits=8             ! Number of iterations for balancing canyon energy budgets (default=6)
real, save         :: tol=0.001            ! Sectant method tolerance for sensible heat flux (default=0.001)
real, save         :: infilalpha=0.5       ! Weighting for dampening factor when calculating internal infiltration
real(kind=8), save :: energytol=0.005_8    ! Tolerance for acceptable energy closure in each timestep
real, save         :: urbtemp=290.         ! reference temperature to improve precision
! physical parameters
real, parameter    :: waterden=1000.       ! water density (kg m^-3)
real, parameter    :: icelambda=2.22       ! conductance of ice (W m^-1 K^-1)
real, parameter    :: aircp=1004.64        ! Heat capacity of dry air (J kg^-1 K^-1)
real, parameter    :: icecp=2100.          ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter    :: grav=9.80616         ! gravity (m s^-2)
real, parameter    :: vkar=0.4             ! von Karman constant
real, parameter    :: lv=2.501e6           ! Latent heat of vaporisation (J kg^-1)
real, parameter    :: lf=3.337e5           ! Latent heat of fusion (J kg^-1)
real, parameter    :: ls=lv+lf             ! Latent heat of sublimation (J kg^-1)
real, parameter    :: pi=3.14159265        ! pi (must be rounded down for shortwave)
real, parameter    :: rd=287.04            ! Gas constant for dry air
real, parameter    :: rv=461.5             ! Gas constant for water vapour
real, parameter    :: sbconst=5.67e-8      ! Stefan-Boltzmann constant
! snow parameters
real, save         :: zosnow=0.001         ! Roughness length for snow (m)
real, save         :: snowemiss=1.         ! snow emissivity
real, save         :: maxsnowalpha=0.85    ! max snow albedo
real, save         :: minsnowalpha=0.5     ! min snow albedo
real, save         :: maxsnowden=300.      ! max snow density (kg m^-3)
real, save         :: minsnowden=100.      ! min snow density (kg m^-3)
! generic urban parameters
real, save         :: refheight=0.6        ! Displacement height as a fraction of building height (Kanda et al 2007)
real, save         :: zomratio=0.10        ! Ratio of roughness length to building height (default=0.1 or 10%)
real, save         :: zocanyon=0.01        ! Roughness length of in-canyon surfaces (m)
real, save         :: zoroof=0.01          ! Roughness length of roof surfaces (m)
real, save         :: maxrfwater=1.        ! Maximum roof water (kg m^-2)
real, save         :: maxrdwater=1.        ! Maximum road water (kg m^-2)
real, save         :: maxrfsn=1.           ! Maximum roof snow (kg m^-2)
real, save         :: maxrdsn=1.           ! Maximum road snow (kg m^-2)
real, save         :: maxvwatf=0.1         ! Factor multiplied to LAI to predict maximum leaf water (kg m^-2)
real, save         :: acfactor=2.          ! Air conditioning inefficiency factor
real, save         :: ac_heatcap=10.       ! Maximum heating/cooling capacity (W m^-3)
real, save         :: ac_coolcap=10.       ! Maximum heating/cooling capacity (W m^-3)
real, save         :: ac_deltat=1.         ! Comfort range for temperatures (+-K)
! atmosphere stability parameters
integer, save      :: icmax=5              ! number of iterations for stability functions (default=5)
real, save         :: a_1=1.
real, save         :: b_1=2./3.
real, save         :: c_1=5.
real, save         :: d_1=0.35
! diurnal cycle profiles (start & end at 1am)
! traffic diurnal cycle weights approximated from Chapman et al., 2016
real, dimension(25), save :: cyc_traf = (/ 0.17, 0.12, 0.12, 0.17, 0.37, 0.88, & 
                                           1.29, 1.48, 1.37, 1.42, 1.5 , 1.52, &
                                           1.50, 1.57, 1.73, 1.84, 1.84, 1.45, &
                                           1.01, 0.77, 0.65, 0.53, 0.41, 0.27, 0.17 /)
! base electricity demand cycle weights from Thatcher (2007), mean for NSW, VIC, QLD, SA
real, dimension(25), save :: cyc_base = (/ 0.92, 0.86, 0.81, 0.78, 0.8 , 0.87, & 
                                           0.98, 1.06, 1.08, 1.09, 1.09, 1.09, &
                                           1.08, 1.08, 1.06, 1.06, 1.08, 1.11, &
                                           1.08, 1.06, 1.03, 1.00, 0.98, 0.95, 0.92 /)
! normalised proportion of heating/cooling appliances in use approximated from Thatcher (2007)
real, dimension(25), save :: cyc_prop = (/ 0.43, 0.40, 0.36, 0.33, 0.31, 0.36, &
                                           0.53, 0.66, 0.70, 0.64, 0.60, 0.60, &
                                           0.64, 0.69 ,0.75, 0.81, 0.88, 1.00, &
                                           0.99, 0.91, 0.81, 0.69, 0.57, 0.43, 0.43 /)
! internal temperature translation cycle from Thatcher (2007) (only used for intairtmeth=0)
real, dimension(25), save :: cyc_tran = (/ -1.09, -1.21, -2.12, -2.77, -3.06, -2.34, &
                                           -0.37,  1.03,  1.88,  2.37,  2.44,  2.26, &
                                            1.93,  1.41,  0.74,  0.16,  0.34,  1.48, & 
                                            1.03,  0.14, -0.74, -1.17, -1.15, -1.34, -1.09/)

end module uclem_parameters
