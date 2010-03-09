MODULE albedo_module
  ! use Goudriaan's radiation scheme for calculate radiation absorbed by canopy 
  ! and soil, treat diffuse, direct separately for three wavebands (nrb),
  ! nrb=1, visible; =2 for nir and  3 for thermal
  ! input variables
  !  fsd: incoming shortwave (0.5 visible, 0.5 nir)
  !  veg%vlai: canopy LAI
  !  veg%xfang: leaf inclination angle distrbution parameter (<0 more vertical)
  !                                                          =0 spherical
  !                                                          >0 more horizontal)
  !  ssoil%albsoilsn: soil+snow albedo 
  !  taul: leaf transmittance
  !  rhol: leaf reflectance
  ! output varibales
  ! rad%albedo
  ! reffdf
  ! reffbm
  USE math_constants
  USE other_constants
  USE define_types
  USE physical_constants
!  USE cable_variables

  IMPLICIT NONE
  ! This module contains the following subroutines:
  PUBLIC albedo
!  PRIVATE spitter ! available only from within this module
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE albedo(ktau)
!  SUBROUTINE albedo(ktau,ssoil, veg, air, met, rad,canopy)
!    TYPE (soil_snow_type),INTENT(INOUT)	:: ssoil
!    TYPE (veg_parameter_type),INTENT(IN):: veg
!    TYPE (air_type),INTENT(IN)	        :: air
!    TYPE (met_type),INTENT(INOUT)	:: met
!    TYPE (radiation_type),INTENT(INOUT)	:: rad
!    TYPE (canopy_type),INTENT(INOUT)    :: canopy
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number

!    REAL(r_1), DIMENSION(nrb) :: c1	! sqrt(1. - taul - refl)
    REAL(r_1), DIMENSION(mp,nrb) :: c1        ! sqrt(1. - taul - refl)

!    REAL(r_1), DIMENSION(mp)  :: cexpkbm ! canopy beam transmittance
!    REAL(r_1), DIMENSION(mp)  :: cexpkdm ! canopy diffuse transmittance
!    REAL(r_1), DIMENSION(mp,nrb) :: extkbm	! modified k beam(6.20)(for leaf scattering)
!    REAL(r_1), DIMENSION(mp,nrb) :: extkdm ! modified k diffuse(6.20)(for leaf scattering)
!    REAL(r_1), DIMENSION(mp)  :: fbeam	! beam fraction
    INTEGER(i_d)	      :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
!    REAL(r_1), DIMENSION(mp,nrb) :: reffbm	! effective conopy beam reflectance
!    REAL(r_1), DIMENSION(mp,nrb) :: reffdf	! effective conopy diffuse reflectance
!    REAL(r_1), DIMENSION(nrb)    :: rhoch ! canopy reflection black horizontal leaves(6.19)
     REAL(r_1), DIMENSION(mp,nrb)    :: rhoch ! canopy reflection black horizontal leaves(6.19)

!    REAL(r_1), DIMENSION(mp,nrb) :: rhocbm	! modified canopy beam reflectance (6.21)
!    REAL(r_1), DIMENSION(mp)     :: xphi1	! leaf angle parmameter 1
!    REAL(r_1), DIMENSION(mp)     :: xphi2	! leaf angle parmameter 2
    LOGICAL, DIMENSION(mp)       :: mask     ! select points for calculation
    INTEGER(i_d) :: k
   
    ! coszen is set during met data read in.

    
    rad%reffbm = ssoil%albsoilsn   ! initialise effective conopy beam reflectance
    rad%reffdf = ssoil%albsoilsn   ! initialise effective conopy beam reflectance

    ! Define beam fraction, fbeam:
!    print *,'rad 1',met%doy,met%coszen,met%fsd,met%fld
!    fbeam = spitter(met%doy, met%coszen, met%fsd)
!    print *,'rad 1',met%coszen,fbeam
!    WHERE (met%coszen <1.0e-2)
!       fbeam = 0.0
!    END WHERE
    ! Define vegetation mask:
    mask = canopy%vlaiw > 1e-2 .AND. met%fsd(:,3) > 1.0e-2
    ! See Sellers 1985, eq.13 (leaf angle parameters):
!    xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
!    xphi2 = 0.877 - (0.877 * 2.0) * xphi1
!    print *,'rad 2',mask,canopy%vlaiw,rad%extkn,veg%xfang,xphi1,xphi2
!    WHERE (canopy%vlaiw > 1e-2)    ! In gridcells where vegetation exists....
!       ! SW beam extinction coefficient ("black" leaves, extinction neglects
       ! leaf SW transmittance and reflectance):
!       rad%extkb = xphi1 / met%coszen + xphi2
!    ELSEWHERE	! i.e. bare soil
!       rad%extkb = 0.5
!    END WHERE
!    WHERE(fbeam < 1.0e-3)
!       rad%extkb=1.
!    END WHERE

!    print *,'albedo',ktau
!    c1 = SQRT(1. - taul - refl)
    c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
    c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
    c1(:,3) = 1.
    ! Define canopy reflection black horizontal leaves(6.19)
    rhoch = (1.0 - c1) / (1.0 + c1)
    ! Update extinction coefficients and fractional transmittance for 
    ! leaf transmittance and reflection (ie. NOT black leaves):
    DO b = 1, 2	! 1 = visible, 2 = nir radiaition
      rad%extkdm(:,b) = rad%extkd * c1(:,b)
       ! Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)
       ! Calculate effective diffuse reflectance (fraction):
       rad%reffdf(:,b) = rad%rhocdf(:,b) + &
            (ssoil%albsoilsn(:,b) - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2
!       print *,'RADIAT',ktau,b,rad%rhocdf(:,b),ssoil%albsoilsn(:,b),cexpkdm,fbeam,reffdf(:,b)
       WHERE (mask) ! i.e. vegetation and sunlight are present
          rad%extkbm(:,b) = rad%extkb * c1(:,b)
          ! Canopy reflection (6.21) beam:
          rad%rhocbm(:,b) = 2.*rad%extkb/(rad%extkb+rad%extkd)*rhoch(:,b)
          ! Canopy beam transmittance (fraction):
          rad%cexpkbm(:,b) = EXP(-rad%extkbm(:,b)*canopy%vlaiw)
          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssoil%albsoilsn(:,b) - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2
       END WHERE
       ! Define albedo:
       rad%albedo(:,b) = (1.0-rad%fbeam(:,b))*rad%reffdf(:,b)+rad%fbeam(:,b)*rad%reffbm(:,b)
    END DO
    rad%reffdf(:,2)=rad%reffdf(:,2)*veg%xalbnir(:)
    rad%reffbm(:,2)=rad%reffbm(:,2)*veg%xalbnir(:)
    rad%albedo(:,2) = rad%albedo(:,2) * veg%xalbnir(:)

!    print 101,ktau,rad%albedo(1,2),rad%albedo(2,2),fbeam,reffbm(:,1),reffdf(:,2),canopy%vlaiw,met%coszen
!    print 101,rad%albedo(1,2),rad%albedo(2,2),fbeam,reffbm(:,1),reffdf(:,2)
!     print 101,ktau,rad%albedo(1,1),rad%albedo(2,1),rad%albedo(1,2),rad%albedo(2,2),&
!                         rad%reffdf(1,1),rad%reffdf(2,1),rad%reffdf(1,2),rad%reffdf(2,2), &
!                         rad%reffbm(1,1),rad%reffbm(2,1),rad%reffbm(1,2),rad%reffbm(2,2),rad%fbeam
101 format(1x,'atm_ph1albedo',i6,12f6.3,2f6.3)
!101 format(1x,'atm_ph2albedo',i6,8f6.3,2f5.2,2f6.3)
  END SUBROUTINE albedo 
END MODULE albedo_module
