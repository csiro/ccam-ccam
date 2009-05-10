MODULE cab_albedo_module
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
  ! rad%reffdf
  ! rad%reffbm
  ! rad%extkdm,rad%extkbm,rad%cexpkdm,rad%cexpkbm,rad%rhocbm,
  USE math_constants
  USE other_constants
  USE define_types
  USE physical_constants
  USE cable_variables

  IMPLICIT NONE
  ! This module contains the following subroutines:
  PUBLIC cab_albedo
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE cab_albedo(istep_cur,dels,ssoil, veg, air, met, rad,  &
                        soil, L_RADUM)
    TYPE (soil_snow_type),INTENT(INOUT)	:: ssoil
    TYPE (veg_parameter_type),INTENT(IN):: veg
    TYPE (air_type),INTENT(IN)	        :: air
    TYPE (met_type),INTENT(INOUT)	:: met
    TYPE (radiation_type),INTENT(INOUT)	:: rad
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil   

    REAL(r_1), INTENT(IN)           :: dels      !integration time step (s)
    INTEGER(i_d), INTENT(IN)        :: istep_cur ! current dt
    LOGICAL, INTENT(IN)             :: L_RADUM   ! true if called from HADGEM
    REAL(r_1), DIMENSION(nrb) :: c1	! sqrt(1. - taul - refl)
    REAL(r_1), DIMENSION(nrb) :: rhoch  ! canopy reflection black horizontal leaves(6.19)
    REAL(r_1), DIMENSION(mp)    :: alv ! Snow albedo for visible
    REAL(r_1), DIMENSION(mp)    :: alir ! Snow albedo for near infra-red
    REAL(r_1), PARAMETER        :: alvo  = 0.95 ! albedo for vis. on a new snow
    REAL(r_1), PARAMETER        :: aliro = 0.65 ! albedo for near-infr. on a new snow
    REAL(r_1), DIMENSION(mp)    :: ar1 ! crystal growth  (-ve)
    REAL(r_1), DIMENSION(mp)    :: ar2 ! freezing of melt water
    REAL(r_1), DIMENSION(mp)    :: ar3
    REAL(r_1), DIMENSION(mp)    :: dnsnow ! new snow albedo
    REAL(r_1), DIMENSION(mp)    :: dtau
    REAL(r_1), DIMENSION(mp)    :: fage !age factor
    REAL(r_1), DIMENSION(mp)    :: fzenm
    REAL(r_1), DIMENSION(mp)    :: sfact
    REAL(r_1), DIMENSION(mp)    :: snr
    REAL(r_1), DIMENSION(mp)    :: snage_UM
    REAL(r_1), DIMENSION(mp)    :: snrat
    REAL(r_1), DIMENSION(mp)    :: talb ! snow albedo
    REAL(r_1), DIMENSION(mp)    :: tmp ! temporary value

    INTEGER(i_d)            :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
    LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation
    INTEGER(i_d)            :: k
   
    ! coszen is set during met data read in.
!==sxy
    !    calculate soil/snow albedo
!    print *,'in cab_alb',soil%albsoil
    sfact = 0.68
    WHERE (soil%albsoil <= 0.14)
       sfact = 0.5
    ELSEWHERE (soil%albsoil > 0.14 .and. soil%albsoil <= 0.20)
       sfact = 0.62
    END WHERE
    ssoil%albsoilsn(:,2) = 2. * soil%albsoil / (1. + sfact)
    ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
    !          new snow albedo (needs osnowd from the previous dt)
!    print *,'in cab_alb2',ssoil%albsoilsn,ssoil%snage
    snrat=0.
    alir =0.
    alv  =0.
!    print *,'in cab_alb3',ssoil%snowd, ssoil%osnowd,ssoil%isflag
    WHERE (ssoil%snowd > 0 .and. .NOT. L_RADUM ) 
       dnsnow = min (1., .1 * max (0., ssoil%snowd - ssoil%osnowd ) ) ! new snow (cm H2O)
       !         Snow age depends on snow crystal growth, freezing of melt water,
       !         accumulation of dirt and amount of new snow.
       tmp = ssoil%isflag * ssoil%tggsn(:,1) + (1 - ssoil%isflag ) * ssoil%tgg(:,1)
       tmp = min (tmp, 273.15)
       ar1 = 5000. * (1. / 273.15 - 1. / tmp) ! crystal growth  (-ve)
       ar2 = 10. * ar1 ! freezing of melt water
       snr = ssoil%snowd / max (ssoil%ssdnn, 100.)
       WHERE (ssoil%snowd < 10. ) snr = ssoil%snowd / max (ssoil%ssdnn, 400.)
          ! fixes for Arctic & Antarctic
       WHERE (soil%isoilm == 9)
          ar3 = .0005
          dnsnow = max (dnsnow, .002) !increase refreshing of snow in Antarctic
          snrat = min (1., snr / (snr + .001) )
       ELSEWHERE
          ! accumulation of dirt
          ar3 = .1
          snrat = min (1., snr / (snr + .01) )
       END WHERE
          dtau = 1.e-6 * (exp(ar1) + exp(ar2) + ar3) * dels 
       WHERE (ssoil%snowd <= 1.0)
          ssoil%snage = 0.
          snage_UM = 0.
       ELSEWHERE
          ssoil%snage = max (0.,(ssoil%snage+dtau)*(1.-dnsnow))
          snage_UM = ssoil%snage
       END WHERE
!       fage = 1. - 1. / (1. + snage_UM ) !age factor
       fage = 1. - 1. / (1. + ssoil%snage ) !age factor
       !
       !       Snow albedo is dependent on zenith angle and  snow age.
       !       albedo zenith dependence
       !       alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
       !               where cs = 0.2, cn = 0.5, b = 2.0
       tmp = max (.17365, met%coszen )
       fzenm = max(merge(0.0, (1. + 1. / 2.) / (1. + 2. * 2. * tmp) - 1. / 2., tmp > 0.5), 0.)
       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo
    ENDWHERE        ! snowd > 0
    WHERE (ssoil%snowd > 0 .and.  L_RADUM )
       dnsnow = min (1., .1 * max (0., ssoil%snowd - ssoil%osnowd ) ) ! new snow (cm H2O)
       !         Snow age depends on snow crystal growth, freezing of melt water,
       !         accumulation of dirt and amount of new snow.
       tmp = ssoil%isflag * ssoil%tggsn(:,1) + (1 - ssoil%isflag ) * ssoil%tgg(:,1)
       tmp = min (tmp, 273.15)
       ar1 = 5000. * (1. / 273.15 - 1. / tmp) ! crystal growth  (-ve)
       ar2 = 10. * ar1 ! freezing of melt water
       snr = ssoil%snowd / max (ssoil%ssdnn, 100.)
       WHERE (ssoil%snowd < 10. ) snr = ssoil%snowd / max (ssoil%ssdnn, 400.)
          ! fixes for Arctic & Antarctic
       WHERE (soil%isoilm == 9)
          ar3 = .0005
          dnsnow = max (dnsnow, .002) !increase refreshing of snow in Antarctic
          snrat = min (1., snr / (snr + .001) )
       ELSEWHERE
          ! accumulation of dirt
          ar3 = .1
          snrat = min (1., snr / (snr + .01) )
       END WHERE
          dtau = 1.e-6 * (exp(ar1) + exp(ar2) + ar3) * dels
       WHERE (ssoil%snowd <= 1.0)
          ssoil%snage = 0.
          snage_UM = 0.
       ELSEWHERE
          snage_UM = ssoil%snage
       END WHERE
       fage = 1. - 1. / (1. + snage_UM ) !age factor
       !
       !       Snow albedo is dependent on zenith angle and  snow age.
       !       albedo zenith dependence
       !       alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
       !               where cs = 0.2, cn = 0.5, b = 2.0
       tmp = max (.17365, met%coszen )
       fzenm = max(merge(0.0, (1. + 1. / 2.) / (1. + 2. * 2. * tmp) - 1. / 2., tmp > 0.5), 0.)
       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo
    ENDWHERE        ! snowd > 0


    ssoil%albsoilsn(:,2) = (1. - snrat) * ssoil%albsoilsn(:,2) + snrat * alir
    ssoil%albsoilsn(:,1) = (1. - snrat) * ssoil%albsoilsn(:,1) + snrat * alv
!    print *,'ALBSOILSN',snrat,ssoil%albsoilsn(:,1),ssoil%albsoilsn(:,2),ssoil%snowd,ssoil%ssdnn
!==sxy
    
    rad%reffbm = ssoil%albsoilsn   ! initialise effective conopy beam reflectance
    rad%reffdf = ssoil%albsoilsn   ! initialise effective conopy beam reflectance
!    print *,'in cab_alb5', ssoil%albsoilsn,rad%reffbm,rad%reffdf
    ! Define vegetation mask:
    mask = veg%vlaiw > 1e-2 .AND. met%fsd(:,3) > 1.0e-2

    c1 = SQRT(1. - taul - refl)
    ! Define canopy reflection black horizontal leaves(6.19)
    rhoch = (1.0 - c1) / (1.0 + c1)
    ! Update extinction coefficients and fractional transmittance for 
    ! leaf transmittance and reflection (ie. NOT black leaves):
    DO b = 1, 2	! 1 = visible, 2 = nir radiaition
      rad%extkdm(:,b) = rad%extkd * c1(b)
       ! Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * veg%vlaiw)
       ! Calculate effective diffuse reflectance (fraction):
       rad%reffdf(:,b) = rad%rhocdf(:,b) + &
            (ssoil%albsoilsn(:,b) - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2
!       print *,'RADIAT',ktau,b,rad%rhocdf(:,b),ssoil%albsoilsn(:,b),cexpkdm,fbeam,reffdf(:,b)
       WHERE (mask) ! i.e. vegetation and sunlight are present
          rad%extkbm(:,b) = rad%extkb * c1(b)
          ! Canopy reflection (6.21) beam:
          rad%rhocbm(:,b) = 2.*rad%extkb/(rad%extkb+rad%extkd)*rhoch(b)
          ! Canopy beam transmittance (fraction):
          rad%cexpkbm(:,b) = EXP(-rad%extkbm(:,b)*veg%vlaiw)
          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssoil%albsoilsn(:,b) - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2
       END WHERE
       ! Define albedo:
       rad%albedo(:,b) = (1.0-rad%fbeam(:,b))*rad%reffdf(:,b)+rad%fbeam(:,b)*rad%reffbm(:,b)
    END DO
!    print *,'in cab_alb8',rad%albedo
!    print 101,istep_cur,rad%albedo(1,1),rad%albedo(2,1),rad%albedo(1,2),rad%albedo(2,2),&
!                         rad%reffdf(1,1),rad%reffdf(2,1),rad%reffdf(1,2),rad%reffdf(2,2), &
!                         rad%reffbm(1,1),rad%reffbm(2,1),rad%reffbm(1,2),rad%reffbm(2,2),rad%fbeam
!101 format(1x,'alb_atmph1',i6,12f6.3,6f6.3)
  END SUBROUTINE cab_albedo 
END MODULE cab_albedo_module
