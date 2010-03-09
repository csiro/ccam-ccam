MODULE radiation_module
  ! use Goudriaan's radiation scheme for calculate radiation absorbed by canopy 
  ! and soil, treat diffuse, direct separately for three wavebands (nrb),
  ! nrb=1, visible; =2 for nir and 3 for thermal
  ! input variables
  !  fsd: incoming shortwave (0.5 visible, 0.5 nir)
  !  fld: incoming longwave radiation
  !  veg%vlai: canopy LAI
  !  veg%xfang: leaf inclination angle distrbution parameter (<0 more vertical)
  !                                                          =0 spherical
  !                                                          >0 more horizontal)
  !  rad%latitude: site latitude in degree (- for SH and + for NH)
  !  soil%ref: soil reflectance
  !  taul: leaf transmittance
  !  rhol: leaf reflectance
  ! output varibales
  !  qcan: absorbed radiation by canopy
  !  qssabs: absorbed radiation (shortwave only) by soil
  !  scalex: scaling for sunlit/shaded leaves
  !  transd: diffuse transmittance of canopy
  !  rad%fvlai: LAI of sunlit/shaded leaves
  !  rad%albedo: surface (canopy+soil) albedo
  !  rad%gradis: isothermal radiative conductance
  USE math_constants
  USE other_constants
  USE define_types
  USE physical_constants
  IMPLICIT NONE
  ! This module contains the following subroutines:
  PRIVATE
  PUBLIC init_radiation, radiation, sinbet ! available outside this module
  PRIVATE spitter ! available only from within this module
CONTAINS
  !--------------------------------------------------------------------------
  SUBROUTINE init_radiation
    REAL(r_1), DIMENSION(mp,nrb) :: c1    ! sqrt(1. - taul - refl)
    REAL(r_1), DIMENSION(3)      :: cos3  ! cos(15 45 75 degrees)
!    REAL(r_1), DIMENSION(nrb)   :: rhoch ! canopy reflect'n black horiz leaves(6.19)
    REAL(r_1), DIMENSION(mp,nrb) :: rhoch ! canopy reflect'n black horiz leaves(6.19)
    REAL(r_1), DIMENSION(mp,3)  :: xvlai2 ! 2D vlai
    REAL(r_1), DIMENSION(mp,3)  :: xk    ! extinction coefficient for beam radiation and black leaves
    REAL(r_1), DIMENSION(mp)    :: xphi1 ! leaf angle parmameter 1
    REAL(r_1), DIMENSION(mp)    :: xphi2 ! leaf angle parmameter 2
    INTEGER(i_d)                :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave

    cos3 = COS(pi180 * (/ 15.0, 45.0, 75.0 /))
    WHERE (canopy%vlaiw > 1e-2)
       ! See Sellers 1985, eq.13 (leaf angle parameters):
       xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
       xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
    END WHERE
    ! 2 dimensional LAI
    xvlai2 = SPREAD(canopy%vlaiw, 2, 3)
!    print *,'initrad',canopy%vlaiw,xvlai2,veg%xfang
!    print *,'initrad1',xphi1,xphi2,cos3,pi180
    WHERE (xvlai2 > 1e-2) ! vegetated
       xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
    ELSEWHERE ! i.e. bare soil
       xk = 0.0		 
    END WHERE
!    print *,'initrad2',xk
    WHERE (canopy%vlaiw > 1e-2) ! vegetated
       ! Extinction coefficient for diffuse radiation for black leaves:
       rad%extkd = -LOG(SUM(SPREAD(gauss_w, 1, mp) * EXP(-xk * xvlai2), 2)) / canopy%vlaiw
    ELSEWHERE ! i.e. bare soil
       rad%extkd = 0.7
    END WHERE
!    where( veg%iveg == 1 ) rad%extkd = rad%extkd * 0.7
   ! Extinction coefficient for leaf nitrogen profile in canopy:
    rad%extkn = rad%extkd * SQRT(1.0 - veg%taul(:,1) - veg%refl(:,1))
    
!   rad%fbeam = spitter(met%doy, met%coszen, met%fsd)
    WHERE (met%coszen <1.0e-2)
       rad%fbeam(:,1) = 0.0
       rad%fbeam(:,2) = 0.0
       rad%fbeam(:,3) = 0.0
    END WHERE
    WHERE (canopy%vlaiw > 1e-2)    ! In gridcells where vegetation exists....
      ! SW beam extinction coefficient ("black" leaves, extinction neglects
      ! leaf SW transmittance and reflectance):
      rad%extkb = xphi1 / met%coszen + xphi2
    ELSEWHERE ! i.e. bare soil
      rad%extkb = 0.5
    END WHERE
!    where( veg%iveg == 1 ) rad%extkb = rad%extkb * 0.7
    WHERE ( abs(rad%extkb - rad%extkd)  < 0.001 )
      rad%extkb = rad%extkd + 0.001
    END WHERE
    WHERE(rad%fbeam(:,3) < 1.0e-4)
      rad%extkb=1.0e5
    END WHERE
!    c1 = SQRT(1. - taul - refl)
    c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
    c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
    c1(:,3) = 1.
    ! Canopy reflection black horiz leaves (eq. 6.19 in Goudriaan and van Laar, 1994):
    rhoch = (1.0 - c1) / (1.0 + c1)
    ! Canopy reflection of diffuse radiation for black leaves:
!    rad%rhocdf = SPREAD(2.0 * gauss_w * rhoch, 1, mp) * xk / (xk + SPREAD(rad%extkd, 2, 3))
    rad%rhocdf = SPREAD(2.0 * gauss_w , 1, mp) * rhoch * xk / (xk + SPREAD(rad%extkd, 2, 3))
!    print *,'rhocdf',rad%rhocdf
    do b=1,nrb
      rad%rhocdf(:,b) = rhoch(:,b) * ( gauss_w(1)*xk(:,1)/( xk(:,1) + rad%extkd(:) ) +  &
                                       gauss_w(2)*xk(:,2)/( xk(:,2) + rad%extkd(:) ) +  &
                                       gauss_w(3)*xk(:,3)/( xk(:,3) + rad%extkd(:) )  )
    enddo
!    print *,'rhocdfnew',rad%rhocdf
!    print *,'initrad4',rad%extkn,veg%taul,veg%refl,c1
!    print *,'initrad5',rhoch
  END SUBROUTINE init_radiation
  !-------------------------------------------------------------------------------
  SUBROUTINE radiation(ktau)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp)  :: cf1	! (1.0 - transb * cexpkdm) / (extkb + extkdm(:,b))
    REAL(r_1), DIMENSION(mp)  :: cf3	! (1.0 - transb * cexpkbm) / (extkb + extkbm(:,b))
    REAL(r_1), DIMENSION(mp)  :: cf2n	! exp(-extkn * vlai) (nitrogen)
    REAL(r_1), DIMENSION(mp)  :: emair	! air emissivity
    REAL(r_1), DIMENSION(mp)  :: flpwb	! black-body long-wave radiation
    REAL(r_1), DIMENSION(mp)  :: flwv	! vegetation long-wave radiation (isothermal)
!    LOGICAL, DIMENSION(mp)    :: ask	! select points for calculation
    INTEGER(i_d)	      :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
    LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation
    INTEGER(i_d) :: k

    ! Define vegetation mask:
    mask = canopy%vlaiw > 1e-2 .AND. met%fsd(:,3) > 1.0e-2
    ! Relative leaf nitrogen concentration within canopy:
    cf2n = EXP(-rad%extkn * canopy%vlaiw)
    rad%transd = 1.0
    WHERE (canopy%vlaiw > 1e-2)    ! In gridcells where vegetation exists....
       ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
       ! leaf SW transmittance and reflectance);
       ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
       rad%transd = EXP(-rad%extkd * canopy%vlaiw)
!    ELSEWHERE	! i.e. bare soil
!       rad%transd = 1.0
    END WHERE
    ! Define fraction of SW beam tranmitted through canopy:
    rad%transb = EXP(-rad%extkb * canopy%vlaiw)
    ! Define longwave from vegetation:
    flpwb = sboltz * (met%tk) ** 4
    flwv = emleaf * flpwb
    ! Combined soil/snow temperature:
    ssoil%tss=(1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    ! Define longwave from soil/snow surface:
    rad%flws = sboltz*emsoil* ssoil%tss **4
    ! Define air emissivity:
    emair = met%fld / flpwb
!    print 21,ktau,emair,met%fld,flpwb,rad%trad,met%tk
21  format(1x,'emair',i5,2f6.2,4f6.0,4f6.1)
    rad%gradis = 0.0 ! initialise radiative conductance
    rad%qcan = 0.0   ! initialise radiation absorbed by canopy
    WHERE (canopy%vlaiw > 1.e-2)
       ! Define radiative conductance (Leuning et al, 1995), eq. D7:
       rad%gradis(:,1) = (4.0 * emleaf / (capp * air%rho)) * flpwb / &
            (met%tvrad) * rad%extkd * &
            ((1.0 - rad%transb * rad%transd) / (rad%extkb + rad%extkd) + &
            (rad%transd - rad%transb) / (rad%extkb - rad%extkd))
       rad%gradis(:,2) = (8.0 * emleaf / (capp * air%rho)) * flpwb / (met%tvrad) * &
            rad%extkd * (1.0 - rad%transd) / rad%extkd - rad%gradis(:,1)
       ! Longwave radiation absorbed by sunlit canopy fraction:
       rad%qcan(:,1,3) = (rad%flws-flwv) *rad%extkd * (rad%transd - rad%transb) / &
            (rad%extkb - rad%extkd) &
!            + (met%fld-emleaf*flpwb) * rad%extkd * (1.0 - rad%transd * rad%transb) / &
            + (emair - emleaf) * rad%extkd*flpwb * (1.0 - rad%transd*rad%transb) / &
            ( rad%extkb + rad%extkd ) 
       ! Longwave radiation absorbed by shaded canopy fraction:
       rad%qcan(:,2,3) = (1.0 - rad%transd) * &
            (rad%flws + met%fld - 2.0 * flwv) - rad%qcan(:,1,3)
!            (rad%flws + flpwb * (emair - 2.0 * emleaf)) - rad%qcan(:,1,3) 
    END WHERE
    ! Convert radiative conductance from m/s to mol/m2/s:
    rad%gradis = SPREAD(air%cmolar, 2, mf)*rad%gradis
    rad%gradis = MAX(1.e-3,rad%gradis)

!    print *,'rad7', rad%gradis

    DO b = 1, 2	! 1 = visible, 2 = nir radiaition

!    print *,'rad8', b,rad%transb,rad%cexpkdm(:,b)
!    print *,'rad8a', b,rad%extkb,rad%extkdm,rad%extkbm

       WHERE (mask) ! i.e. vegetation and sunlight are present
          cf1 = (1.0 - rad%transb * rad%cexpkdm(:,b)) / (rad%extkb + rad%extkdm(:,b))
          cf3 = (1.0 - rad%transb * rad%cexpkbm(:,b)) / (rad%extkb + rad%extkbm(:,b))
          ! Define qcan for short wave (par, nir) for sunlit leaf:
          rad%qcan(:,1,b) = met%fsd(:,b) * ( &                ! scale to real sunlit flux
               (1.0-rad%fbeam(:,b))*(1.0-rad%reffdf(:,b))*rad%extkdm(:,b)*cf1 &
               + rad%fbeam(:,b)*(1.0-rad%reffbm(:,b))*rad%extkbm(:,b)*cf3 &
               + rad%fbeam(:,b)*(1.0-veg%taul(:,b)-veg%refl(:,b))*rad%extkb &
               * ((1-rad%transb)/rad%extkb - (1-rad%transb**2)/(rad%extkb+rad%extkb)))
          ! Define qcan for short wave (par, nir) for shaded leaf:
          rad%qcan(:,2,b) = met%fsd(:,b) * ( &                ! scale to real shaded flux
               (1.0-rad%fbeam(:,b))*(1.0-rad%reffdf(:,b))*rad%extkdm(:,b)* &
               ((1.0 - rad%cexpkdm(:,b)) / rad%extkdm(:,b) - cf1) &
               + rad%fbeam(:,b)*(1.-rad%reffbm(:,b))*rad%extkbm(:,b) &
               * ((1.0 - rad%cexpkbm(:,b)) / rad%extkbm(:,b) - cf3) &
               - rad%fbeam(:,b)*(1.0-veg%taul(:,b)-veg%refl(:,b))*rad%extkb &
               * ((1-rad%transb)/rad%extkb - (1-rad%transb**2)/(rad%extkb+rad%extkb)))
       END WHERE
    END DO

!    print *,'rad9a',rad%qcan
!    print 101,ktau,rad%albedo(1,2),rad%albedo(2,2),fbeam,reffbm(:,1),reffdf(:,2),canopy%vlaiw,met%coszen
101 format(1x,'radalbedo',i6,8f6.3,2f5.2,2f6.3)

!les    ! Define IR albedo - CURRENTLY NOT USED elsewhere
!    rad%albedo(:,3) = 0.05

    WHERE (mask) ! i.e. vegetation and sunlight are present
       ! Calculate shortwave radiation absorbed by soil:
       ! (av. of transmitted NIR and PAR through canopy)*SWdown
       rad%qssabs =  met%fsd(:,1) * ( &
            rad%fbeam(:,1)*(1.-rad%reffbm(:,1))*EXP(-rad%extkbm(:,1)*canopy%vlaiw) &
            +(1.-rad%fbeam(:,1))*(1.-rad%reffdf(:,1))*EXP(-rad%extkdm(:,1)*canopy%vlaiw)) + &
            met%fsd(:,2)*(                                                               &
            rad%fbeam(:,2)*(1.-rad%reffbm(:,2))*rad%cexpkbm(:,2) + &
            (1.-rad%fbeam(:,2))*(1.-rad%reffdf(:,2))*rad%cexpkdm(:,2))
       ! Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C:
       rad%scalex(:,1) = (1.0 - rad%transb * cf2n) / (rad%extkb + rad%extkn)
       ! Leaf area index of big leaf, sunlit, shaded, respectively:
       rad%fvlai(:,1) = (1.0 - rad%transb) / rad%extkb
       rad%fvlai(:,2) = canopy%vlaiw - rad%fvlai(:,1)
    ELSEWHERE ! i.e. either vegetation or sunlight are NOT present
       ! Shortwave absorbed by soil/snow surface:
       rad%qssabs = (1.0 - (0.5 *(ssoil%albsoilsn(:,1) + ssoil%albsoilsn(:,2)))) * met%fsd(:,3)
       rad%scalex(:,1) = 0.0
       rad%fvlai(:,1) = 0.0
       rad%fvlai(:,2) = canopy%vlaiw
    END WHERE
    rad%scalex(:,2) = (1.0 - cf2n) / rad%extkn - rad%scalex(:,1)
    ! Total energy absorbed by canopy:
    rad%rniso = sum(rad%qcan, 3)

!    print *,'rad10',rad%qssabs,rad%rniso
  END SUBROUTINE radiation
  !-----------------------------------------------------------------
  ELEMENTAL FUNCTION sinbet(doy,xslat,hod) RESULT(z)
    ! calculate sin(bet), bet = elevation angle of sun
    ! calculations according to goudriaan & van laar 1994 p30
    REAL(r_1), INTENT(IN)	:: doy		! day of year
    REAL(r_1), INTENT(IN)	:: xslat	! latitude (degrees north)
    REAL(r_1), INTENT(IN)	:: hod		! hour of day
    REAL(r_1)			:: sindec	! sine of maximum declination
    REAL(r_1)			:: z		! result
    sindec = -SIN(23.45 * pi180) * COS(two_pi * (doy + 10.0) / 365.0)
    z = MAX( &
         SIN(pi180 * xslat) * sindec + COS(pi180 * xslat) * SQRT(1. - sindec * sindec) &
         * COS(pi_c * (hod - 12.) / 12.), 1e-8)
!    print *,'from sinbet',doy,xslat,hod,z
  END FUNCTION sinbet
  !-------------------------------------------------------------------
  FUNCTION spitter(doy, coszen, fsd) RESULT(fbeam)
    ! Calculate beam fraction
    ! See spitters et al. 1986, agric. for meteorol., 38:217-229
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: doy	! day of year
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: coszen ! cos(zenith angle of sun)
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: fsd	! short wave down (positive) w/m^2
    REAL(r_1), DIMENSION(mp) :: fbeam	! beam fraction (result)
    REAL(r_1), PARAMETER :: solcon = 1370.0
    REAL(r_1), DIMENSION(mp) :: tmpr !
    REAL(r_1), DIMENSION(mp) :: tmpk !
    REAL(r_1), DIMENSION(mp) :: tmprat !
!    print *,'spitter',doy,coszen,fsd,two_pi
    fbeam = 0.0
    tmpr = 0.847 + coszen * (1.04 * coszen - 1.61)
    tmpk = (1.47 - tmpr) / 1.66
    WHERE (coszen > 1.0e-10 .AND. fsd > 10.0)
       tmprat = fsd / (solcon * (1.0 + 0.033 * COS(two_pi * (doy-10.0) / 365.0)) * coszen)
    ELSEWHERE
       tmprat = 0.0
    END WHERE
!    print *,'tmprat',tmprat
    WHERE (tmprat > 0.22) fbeam = 6.4 * (tmprat - 0.22) ** 2
    WHERE (tmprat > 0.35) fbeam = MIN(1.66 * tmprat - 0.4728, 1.0)
    WHERE (tmprat > tmpk) fbeam = MAX(1.0 - tmpr, 0.0)
  END FUNCTION spitter
END MODULE radiation_module
