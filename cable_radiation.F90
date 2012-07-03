
MODULE radiation_module
  USE math_constants
  USE other_constants
  USE define_types
  USE define_dimensions
  USE physical_constants
  IMPLICIT NONE

  PRIVATE
  PUBLIC init_radiation, radiation, sinbet 

CONTAINS

   SUBROUTINE init_radiation(met,rad,veg,canopy)
      use cable_common_module
      implicit none
      TYPE (veg_parameter_type), INTENT(INout) :: veg
      TYPE (radiation_type), INTENT(INOUT) :: rad
      TYPE (met_type),INTENT(INOUT)        :: met
      TYPE (canopy_type),INTENT(IN)                :: canopy
      REAL(r_1), DIMENSION(nrb)     :: cos3  ! cos(15 45 75 degrees)
      REAL(r_1), DIMENSION(mp,nrb)  :: xvlai2 ! 2D vlai
      REAL(r_1), DIMENSION(mp,nrb)  :: xk    ! extinct. coef.for beam rad. and black leaves
      REAL(r_1), DIMENSION(mp)    :: xphi1 ! leaf angle parmameter 1
      REAL(r_1), DIMENSION(mp)    :: xphi2 ! leaf angle parmameter 2
      INTEGER(i_d) :: ictr
      LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation
      real, dimension(mp,nrb) :: c1, rhoch ! MJT suggestion

         cos3 = COS(pi180 * (/ 15.0, 45.0, 75.0 /))

         ! See Sellers 1985, eq.13 (leaf angle parameters):
         WHERE (canopy%vlaiw > LAI_THRESH)
            xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
            xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
         END WHERE
         ! 2 dimensional LAI
         xvlai2 = SPREAD(canopy%vlaiw, 2, 3)

         ! Extinction coefficient for beam radiation and black leaves;
         ! eq. B6, Wang and Leuning, 1998
         WHERE (xvlai2 > LAI_THRESH) ! vegetated
            xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
         ELSEWHERE ! i.e. bare soil
            xk = 0.0          
         END WHERE
     
         WHERE (canopy%vlaiw > LAI_THRESH ) ! vegetated
            ! Extinction coefficient for diffuse radiation for black leaves:
            rad%extkd = -LOG(SUM(SPREAD(gauss_w, 1, mp) * EXP(-xk * xvlai2), 2)) / canopy%vlaiw
         ELSEWHERE ! i.e. bare soil
            rad%extkd = 0.7
         END WHERE

         mask = canopy%vlaiw > LAI_THRESH  .AND.                                 &
               ( met%fsd(:,1) + met%fsd(:,2) ) > RAD_THRESH

         call calc_rhoch( veg, c1, rhoch )

         ! Canopy reflection of diffuse radiation for black leaves:
         do ictr=1,nrb
           rad%rhocdf(:,ictr) = rhoch(:,ictr) * ( gauss_w(1)*xk(:,1)/( xk(:,1) + rad%extkd(:) ) +  &
                                            gauss_w(2)*xk(:,2)/( xk(:,2) + rad%extkd(:) ) +  &
                                            gauss_w(3)*xk(:,3)/( xk(:,3) + rad%extkd(:) )  )
         enddo
         
         ! MJT bug fix for CCAM -------------------------------------------------------------------
         !if( .NOT. cable_runtime%um) then
         !   ! Define beam fraction, fbeam:
         !   rad%fbeam(:,1) = spitter(met%doy, met%coszen, met%fsd(:,1))
         !   rad%fbeam(:,2) = spitter(met%doy, met%coszen, met%fsd(:,2))
         !   ! coszen is set during met data read in.
         !   WHERE (met%coszen <1.0e-2)
         !      rad%fbeam(:,1) = 0.0
         !      rad%fbeam(:,2) = 0.0
         !   END WHERE
         !endif
         ! ----------------------------------------------------------------------------------------

   WHERE (canopy%vlaiw > LAI_THRESH)    ! In gridcells where vegetation exists....
      ! SW beam extinction coefficient ("black" leaves, extinction neglects
      ! leaf SW transmittance and reflectance):
      rad%extkb = xphi1 / met%coszen + xphi2
   ELSEWHERE ! i.e. bare soil
      rad%extkb = 0.5
   END WHERE
   
   WHERE ( abs(rad%extkb - rad%extkd)  < 0.001 )
      rad%extkb = rad%extkd + 0.001
   END WHERE
   
   WHERE(rad%fbeam(:,1) < RAD_THRESH )
      rad%extkb=30.0         ! keep cexpkbm within real*4 range (BP jul2010)
   END WHERE
   
      return
   END SUBROUTINE init_radiation



   SUBROUTINE radiation(bal, soil, ssoil, veg, air, met, rad, canopy, dels)
      use cable_common_module, only : cable_runtime, cable_user
      use cable_diag_module, only : cable_stat
      implicit none
      TYPE (soil_parameter_type),INTENT(IN)               :: soil
      TYPE (soil_snow_type),INTENT(INOUT) :: ssoil
      TYPE (veg_parameter_type),INTENT(IN)                :: veg
      TYPE (canopy_type),INTENT(IN)                :: canopy
      TYPE (air_type),INTENT(IN)          :: air
      TYPE (met_type),INTENT(INOUT)       :: met
      TYPE (radiation_type),INTENT(INOUT) :: rad
      TYPE (balances_type),INTENT(INOUT)  :: bal
      real(r_1), intent(in)               :: dels ! integration time setp (s)

  
      REAL(r_1), DIMENSION(mp)  :: cf1    ! (1.0 - rad%transb * cexpkdm) / (extkb + extkdm(:,b))
      REAL(r_1), DIMENSION(mp)  :: cf3    ! (1.0 - rad%transb * cexpkbm) / (extkb + extkbm(:,b))
      REAL(r_1), DIMENSION(mp)  :: cf2n   ! exp(-extkn * vlai) (nitrogen)
      REAL(r_1), DIMENSION(mp)  :: emair  ! air emissivity
      REAL(r_1), DIMENSION(mp)  :: flpwb  ! black-body long-wave radiation
      REAL(r_1), DIMENSION(mp)  :: flwv   ! vegetation long-wave radiation (isothermal)
      INTEGER(i_d)              :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
      LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation
      REAL(r_1), DIMENSION(mp)  :: xx1,tssp   ! vegetation long-wave radiation (isothermal)
      integer, save :: call_number =0
      REAL(r_2), DIMENSION(mp)  :: dummy2
      REAL(r_2), DIMENSION(mp)  :: dummy
      REAL(r_1) s1,s2,s3,step


      call_number = call_number + 1
  
      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('radiation')
    ! Define vegetation mask:
    mask = canopy%vlaiw > LAI_THRESH .AND. ( met%fsd(:,1)+met%fsd(:,2) ) > RAD_THRESH 
    ! Relative leaf nitrogen concentration within canopy:
    ! rad%extkn renamed veg%extkn
!jhan:Evaa made this change Aug-Jan after agreed to mv to veg%
!jh        where( veg%iveg == 2 ) rad%extkn = 0.01

    cf2n = EXP(-veg%extkn * canopy%vlaiw)
   
         rad%transd = 1.0
         WHERE (canopy%vlaiw > LAI_THRESH )    ! In gridcells where vegetation exists....
            ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
            ! leaf SW transmittance and reflectance);
            ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
            rad%transd = EXP(-rad%extkd * canopy%vlaiw)
         END WHERE

    ! Define fraction of SW beam tranmitted through canopy:
    dummy2 = -rad%extkb * canopy%vlaiw
    dummy = EXP(dummy2)
    rad%transb = REAL(dummy, r_1)
!    rad%transb = EXP(-rad%extkb * canopy%vlaiw)

    ! Define longwave from vegetation:
    !jhan: should we use tk or tvrad here    
   !jhanflpwb = sboltz * (met%tvrad) ** 4   ! YP Nov2009 (fix cold bias problem)
    flpwb = sboltz * (met%tk) ** 4
    flwv = emleaf * flpwb

    rad%flws = sboltz*emsoil* ssoil%tss **4
    ! Define air emissivity:
    emair = met%fld / flpwb
    
    rad%gradis = 0.0 ! initialise radiative conductance
    rad%qcan = 0.0   ! initialise radiation absorbed by canopy
    
    WHERE (canopy%vlaiw > LAI_THRESH )
       ! Define radiative conductance (Leuning et al, 1995), eq. D7:
       rad%gradis(:,1) = (4.0 * emleaf / (capp * air%rho)) * flpwb &
            & / (met%tk) * rad%extkd &
            & * ((1.0 - rad%transb * rad%transd) / (rad%extkb + rad%extkd) &
            & + (rad%transd - rad%transb) / (rad%extkb - rad%extkd))
       rad%gradis(:,2) = (8.0*emleaf/(capp*air%rho)) * flpwb / (met%tk) &
            & * rad%extkd * (1.0 - rad%transd) / rad%extkd - rad%gradis(:,1)
       ! Longwave radiation absorbed by sunlit canopy fraction:
       rad%qcan(:,1,3) = (rad%flws-flwv) *rad%extkd * (rad%transd - rad%transb) &
            & / (rad%extkb - rad%extkd) &
            & + (emair-emleaf) * rad%extkd * flpwb * (1.0-rad%transd*rad%transb) &
            & / ( rad%extkb + rad%extkd) 
       ! Longwave radiation absorbed by shaded canopy fraction:
       rad%qcan(:,2,3) = (1.0 - rad%transd) * &
            & (rad%flws + met%fld - 2.0 * flwv) - rad%qcan(:,1,3)
    END WHERE
    ! Convert radiative conductance from m/s to mol/m2/s:
    rad%gradis=SPREAD(air%cmolar, 2, mf)*rad%gradis
    rad%gradis = MAX(1.0e-3_r_2,rad%gradis)

         ! Update extinction coefficients and fractional transmittance for 
         ! leaf transmittance and reflection (ie. NOT black leaves):
         ! Define qcan for short wave (par, nir) for sunlit leaf:
         !jhan:prev. offline//Mk3l used met%fsd(:), hence assumed 1/2 met%fsd in each "b" calc here.
         !UM recieves met%fsd(:,b) forcing. assumed for offline that USED met%fsd(:,b) = 1/2* INPUT met%fsd
         DO b = 1, 2 ! 1 = visible, 2 = nir radiaition
            WHERE (mask) ! i.e. vegetation and sunlight are present
               cf1 = (1.0 - rad%transb * rad%cexpkdm(:,b)) / (rad%extkb + rad%extkdm(:,b))
               cf3 = (1.0 - rad%transb * rad%cexpkbm(:,b)) / (rad%extkb + rad%extkbm(:,b))
               rad%qcan(:,1,b) = met%fsd(:,b) * ( &          ! scale to real sunlit flux
                     (1.0-rad%fbeam(:,b))*(1.0-rad%reffdf(:,b))*rad%extkdm(:,b)*cf1 &
                     + rad%fbeam(:,b)*(1.0-rad%reffbm(:,b))*rad%extkbm(:,b)*cf3 &
                     + rad%fbeam(:,b)*(1.0-veg%taul(:,b)-veg%refl(:,b))*rad%extkb &
                     * ((1-rad%transb)/rad%extkb - (1-rad%transb**2)/(rad%extkb+rad%extkb)))
               ! Define qcan for short wave (par, nir) for shaded leaf:
               rad%qcan(:,2,b) = met%fsd(:,b) * ( &          ! scale to real shaded flux
                     (1.0-rad%fbeam(:,b))*(1.0-rad%reffdf(:,b))*rad%extkdm(:,b)* &
                     ((1.0 - rad%cexpkdm(:,b)) / rad%extkdm(:,b) - cf1) &
                     + rad%fbeam(:,b)*(1.-rad%reffbm(:,b))*rad%extkbm(:,b) &
                     * ((1.0 - rad%cexpkbm(:,b)) / rad%extkbm(:,b) - cf3) &
                     - rad%fbeam(:,b)*(1.0-veg%taul(:,b)-veg%refl(:,b))*rad%extkb &
                     * ((1-rad%transb)/rad%extkb - (1-rad%transb**2)/(rad%extkb+rad%extkb)))
       END WHERE
    END DO
    
    rad%qssabs = 0.
    
         WHERE (mask) ! i.e. vegetation and sunlight are present
            ! Calculate shortwave radiation absorbed by soil:
            ! (av. of transmitted NIR and PAR through canopy)*SWdown
            rad%qssabs = met%fsd(:,1) * (                                                 &
               rad%fbeam(:,1)*(1.-rad%reffbm(:,1))*EXP(-rad%extkbm(:,1)*canopy%vlaiw) +      &
               (1.-rad%fbeam(:,1))*(1.-rad%reffdf(:,1))*EXP(-rad%extkdm(:,1)*canopy%vlaiw) ) &
              + met%fsd(:,2)*(                                                              &
                rad%fbeam(:,2)*(1.-rad%reffbm(:,2))*rad%cexpkbm(:,2) +                    &
                (1.-rad%fbeam(:,2))*(1.-rad%reffdf(:,2))*rad%cexpkdm(:,2) )
     
            ! Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C:
            rad%scalex(:,1) = (1.0 - rad%transb * cf2n) / (rad%extkb + veg%extkn)
            ! Leaf area index of big leaf, sunlit, shaded, respectively:
            rad%fvlai(:,1) = (1.0 - rad%transb) / rad%extkb
            rad%fvlai(:,2) = canopy%vlaiw - rad%fvlai(:,1)
         ELSEWHERE ! i.e. either vegetation or sunlight are NOT present
            ! Shortwave absorbed by soil/snow surface:
            rad%qssabs = (1.0 - ssoil%albsoilsn(:,1))*met%fsd(:,1)+                         &
                         (1.0 - ssoil%albsoilsn(:,2))*met%fsd(:,2)
            rad%scalex(:,1) = 0.0
            rad%fvlai(:,1) = 0.0
            rad%fvlai(:,2) = canopy%vlaiw
         END WHERE

    rad%scalex(:,2) = (1.0 - cf2n) / veg%extkn - rad%scalex(:,1)
    ! Total energy absorbed by canopy:
    rad%rniso = SUM(rad%qcan, 3)
    
  END SUBROUTINE radiation

   subroutine calc_rhoch(veg,c1,rhoch) 
      use define_types
      use other_constants
      use cable_common_module, only : cable_runtime   
      implicit none
      type (veg_parameter_type), intent(inout) :: veg
      real, intent(inout), dimension(:,:) :: c1, rhoch
         
!jhan:UM uses rad%extkn instead of veg%extkn, which should be read from par. file anyway
!jhan:cahnge Mk3l to read veg%taul like UM 
         ! MJT bug fix ----------------------------------------------------------------------------
         !if( .NOT. cable_runtime%um) then
         !   veg%taul(:,1) = taul(1)
         !   veg%taul(:,2) = taul(2)
         !   veg%refl(:,1) = refl(1) 
         !   veg%refl(:,2) = refl(2) 
         !endif
         !-----------------------------------------------------------------------------------------                  
         c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
         c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
         c1(:,3) = 1.
          
         ! Canopy reflection black horiz leaves (eq. 6.19 in Goudriaan and van Laar, 1994):
         rhoch = (1.0 - c1) / (1.0 + c1)
      return
   end subroutine calc_rhoch 


  !-----------------------------------------------------------------
  ELEMENTAL FUNCTION sinbet(doy,xslat,hod) RESULT(z)
    ! calculate sin(bet), bet = elevation angle of sun
    ! calculations according to goudriaan & van laar 1994 p30
    REAL(r_1), INTENT(IN) :: doy    ! day of year
    REAL(r_1), INTENT(IN) :: xslat  ! latitude (degrees north)
    REAL(r_1), INTENT(IN) :: hod    ! hour of day
    REAL(r_1)             :: sindec ! sine of maximum declination
    REAL(r_1)             :: z      ! result
    sindec = -SIN(23.45 * pi180) * COS(two_pi * (doy + 10.0) / 365.0)
    z = MAX( SIN(pi180 * xslat) * sindec &
      & + COS(pi180 * xslat) * SQRT(1. - sindec * sindec) &
      & * COS(pi_c * (hod - 12.0) / 12.0), 1e-8 )
  END FUNCTION sinbet
  !-------------------------------------------------------------------
   FUNCTION spitter(doy, coszen, fsd) RESULT(fbeam)
    ! Calculate beam fraction
    ! See spitters et al. 1986, agric. for meteorol., 38:217-229
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: doy ! day of year
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: coszen ! cos(zenith angle of sun)
    REAL(r_1), DIMENSION(mp), INTENT(IN) :: fsd ! short wave down (positive) w/m^2
    REAL(r_1), DIMENSION(mp) :: fbeam   ! beam fraction (result)
    REAL(r_1), PARAMETER :: solcon = 1370.0
    REAL(r_1), DIMENSION(mp) :: tmpr !                                      
    REAL(r_1), DIMENSION(mp) :: tmpk !                                
    REAL(r_1), DIMENSION(mp) :: tmprat !                               
    fbeam = 0.0
    tmpr = 0.847 + coszen * (1.04 * coszen - 1.61)
    tmpk = (1.47 - tmpr) / 1.66
    WHERE (coszen > 1.0e-10 .AND. fsd > 10.0)
       tmprat = fsd / (solcon * (1.0 + 0.033 * COS(two_pi * (doy-10.0) / 365.0)) * coszen)
    ELSEWHERE
      tmprat = 0.0
    END WHERE
    WHERE (tmprat > 0.22) fbeam = 6.4 * (tmprat - 0.22) ** 2
    WHERE (tmprat > 0.35) fbeam = MIN(1.66 * tmprat - 0.4728, 1.0)
    WHERE (tmprat > tmpk) fbeam = MAX(1.0 - tmpr, 0.0)
  END FUNCTION spitter
END MODULE radiation_module
