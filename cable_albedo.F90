
MODULE albedo_module
   IMPLICIT NONE
   private
   PUBLIC surface_albedo
   CONTAINS

   SUBROUTINE surface_albedo(ssoil, veg, met, rad, soil, canopy)
      use cable_common_module
      use define_types
      use define_dimensions
      use other_constants, only : LAI_THRESH, RAD_THRESH 
      use cable_diag_module, only : cable_stat
      implicit none
      TYPE (soil_snow_type),INTENT(INOUT) :: ssoil
      TYPE (veg_parameter_type),INTENT(INout):: veg
      TYPE (canopy_type),INTENT(IN)          :: canopy
      TYPE (met_type),INTENT(INOUT)       :: met
      TYPE (radiation_type),INTENT(INOUT) :: rad
      TYPE(soil_parameter_type), INTENT(INOUT) :: soil   

      INTEGER(i_d)            :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
      LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation
      real, dimension(:,:), allocatable, save :: c1, rhoch
      REAL(r_2), DIMENSION(mp)  :: dummy2
      REAL(r_2), DIMENSION(mp)  :: dummy

      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('surface_albedo')
      
      if (.not. allocated(c1)) &
         allocate( c1(mp,nrb), rhoch(mp,nrb) )

      if(cable_runtime%um) &
         call surface_albedosn(ssoil, veg, met, soil)

      !jhan:Eva has added this
      WHERE (soil%isoilm == 9)          ! use dry snow albedo
       ssoil%albsoilsn(:,2) = 0.82
        ssoil%albsoilsn(:,1) = 0.82
      END WHERE
      
      rad%cexpkbm = 0.0
      rad%extkbm  = 0.0
      rad%rhocbm      = 0.0
 
      ! Initialise effective conopy beam reflectance:
      rad%reffbm = ssoil%albsoilsn
      rad%reffdf = ssoil%albsoilsn
      rad%albedo = ssoil%albsoilsn
 
      ! Define vegetation mask:
      mask = canopy%vlaiw > LAI_THRESH .AND. ( met%fsd(:,1) + met%fsd(:,2) ) > RAD_THRESH     

      call calc_rhoch( veg, c1, rhoch )

      ! Update extinction coefficients and fractional transmittance for 
      ! leaf transmittance and reflection (ie. NOT black leaves):
      !---1 = visible, 2 = nir radiaition
      DO b = 1, 2        
         
         rad%extkdm(:,b) = rad%extkd * c1(:,b)
         !--Define canopy diffuse transmittance (fraction):
         rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)
 
         !---Calculate effective diffuse reflectance (fraction):
         where( canopy%vlaiw > 1e-2 ) 
            rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssoil%albsoilsn(:,b) &
                  - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2
         END WHERE
         
         !---where vegetated and sunlit 
         WHERE (mask)                
            rad%extkbm(:,b) = rad%extkb * c1(:,b)
            ! Canopy reflection (6.21) beam:
            !jhan:BP uses local var rhocbm here not rad%
            rad%rhocbm(:,b) = 2.*rad%extkb/(rad%extkb+rad%extkd)*rhoch(:,b)
            ! Canopy beam transmittance (fraction):
            dummy2 = -rad%extkbm(:,b)*canopy%vlaiw
            dummy  = EXP(dummy2)
            rad%cexpkbm(:,b) = REAL(dummy, r_1)
!            rad%cexpkbm(:,b) = EXP(-rad%extkbm(:,b)*canopy%vlaiw)
            ! Calculate effective beam reflectance (fraction):
            rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssoil%albsoilsn(:,b) &
                  - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2
         END WHERE
 
         ! Define albedo:
         where( canopy%vlaiw> LAI_THRESH ) 
            rad%albedo(:,b) = (1.-rad%fbeam(:,b))*rad%reffdf(:,b) + &
                  rad%fbeam(:,b)*rad%reffbm(:,b)
         END WHERE
          
      END DO

      return
   END SUBROUTINE surface_albedo 

   !jhan:subr was reintroduced here to temporarily resolve issue when 
   !creating libcable.a 
   subroutine calc_rhoch(veg,c1,rhoch) 
      use define_types
      use other_constants
      implicit none
      type (veg_parameter_type), intent(inout) :: veg
      real, intent(inout), dimension(:,:) :: c1, rhoch

!jhan:UM uses rad%extkn instead of veg%extkn, which should be read from par. file anyway
!jhan:cahnge Mk3l to read veg%taul like UM 
!#ifndef ONLINE_UM
!         veg%taul(:,1) = taul(1)
!         veg%taul(:,2) = taul(2)
!         veg%refl(:,1) = refl(1) 
!         veg%refl(:,2) = refl(2) 
!#endif                  
         c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
         c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
         c1(:,3) = 1.
          
         ! Canopy reflection black horiz leaves (eq. 6.19 in Goudriaan and van Laar, 1994):
         rhoch = (1.0 - c1) / (1.0 + c1)
      return
   end subroutine calc_rhoch 

   SUBROUTINE surface_albedosn(ssoil, veg, met, soil)
      use define_types
      use define_dimensions
      use cable_common_module
      implicit none
      TYPE (soil_snow_type),INTENT(INOUT) :: ssoil
      TYPE (veg_parameter_type),INTENT(INout):: veg
      TYPE (met_type),INTENT(INOUT)       :: met
      TYPE(soil_parameter_type), INTENT(INOUT) :: soil   
      REAL(r_1), DIMENSION(mp)    :: alv ! Snow albedo for visible
      REAL(r_1), DIMENSION(mp)    :: alir ! Snow albedo for near infra-red
      REAL(r_1), PARAMETER        :: alvo  = 0.95 ! albedo for vis. on a new snow
      REAL(r_1), PARAMETER        :: aliro = 0.70 ! albedo for near-infr. on a new snow
      REAL(r_1), DIMENSION(mp)    :: ar1 ! crystal growth  (-ve)
      REAL(r_1), DIMENSION(mp)    :: ar2 ! freezing of melt water
      REAL(r_1), DIMENSION(mp)    :: ar3
      REAL(r_1), DIMENSION(mp)    :: dnsnow ! new snow albedo
      REAL(r_1), DIMENSION(mp)    :: dtau
      REAL(r_1), DIMENSION(mp)    :: fage !age factor
      REAL(r_1), DIMENSION(mp)    :: fzenm
      REAL(r_1), DIMENSION(mp)    :: sfact
      REAL(r_1), DIMENSION(mp)    :: snr
      REAL(r_1), DIMENSION(mp)    :: snrat
      REAL(r_1), DIMENSION(mp)    :: talb ! snow albedo
      REAL(r_1), DIMENSION(mp)    :: tmp ! temporary value
      INTEGER(i_d)            :: k,i,j,l,l1,l2

!jhan:incurrent code in coupled model soil%albsoil(mp)    
    soil%albsoilf = soil%albsoil(:,1)
    where( veg%iveg == 16 )
      soil%albsoilf = -0.022*( min(275., max(260.,met%tk) ) - 260.) + 0.45
    end where
!jhan:Eva has added this
    where(ssoil%snowd > 1. .and. veg%iveg == 16 ) soil%albsoilf = 0.85
    sfact = 0.68
    WHERE (soil%albsoilf <= 0.14)
       sfact = 0.5
    ELSEWHERE (soil%albsoilf > 0.14 .and. soil%albsoilf <= 0.20)
       sfact = 0.62
    END WHERE
    ssoil%albsoilsn(:,2) = 2. * soil%albsoilf / (1. + sfact)
    ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
    snrat=0.
!jhan:Eva has added this
!.6 and .9 -> 0.
    alir =0.
    alv  =0.
!jhan:Eva has added this 0 ->1
    WHERE (ssoil%snowd > 1. .and. .NOT. cable_runtime%um_radiation) 
       dnsnow = min (1., .1 * max (0., ssoil%snowd - ssoil%osnowd ) ) ! new snow (cm H2O)
       !         Snow age depends on snow crystal growth, freezing of melt water,
       !         accumulation of dirt and amount of new snow.
       tmp = ssoil%isflag * ssoil%tggsn(:,1) + (1 - ssoil%isflag ) * ssoil%tgg(:,1)
       tmp = min (tmp, 273.15)
       ar1 = 5000. * (1. / 273.15 - 1. / tmp) ! crystal growth  (-ve)
       ar2 = 10. * ar1 ! freezing of melt water
!jhan:Eva has added this 300 -> 200
       snr = ssoil%snowd / max (ssoil%ssdnn, 200.)
       WHERE (soil%isoilm == 9)
          ar3 = .0000001
          dnsnow = max (dnsnow, .5) !increase refreshing of snow in Antarctic
  !jhan:Eva has added this
            dnsnow = 1.0
          snrat = 1.
       ELSEWHERE
          ! accumulation of dirt
  !jhan:Eva has added this
          !ar3 = .2
          ar3 = .1
          snrat = min (1., snr / (snr + .1) )    ! snow covered fraction of the grid
       END WHERE
       dtau = 1.e-6 * (exp(ar1) + exp(ar2) + ar3) * kwidth_gl 
       WHERE (ssoil%snowd <= 1.0)
          ssoil%snage = 0.
       ELSEWHERE
          ssoil%snage = max (0.,(ssoil%snage+dtau)*(1.-dnsnow))
       END WHERE
       fage = 1. - 1. / (1. + ssoil%snage ) !age factor

       tmp = max (.17365, met%coszen )
       fzenm = max(merge(0.0, (1. + 1./2.)/(1. + 2.*2.*tmp) - 1./2., tmp > 0.5), 0.)
       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
       WHERE (soil%isoilm == 9)             ! use dry snow albedo for pernament land ice
  !jhan:Eva has channged this - needs to be sym link
         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
  !jhan:Eva has channged this - needs to be sym link
         tmp = 0.75 * (1. - .5 * fage)
       END WHERE
       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo
    ENDWHERE        ! snowd > 0
   !   wnhen it is called from  cable_RADum no need to recalculate snage 
  !jhan:Eva has channged this - needs to be sym link
    WHERE (ssoil%snowd > 1 .and. cable_runtime%um_radiation )
       snr = ssoil%snowd / max (ssoil%ssdnn, 200.)
       WHERE (soil%isoilm == 9)
          snrat = 1.
       ELSEWHERE
          snrat = min (1., snr / (snr + .1) )
       END WHERE
       fage = 1. - 1. / (1. + ssoil%snage ) !age factor
       tmp = max (.17365, met%coszen )
       fzenm = max(merge(0.0, (1. + 1./2.)/(1. + 2.*2.*tmp) - 1./2., tmp > 0.5), 0.)
       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
       WHERE (soil%isoilm == 9)          ! use dry snow albedo
  !jhan:Eva has channged this - needs to be sym link
         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
  !jhan:Eva has channged this - needs to be sym link
         tmp = 0.75 * (1. - .5 * fage)
       END WHERE
       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo
    ENDWHERE        ! snowd > 0
    ssoil%albsoilsn(:,2) = min(aliro,(1. - snrat) * ssoil%albsoilsn(:,2) + snrat * alir)
    ssoil%albsoilsn(:,1) = min(alvo,(1. - snrat) * ssoil%albsoilsn(:,1) + snrat * alv)

   end subroutine surface_albedosn

END MODULE albedo_module
