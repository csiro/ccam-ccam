MODULE roughness_module
  USE physical_constants
  USE define_types
  USE define_dimensions
  IMPLICIT NONE
  PRIVATE
  PUBLIC ruff_resist
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE ruff_resist(veg, rough, ssoil, soil, met, canopy)
    ! m.r. raupach, 24-oct-92
    ! see: Raupach, 1992, BLM 60 375-395
    !      MRR notes "Simplified wind model for canopy", 23-oct-92
    !      MRR draft paper "Simplified expressions...", dec-92
    ! modified to include resistance calculations by Ray leuning 19 Jun 1998  
      use cable_common_module, only : cable_runtime, cable_user
      use cable_diag_module, only : cable_stat
    implicit none
    TYPE (veg_parameter_type), INTENT(INOUT)       :: veg
    TYPE (soil_snow_type), INTENT(IN)      :: ssoil
    TYPE (soil_parameter_type), INTENT(IN) :: soil  ! soil parameters
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    REAL(r_1), DIMENSION(mp) :: xx ! =ccd*LAI; working variable 
    REAL(r_1), DIMENSION(mp) :: dh ! d/h where d is zero-plane displacement
    REAL(r_1), DIMENSION(mp) :: hmax ! maximum height of canopy from
                                     ! tiles belonging to the same grid
 
      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('ruff_resist')

    ! Set canopy height above snow level:
!jhan:Eva has changed 0.01 to 0.001
    rough%hruff = MAX(0.01,veg%hc-1.2*ssoil%snowd/MAX(ssoil%ssdnn,100.)) 
    ! maximum height of canopy from tiles belonging to the same grid
    hmax = rough%hruff_grmx
    ! LAI decreases due to snow and vegetation fraction:
    !jh:use this in UM
!jhan:Eva still using veg%
    canopy%vlaiw = veg%vlai * rough%hruff/MAX(0.01,veg%hc)
    canopy%rghlai = canopy%vlaiw
    where(ssoil%snowd.lt.0.001.and.veg%iveg.ne.1) canopy%rghlai = min( 3.,canopy%vlaiw)
!            (veg%iveg.ne.1.and.veg%iveg.ne.6.and.veg%iveg.ne.9)) &

    ! Roughness length of bare soil (m):
    rough%z0soil = 1.e-6
    !rough%z0soil = MIN(0.001,MAX(0.0011*exp(-canopy%vlaiw),1.e-6))
    !rough%z0soilsn = MAX(MIN(-7.5e-6*(0.01*MIN(ssoil%snowd,20.))+ &
    rough%z0soilsn = max(rough%z0soil-0.5e-7*min(ssoil%snowd,20.),0.1e-7)

!jhan:Eva adds this cond.
    !jh:WHERE( soil%isoilm == 9 ) rough%z0soilsn = 1.e-5
    
    WHERE (canopy%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn) ! BARE SOIL SURFACE
       rough%z0m = rough%z0soilsn
       rough%hruff = 0.0
       rough%rt0us = 0.0  
       rough%disp = 0.0
    ! Reference height zref is height above the displacement height
       !rough%zref_uv = max(3.5,rough%za_uv - rough%disp + hmax -  &
       !                                      min(1.,ssoil%snowd/max(ssoil%ssdnn,100.)))
       !rough%zref_tq = max(3.5,rough%za_tq - rough%disp + hmax - &
       !                                      min(1.,ssoil%snowd/max(ssoil%ssdnn,100.)))
!       rough%zref_uv = max(3.5,rough%za_uv + rough%hruff)
       rough%zref_uv = max(3.5,rough%za_uv )
       rough%zref_tq = max(3.5,rough%za_tq )

       rough%zruffs = 0.0
       rough%rt1usa = 0.0 
       rough%rt1usb = 0.0
       ! Friction velocity/windspeed at canopy height
       ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
       ! (usuhm set in physical_constants module):
       rough%usuh = MIN(SQRT(csd+crd*(canopy%vlaiw*0.5)), usuhm)
       ! xx is ccd (see physical_constants) by LAI
       xx = SQRT(ccd*MAX((canopy%vlaiw*0.5),0.0005))
       ! Displacement height/canopy height:
       ! eq.8 Raupach 1994, BLM, vol 71, p211-216
       dh = 1.0 - (1.0 - EXP(-xx))/xx
       ! Extinction coefficient for wind profile in canopy:
       ! eq. 3.14, SCAM manual (CSIRO tech report 132)
       rough%coexp = rough%usuh / (vonk*ccw_c*(1.0 - dh))
    ELSEWHERE ! VEGETATED SURFACE
       ! Friction velocity/windspeed at canopy height
       ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
       ! (usuhm set in physical_constants module):
!jhan:Eva uses         
        !rough%usuh = MIN(SQRT(csd+crd*(canopy%vlaiw*0.5)), usuhm)
        rough%usuh = MIN(SQRT(csd+crd*(canopy%rghlai*0.5)), usuhm)
       ! xx is ccd (see physical_constants) by LAI:
!jhan:Eva uses         
       !xx = SQRT(ccd*MAX((canopy%vlaiw*0.5),0.0005))
       xx = SQRT(ccd*MAX((canopy%rghlai*0.5),0.0005))
       ! eq.8 Raupach 1994, BLM, vol 71, p211-216:
       dh = 1.0 - (1.0 - EXP(-xx))/xx
       ! Calculate zero-plane displacement:
       rough%disp = dh*rough%hruff
    ! Reference height zref is height above the displacement height
       !rough%zref_uv = max(3.5,rough%za_uv -  & 
       !                max(rough%disp,min(1.,ssoil%snowd/max(ssoil%ssdnn,100.))) + hmax)
       !rough%zref_tq = max(3.5,rough%za_tq -  &
       !                max(rough%disp,min(1.,ssoil%snowd/max(ssoil%ssdnn,100.))) + hmax)
       rough%zref_uv = max(3.5,rough%za_uv )
       rough%zref_tq = max(3.5,rough%za_tq )
!       rough%zref_uv = max(3.5,rough%za_uv + rough%hruff)
!       rough%zref_tq = max(3.5,rough%za_tq + rough%hruff)
       ! Calcualte roughness length:
       rough%z0m = ( (1.0 - dh) * EXP( LOG(ccw_c) - 1. + 1./ccw_c &
            & - vonk/rough%usuh ) ) * rough%hruff
       ! find coexp: see notes "simplified wind model ..." eq 34a
       ! Extinction coefficient for wind profile in canopy:
       ! eq. 3.14, SCAM manual (CSIRO tech report 132)
       rough%coexp = rough%usuh / (vonk*ccw_c*(1.0 - dh))
!jhan:Eva uses         
       !jh rough%term2  = EXP(2*csw*min(veg%vlaiw,3.)*(1-rough%disp/rough%hruff))
       !jh rough%term3  = a33**2*ctl*2*csw*min(veg%vlaiw,3.)

       rough%term2  = EXP(2*csw*canopy%rghlai*(1-rough%disp/rough%hruff))
       rough%term3  = a33**2*ctl*2*csw*canopy%rghlai
       rough%term5  = MAX((2./3.)*rough%hruff/rough%disp, 1.0)
       rough%term6 =  exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
       ! eq. 3.54, SCAM manual (CSIRO tech report 132)
       rough%rt0us  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!jhan:Eva uses         
            !jh + (1-zdlin))*(EXP(2*csw*min(veg%vlaiw,3.)) - rough%term2)/rough%term3  ! &
            + (1-zdlin))*(EXP(2*csw*canopy%rghlai) - rough%term2)/rough%term3  ! &
!              / rough%term6
       !        rt1 = turbulent resistance from canopy (z1 = disp) to
       !        reference level zref (from disp as origin). Normalisation:
       !        rt1us = us*rt1 = rt1usa + rt1usb + rt1usc
       !        with term a = resistance from disp to hruf
       !        term b = resistance from hruf to zruffs (or zref if zref<zruffs)
       !        term c = resistance from zruffs to zref (= 0 if zref<zruffs)
       !        where zruffs = SCALAR roughness sublayer depth (ground=origin)
       !        xzrufs = xdisp + xhruf*a33**2*ctl/vonk
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
       rough%zruffs = rough%disp + rough%hruff*a33**2*ctl/vonk/rough%term5
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
       rough%rt1usa = rough%term5*(rough%term2 - 1.0)/rough%term3
       rough%rt1usb = rough%term5*(MIN(rough%zref_tq+rough%disp,rough%zruffs) - rough%hruff)/ &
            (a33**2*ctl*rough%hruff)
       rough%rt1usb = MAX(rough%rt1usb,0.0)       ! in case zrufs < rough%hruff
    END WHERE
  END SUBROUTINE ruff_resist
END MODULE roughness_module
