MODULE roughness_module
  USE physical_constants
  USE define_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC ruff_resist
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE ruff_resist
!  SUBROUTINE ruff_resist(veg, rough, ssoil,canopy)
    ! m.r. raupach, 24-oct-92
    ! see: Raupach, 1992, BLM 60 375-395
    !      MRR notes "Simplified wind model for canopy", 23-oct-92
    !      MRR draft paper "Simplified expressions...", dec-92
    ! modified to include resistance calculations by Ray leuning 19 Jun 1998  
!    TYPE (veg_parameter_type), INTENT(INOUT)	   :: veg
!    TYPE (soil_snow_type), INTENT(IN)	   :: ssoil
!    TYPE (roughness_type), INTENT(INOUT) :: rough
!    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    REAL(r_1), DIMENSION(mp)	   :: xx ! =ccd*LAI; working variable
    REAL(r_1), DIMENSION(mp)	   :: dh ! d/h where d is zero-plane displacement
    REAL(r_1), DIMENSION(mp)	     :: hmax ! maximum height of canopy from
                                           ! tiles belonging to the same grid
    REAL(r_1), DIMENSION(mp)       :: rghlai
    ! Set canopy height above snow level:
!    print *,'ruff_resist 1',veg%hc,ssoil%snowd,ssoil%ssdnn
    rough%hruff= MAX(0.01,veg%hc-1.2*ssoil%snowd/max(ssoil%ssdnn,100.)) 
    ! maximum height of canopy from tiles belonging to the same grid
    hmax = rough%hruff_grmx
    ! LAI decreases due to snow and vegetation fraction:
    canopy%vlaiw = veg%vlai * rough%hruff/MAX(0.01,veg%hc)
    rghlai = canopy%vlaiw
    where(ssoil%snowd.lt.0.001) rghlai = min( 3.,canopy%vlaiw)
    ! Roughness length of bare soil (m):
    rough%z0soil = 1.e-6
    where(ssoil%snowd.lt.0.001) rough%z0soil = 0.76e-6 + 1.e-4*exp( -(6.-min(canopy%vlaiw, 6.0)) )
    
    rough%z0soilsn = max(rough%z0soil-0.5e-7*min(ssoil%snowd,20.),0.1e-7)
    WHERE( soil%isoilm == 9 ) rough%z0soilsn = 1.e-5    
    
    WHERE (canopy%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn) ! i.e. BARE SOIL SURFACE
       rough%z0m = rough%z0soilsn
       rough%hruff = 0.0
       rough%rt0us = 0.0  
       rough%disp  = 0.0
       ! Reference height 
       rough%zref_uv = max(3.5,rough%za_uv - rough%disp + hmax)
       rough%zref_tq = max(3.5,rough%za_tq - rough%disp + hmax)
!       rough%zref_uv = max(3.5,rough%za_uv + rough%hruff)
!       rough%zref_tq = max(3.5,rough%za_tq + rough%hruff)
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
       rough%usuh = MIN(SQRT(csd+crd*(rghlai*0.5)), usuhm)
       ! xx is ccd (see physical_constants) by LAI:
       xx = SQRT(ccd*MAX((rghlai*0.5),0.0005))
       ! eq.8 Raupach 1994, BLM, vol 71, p211-216:
       dh = 1.0 - (1.0 - EXP(-xx))/xx
       ! Calculate zero-plane displacement:
       rough%disp = dh*rough%hruff
       ! Reference height zref is height above the displacement height
       rough%zref_uv = max(3.5,rough%za_uv - rough%disp + hmax)
       rough%zref_tq = max(3.5,rough%za_tq - rough%disp + hmax)
!       rough%zref_uv = max(3.5,rough%za_uv + rough%hruff)
!       rough%zref_tq = max(3.5,rough%za_tq + rough%hruff)
       ! Calcualte roughness length:
       rough%z0m = ((1. - dh)*EXP(LOG(ccw_c)-1. + 1./ccw_c - vonk/rough%usuh))*rough%hruff
       !	find coexp: see notes "simplified wind model ..." eq 34a
       ! Extinction coefficient for wind profile in canopy:
       ! eq. 3.14, SCAM manual (CSIRO tech report 132)
       rough%coexp = rough%usuh / (vonk*ccw_c*(1.0 - dh))
       !	rt0 = turbulent resistance from soil (z0 = 0) to canopy
       !	(z1 = zero-plane displacement), normalized as rt0us=rt0*us
       !	nb: rough%term4 added 13-sep-95 to make TL proportional to z near ground.
       !	nb: rough%term5 added 03-oct-96 to account for sparse canopies. Constant
       !	length scale ctl*hruf replaced by ctl*(3/2)*disp, taking effect
       !	when disp<(2/3)*hruf, or total LAI < 1.11. Otherwise, rough%term5=1.
       rough%term2  = EXP(2*csw*rghlai*(1-rough%disp/rough%hruff))
       rough%term3  = a33**2*ctl*2*csw*rghlai
       rough%term5  = MAX(2./3.*rough%hruff/rough%disp, 1.)
       rough%term6  = exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
       ! eq. 3.54, SCAM manual (CSIRO tech report 132)
       rough%rt0us  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) + (1.-zdlin))* &
                      (EXP(2*csw*rghlai) - rough%term2)/rough%term3  
!              / rough%term6
       !	rt1 = turbulent resistance from canopy (z1 = disp) to
       !	reference level zref (from disp as origin). Normalisation:
       !	rt1us = us*rt1 = rt1usa + rt1usb + rt1usc
       !	with term a = resistance from disp to hruf
       !	term b = resistance from hruf to zruffs (or zref if zref<zruffs)
       !	term c = resistance from zruffs to zref (= 0 if zref<zruffs)
       !	where zruffs = SCALAR roughness sublayer depth (ground=origin)
       !	xzrufs = xdisp + xhruf*a33**2*ctl/vonk
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
       rough%zruffs = rough%disp + rough%hruff*a33**2*ctl/vonk/rough%term5
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
       rough%rt1usa = rough%term5*(rough%term2 - 1.)/rough%term3
       rough%rt1usb = rough%term5*(MIN(rough%zref_tq+rough%disp,rough%zruffs) - rough%hruff)/ &
                      (a33**2*ctl*rough%hruff)
       rough%rt1usb = MAX(rough%rt1usb,0.)       ! in case zrufs < rough%hruff
    END WHERE
!   print *,'z0m',rough%z0m,canopy%vlaiw,rough%hruff,hmax
!   print *,'ruff_resEND',rough%rt1usa,rough%rt1usb,rough%disp,rough%zruffs
  END SUBROUTINE ruff_resist

END MODULE roughness_module
