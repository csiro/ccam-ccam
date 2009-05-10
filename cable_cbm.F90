MODULE cbm_module
  USE canopy_module
  USE cab_albedo_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC cbm 
CONTAINS
  SUBROUTINE cbm(ktau, kend,  dels, air, bgc, canopy, met, bal, &
             rad, rough, soil, ssoil, sum_flux, veg, L_EXPLICIT, L_HYD )
!            ftlt1,ftlt2,fqwt1,fqwt2)
    USE carbon_module
    USE soil_snow_module
    USE define_types
    USE physical_constants
    USE roughness_module
    USE radiation_module
    USE air_module
    USE albedo_module

    INTEGER(i_d) :: k  
    INTEGER(i_d), INTENT(IN)		:: ktau ! integration step number
!    INTEGER(i_d), INTENT(IN)	       	:: kstart ! starting value of ktau
    INTEGER(i_d), INTENT(IN)	       	:: kend ! total # timesteps in run
!    INTEGER(i_d), INTENT(IN)	       	:: ktauyear
    REAL(r_1), INTENT(IN)		:: dels ! time setp size (s)
    LOGICAL, INTENT(IN)                 :: L_EXPLICIT 
    LOGICAL, INTENT(IN)                 :: L_HYD 
    TYPE (air_type), INTENT(INOUT)	:: air
    TYPE (bgc_pool_type), INTENT(INOUT)	:: bgc	
    TYPE (canopy_type), INTENT(INOUT)	:: canopy
    TYPE (met_type), INTENT(INOUT) 	:: met
    TYPE (balances_type), INTENT(INOUT) 	:: bal
    TYPE (radiation_type), INTENT(INOUT) 	:: rad
    TYPE (roughness_type), INTENT(INOUT) 	:: rough
    TYPE (soil_parameter_type), INTENT(INOUT)	:: soil	
    TYPE (soil_snow_type), INTENT(INOUT)	:: ssoil
    TYPE (sum_flux_type), INTENT(INOUT)	:: sum_flux
    TYPE (veg_parameter_type), INTENT(INOUT)	:: veg	
!    REAL(r_1), INTENT(IN)		:: ftlt1 
!    REAL(r_1), INTENT(IN)		:: ftlt2 
!    REAL(r_1), INTENT(IN)		:: fqwt1 
!    REAL(r_1), INTENT(IN)		:: fqwt2 
    LOGICAL                 :: L_RADUM

    veg%meth = 1
    L_RADUM = .false.

!    print *,' cable_cbm is used',ms,ktau,veg%vlai,ssoil%wb,ssoil%tgg

!    IF(  L_EXPLICIT .and. .NOT. L_HYD ) CALL ruff_resist(veg, rough, ssoil)
    IF( .NOT. L_HYD ) CALL ruff_resist(veg, rough, ssoil)

!    print *,'bef define air'
    IF( .NOT. L_HYD ) CALL define_air (met, air)

!    print *,'bef init_rad'
!    IF(  L_EXPLICIT .and. .NOT. L_HYD ) CALL init_radiation(rad,veg)  
    IF( .NOT. L_HYD ) CALL init_radiation(met,rad,veg) ! need to be called at every dt

    if( ktau .lt. 2 .and. .NOT. L_HYD) CALL albedo(ktau,ssoil, veg, air, met, rad)

!    print *,'bef cab_alb'
    IF( L_EXPLICIT .and. .NOT. L_HYD ) &   ! call albedo once only
    CALL cab_albedo(ktau, dels, ssoil, veg, air, met, rad, soil, L_RADUM)

!    print *,'befcan',canopy%zetar(1,:)
    ! Calculate canopy variables:
!    IF(  L_EXPLICIT .and. .NOT. L_HYD ) CALL define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg,bgc,canopy)
     IF(.NOT. L_HYD ) CALL define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg,bgc,canopy,L_EXPLICIT)

!      print 101,ktau,rad%albedo(1,1),rad%albedo(2,1),rad%albedo(1,2),rad%albedo(2,2), met%coszen(1),met%fsd(1), &
!     IF(  L_EXPLICIT ) print 101,ktau,rad%albedo(1,1),rad%albedo(2,1),rad%albedo(1,2),rad%albedo(2,2), met%coszen(1),met%fsd(1), &
! rad%qcan(1,1,1)+rad%qcan(1,2,1)+rad%qcan(1,1,2)+rad%qcan(1,2,2)+rad%qssabs(1), &
! rad%qcan(2,1,1)+rad%qcan(2,2,1)+rad%qcan(2,1,2)+rad%qcan(2,2,2)+rad%qssabs(2), &
! ssoil%snowd
!101 format(1x,'ALBEDO',i5,4f6.3,1x,f7.4,f7.0,1x,2f7.0,1x,2f7.2)

#if defined(SCMA)
    IF(L_EXPLICIT .and. .NOT. L_HYD ) print 101,ktau,rad%albedo(1,1),rad%albedo(2,1),rad%albedo(1,2),rad%albedo(2,2),&
                             rad%reffdf(1,1),rad%reffdf(2,1),rad%reffdf(1,2),rad%reffdf(2,2),     &
                             rad%reffbm(1,1),rad%reffbm(2,1),rad%reffbm(1,2),rad%reffbm(2,2),rad%fbeam
     IF(L_EXPLICIT .and. .NOT. L_HYD ) print 209,ktau,                                    &
   (rad%qcan(1,1,1)+rad%qcan(1,2,1)+rad%qcan(1,1,2)+rad%qcan(1,2,2)+rad%qssabs(1)), &
   (rad%qcan(2,1,1)+rad%qcan(2,2,1)+rad%qcan(2,1,2)+rad%qcan(2,2,2)+rad%qssabs(2)), &                  
    met%fsd

101 format(1x'alb_cab',i5,18f6.3)
209 format(1x'rad%qcan',i5,8f7.0)

     IF( L_EXPLICIT .and. .NOT. L_HYD ) print 36,ktau,niter,canopy%zetar(1,:),&
          canopy%zetar(2,:),met%tk(1),canopy%tv(1),met%ua(1),canopy%us
36     format(1x,'ZETAR1',i5,i2,4f7.2,3x,4f7.2,1x,2f6.1,1x,f5.1,2f7.4)
#endif


    ! Calculate soil and snow variables:
!     print *,'canopy%ga',canopy%ga,canopy%fev,canopy%fhv,canopy%fevc,canopy%fevw
!     print 133,ktau,canopy%fes,canopy%fhs,canopy%fev,canopy%fevc,canopy%fevw,canopy%fhv
133 format(1x,'fluxess bef soilsnow',i5,12f6.0)
!     print *,'cbm L_EXPLICIT vlai',L_EXPLICIT,veg%vlai

!     IF( .not. L_EXPLICIT ) &
     
!     IF( L_HYD ) &
     IF( L_HYD ) THEN
          ssoil%owetfac = ssoil%wetfac
          CALL soil_snow(dels, ktau, soil, ssoil, canopy, met, bal)
     ENDIF

!     print 131,ktau,ssoil%tgg
!131  format(1x,' after soilsnow ssoil%tgg',i5,10f7.2)
!     print 132,ktau,ssoil%wb
!132   format(1x,'ssoil%wb',i5,10f6.3)
!     print 134,ktau,canopy%tv,met%tk
!134   format(1x,'tv,tk',i5,10f6.1)

    !	need to adjust fe after soilsnow
    canopy%fev	= canopy%fevc + canopy%fevw
    ! Calculate total latent heat flux:
    canopy%fe = canopy%fev + canopy%fes
    ! Calculate net radiation absorbed by soil + veg
    canopy%rnet = canopy%fns + canopy%fnv
!    print 201,ktau,rad%rniso,canopy%fe,canopy%fh,canopy%fev,canopy%fevc,canopy%fevw,canopy%fhv
201 format('SEB1',i5,15f7.1)
!    print 202,ktau,met%tk(1),met%fsd(1),met%fld(1),rad%qssabs,canopy%fes,canopy%fhs,canopy%ga
!202  format('SEB2',i5,f6.1,f6.0,f5.0,15f6.0)
!    print 203,ktau,canopy%fes,canopy%fhs,canopy%fev,canopy%fevc,canopy%fevw,canopy%fhv
!203  format(1x,'fluxess',i5,12f6.0)
!    print 204,ktau,SUM(rad%qcan(:,:,1),2),SUM(rad%qcan(:,:,2),2),rad%qssabs, &
!    met%fld,sboltz*emleaf*canopy%tv**4*(1-rad%transd),rad%flws*rad%transd,  &
!    canopy%fev,canopy%fes,canopy%fh ,canopy%ga
!204  format(1x,'SEB3',i5,22f6.0)
!    print 205,ktau,SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs+ &
!    met%fld-sboltz*emleaf*canopy%tv**4*(1-rad%transd)-rad%flws*rad%transd-  &
!    canopy%fev-canopy%fes-canopy%fh -canopy%ga
!205  format(1x,'SEB4',i5,22f6.0)
!    print 105,ktau,canopy%fe,canopy%fh,canopy%rnet
!105  format(1x,'grid fluxes',i5,12f6.0)
    ! Calculate radiative/skin temperature:
!    rad%trad = ( (1.-rad%transd)*canopy%tv**4 + rad%transd * ssoil%tss**4 )**0.25
!    print *,'trad',ktau,rad%trad,rad%transd,veg%vlai,canopy%fev,canopy%fevc
    
!    sum_flux%sumpn = sum_flux%sumpn+canopy%fpn*dels
!    sum_flux%sumrp = sum_flux%sumrp+canopy%frp*dels
!    sum_flux%sumrpw = sum_flux%sumrpw+canopy%frpw*dels
!    sum_flux%sumrpr = sum_flux%sumrpr+canopy%frpr*dels
!    sum_flux%sumrd = sum_flux%sumrd+canopy%frday*dels
!    sum_flux%dsumpn = sum_flux%dsumpn+canopy%fpn*dels
!    sum_flux%dsumrp = sum_flux%dsumrp+canopy%frp*dels
!    sum_flux%dsumrd = sum_flux%dsumrd+canopy%frday*dels
!    rad%flws = sboltz*emsoil* ssoil%tss **4
!    
!    CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)
!
!    CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)
!
!     print *,'bef call to coeftest',canopy%us
!    IF( L_EXPLICIT ) &
!    CALL coeftest(rad,rough,air,met,dels,ssoil,soil,canopy,ftlt1,ftlt2,fqwt1,fqwt2)
!
!    ! canopy%frs set in soilcarb
!    sum_flux%sumrs = sum_flux%sumrs+canopy%frs*dels
!    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
  
END SUBROUTINE cbm

END MODULE cbm_module
