MODULE canopy_module
  USE photosynthetic_constants
  USE radiation_module
  USE roughness_module
  USE air_module
  USE define_types
  USE physical_constants
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_canopy, sinbet
!  PUBLIC define_canopy, sinbet, coeftest
CONTAINS
  
  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg,bgc,canopy,L_EXPLICIT)
    TYPE (balances_type),INTENT(INOUT)  :: bal
    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)	:: air
    TYPE (met_type), INTENT(INOUT)	:: met
    REAL(r_1), INTENT(IN)		:: dels ! integration time setp (s)
    TYPE (soil_snow_type), INTENT(INOUT):: ssoil
    TYPE (bgc_pool_type),INTENT(IN)	:: bgc
    TYPE (soil_parameter_type), INTENT(INOUT)	:: soil
    TYPE (veg_parameter_type), INTENT(INOUT)	:: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    LOGICAL, INTENT(IN)                 :: L_EXPLICIT
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp,mf)	        :: abs_deltlf ! ABS(deltlf)
    REAL(r_1), DIMENSION(mp,mf,3)	:: ancj ! soln to quad eqn
    REAL(r_1), DIMENSION(mp,mf)		:: anx ! net photos. prev iteration
    REAL(r_1), DIMENSION(mp,mf)		:: an_y ! net photosynthesis soln
    REAL(r_1), DIMENSION(mp)		:: avgtrs !root weighted mean soil temperature
    REAL(r_1), DIMENSION(mp)		:: avgwrs !root weighted mean soil moisture
    REAL(r_1), DIMENSION(mp,mf)		:: ca2	 ! 2D CO2 concentration
    REAL(r_1), DIMENSION(mp)		:: cansat ! max canopy intercept. (mm)
    REAL(r_1), DIMENSION(mp,mf,3)	:: ci ! intercellular CO2 conc.
    REAL(r_1), PARAMETER		:: co2cp3=0.0 ! CO2 compensation pt C3
    REAL(r_1), DIMENSION(mp,mf,3)	:: coef0 ! CO2 comp. pt coeff 1
    REAL(r_1), DIMENSION(mp,mf,3)	:: coef1 ! " 2
    REAL(r_1), DIMENSION(mp,mf,3)	:: coef2 ! " 3
    REAL(r_1), DIMENSION(mp,mf)		:: conkct ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp,mf)		:: conkot ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp,mf)		:: csx ! leaf surface CO2 concentration
    REAL(r_1), DIMENSION(mp,mf,3)	:: cx  ! "d_{3}" in Wang and Leuning, 1998, appendix E
    REAL(r_1), DIMENSION(mp,mf)		:: da2 ! 2D sat vap pres deficit
    REAL(r_1), DIMENSION(mp,mf)		:: dva2 ! 2D in canopy sat vap pres deficit
    REAL(r_1), DIMENSION(mp,mf,3)	:: delcx ! discriminant  in quadratic in eq. E7 Wang and Leuning, 1998
    REAL(r_1), DIMENSION(mp,mf)		:: deltlf ! deltlfy of prev iter.
    REAL(r_1), DIMENSION(mp,mf)		:: deltlfy ! del temp successive iteration
    REAL(r_1), DIMENSION(mp)		:: dq ! sat spec hum diff.
    REAL(r_1), DIMENSION(mp,mf)		:: dsatdk2	! 2D dsatdk
    REAL(r_1), DIMENSION(mp,mf)		:: dsx ! leaf surface vpd
    REAL(r_1), DIMENSION(mp,mf)		:: ecx ! lat. hflux big leaf
    REAL(r_1), DIMENSION(mp,mf)		:: ejmax2 ! jmax of big leaf
    REAL(r_1), DIMENSION(mp,mf)		:: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_1), DIMENSION(mp,mf)		:: ecy ! lat heat fl dry big leaf
    REAL(r_1), DIMENSION(mp,mf)		:: frac42	! 2D frac4
    REAL(r_1), DIMENSION(mp)		:: fwsoil ! soil water modifier of stom. cond.
    REAL(r_1), DIMENSION(mp)		:: gaw ! aerodynamic conduct. for water
    REAL(r_1), DIMENSION(mp,mf)		:: gaw2	! 2D gaw
    REAL(r_1), DIMENSION(mp,mf)		:: gbhf ! freeConvectionBndLayerConductance mol/m2/s
    REAL(r_1), DIMENSION(mp,mf)		:: gbhu ! forcedConvectionBoundaryLayerConductance
    REAL(r_1), DIMENSION(mp)		:: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp,mf)		:: gras ! Grashof coeff
    REAL(r_1), DIMENSION(mp,mf)		:: gswmin ! min stomatal conductance
!    REAL(r_1), DIMENSION(mp,mf)	:: gswx ! stom cond for water
    REAL(r_1), DIMENSION(mp,mf)		:: gw  ! cond for water for a dry canopy
    REAL(r_1), DIMENSION(mp,mf)		:: gh  ! cond for heat for a dry canopy
    REAL(r_1), DIMENSION(mp,mf)		:: ghr ! dry canopy cond for heat & thermal radiat'n
    REAL(r_1), DIMENSION(mp)		:: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp)		:: ghwet  ! cond for heat for a wet canopy
    REAL(r_1), DIMENSION(mp)		:: ghrwet ! wet canopy cond: heat & thermal radiat'n
    REAL(r_1), DIMENSION(mp,mf)		:: hcx ! sens heat fl big leaf prev iteration
    REAL(r_1), DIMENSION(mp,mf)		:: hcy ! veg. sens heat
    INTEGER(i_d)			:: iter ! iteration #
    INTEGER(i_d)			:: iterplus !
    INTEGER(i_d)			:: k		! interation count
    INTEGER(i_d)			:: kk		! interation count
    REAL(r_1), DIMENSION(mp,mf)		:: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp,mf)		:: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_1), DIMENSION(mp,mf)		:: rdy ! daytime leaf resp rate
    REAL(r_1), DIMENSION(mp,mf)		:: rnx ! net rad prev timestep
    REAL(r_1), DIMENSION(mp,mf)		:: rny ! net rad
    REAL(r_1), DIMENSION(mp)		:: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp)		:: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp)		:: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL(r_1), DIMENSION(mp)		:: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp,mf)		:: tair2 ! 2D tair
    REAL(r_1), DIMENSION(mp,mf)		:: tvair2 ! 2D tair
    REAL(r_1), DIMENSION(mp,mf)		:: tdiff ! leaf air temp diff.
    REAL(r_1), DIMENSION(mp,mf)		:: tlfx ! leaf temp prev. iteration
    REAL(r_1), DIMENSION(mp,mf)		:: tlfxx ! leaf temperature of current iteration
    REAL(r_1), DIMENSION(mp,mf)		:: tlfy ! leaf temp
    REAL(r_1), DIMENSION(mp,mf)		:: vcmax2 ! vcmax big leaf
    REAL(r_1), DIMENSION(mp,mf)		:: vcmxt3 ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp,mf)		:: vcmxt4 ! vcmax big leaf C4
    REAL(r_1), DIMENSION(mp,mf,2)	:: vx3 ! carboxylation C3 plants
    REAL(r_1), DIMENSION(mp,mf,2)	:: vx4 ! carboxylation C4 plants
!    REAL(r_1), DIMENSION(mp)		:: wetfac ! degree of soil water limitation on stage 2 soil evaporation
    REAL(r_1), DIMENSION(mp,mf)		:: xdleaf2	! 2D dleaf
    REAL(r_1), DIMENSION(mp,mf)		:: xleuning ! leuning stomatal coeff
!    REAL(r_1), DIMENSION(mp,niter)	:: zetar ! stability correction
    REAL(r_1), PARAMETER		:: jtomol = 4.6e-6 ! Conversion from Joule to Mol for light
    REAL(r_1), PARAMETER		:: EHaVc  = 73637.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EHdVc  = 149252.0 !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EntropVc = 486.0  !J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER		:: xVccoef = 1.17461 !derived parameter
    !	xvccoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    REAL(r_1), PARAMETER		:: EHaJx  = 50300.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EHdJx  = 152044.0	!J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EntropJx = 495.0	!J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER		:: xjxcoef = 1.16715	!derived parameter
    REAL(r_1), PARAMETER		:: effc4 = 4000.0  !Vc=effc4*Ci*Vcmax (see
    ! Bonan,LSM version 1.0, p106)
!    REAL(r_1), DIMENSION(mp)		:: oldcansto ! prev t step canopy storage
    REAL(r_1), DIMENSION(mp)		:: cc ! limitation term for canopy interception per timestep		   
    REAL(r_1), DIMENSION(mp)		:: ccfevw ! limitation term for wet canopy evaporation rate  
    REAL(r_1), DIMENSION(mp)		:: denom ! denominator in calculating screen temperature, humidity etc
    REAL(r_1), DIMENSION(mp)		:: tstar ! 
    REAL(r_1), DIMENSION(mp)		:: zscrn !
    REAL(r_1), DIMENSION(mp)		:: qstar !
    REAL(r_1), DIMENSION(mp)		:: rsts  !
    REAL(r_1), DIMENSION(mp)		:: qsurf !
    REAL(r_1), DIMENSION(mp)		:: qtgnet !
    REAL(r_1), DIMENSION(mp)		:: evapfb !
    REAL(r_1), DIMENSION(mp,ms)		:: evapfbl !
    REAL(r_1), DIMENSION(mp)		:: phenps ! Leaf phenology influence on vcmax and jmax
    REAL(r_1), DIMENSION(mp)		:: poolcoef1 ! leaf carbon turnover rate * leaf pool size
    REAL(r_1), DIMENSION(mp)		:: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp)		:: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp)		:: rbw  ! leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp)		:: rrbw ! recipr. leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp)		:: rsw  ! stomatal resistance for water
    REAL(r_1), DIMENSION(mp)		:: rrsw ! recipr. stomatal resistance for water
    REAL(r_1), DIMENSION(mp)		:: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)		:: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)		:: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)		:: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)		:: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)		:: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp,mf)		:: tss2 ! 2D soil/snow temperature
    REAL(r_1), DIMENSION(mp)		:: tss4 ! soil/snow temperature**4
    REAL(r_1), DIMENSION(mp)		:: sss ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)		:: cc1 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)		:: cc2 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)		:: cc1T ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)		:: cc2T ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)		:: qstvair ! sat spec hunidity at leaf temperature
    REAL(r_1), DIMENSION(mp)		:: qstss ! sat spec hunidity at soil/snow temperature
    REAL(r_1), DIMENSION(mp)		:: xx ! delta-type function for sparse canopy limit, p20 SCAM manual
    REAL(r_2), DIMENSION(mp)		:: dxx 
    REAL(r_1), DIMENSION(mp,mf)		:: temp ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp,mf)         :: deltecy ! YP & Mao (jun08)

!%% changes by Ashok Luhar (low wind speed)
    REAL(r_1), PARAMETER		:: alpha1=4.0
    REAL(r_1), PARAMETER		:: beta1=0.5
    REAL(r_1), PARAMETER		:: gamma1=0.3
    REAL(r_1), DIMENSION(mp)	        :: zeta1
    REAL(r_1), DIMENSION(mp)	        :: zeta2
    REAL(r_1), DIMENSION(mp)	        :: usA
!%%

    !
    !	xjxcoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    ! 1-oct-2002 ypw: to keep the unit consistent for resistance or conductance
    ! s/m for r; mol/m2/s for g, and always use g where appropriate
    ! replace rh, rhr, rw  with ghdry/ghwet,ghrdry, ghrwet, gwdry, gwwet

    ! Set surface water vapour pressure deficit:
!    print *,'define_can',met%tc,met%pmb,met%qv
!    print *,'df-c ssoil%wb',ssoil%wb
!    print *,'df-c ssoil%tgg',ssoil%tgg
!    print *,'define_can 1',rmair,rmh2o
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    met%da = (qsatf(met%tc,met%pmb) - met%qv ) * rmair/rmh2o * met%pmb * 100.
    ! Soil water limitation on stomatal conductance:
    rwater = MAX(1.0e-4_r_2, &
         SUM(soil%froot * MIN(1.0_r_2,ssoil%wb - SPREAD(soil%swilt, 2, ms)),2) &
         /(soil%sfc-soil%swilt))
!    print *,'define_can 5',rwater,soil%froot(1,:),ktau,soil%froot(2,:)
    ! construct function to limit stom cond for soil water availability
    fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * veg%vlaiw
    ! Leaf phenology influence on vcmax and jmax
    !phenps = max (1.0e-4, MIN(1.,1. - ( (veg%tmaxvj - ssoil%tgg(:,4)+tfrz)/ &
    !     (veg%tmaxvj - veg%tminvj) )**2 ) )
    !WHERE ( ssoil%tgg(:,4) < (veg%tminvj + tfrz) ) phenps = 0.0
    !WHERE ( ssoil%tgg(:,4) > (veg%tmaxvj + tfrz) ) phenps = 1.0
    phenps = 1. ! MJT fix from Eva
    ! Set previous time step canopy water storage:
    canopy%oldcansto=canopy%cansto
    ! Rainfall variable is limited so canopy interception is limited,
    ! used to stabilise latent fluxes.
!    cc =min(met%precip, 4./(1440./(dels/60.)))! to avoid canopy temp. oscillations
!    cc =min(met%precip-met%precip_s, 4./(1440./(dels/60.)))! to avoid canopy temp. oscillations
    cc =min( met%precip-met%precip_s, 4./( 1440./( min(dels,1800.)/60.) ) )! to avoid canopy temp. oscillations

    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - canopy%cansto,0.0), cc), 0.0, &
         cc > 0.0  )
!         cc > 0.0  .AND. met%tk > tfrz)
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_s+MIN(met%precip-met%precip_s, &
                     MAX(0.0, met%precip-met%precip_s - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint

!    wetfac = MAX(0.0_r_2, MIN(1.0_r_2, &
    ssoil%wetfac = MAX(0.0_r_2, MIN(1.0_r_2, &
         (ssoil%wb(:,1) - soil%swilt) / (soil%sfc - soil%swilt)))
    ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)

    ! Temporay fixer for accounting of reduction of soil evaporation due to freezing
    where ( ssoil%wbice(:,1) > 0. )
       ! Prevents divide by zero at glaciated points where wb and wbice=0.
       ssoil%wetfac = ssoil%wetfac * ( 1.0 - ssoil%wbice(:,1)/ssoil%wb(:,1) )**2
!       wetfac = wetfac * ( 1.0 - ssoil%wbice(:,1)/ssoil%wb(:,1) )**2
    endwhere

!    print *,'define_can 8',canopy%wcint,cansat,canopy%cansto
!    print *,'define_can 9',canopy%through,canopy%cansto,ssoil%wetfac
!    print *,'zeta0',zeta0,zetpos

    canopy%zetar(:,1) = zeta0 ! stability correction terms
    canopy%zetar(:,2) = zetpos + 1 
    
!    print *,'zetar',canopy%zetar

    xdleaf2 = SPREAD(veg%dleaf, 2, mf) ! characteristic leaf length
    dsatdk2 = SPREAD(air%dsatdk, 2, mf)! deriv of sat vap pressure wrt temp
    ca2 = SPREAD(met%ca, 2, mf)        ! CO2 concentration
    csx = ca2                     ! initialise leaf surface CO2 concentration
    da2 = SPREAD(met%da, 2, mf)   ! water vapour pressure deficit
    dsx = da2                     ! init. leaf surface vpd
    tair2 = SPREAD(met%tc, 2, mf) ! air temp in C
    tss2 = SPREAD(ssoil%tss, 2, mf) ! 2D soil/snow 1st layer temperature
    ejmax2 = SPREAD(veg%ejmax*phenps, 2,mf) !max. pot. electr transp. rate top leaf(mol/m2s)
    vcmax2 = SPREAD(veg%vcmax*phenps, 2,mf) !max. RuBP carboxylsation rate top leaf(mol/m2s)
!    print *,'vcmax,ejmax',veg%vcmax,veg%ejmax
!    print *,'define_can 10',vcmax2,tss2,tair2,ca2
!    print *,'define_can 11',canopy%zetar(:,1),canopy%zetar(:,2)
!    print *,'define_can 12',xdleaf2,dsatdk2,ca2

    tlfx = tair2  ! initialise leaf temp iteration variable
    tlfy = tair2  ! initialise current leaf temp
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    ! weight min stomatal conductance by C3 an C4 plant fractions
    rdy = 0.0       ! init daytime leaf respiration rate
    rdx = 0.0       ! init daytime leaf respiration rate
    an_y = 0.0      ! init current estimate net photos.
!    gswx =1e-3     ! default stomatal conuctance
    canopy%gswx = 1e-3     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    anx = 0.0       ! init net photos. iteration memory variable
    ancj = 0.0    
    psycst = SPREAD(air%psyc, 2, mf) ! modified psyc. constant
    ! add on by ypw 1-oct-2002
    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    gwwet = 1.0e-3
    ghwet = 1.0e-3
    ghrwet= 1.0e-3
    ! Initialise in-canopy temperatures and humidity:
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    ortsoil = ssoil%rtsoil
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    tss4 = ssoil%tss**4
    DO iter = 1, niter
       CALL define_air (met, air)
       psycst = SPREAD(air%psyc, 2, mf)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)
       CALL radiation(ktau,bal, soil, ssoil, veg, air, met, rad)
       hcx = 0.0       ! init sens heat iteration memory variable
       ecx = rad%rniso ! init lat heat iteration memory variable
       rnx = rad%rniso ! init net rad iteration memory variable
       rny = rad%rniso ! init current estimate net rad
       hcy = 0.0       ! init current estimate lat heat
       ecy = rny - hcy ! init current estimate lat heat
       
       gswmin = rad%scalex * (gsw03 * (1. - frac42) + gsw04 * frac42)
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX(1.e-6, &
            vonk * MAX(met%ua,umin) / ( &
            LOG(rough%zref_uv / rough%z0m) - &
            psim(canopy%zetar(:,iter)) + &
            psim(canopy%zetar(:,iter) * rough%z0m / rough%zref_uv) ))
!%%change by Ashok Luhar - low wind formulation
            usA = 0.0
        where (canopy%zetar(:,iter) > 0.7)
            zeta1=canopy%zetar(:,iter) * rough%z0m / rough%zref_uv
!            usA = MAX(1.e-6, &
            canopy%us = MAX(1.e-6, &
            vonk * MAX(met%ua,umin) / ( &
            alpha1* ((canopy%zetar(:,iter)**beta1*  &
               (1.0+gamma1*canopy%zetar(:,iter)**(1.0-beta1)))  &       
             - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
        endwhere         
!%%
!        print 191,iter,ktau,met%ua(1),canopy%us,usA, &
!                  canopy%zetar(1,iter),canopy%zetar(2,iter), &
!                  rough%zref_uv(1),rough%zref_uv(2),met%tc,canopy%tv 
!191     format(1x,'windspeed',i2,i6,f6.1,1x,4f7.4,2x,2f8.3,2f6.0,9f7.1)
!        print *,'canopy%us',canopy%us,UMIN
       ! Turbulent aerodynamic resistance from roughness sublayer depth to reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
!              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * (LOG(rough%zref_tq/MAX(rough%zruffs-rough%disp, rough%z0soilsn)) &
        - psis( canopy%zetar(:,iter) ) &
        + psis( canopy%zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn))/rough%zref_tq ) &
          )/vonk

       ! rt0 = turbulent resistance from soil to canopy:
!!$       ! correction  by Ian Harman to rough%rt0us = f( canopy%zetar )
!!$       WHERE (veg%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn)
!!$       rough%rt0us  = 0.0
!!$       rt0old  = 0.0
!!$       ELSEWHERE
!!$!       rough%term6 =  exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
!!$       rt0old  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$            + (1-zdlin))*(EXP(2*csw*veg%vlaiw) - rough%term2)/rough%term3
!!$       rough%rt0us  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$!            - psis( canopy%zetar(:,iter) * rough%disp/rough%zref/rough%term6)  &
!!$!           + psis( canopy%zetar(:,iter) * rough%z0soilsn/rough%zref/rough%term6) &
!!$            + (1-zdlin))*(EXP(2*csw*veg%vlaiw) - rough%term2)/rough%term3 &
!!$              / rough%term6
!!$       ENDWHERE
!!$       rt0old = rt0old / canopy%us
!!$       rt0 = max(5.,rough%rt0us / canopy%us)
       rt0 = rough%rt0us / canopy%us
!       print *,'rt1usc',rt1usc,rt0,rough%rt0us
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)

!       print *,'rt1',rough%rt1
       WHERE (ssoil%snowd > 0.1)
          ssoil%wetfac = 1.
!          wetfac = 1.
       END WHERE
       ssoil%rtsoil = rt0 + rough%rt1*(0.5+sign(0.5,0.02-veg%vlaiw)) 
       ssoil%rtsoil = max(25.,ssoil%rtsoil)   
!       print *,'rtsoil',ssoil%rtsoil
       WHERE ( ssoil%rtsoil .GT. 2.* ortsoil .OR. ssoil%rtsoil .LT. 0.5*ortsoil )
          ssoil%rtsoil = MAX(25.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE
!       print *,'ssoil%rtsoil',ssoil%rtsoil
       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf *	&
            (canopy%us / MIN(rough%usuh, 0.2) * &
            veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
       ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
!                                gbhu corrected by F.Cruz & A.Pitman on 13/03/07
       gbhu(:,1) = gbvtop*(1.0-EXP(-veg%vlaiw*(0.5*rough%coexp+rad%extkb))) / &
                                               (rad%extkb+0.5*rough%coexp)
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
                            (1.0-EXP(-0.5*rough%coexp*veg%vlaiw))-gbhu(:,1)
!       print *,'gbvtop,gbhu',gbvtop,gbhu
!       print *,'rough%rt1',rough%rt1
       ! Aerodynamic conductance:
       gaw = air%cmolar / rough%rt1
       WHERE (veg%meth > 0 )
          gaw=100000.0
       END WHERE
       gaw2 = SPREAD(gaw, 2, mf)
       abs_deltlf = 999.0
       deltlfy = 999.0
       ! Initialise, over each gridpoint, sunlit and shaded leaves:
       DO k=1,mp
          DO kk=1,mf
             IF(rad%fvlai(k,kk) <=1.0e-2) THEN
                abs_deltlf(k,kk)=0.0
                hcx(k,kk) = 0.0 ! intialise
                ecx(k,kk) = 0.0 ! intialise
                anx(k,kk) = 0.0 ! intialise
                rnx(k,kk) = 0.0 ! intialise
                rny(k,kk) = rnx(k,kk) ! store initial values
                hcy(k,kk) = hcx(k,kk) ! store initial values
                ecy(k,kk) = ecx(k,kk) ! store initial values
                rdy(k,kk) = rdx(k,kk) ! store initial values
                an_y(k,kk) = anx(k,kk) ! store initial values
             END IF
          ENDDO
       ENDDO
       deltlfy = abs_deltlf
       k = 0
       DO WHILE (ANY(abs_deltlf > 0.1)  .AND.  k < maxiter)
          k = k + 1
!          print *,'DO WHILE',k,abs_deltlf
          ! Where vegetation and no convergence...
          WHERE (rad%fvlai > 1e-2 .and. abs_deltlf > 0.1)
             ! Grashof number (Leuning et al, 1995) eq E4:
             gras = max(1.0e-6,1.595E8*ABS(tlfx-tair2)*(xdleaf2**3.))
             ! See Appendix E in (Leuning et al, 1995):
             gbhf = rad%fvlai*SPREAD(air%cmolar, 2, mf)*0.5*dheat *(gras**0.25)/xdleaf2
             ! Conductance for heat:
             gh = 1.0/(MIN(1e3, SPREAD(1.0/gaw, 2, mf) + 0.5/(gbhu+gbhf)))
             ! Conductance for heat and longwave radiation:
             ghr = rad%gradis+gh
             temp =  xvcmxt3(tlfx+tfrz)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Vcmax for C3 plants:
             vcmxt3 = (1.0-frac42)*vcmax2*rad%scalex *temp
             temp=  xvcmxt4(tlfx)
             ! Temperature of Vcmax for C4 plants (Collatz et al 1989):
             vcmxt4 = frac42 * vcmax2 * rad%scalex * temp
             temp= xejmxt3(tlfx+tfrz)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Jmax for C3 plants:
             ejmxt3 = (1.0-frac42) * ejmax2 * rad%scalex * temp
             ! Difference between leaf temperature and reference temperature:
             tdiff  =  tlfx+tfrz-trefk
             ! Michaelis menten constant of Rubisco for CO2:
             conkct = conkc0*EXP((ekc/(rgas* trefk)) *(1.-trefk/(tlfx+tfrz)))
             ! Michaelis menten constant of Rubisco for oxygen:
             conkot = conko0*EXP((eko/(rgas* trefk)) *(1.-trefk/(tlfx+tfrz)))
             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx(:,:,1) = conkct*(1.0+0.21/conkot)
             cx(:,:,2) = 2.0* gam0*(1.+gam1*tdiff + gam2*tdiff*tdiff) !gamma*
             ! All equations below in appendix E in Wang and Leuning 1998 are
             ! for calculating anx, csx and gswx for Rubisco limited, RuBP limited,
             ! sink limited
             vx3(:,:,1) = vcmxt3
             vx4(:,:,1) = vcmxt4
             temp = rad%qcan(:,:,1)*jtomol*(1.0-frac42)
             vx3(:,:,2) = ej3x(temp,ejmxt3)
             temp = frac42*rad%qcan(:,:,1)*jtomol
             vx4(:,:,2) = ej4x(temp,vcmxt4)
             rdx = (cfrd3*vcmxt3+cfrd4*vcmxt4)*SPREAD(fwsoil, 2, mf)
             xleuning = (1.0-frac42)*a1c3/(1.0+dsx/d0c3) +frac42*a1c4/(1.0+dsx/d0c4)
             xleuning = xleuning * SPREAD(fwsoil, 2, mf) / (csx-co2cp3)
             ! Rubisco limited:
             coef2(:,:,1) = gswmin/rgswc+xleuning *(vx3(:,:,1)-(rdx-vx4(:,:,1)))
             coef1(:,:,1) = (1.0-csx*xleuning) *(vx3(:,:,1)+vx4(:,:,1)-rdx)	&
                  +(gswmin/rgswc)*(cx(:,:,1)-csx) -xleuning*(vx3(:,:,1)*cx(:,:,2)/2.0 &
                  +cx(:,:,1)*(rdx-vx4(:,:,1)))
             coef0(:,:,1) = -(1.0-csx*xleuning) *(vx3(:,:,1)*cx(:,:,2)/2.0  &
                  +cx(:,:,1)*(rdx-vx4(:,:,1))) -(gswmin/rgswc)*cx(:,:,1)*csx
             ! Discriminant in quadratic in eq. E7 Wang and Leuning, 1998
             delcx(:,:,1) = coef1(:,:,1)**2 -4.0*coef0(:,:,1)*coef2(:,:,1)
             ci(:,:,1) = (-coef1(:,:,1)+SQRT(MAX(0.0,delcx(:,:,1)))) /(2.0*coef2(:,:,1))
             ci(:,:,1) = MAX(0.0,ci(:,:,1))
             ancj(:,:,1) = vx3(:,:,1)*(ci(:,:,1)-cx(:,:,2)/2.0)	&
                  / (ci(:,:,1) + cx(:,:,1)) + vx4(:,:,1) - rdx
             ! RuBP limited:
             coef2(:,:,2) = gswmin/rgswc+xleuning *(vx3(:,:,2)-(rdx-vx4(:,:,2)))
             coef1(:,:,2) = (1.0-csx*xleuning) *(vx3(:,:,2)+vx4(:,:,2)-rdx)	&
                  +(gswmin/rgswc)*(cx(:,:,2)-csx) -xleuning*(vx3(:,:,2)*cx(:,:,2)/2.0 &
                  +cx(:,:,2)*(rdx-vx4(:,:,2)))
             coef0(:,:,2) = -(1.0-csx*xleuning) *(vx3(:,:,2)*cx(:,:,2)/2.0  &
                  +cx(:,:,2)*(rdx-vx4(:,:,2))) -(gswmin/rgswc)*cx(:,:,2)*csx
             delcx(:,:,2) = coef1(:,:,2)**2 -4.0*coef0(:,:,2)*coef2(:,:,2)
             ci(:,:,2) = (-coef1(:,:,2)+SQRT(MAX(0.0,delcx(:,:,2)))) /(2.0*coef2(:,:,2))
             ci(:,:,2) = MAX(0.0,ci(:,:,2))
             ancj(:,:,2) = vx3(:,:,2)*(ci(:,:,2)-cx(:,:,2)/2.0)	&
                  /(ci(:,:,2)+cx(:,:,2)) +vx4(:,:,2)-rdx
             ! Sink limited:
             coef2(:,:,3) = xleuning
             coef1(:,:,3) = gswmin/rgswc + xleuning * (rdx - 0.5*vcmxt3)  +  &
                  effc4 * vcmxt4 - xleuning * csx * effc4 *vcmxt4
             coef0(:,:,3) = -(gswmin/rgswc)*csx *effc4*vcmxt4 +	&
                  (rdx -0.5*vcmxt3)*gswmin/rgswc
             delcx(:,:,3) = coef1(:,:,3)**2 -4.0*coef0(:,:,3)*coef2(:,:,3)
             ancj(:,:,3)  = (-coef1(:,:,3)+SQRT(MAX(0.0,delcx(:,:,3)))) &
                  /(2.0*coef2(:,:,3))
             anx = MIN(ancj(:,:,1),ancj(:,:,2),ancj(:,:,3))
             csx = ca2 - anx * (1.0/gaw2+rgbwc/(gbhu + gbhf))
             canopy%gswx = gswmin+MAX(0.0,rgswc*xleuning *anx)
             ! Recalculate conductance for water:
             gw = 1.0/(1.0/canopy%gswx + 1.0/(1.075*(gbhu+gbhf)) + SPREAD(1.0/gaw, 2, mf))
             ! Modified psychrometric constant (Monteith and Unsworth, 1990)
             psycst = SPREAD(air%psyc, 2, mf) *ghr/gw
             ! Store leaf temperature:
             tlfxx = tlfx
             ! Update canopy latent heat flux:
             ecx = (dsatdk2*rad%rniso +capp*rmair*da2*ghr) /(dsatdk2+psycst)
             ! Update canopy sensible heat flux:
             hcx = (rad%rniso-ecx)*gh/ghr
             ! Update leaf temperature:
             tlfx=tair2+hcx/(capp*rmair*gh)
             !             tlfx=max(tlfx,min(tss2,tair2)-3.)
             !             tlfx=min(tlfx,max(tss2,tair2)+3.)
             ! Update net radiation for canopy:
             rnx = rad%rniso-capp*rmair*(tlfx -tair2)*rad%gradis
             ! Update leaf surface vapour pressure deficit:
            ! dsx = ecx*100.0* SPREAD(met%pmb, 2, mf) /(canopy%gswx*rmh2o*SPREAD(air%rlam, 2, mf))
             dsx = da2 + dsatdk2 * (tlfx-tair2)
             ! Store change in leaf temperature between successive iterations:
             deltlf = tlfxx-tlfx
             abs_deltlf = ABS(deltlf)
          END WHERE
!          print *,'def-can loop ecx',ecx,rad%rniso,dsatdk2,psycst
!          print *,'def-can loop tlx',tlfx
!          print *,'def-can loop rniso',rad%rniso
          ! Where leaf temp change b/w iterations is significant, and difference is 
          ! smaller than the previous iteration, store results:
          WHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
             deltlfy = deltlf
             tlfy = tlfx
             rny = rnx
             hcy = hcx
             ecy = ecx
             rdy = rdx
             an_y = anx
          END WHERE
          WHERE (abs_deltlf > 0.1)
             !        after four iteration, take the mean value of current and previous estimates
             !        as the next estimate of leaf temperature, to avoid oscillation
             tlfx = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx + &
                  (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx
          END WHERE
          IF(k==1) THEN
             !        taken the first iterated estimates as the defaults
             tlfy = tlfx
             rny = rnx
             hcy = hcx
             ecy = ecx
             rdy = rdx
             an_y = anx
          END IF
       END DO  ! DO WHILE (ANY(abs_deltlf > 0.1)	.AND.  k < maxiter)

       ! VEGETATION SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
       ! calculate total thermal resistance, rthv in s/m
       ghwet=gaw*SUM(gbhu,2)/(0.5*gaw+SUM(gbhu,2))
       gwwet=1.075*SUM(gbhu,2)*gaw/(1.075*SUM(gbhu,2)+gaw)
       WHERE (veg%meth > 0 )
          ghwet=2.*SUM(gbhu,2)
          gwwet=1.075*SUM(gbhu,2)
       END WHERE
       ghrwet=SUM(rad%gradis,2)+ghwet
       ! Calculate fraction of canopy which is wet:
       veg%fwet   = max(0.0,min(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))
       ! Calculate lat heat from wet canopy, may be neg. if dew onto wet canopy
       ! to avoid excessive evaporation:
       ccfevw = MIN(canopy%cansto*air%rlam/dels, 2./(1440./(dels/60.))*air%rlam)
!                             potential evaporation from vegetation
       canopy%fevw_pot = (air%dsatdk*SUM(rad%rniso,2)+ &
          capp*rmair*met%da*ghrwet) /(air%dsatdk+air%psyc*ghrwet/gwwet)
       
       canopy%fevw = MIN(veg%fwet*((air%dsatdk*SUM(rad%rniso,2)+ &
            capp*rmair*met%da*ghrwet) &
            /(air%dsatdk+air%psyc*ghrwet/gwwet)), ccfevw)
       ! Calculate sens heat from wet canopy:
       canopy%fhvw = (veg%fwet*SUM(rad%rniso,2)-canopy%fevw)*ghwet/ghrwet
       ! Calculate (dry) canopy transpiration, may be negative if dew
       canopy%fevc = (1.0 - veg%fwet) * sum(ecy,2)
!       print *,'define_can 13 0',canopy%fevc,veg%fwet,sum(ecy,2)
!       canopy%fevc = max(0.,canopy%fevc)
!sxy       evapfb = 0.
       evapfbl = 0.
       DO k = 1,ms
          WHERE (canopy%fevc > 0.)
             evapfb =canopy%fevc  * dels/air%rlam ! convert to mm/dt
	     dxx = evapfb*soil%froot(:,k)
!            evapfbl(:,k) =min(evapfb*soil%froot(:,k),max(0._r_2,min(ssoil%wb(:,k)-soil%swilt, & 
             evapfbl(:,k) =min(dxx,max(0._r_2,min(ssoil%wb(:,k)-soil%swilt, & 
                  ssoil%wb(:,k)-1.05*ssoil%wbice(:,k)))*soil%zse(k)*1000.)
          END WHERE
!         print *,'def_can evap',k,canopy%fevc,evapfb,air%rlam,dels,evapfbl(:,k)
       END DO
!          print *,'def_can evap l',canopy%fevc,evapfb,air%rlam,dels,evapfbl
!sxy       WHERE (evapfb > 0 )
!sxy        canopy%fevc=sum(evapfbl,2)*air%rlam/dels
!        canopy%fevc=(evapfbl(:,1)+evapfbl(:,2)+evapfbl(:,3)+ &
!           evapfbl(:,4)+evapfbl(:,5)+evapfbl(:,6))*air%rlam/dels
!sxy       END WHERE
        WHERE(SUM(ecy,2)>0.0)
          canopy%fevc=SUM(evapfbl,2)*air%rlam/dels
        END WHERE 

!          print *,'def_can evap l2',canopy%fevc,evapfb,air%rlam,dels,evapfbl
       ! Calculate latent heat from vegetation:
       canopy%fev = canopy%fevc + canopy%fevw
       ! recalculate for checking energy balance (YP & Mao, jun08)
       deltecy(:,1) = (canopy%fevc/(max((1.0-veg%fwet),1.0e-10)))*rad%fvlai(:,1) &
                      /(rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
       deltecy(:,2) = (canopy%fevc/(max((1.0-veg%fwet),1.0e-10)))*rad%fvlai(:,2) &
                      /(rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
       ecy(:,:)     = deltecy(:,:)
       hcy(:,:)     = (rad%rniso(:,:) - ecy(:,:))*gh(:,:)/ghr(:,:)
       tvair2 = SPREAD(met%tvair-tfrz, 2, mf) ! within-canopy air temp in C
       tlfy(:,:)    = tvair2(:,:)+REAL(hcy(:,:),r_1)/(capp*rmair*gh(:,:))
       rny(:,:)     = rad%rniso(:,:) - capp*rmair * (tlfy(:,:) &
                    - tvair2(:,:)) * rad%gradis(:,:)
       ! Calculate sensible heat from vegetation:
       canopy%fhv = (1.0 - veg%fwet) * sum(hcy,2)  + canopy%fhvw
       ! Calculate net rad absorbed by canopy:
       canopy%fnv = (1.0-veg%fwet)*SUM(rny,2)+canopy%fevw+canopy%fhvw
!    print *,'define_can 13',canopy%fevc,canopy%fevw,canopy%fev
       ! canopy radiative temperature is calculated based on long-wave radiation balance
       ! Q_lw=Q_lw_iso - (1.0-fwet)*SUM(capp*rmair*(tlfy-tair)*gri - canopy%fhvw*gr/ghw
       ! Q_lw=(1-transd)*(L_soil+L_sky-2.0*L_can)
       ! therefore
       ! Q_lw_iso-Q_lw=2(1-transd)*emleaf*(Tv^4-Tc^4)
       !	    rad%lwabv = (1.0-veg%fwet)*(capp*rmair*(tlfy(:,1) - met%tc)*rad%gradis(:,1) &
       !		 +capp*rmair*(tlfy(:,2) - met%tc)*rad%gradis(:,2)) &
       !		 + canopy%fhvw*SUM(rad%gradis,2)/ghwet
       !YP & Mao (jun08) replaced met%tk with tvair2
       rad%lwabv = (1.0-veg%fwet)*(capp*rmair*(tlfy(:,1) - &
            tvair2(:,1))*rad%gradis(:,1) &
            +capp*rmair*(tlfy(:,2) - tvair2(:,2))*rad%gradis(:,2)) &
            + canopy%fhvw*SUM(rad%gradis,2)/ghwet
       ! add if condition here to avoid dividing by zero ie when rad%transd=1.0 Ypw:24-02-2003
       WHERE (veg%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
       !YP & Mao (jun08) replaced met%tk with met%tvair 
       !   canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tk**4)**0.25
          canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tvair**4)**0.25
       ELSEWHERE ! sparse canopy
          !canopy%tv = met%tk
          canopy%tv = met%tvair
       END WHERE
       ! Calculate ground heat flux:
!       canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
       ! Saturation specific humidity at soil/snow surface temperature:
       qstss = qsatf((ssoil%tss - tfrz),met%pmb)
       ! Spec hum deficit at soil/snow surface:
       dq = qstss - met%qv
       !			       excessive dew over snow area
       WHERE (ssoil%snowd > 0.1)
          dq = max( -0.1e-3, dq)
       END WHERE
       ! Calculate net rad to soil:
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
            sboltz*canopy%tv**4 - emsoil*sboltz* tss4
!       ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
       ! Penman-Monteith formula
          sss=air%dsatdk
          cc1=sss/(sss+air%psyc )
          cc2=air%psyc /(sss+air%psyc )
!          ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) +  &
          ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
          cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil 
          cc1T = cc1 * (canopy%fns - canopy%ga)
          cc2T = cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
       ! Soil latent heat:
       canopy%fes= ssoil%wetfac * ssoil%potev
!       canopy%fes= wetfac * ssoil%potev
       WHERE (ssoil%snowd < 0.1 .and. canopy%fes .gt. 0. )
        canopy%fes= min(canopy%fes,max(0._r_2,(ssoil%wb(:,1)-soil%swilt))*soil%zse(1) &
               * 1000. * air%rlam / dels)
        canopy%fes = min(canopy%fes,(ssoil%wb(:,1)-ssoil%wbice(:,1))* soil%zse(1) &
               * 1000. * air%rlam / dels)
       END WHERE
       ssoil%cls=1.
       WHERE (ssoil%snowd >= 0.1)
          ssoil%cls = 1.1335
          canopy%fes= min(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
!          canopy%fes= min(wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
       END WHERE
       ! Calculate soil sensible heat:
       canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
       ! Calculate ground heat flux:
       canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
!       print *,'cdefine',canopy%fns,canopy%fhs,canopy%fes*ssoil%cls,ssoil%tss,met%tk,ssoil%rtsoil
       ! Calculate total latent heat:
       canopy%fe = canopy%fev + canopy%fes
       ! Calculate total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs
       ! Initialise in-canopy temperature and humidity:
       met%tvair = met%tk
       met%qvair = met%qv
!       print *,'canopy_module fe',canopy%fe,canopy%fev,canopy%fevc,canopy%fevw,canopy%fh
!       print 106,iter,ktau,sum(gbhu+gbhf,2),rbw,sum(canopy%gswx,2),rsw
!106    format(1x,'rbw,rsw',i3,i6,10f10.4)

      WHERE (veg%meth > 0 .and. veg%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn) 
          !      use the dispersion matrix (DM) to find the air temperature and specific humidity 
          !      (Raupach, Finkele and Zhang 1997, pp 17)
          ! leaf boundary layer resistance for water
          rbw = air%cmolar/sum(gbhu+gbhf,2)
          rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT 
          ! leaf stomatal resistance for water
          rsw = air%cmolar/sum(canopy%gswx,2)
          rrsw = sum(canopy%gswx,2)/air%cmolar ! MJT
          ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmah = (rt0+rough%rt1)*((1.+air%epsi)/rsw +1.0/rbw) &
!               + air%epsi * (rt0*rough%rt1)/(rbw*rsw)
          dmah = (rt0+rough%rt1)*((1.+air%epsi)*rrsw + rrbw) &
               + air%epsi * (rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmbh = (-air%rlam/capp)*(rt0*rough%rt1)/(rbw*rsw)
          dmbh = (-air%rlam/capp)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmch = ((1.+air%epsi)/rsw +1.0/rbw)*rt0*rough%rt1* &
          dmch = ((1.+air%epsi)*rrsw + rrbw)*rt0*rough%rt1* &
               (canopy%fhv + canopy%fhs)/(air%rho*capp)
          ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)/(rbw*rsw)
          dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmbe = (rt0+ssoil%wetfac*rough%rt1)*((1.+air%epsi)/rsw +1.0/rbw)+(rt0*rough%rt1)/(rbw*rsw)
          dmbe = (rt0+ssoil%wetfac*rough%rt1)* &
!          dmbe = (rt0+wetfac*rough%rt1)* &
                 ((1.+air%epsi)*rrsw + rrbw)+(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmce = ((1.+air%epsi)/rsw +1.0/rbw)*rt0*rough%rt1*
          dmce = ((1.+air%epsi)*rrsw + rrbw)*rt0*rough%rt1* &
                 (canopy%fev + canopy%fes)/(air%rho*air%rlam)
          ! Within canopy air temperature:
          met%tvair = met%tk + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          met%tvair = max(met%tvair , min( ssoil%tss, met%tk) - 5.0)
          met%tvair = min(met%tvair , max( ssoil%tss, met%tk) + 5.0)	  
          ! Within canopy specific humidity:
          met%qvair = met%qv + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          met%qvair = max(0.0,met%qvair)
       END WHERE
       where (met%tvair.lt.min(met%tk,ssoil%tss)-4.9) ! MJT
         met%qvair=met%qv
       end where
       where (met%tvair.gt.max(met%tk,ssoil%tss)+4.9) ! MJT
         met%qvair=met%qv
       end where
       ! Saturated specific humidity in canopy:
       qstvair = qsatf((met%tvair-tfrz),met%pmb)
       ! Saturated vapour pressure deficit in canopy:
       met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.
       ! 2 Dim saturated vapour pressure deficit in canopy:
       dva2 = SPREAD(met%dva, 2, mf)
       ! 2 dim Within canopy air temperature in degrees C:
       tvair2 = SPREAD(met%tvair-tfrz, 2, mf)
       ! Set radiative temperature as within canopy air temp:
       met%tvrad = met%tvair
       ! recalculate using canopy within temperature
       !     where (veg%meth .eq. 0 )
! EAK  recalculate air%dsatdk 
       CALL define_air (met, air)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

       WHERE (rad%fvlai > 1e-2) ! where LAI of sunlit or shaded leaf is significant:
          ! Recalculate fluxes and leaf temperature using within canopy air vpd:
          ecy = (dsatdk2*rad%rniso +capp*rmair*dva2*ghr) /(dsatdk2+psycst)
          hcy = (rad%rniso-ecy)*gh/ghr
          !             tlfx=tvair+hcx/(capp*rmair*gh)
          tlfy=tvair2+hcy/(capp*rmair*gh)
          ! YP & Mao (jun08) added rny calculation here
          rny = rad%rniso - capp*rmair * (tlfy - tvair2) * rad%gradis
       END WHERE
       WHERE (veg%meth > 0 )
          canopy%fevc = (1.0 - veg%fwet) * sum(ecy,2)
          !WHERE (canopy%fevc > 0.)
          WHERE (SUM(ecy,2) > 0.0) 
             evapfb =canopy%fevc  * dels/air%rlam             ! convert to mm/dt
             ! Calcualte contribution by different soil layers to canopy transpiration:
	     dxx=evapfb*soil%froot(:,1)
!            evapfbl(:,1) =min(evapfb*soil%froot(:,1),max(0._r_2,min(ssoil%wb(:,1)-soil%swilt, &
             evapfbl(:,1) =min(dxx,max(0._r_2,min(ssoil%wb(:,1)-soil%swilt, &
                  ssoil%wb(:,1)-1.05*ssoil%wbice(:,1)))*soil%zse(1)*1000.)
	     dxx=evapfb*soil%froot(:,2)
!            evapfbl(:,2) =min(evapfb*soil%froot(:,2),max(0._r_2,min(ssoil%wb(:,2)-soil%swilt, &
             evapfbl(:,2) =min(dxx,max(0._r_2,min(ssoil%wb(:,2)-soil%swilt, &
                  ssoil%wb(:,2)-1.05*ssoil%wbice(:,2)))*soil%zse(2)*1000.)
	     dxx=evapfb*soil%froot(:,3)
!            evapfbl(:,3) =min(evapfb*soil%froot(:,3),max(0._r_2,min(ssoil%wb(:,3)-soil%swilt, &
             evapfbl(:,3) =min(dxx,max(0._r_2,min(ssoil%wb(:,3)-soil%swilt, &
                  ssoil%wb(:,3)-1.05*ssoil%wbice(:,3)))*soil%zse(3)*1000.)
	     dxx=evapfb*soil%froot(:,4)
!            evapfbl(:,4) =min(evapfb*soil%froot(:,4),max(0._r_2,min(ssoil%wb(:,4)-soil%swilt, &
             evapfbl(:,4) =min(dxx,max(0._r_2,min(ssoil%wb(:,4)-soil%swilt, &
                  ssoil%wb(:,4)-1.05*ssoil%wbice(:,4)))*soil%zse(4)*1000.)
	     dxx=evapfb*soil%froot(:,5)
!            evapfbl(:,5) =min(evapfb*soil%froot(:,5),max(0._r_2,min(ssoil%wb(:,5)-soil%swilt, &
             evapfbl(:,5) =min(dxx,max(0._r_2,min(ssoil%wb(:,5)-soil%swilt, &
                  ssoil%wb(:,5)-1.05*ssoil%wbice(:,5)))*soil%zse(5)*1000.)
	     dxx=evapfb*soil%froot(:,6)
!            evapfbl(:,6) =min(evapfb*soil%froot(:,6),max(0._r_2,min(ssoil%wb(:,6)-soil%swilt, &
             evapfbl(:,6) =min(dxx,max(0._r_2,min(ssoil%wb(:,6)-soil%swilt, &
                  ssoil%wb(:,6)-1.05*ssoil%wbice(:,6)))*soil%zse(6)*1000.)
            ! fevc recalculated within WHERE construct (YP & Mao jun08)
            ! so that negative fevc values can be retained
              canopy%fevc=sum(evapfbl,2)*air%rlam/dels
          END WHERE
!sxy          WHERE (evapfb > 0 )
!             canopy%fevc=(evapfbl(:,1)+evapfbl(:,2)+evapfbl(:,3)+ &
!             evapfbl(:,4)+evapfbl(:,5)+evapfbl(:,6))*air%rlam/dels
!sxy             canopy%fevc=sum(evapfbl,2)*air%rlam/dels
!sxy          END WHERE

          ! Set total vegetation latent heat:
          canopy%fev  = canopy%fevc + canopy%fevw
          ! recalculate for checking energy balance (YP & Mao, jun08)
          deltecy(:,1) = (canopy%fevc/(max((1.0-veg%fwet),1.0e-10)))*rad%fvlai(:,1) &
                      /(rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
          deltecy(:,2) = (canopy%fevc/(max((1.0-veg%fwet),1.0e-10)))*rad%fvlai(:,2) &
                      /(rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
          ecy(:,1)     = deltecy(:,1)
          ecy(:,2)     = deltecy(:,2)
          hcy(:,1)     = (rad%rniso(:,1) - ecy(:,1))*gh(:,1)/ghr(:,1)
          hcy(:,2)     = (rad%rniso(:,2) - ecy(:,2))*gh(:,2)/ghr(:,2)
          tlfy(:,1)    = tvair2(:,1)+REAL(hcy(:,1),r_1)/(capp*rmair*gh(:,1))
          tlfy(:,2)    = tvair2(:,2)+REAL(hcy(:,2),r_1)/(capp*rmair*gh(:,2))
          rny(:,1)     = rad%rniso(:,1) - capp*rmair * (tlfy(:,1) &
                       - tvair2(:,1)) * rad%gradis(:,1)
          rny(:,2)     = rad%rniso(:,2) - capp*rmair * (tlfy(:,2) &
                       - tvair2(:,2)) * rad%gradis(:,2)
          ! Set total vegetation sensible heat:
          canopy%fhv  = (1.0 - veg%fwet) * sum(hcy,2)  + canopy%fhvw
          ! Longwave absorbed by vegetation:
          ! YP & Mao (jun08) replaced met%tvair with tvair2
          rad%lwabv = (1.0-veg%fwet)*(capp*rmair*(tlfy(:,1) - tvair2(:,1))*rad%gradis(:,1) &
               +capp*rmair*(tlfy(:,2) - tvair2(:,2))*rad%gradis(:,2)) &
               + canopy%fhvw*SUM(rad%gradis,2)/ghwet
          ! Set canopy temperature:
          WHERE (rad%transd <= 0.98)
             canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf) + met%tvair**4)**0.25
          ELSEWHERE
             ! sparse canopy 
             canopy%tv = met%tvair
          END WHERE
          ! Ground heat flux:
!          canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
          dq = qstss - met%qvair
          WHERE (ssoil%snowd > 0.1)
             dq = max( -0.1e-3, dq)
          END WHERE
          ! Net radiation absorbed by soil: 
          canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
                         sboltz*canopy%tv**4 - emsoil*sboltz* tss4
!          ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
          sss=air%dsatdk
          cc1=sss/(sss+air%psyc )
          cc2=air%psyc /(sss+air%psyc )
!          ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) +  &
          ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
               cc2 * air%rho * air%rlam*(qsatf((met%tvair-tfrz),met%pmb) - met%qvair)/ssoil%rtsoil
!          cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
!          cc1T = cc1 * (canopy%fns - canopy%ghflux)
!          cc2T = cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
       ! Soil latent heat:
          canopy%fes= ssoil%wetfac * ssoil%potev
!          canopy%fes= wetfac * ssoil%potev
          WHERE (ssoil%snowd < 0.1 .and. canopy%fes .gt. 0. )
             canopy%fes= min(canopy%fes,max(0._r_2,(ssoil%wb(:,1)-soil%swilt))* soil%zse(1) &
                  * 1000. * air%rlam / dels)
             canopy%fes = min(canopy%fes,(ssoil%wb(:,1)-ssoil%wbice(:,1)) * soil%zse(1) &
                  * 1000. * air%rlam / dels)
          END WHERE
          ssoil%cls=1.
          WHERE (ssoil%snowd >= 0.1)
             ssoil%cls = 1.1335
             canopy%fes= min(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
!             canopy%fes= min(wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
          END WHERE
          ! Soil sensible heat:
          canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
          canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
          ! Set total latent heat:
          canopy%fe = canopy%fev + canopy%fes
          WHERE (ssoil%potev .ge. 0.)
             ssoil%potev = max(0.001,ssoil%potev)
          ELSEWHERE
             ssoil%potev = min(-0.002,ssoil%potev)
          END WHERE
          WHERE (canopy%fevw_pot .ge. 0.)
             canopy%fevw_pot = max(0.001,canopy%fevw_pot)
          ELSEWHERE
             canopy%fevw_pot = min(-0.002,canopy%fevw_pot)
          END WHERE
!          canopy%rnet = (1.-rad%transd)*canopy%fnv + rad%transd*canopy%fns  
          canopy%rnet = canopy%fnv + canopy%fns  
          canopy%epot = ((1.-rad%transd)*canopy%fevw_pot + rad%transd*ssoil%potev) * dels/air%rlam  ! convert to mm/day
          canopy%wetfac_cs = min(1.0,canopy%fe / (canopy%fevw_pot + ssoil%potev))
          WHERE ( canopy%wetfac_cs .le. 0. )  &
           canopy%wetfac_cs = max(0.,min(1., &
                      max(canopy%fev/canopy%fevw_pot,canopy%fes/ssoil%potev)))
          ! Set total sensible heat:
          canopy%fh = canopy%fhv + canopy%fhs
       END WHERE ! veg%meth > 0

!       print *,'NETRADIATION',canopy%rnet, canopy%fnv,canopy%fns
!       print *,'canopy%epot',canopy%fevw_pot,ssoil%potev,dels/air%rlam, &
!        rad%transd,canopy%epot,canopy%epot*air%rlam/dels,canopy%fev,canopy%fes, &
!        met%qv,met%qvair,met%tvair,met%tk,canopy%tv,ssoil%tss

!      print *,'cdefine2',canopy%fns,canopy%fhs,canopy%fes*ssoil%cls,ssoil%tss

!      print *,'canopy_module fe 2',veg%meth,canopy%fevc,canopy%fevw,canopy%fev,canopy%fe,canopy%fhv,canopy%fhvw
       ! monin-obukhov stability parameter zetar=zref/l
       !	recompute zetar for the next iteration, except on last iteration
       IF (iter < niter) THEN ! dont compute zetar on the last iter
          iterplus = max(iter+1,2)
          canopy%zetar(:,iterplus) = -(vonk*grav*rough%zref_tq*(canopy%fh+0.07*canopy%fe))/ &
               (air%rho*capp*met%tk*canopy%us**3)
          ! case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
          IF (niter == 2) THEN
             canopy%zetar(:,2) = zetmul * canopy%zetar(:,2)
             WHERE (met%fsd(:,3) ==  0.0)
                canopy%zetar(:,2) = 0.5 * canopy%zetar(:,2)
             END WHERE
          END IF
          !	constrain zeta to zetpos and zetneg (set in param0)
          canopy%zetar(:,iterplus) = min(zetpos,canopy%zetar(:,iterplus))	 ! zetar too +
          canopy%zetar(:,iterplus) = max(zetneg,canopy%zetar(:,iterplus))	 ! zetar too -
       END IF
    END DO	     ! do iter = 1, niter
     !   check with Ying Ping ******
    canopy%gswx_T = rad%fvlai(:,1)/max(0.01,veg%vlaiw(:))*canopy%gswx(:,1) &
           + rad%fvlai(:,2)/max(0.01,veg%vlaiw(:))*canopy%gswx(:,2)
!    print 36,ktau,niter,canopy%zetar(1,:),canopy%zetar(2,:),met%tk(1),canopy%tv(1),met%ua(1),canopy%us
!36     format(1x,'ZETAR1',i5,i2,4f7.2,3x,4f7.2,1x,2f6.1,1x,f5.1,2f7.4)

    canopy%cduv = canopy%us * canopy%us / (max(met%ua,umin))**2

! screen temp., windspeed and relative humidity at 1.5m

    tstar = - canopy%fh / ( air%rho*capp*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
    zscrn = max(rough%z0m,1.5-rough%disp)
!    denom = ( log(rough%zref_tq/zscrn)- psim(canopy%zetar(:,iterplus)) + &
!         psim(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk
    denom = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) + &
         psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk

!%% change by Ashok Luhar
    where (canopy%zetar(:,iterplus) > 0.7)
!            zeta2(:)=zetar(:,iterplus) * zscrn / rough%zref
            zeta2=canopy%zetar(:,iterplus) * zscrn / rough%zref_tq
            denom =alpha1* ((canopy%zetar(:,iterplus)**beta1* &
               (1.0+gamma1*canopy%zetar(:,iterplus)**(1.0-beta1)))  &       
             - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk 
     endwhere
!%%

    ! Calculate screen temperature:
    canopy%tscrn = met%tc - tstar * denom
    ! Calculate radiative/skin temperature:
    rad%trad = ( (1.-rad%transd)*canopy%tv**4 + rad%transd * ssoil%tss**4 )**0.25
!   temporary fix for screen temperature to stay in the range
!    where ( met%tc .gt. rad%trad-tfrz )  &  ! night time
!           canopy%tscrn = min(met%tc,max( rad%trad-tfrz,canopy%tscrn))
!    where ( met%tc .le. rad%trad-tfrz )  &  ! day time
!           canopy%tscrn = min(rad%trad-tfrz,max( met%tc,canopy%tscrn))
    rsts = qsatf(canopy%tscrn, met%pmb)
    qtgnet = rsts * ssoil%wetfac - met%qv
!    qtgnet = rsts * wetfac - met%qv
    WHERE (qtgnet .gt. 0. )
       qsurf = rsts * ssoil%wetfac
!       qsurf = rsts * wetfac
    ELSEWHERE
       qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
!       qsurf = 0.1*rsts*wetfac + 0.9*met%qv
    END WHERE
!    canopy%qscrn = qsurf + qstar * denom
!    canopy%qscrn = (met%qv - qstar * denom)/qsatf(canopy%tscrn,met%pmb)
    canopy%qscrn = met%qv - qstar * denom
!   temporary fix for screen humidity to stay in the range
!    where ( met%qv .gt. qsurf )  &  
!           canopy%qscrn = min(met%qv,max( qsurf,canopy%qscrn))
!    where ( met%qv .le. qsurf )  &  
!           canopy%qscrn = min(qsurf,max( met%qv,canopy%qscrn))
!#if defined(SCMA)
!    IF(L_EXPLICIT)THEN
!      print *,'rough%disp',rough%disp,zscrn,canopy%zetar(:,iterplus),rough%z0m 
!      print 88,ktau,canopy%qscrn,qsurf,met%qv,qstar,denom,ssoil%wb(1,1), &
!             ssoil%wb(2,1),canopy%fe 
!      print 89,ktau,canopy%tscrn,rad%trad-273.16,met%tc,ssoil%tss-273.16,canopy%tv-273.16,tstar
!      print 90,ktau,canopy%fe,canopy%fes,canopy%fev,canopy%fevc,canopy%fevw,canopy%fh,canopy%fhs,&
!                    canopy%fhv,canopy%fhvw,canopy%us,met%ua/10.
!      
!    ENDIF
!#endif

!88  format(x,'screenvarqv',i5,2f8.5,2x,2f8.5,2x,2f8.5,2x,2f8.5,x,2f5.1,2x,2f5.2,2f5.0)
!89  format(1x,'screenvartsc',i5,1x,2f7.1,1x,2f7.1,1x,2f7.1,4f7.1,2x,2f6.2)
!90  format(1x,'canopy ls',i5,18f8.3,1x,4f7.3)

!    uscrn not required
!    canopy%uscrn = max(0., max(met%ua,umin) - canopy%us * denom ) !at present incorrect

    avgwrs = sum(soil%froot * ssoil%wb,2)
    avgtrs = max(0.0,sum(soil%froot * ssoil%tgg,2)-tfrz)
    poolcoef1=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) - &
         bgc%ratecp(1)*bgc%cplant(:,1))
    poolcoef1w=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))
    poolcoef1r=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))
    ! Carbon uptake from photosynthesis: 
    canopy%frp = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
         * poolcoef1 /(365.0*24.0*3600.0)		 ! 24/05
    canopy%frpw = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
         * poolcoef1w /(365.0*24.0*3600.0)		  ! 24/05
    canopy%frpr = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
         * poolcoef1r /(365.0*24.0*3600.0)		  ! 24/05

    ! This section to be updated as part of carbon module upgrade;
    ! frs is currently calculated in carbon module.
    canopy%frs  = rsoil(soil%rs20, avgwrs, avgtrs)
    canopy%frs  = canopy%frs &
         * sum(spread(bgc%ratecs,1, mp) * bgc%csoil,2)	&
         /(365.0*24.0*3600.0)		     !convert 1/year to 1/second
    WHERE (ssoil%snowd > 1.)
       canopy%frs	= canopy%frs / min(100.,ssoil%snowd)
    END WHERE
    canopy%frday = 12.0 * sum(rdy, 2)
    canopy%fpn = -12.0 * sum(an_y, 2)
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp

!    print 112,ktau,(canopy%frs(k),k=401,403),(avgwrs(k),k=401,403),(avgtrs(k),k=401,403) 
!112 format('soilfrs',i5,3f10.7,1x,3f6.2,3f6.1)
!    print 113,ktau,(canopy%fpn(k),k=401,403),(canopy%frp(k),k=401,403), &
!                   (canopy%frday(k),k=401,403),(canopy%frs(k),k=401,403), &
!                   (canopy%fnee(k),k=401,403) 
!113 format('Cfluxes',i5,15f8.5)
!    print 114,ktau,(canopy%fnpp(k),k=401,410)
!114  format('nppfluxes',i5,10f10.6)
    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
         dels * 1.0e3 / (rhow*air%rlam)
    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm
    ! Calculate canopy water storage excess:
    canopy%spill=max(0.,min(0.2*canopy%cansto,max(0.0, canopy%cansto-cansat)))
    ! Move excess canopy water to throughfall:
    canopy%through = canopy%through + canopy%spill
    ! Initialise 'throughfall to soil' as 'throughfall from canopy'; snow may absorb
    canopy%precis = max(0.,canopy%through)
    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill
    ! Modify canopy water storage for evaporation:
    canopy%cansto = max(canopy%cansto-max(0.0,canopy%fevw)*dels*1.0e3/ &
         (rhow*air%rlam), 0.0)
    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-canopy%oldcansto
    ! calculate dgdtg, derivative of ghflux
    ssoil%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssoil%tss  ! d(canopy%fns)/d(ssoil%tgg)
    ssoil%dfh_dtg = air%rho*capp/ssoil%rtsoil	   ! d(canopy%fhs)/d(ssoil%tgg)
    ssoil%dfe_ddq = ssoil%wetfac*air%rho*air%rlam/ssoil%rtsoil	! d(canopy%fes)/d(dq)
!    ssoil%dfe_ddq = wetfac*air%rho*air%rlam/ssoil%rtsoil	! d(canopy%fes)/d(dq)
    ssoil%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
         /((tetenc+ssoil%tss-tfrz)**2)*exp(tetenb*(ssoil%tss-tfrz)/(tetenc+ssoil%tss-tfrz))
    canopy%dgdtg = ssoil%dfn_dtg - ssoil%dfh_dtg - ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg
    !ypw: energy balance of the dry canopy
    !    bal%drybal=ecy(:,1)+ecy(:,2)+hcy(:,1)+hcy(:,2)-rad%rniso(:,1)-rad%rniso(:,2) &
!         +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*rad%gradis(:,1) &
!         +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*rad%gradis(:,2)
    bal%drybal=ecy(:,1)+hcy(:,1)-rad%rniso(:,1) &
         +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*rad%gradis(:,1)

    bal%drybal=ecy(:,2)+hcy(:,2)-rad%rniso(:,2) &
         +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*rad%gradis(:,2)
    !ypw: energy balance of the wet canopy
    bal%wetbal=canopy%fevw+canopy%fhvw-(rad%rniso(:,1)+rad%rniso(:,2))*veg%fwet &
         +canopy%fhvw*(rad%gradis(:,1)+rad%gradis(:,2))/ghwet

  CONTAINS
    !--------------------------------------------------------------------------
    ELEMENTAL FUNCTION qsatf(tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1)		  :: r    ! result; sat sp humidity
      r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END FUNCTION qsatf
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej3x(parx,x) result(z)
      REAL(r_1), INTENT(IN)	:: parx
      REAL(r_1), INTENT(IN)	:: x
      REAL(r_1)			:: z
      z = max(0.0, &
           0.25*((alpha3*parx+x-sqrt((alpha3*parx+x)**2 - &
           4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
    END FUNCTION ej3x
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej4x(parx,x) result(z)
      REAL(r_1), INTENT(IN)	:: parx
      REAL(r_1), INTENT(IN)	:: x
      REAL(r_1)			:: z
      z = max(0.0, &
           (alpha4*parx+x-sqrt((alpha4*parx+x)**2 - &
           4.0*convx4*alpha4*parx*x))/(2.0*convx4))
    END FUNCTION ej4x
    !---------------------------------------------------------
    ! Explicit array dimensions as temporary work around for NEC inlining problem
    FUNCTION xvcmxt4(x) result(z)
      REAL(r_1), PARAMETER	:: q10c4 = 2.0
      REAL(r_1), DIMENSION(mp,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp,mf)			:: z
      z = q10c4 ** (0.1 * x - 2.5) / &
           ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
    END FUNCTION xvcmxt4
    !---------------------------------------------------------
    FUNCTION xvcmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for vcmax for c3 plants
      REAL(r_1), DIMENSION(mp,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp,mf)		:: xvcnum
      REAL(r_1), DIMENSION(mp,mf)		:: xvcden
      REAL(r_1), DIMENSION(mp,mf)		:: z
      xvcnum=xvccoef*exp((ehavc/(rgas*trefk))*(1.-trefk/x))
      xvcden=1.0+exp((entropvc*x-ehdvc)/(rgas*x))
      z = max(0.0,xvcnum/xvcden)
    END FUNCTION xvcmxt3
    !---------------------------------------------------------
    FUNCTION xejmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants
      REAL(r_1), DIMENSION(mp,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp,mf)		:: xjxnum
      REAL(r_1), DIMENSION(mp,mf)		:: xjxden
      REAL(r_1), DIMENSION(mp,mf)		:: z
      xjxnum=xjxcoef*exp((ehajx/(rgas*trefk))*(1.-trefk/x))
      xjxden=1.0+exp((entropjx*x-ehdjx)/(rgas*x))
      z = max(0.0, xjxnum/xjxden)
    END FUNCTION xejmxt3
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psim(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psim(z/l) (z/l=zeta)
      ! for momentum, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      USE math_constants
      REAL(r_1), INTENT(IN)	:: zeta
      REAL(r_1)			:: r
      REAL(r_1)			:: x
      REAL(r_1), PARAMETER	:: gu = 16.0
      REAL(r_1), PARAMETER	:: gs = 5.0

      REAL(r_1)                 :: z 
      REAL(r_1)                 :: stable 
      REAL(r_1)                 :: unstable 
!      x = (1.0 + gu*abs(zeta))**0.25
!      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
!           + pi_c*0.5, -gs*zeta, zeta < 0.0)
             z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
      stable = -gs*zeta
      x      = (1.0 + gu*abs(zeta))**0.25
      unstable = alog((1.0+x*x)*(1.0+x)**2/8) - 2.0*atan(x) + pi_c*0.5
      r   = z*stable + (1.0-z)*unstable
       
    END FUNCTION psim
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psis(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      REAL(r_1), INTENT(IN)	:: zeta
      REAL(r_1)			:: r
      REAL(r_1), PARAMETER	:: gu = 16.0
      REAL(r_1), PARAMETER	:: gs = 5.0

      REAL(r_1)                 :: z
      REAL(r_1)                 :: y
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable

!      r = merge(2.0 * log((1.0 + sqrt(1.0 + gu * abs(zeta))) * 0.5), &
!           - gs * zeta, zeta < 0.0)

      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
      stable = -gs*zeta
      y      = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
      REAL(r_1), INTENT(IN)	:: rpconst
      REAL(r_1), INTENT(IN)	:: rpcoef
      REAL(r_1), INTENT(IN)	:: tair
      REAL(r_1)			:: z
      z = rpconst * exp(rpcoef * tair)
    END FUNCTION rplant
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rsoil(rsconst, avgwrs, avgtrs) result(z)
      REAL(r_1), INTENT(IN)	:: rsconst
      REAL(r_1), INTENT(IN)	:: avgwrs
      REAL(r_1), INTENT(IN)	:: avgtrs
      REAL(r_1)			:: z
      z = rsconst * min(1.0, max(0.0, min(&
           -0.0178+0.2883*avgwrs+5.0176*avgwrs*avgwrs-4.5128*avgwrs*avgwrs*avgwrs, &
           0.3320+22.6726*exp(-5.8184*avgwrs)))) &
           * min(1.0, max(0.0, min( 0.0104*(avgtrs**1.3053), 5.5956-0.1189*avgtrs)))
    END FUNCTION rsoil
    !---------------------------------------------------------
  END SUBROUTINE define_canopy
!  SUBROUTINE coeftest(rad,rough,air,met,dels,ssoil,soil,canopy, &
!             ftlt1,ftlt2,fqwt1,fqwt2)
!
!    TYPE (radiation_type), INTENT(INOUT):: rad
!    TYPE (roughness_type), INTENT(INOUT):: rough
!    TYPE (air_type), INTENT(INOUT)      :: air
!    TYPE (met_type), INTENT(INOUT)      :: met
!    REAL(r_1), INTENT(IN)               :: dels ! integration time setp (s)
!    TYPE (soil_snow_type), INTENT(INOUT):: ssoil
!    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
!    TYPE (canopy_type), INTENT(INOUT)   :: canopy
!!    REAL(r_1), DIMENSION(mp,niter)	:: zetar ! stability correction
!    INTEGER(i_d)                        :: iter ! iteration #
!    INTEGER(i_d)                        :: iterplus !
!    REAL(r_1), INTENT(IN)               :: ftlt1
!    REAL(r_1), INTENT(IN)               :: ftlt2
!    REAL(r_1), INTENT(IN)               :: fqwt1
!    REAL(r_1), INTENT(IN)               :: fqwt2
!    REAL(r_1), DIMENSION(mp)            :: ftl
!    REAL(r_1), DIMENSION(mp)            :: fqw
!    REAL(r_1), DIMENSION(mp)            :: usnew 
!!
!!    print *,'in testcoef',ftlt1,ftlt2,fqwt1,fqwt2
!    ftl(1) = ftlt1 * capp
!    ftl(2) = ftlt2 * capp
!    fqw(1) = fqwt1 *  air%rlam(1)/dels
!    fqw(2) = fqwt2 *  air%rlam(1)/dels
!
!!   to test calculation od u* with prescribed fluxes
!    canopy%zetar(:,1) = zeta0 ! stability correction terms
!    canopy%zetar(:,2) = zetpos + 1
!    DO iter = 1, niter
!!       canopy%us = MAX(1.e-6, &
!       usnew = MAX(1.e-6, &
!            vonk * MAX(met%ua,umin) / ( &
!            LOG(rough%zref_uv / rough%z0m) - &
!            psimn(canopy%zetar(:,iter)) + &
!            psimn(canopy%zetar(:,iter) * rough%z0m / rough%zref_uv) ))
!!
!       ! monin-obukhov stability parameter zetar=zref/l
!       !        recompute zetar for the next iteration, except on last iteration
!       IF (iter < niter) THEN ! dont compute zetar on the last iter
!          iterplus = max(iter+1,2)
!          canopy%zetar(:,iterplus) = -(vonk*grav*rough%zref_tq*(ftl+0.07*fqw))/ &
!               (air%rho*capp*met%tk*canopy%us**3)
!          ! case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
!          IF (niter == 2) THEN
!             canopy%zetar(:,2) = zetmul * canopy%zetar(:,2)
!             WHERE (met%fsd(:,3) ==  0.0)
!                canopy%zetar(:,2) = 0.5 * canopy%zetar(:,2)
!             END WHERE
!          END IF
!          !     constrain zeta to zetpos and zetneg (set in param0)
!          canopy%zetar(:,iterplus) = min(zetpos,canopy%zetar(:,iterplus))      ! zetar too +
!          canopy%zetar(:,iterplus) = max(zetneg,canopy%zetar(:,iterplus))      ! zetar too -
!       END IF
!    END DO           ! do iter = 1, niter
!! print *,'in testcoef1',iter,usnew,canopy%zetar(:,iterplus)
!! print *,'in testcoef1',iter,canopy%us,canopy%zetar(:,iterplus),usnew
!  CONTAINS
!    !--------------------------------------------------------------------------
!    !---------------------------------------------------------
!    ELEMENTAL FUNCTION psimn(zeta) result(r)
!      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
!      ! computes integrated stability function psim(z/l) (z/l=zeta)
!      ! for momentum, using the businger-dyer form for unstable cases
!      ! and the webb form for stable cases. see paulson (1970).
!      USE math_constants
!      REAL(r_1), INTENT(IN)     :: zeta
!      REAL(r_1)                 :: r
!      REAL(r_1)                 :: x
!      REAL(r_1), PARAMETER      :: gu = 16.0
!      REAL(r_1), PARAMETER      :: gs = 5.0
!      x = (1.0 + gu*abs(zeta))**0.25
!      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
!           + pi_c*0.5, -gs*zeta, zeta < 0.0)
!    END FUNCTION psimn
!    !---------------------------------------------------------
!  END SUBROUTINE  coeftest

END MODULE canopy_module
