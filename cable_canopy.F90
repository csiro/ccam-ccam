! cable_canopy.f90
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach,
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file has the canopy_module containing subroutines define_canopy, dryLeaf
! and wetLeaf.
! The functions included are:
!   qsatf,
!   ej3x,
!   ej4x,
!   xvcmxt4,
!   xvcmxt3,
!   xejmxt3,
!   psim,
!   psis, and
!   photosynthesis.
!
MODULE canopy_module
  USE photosynthetic_constants
  USE radiation_module
  USE roughness_module
  USE air_module
  USE define_types
  USE physical_constants
  IMPLICIT NONE
  REAL(r_1), DIMENSION(:), POINTER :: cansat ! max canopy intercept. (mm)
  REAL(r_1), DIMENSION(:), POINTER :: ghwet  ! cond for heat for a wet canopy
  REAL(r_1), DIMENSION(:), POINTER :: dsx ! leaf surface vpd
  REAL(r_1), DIMENSION(:), POINTER :: fwsoil ! soil water modifier of stom. cond
  REAL(r_1), DIMENSION(:), POINTER :: tlfx ! leaf temp prev. iter (K)
  REAL(r_1), DIMENSION(:), POINTER :: tlfy ! leaf temp (K)
  REAL(r_2), DIMENSION(:), POINTER :: ecy ! lat heat fl dry big leaf
  REAL(r_2), DIMENSION(:), POINTER :: hcy ! veg. sens heat
  REAL(r_2), DIMENSION(:), POINTER :: rny ! net rad
  REAL(r_1), DIMENSION(:,:), POINTER :: gbhu ! forcedConvectionBndryLayerCond
  REAL(r_1), DIMENSION(:,:), POINTER :: gbhf ! freeConvectionBndryLayerCond
  REAL(r_1), DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance
!  REAL(r_1), DIMENSION(:,:), POINTER :: gswx ! stom cond for water
  REAL(r_1), DIMENSION(:,:), POINTER :: csx ! leaf surface CO2 concentration
  PRIVATE
  PUBLIC define_canopy, sinbet
!  PUBLIC define_canopy, sinbet, coeftest
CONTAINS

  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy(ktau,dels,L_EXPLICIT)
    REAL(r_1), INTENT(IN)		:: dels ! integration time setp (s)
    LOGICAL, INTENT(IN)                 :: L_EXPLICIT
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp)		:: avgtrs !root weighted mean soil temperature
    REAL(r_1), DIMENSION(mp)		:: avgwrs !root weighted mean soil moisture
    REAL(r_1), DIMENSION(mp)		:: dq ! sat spec hum diff.
    REAL(r_1), DIMENSION(mp,mf)		:: dsatdk2	! 2D dsatdk
    REAL(r_2), DIMENSION(mp,mf)		:: ecx ! lat. hflux big leaf
    REAL(r_1), DIMENSION(mp,mf)		:: frac42	! 2D frac4
    REAL(r_1), DIMENSION(mp)		:: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp,mf)		:: gw  ! cond for water for a dry canopy
    REAL(r_1), DIMENSION(mp,mf)		:: gh  ! cond for heat for a dry canopy
    REAL(r_2), DIMENSION(mp,mf)		:: hcx ! sens heat fl big leaf prev iteration

    INTEGER(i_d)			:: iter ! iteration #
    INTEGER(i_d)			:: iterplus !
    INTEGER(i_d)			:: k		! interation count
    INTEGER(i_d)			:: kk		! interation count
    REAL(r_1), DIMENSION(mp,mf)		:: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp)		:: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp)		:: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp)		:: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL(r_1), DIMENSION(mp)		:: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp)		:: cc ! limitation term for canopy interception per timestep		   
    REAL(r_1), DIMENSION(mp)		:: ccfevw ! limitation term for wet canopy evaporation rate  
    REAL(r_1), DIMENSION(mp)		:: denom ! denominator in calculating screen temperature, humidity etc
!    REAL(r_1), DIMENSION(mp)            :: denom1 ! denominator in calculating screen temperature, humidity etc
!    REAL(r_1), DIMENSION(mp)            :: denom2 ! denominator in calculating screen temperature, humidity etc
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
    REAL(r_1), DIMENSION(mp,mf)		:: temp ! vcmax big leaf C3

!%% changes by Ashok Luhar (low wind speed)
    REAL(r_1), PARAMETER		:: alpha1=4.0
    REAL(r_1), PARAMETER		:: beta1=0.5
    REAL(r_1), PARAMETER		:: gamma1=0.3
    REAL(r_1), DIMENSION(mp)	        :: zeta1
    REAL(r_1), DIMENSION(mp)	        :: zeta2
    REAL(r_1), DIMENSION(mp)	        :: usA
!
!   temporary working variables
    REAL(r_1), DIMENSION(mp,3)	        :: xi 
    REAL(r_1), DIMENSION(mp,3)	        :: ti 
    REAL(r_1), DIMENSION(mp,3)	        :: si 
    REAL(r_1), DIMENSION(mp)	        :: r_sc 
    REAL(r_1), DIMENSION(mp)	        :: qscrn1,qscrn2,qscrn3 
    REAL(r_1), DIMENSION(mp)	        :: term1 
    REAL(r_1), DIMENSION(mp)	        :: term2 
    REAL(r_1), DIMENSION(mp)	        :: term3 
    REAL(r_1), DIMENSION(mp)	        :: term5 
    REAL(r_1), DIMENSION(mp)	        :: beta_sc 
    REAL(r_1), DIMENSION(mp)	        :: zscrn_10m
    REAL(r_1), DIMENSION(mp)	        :: ts1v 
    REAL(r_1), DIMENSION(mp)	        :: ts2v 
    REAL(r_1), DIMENSION(mp)	        :: ts3v 
    REAL(r_1), DIMENSION(mp)	        :: pnt1 
    REAL(r_1), DIMENSION(mp)	        :: pnt2 
    REAL(r_1), DIMENSION(mp)	        :: zscl 
    REAL(r_1), DIMENSION(mp)	        :: zscl_10m,zscl_scrn
    REAL(r_1), DIMENSION(mp)	        :: yy1,yy2,yy3,yy4 


    ALLOCATE(cansat(mp),ghwet(mp),gbhu(mp,mf))
    ALLOCATE(dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp))
    ALLOCATE(ecy(mp), hcy(mp), rny(mp))
    ALLOCATE(gbhf(mp,mf), gswmin(mp,mf), csx(mp,mf))
!    ALLOCATE(gbhf(mp,mf), gswmin(mp,mf), gswx(mp,mf), csx(mp,mf))

    ! Set surface water vapour pressure deficit:

    met%da = (qsatf(met%tc,met%pmb) - met%qv ) * rmair/rmh2o * met%pmb * 100.
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    ! Soil water limitation on stomatal conductance:
!   print*, "MRD RWATER", minval(soil%sfc), minval(soil%swilt), &
!             minloc(soil%sfc), minloc(soil%swilt)
    rwater = MAX(1.0e-4_r_2, &
         SUM(veg%froot * MIN(1.0_r_2,ssoil%wb - SPREAD(soil%swilt, 2, ms)),2) &
         /(soil%sfc-soil%swilt))
    ! construct function to limit stom cond for soil water availability
    fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw
    ! Leaf phenology influence on vcmax and jmax
    phenps = 1.0
!   WHERE(veg%deciduous)
!    phenps = max (1.0e-4, MIN(1.,1. - ( (veg%tmaxvj - ssoil%tgg(:,4)+tfrz)/ &
!         (veg%tmaxvj - veg%tminvj) )**2 ) )
!    WHERE ( ssoil%tgg(:,4) < (veg%tminvj + tfrz) ) phenps = 0.0
!    WHERE ( ssoil%tgg(:,4) > (veg%tmaxvj + tfrz) ) phenps = 1.0
!   ELSEWHERE
!    phenps = 1.0
!   END WHERE 

    ! Set previous time step canopy water storage:
    canopy%oldcansto=canopy%cansto
    ! canopy intercepted rainfall rate is limited to avoid excessive direct 
    ! canopy evaporation, modified further by timestep requirement (EAK aug08)
    cc =min( met%precip-met%precip_s, 4./( 1440./( min(dels,1800.)/60.) ) )

    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - canopy%cansto,0.0), cc), 0.0, &
         cc > 0.0  )
!         cc > 0.0  .AND. met%tk > tfrz)
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_s+MIN(met%precip-met%precip_s, &
                     MAX(0.0, met%precip-met%precip_s - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint
    ! Calculate fraction of canopy which is wet:
    canopy%fwet   = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))

    ssoil%wetfac = MAX(0.0_r_2, MIN(1.0_r_2, &
         (ssoil%wb(:,1) - soil%swilt) / (soil%sfc - soil%swilt)))

    ! owetfac introduced to reduce sharp changes in soil moisture in dry regions,
    ! especially in offline runs where there may be discrepancies between
    ! timing of precip and temperature change (EAK apr2009)
    ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)

    ! Temporay fixer for accounting of reduction of soil evaporation due to freezing
    ! ssoil%wetfac = wetfac * ( 1.0 - ssoil%wbice(:,1)/ssoil%wb(:,1) )**2
    where ( ssoil%wbice(:,1) > 0. )
       ! Prevents divide by zero at glaciated points where wb and wbice=0.
       ssoil%wetfac = ssoil%wetfac * ( 1.0 - ssoil%wbice(:,1)/ssoil%wb(:,1) )**2
    endwhere

    canopy%zetar(:,1) = zeta0         ! stability correction terms
    canopy%zetar(:,2) = zetpos + 1 

    ! weight min stomatal conductance by C3 an C4 plant fractions
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
!    gswx =1e-3      ! default stomatal conuctance
    canopy%gswx = 1e-3     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    gw = 1.0e-3     ! default values of conductance
    gh = 1.0e-3
    ghwet = 1.0e-3
    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    ortsoil = ssoil%rtsoil
    tss4 = ssoil%tss**4
!    denom1 = 0.0
!    denom2 = 0.0
    qstvair = qsatf((met%tvair-tfrz),met%pmb)
    met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.0
    dsx = met%dva     ! init. leaf surface vpd
    tlfx = met%tvair  ! initialise leaf temp iteration memory variable (K)
    tlfy = met%tvair  ! initialise current leaf temp (K)

    DO iter = 1, niter

       CALL define_air

       psycst = SPREAD(air%psyc, 2, mf)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

       CALL radiation(ktau)

       gswmin = max(1.e-6,rad%scalex * (gsw03 * (1. - frac42) + gsw04 * frac42))

       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX(1.e-6, &
            vonk * MAX(met%ua,umin) / ( &
            LOG(rough%zref_uv / rough%z0m) - &
            psim(canopy%zetar(:,iter)) + &
            psim(canopy%zetar(:,iter) * rough%z0m / rough%zref_uv) ))

!        change by Ashok Luhar - low wind formulation
!            usA = 0.0
!        where (canopy%zetar(:,iter) > 0.7)
!            zeta1=canopy%zetar(:,iter) * rough%z0m / rough%zref_uv
!!            usA = MAX(1.e-6, &
!            canopy%us = MAX(1.e-6, &
!            vonk * MAX(met%ua,umin) / ( &
!            alpha1* ((canopy%zetar(:,iter)**beta1*  &
!               (1.0+gamma1*canopy%zetar(:,iter)**(1.0-beta1)))  &       
!             - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
!        endwhere         
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
        + psis( canopy%zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn)) /  &
          rough%zref_tq ) )/vonk

       ! rt0 = turbulent resistance from soil to canopy:

!!$       ! correction  by Ian Harman to rough%rt0us = f( canopy%zetar )
!!$       WHERE (canopy%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn)
!!$       rough%rt0us  = 0.0
!!$       rt0old  = 0.0
!!$       ELSEWHERE
!!$!       rough%term6 =  exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
!!$       rt0old  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$            + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3
!!$       rough%rt0us  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$!            - psis( canopy%zetar(:,iter) * rough%disp/rough%zref/rough%term6)  &
!!$!           + psis( canopy%zetar(:,iter) * rough%z0soilsn/rough%zref/rough%term6) &
!!$            + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3 &
!!$              / rough%term6
!!$       ENDWHERE
!!$       rt0old = rt0old / canopy%us
!!$       rt0 = max(5.,rough%rt0us / canopy%us)

       rt0 = rough%rt0us / canopy%us
       !  print *,'rt1usc',rt1usc,rt0,rough%rt0us
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)

       WHERE (ssoil%snowd > 0.1)
          ssoil%wetfac = 1.
       END WHERE
       ssoil%rtsoil = rt0 + rough%rt1*(0.5+sign(0.5,0.01-canopy%vlaiw)) 
       ssoil%rtsoil = max(25.,ssoil%rtsoil)   
!       print *,'resistances',rough%rt1usa/canopy%us,rough%rt1usb/canopy%us,rt1usc/canopy%us,rough%rt1, &
!                             ssoil%rtsoil,rt0
       WHERE ( ssoil%rtsoil .GT. 2.* ortsoil .OR. ssoil%rtsoil .LT. 0.5*ortsoil )
          ssoil%rtsoil = MAX(25.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE
       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf *	&
            (canopy%us / MIN(rough%usuh, 0.2) * &
            veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
       ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
       gbhu(:,1) = gbvtop*(1.0-EXP(-canopy%vlaiw*(0.5*rough%coexp+rad%extkb))) / &
                                               (rad%extkb+0.5*rough%coexp)
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
                            (1.0-EXP(-0.5*rough%coexp*canopy%vlaiw))-gbhu(:,1)

       rny = SUM(rad%rniso,2) ! init current estimate net rad
       hcy = 0.0              ! init current estimate lat heat
       ecy = rny - hcy        ! init current estimate lat heat

       CALL dryLeaf(dels,phenps)

       CALL wetLeaf(dels)

       ! Calculate latent heat from vegetation:
       canopy%fev = REAL(canopy%fevc,r_1) + canopy%fevw
       ! Calculate sensible heat from vegetation:
       canopy%fhv = (1.0 - canopy%fwet) * REAL(hcy,r_1)  + canopy%fhvw
       ! Calculate net rad absorbed by canopy:
       canopy%fnv = (1.0-canopy%fwet)*REAL(rny,r_1)+canopy%fevw+canopy%fhvw

       ! canopy radiative temperature is calculated based on long-wave radiation balance
       ! Q_lw=Q_lw_iso - (1.0-canopy%fwet)*SUM(capp*rmair*(tlfy-tair)*gri - canopy%fhvw*gr/ghw
       ! Q_lw=(1-transd)*(L_soil+L_sky-2.0*L_can)
       ! therefore
       ! Q_lw_iso-Q_lw=2(1-transd)*emleaf*(Tv^4-Tc^4)
       !	    rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) - met%tc)*rad%gradis(:,1) &
       !		 +capp*rmair*(tlfy(:,2) - met%tc)*rad%gradis(:,2)) &
       !		 + canopy%fhvw*SUM(rad%gradis,2)/ghwet
!       rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) - &
!            tvair2(:,1))*rad%gradis(:,1) &
!            +capp*rmair*(tlfy(:,2) - tvair2(:,2))*rad%gradis(:,2)) &
!            + canopy%fhvw*SUM(rad%gradis,2)/ghwet
!       rad%lwabv = (1.0-canopy%fwet) * capp * rmair * (tlfy-met%tvair) &
!                 * SUM(rad%gradis,2) &
!                 + canopy%fhvw * SUM(rad%gradis,2) / ghwet
!       rad%lwabv = capp * rmair * (tlfy-met%tk)* SUM(rad%gradis,2)
       ! Calculate canopy temperature
       ! add if condition here to avoid dividing by zero ie when rad%transd=1.0 Ypw:24-02-2003
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          rad%lwabv = capp * rmair * (tlfy-met%tk)* SUM(rad%gradis,2) + &
                 canopy%fhvw*SUM(rad%gradis,2)/max(0.001,ghwet)
!          canopy%tv = ( rad%lwabv/(2.0*(1.0-rad%transd)*sboltz*emleaf) + met%tvair**4)**0.25
          canopy%tv = max( rad%lwabv/(2.*(1.-rad%transd)*sboltz*emleaf)  + met%tk**4,0. )**0.25 ! MJT
       ELSEWHERE ! sparse canopy
!          canopy%tv = met%tvair
          canopy%tv = met%tk
       END WHERE
       where (canopy%tv.lt.met%tk-50.) ! MJT
         canopy%tv=met%tk              ! MJT
       end where                       ! MJT

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
       ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
       ! Penman-Monteith formula
!          sss=air%dsatdk
!          cc1=sss/(sss+air%psyc )
!          cc2=air%psyc /(sss+air%psyc )
!!          ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) +  &
!          ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
!          cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil 
!          cc1T = cc1 * (canopy%fns - canopy%ga)
!          cc2T = cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
!       ! Soil latent heat:
       canopy%fes= ssoil%wetfac * ssoil%potev
       WHERE (ssoil%snowd < 0.1 .and. canopy%fes > 0. )
        ! Reduce for wilting point limitation:
        canopy%fes= min(canopy%fes,MAX(0._r_2,(ssoil%wb(:,1)-soil%swilt))*soil%zse(1) &
               * 1000. * air%rlam / dels)
        ! Reduce for soil ice limitation:
        canopy%fes = MIN(canopy%fes,(ssoil%wb(:,1)-ssoil%wbice(:,1))* soil%zse(1) &
               * 1000. * air%rlam / dels)
       END WHERE
       ssoil%cls=1.
       WHERE (ssoil%snowd >= 0.1)
          ssoil%cls = 1.1335
          canopy%fes= MIN(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
       END WHERE
       ! Calculate soil sensible heat:
       canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
       ! Calculate ground heat flux:
       canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
       ! Calculate total latent heat:
       canopy%fe = canopy%fev + canopy%fes
       ! Calculate total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs
       ! Initialise in-canopy temperature and humidity:
!       met%tvair = met%tk
!       met%qvair = met%qv

      WHERE (veg%meth > 0 .and. canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn) 
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
!          dmbe = (rt0+wetfac*rough%rt1)*((1.+air%epsi)/rsw +1.0/rbw)+(rt0*rough%rt1)/(rbw*rsw)
          dmbe = (rt0+ssoil%wetfac*rough%rt1)* &
                 ((1.+air%epsi)*rrsw + rrbw)+(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
!          dmce = ((1.+air%epsi)/rsw +1.0/rbw)*rt0*rough%rt1*
          dmce = ((1.+air%epsi)*rrsw + rrbw)*rt0*rough%rt1* &
                 (canopy%fev + canopy%fes)/(air%rho*air%rlam)
          ! Within canopy air temperature:
          met%tvair = met%tk + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          ! tvair clobbered at present
          met%tvair = max(met%tvair , min( ssoil%tss, met%tk) - 5.0) ! MJT fix from Eva
          met%tvair = min(met%tvair , max( ssoil%tss, met%tk) + 5.0) ! MJT fix from Eva 
          ! Within canopy specific humidity:
          met%qvair = met%qv + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          met%qvair = min(max(0.0,met%qvair),1.) ! MJT

!          denom1 = dmbe*dmch-dmbh*dmce
!          denom2 = dmah*dmbe-dmae*dmbh+1.0e-12

       END WHERE
!     IF( .not. L_EXPLICIT) THEN 
!    print 93, iter, ktau, denom1, denom2, dmbe, dmch, dmbh, dmce, dmbe*dmch, dmbh*dmce, dmah, dmae, dmah*dmbe, dmae*dmbh, rt0, rough%rt1, air%epsi, rrsw, rrbw, air%rlam, capp, canopy%fhv, canopy%fhs, air%rho, ssoil%wetfac, canopy%fev, canopy%fes, met%tvair, met%qvair, met%tk, met%qv, gbhu, gbhf, air%cmolar, canopy%fwet, hcy, canopy%fhvw, rough%rt1usa, rough%rt1usb, rt1usc, canopy%us, rough%rt0us, ssoil%tss, ssoil%rtsoil  
!,canopy%gswx 
!93 format('denom#2',i2,i6,122(2x,f15.4))
!     END IF

       ! Saturated specific humidity in canopy:
       qstvair = qsatf((met%tvair-tfrz),met%pmb)
       ! Saturated vapour pressure deficit in canopy:
       met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.

       ! Set radiative temperature as within canopy air temp:
       met%tvrad = met%tvair

       ! recalculate using canopy within temperature
       CALL define_air

       ! Ground heat flux:
!        canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
       dq = qstss - met%qvair
       WHERE (ssoil%snowd > 0.1)
          dq = max( -0.1e-3, dq)
       END WHERE
       ! Net radiation absorbed by soil: 
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
                         sboltz*canopy%tv**4 - emsoil*sboltz* tss4
       ! method alternative to P-M formula
       ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
!          sss=air%dsatdk
!          cc1=sss/(sss+air%psyc )
!          cc2=air%psyc /(sss+air%psyc )
!          ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
!               cc2 * air%rho * air%rlam*(qsatf((met%tvair-tfrz),met%pmb) - met%qvair)/ssoil%rtsoil
!!          cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
!!          cc1T = cc1 * (canopy%fns - canopy%ga)
!!          cc2T = cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
       ! Soil latent heat:
       canopy%fes= ssoil%wetfac * ssoil%potev
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
       canopy%rnet = canopy%fnv + canopy%fns  
       canopy%epot = ((1.-rad%transd)*canopy%fevw_pot + rad%transd*ssoil%potev) * dels/air%rlam  ! convert to mm/day
       canopy%wetfac_cs = min(1.0,canopy%fe / (canopy%fevw_pot + ssoil%potev))
       WHERE ( canopy%wetfac_cs .le. 0. )  &
         canopy%wetfac_cs = max(0.,min(1., &
                 max(canopy%fev/canopy%fevw_pot,canopy%fes/ssoil%potev)))
       ! Set total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs

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
          
          !canopy%zetar(:,iterplus) = 0.7*canopy%zetar(:,iterplus)+0.3*canopy%zetar(:,iterplus-1) ! MJT suggestion
       END IF
       !print *,"zetar ",iter,canopy%zetar(1301,iter),maxval(canopy%fpn),maxloc(canopy%fpn)

    END DO	     ! do iter = 1, niter

    canopy%gswx_T = rad%fvlai(:,1)/max(0.01,canopy%vlaiw(:))*canopy%gswx(:,1) &
           + rad%fvlai(:,2)/max(0.01,canopy%vlaiw(:))*canopy%gswx(:,2)
!    print 36,ktau,niter,canopy%zetar(1,:),canopy%zetar(2,:),met%tk(1),canopy%tv(1), &
!           met%ua(1),canopy%us
36     format(1x,'ZETAR1',i5,i2,4f7.2,3x,4f7.2,1x,2f6.1,1x,f5.1,2f7.4)

    canopy%cduv = canopy%us * canopy%us / (max(met%ua,umin))**2
    canopy%cdtq = canopy%cduv *(LOG(rough%zref_uv / rough%z0m) -          &
      psim( canopy%zetar(:,niter) * rough%zref_uv/rough%zref_tq )) /      &
      (LOG( rough%zref_uv /(.1*rough%z0m) ) - psis(canopy%zetar(:,niter)) )


    ! Calculate screen temperature:
    ! 1) original method from SCAM

    ! screen temp., windspeed and relative humidity at 1.8m
    tstar = -canopy%fh / ( air%rho*capp*canopy%us)
    qstar = -canopy%fe / ( air%rho*air%rlam *canopy%us)
!   zscrn = max(rough%z0m,1.5-rough%disp)
    zscrn = max(rough%z0m,1.8-rough%disp)
    denom = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) + &
         psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk

!%% change by Ashok Luhar
!    where (canopy%zetar(:,iterplus) > 0.7)
!            zeta2=canopy%zetar(:,iterplus) * zscrn / rough%zref_tq
!            denom =alpha1* ((canopy%zetar(:,iterplus)**beta1* &
!               (1.0+gamma1*canopy%zetar(:,iterplus)**(1.0-beta1)))  &       
!             - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk 
!     endwhere

    canopy%tscrn1 = met%tc - tstar * denom
    canopy%tscrn2 = canopy%tscrn1

    ! 2) fitting polynomial 
    where( rough%hruff/2.0 .gt. 1.8 )
      xi(:,1) = 0.
      xi(:,2) = rough%hruff/2.0
      xi(:,3) = rough%zref_tq
      ti(:,1) = ssoil%tss
      ti(:,2) = met%tvair
      ti(:,3) = met%tk
      si(:,1) = (1.8 - xi(:,2)) / ( xi(:,1) - xi(:,2)) * ( 1.8  - xi(:,3)) / (xi(:,1) - xi(:,3))
      si(:,2) = (1.8 - xi(:,1)) / ( xi(:,2) - xi(:,1)) * ( 1.8  - xi(:,3)) / (xi(:,2) - xi(:,3))
      si(:,3) = (1.8 - xi(:,1)) / ( xi(:,3) - xi(:,1)) * ( 1.8  - xi(:,2)) / (xi(:,3) - xi(:,2))
      canopy%tscrn2 = ti(:,1) * si(:,1) + ti(:,2) * si(:,2) + ti(:,3) * si(:,3) - 273.16
    endwhere 
!    print *,'canopy%tscrn2',canopy%tscrn2
    ! 3) resistance method by Ian Harman 
    term1=0.
    term2=0.
    term5=0.
    zscl = max(rough%z0soilsn,1.8)
    where ( rough%hruff  > 0.0 .and. rough%disp  > 0.0 )
       term1 = EXP(2*csw*canopy%vlaiw*(1-zscl/rough%hruff))
       term2 = EXP(2*csw*canopy%vlaiw*(1-rough%disp/rough%hruff))
       term5 = MAX(2./3.*rough%hruff/rough%disp, 1.)
    endwhere
    term3 = a33**2*ctl*2*csw*canopy%vlaiw
    pnt1 = 0
    where( zscl < rough%disp ) 
        r_sc = term5 * LOG(zscl/rough%z0soilsn) * ( exp(2*csw*canopy%vlaiw) - term1 ) / term3
        pnt1=1
    elsewhere ( rough%disp <= zscl .and. zscl < rough%hruff )
        r_sc = rough%rt0us + term5 * ( term2 - term1 ) / term3
        pnt1=2
    elsewhere ( rough%hruff <= zscl .and. zscl <  rough%zruffs )
        r_sc = rough%rt0us + rough%rt1usa + term5 * ( zscl - rough%hruff ) /  &
                                                          (a33**2*ctl*rough%hruff)
        pnt1=3
    elsewhere (zscl >= rough%zruffs ) 
        r_sc = rough%rt0us + rough%rt1usa + rough%rt1usb +  &
              ( log( (zscl - rough%disp)/MAX(rough%zruffs-rough%disp, rough%z0soilsn) ) &
                - psis( (zscl - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_tq )  &
                + psis( (rough%zruffs - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_tq )  &
              ) / vonk  
        pnt1=4
    endwhere 
    canopy%tscrn3 = ssoil%tss + (met%tk - ssoil%tss) * r_sc /                            &
        (rough%rt0us + rough%rt1usa + rough%rt1usb + rt1usc)  - tfrz
    canopy%tscrn = canopy%tscrn3

!   calculating qscrn
    rsts = qsatf(canopy%tscrn1, met%pmb)  ! first tscrn1 used
    qtgnet = rsts * ssoil%wetfac - met%qv
    WHERE (qtgnet .gt. 0. )
       qsurf = rsts * ssoil%wetfac
    ELSEWHERE
       qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
    END WHERE
    qscrn1 =  met%qv - qstar * denom

    rsts = qsatf(canopy%tscrn, met%pmb)  ! I.Harman method
    qtgnet = rsts * ssoil%wetfac - met%qv
    WHERE (qtgnet .gt. 0. )
       qsurf = rsts * ssoil%wetfac
    ELSEWHERE
       qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
    END WHERE
    qscrn3 =  qsurf + (met%qv - qsurf) * r_sc /                            &
        (rough%rt0us + rough%rt1usa + rough%rt1usb + rt1usc) 
    canopy%qscrn = qscrn3

!    print 261,ktau,canopy%tscrn3(1),canopy%tscrn1(1),ssoil%tss(1),met%tk(1),r_sc(1),pnt1(1),rough%rt0us(1), &
!                     rough%rt1usa(1),rough%rt1usb(1),rt1usc(1),canopy%us(1), &
!                    (rough%rt0us(1)+rough%rt1usa(1)+rough%rt1usb(1)+rt1usc(1))/canopy%us(1), &
!                     term1(1),term2(1),term3(1),term5(1),yy1(1),yy2(1),yy3(1),rough%hruff(1), &
!                     rough%disp(1)
261 format(x,'tscrn3new',i5,4f6.1,f8.3,f4.1,x,4f6.1,x,f7.4,x,f7.2,4f6.2,x,3f7.2,x,2f7.4)
    where ( rough%hruff > 0.01 ) &
       beta_sc = vonk /                                                                     &
           ( LOG( (rough%hruff - rough%disp) / rough%z0m )                               &
                 - psim( (rough%hruff - rough%disp) / (rough%zref_uv/canopy%zetar(:,iterplus)) )  &
                 + psim(  rough%z0m / (rough%zref_uv/canopy%zetar(:,iterplus)) )                  &
            ) 
    zscl_10m = max(rough%z0soilsn,10.)
    pnt2 = 0
    where ( rough%hruff <= 1.e-2 ) 
       canopy%ua_10m = canopy%us/vonk *                                                     &
           ( LOG( zscl_10m  / rough%z0m )                                  &
                 - psim( zscl_10m * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim( rough%z0m * canopy%zetar(:,iterplus) / rough%zref_uv )  &
            )
       pnt2 = 1
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff > zscl_10m ) 
       canopy%ua_10m = canopy%us *                                                       &
                       exp( (zscl_10m-rough%hruff)*beta_sc/(rough%disp*vonk) ) / beta_sc
       pnt2 = 2
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff <= zscl_10m ) 
       canopy%ua_10m = canopy%us/vonk *                                                     &
           ( LOG( (zscl_10m - rough%disp) / rough%z0m )                                  &
                 - psim( (zscl_10m - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim(  rough%z0m * canopy%zetar(:,iterplus) / rough%zref_uv )                          &
            ) 
       pnt2 = 3
    endwhere
!   calculation of uscrn at 1.8m
    zscl_scrn = max(rough%z0soilsn,1.8)
    canopy%uscrn1 = max(0., max(met%ua,umin) - canopy%us * denom ) ! for bare ground only
    canopy%uscrn = canopy%uscrn1 ! for bare ground only
    
    where ( rough%hruff <= 1.e-2 )
       canopy%uscrn = canopy%us/vonk *                                                     &
           ( LOG( zscl_scrn  / rough%z0m )                                  &
                 - psim( zscl_scrn * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim( rough%z0m * canopy%zetar(:,iterplus)/ rough%zref_uv )    &
            )
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff > zscl_scrn )
       canopy%uscrn = canopy%us *                                                       &
                       exp( (zscl_scrn-rough%hruff)*beta_sc/(rough%disp*vonk) ) / beta_sc
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff <= zscl_scrn )
       canopy%uscrn = canopy%us/vonk *                                                     &
           ( LOG( (zscl_scrn - rough%disp) / rough%z0m )                                  &
                 - psim( (zscl_scrn - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim(  rough%z0m * canopy%zetar(:,iterplus) / rough%zref_uv )                          &
            )
    endwhere

!    print *,'laihc',canopy%vlaiw,veg%hc,rough%hruff
!    if( .not. L_EXPLICIT ) then
!    print 219,ktau,met%ua(1),canopy%us(1),canopy%ua_10m(1),zscl_10m(1),       &
!                   rough%hruff(1),rough%disp(1),rough%z0m(1),beta_sc(1),pnt2(1), &
!                   met%ua(2),canopy%us(2),canopy%ua_10m(2),zscl_10m(2),       &
!                   rough%hruff(2),rough%disp(2),rough%z0m(2),                &
!                   met%tk(2)-tfrz,pnt2(2),ts1v(2),ts2v(2),ts3v(2),canopy%zetar(:,iterplus)
219 format(x,'newtu10',i4,f6.2,f6.3,f6.2,x,f5.1,2f5.1,x,f5.3,f5.2,f3.0,3x,  &
                          f5.2,f6.3,f6.2,x,f5.1,2f8.5,x,f9.7,f5.1,f3.0,x,3f6.3,2f6.2)
!    print 115,ktau,canopy%tscrn1(1),canopy%tscrn2(1),canopy%tscrn3(1), &
!                   met%tk(1)-tfrz,met%tvair(1)-tfrz,ssoil%tss(1)-tfrz,      &                     
!                   zscrn(1),rough%disp(1),rough%hruff(1),rough%zruffs(1),canopy%vlaiw(1),pnt1(1),   & 
!                   canopy%tscrn1(2),canopy%tscrn2(2),canopy%tscrn3(2), &
!                   met%tk(2)-tfrz,met%tvair(2)-tfrz,ssoil%tss(2)-tfrz,      &
!                   zscrn(2),rough%disp(2),rough%hruff(2),rough%zruffs(2),canopy%vlaiw(2),pnt1(2), &
!                   met%fsd(1,3),met%ua(1)
!                   met%ua(1),canopy%us,canopy%ua_10m,zscrn,zscrn_10m,       &
!                   rough%hruff,rough%disp,rough%z0m
115 format(x,'newtscrn',i4,3f5.1,x,3f5.1,x,4f6.2,x,f4.2,f3.0,3x,3f5.1,x,3f5.1,4f5.2,x,f4.2,f3.0,f5.0,f6.2)
!           f5.1,x,2f6.3,2f8.2,x,4f5.1,x,2f6.2,x,2f5.1,2f9.6) 
!    print 117,ktau,canopy%tscrn(1),canopy%vlaiw(1),rough%hruff(1),rough%disp(1), &
!              xi(1,1),xi(1,2),xi(1,3), &
!              ti(1,1),ti(1,2),ti(1,3),rough%hruff_grmx(1),rough%zref_tq(1),rough%za_tq(1)
117 format(1x,'tscrnew',i5,f7.1,f5.2,x,2f6.1,x,3f6.1,x,3f6.1,2x,3f7.1)
!    print 118,ktau,canopy%tscrn(2),canopy%vlaiw(2),rough%hruff(2),rough%disp(2), &
!              xi(2,1),xi(2,2),xi(2,3), &
!              ti(2,1),ti(2,2),ti(2,3),rough%hruff_grmx(2),rough%zref_tq(2)
118 format(1x,'tscrnbg',i5,f7.1,f5.2,x,2f6.1,x,3f6.1,x,3f6.1,2x,2f7.1)
!    endif ! if( .not. L_EXPLICIT )
       
    ! Calculate radiative/skin temperature;  at this stage old soil temperature is used
!    rad%trad = ( (1.-rad%transd)*canopy%tv**4 + rad%transd * ssoil%tss**4 )**0.25

#if defined(SCMA)
    IF(L_EXPLICIT)THEN
      print 87,ktau,met%tk(1),met%tvair(1),ssoil%tss(1),canopy%tv(1),met%fsd(1,3),met%fld(1), &  !  1-7
               met%ua(1),canopy%us(1),rough%disp(1),zscrn(1),canopy%zetar(1,iterplus), &         !  8-12
               rough%z0m(1),canopy%fe(1),canopy%fh(1),rad%trad(1),canopy%qscrn(1),            &  ! 13-17
               canopy%ga(1),canopy%tscrn1(1)+273.16,canopy%tscrn(1)+273.16,qscrn1(1),qscrn3(1),met%qv(1),qsurf(1)
87    format(1x,'rough%dispveg',i4,4f6.1,x,2f5.0,x,f5.2,f6.3,x,2f5.1,2x,f7.3,x,f9.6,x,2f7.1, &
                 f7.1,f9.5,3f6.1,4f8.4)
      print 871,ktau,met%tk(2),met%tvair(2),ssoil%tss(2),canopy%tv(2),met%fsd(2,3),met%fld(2), &
               met%ua(2),canopy%us(2),rough%disp(2),zscrn(2),canopy%zetar(2,iterplus), &
               rough%z0m(2),canopy%fe(2),canopy%fh(2),rad%trad(2),canopy%qscrn(2), &
               canopy%ga(2),canopy%tscrn1(2)+273.16,canopy%tscrn(2)+273.16,qscrn1(2),qscrn3(2),met%qv(2),qsurf(2)
871    format(1x,'rough%dispbgr',i4,4f6.1,x,2f5.0,x,f5.2,f6.3,x,2f5.1,2x,f7.3,x,f9.6,x,2f7.1, &
                 f7.1,f9.5,3f6.1,4f8.4)
      print 98,ktau,met%tk(1),met%tvair(1),ssoil%tss(1),canopy%tv(1),met%fsd(1,3),met%fld(1), &  !  1-7
               met%ua(1),canopy%us(1),rough%disp(1),zscrn(1),canopy%zetar(1,iterplus), &         !  8-12
               rough%z0m(1),canopy%fe(1),canopy%fh(1),canopy%tscrn1(1)+273.16,canopy%qscrn(1), &  ! 13-17
               canopy%ga(1),canopy%tscrn2(1)+273.16,canopy%tscrn3(1)+tfrz,veg%hc(1), &          ! 18-21
               canopy%vlaiw(1),rough%zref_tq(1),rough%zref_uv(1),canopy%ua_10m(1),rough%zruffs(1), & ! 22-26  
               zscl(1),qscrn1(1),qscrn3(1),met%qv(1),qsurf(1), canopy%uscrn1(1),canopy%uscrn(1)
98    format(1x,'tscu10veg',i4,4f6.1,x,2f5.0,x,f5.2,f6.3,x,2f5.1,2x,f7.3,x,f9.6,x,2f7.1, &
                 f7.1,f9.5,f7.1,2f6.1,f6.1,f5.1,2f5.1,x,f6.2,f6.2,f6.2,4f8.4,2f7.3)
      print 898,ktau,met%tk(2),met%tvair(2),ssoil%tss(2),canopy%tv(2),met%fsd(1,3),met%fld(2), &  !  1-7
               met%ua(2),canopy%us(2),rough%disp(2),zscrn(2),canopy%zetar(2,iterplus), &         !  8-12
               rough%z0m(2),canopy%fe(2),canopy%fh(2),canopy%tscrn1(2)+273.16,canopy%qscrn(2), &  ! 13-17
               canopy%ga(2),canopy%tscrn2(2)+273.16,canopy%tscrn3(2)+tfrz,veg%hc(2), &          ! 18-21
               canopy%vlaiw(2),rough%zref_tq(2),rough%zref_uv(2),canopy%ua_10m(2),rough%zruffs(2), & ! 22-26  
               zscl(2),qscrn1(2),qscrn3(2),met%qv(2),qsurf(2),canopy%uscrn1(2),canopy%uscrn(2)
898    format(1x,'tscu10bg',i4,4f6.1,x,2f5.0,x,f5.2,f6.3,x,2f5.1,2x,f7.3,x,f9.6,x,2f7.1, &
                 f7.1,f9.5,f7.1,2f6.1,f6.1,f5.1,2f5.1,x,f6.2,f6.2,f6.2,4f8.4,2f7.3)
!

!
!      print 88,ktau,canopy%qscrn,qsurf,met%qv,qstar,denom,ssoil%wb(1,1), &
!             ssoil%wb(2,1),canopy%fe 
!      print 89,ktau,canopy%tscrn,rad%trad-273.16,met%tc,ssoil%tss-273.16,canopy%tv-273.16,tstar
!      met%tk(1),rad%rniso(1,:),-capp*rmair*(tlfx -tair2)*rad%gradis,canopy%tv(1),canopy%fns(1),canopy%fh(1), &
      print 199,ktau,canopy%fnv(1),canopy%fhv(1),canopy%fev(1),canopy%fevc(1),canopy%fevw(1), &
                canopy%fnv(1)-canopy%fhv(1)-canopy%fev(1),SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2), &
      met%tk(1),rad%rniso(1,:),canopy%tv(1),canopy%fns(1),canopy%fh(1), &
                    canopy%fes(1),canopy%fns(1)-canopy%fes(1)-canopy%fhs(1),canopy%ga(1),canopy%ghflux(1), &
                    met%fsd(1,3),rad%qssabs(1),ssoil%tss(1),SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs(1), &
                    rad%albedo(1,1:2),met%fld(1)
199 format(1x,'checkSEB',i4,5f5.1,f6.2,2x,2f5.1,f6.1,2f5.1,2f6.1,2x,3f5.1,x,2f6.2,2f6.0,f6.1,2f6.0,2f6.3,f5.0)

!      print 90,ktau,canopy%fe,canopy%fes,canopy%fev,canopy%fevc,canopy%fevw,canopy%fh,canopy%fhs,&
!                    canopy%fhv,canopy%fhvw,canopy%us,met%ua/10.
      
    ENDIF
#endif

88  format(x,'screenvarqv',i5,2f8.5,2x,2f8.5,2x,2f8.5,2x,2f8.5,x,2f5.1,2x,2f5.2,2f5.0)
89  format(1x,'screenvartsc',i5,1x,2f7.1,1x,2f7.1,1x,2f7.1,4f7.1,2x,2f6.2)
90  format(1x,'canopy ls',i5,18f8.3,1x,4f7.3)


    avgwrs = sum(veg%froot * ssoil%wb,2)
    avgtrs = max(0.0,sum(veg%froot * ssoil%tgg,2)-tfrz)
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
!les    !frs is currently calculated in carbon module.
    canopy%frs  = rsoil(soil%rs20, avgwrs, avgtrs)
    canopy%frs  = canopy%frs &
         * sum(spread(bgc%ratecs,1, mp) * bgc%csoil,2)	&
         /(365.0*24.0*3600.0)		     !convert 1/year to 1/second
    WHERE (ssoil%snowd > 1.)
       canopy%frs	= canopy%frs / min(100.,ssoil%snowd)
    END WHERE
!    canopy%frday = 12.0 * sum(rdy, 2)
!    canopy%fpn = -12.0 * sum(an_y, 2)
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp

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
!les canopy%precis = canopy%through
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

!    bal%drybal=REAL(ecy(:,1)+hcy(:,1),r_1)-rad%rniso(:,1) &
!         +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*rad%gradis(:,1)
!
!    bal%drybal=REAL(ecy(:,2)+hcy(:,2),r_1)-rad%rniso(:,2) &
!         +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*rad%gradis(:,2)
    bal%drybal=REAL(ecy+hcy,r_1)-SUM(rad%rniso,2) &
         +capp*rmair*(tlfy-met%tvair)*SUM(rad%gradis,2)
    !ypw: energy balance of the wet canopy
    bal%wetbal=canopy%fevw+canopy%fhvw-SUM(rad%rniso,2)*canopy%fwet &
         +canopy%fhvw*SUM(rad%gradis,2)/MAX(0.001,ghwet)

    DEALLOCATE(cansat,ghwet,gbhu)
    DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
    DEALLOCATE(ecy, hcy, rny)
    DEALLOCATE(gbhf, gswmin, csx)
!    DEALLOCATE(gbhf, gswmin, gswx, csx)

  END SUBROUTINE define_canopy

!  CONTAINS
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
      REAL(r_1), DIMENSION(mp), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp)			:: z
      z = q10c4 ** (0.1 * x - 2.5) / &
           ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
    END FUNCTION xvcmxt4
    !---------------------------------------------------------
    FUNCTION xvcmxt3(x) RESULT(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for vcmax for c3 plants
      REAL(r_1), DIMENSION(mp), INTENT(IN) :: x
      REAL(r_1), PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
      REAL(r_1), PARAMETER  :: xVccoef = 1.17461 ! derived parameter
                        ! xVccoef=1.0+exp((EntropJx*TrefK-EHdJx)/(Rconst*TrefK))
      REAL(r_1), DIMENSION(mp)  :: xvcnum
      REAL(r_1), DIMENSION(mp)  :: xvcden
      REAL(r_1), DIMENSION(mp)  :: z
      xvcnum=xVccoef*exp((EHaVc/(rgas*trefk))*(1.0-trefk/x))
      xvcden=1.0+exp((EntropVc*x-EHdVc)/(rgas*x))
      z = MAX(0.0,xvcnum/xvcden)
    END FUNCTION xvcmxt3
    !---------------------------------------------------------
    FUNCTION xejmxt3(x) RESULT(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants
      REAL(r_1), DIMENSION(mp), INTENT(IN) :: x
      REAL(r_1), PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
      REAL(r_1), PARAMETER  :: xjxcoef = 1.16715 ! derived parameter
      REAL(r_1), DIMENSION(mp)  :: xjxnum
      REAL(r_1), DIMENSION(mp)  :: xjxden
      REAL(r_1), DIMENSION(mp)  :: z
      xjxnum=xjxcoef*exp((EHaJx/(rgas*trefk))*(1.0-trefk/x))
      xjxden=1.0+exp((EntropJx*x-EHdJx)/(rgas*x))
      z = MAX(0.0, xjxnum/xjxden)
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
  SUBROUTINE dryLeaf(dels,phenps)
    REAL(r_1), INTENT(IN)     :: dels ! integration time step (s)
    REAL(r_1), PARAMETER  :: co2cp3 = 0.0 ! CO2 compensation pt C3
    REAL(r_1), PARAMETER  :: jtomol = 4.6e-6 ! Convert from J to Mol for light
    REAL(r_1), DIMENSION(mp)  :: conkct ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp)  :: conkot ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp)  :: cx1  ! "d_{3}" in Wang and Leuning,
    REAL(r_1), DIMENSION(mp)  :: cx2  !     1998, appendix E
    REAL(r_1), DIMENSION(mp)  :: tdiff ! leaf air temp diff.
    REAL(r_1), DIMENSION(mp)  :: tlfxx ! leaf temp of current iteration (K)
    REAL(r_1), DIMENSION(mp)  :: abs_deltlf ! ABS(deltlf)
    REAL(r_1), DIMENSION(mp)  :: deltlf ! deltlfy of prev iter.
    REAL(r_1), DIMENSION(mp)  :: deltlfy ! del temp successive iter.
    REAL(r_1), DIMENSION(mp)  :: gras ! Grashof coeff
    REAL(r_1), DIMENSION(mp)  :: evapfb !
    REAL(r_1), DIMENSION(mp)  :: temp
    REAL(r_2), DIMENSION(mp)  :: ecx ! lat. hflux big leaf
    REAL(r_2), DIMENSION(mp)  :: hcx ! sens heat fl big leaf prev iteration
    REAL(r_2), DIMENSION(mp)  :: rnx ! net rad prev timestep
    REAL(r_1), DIMENSION(mp,ms)  :: evapfbl !
    REAL(r_1), DIMENSION(mp,mf)  :: gw  ! cond for water for a dry canopy
    REAL(r_1), DIMENSION(mp,mf)  :: gh  ! cond for heat for a dry canopy
    REAL(r_1), DIMENSION(mp,mf)  :: ghr ! dry canopy cond for heat & thermal rad
    REAL(r_1), DIMENSION(mp,mf)  :: anx ! net photos. prev iteration
    REAL(r_1), DIMENSION(mp,mf)  :: an_y ! net photosynthesis soln
    REAL(r_1), DIMENSION(mp,mf)  :: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_1), DIMENSION(mp,mf)  :: rdy ! daytime leaf resp rate
    REAL(r_1), DIMENSION(mp,mf)  :: ejmax2 ! jmax of big leaf
    REAL(r_1), DIMENSION(mp,mf)  :: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_1), DIMENSION(mp,mf)  :: vcmxt3 ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp,mf)  :: vcmxt4 ! vcmax big leaf C4
    REAL(r_1), DIMENSION(mp,mf)  :: vx3 ! carboxylation C3 plants
    REAL(r_1), DIMENSION(mp,mf)  :: vx4 ! carboxylation C4 plants
    REAL(r_1), DIMENSION(mp,mf)  :: xleuning ! leuning stomatal coeff
    REAL(r_1), DIMENSION(mp,mf)  :: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp,mf)  :: temp2
    REAL(r_1), dimension(mp)     :: phenps
    INTEGER(i_d)   :: k  ! iteration count
    INTEGER(i_d)   :: kk  ! iteration count

    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    rdx = 0.0
    rnx = SUM(rad%rniso,2)
    abs_deltlf = 999.0
    DO kk=1,mp
      IF(canopy%vlaiw(kk) <=1.0e-2) THEN
          abs_deltlf(kk)=0.0
          hcx(kk) = 0.0 ! intialise
          ecx(kk) = 0.0 ! intialise
          anx(kk,:) = 0.0 ! intialise
          rny(kk) = rnx(kk) ! store initial values
          hcy(kk) = hcx(kk) ! store initial values
          ecy(kk) = ecx(kk) ! store initial values
          rdy(kk,:) = rdx(kk,:) ! store initial values
          an_y(kk,:) = anx(kk,:) ! store initial values
       END IF
    ENDDO
    deltlfy = abs_deltlf
    k = 0
    DO WHILE (ANY(abs_deltlf > 0.1)  .AND.  k < maxiter)
       k = k + 1
       ! Grashof number (Leuning et al, 1995) eq E4:
       gras = MAX(1.0e-6, 1.595E8*ABS(tlfx-met%tvair)*(veg%dleaf**3.0))
       ! See Appendix E in (Leuning et al, 1995):
       gbhf(:,1) = rad%fvlai(:,1) * air%cmolar * 0.5*dheat *(gras(:)**0.25) &
                                  / veg%dleaf
       gbhf(:,2) = rad%fvlai(:,2) * air%cmolar * 0.5*dheat *(gras(:)**0.25) &
                                  / veg%dleaf
       gbhf = max(1.e-6,gbhf)
       ! Conductance for heat:
       gh(:,:) = 2.0 * (gbhu(:,:) + gbhf(:,:))
       ! Conductance for heat and longwave radiation:
       ghr(:,:) = rad%gradis(:,:)+gh(:,:)
       !  Leuning 2002 (P C & E) equation for temperature response
       !  used for Vcmax for C3 plants:
       temp =  phenps * xvcmxt3(tlfx) * veg%vcmax * (1.0-veg%frac4)
       vcmxt3 =  rad%scalex * SPREAD(temp, 2, mf)
       ! Temperature response of Vcmax for C4 plants (Collatz et al 1989):
       temp = phenps * xvcmxt4(tlfx-tfrz) * veg%vcmax * veg%frac4
       vcmxt4 =  rad%scalex * SPREAD(temp, 2, mf)
       !  Leuning 2002 (P C & E) equation for temperature response
       !  used for Jmax for C3 plants:
       temp = phenps * xejmxt3(tlfx) * veg%ejmax * (1.0-veg%frac4)
       ejmxt3 =  rad%scalex * SPREAD(temp, 2, mf)
       ! Difference between leaf temperature and reference temperature:
       tdiff = tlfx - trefk
       ! Michaelis menten constant of Rubisco for CO2:
       conkct = conkc0 * EXP((ekc/(rgas*trefk)) * (1.0-trefk/tlfx))
       ! Michaelis menten constant of Rubisco for oxygen:
       conkot = conko0 * EXP((eko/(rgas*trefk)) * (1.0-trefk/tlfx))
       ! "d_{3}" in Wang and Leuning, 1998, appendix E:
       cx1 = conkct * (1.0+0.21/conkot)
       cx2 = 2.0 * gam0 * (1.0 + gam1*tdiff + gam2*tdiff*tdiff)
       ! All equations below in appendix E in Wang and Leuning 1998 are
       ! for calculating anx, csx and gswx for Rubisco limited, RuBP limited,
       ! sink limited
       temp2(:,1) = rad%qcan(:,1,1) * jtomol * (1.0-veg%frac4)
       temp2(:,2) = rad%qcan(:,2,1) * jtomol * (1.0-veg%frac4)
       vx3(:,1)  = ej3x(temp2(:,1),ejmxt3(:,1))
       vx3(:,2)  = ej3x(temp2(:,2),ejmxt3(:,2))


       temp2(:,1) = rad%qcan(:,1,1) * jtomol * veg%frac4
       temp2(:,2) = rad%qcan(:,2,1) * jtomol * veg%frac4
       vx4(:,1)  = ej4x(temp2(:,1),vcmxt4(:,1))
       vx4(:,2)  = ej4x(temp2(:,2),vcmxt4(:,2))

       rdx(:,1) = (cfrd3*vcmxt3(:,1) + cfrd4*vcmxt4(:,1))*fwsoil
       rdx(:,2) = (cfrd3*vcmxt3(:,2) + cfrd4*vcmxt4(:,2))*fwsoil
       xleuning(:,1) = (fwsoil / (csx(:,1)-co2cp3))  &
                     * ((1.0-veg%frac4) * a1c3 / (1.0+dsx/d0c3) &
                         + veg%frac4    * a1c4 / (1.0+dsx/d0c4))
       xleuning(:,2) = (fwsoil / (csx(:,2)-co2cp3))  &
                     * ((1.0-veg%frac4) * a1c3 / (1.0+dsx/d0c3) &
                         + veg%frac4    * a1c4 / (1.0+dsx/d0c4))

!       print *,'befph',csx(:,1),cx1(:),cx2(:),gswmin(:,1),rdx(:,1),&
!                     vcmxt3(:,1),vcmxt4(:,1),vx3(:,1),vx4(:,1),xleuning(:,1)
       anx(:,1) = photosynthesis(csx(:,1),cx1(:),cx2(:),gswmin(:,1),rdx(:,1),&
                     vcmxt3(:,1),vcmxt4(:,1),vx3(:,1),vx4(:,1),xleuning(:,1))

       anx(:,2) = photosynthesis(csx(:,2),cx1(:),cx2(:),gswmin(:,2),rdx(:,2), &
                     vcmxt3(:,2),vcmxt4(:,2),vx3(:,2),vx4(:,2),xleuning(:,2))

!       write(79,*) 'within while loop=',anx(:,1)/1.2e-5,rdx(:,1)/1.2e-5,anx(:,2)/1.2e-5,rdx(:,2)/1.2e-5
       csx = SPREAD(met%ca, 2, mf) - rgbwc * anx / (gbhu + gbhf)
       !canopy%gswx = gswmin + MAX(0.0, rgswc*xleuning*anx)
       canopy%gswx = MAX(1.e-3, gswmin + MAX(0.0, rgswc*xleuning*anx))
!       print *,'in dryL ,gbhu1',gbhu
!       print *,'in dryL ,gbhu2',canopy%gswx
!       print *,'in dryL ,gbhu3',gbhf,csx
!       print *,'in dryL ,gbhu4',csx
!       print *,'in dryL ,anx',anx
!       print *,'in dryL ,gbhf',gbhf
!       print *,'in dryL ,gw',gw

       ! Recalculate conductance for water:
       gw = 1.0/(1.0/canopy%gswx + 1.0/(1.075*(gbhu+gbhf)))
       gw = max(gw, 0.000001) 
       ! Modified psychrometric constant (Monteith and Unsworth, 1990)
!       print *,'gw',gw
!       print *,'ghr',ghr
!       print *,'psyc',air%psyc
!       print *,'tlfx',tlfx
       psycst = SPREAD(air%psyc, 2, mf) *ghr/gw
       ! Store leaf temperature:
       tlfxx = tlfx
       ! Update canopy latent heat flux:
!       ecx(:) = (air%dsatdk(:)*rad%rniso(:,1) &
!                 + capp*rmair*met%dva(:)*ghr(:,1)) &
!                 / (air%dsatdk(:)+psycst(:,1)) &
!             + (air%dsatdk(:)*rad%rniso(:,2) &
!                 + capp*rmair*met%dva(:)*ghr(:,2)) &
!                 / (air%dsatdk(:)+psycst(:,2))


       ecx(:) = (air%dsatdk(:)*(rad%rniso(:,1)- capp*rmair*(met%tvair(:)-met%tk(:))*rad%gradis(:,1)) &
                 + capp*rmair*met%dva(:)*ghr(:,1)) &
                 / (air%dsatdk(:)+psycst(:,1)) &
             + (air%dsatdk(:)*(rad%rniso(:,2)- capp*rmair*(met%tvair(:)-met%tk(:))*rad%gradis(:,2)) &
                 + capp*rmair*met%dva(:)*ghr(:,2)) &
                 / (air%dsatdk(:)+psycst(:,2))


       evapfbl = 0.0
       evapfb = (1.0-canopy%fwet) * REAL(ecx,r_1) * dels/air%rlam
                      ! convert W/m2 to mm/dt using dels/air%rlam
       DO kk = 1,ms
          WHERE (ecx > 0.0.and.canopy%fwet < 1.0)
             evapfbl(:,kk) = MIN(evapfb*veg%froot(:,kk), &
                  MAX(0.0,MIN(REAL(ssoil%wb(:,kk),r_1)-soil%swilt, &
                  REAL(ssoil%wb(:,kk)-1.05*ssoil%wbice(:,kk),r_1))) &
                  * soil%zse(kk)*1000.0)
             canopy%fevc = 0.0
          END WHERE
       END DO
       canopy%fevc=canopy%fevc+ SUM(evapfbl,2)*air%rlam/dels
       WHERE (ecx > 0.0.and.canopy%fwet < 1.0) &
          ecx = canopy%fevc /(1.0-canopy%fwet)

       ! Update canopy sensible heat flux:
       hcx = (SUM(rad%rniso,2)-ecx-capp*rmair*(met%tvair-met%tk)*SUM(rad%gradis,2))*SUM(gh,2)/SUM(ghr,2)
       ! Update leaf temperature:
       tlfx=met%tvair+REAL(hcx,r_1)/(capp*rmair*SUM(gh,2))
       ! Update net radiation for canopy:
!       rnx = SUM(rad%rniso,2) - capp*rmair*(tlfx-met%tvair)*SUM(rad%gradis,2)
       rnx = SUM(rad%rniso,2) - capp*rmair*(tlfx-met%tk)*SUM(rad%gradis,2)
       ! Update leaf surface vapour pressure deficit:
       dsx = met%dva + air%dsatdk * (tlfx-met%tvair)
       ! Store change in leaf temperature between successive iterations:
       deltlf = tlfxx-tlfx
       abs_deltlf = ABS(deltlf)
       ! Where leaf temp change b/w iterations is significant, and
       ! difference is smaller than the previous iteration, store results:
       WHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
          deltlfy = deltlf
          tlfy = tlfx
          rny = rnx
          hcy = hcx
          ecy = ecx
          rdy(:,1) = rdx(:,1)
          rdy(:,2) = rdx(:,2)
          an_y(:,1) = anx(:,1)
          an_y(:,2) = anx(:,2)
       END WHERE
       WHERE (abs_deltlf > 0.1)
       ! after 4 iterations, take mean value of current & previous estimates
       ! as the next estimate of leaf temperature, to avoid oscillation
          tlfx = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx + &
               (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx
       END WHERE
       IF(k==1) THEN
       ! take the first iterated estimates as the defaults
          tlfy = tlfx
          rny = rnx
          hcy = hcx
          ecy = ecx
          rdy = rdx
          an_y = anx
       END IF
       !print *,"fpn ",k,-12.*sum(an_y(1301,:)),tlfy(1301)
    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < maxiter)
    ! dry canopy flux
    canopy%fevc = (1.0-canopy%fwet) * ecy

    canopy%frday = 12.0 * SUM(rdy, 2)
    canopy%fpn = -12.0 * SUM(an_y, 2)

  END SUBROUTINE dryLeaf
    !---------------------------------------------------------
    FUNCTION photosynthesis(csxz,cx1z,cx2z,gswminz,rdxz,vcmxt3z,vcmxt4z, &
                            vx3z,vx4z,xleuningz) RESULT(z)

    ! inputs:
    REAL(r_1), DIMENSION(mp)  :: csxz
    REAL(r_1), DIMENSION(mp)  :: cx1z
    REAL(r_1), DIMENSION(mp)  :: cx2z
    REAL(r_1), DIMENSION(mp)  :: gswminz
    REAL(r_1), DIMENSION(mp)  :: rdxz
    REAL(r_1), DIMENSION(mp)  :: vcmxt3z
    REAL(r_1), DIMENSION(mp)  :: vcmxt4z
    REAL(r_1), DIMENSION(mp)  :: vx4z
    REAL(r_1), DIMENSION(mp)  :: vx3z
    REAL(r_1), DIMENSION(mp)  :: xleuningz
    !local variables
    REAL(r_2), DIMENSION(mp)  :: coef0z,coef1z,coef2z
    REAL(r_2), DIMENSION(mp)  :: ciz,delcxz
    REAL(r_2), DIMENSION(mp)  :: anrubiscoz,anrubpz,ansinkz
    REAL(r_1), DIMENSION(mp)  :: z
    REAL(r_1), PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
                                             ! Bonan,LSM version 1.0, p106)
    ! rgswc inherited from canopy_module's USE photosynthetic_constants
    ! mp    inherited from canopy_module
!    print *,'inph',mp,csxz,cx1z,cx2z,gswminz,rdxz,vcmxt3z,vcmxt4z, &
!                            vx3z,vx4z,xleuningz

   ! Rubisco limited:
!     print *, 'coef2zd',gswminz,rgswc,xleuningz,vcmxt3z,rdxz,vcmxt4z

     coef2z = gswminz/rgswc+xleuningz *(vcmxt3z-(rdxz-vcmxt4z))
     coef1z = (1.0-csxz*xleuningz) *(vcmxt3z+vcmxt4z-rdxz)      &
                   + (gswminz/rgswc)*(cx1z-csxz) -xleuningz*(vcmxt3z*cx2z/2.0 &
                   + cx1z*(rdxz-vcmxt4z))
     coef0z = -(1.0-csxz*xleuningz) *(vcmxt3z*cx2z/2.0  &
                   + cx1z*(rdxz-vcmxt4z)) -(gswminz/rgswc)*cx1z*csxz
!     delcxz = coef1z**2 -4.0*coef0z*coef2z
!     ciz    = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*max(1.e-6,coef2z))
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) < 1.e-9)
       ! no solution, give it a huge number
       ciz = 99999.0
     END WHERE
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) >= 1.e-9)
       ! solve linearly
       ciz = -1.0 * coef0z / coef1z
     END WHERE
     WHERE (ABS(coef2z) >= 1.e-9)
       ! solve quadratic 
       delcxz = coef1z**2 -4.0*coef0z*coef2z
       ciz    = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*coef2z)
     END WHERE
     ciz    = MAX(0.0_r_2,ciz)
     anrubiscoz = vcmxt3z*(ciz-cx2z/2.0) / (ciz + cx1z) + vcmxt4z - rdxz
          
!     print *,'anr',anrubiscoz
   ! RuBP limited:
     coef2z = gswminz/rgswc+xleuningz *(vx3z-(rdxz-vx4z))
     coef1z = (1.0-csxz*xleuningz) *(vx3z+vx4z-rdxz)    &
                   + (gswminz/rgswc)*(cx2z-csxz) -xleuningz*(vx3z*cx2z/2.0 &
                   + cx2z*(rdxz-vx4z))
     coef0z = -(1.0-csxz*xleuningz) *(vx3z*cx2z/2.0  &
                   + cx2z*(rdxz-vx4z)) -(gswminz/rgswc)*cx2z*csxz
!     print *,'coef012',coef2z,coef1z,coef0z
     !delcxz = coef1z**2 -4.0*coef0z*coef2z
     !ciz    = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*max(1.e-6,coef2z))
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) < 1.e-9)
       ! no solution, give it a huge number
       ciz = 99999.0
     END WHERE
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) >= 1.e-9)
       ! solve linearly
       ciz = -1.0 * coef0z / coef1z
     END WHERE
     WHERE (ABS(coef2z) >= 1.e-9)
       ! solve quadratic 
       delcxz = coef1z**2 -4.0*coef0z*coef2z
       ciz    = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*coef2z)
     END WHERE
     ciz    = MAX(0.0_r_2,ciz)
     anrubpz  = vx3z*(ciz-cx2z/2.0) /(ciz+cx2z) +vx4z-rdxz
!     print *,'ciz and anr',ciz,anrubpz
   ! Sink limited:
     coef2z = xleuningz
     coef1z = gswminz/rgswc + xleuningz * (rdxz - 0.5*vcmxt3z)  +  &
                     effc4 * vcmxt4z - xleuningz * csxz * effc4 *vcmxt4z
     coef0z = -(gswminz/rgswc)*csxz *effc4*vcmxt4z +    &
                    (rdxz -0.5*vcmxt3z)*gswminz/rgswc
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) < 1.e-9)
       ! no solution, give it a huge number
       ciz = 99999.0
     END WHERE
     WHERE (ABS(coef2z) < 1.e-9 .AND. ABS(coef1z) >= 1.e-9)
       ! solve linearly
       ciz = -1.0 * coef0z / coef1z
     END WHERE
     WHERE (ABS(coef2z) >= 1.e-9)
       ! solve quadratic 
       delcxz = coef1z**2 -4.0*coef0z*coef2z
       ciz    = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*coef2z)
     END WHERE
     ciz    = MAX(0.0_r_2,ciz)
     ansinkz = ciz
     !ansinkz  = (-coef1z+SQRT(MAX(0.0_r_2,delcxz))) /(2.0*max(1.e-6,coef2z))
!     print *,'sinl',coef2z,coef0z,delcxz,ansinkz
   ! minimal of three limited rates
     z    = REAL(MIN(anrubiscoz,anrubpz,ansinkz),r_1)
     !print *,"rubiscoz,rubpz,sinkz ",anrubiscoz(1301),anrubpz(1301),ansinkz(1301)
!     print *,'z',z
    END FUNCTION photosynthesis

    !---------------------------------------------------------

  SUBROUTINE wetLeaf(dels)
   ! assuming the temperature of wet leaf is equal that of dry leaf ="tlfy"
    REAL(r_1), INTENT(IN)     :: dels ! integration time step (s)
    REAL(r_1), DIMENSION(mp)  :: ccfevw ! limitation term for
                                        ! wet canopy evaporation rate
    REAL(r_1), DIMENSION(mp)  :: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp)  :: ghrwet ! wet canopy cond: heat & thermal rad
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    WHERE (canopy%vlaiw > 0.01)
    ! VEG SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
    ! calculate total thermal resistance, rthv in s/m
       ghwet = 2.0   * SUM((gbhu+gbhf),2)
       gwwet = 1.075 * SUM((gbhu+gbhf),2)
       ghrwet = SUM(rad%gradis,2) + ghwet
       ! Calculate fraction of canopy which is wet:
       canopy%fwet = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))
       ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
       ! to avoid excessive evaporation:
       ccfevw = MIN(canopy%cansto * air%rlam / dels, &
                    2.0 / (1440.0 / (dels/60.0)) * air%rlam)

       canopy%fevw = MIN(canopy%fwet * (air%dsatdk * (SUM(rad%rniso,2)- capp*rmair*(met%tvair(:)-met%tk(:))*sum(rad%gradis,2)) &
                                         + capp*rmair*met%dva*ghrwet) &
                         / (air%dsatdk+air%psyc*ghrwet/gwwet), ccfevw)

       ! Calculate sens heat from wet canopy:
       !canopy%fhvw = (canopy%fwet*(SUM(rad%rniso,2)-capp*rmair*(met%tvair(:)-met%tk(:))*sum(rad%gradis,2))-canopy%fevw)*ghwet/ghrwet

       canopy%fhvw = canopy%fwet*(SUM(rad%rniso,2)-capp*rmair*(tlfy-met%tk(:))*sum(rad%gradis,2)) &
                    -canopy%fevw
    END WHERE
  END SUBROUTINE wetLeaf
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
