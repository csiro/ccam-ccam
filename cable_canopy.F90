#include "cable_directives.h"

MODULE canopy_module
  USE photosynthetic_constants
  USE radiation_module
  USE roughness_module
  USE air_module
  USE define_types
  USE physical_constants
  USE cable_common_module

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
  REAL(r_2), DIMENSION(:,:), POINTER :: csx ! leaf surface CO2 concentration
  PRIVATE
  PUBLIC define_canopy, sinbet
CONTAINS
  
  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg,bgc,canopy,L_EXPLICIT)
    TYPE (balances_type),INTENT(INOUT)  :: bal
    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (met_type), INTENT(INOUT)      :: met
    REAL(r_1), INTENT(IN)               :: dels ! integration time setp (s)
    TYPE (soil_snow_type), INTENT(INOUT):: ssoil
    TYPE (bgc_pool_type),INTENT(IN)     :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    LOGICAL, INTENT(IN)                 :: L_EXPLICIT
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp)            :: avgtrs !root weighted mean soil temperature
    REAL(r_1), DIMENSION(mp)            :: avgwrs !root weighted mean soil moisture
    REAL(r_1), DIMENSION(mp)            :: dq ! sat spec hum diff.
    REAL(r_1), DIMENSION(mp,mf)         :: frac42       ! 2D frac4
    REAL(r_1), DIMENSION(mp)            :: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp,mf)         :: gw  ! cond for water for a dry canopy

    INTEGER(i_d)                        :: iter ! iteration #
    INTEGER(i_d)                        :: iterplus !
    INTEGER(i_d)                        :: k            ! interation count
    INTEGER(i_d)                        :: kk           ! interation count
    REAL(r_1), DIMENSION(mp)            :: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp)            :: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp)            :: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL(r_1), DIMENSION(mp)            :: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp)            :: rwaters !EAK,09/10
    REAL(r_1), DIMENSION(mp)            :: cc ! limitation term for canopy interception per timestep               
    REAL(r_1), DIMENSION(mp)            :: ccfevw ! limitation term for wet canopy evaporation rate  
    REAL(r_1), DIMENSION(mp)            :: pwet
    REAL(r_1), DIMENSION(mp)            :: denom ! denominator in calculating screen temperature, humidity etc
    REAL(r_1), DIMENSION(mp)            :: tstar ! 
    REAL(r_1), DIMENSION(mp)            :: zscrn !
    REAL(r_1), DIMENSION(mp)            :: qstar !
    REAL(r_1), DIMENSION(mp)            :: rsts  !
    REAL(r_1), DIMENSION(mp)            :: qsurf !
    REAL(r_1), DIMENSION(mp)            :: qtgnet !
    REAL(r_1), DIMENSION(mp)            :: phenps ! Leaf phenology influence on vcmax and jmax
    REAL(r_1), DIMENSION(mp)            :: poolcoef1 ! leaf carbon turnover rate * leaf pool size
    REAL(r_1), DIMENSION(mp)            :: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp)            :: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp)            :: rbw  ! leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp)            :: rrbw ! recipr. leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp)            :: rsw  ! stomatal resistance for water
    REAL(r_1), DIMENSION(mp)            :: rrsw ! recipr. stomatal resistance for water
    REAL(r_1), DIMENSION(mp)            :: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp)            :: tss4 ! soil/snow temperature**4
    REAL(r_1), DIMENSION(mp)            :: sss ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)            :: cc1 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)            :: cc2 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)            :: cc1T ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)            :: cc2T ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp)            :: qstvair ! sat spec hunidity at leaf temperature
!    REAL(r_1), DIMENSION(mp)            :: qstss ! sat spec hunidity at soil/snow temperature
    REAL(r_1), DIMENSION(mp)            :: xx ! delta-type function for sparse canopy limit, p20 SCAM manual

!%% changes by Ashok Luhar (low wind speed)
    REAL(r_1), PARAMETER                :: alpha1=4.0
    REAL(r_1), PARAMETER                :: beta1=0.5
    REAL(r_1), PARAMETER                :: gamma1=0.3
    REAL(r_1), DIMENSION(mp)            :: zeta1
    REAL(r_1), DIMENSION(mp)            :: zeta2
    REAL(r_1), DIMENSION(mp)            :: usA
    REAL(r_1), DIMENSION(mp,3)          :: xi
    REAL(r_1), DIMENSION(mp,3)          :: ti
    REAL(r_1), DIMENSION(mp,3)          :: si
    REAL(r_1), DIMENSION(mp)            :: term1
    REAL(r_1), DIMENSION(mp)            :: term2
    REAL(r_1), DIMENSION(mp)            :: term3
    REAL(r_1), DIMENSION(mp)            :: term5
    REAL(r_1), DIMENSION(mp)            :: r_sc
    REAL(r_1), DIMENSION(mp)            :: beta_sc
    REAL(r_1), DIMENSION(mp)            :: zscrn_10m
    REAL(r_1), DIMENSION(mp)            :: zscl
    REAL(r_1), DIMENSION(mp)            :: zscl_scrn
    REAL(r_1), DIMENSION(mp)       :: rghlai

    INTEGER(i_d) :: idjd1,idjd2,idjd3
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    CHARACTER(len=3) :: DIAG_SOIL_RESP = FDIAG_SOIL_RESP
    integer, dimension(1) :: temp1, temp2

    ALLOCATE(cansat(mp),ghwet(mp),gbhu(mp,mf))
    ALLOCATE(dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp))
    ALLOCATE(ecy(mp), hcy(mp), rny(mp))
    ALLOCATE(gbhf(mp,mf), gswmin(mp,mf), csx(mp,mf))
      
    ghwet = 1.e-3
    idjd1 = 3419
    idjd2 = 6488
    idjd3 = 9985

     if( ktau_gl .le. 1) canopy%cansto=0.0
    !   xjxcoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    ! 1-oct-2002 ypw: to keep the unit consistent for resistance or conductance
    ! s/m for r; mol/m2/s for g, and always use g where appropriate
    ! replace rh, rhr, rw  with ghdry/ghwet,ghrdry, ghrwet, gwdry, gwwet
    ! Set surface water vapour pressure deficit:

    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    met%da = (qsatf(met%tc,met%pmb) - met%qv ) * rmair/rmh2o * met%pmb * 100.

    ! Soil water limitation on stomatal conductance:
    !EAK, 09/10 - replace linear approx by polynomial fitting
    rwater = MAX(1.0e-4_r_2, &
         SUM(soil%froot * MAX(0.024,MIN(1.0_r_2,ssoil%wb - &
              SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))
    !fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))

    rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)
    xi(:,1) = soil%swilt
    xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
    xi(:,3) = soil%sfc
    ti(:,1) = 0.
    ti(:,2) = 0.9
    ti(:,3) = 1.0
    si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) *  &
              (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))
    si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) *  &
              (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))
    si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) *  &
              (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))
    fwsoil = 1.
    where (rwater < soil%sfc - 0.02)
        fwsoil = max(0.,min(1., ti(:,1)*si(:,1) + &
                        ti(:,2)*si(:,2) + ti(:,3)*si(:,3)))
    endwhere

    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * veg%vlaiw
    phenps = 1.0  !EAK, 09/10

    ! Leaf phenology influence on vcmax and jmax
    !phenps = max (1.0e-4, MIN(1.,1. - ( (veg%tmaxvj - ssoil%tgg(:,4)+tfrz)/ &
    !(veg%tmaxvj - veg%tminvj) )**2 ) )
    !WHERE ( ssoil%tgg(:,4) < (veg%tminvj + tfrz) ) phenps = 0.0
    !WHERE ( ssoil%tgg(:,4) > (veg%tmaxvj + tfrz) ) phenps = 1.0

    ! Set previous time step canopy water storage:
    ! canopy%oldcansto=canopy%cansto
    ! canopy intercepted rainfall rate is limited to avoid excessive direct
    ! canopy evaporation, modified further by timestep requirement (EAK aug08)
    cc =min( met%precip-met%precip_sn, 4./( 1440./( min(dels,1800.)/60.) ) )

    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - canopy%cansto,0.0), cc), 0.0, &
         cc > 0.0  .AND. met%tk > tfrz)  !EAK, 09/10
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_sn+MIN(met%precip-met%precip_sn, &
                     MAX(0.0, met%precip-met%precip_sn - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint

    ! Calculate fraction of canopy which is wet:
    !veg%fwet   = max(0.0,min(0.9,0.8*canopy%cansto/MAX(cansat,0.01)))
    veg%fwet   = max(0.0,min(0.9,0.8*canopy%cansto/MAX(cansat,0.1)))

    !EAK, 09/10 - used new rwaters calculation
    ssoil%wetfac = MAX(0.0_r_2, MIN(1.0_r_2, &      
       (ssoil%wb(:,1) - soil%swilt) / (soil%sfc - soil%swilt))) 
       !(ssoil%wb(:,1) - soil%swilt/3.0) / (soil%sfc - soil%swilt))) 
!    rwaters = MAX(0.0_r_2,MIN(1.0_r_2, &
!         (ssoil%wb(:,1) - soil%swilt/2.0) / (soil%sfc - soil%swilt/2.)))

!    rwaters = soil%swilt/2. + rwaters * (soil%sfc-soil%swilt/2.)
!    xi(:,1) = soil%swilt/2.
!    xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
!    xi(:,3) = soil%sfc
!    ti(:,1) = 0.
!    ti(:,2) = 0.8
!    ti(:,3) = 1.0
!    si(:,1) = (rwaters - xi(:,2)) / (xi(:,1) - xi(:,2)) *  &
!              (rwaters - xi(:,3)) / (xi(:,1) - xi(:,3))
!    si(:,2) = (rwaters - xi(:,1)) / (xi(:,2) - xi(:,1)) *  &
!              (rwaters - xi(:,3)) / (xi(:,2) - xi(:,3))
!    si(:,3) = (rwaters - xi(:,1)) / (xi(:,3) - xi(:,1)) *  &
!              (rwaters - xi(:,2)) / (xi(:,3) - xi(:,2))
!    ssoil%wetfac = 1.
!    where (rwaters < soil%sfc - 0.02)
!        ssoil%wetfac = max(0.,min(1., ti(:,1)*si(:,1) + &
!                      ti(:,2)*si(:,2) + ti(:,3)*si(:,3)))
!    endwhere
!    where (rwaters < soil%swilt/2.) ssoil%wetfac = 0.0

    ! Temporay fixer for accounting of reduction of soil evaporation 
    ! due to freezing
    where ( ssoil%wbice(:,1) > 0. )
       ssoil%wetfac = ssoil%wetfac*max(0.5,(1.-(ssoil%wbice(:,1)/ssoil%wb(:,1))**2))
    endwhere
    ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)

    canopy%zetar(:,1) = zeta0 ! stability correction terms
    canopy%zetar(:,2) = zetpos + 1 

    ! weight min stomatal conductance by C3 an C4 plant fractions
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants

    canopy%gswx = 1e-3     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    gw = 1.0e-3 ! default values of conductance

    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    canopy%tv = met%tvair
    canopy%fevw_pot = 0.0

    CALL define_air (met, air)

    qstvair = qsatf((met%tvair-tfrz),met%pmb)
    met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.0
    dsx = met%dva     ! init. leaf surface vpd
    ortsoil = ssoil%rtsoil
    canopy%fes = 0.
    canopy%fess = 0.
    canopy%fesp = 0.
    canopy%fevw = 0.
    canopy%fevc = 0.
    ssoil%evapfbl = 0.0
    tss4 = ssoil%tss**4
    tlfx  = met%tk 
    tlfy  = met%tk 
    rt0 = 10.
    rough%rt1 = 10.

    CALL radiation(ktau,bal, soil, ssoil, veg, air, met, rad)

    DO iter = 1, niter

!       CALL define_air (met, air)

!       CALL radiation(ktau,bal, soil, ssoil, veg, air, met, rad)

       gswmin = max(1.e-6,rad%scalex * (gsw03 * (1. - frac42) + gsw04 * frac42))

       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX(1.e-6, &
            vonk * MAX(met%ua,umin) / ( &
            LOG(rough%zref_uv / rough%z0m) - &
            psim(canopy%zetar(:,iter)) + &
            psim(canopy%zetar(:,iter) * rough%z0m / rough%zref_uv) ))

!%%change by Ashok Luhar - low wind formulation
!            usA = 0.0
!        where (canopy%zetar(:,iter) > 0.7 .and. ssoil%snowd < 0.01 )
!            zeta1=canopy%zetar(:,iter) * rough%z0m / rough%zref_uv
!!            usA = MAX(1.e-6, &
!            canopy%us = MAX(1.e-6, &
!            vonk * MAX(met%ua,umin) / ( &
!            alpha1* ((canopy%zetar(:,iter)**beta1*  &
!               (1.0+gamma1*canopy%zetar(:,iter)**(1.0-beta1)))  &       
!             - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
!        endwhere         
!%%

       ! Turbulent aerodynamic resistance from roughness sublayer depth to reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
!              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * (LOG(rough%zref_tq/MAX(rough%zruffs-rough%disp, rough%z0soilsn)) &
        - psis( canopy%zetar(:,iter) ) &
        + psis( canopy%zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn))/rough%zref_tq ) &
          )/vonk

       !rt0 = max(25.,rough%rt0us / canopy%us)
       rt0 = max(5.,rough%rt0us / canopy%us)
!        rt0 = min(rt0,1000.0)

       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       !rough%rt1 = max(25.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
!        rough%rt1 = min(rough%rt1,1000.0)

       WHERE (ssoil%snowd > 0.1)
          ssoil%wetfac = 1.00
       END WHERE
       ssoil%rtsoil = rt0 + rough%rt1*(0.5+sign(0.5,0.01-veg%vlaiw)) 
       ssoil%rtsoil = max(5.,ssoil%rtsoil)   
!       if( ktau <= 12 ) ssoil%rtsoil = max(25.,ssoil%rtsoil) 
!       ssoil%rtsoil = max(1.,ssoil%rtsoil)   

       WHERE ( ssoil%rtsoil .GT. 2.* ortsoil .OR. ssoil%rtsoil .LT. 0.5*ortsoil )
          ssoil%rtsoil = MAX(5.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf *       &
            (canopy%us / MIN(rough%usuh, 0.2) * &
            veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
       ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
       gbhu(:,1) = gbvtop*(1.0-EXP(-veg%vlaiw*(0.5*rough%coexp+rad%extkb))) / &
                                               (rad%extkb+0.5*rough%coexp)
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
                            (1.0-EXP(-0.5*rough%coexp*veg%vlaiw))-gbhu(:,1)
       rny = SUM(rad%rniso,2) ! init current estimate net rad
       hcy = 0.0              ! init current estimate lat heat
       ecy = rny - hcy        ! init current estimate lat heat

       CALL dryLeaf(dels,phenps,rad,rough,air,met,veg,canopy)

       CALL wetLeaf(dels,rad,rough,air,met,veg,canopy)

!       if ((minval(tlfy).lt.200.or.maxval(tlfy).gt.360).and.veg%vlaiw > 0.01) then
!           write(6,*) 'rmlcheck tlfy: ',minval(tlfy),maxval(tlfy)
!           temp1 = minloc(tlfy); temp2 = maxloc(tlfy)
!           write(6,*) 'rmlcheck: min/max loc',minloc(tlfy),maxloc(tlfy)
!           write(6,*) 'rmlcheck: met/tvair',met%tvair(temp1(1)), met%tvair(temp2(1))
!           write(6,*) 'rmlcheck: veg type',veg%iveg(temp1(1)), veg%iveg(temp2(1))
!           write(6,*) 'rmlcheck: lai ',veg%vlaiw(temp1(1)), veg%vlaiw(temp2(1))
!       stop 'Error in cable_canopy: tlf'
!       endif


       ! Calculate latent heat from vegetation:

       canopy%fev = canopy%fevc + canopy%fevw
       ! recalculate for checking energy balance (YP & Mao, jun08)
       ! Calculate sensible heat from vegetation:
       canopy%fhv = (1.0 - veg%fwet) *  REAL(hcy,r_1) + canopy%fhvw
       ! Calculate net rad absorbed by canopy:
       canopy%fnv = (1.0-veg%fwet)*REAL(rny,r_1)+canopy%fevw+canopy%fhvw
       ! canopy radiative temperature is calculated based on long-wave radiation balance
       WHERE (veg%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          rad%lwabv = capp * rmair * (tlfy-met%tk)* SUM(rad%gradis,2) + &
                      canopy%fhvw*SUM(rad%gradis,2)/max(0.001,ghwet)

          !canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tk**4)**0.25
          canopy%tv = max(rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tk**4,0.)**0.25 ! MJT
       ELSEWHERE ! sparse canopy
           canopy%tv = met%tk
       END WHERE
       where (canopy%tv.lt.met%tk-50.) ! MJT
         canopy%tv=met%tk              ! MJT
       end where                       ! MJT

       ! Calculate ground heat flux:
       ! canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
       ! Saturation specific humidity at soil/snow surface temperature:
       ssoil%qstss = qsatf((ssoil%tss - tfrz),met%pmb)
       ! Spec hum deficit at soil/snow surface:
       dq = ssoil%qstss - met%qv
       !                               excessive dew over snow area
       WHERE (ssoil%snowd > 1.0)
          dq = max( -0.1e-3, dq)
       END WHERE
       ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
       ! Calculate net rad to soil:
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
            sboltz*canopy%tv**4 - emsoil*sboltz* tss4
       ! Penman-Monteith formula
!          sss=air%dsatdk
!          cc1=sss/(sss+air%psyc )
!          cc2=air%psyc /(sss+air%psyc )
!          !ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) +  &
!          ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
!          cc2 * air%rho * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil 

       ! Soil latent heat:
       canopy%fess= ssoil%wetfac * ssoil%potev
       ! Reduce soil evap due to presence of puddle
       pwet = max(0.,min(0.2,ssoil%pudsto/max(1.,ssoil%pudsmx)))
       canopy%fess = canopy%fess * (1.-pwet)

       WHERE (ssoil%snowd < 0.1 .and. canopy%fess .gt. 0. )
            !allow the top soil moistur below wilting i.e. soil%swilt/3.
        !canopy%fess= min(canopy%fess,(max(0._r_2,ssoil%wb(:,1)-soil%swilt)*soil%zse(1) &  
        !       *1000. - ssoil%evapfbl(:,1))*air%rlam / dels)                                                 
        canopy%fess= min(canopy%fess,max(0._r_2,((ssoil%wb(:,1)-soil%swilt)*soil%zse(1) &
                         *1000.) - ssoil%evapfbl(:,1) )*air%rlam / dels)

        canopy%fess = min(canopy%fess,(ssoil%wb(:,1)-ssoil%wbice(:,1))* soil%zse(1) &
               * 1000. * air%rlam / dels)
       END WHERE
       ssoil%cls=1.
       WHERE (ssoil%snowd >= 0.1)
          ssoil%cls = 1.1335
          canopy%fess= min(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
       END WHERE
       ! Calculate soil sensible heat:
       canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
       ! Evaporation form soil puddle
       canopy%fesp = min(ssoil%pudsto/dels*air%rlam,max(pwet*ssoil%potev,0.))
       canopy%fes = canopy%fess + canopy%fesp

       ! temporary fix to Jhan's code by EAK
!       where( canopy%fhs <= -80.0 .and. ktau <= 12 ) canopy%fhs = -80.0
!       where( canopy%fhs > 300.0  .and. ktau <= 12 ) canopy%fhs = 300.0
!       where( canopy%fhv <= -80.0 .and. ktau <= 12 ) canopy%fhv = -80.0
!       where( canopy%fhv > 300.0  .and. ktau <= 12 ) canopy%fhv = 300.0
!       where( canopy%fes <= -10.0 .and. ktau <= 12 ) canopy%fes = -10.0
!       where( canopy%fes > 300.0  .and. ktau <= 12 ) canopy%fes = 300.0
!       where( canopy%fev <= -10.0 .and. ktau <= 12 ) canopy%fev = -10.0
       ! Calculate total latent heat:
       canopy%fe = canopy%fev + canopy%fes
       ! Calculate total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs
       ! Calculate ground heat flux:
       canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
       ! Initialise in-canopy temperature and humidity:
       met%tvair = met%tk
       met%qvair = met%qv

       WHERE (veg%meth > 0 .and. veg%vlaiw > 0.01 .and. &
                 rough%hruff > rough%z0soilsn) 
          !   use the dispersion matrix (DM) to find the air temperature 
          !   and specific humidity 
          !   (Raupach, Finkele and Zhang 1997, pp 17)
          ! leaf boundary layer resistance for water
          rbw = air%cmolar/sum(gbhu+gbhf,2)
          rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT 
          ! leaf stomatal resistance for water
          rsw = air%cmolar/sum(canopy%gswx,2)
          rrsw = sum(canopy%gswx,2)/air%cmolar ! MJT
          ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmah = (rt0+rough%rt1)*((1.+air%epsi)*rrsw + rrbw) &
               + air%epsi * (rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbh = (-air%rlam/capp)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmch = ((1.+air%epsi)*rrsw + rrbw)*rt0*rough%rt1* &
               (canopy%fhv + canopy%fhs)/(air%rho*capp)
          ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbe = (rt0+ssoil%wetfac*rough%rt1)* &
                 ((1.+air%epsi)*rrsw + rrbw)+(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmce = ((1.+air%epsi)*rrsw + rrbw)*rt0*rough%rt1* &
                 (canopy%fev + canopy%fes)/(air%rho*air%rlam)
          ! Within canopy air temperature:
          met%tvair = met%tk + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)

! tvair clobbered at present
          met%tvair = max(met%tvair , min( ssoil%tss, met%tk) - 5.0)
          met%tvair = min(met%tvair , max( ssoil%tss, met%tk) + 5.0)

          ! recalculate using canopy within temperature
          !     where (veg%meth .eq. 0 )
          met%qvair = met%qv + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          met%qvair = max(0.0,met%qvair)
! qvair clobbered at present
          met%qvair =  max(met%qvair ,min( ssoil%qstss, met%qv))
          met%qvair =  min(met%qvair ,max(ssoil%qstss, met%qv))

          ! Saturated specific humidity in canopy:
          qstvair = qsatf((met%tvair-tfrz),met%pmb)
          ! Saturated vapour pressure deficit in canopy:
          met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.


       END WHERE ! (veg%meth > 0 .and. veg%vlaiw > .01 .and. rough%hruff > rough%z0soilsn)

!       CALL define_air (met, air)

       ! Set radiative temperature as within canopy air temp:
!       met%tvrad = met%tvair      ! corrected 9/10/09

       ! Ground heat flux:
       !canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
       !

       dq = ssoil%qstss - met%qvair
       WHERE (ssoil%snowd > 1.0)
          dq = max( -0.1e-3, dq)
       END WHERE
       ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
       ! Net radiation absorbed by soil: 
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
                    sboltz*canopy%tv**4 - emsoil*sboltz* tss4

!       sss=air%dsatdk
!       cc1=sss/(sss+air%psyc )
!       cc2=air%psyc /(sss+air%psyc )
!       ssoil%potev = cc1 * (canopy%fns - canopy%ga) +  &
!          cc2 * air%rho * air%rlam*(qsatf((met%tvair-tfrz),met%pmb) - &
!          met%qvair)/ssoil%rtsoil
       ! Soil latent heat:
          canopy%fess= ssoil%wetfac * ssoil%potev
       ! canopy%fess = canopy%fess * (1-ssoil%pudsto/ssoil%pudsmx)
          canopy%fess = canopy%fess * (1-pwet)

          WHERE (ssoil%snowd < 0.1 .and. canopy%fess .gt. 0. )
            !allow the top soil moisture below wilting i.e. soil%swilt/3.
             !canopy%fess= min(canopy%fess,max(0._r_2,ssoil%wb(:,1)-soil%swilt/3.0)*soil%zse(1) &
             !  * 1000. * air%rlam / dels)
!             canopy%fess= min(canopy%fess,max(0._r_2,ssoil%wb(:,1)-soil%swilt)*soil%zse(1) &  
!               * 1000. * air%rlam / dels)                                                 
             canopy%fess= min(canopy%fess,max(0._r_2,((ssoil%wb(:,1)-soil%swilt)*soil%zse(1) &
                                  *1000.) - ssoil%evapfbl(:,1) )*air%rlam / dels)
         
             canopy%fess = min(canopy%fess,(ssoil%wb(:,1)-ssoil%wbice(:,1)) * soil%zse(1) &
                  * 1000. * air%rlam / dels)
          END WHERE
          ssoil%cls=1.
          WHERE (ssoil%snowd >= 0.1)
             ssoil%cls = 1.1335
             canopy%fess= min(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
          END WHERE

          ! Soil sensible heat:
          canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
          !EAK, 09/10 - commented out sets
          !where( canopy%fhs <= -80.0 .and. ktau <= 12 ) canopy%fhs = -80.0
          !where( canopy%fhs >  300.0 .and. ktau <= 12 ) canopy%fhs = 300.0
          !where( canopy%fhv <= -80.0 .and. ktau <= 12 ) canopy%fhv = -80.0
          !where( canopy%fhv >  300.0 .and. ktau <= 12 ) canopy%fhv = 300.0
          !where( canopy%fes <= -10.0 .and. ktau <= 12 ) canopy%fes = -10.0
          !where( canopy%fes >  300.0 .and. ktau <= 12 ) canopy%fes = 300.0
          !where( canopy%fev <= -10.0 .and. ktau <= 12 ) canopy%fev = -10.0
          !where( canopy%fev >  300.0 .and. ktau <= 12 ) canopy%fev = 300.0
!          where( canopy%fhs <= -80.0 .and. ktau <= 20 ) canopy%fhs = -80.0
         ! Evaporation form soil puddle already done
!1          canopy%fesp = min(ssoil%pudsto/dels*air%rlam,max(ssoil%pudsto/ssoil%pudsmx*ssoil%potev,0.))
          canopy%fesp = min(ssoil%pudsto/dels*air%rlam,max(pwet*ssoil%potev,0.))
          canopy%fes = canopy%fess + canopy%fesp
                                         
          ! Set total latent heat:
          canopy%fe = canopy%fev + canopy%fes

          ! Set total sensible heat:
          canopy%fh = canopy%fhv + canopy%fhs

          canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls

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
          canopy%epot = ((1.-rad%transd)*canopy%fevw_pot + rad%transd*ssoil%potev) * dels/air%rlam  
                      ! convert to mm/day
          canopy%wetfac_cs = min(1.0,canopy%fe / (canopy%fevw_pot + ssoil%potev))
          WHERE ( canopy%wetfac_cs .le. 0. )  &
           canopy%wetfac_cs = max(0.,min(1., &
                      max(canopy%fev/canopy%fevw_pot,canopy%fes/ssoil%potev)))

       ! monin-obukhov stability parameter zetar=zref/l
       !        recompute zetar for the next iteration, except on last iteration
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
          !     constrain zeta to zetpos and zetneg (set in param0)
          canopy%zetar(:,iterplus) = min(zetpos,canopy%zetar(:,iterplus))        ! zetar too +
          canopy%zetar(:,iterplus) = max(zetneg,canopy%zetar(:,iterplus))        ! zetar too -
       END IF ! (iter < niter) 
    END DO           ! do iter = 1, niter

    canopy%gswx_T = rad%fvlai(:,1)/max(0.01,veg%vlaiw(:))*canopy%gswx(:,1) &
           + rad%fvlai(:,2)/max(0.01,veg%vlaiw(:))*canopy%gswx(:,2)
    canopy%gswx_T = max(1.e-05,canopy%gswx_T )
    canopy%gs_vs = canopy%gswx_T + rad%transd * (0.01*ssoil%wb(:,1)/soil%sfc)**2
    where ( soil%isoilm == 9 ) canopy%gs_vs = 1.e6


    canopy%cduv = canopy%us * canopy%us / (max(met%ua,umin))**2
    canopy%cdtq = canopy%cduv *(LOG(rough%zref_uv / rough%z0m) -          &
      psim( canopy%zetar(:,niter) * rough%zref_uv/rough%zref_tq )) /      &
     (LOG( rough%zref_uv /(0.1*rough%z0m) ) - psis(canopy%zetar(:,niter)) )

    ! Calculate screen temperature:
    ! 1) original method from SCAM

   ! screen temp., windspeed and relative humidity at 1.5m
   ! screen temp., windspeed and relative humidity at 2.0m
    tstar = - canopy%fh / ( air%rho*capp*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
    zscrn = max(rough%z0m,2.0-rough%disp)
    denom = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) + &
         psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk

!%% change by Ashok Luhar
!    where (canopy%zetar(:,iterplus) > 0.7  .and. ssoil%snowd < 0.01)
!!            zeta2(:)=zetar(:,iterplus) * zscrn / rough%zref
!            zeta2=canopy%zetar(:,iterplus) * zscrn / rough%zref_tq
!            denom =alpha1* ((canopy%zetar(:,iterplus)**beta1* &
!               (1.0+gamma1*canopy%zetar(:,iterplus)**(1.0-beta1)))  &       
!             - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk 
!     endwhere
!%%
    ! Calculate screen temperature:
    canopy%tscrn = met%tc - tstar * denom

    ! Calculate radiative/skin temperature; at this stage old soil temperature is used
    rad%trad = ( (1.-rad%transd)*canopy%tv**4 + rad%transd * ssoil%tss**4 )**0.25

! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman

     rghlai = veg%vlaiw
!     where(ssoil%snowd.lt.0.001) rghlai = min( 3.,veg%vlaiw)
     where(ssoil%snowd.lt.0.001.and.veg%iveg .ne. 1) rghlai = min(3.,veg%vlaiw)
     term1=0.
     term2=0.
     term5=0.
     term3 = 0. ! Work around for Intel compiler problem with nested wheres
     zscl = 0.
     where ( veg%vlaiw > 0.01 .and. rough%hruff > 0.01)
        zscl = max(rough%z0soilsn,2.0)
        where ( rough%hruff  > 0.0 .and. rough%disp  > 0.0 )
           term1 = EXP(2*csw*rghlai*(1-zscl/rough%hruff))
           term2 = EXP(2*csw*rghlai*(1-rough%disp/rough%hruff))
           term5 = MAX(2./3.*rough%hruff/rough%disp, 1.)
        endwhere
        term3 = a33**2*ctl*2*csw*rghlai
        where( zscl < rough%disp )
            r_sc = term5 * LOG(zscl/rough%z0soilsn) * ( exp(2*csw*rghlai) - term1 ) / term3
        elsewhere ( rough%disp <= zscl .and. zscl < rough%hruff )
            r_sc = rough%rt0us + term5 * ( term2 - term1 ) / term3
        elsewhere ( rough%hruff <= zscl .and. zscl <  rough%zruffs )
            r_sc = rough%rt0us + rough%rt1usa + term5 * ( zscl - rough%hruff ) /  &
                                                          (a33**2*ctl*rough%hruff)
        elsewhere (zscl >= rough%zruffs )
            r_sc = rough%rt0us + rough%rt1usa + rough%rt1usb +  &
              ( log( (zscl - rough%disp)/MAX(rough%zruffs-rough%disp, rough%z0soilsn) ) &
                - psis( (zscl - rough%disp) / (rough%zref_tq/canopy%zetar(:,iterplus)) )  &
                + psis( (rough%zruffs - rough%disp) / (rough%zref_tq/canopy%zetar(:,iterplus)) )  &
              ) / vonk
        endwhere

        canopy%tscrn = ssoil%tss + (met%tk - ssoil%tss) * min(1.,r_sc /                            &
            max(1.,rough%rt0us + rough%rt1usa + rough%rt1usb + rt1usc)) - tfrz   ! in deg C
     endwhere         
     rsts = qsatf(canopy%tscrn, met%pmb)
     qtgnet = rsts * ssoil%wetfac - met%qv
     WHERE (qtgnet .gt. 0. )
        qsurf = rsts * ssoil%wetfac
     ELSEWHERE
        qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
     END WHERE
     canopy%qscrn = met%qv - qstar * denom
     where ( veg%vlaiw > 0.01 .and. rough%hruff > 0.01)
        canopy%qscrn =  qsurf + (met%qv - qsurf) * min(1.,r_sc /                            &
                            max(1.,rough%rt0us + rough%rt1usa + rough%rt1usb + rt1usc))
     endwhere

    ! MJT suggestion
!   calculation of ua_10 at 10m
    where ( rough%hruff > 0.01 ) &
       beta_sc = vonk /                                                                     &
           ( LOG( (rough%hruff - rough%disp) / rough%z0m )                               &
                 - psim( (rough%hruff - rough%disp) / (rough%zref_uv/canopy%zetar(:,iterplus)) )  &
                 + psim(  rough%z0m / (rough%zref_uv/canopy%zetar(:,iterplus)) )                  &
            ) 
    zscl_scrn = max(rough%z0soilsn,10.)
    where ( rough%hruff <= 1.e-2 )
       canopy%ua_10m = canopy%us/vonk *                                           &
           ( LOG( zscl_scrn  / rough%z0m )                                        &
                 - psim( zscl_scrn * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim( rough%z0m * canopy%zetar(:,iterplus)/ rough%zref_uv )    &
            )
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff > zscl_scrn )
       canopy%ua_10m = canopy%us *                                                       &
                       exp( (zscl_scrn-rough%hruff)*beta_sc/(rough%disp*vonk) ) / beta_sc
    elsewhere
       canopy%ua_10m = canopy%us/vonk *                                                          &
           ( LOG( (zscl_scrn - rough%disp) / rough%z0m )                                         &
                 - psim( (zscl_scrn - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim(  rough%z0m * canopy%zetar(:,iterplus) / rough%zref_uv )                 &
            )
    endwhere    
!   calculation of uscrn at 1.8m
    zscl_scrn = max(rough%z0soilsn,1.8)
    where ( rough%hruff <= 1.e-2 )
       canopy%uscrn = canopy%us/vonk *                                            &
           ( LOG( zscl_scrn  / rough%z0m )                                        &
                 - psim( zscl_scrn * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim( rough%z0m * canopy%zetar(:,iterplus)/ rough%zref_uv )    &
            )
    elsewhere ( rough%hruff > 1.e-2 .and.  rough%hruff > zscl_scrn )
       canopy%uscrn = canopy%us *                                                       &
                       exp( (zscl_scrn-rough%hruff)*beta_sc/(rough%disp*vonk) ) / beta_sc
    elsewhere
       canopy%uscrn = canopy%us/vonk *                                                           &
           ( LOG( (zscl_scrn - rough%disp) / rough%z0m )                                         &
                 - psim( (zscl_scrn - rough%disp) * canopy%zetar(:,iterplus) / rough%zref_uv )   &
                 + psim(  rough%z0m * canopy%zetar(:,iterplus) / rough%zref_uv )                 &
            )
    endwhere    


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
         * poolcoef1 /(365.0*24.0*3600.0)                ! 24/05
    canopy%frpw = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
         * poolcoef1w /(365.0*24.0*3600.0)                ! 24/05
    canopy%frpr = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
         * poolcoef1r /(365.0*24.0*3600.0)                ! 24/05

    ! This section to be updated as part of carbon module upgrade;
    ! frs is currently calculated in carbon module.

    !---?use? diagnostics for soil respiration
    !---def. DIAG_SOIL_RESP cable_directives
    if(DIAG_SOIL_RESP == 'off') then
       canopy%frs = rsoil(soil%rs20,avgwrs,avgtrs)
       canopy%frs = canopy%frs &
          * sum(spread(bgc%ratecs,1,mp) * bgc%csoil,2) &
          /(365.0*24.0*3600.0)   !convert 1/year to 1/second
    endif

    WHERE (ssoil%snowd > 1.)
       canopy%frs       = canopy%frs / max(0.001,min(100.,ssoil%snowd))
    END WHERE
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp

    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
         dels * 1.0e3 / (rhow*air%rlam)
    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm
    ! Modify canopy water storage for evaporation:
    canopy%cansto = max(canopy%cansto-max(0.0,canopy%fevw)*dels*1.0e3/ &
         (rhow*air%rlam), 0.0)
    ! Calculate canopy water storage excess:
    !canopy%spill=max(0.,min(0.2*canopy%cansto,max(0.0, canopy%cansto-cansat)))
    canopy%spill=max(0.0, canopy%cansto-cansat)
    ! Move excess canopy water to throughfall:
    canopy%through = canopy%through + canopy%spill
    ! Initialise 'throughfall to soil' as 'throughfall from canopy'; snow may absorb
    canopy%precis = max(0.,canopy%through)
    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill
    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-canopy%oldcansto
    ! calculate dgdtg, derivative of ghflux
    ssoil%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssoil%tss  ! d(canopy%fns)/d(ssoil%tgg)
    ssoil%dfh_dtg = air%rho*capp/ssoil%rtsoil      ! d(canopy%fhs)/d(ssoil%tgg)
    ssoil%dfe_ddq = ssoil%wetfac*air%rho*air%rlam/ssoil%rtsoil  ! d(canopy%fes)/d(dq)
    ssoil%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
         /((tetenc+ssoil%tss-tfrz)**2)*exp(tetenb*(ssoil%tss-tfrz)/(tetenc+ssoil%tss-tfrz))
    canopy%dgdtg = ssoil%dfn_dtg - ssoil%dfh_dtg - ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg

    bal%drybal=REAL(ecy+hcy,r_1)-SUM(rad%rniso,2) &
         +capp*rmair*(tlfy-met%tvair)*SUM(rad%gradis,2)
    !ypw: energy balance of the wet canopy
    bal%wetbal=canopy%fevw+canopy%fhvw-SUM(rad%rniso,2)*veg%fwet &
         +canopy%fhvw*SUM(rad%gradis,2)/MAX(0.001,ghwet)

    DEALLOCATE(cansat,ghwet,gbhu)
    DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
    DEALLOCATE(ecy, hcy, rny)
    DEALLOCATE(gbhf, gswmin, csx)

  CONTAINS
    !--------------------------------------------------------------------------
    ELEMENTAL FUNCTION qsatf(tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1)           :: r    ! result; sat sp humidity
      r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END FUNCTION qsatf
    !---------------------------------------------------------
    FUNCTION ej3x(parx,x) result(z)
      REAL(r_1), INTENT(IN)     :: parx
      REAL(r_1), INTENT(IN)     :: x
      REAL(r_1)                 :: z
      z = max(0.0, &
           0.25*((alpha3*parx+x-sqrt((alpha3*parx+x)**2 - &
           4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
    END FUNCTION ej3x
    !---------------------------------------------------------
    FUNCTION ej4x(parx,x) result(z)
      REAL(r_1), INTENT(IN)     :: parx
      REAL(r_1), INTENT(IN)     :: x
      REAL(r_1)                 :: z
      z = max(0.0, &
           (alpha4*parx+x-sqrt((alpha4*parx+x)**2 - &
           4.0*convx4*alpha4*parx*x))/(2.0*convx4))
    END FUNCTION ej4x
    !---------------------------------------------------------
    ! Explicit array dimensions as temporary work around for NEC inlining problem
    FUNCTION xvcmxt4(x) result(z)
      REAL(r_1), PARAMETER      :: q10c4 = 2.0

      ! modifying input to single real - kdcorbin, 09/10
      !REAL(r_1), DIMENSION(mp), INTENT(IN)   :: x
      !REAL(r_1), DIMENSION(mp)                  :: z
      REAL(r_1), INTENT(IN) :: x
      REAL(r_1) :: z

      z = q10c4 ** (0.1 * x - 2.5) / &
           ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
    END FUNCTION xvcmxt4
    !---------------------------------------------------------
    FUNCTION xvcmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for vcmax for c3 plants

      ! modifying input to single real - kdcorbin, 09/10
      !REAL(r_1), DIMENSION(mp), INTENT(IN)   :: x
      !REAL(r_1), DIMENSION(mp)               :: xvcnum
      !REAL(r_1), DIMENSION(mp)               :: xvcden
      !REAL(r_1), DIMENSION(mp)               :: z

      REAL(r_1), INTENT(IN) :: x
      REAL(r_1) :: xvcnum,xvcden,z

      REAL(r_1), PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
      REAL(r_1), PARAMETER  :: xVccoef = 1.17461 ! derived parameter
                        ! xVccoef=1.0+exp((EntropJx*TrefK-EHdJx)/(Rconst*TrefK))

      xvcnum=xvccoef*exp((ehavc/(rgas*trefk))*(1.-trefk/x))
      xvcden=1.0+exp((entropvc*x-ehdvc)/(rgas*x))
      z = max(0.0,xvcnum/xvcden)
    END FUNCTION xvcmxt3
    !---------------------------------------------------------
    FUNCTION xejmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants

      ! modifying input to single real - kdcorbin, 09/10
      !REAL(r_1), DIMENSION(mp), INTENT(IN)   :: x
      !REAL(r_1), DIMENSION(mp)               :: xjxnum
      !REAL(r_1), DIMENSION(mp)               :: xjxden
      !REAL(r_1), DIMENSION(mp)               :: z

      REAL(r_1), INTENT(IN) :: x
      REAL(r_1) :: xjxnum,xjxden,z   

      REAL(r_1), PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
      REAL(r_1), PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
      REAL(r_1), PARAMETER  :: xjxcoef = 1.16715 ! derived parameter

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
      REAL(r_1), INTENT(IN)     :: zeta
      REAL(r_1)                 :: r
      REAL(r_1)                 :: x
      REAL(r_1), PARAMETER      :: gu = 16.0
      REAL(r_1), PARAMETER      :: gs = 5.0

      REAL(r_1)                 :: z 
      REAL(r_1)                 :: stable 
      REAL(r_1)                 :: unstable 

      z = 0.5 + sign(0.5,zeta) !z=1 in stable, 0 in unstable
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
      REAL(r_1), INTENT(IN)     :: zeta
      REAL(r_1)                 :: r
      REAL(r_1), PARAMETER      :: gu = 16.0
      REAL(r_1), PARAMETER      :: gs = 5.0

      REAL(r_1)                 :: z
      REAL(r_1)                 :: y
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable

      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
      stable = -gs*zeta
      y      = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
      REAL(r_1), INTENT(IN)     :: rpconst
      REAL(r_1), INTENT(IN)     :: rpcoef
      REAL(r_1), INTENT(IN)     :: tair
      REAL(r_1)                 :: z
      z = rpconst * exp(rpcoef * tair)
    END FUNCTION rplant
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rsoil(rsconst, avgwrs, avgtrs) result(z)
      REAL(r_1), INTENT(IN)     :: rsconst
      REAL(r_1), INTENT(IN)     :: avgwrs
      REAL(r_1), INTENT(IN)     :: avgtrs
      REAL(r_1)                 :: z
      z = rsconst * min(1.0, max(0.0, min(&
           -0.0178+0.2883*avgwrs+5.0176*avgwrs*avgwrs-4.5128*avgwrs*avgwrs*avgwrs, &
           0.3320+22.6726*exp(-5.8184*avgwrs)))) &
           * min(1.0, max(0.0, min( 0.0104*(avgtrs**1.3053), 5.5956-0.1189*avgtrs)))
    END FUNCTION rsoil
    !---------------------------------------------------------
  SUBROUTINE dryLeaf(dels,phenps,rad,rough,air,met,veg,canopy)
    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy

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
!    REAL(r_1), DIMENSION(mp,ms)  :: evapfbl !
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
    REAL(r_1), DIMENSION(mp,ms)  :: oldevapfbl

    INTEGER(i_d)   :: i,j !loop counts - kdcorbin, 09/10
    REAL, PARAMETER :: min_lai=0.001 !kdcorbin, 09/10

    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    rdx = 0.0

    !kdcorbin, 08/10
    csx = SPREAD(met%ca,2,mf)
    anx = 0.0

    rnx = SUM(rad%rniso,2)
    abs_deltlf = 999.0
    !canopy%fevc = 0.0
    !ssoil%evapfbl = 0.0
    oldevapfbl = 0.0

    DO kk=1,mp
      IF(veg%vlaiw(kk) <= min_lai) THEN
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

    !kdcorbin, 08/10 - doing all points all the time
    DO WHILE (k < maxiter)
       k = k + 1

       !kdcorbin, 09/10 - put calculations in do/if statements
       DO i=1,mp
          IF (veg%vlaiw(i) > min_lai .AND. abs_deltlf(i) > 0.1) Then

           ! Grashof number (Leuning et al, 1995) eq E4:
           gras(i) = MAX(1.0e-6, &
               1.595E8*ABS(tlfx(i)-met%tvair(i))*(veg%dleaf(i)**3.0))
           ! See Appendix E in (Leuning et al, 1995):
           gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*dheat &
                *(gras(i)**0.25) / veg%dleaf(i)
           gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5*dheat &
                *(gras(i)**0.25) / veg%dleaf(i)
           gbhf(i,:) = max(1.e-6,gbhf(i,:))

           ! Conductance for heat:
           gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))

           ! Conductance for heat and longwave radiation:
           ghr(i,:) = rad%gradis(i,:)+gh(i,:)

           ! Leuning 2002 (P C & E) equation for temperature response
           ! used for Vcmax for C3 plants:
           temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
           vcmxt3(i,1) = rad%scalex(i,1) * temp(i)
           vcmxt3(i,2) = rad%scalex(i,2) * temp(i)

           ! Temperature response of Vcmax for C4 plants (Collatz et al 1989):
           temp(i) = xvcmxt4(tlfx(i)-tfrz) * veg%vcmax(i) * veg%frac4(i)
           vcmxt4(i,1) = rad%scalex(i,1) * temp(i)
           vcmxt4(i,2) = rad%scalex(i,2) * temp(i)

           ! Leuning 2002 (P C & E) equation for temperature response
           ! used for Jmax for C3 plants:
           temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
           ejmxt3(i,1) = rad%scalex(i,1) * temp(i)
           ejmxt3(i,2) = rad%scalex(i,2) * temp(i)

           ! Difference between leaf temperature and reference temperature:
           tdiff(i) = tlfx(i) - trefk
           ! Michaelis menten constant of Rubisco for CO2:
           conkct(i) = conkc0 * EXP((ekc/(rgas*trefk)) * (1.0-trefk/tlfx(i)))
           ! Michaelis menten constant of Rubisco for oxygen:
           conkot(i) = conko0 * EXP((eko/(rgas*trefk)) * (1.0-trefk/tlfx(i)))

           ! Store leaf temperature
           tlfxx(i) = tlfx(i)

           ! "d_{3}" in Wang and Leuning, 1998, appendix E:
           cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
           cx2(i) = 2.0 * gam0 * (1.0 + gam1*tdiff(i) + gam2*tdiff(i)*tdiff(i))

           ! All equations below in appendix E in Wang and Leuning 1998 are
           ! for calculating anx, csx and gswx for Rubisco limited,
           ! RuBP limited, sink limited
           temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
           temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
           vx3(i,1)  = ej3x(temp2(i,1),ejmxt3(i,1))
           vx3(i,2)  = ej3x(temp2(i,2),ejmxt3(i,2))

           temp2(i,1) = rad%qcan(i,1,1) * jtomol * veg%frac4(i)
           temp2(i,2) = rad%qcan(i,2,1) * jtomol * veg%frac4(i)
           vx4(i,1)  = ej4x(temp2(i,1),vcmxt4(i,1))
           vx4(i,2)  = ej4x(temp2(i,2),vcmxt4(i,2))

           rdx(i,1) = (cfrd3*vcmxt3(i,1) + cfrd4*vcmxt4(i,1))*fwsoil(i)
           rdx(i,2) = (cfrd3*vcmxt3(i,2) + cfrd4*vcmxt4(i,2))*fwsoil(i)
           xleuning(i,1) = (fwsoil(i) / (csx(i,1)-co2cp3))  &
                     * ((1.0-veg%frac4(i)) * a1c3 / (1.0+dsx(i)/d0c3) &
                         + veg%frac4(i)    * a1c4 / (1.0+dsx(i)/d0c4))
           xleuning(i,2) = (fwsoil(i) / (csx(i,2)-co2cp3))  &
                     * ((1.0-veg%frac4(i)) * a1c3 / (1.0+dsx(i)/d0c3) &
                         + veg%frac4(i)    * a1c4 / (1.0+dsx(i)/d0c4))
       ENDIF
      ENDDO !i=1,mp

      !kdcorbin, 09/10 - put in photosynthesis subroutine to replace function
      call photosynthesis(csx(:,:),SPREAD(cx1(:),2,mf), &
                SPREAD(cx2(:),2,mf),gswmin(:,:),rdx(:,:), &
                vcmxt3(:,:),vcmxt4(:,:),vx3(:,:),vx4(:,:),xleuning(:,:), &
                rad%fvlai(:,:),SPREAD(abs_deltlf,2,mf),anx(:,:))

    !kdcorbin, 09/10 - put calculations in do/if statements
    DO i=1,mp
         IF (veg%vlaiw(i) > min_lai .AND. abs_deltlf(i) > 0.1) Then
           !kdcorbin, 10/10 - added second do and test
            DO kk=1,mf
               !if (rad%fvlai(i,kk) > min_lai) Then
                 csx(i,kk) = met%ca(i) - rgbwc*anx(i,kk) / &
                           (gbhu(i,kk) + gbhf(i,kk))
                 csx(i,kk) = met%ca(i) - rgbwc*anx(i,kk) / &
                           (gbhu(i,kk) + gbhf(i,kk))
                 csx(i,kk) = MAX(1.0e-4,csx(i,kk))
                 canopy%gswx(i,kk) = MAX(1.e-3, gswmin(i,kk) + &
                           MAX(0.0,rgswc*xleuning(i,kk)*anx(i,kk)))
 
                 !Recalculate conductance for water:
                 gw(i,kk) = 1.0/(1.0/canopy%gswx(i,kk) + &
                          1.0/(1.075*(gbhu(i,kk)+gbhf(i,kk))))
                 gw(i,kk) = MAX(gw(i,kk),0.00001)
           
                 !Modified psychrometric constant (Monteith and Unsworth, 1990)
                 psycst(i,kk) = air%psyc(i)*REAL(ghr(i,kk)/gw(i,kk),r_1)
               !endif
             ENDDO  !kk=1,mf

             ecx(i) = (air%dsatdk(i)*(rad%rniso(i,1) &
                 - capp*rmair*(met%tvair(i)-met%tk(i)) &
                 * rad%gradis(i,1)) + capp*rmair*met%dva(i)*ghr(i,1)) &
                 / (air%dsatdk(i)+psycst(i,1)) &
                 + (air%dsatdk(i)*(rad%rniso(i,2) &
                 -  capp*rmair*(met%tvair(i)-met%tk(i))*rad%gradis(i,2)) &
                 + capp*rmair*met%dva(i)*ghr(i,2)) &
                 / (air%dsatdk(i)+psycst(i,2))

      
             IF (ecx(i) > 0.0 .AND. veg%fwet(i) < 1.0) Then
                evapfb(i) = (1.0-veg%fwet(i))*REAL(ecx(i),r_1)*dels/air%rlam(i)
                DO kk = 1,ms
                   ssoil%evapfbl(i,kk) = MIN(evapfb(i)*soil%froot(i,kk), &
                      MAX(0.0,REAL(ssoil%wb(i,kk),r_1) - 1.1*soil%swilt(i)) &
                           * soil%zse(kk)*1000.0)
!                   soil%vapfbl(i,kk) = MIN(evapfb(i)*soil%froot(i,kk), &
!                      MAX(0.0,MIN(REAL(ssoil%wb(i,kk),r_1) - 1.1*soil%swilt(i), &
!                            REAL(ssoil%wb(i,kk) - 1.05*ssoil%wbice(i,kk),r_1))) &
!                           * soil%zse(kk)*1000.0)
                 ENDDO

                 canopy%fevc(i) = SUM(ssoil%evapfbl(i,:))*air%rlam(i)/dels
          
                 ecx(i) = canopy%fevc(i) / (1.0-veg%fwet(i))

             ENDIF

             ! Update canopy sensible heat flux:
             hcx(i) = (SUM(rad%rniso(i,:))-ecx(i) &
                - capp*rmair*(met%tvair(i)-met%tk(i))  &
                * SUM(rad%gradis(i,:)))    &
                * SUM(gh(i,:))/ SUM(ghr(i,:))

             ! Update leaf temperature:
             tlfx(i)=met%tvair(i)+REAL(hcx(i),r_1)/(capp*rmair*SUM(gh(i,:)))
       
             ! Update net radiation for canopy:
             rnx(i) = SUM(rad%rniso(i,:)) - &
                      capp*rmair*(tlfx(i)-met%tk(i))*  &
                      SUM(rad%gradis(i,:))

             ! Update leaf surface vapour pressure deficit:
             dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

             ! Store change in leaf temperature between successive iterations:
             deltlf(i) = tlfxx(i)-tlfx(i)
             abs_deltlf(i) = ABS(deltlf(i))

           ENDIF !lai/abs_deltlf
        ENDDO !i=1,mp

       ! Where leaf temp change b/w iterations is significant, and
       ! difference is smaller than the previous iteration, store results:
!       WHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
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
          !oldevapfbl(:,:) = ssoil%evapfbl(:,:)
          oldevapfbl(:,1) = ssoil%evapfbl(:,1)
          oldevapfbl(:,2) = ssoil%evapfbl(:,2)
          oldevapfbl(:,3) = ssoil%evapfbl(:,3)
          oldevapfbl(:,4) = ssoil%evapfbl(:,4)
          oldevapfbl(:,5) = ssoil%evapfbl(:,5)
          oldevapfbl(:,6) = ssoil%evapfbl(:,6)
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
          oldevapfbl(:,:) = ssoil%evapfbl(:,:)
       END IF
    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < maxiter)
    ! dry canopy flux
    canopy%fevc = (1.0-veg%fwet) * ecy

    canopy%frday = 12.0 * SUM(rdy, 2)
    canopy%fpn = -12.0 * SUM(an_y, 2)
    ssoil%evapfbl(:,:) = oldevapfbl(:,:)


  END SUBROUTINE dryLeaf
    !---------------------------------------------------------
  SUBROUTINE photosynthesis(csxz,cx1z,cx2z,gswminz,rdxz,vcmxt3z,vcmxt4z, &
                    vx3z,vx4z,xleuningz,vlaiz,deltlfz,anxz)
    ! inputs:
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: cx1z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: cx2z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: gswminz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: rdxz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vcmxt3z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vcmxt4z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vx4z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vx3z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: xleuningz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vlaiz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: deltlfz
    REAL(r_1), DIMENSION(mp,mf), INTENT(INOUT) :: anxz
    !local variables
    REAL(r_2), DIMENSION(mp,mf) :: coef0z,coef1z,coef2z
    REAL(r_2), DIMENSION(mp,mf) :: ciz,delcxz
    REAL(r_2), DIMENSION(mp,mf) :: anrubiscoz,anrubpz,ansinkz
    REAL(r_1), PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
                                             ! Bonan,LSM version 1.0, p106)
    INTEGER :: i,j   !kdcorbin, 09/10
    REAL, PARAMETER :: min_lai=0.001  !kdcorbin, 09/10

    ! rgswc - inherited from canopy_module's USE photosynthetic_constants
    ! mp - inherited from canopy_module
    ! mf - inherited from canopy_module
   
   DO i=1,mp
      !if (SUM(vlaiz(1,:)) > min_lai) Then  !kdcorbin, 10/10
         DO j=1,mf
            IF (vlaiz(i,j) .gt. min_lai .AND. deltlfz(i,j) .gt. 0.1) Then

    ! Rubisco limited:
     coef2z(i,j) = gswminz(i,j)/rgswc+xleuningz(i,j) * &
                   (vcmxt3z(i,j)-(rdxz(i,j)-vcmxt4z(i,j)))
     coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) * &
                   (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j)) &
                   + (gswminz(i,j)/rgswc)*(cx1z(i,j)-csxz(i,j)) &
                   - xleuningz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0 &
                   + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j)))
     coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) * &
                    (vcmxt3z(i,j)*cx2z(i,j)/2.0  &
                   + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j))) &
                   -(gswminz(i,j)/rgswc)*cx1z(i,j)*csxz(i,j)

     !kdcorbin,09/10 - new calculations
     IF (ABS(coef2z(i,j)) .gt. 1.0e-9 .AND. &
           ABS(coef1z(i,j)) .lt. 1.0e-9) Then
       ! no solution, give it a huge number
       ciz(i,j) = 99999.0 ! quadratic below cannot handle zero denominator
       anrubiscoz(i,j) = 99999.0    !should ciz=0 and anrubiscoz calculated?
     ENDIF

     ! solve linearly
     IF (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1e-9) Then
       ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)   ! same reason as above
       ciz(i,j)    = MAX(0.0_r_2,ciz(i,j))
       anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
                         (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
     ENDIF

     ! solve quadratic (only take the more positive solution)
     IF (ABS(coef2z(i,j)) >= 1.e-9) Then
       delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
       ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                      /(2.0*coef2z(i,j))
       ciz(i,j) = MAX(0.0_r_2,ciz(i,j))   ! must be positive, why?
       anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
              (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
     ENDIF

   ! RuBP limited:
     coef2z(i,j) = gswminz(i,j)/rgswc+xleuningz(i,j) &
                    *(vx3z(i,j)-(rdxz(i,j)-vx4z(i,j)))
     coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) * &
                   (vx3z(i,j)+vx4z(i,j)-rdxz(i,j))    &
                   + (gswminz(i,j)/rgswc)* &
                     (cx2z(i,j)-csxz(i,j))-xleuningz(i,j)  &
                   *(vx3z(i,j)*cx2z(i,j)/2.0 + cx2z(i,j)* &
                     (rdxz(i,j)-vx4z(i,j)))
     coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) * &
                     (vx3z(i,j)*cx2z(i,j)/2.0  &
                   + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j))) &
                   -(gswminz(i,j)/rgswc)*cx2z(i,j)*csxz(i,j)

     !kdcorbin, 09/10 - new calculations
     ! no solution, give it a huge number
     IF (ABS(coef2z(i,j)) < 1.0e-9 .AND. ABS(coef1z(i,j)) < 1.0e-9) Then
       ciz(i,j) = 99999.0
       anrubpz(i,j)  = 99999.0
     ENDIF
     ! solve linearly
     IF (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1.e-9) Then
        ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
        ciz(i,j) = MAX(0.0_r_2,ciz(i,j))
        anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
                (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
     ENDIF
     ! solve quadratic (only take the more positive solution)
     IF (ABS(coef2z(i,j)) >= 1.e-9) Then
         delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
         ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                      /(2.0*coef2z(i,j))
         ciz(i,j) = MAX(0.0_r_2,ciz(i,j)) 
         anrubpz(i,j)  = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
               (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
     ENDIF
       
   ! Sink limited:
     coef2z(i,j) = xleuningz(i,j)
     coef1z(i,j) = gswminz(i,j)/rgswc + xleuningz(i,j) &
                     * (rdxz(i,j) - 0.5*vcmxt3z(i,j))  &
                     + effc4 * vcmxt4z(i,j) - xleuningz(i,j) &
                     * csxz(i,j) * effc4 * vcmxt4z(i,j)
     coef0z(i,j) = -(gswminz(i,j)/rgswc)*csxz(i,j)*effc4*vcmxt4z(i,j) + &
                    (rdxz(i,j) -0.5*vcmxt3z(i,j))*gswminz(i,j)/rgswc

     ! no solution, give it a huge number
     IF (ABS(coef2z(i,j)) < 1.0e-9 .AND. ABS(coef1z(i,j)) < 1.0e-9) Then
       ciz(i,j) = 99999.0
       ansinkz(i,j)  = 99999.0
     ENDIF

     ! solve linearly
     IF (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1.e-9) Then
        ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
        ansinkz(i,j)  = ciz(i,j)
     ENDIF
     ! solve quadratic (only take the more positive solution)
     IF (ABS(coef2z(i,j)) >= 1.e-9) Then
        delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
        ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                     /(2.0*coef2z(i,j))
        ansinkz(i,j) = ciz(i,j)
     ENDIF
       
   ! minimal of three limited rates
     anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
     ENDIF
     ENDDO
    !ENDIF  !sum(vlaiz > min_lai)
   ENDDO
     
  END SUBROUTINE photosynthesis
    !---------------------------------------------------------

  SUBROUTINE wetLeaf(dels,rad,rough,air,met,veg,canopy)

    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy

   ! assuming the temperature of wet leaf is equal that of dry leaf ="tlfy"
    REAL(r_1), INTENT(IN)     :: dels ! integration time step (s)
    REAL(r_1), DIMENSION(mp)  :: ccfevw ! limitation term for
                                        ! wet canopy evaporation rate
    REAL(r_1), DIMENSION(mp)  :: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp)  :: ghrwet ! wet canopy cond: heat & thermal rad

    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    ghwet = 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    canopy%fevw_pot = 0.0
    WHERE (veg%vlaiw > 0.01)
    ! VEG SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
    ! calculate total thermal resistance, rthv in s/m
       ghwet = 2.0   * SUM((gbhu+gbhf),2)
       gwwet = 1.075 * SUM((gbhu+gbhf),2)
       ghrwet = SUM(rad%gradis,2) + ghwet
       ! Calculate fraction of canopy which is wet:
       veg%fwet = MAX(0.0,MIN(0.9,0.8*canopy%cansto/MAX(cansat,0.01)))
       ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
       ! to avoid excessive evaporation:
       ccfevw = MIN(canopy%cansto * air%rlam / dels, &
                    2.0 / (1440.0 / (dels/60.0)) * air%rlam)

       canopy%fevw = MIN(veg%fwet * (air%dsatdk*(SUM(rad%rniso,2)- &
                     capp*rmair*(met%tvair(:)-met%tk(:))*sum(rad%gradis,2)) &
                     + capp*rmair*met%dva*ghrwet) &
                         / (air%dsatdk+air%psyc*ghrwet/gwwet), ccfevw)
       canopy%fevw_pot = (air%dsatdk*(SUM(rad%rniso,2)- &
                     capp*rmair*(met%tvair(:)-met%tk(:))*sum(rad%gradis,2)) &
                     + capp*rmair*met%dva*ghrwet) &
                         / (air%dsatdk+air%psyc*ghrwet/gwwet)

       ! Calculate sens heat from wet canopy:
       !canopy%fhvw = (veg%fwet*(SUM(rad%rniso,2)-capp*rmair *              &
       !(met%tvair(:)-met%tk(:))*sum(rad%gradis,2))-canopy%fevw)*ghwet/ghrwet

       canopy%fhvw = veg%fwet*(SUM(rad%rniso,2)-capp*rmair*(tlfy-met%tk(:))* &
                     sum(rad%gradis,2)) - canopy%fevw
    END WHERE
  END SUBROUTINE wetLeaf

  END SUBROUTINE define_canopy

END MODULE canopy_module
