!!$ cable_checks.f90
!!$
!!$ Energy and mass conservation routines; acceptable ranges for i/o  
!!$ variables in offline netcdf driver.
!!$
!!$ Eva Kowalczyk, Gab Abramowitz 2006 CSIRO Marine and Atmospheric
!!$ Research/ University of New South Wales; gabsun@gmail.com
!!$
!!$ ranges_type below sets the acceptable ranges for all variables
!!$ coming in or going out of the offline netcdf driver. The mass_balance
!!$ and energy_balance subroutines calculate cumulative and per-timestep 
!!$ balances, as well as allow user to scrutinise balances in
!!$ particular sections of the code - largely for diagnostics/fault finding.
!!$ rh_sh - converts relative to sensible humidity if met file units require it

MODULE checks_module
  USE canopy_module
  USE define_types
  USE physical_constants
  IMPLICIT NONE
  TYPE units_type
     CHARACTER(LEN=1) :: Rainf ! 's' (mm/s) or 'h' (mm/h)
     CHARACTER(LEN=1) :: PSurf  ! 'h'(hPa or mbar) or 'P'(Pa)
     CHARACTER(LEN=1) :: Tair  ! 'C' or 'K'
     CHARACTER(LEN=1) :: Qair  ! '%' or 'g' (spec hum)
     CHARACTER(LEN=1) :: CO2air ! 'p' (ppmv)
     CHARACTER(LEN=1) :: Wind ! 'm'(m/s)
  END TYPE units_type
  TYPE(units_type) :: units
  TYPE ranges_type 
     REAL(r_1), DIMENSION(2) :: nav_lon = (/-360.0,360.0/)   
     REAL(r_1), DIMENSION(2) :: nav_lat = (/-90.0,90.0/)   
     REAL(r_1), DIMENSION(2) :: time     
     REAL(r_1), DIMENSION(2) :: timestp      
     ! possible forcing variables for CABLE
     REAL(r_1), DIMENSION(2) :: SWdown = (/0.0,1360.0/)  ! W/m^2
     REAL(r_1), DIMENSION(2) :: LWdown = (/0.0,750.0/)   ! W/m^2
     REAL(r_1), DIMENSION(2) :: Rainf = (/0.0,0.03/)     ! mm/s
     REAL(r_1), DIMENSION(2) :: Snowf = (/0.0,0.0085/)   ! mm/s
     REAL(r_1), DIMENSION(2) :: PSurf = (/500.0,1100.0/) ! mbar/hPa
     REAL(r_1), DIMENSION(2) :: Tair = (/200.0,333.0/)   ! K
     REAL(r_1), DIMENSION(2) :: Qair = (/0.0,0.04/)      ! g/g
     REAL(r_1), DIMENSION(2) :: CO2air = (/160.0,2000.0/)! ppmv   
     REAL(r_1), DIMENSION(2) :: Wind = (/0.0,75.0/)      ! m/s
     REAL(r_1), DIMENSION(2) :: Wind_N = (/-75.0,75.0/)  ! m/s
     REAL(r_1), DIMENSION(2) :: Wind_E = (/-75.0,75.0/)  ! m/s
     ! possible output variables
     REAL(r_1), DIMENSION(2) :: Qh = (/-1000.0,1000.0/)    ! W/m^2
     REAL(r_1), DIMENSION(2) :: Qle = (/-1000.0,1000.0/)   ! W/m^2
     REAL(r_1), DIMENSION(2) :: Qg = (/-1000.0,1000.0/)    ! W/m^2   
     REAL(r_1), DIMENSION(2) :: SWnet = (/0.0,1250.0/)     ! W/m^2
     REAL(r_1), DIMENSION(2) :: LWnet = (/-500.0,510.0/)   ! W/m^2 
     REAL(r_1), DIMENSION(2) :: Rnet = (/-500.0,1250.0/)   ! W/m^2 
     REAL(r_1), DIMENSION(2) :: Evap = (/-0.0003,0.00035/)      
     REAL(r_1), DIMENSION(2) :: Ewater = (/-0.0003,0.0003/)
     REAL(r_1), DIMENSION(2) :: ESoil = (/-0.0003,0.0003/)     
     REAL(r_1), DIMENSION(2) :: TVeg = (/-0.0003,0.0003/)    
     REAL(r_1), DIMENSION(2) :: ECanop = (/-0.0003,0.0003/)   
     REAL(r_1), DIMENSION(2) :: PotEvap = (/-0.0006,0.0006/)     
     REAL(r_1), DIMENSION(2) :: ACond = (/0.0,1.0/)    
     REAL(r_1), DIMENSION(2) :: SoilWet = (/-0.4,1.2/) 
     REAL(r_1), DIMENSION(2) :: Albedo = (/0.0,1.0/)    
     REAL(r_1), DIMENSION(2) :: VegT = (/213.0,333.0/)     
     REAL(r_1), DIMENSION(2) :: SoilTemp = (/213.0,343.0/)   
     REAL(r_1), DIMENSION(2) :: SoilMoist = (/0.0,2000.0/) 
     REAL(r_1), DIMENSION(2) :: Qs = (/0.0,5.0/)
     REAL(r_1), DIMENSION(2) :: Qsb = (/0.0,5.0/)
     REAL(r_1), DIMENSION(2) :: DelSoilMoist  = (/-2000.0,2000.0/) 
     REAL(r_1), DIMENSION(2) :: DelSWE  = (/-2000.0,2000.0/)       
     REAL(r_1), DIMENSION(2) :: DelIntercept = (/-100.0,100.0/)  
     REAL(r_1), DIMENSION(2) :: SnowT  = (/213.0,280.0/)        
     REAL(r_1), DIMENSION(2) :: BaresoilT = (/213.0,343.0/)     
     REAL(r_1), DIMENSION(2) :: AvgSurfT = (/213.0,333.0/)      
     REAL(r_1), DIMENSION(2) :: RadT = (/213.0,353.0/)         
     REAL(r_1), DIMENSION(2) :: SWE = (/0.0,2000.0/)           
     REAL(r_1), DIMENSION(2) :: RootMoist = (/0.0,2000.0/)     
     REAL(r_1), DIMENSION(2) :: CanopInt = (/0.0,100.0/)  
     REAL(r_1), DIMENSION(2) :: NEE = (/-70.0,50.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: NPP = (/-20.0,50.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: GPP = (/-10.0,50.0/) ! umol/m2/s 
     REAL(r_1), DIMENSION(2) :: AutoResp = (/-50.0,20.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: HeteroResp = (/-50.0,20.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: HSoil = (/-1000.0,1000.0/) 
     REAL(r_1), DIMENSION(2) :: HVeg = (/-1000.0,1000.0/)
     REAL(r_1), DIMENSION(2) :: SnowDepth = (/0.0,5.0/)
     ! parameters:
     REAL(r_1), DIMENSION(2) :: bch = (/2.0,15.0/)  
     REAL(r_1), DIMENSION(2) :: latitude = (/-90.0,90.0/)          
     REAL(r_1), DIMENSION(2) :: clay = (/0.0,1.0/)  
     REAL(r_1), DIMENSION(2) :: css = (/700.0,2000.0/)         
     REAL(r_1), DIMENSION(2) :: rhosoil = (/300.0,3000.0/)    
     REAL(r_1), DIMENSION(2) :: hyds = (/1.0E-6,2.0E-4/)
     REAL(r_1), DIMENSION(2) :: rs20 = (/0.1,10.0/)
     REAL(r_1), DIMENSION(2) :: sand = (/0.0,1.0/)      
     REAL(r_1), DIMENSION(2) :: sfc = (/0.1,0.5/)        
     REAL(r_1), DIMENSION(2) :: silt = (/0.0,1.0/)
     REAL(r_1), DIMENSION(2) :: ssat = (/0.35,0.5/)      
     REAL(r_1), DIMENSION(2) :: sucs = (/-0.8,-0.03/)       
     REAL(r_1), DIMENSION(2) :: swilt = (/0.05,0.4/)
     REAL(r_1), DIMENSION(2) :: froot = (/0.0,1.0/) 
     REAL(r_1), DIMENSION(2) :: zse = (/0.0,5.0/)  
     REAL(r_1), DIMENSION(2) :: canst1 = (/0.05,0.15/)     
     REAL(r_1), DIMENSION(2) :: dleaf = (/0.005,0.2/)      
     REAL(r_1), DIMENSION(2) :: ejmax = (/1.0E-5,3.0E-4/) 
     REAL(r_1), DIMENSION(2) :: frac4 = (/0.0,1.0/)
     REAL(r_1), DIMENSION(2) :: hc = (/0.0,100.0/)         
     REAL(r_1), DIMENSION(2) :: lai = (/0.0,8.0/)
     REAL(r_1), DIMENSION(2) :: rp20 = (/0.1,10.0/)       
     REAL(r_1), DIMENSION(2) :: rpcoef = (/0.05,1.5/)
     REAL(r_1), DIMENSION(2) :: shelrb = (/1.0,3.0/)     
     REAL(r_1), DIMENSION(2) :: vcmax = (/5.0E-6,1.5E-4/)      
     REAL(r_1), DIMENSION(2) :: xfang = (/-1.0,1.0/)      
     REAL(r_1), DIMENSION(2) :: ratecp = (/0.5,3.0/)    
     REAL(r_1), DIMENSION(2) :: ratecs = (/0.5,3.0/)   
     REAL(r_1), DIMENSION(2) :: refsbare = (/0.0,0.5/) 
     REAL(r_1), DIMENSION(2) :: taul = (/0.0,0.3/)     
     REAL(r_1), DIMENSION(2) :: refl = (/0.0,0.5/)    
     REAL(r_1), DIMENSION(2) :: tauw = (/0.0,0.1/)     
     REAL(r_1), DIMENSION(2) :: refw = (/0.0,0.5/)    
     REAL(r_1), DIMENSION(2) :: tminvj = (/-10.0,10.0/)   
     REAL(r_1), DIMENSION(2) :: tmaxvj = (/-5.0,15.0/)  
     REAL(r_1), DIMENSION(2) :: veg_class = (/1.0,20.0/)  
     REAL(r_1), DIMENSION(2) :: soil_class = (/1.0,20.0/)  
  END TYPE ranges_type
  TYPE(ranges_type),SAVE :: ranges
CONTAINS
  SUBROUTINE mass_balance(ktau,dels,ssoil,soil,canopy,met,air,bal)
    INTEGER(i_d), INTENT(IN)          :: ktau ! time step
    REAL(r_1),INTENT(IN)              :: dels ! time step size
    TYPE (soil_snow_type),INTENT(IN) :: ssoil ! soil data
    TYPE (soil_parameter_type),INTENT(IN) :: soil ! soil data
    TYPE (canopy_type),INTENT(IN):: canopy ! canopy variable data
    TYPE(met_type),INTENT(IN) :: met  ! met data
    TYPE (air_type),INTENT(IN) 	:: air
    REAL(r_2), DIMENSION(:,:,:),POINTER, SAVE :: bwb ! volumetric soil moisture
    REAL(r_2), DIMENSION(mp) :: delwb ! change in soilmoisture b/w tsteps
    REAL(r_1), DIMENSION(mp) :: canopy_wbal ! canopy water balance
    TYPE (balances_type),INTENT(INOUT):: bal 
    INTEGER :: j, k ! do loop counter
    
    IF(ktau==1) THEN
       ALLOCATE( bwb(mp,ms,2) )
       bwb(:,:,1)=ssoil%wb ! initial vlaue of soil moisture
    ELSE
       ! Calculate change in soil moisture b/w timesteps:
       IF(MOD(REAL(ktau),2.0)==1.0) THEN         ! if odd timestep
          bwb(:,:,1)=ssoil%wb
          DO k=1,mp           ! current smoist - prev tstep smoist
             delwb(k)=SUM((bwb(k,:,1)-(bwb(k,:,2)))*soil%zse)*1000.0
          END DO
       ELSE IF(MOD(REAL(ktau),2.0)==0.0) THEN    ! if even timestep
          bwb(:,:,2)=ssoil%wb
          DO k=1,mp           !  current smoist - prev tstep smoist
             delwb(k)=SUM((bwb(k,:,2)-(bwb(k,:,1)))*soil%zse)*1000.0
          END DO
       END IF
    END IF
    ! net water into soil (precip-(change in canopy water storage) 
       !  - (change in snow depth) - (surface runoff) - (deep drainage)
       !  - (evaporated water from vegetation and soil(excluding fevw, since
       !      it's included in change in canopy storage calculation))
       bal%wbal = REAL(met%precip - canopy%delwc - ssoil%snowd+ssoil%osnowd & 
            -ssoil%rnof1-ssoil%rnof2-(canopy%fevw+canopy%fevc + &
            canopy%fes/ssoil%cls)*dels/air%rlam - delwb,r_1)

       ! Canopy water balance: precip-change.can.storage-throughfall-evap+dew
       canopy_wbal=REAL(met%precip-canopy%delwc-canopy%through - &
            (canopy%fevw+MIN(canopy%fevc,0.0))*dels/air%rlam,r_1)
       
       IF(ktau>10) THEN
          DO j=1,mp	
             IF(ABS(canopy_wbal(j))>1e-4) THEN
                WRITE(*,*) 'Imbalance: ',canopy_wbal(j)
                WRITE(*,*) 'Timestep:',ktau, 'land point #:',j,'In mm:'
                WRITE(*,*) 'Precipitation:',met%precip(j),'canopy interception',&
                     canopy%wcint(j)
                WRITE(*,*) 'change in canopy water store:', canopy%delwc(j)
                WRITE(*,*) 'throughfall',canopy%through(j),'latent from wet canopy', &
                     canopy%fevw(j)*dels/air%rlam(j)
                WRITE(*,*) 'dew',MIN(canopy%fevc(j),0.0)*dels/air%rlam(j)
                WRITE(*,*) ''
                WRITE(*,*) 'Non-precip:', canopy%delwc(j)+canopy%through(j)+ & 
                     (canopy%fevw(j)+MIN(canopy%fevc(j),0.0))*dels/air%rlam(j)
                STOP 'Water balance failure within canopy.'
             END IF
          END DO
          ! Add current water imbalance to total imbalance (method 1 for water balance):
          bal%wbal_tot = bal%wbal_tot + bal%wbal
          ! Add to accumulation variables:
          bal%precip_tot = bal%precip_tot + met%precip
          bal%rnoff_tot = bal%rnoff_tot + ssoil%rnof1 + ssoil%rnof2
          bal%evap_tot = bal%evap_tot + &
               (canopy%fev+canopy%fes/ssoil%cls) * dels/air%rlam
          
       END IF
    
  END SUBROUTINE mass_balance
 !==================================================================================
  SUBROUTINE energy_balance(ktau,dels,met,rad,canopy,bal,ssoil)
    USE physical_constants
    INTEGER(i_d), INTENT(IN)     :: ktau ! time step
    REAL(r_1),INTENT(IN)         :: dels ! time step size
    TYPE (canopy_type),INTENT(IN):: canopy ! canopy variable data
    TYPE(met_type),INTENT(IN)    :: met  ! met data
    TYPE(radiation_type),INTENT(IN)   :: rad  ! met data
    TYPE (balances_type),INTENT(INOUT):: bal 
    TYPE (soil_snow_type),INTENT(IN)  :: ssoil ! soil data
    INTEGER :: j ! do loop counter

    ! SW absorbed + LW absorbed - (LH+SH+ghflux) should = 0
    bal%ebal = SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs &
         +met%fld-sboltz*emleaf*canopy%tv**4*(1-rad%transd)- &
         rad%flws*rad%transd -canopy%fev-canopy%fes/ssoil%cls &
         -canopy%fh -canopy%ga
    ! Add to cumulative balance:
    bal%ebal_tot = bal%ebal_tot + bal%ebal

  END SUBROUTINE energy_balance
!===============================================================================
     SUBROUTINE units_in(met,rad,dels)
       ! Changes units for CABLE text driver, if required, and 
       ! initialises several met variables.
       TYPE (met_type), INTENT(INOUT) :: met
       TYPE (radiation_type),INTENT(IN) :: rad
       REAL(r_1),INTENT(IN) :: dels
       REAL(r_1) :: temp_hum ! temporary humidity variable
       INTEGER :: i ! do loop counter
       
       ! Adjust rainfall units:
       IF(units%Rainf=='s') THEN
          met%precip = met%precip*dels
       ELSE IF(units%Rainf=='h') THEN
          met%precip = met%precip*dels/3600.0
       ELSE IF(units%Rainf=='t') THEN
          ! no change required
       ELSE
          WRITE(*,*) 'unknown rainfall units'
       END IF
       ! Adjust surface pressure units:
       IF(units%PSurf=='P') THEN
          met%pmb = met%pmb*0.01
       ELSE IF(units%PSurf=='h') THEN
          ! no change required
       ELSE
          WRITE(*,*) 'unknown pressure units'
       END IF
       ! Adjust air temperature units:
       IF(units%Tair=='K') THEN
          met%tc = met%tc - tfrz
       ELSE IF(units%Tair=='C') THEN
          ! no change required
       ELSE
          WRITE(*,*) 'unknown temperature units'
       END IF
       ! Create Kelvin temperature:
       met%tk = met%tc + tfrz
       ! Adjust humidity units:
       IF(units%Qair=='%') THEN
          DO i=1,mp
             temp_hum = met%qv(i)
             CALL rh_sh(temp_hum, met%tk(i), &
                  met%pmb(i),met%qv(i))
          END DO
       ELSE IF(units%Qair=='g') THEN
          ! no change required
       ELSE
          WRITE(*,*) 'unknown humidity units'
       END IF
       ! Adjust CO2 concentration units:
       IF(units%CO2air=='p') THEN
          met%ca = met%ca/1000000.0
       ELSE IF(units%CO2air=='m') THEN
          ! no change required
       ELSE
          WRITE(*,*) 'unknown temperature units'
       END IF
       
       ! Initialise other met variables:
       WHERE(met%tc<0.0) 
          met%precips = met%precip
       ELSEWHERE
          met%precips=0.0
       END WHERE
       met%tvair = met%tk 
       met%tvrad = met%tk 
       met%coszen = sinbet(met%doy, rad%latitude, met%hod)
       
       ! Check ranges are okay:
       ! Multiply acceptable Rainf ranges by time step size:
       ranges%Rainf = ranges%Rainf*dels ! range therefore depends on dels
       IF(ANY(met%fsd<ranges%SWdown(1)).OR.ANY(met%fsd>ranges%SWdown(2))) THEN
          WRITE(*,*)'SWdown out of specified ranges! Check units.'
          STOP
       ELSE IF(ANY(met%fld<ranges%LWdown(1)).OR.ANY(met%fld>ranges%LWdown(2)))THEN
          WRITE(*,*) 'LWdown out of specified ranges! Check units.'
          STOP
       ELSE IF(ANY(met%qv<ranges%Qair(1)).OR.ANY(met%qv>ranges%Qair(2))) THEN
          WRITE(*,*) 'Qair out of specified ranges! Check units.'
          STOP
       ELSE IF(ANY(met%precip<ranges%Rainf(1)).OR.ANY(met%precip>ranges%Rainf(2))) THEN
          WRITE(*,*) 'Rainf out of specified ranges! Check units'
          STOP
       ELSE IF(ANY(met%ua<ranges%Wind(1)).OR.ANY(met%ua>ranges%Wind(2))) THEN
          WRITE(*,*) 'Wind out of specified ranges! Check units'
          STOP
       ELSE IF(ANY(met%tk<ranges%Tair(1)).OR.ANY(met%tk>ranges%Tair(2))) THEN        
          WRITE(*,*)  'Tair out of specified ranges! Check units.'
          STOP
       ELSE IF(ANY(met%pmb<ranges%PSurf(1)).OR.ANY(met%pmb>ranges%PSurf(2))) THEN
          WRITE(*,*) 'PSurf out of specified ranges! Check units'
          STOP
       END IF

     END SUBROUTINE units_in

!===========================================================================
     SUBROUTINE rh_sh (relHum,tk,psurf,specHum)
       ! Converts relative humidity to specific humidity
       REAL, INTENT (IN)  :: psurf  ! surface pressure (hPa)
       REAL, INTENT (IN)  :: relHum ! relative humidity (%)
       REAL, INTENT (OUT) :: specHum ! specific humidity (kg/kg)
       REAL, INTENT (IN)  :: tk     ! air temp (K) 
       REAL :: es ! saturation vapour pressure
       REAL :: ws ! specific humidity at saturation
       es = svp (tk) ! saturation vapour pressure
       ws = 0.622 * es / (psurf - es) ! specific humidity at saturation
       specHum = (relHum/100.0) * ws ! specific humidity
     END SUBROUTINE rh_sh
!------------------------------------------------------------------------------------
     FUNCTION svp(tk) RESULT (F_Result)
       ! Calculate saturation vapour pressure
       REAL :: eilog
       REAL :: ewlog, ewlog2, ewlog3, ewlog4
       REAL :: F_Result
       REAL :: temp, tk
       REAL :: toot, toto, tsot
       temp = tk - 273.155
       IF (temp < -20.) THEN
          ! ice saturation
          toot = 273.16 / tk
          toto = 1. / toot
          eilog =   -9.09718 * (toot-1)                                      &
               -  3.56654 * (LOG (toot) / LOG (10.))                      &
               +  .876793 * (1-toto)                                      &
               +  (LOG (6.1071) / LOG (10.))
          F_Result = 10.**eilog
       ELSE
          tsot = 373.16 / tk
          ewlog = -7.90298 * (tsot-1) + 5.02808 * (LOG (tsot) / LOG (10.))
          ewlog2 =   ewlog                                                   &
               - 1.3816e-07 * (10**(11.344 * (1 - (1/tsot))) - 1)
          ewlog3 =    ewlog2                                                 &
               + .0081328 * (10**(-3.49149 * (tsot-1)) - 1)
          ewlog4 = ewlog3 + (LOG (1013.246) / LOG (10.))
          F_Result = 10.**ewlog4
       END IF
     END FUNCTION svp
!============================================================================
END MODULE checks_module
