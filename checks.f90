MODULE checks_module
  USE canopy_module
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
     ! possible output/assimilatable variables
     REAL(r_1), DIMENSION(2) :: Qh = (/-1000.0,1000.0/)    ! W/m^2
     REAL(r_1), DIMENSION(2) :: Qle = (/-1000.0,1000.0/)   ! W/m^2
     REAL(r_1), DIMENSION(2) :: Qg = (/-1000.0,1000.0/)    ! W/m^2   
     REAL(r_1), DIMENSION(2) :: SWnet = (/0.0,1250.0/)     ! W/m^2
     REAL(r_1), DIMENSION(2) :: LWnet = (/-500.0,510.0/)   ! W/m^2 
     REAL(r_1), DIMENSION(2) :: Evap = (/-0.0003,0.00035/)      
     REAL(r_1), DIMENSION(2) :: Ewater = (/-0.0003,0.0003/)
     REAL(r_1), DIMENSION(2) :: ESoil = (/-0.0003,0.0003/)     
     REAL(r_1), DIMENSION(2) :: Tveg = (/-0.0003,0.0003/)    
     REAL(r_1), DIMENSION(2) :: Ecanop = (/-0.0003,0.0003/)   
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
     REAL(r_1), DIMENSION(2) :: NEE = (/-50.0,20.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: NPP = (/-20.0,50.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: GPP = (/-10.0,50.0/) ! umol/m2/s 
     REAL(r_1), DIMENSION(2) :: AutoResp = (/-50.0,20.0/) ! umol/m2/s
     REAL(r_1), DIMENSION(2) :: HeteroResp = (/-50.0,20.0/) ! umol/m2/s   
     ! parameters:
     REAL(r_1), DIMENSION(2) :: bch = (/2.0,15.0/)  
     REAL(r_1), DIMENSION(2) :: latitude = (/-90.0,90.0/)   
     REAL(r_1), DIMENSION(2) :: c3 = (/0.0,1.0/)         
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
     REAL(r_1), DIMENSION(2) :: rpcoef = (/0.8,1.5/)
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
    REAL(r_1), DIMENSION(:,:,:),allocatable, SAVE :: bwb ! volumetric soil moisture
    REAL(r_1), DIMENSION(mp) :: delwb ! change in soilmoisture b/w tsteps
    REAL(r_1), DIMENSION(mp) :: into_soil ! moisture into soil
    REAL(r_1), DIMENSION(mp) :: canopy_wbal ! canopy water balance
    TYPE (balances_type),INTENT(INOUT):: bal 
    INTEGER :: k ! do loop counter
    
    IF(ktau==1) THEN
       allocate ( bwb(mp,ms,2) )
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
       bal%wbal = met%precip - canopy%delwc - ssoil%snowd+ssoil%osnowd & 
            -ssoil%rnof1-ssoil%rnof2-(canopy%fevw+canopy%fevc + &
            canopy%fes/ssoil%cls)*dels/air%rlam - delwb

       ! Canopy water balance: precip-change.can.storage-throughfall-evap+dew
       canopy_wbal=met%precip-canopy%delwc-canopy%through - &
            (canopy%fevw+MIN(canopy%fevc,0.0))*dels/air%rlam

       ! Water into soil: (canopy throughfall) - (change in snow water) -
       ! (runoff and deep drainage) - (dry canopy transpiration) - (soil latent)
       into_soil = canopy%through - ssoil%snowd+ssoil%osnowd & 
            -ssoil%rnof1-ssoil%rnof2-(MAX(canopy%fevc,0.0) + canopy%fes)* & 
            dels/air%rlam
       
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
                CALL abort('Water balance failure within canopy.')
             END IF
!!$             IF(ABS(bal%wbal(j))>=5e-4) THEN
!!$                WRITE(*,*) '****',(bal%wbal), 'Excess in water balance:'
!!$                WRITE(*,*) 'Timestep:',ktau, 'land point #:',j
!!$                print*, 'net water into soil:',into_soil(j),'soil moisture change:',delwb(j)
!!$                print*, 'canopy water balance:',canopy_wbal(j)
!!$                print*,'fev:',canopy%fev(j)*dels/air%rlam(j),'fevc:',canopy%fevc(j) &
!!$                     *dels/air%rlam(j),'fevw:',canopy%fevw(j)*dels/air%rlam(j)
!!$                print*, 'fes:',canopy%fes(j)*dels/air%rlam(j),'dew:',canopy%dewmm(j)
!!$                print*, 'airT:',met%tc(j),'fe:',(canopy%fev(j)+canopy%fes(j))* &
!!$                     dels/air%rlam(j)
!!$                print*, 'fes:',canopy%fes(j)*dels/air%rlam(j)
!!$                print*, 'precip', met%precip(j),'canopy throughfall',canopy%through(j), &
!!$                     'delwc',canopy%delwc(j)
!!$                print*, 'osd:',ssoil%osnowd(j), 'snowd:',ssoil%snowd(j), &
!!$                     'Diff:',ssoil%snowd(j)-ssoil%osnowd(j)
!!$                print*, 'runoff:',ssoil%rnof1(j),'deepd:',ssoil%rnof2(j)
!!$                CALL abort('water balance failure.')
!!$             END IF
!!$             ELSE
!!$                IF(water_dump) THEN
!!$                   ! Dump excess water into bottom layer soil:
!!$                   ssoil%wb(j,6)=ssoil%wb(j,6)+(delwat(j)-wbal(j))/(soil%zse(6)*1000)
!!$    
!!$                   ! Recalculate this timestep's bwb:
!!$                   IF(MOD(REAL(ktau),2.0)==1.0) THEN         ! if odd timestep
!!$                      bwb(:,:,1)=ssoil%wb
!!$                      do k=1,mp           ! wbal = current smoist - prev tstep smoist
!!$                         wbal(k)=SUM((bwb(k,:,1)-(bwb(k,:,2)))*soil%zse)*1000
!!$                      end do
!!$                   ELSE IF(MOD(REAL(ktau),2.0)==0.0) THEN    ! if even timestep
!!$                      bwb(:,:,2)=ssoil%wb
!!$                      do k=1,mp           ! wbal = current smoist - prev tstep smoist
!!$                         wbal(k)=SUM((bwb(k,:,2)-(bwb(k,:,1)))*soil%zse)*1000
!!$                      end do
!!$                   END IF
!!$                   ! Check dump okay:
!!$                   IF(ABS((delwat(j)-wbal(j)))>5e-5) CALL abort('Error in water dump (checks.f90)')
!!$                   ! Update soilsnow soil moisture variable:
!!$                   wb(:,6)=UNPACK(ssoil%wb(:,6),land,wb(:,6))
!!$                   
!!$                !   IF(ktau>38600) print*, ktau, 'post', delwat(j)-wbal(j),'wb',wb
!!$
!!$                END IF
            
          END DO
          ! Add current water imbalance to total imbalance:
          bal%wbal_tot = bal%wbal_tot + bal%wbal
          ! Add to accumulation variables:
          bal%precip_tot = bal%precip_tot + met%precip
          bal%rnoff_tot = bal%rnoff_tot + ssoil%rnof1 + ssoil%rnof2
          bal%evap_tot = bal%evap_tot + &
               (canopy%fev+canopy%fes/ssoil%cls) * dels/air%rlam
       END IF
    
  END SUBROUTINE mass_balance
  
  SUBROUTINE energy_balance(ktau,dels,met,rad,canopy,bal,ssoil)
    INTEGER(i_d), INTENT(IN)     :: ktau ! time step
    REAL(r_1),INTENT(IN)         :: dels ! time step size
    TYPE (canopy_type),INTENT(IN):: canopy ! canopy variable data
    TYPE(met_type),INTENT(IN)    :: met  ! met data
    TYPE(radiation_type),INTENT(IN)   :: rad  ! met data
    TYPE (balances_type),INTENT(INOUT):: bal 
    TYPE (soil_snow_type),INTENT(IN)  :: ssoil ! soil data
    REAL(r_1), DIMENSION(mp)     :: ebalcan,ebalcan2 ! energy balances
    REAL(r_1), DIMENSION(mp)     :: SWbal ! SW rad balance
    INTEGER :: j ! do loop counter

    ! SW absorbed + LW absorbed - (LH+SH+ghflux) should = 0
    bal%ebal = SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs &
         +met%fld-sboltz*emleaf*canopy%tv**4*(1-rad%transd)- &
        ! rad%flws*rad%transd - canopy%fe -canopy%fh -canopy%ghflux   
      !   rad%flws*rad%transd - canopy%fe -canopy%fh -canopy%ga
         rad%flws*rad%transd -canopy%fev-canopy%fes/ssoil%cls &
         -canopy%fh -canopy%ga

    ebalcan2=bal%drybal+bal%wetbal

    ! distribution of fnv amongst fluxes
    ebalcan=canopy%fnv-canopy%fhv-canopy%fev

    ! will not hold when fes updated for snow evap
   ! ebalsoil=canopy%fns-canopy%fhs-canopy%fes-canopy%ghflux

    ! SW absorbed by canopy + abs by soil - (SWin-SW reflected) should =0
    SWbal=SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs &
         - met%fsd*(1-(rad%albedo(:,1)+rad%albedo(:,2))/2)

    ! Stop if imbalances.
    IF(ktau>10) THEN
          DO j=1,mp	
!!$            IF(ABS(ebalcan(j))>1e-4) CALL abort('canopy energy partitioning failure.')
!!$             !  IF(ABS(ebalsoil(j))>1e-4) CALL abort('soil energy balance failure.')
!!$             IF(ABS(SWbal(j))>1e-4) THEN
!!$                print*, 'ktau:',ktau, 'mp:',j
!!$                print*, 'fsd:',met%fsd(1),'albedo:',rad%albedo(1,:)
!!$                print*, 'qcanSW:',rad%qcan(1,:,1),rad%qcan(1,:,2)
!!$                print*, 'qssabs', rad%qssabs
!!$                CALL abort('SW balance failure')
!!$             END IF
!!$             IF(ABS(ebalcan2(j))>1e-4) THEN
!!$                print*, 'ktau:',ktau, 'mp:',j
!!$                CALL abort('canopy energy balance failure')
!!$             END IF
!!$             IF(ABS(bal%ebal(j))>1e-3)  THEN
!!$                WRITE(*,*) '**** Energy balance problem:'
!!$                WRITE(*,*) 'Timestep:',ktau, 'land point #:',j
!!$                WRITE(*,*) 'enegry balance:', bal%ebal(j)
!!$                WRITE(*,*) 'SW absorbed:',SUM(rad%qcan(j,:,1))+ &
!!$                     SUM(rad%qcan(j,:,2))+rad%qssabs(j)
!!$                WRITE(*,*) 'LW absorbed:',met%fld(1) - sboltz*emleaf* &
!!$                     canopy%tv**4*(1-rad%transd) - rad%flws*rad%transd
!!$                WRITE(*,*) 'sum of energy fluxes:', canopy%fe+canopy%fh+canopy%ga
!!$               
!!$                CALL abort('General energy imbalance ')
!!$             END IF

!!$             ELSE
!!$                IF(energy_dump) THEN
!!$                   ! Add energy to soil sensible heat:
!!$                   canopy%fh(j) = canopy%fh(j) + bal%ebal(j)
!!$                   canopy%fhs(j) = canopy%fhs(j) + bal%ebal(j)
!!$                   ! Recalculate energy balance:
!!$                   bal%ebal = met%fsd*(1-(rad%albedo(:,1)+rad%albedo(:,2))/2) &
!!$                        +met%fld(1)-sboltz*emleaf*canopy%tv**4*(1-rad%transd) &
!!$                        -rad%flws*rad%transd &
!!$                        -canopy%fe -canopy%fh -canopy%ghflux   
!!$                   ! Check okay
!!$                   IF(ABS(bal%ebal(j))>5e-4) CALL abort('Ebal dump failed!')
!!$                END IF
            
          END DO
          bal%ebal_tot = bal%ebal_tot + bal%ebal
       END IF
     END SUBROUTINE energy_balance

END MODULE checks_module
