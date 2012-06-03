
module canopy_module
   implicit none
   private
   public define_canopy

   contains
 
   subroutine define_canopy(bal,rad,rough,air,met,dels,ssoil,soil,veg, canopy)
      use define_types
      use define_dimensions
      use radiation_module
      use air_module
      use other_constants
      use physical_constants
      use photosynthetic_constants
      use cable_common_module   
      use cable_diag_module, only : cable_stat

      implicit none
      type (balances_type),intent(inout)  :: bal
      type (radiation_type), intent(inout):: rad
      type (roughness_type), intent(inout):: rough
      type (air_type), intent(inout)      :: air
      type (met_type), intent(inout)      :: met
      type (soil_snow_type), intent(inout):: ssoil
      type (soil_parameter_type), intent(inout)   :: soil
      type (veg_parameter_type), intent(inout)    :: veg
      type (canopy_type), intent(inout)   :: canopy
      real(r_1), intent(in)               :: dels ! integration time setp (s)
      
      !jhan: can be alloc/dealloc
      REAL(r_2), DIMENSION(mp)            :: gbvtop ! bnd layer cond. top leaf
      REAL(r_1), DIMENSION(mp)            :: rt0 ! turbulent resistance
      REAL(r_1), DIMENSION(mp)            :: ortsoil ! turb. resist. prev t-step
      REAL(r_1), DIMENSION(mp)            :: rt1usc ! eq. 3.53, SCAM manual, 1997
      REAL(r_1), DIMENSION(mp) :: tstar ! 
      REAL(r_1), DIMENSION(mp) :: zscrn !
      REAL(r_1), DIMENSION(mp) :: qstar !
      REAL(r_1), DIMENSION(mp) :: rsts  !
      REAL(r_1), DIMENSION(mp) :: qsurf !
      REAL(r_1), DIMENSION(mp) :: qtgnet !
      REAL(r_1), DIMENSION(mp) :: tss4 ! soil/snow temperature**4
      REAL(r_1), DIMENSION(mp) :: qstvair ! sat spec humidity at leaf temperature
      INTEGER(i_d) :: iter ! iteration #
      INTEGER(i_d) :: iterplus !

      REAL(r_1), DIMENSION(mp) :: xx ! delta-type func for sparse canopy limit, p20 SCAM manual
      REAL(r_1), DIMENSION(mp)  :: term1, term2, term3, term5
      REAL(r_1), DIMENSION(mp)  :: r_sc
      REAL(r_1), DIMENSION(mp)  :: zscl

      real(r_1), dimension(:), pointer :: cansat ! max canopy intercept. (mm)
      real(r_1), dimension(:), pointer :: dsx ! leaf surface vpd
      real(r_1), dimension(:), pointer :: fwsoil ! soil water modifier of stom. cond
      real(r_1), dimension(:), pointer :: tlfx ! leaf temp prev. iter (K)
      real(r_1), dimension(:), pointer :: tlfy ! leaf temp (K)
      real(r_2), dimension(:), pointer :: ecy ! lat heat fl dry big leaf
      real(r_2), dimension(:), pointer :: hcy ! veg. sens heat
      real(r_2), dimension(:), pointer :: rny ! net rad
      real(r_2), dimension(:,:), pointer :: gbhu ! forcedConvectionBndryLayerCond
      real(r_2), dimension(:,:), pointer :: gbhf ! freeConvectionBndryLayerCond
      real(r_2), dimension(:,:), pointer :: csx ! leaf surface CO2 concentration
      real(r_2), dimension(:), pointer :: ghwet  ! cond for heat for a wet canopy
      REAL(r_1), DIMENSION(mp)            :: pwet
      REAL(r_1), DIMENSION(mp)            :: dq ! sat sp
      REAL(r_1), DIMENSION(mp)            :: xx1
      REAL(r_1), DIMENSION(mp)  :: sum_rad_rniso 

      !jhan:local temporary vars
      real(r_1), dimension(mp)  :: ftemp,z_eff,psim_arg, psim_1, psim_2, rlower_limit  
      real  :: rt_min
      real(r_1), dimension(mp)  :: sum_rad_gradis 
      
      !  rml temporary arrays to write out air temp when leaf temp goes crazy
      integer, dimension(1) :: temp1, temp2
      integer :: j
      integer, save :: call_number =0
      call_number = call_number + 1
               
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('define_canopy')
      
         !jhan:in UM interface. Mk3l to follow
         if( .NOT. cable_runtime%um) &
            canopy%cansto =  canopy%oldcansto

         allocate(cansat(mp),gbhu(mp,mf))
         allocate(dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp))
         allocate(ecy(mp), hcy(mp), rny(mp))
         allocate(gbhf(mp,mf), csx(mp,mf))
         allocate(ghwet(mp))

         !jhan:never used
         !phenps = 1.0  

         ! BATS-type canopy saturation proportional to LAI:
         cansat = veg%canst1 * canopy%vlaiw
         !cansat = 0.

         !---compute surface wetness factor, update cansto, through
         call surf_wetness_fact( cansat, canopy, ssoil,veg,met, soil, dels )

         canopy%fevw_pot = 0.0
         canopy%gswx = 1e-3     ! default stomatal conuctance 
         gbhf = 1e-3     ! default free convection boundary layer conductance
         gbhu = 1e-3     ! default forced convection boundary layer conductance
         ssoil%evapfbl = 0.0

     
         ! Initialise in-canopy temperatures and humidity:
         csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
         met%tvair = met%tk
         met%qvair = met%qv
         canopy%tv = met%tvair
     
         CALL define_air (met, air)
         
         call qsatfjh(qstvair,met%tvair-tfrz,met%pmb)

         met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.0
         dsx = met%dva     ! init. leaf surface vpd
         
         tlfx = met%tk  ! initialise leaf temp iteration memory variable (K)
         tlfy = met%tk  ! initialise current leaf temp (K)
         
         ortsoil = ssoil%rtsoil
         ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
         tss4 = ssoil%tss**4
         canopy%fes = 0.
         canopy%fess = 0.
         canopy%fesp = 0.
         ssoil%potev = 0.
         canopy%fevw_pot = 0.
    
         CALL radiation(bal, soil, ssoil, veg, air, met, rad,canopy, dels)
    
     
         !jhan:mk3l to make global as well
         canopy%zetar(:,1) = zeta0 ! stability correction terms
         canopy%zetar(:,2) = zetpos + 1 

         DO iter = 1, niter

            ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
            ! resistances rt0, rt1 (elements of dispersion matrix):
            ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
            
            call comp_friction_vel()

            ! Turbulent aerodynamic resistance from roughness sublayer depth to reference height,
            ! x=1 if zref+disp>zruffs, 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
            ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
            xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
            !              correction  by Ian Harman to the 2nd psis term
            
            rt1usc = xx * (LOG(rough%zref_tq/MAX(rough%zruffs-rough%disp, rough%z0soilsn)) &
                  - psis( canopy%zetar(:,iter) ) &
                  + psis( canopy%zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn))/rough%zref_tq ) &
                  )/vonk
            
            rt_min = 5.      
            rt0 = max(rt_min,rough%rt0us / canopy%us)
     
            ! Aerodynamic resistance (sum 3 height integrals)/us
            ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
            
            rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
    
            !jhan:Eva uses .9 -> 1.
            !jhan:in latest vn doesn't use this block at all
            !WhHERE (ssoil%snowd > 0.1)
            !   ssoil%wetfac = 0.90
            !END WHhERE
            
            do j=1,mp
               if(canopy%vlaiw(j) > LAI_THRESH) then
                 ssoil%rtsoil(j) = rt0(j)
               ELSE
                 ssoil%rtsoil(j) = rt0(j) + rough%rt1(j)
               ENDif 
            enddo 
     
            ssoil%rtsoil = max(rt_min,ssoil%rtsoil)   
            
            do j=1,mp
               if( ssoil%rtsoil(j) > 2.*ortsoil(j) .OR. ssoil%rtsoil(j) < 0.5*ortsoil(j) ) &
                  ssoil%rtsoil(j) = MAX(rt_min,0.5*(ssoil%rtsoil(j) + ortsoil(j)))
            enddo 
     
            ! Vegetation boundary-layer conductance (mol/m2/s)
            ! prandt = kinematic viscosity/molecular diffusivity
            ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
            do j=1,mp
               if(canopy%vlaiw(j) > lai_thresh) then
                  gbvtop(j) = air%cmolar(j)*apol * air%visc(j) / prandt / veg%dleaf(j) *       &
                        (canopy%us(j) / MAX(rough%usuh(j),1.e-6) * &
                        veg%dleaf(j) / air%visc(j) )**0.5 * prandt**(1.0/3.0) / veg%shelrb(j)
                  gbvtop(j) = MAX (0.05,gbvtop(j) )      ! for testing (BP aug2010)
                  
                  ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
                  gbhu(j,1) = gbvtop(j)*(1.0-EXP(-canopy%vlaiw(j)*(0.5*rough%coexp(j)+rad%extkb(j) ))) / &
                        (rad%extkb(j)+0.5*rough%coexp(j))
                  gbhu(j,2) = (2.0/rough%coexp(j))*gbvtop(j)*  &
                        (1.0-EXP(-0.5*rough%coexp(j)*canopy%vlaiw(j)))-gbhu(j,1)
               endif 
            enddo 
     
            rny = SUM(rad%rniso,2) ! init current estimate net rad
            hcy = 0.0              ! init current estimate lat heat
            ecy = rny - hcy        ! init current estimate lat heat
     
            sum_rad_rniso = SUM(rad%rniso,2)
            CALL dryLeaf(dels,rad,rough,air,met,veg,canopy, soil,ssoil, dsx,fwsoil,    &
                  tlfx,tlfy, ecy,hcy, rny, gbhu, gbhf, csx, cansat, ghwet,  iter )
     
            CALL wetLeaf(dels,rad,rough,air,met,veg,canopy,cansat, tlfy, gbhu,gbhf, ghwet )

!  rml 8/2/11 check leaf temperature
! need an lai check too since Antarctica can get this cold
!           if (minval(tlfy).lt.200.or.maxval(tlfy).gt.400) then
!             write(6,*) 'rmlcheck tlfy: ',minval(tlfy),maxval(tlfy)
!             temp1 = minloc(tlfy); temp2 = maxloc(tlfy)
!             write(6,*) 'rmlcheck: min/max loc',minloc(tlfy),maxloc(tlfy)
!             write(6,*) 'rmlcheck: met/tvair',met%tvair(temp1(1)), met%tvair(temp2(1))
!             write(6,*) 'rmlcheck: veg type',veg%iveg(temp1(1)), veg%iveg(temp2(1))
!             write(6,*) 'rmlcheck: lai ',canopy%vlaiw(temp1(1)), canopy%vlaiw(temp2(1))
!             stop 'Error in cable_canopy: tlf'
!           endif
     
            ! Calculate latent heat from vegetation:
            ! Calculate sensible heat from vegetation:
            ! Calculate net rad absorbed by canopy:
            canopy%fev = REAL(canopy%fevc + canopy%fevw,r_1)
            ftemp = (1.0 - canopy%fwet) *  REAL(hcy,r_1) + canopy%fhvw
            canopy%fhv = real(ftemp,r_1) 
            ftemp= (1.0-canopy%fwet)*REAL(rny,r_1)+canopy%fevw+canopy%fhvw
            canopy%fnv = real(ftemp,r_1)
     
            ! canopy radiative temperature is calculated based on long-wave radiation balance
            sum_rad_gradis = SUM(rad%gradis,2)
            do j=1,mp
               if ( canopy%vlaiw(j) > LAI_THRESH .and. rough%hruff(j) > rough%z0soilsn(j) ) then
                  rad%lwabv(j) = capp * rmair * (tlfy(j)-met%tk(j))* sum_rad_gradis(j) 

                  canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))*sboltz*emleaf)+met%tk(j)**4)**0.25
                  !canopy%tv(j) = ( (rad%qcan(j,1,1)+rad%qcan(j,1,2) + rad%qcan(j,2,1)+rad%qcan(j,2,2) &
                  !    + (1.0-rad%transd(j))*(rad%flws(j) + met%fld(j)) - canopy%fev(j) - canopy%fhv(j) ) / &
                  !                                  (2.*(1.0-rad%transd(j))*sboltz * emleaf ) ) ** 0.25
               ELSE! sparse canopy
                 canopy%tv(j) = met%tk(j)
               ENDif 
            enddo 
           

            ! Calculate net rad to soil:
            canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
                  sboltz*canopy%tv**4 - emsoil*sboltz* tss4


            !jhan:?? to do here. UM uses ssoil%qstss in calc of met%qvair as well 
            !jhan:ssoil%qstss include in mk3l?
            ! Saturation specific humidity at soil/snow surface temperature:
    
            !call qsatfjh(ssoil%qstss,ssoil%tss-tfrz,met%pmb)
            call qsatfjh(ssoil%qstss,ssoil%tss-tfrz,met%pmb)

            !if(cable_runtime%um) then
            dq = ssoil%qstss - met%qv
            ssoil%potev =  Humidity_deficit_method(dq,ssoil%qstss ) 
            !else
            !   !--- uses %ga from previous timestep    
            !   ssoil%potev =  Penman_Monteith(canopy%ga) 
            !endif

            ! Soil latent heat:
            call latent_heat_flux()
  
            ! Calculate soil sensible heat:
            !canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
            canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil


            call within_canopy( gbhu, gbhf )

         
            !if(cable_runtime%um) then
            call qsatfjh(ssoil%qstss,ssoil%tss-tfrz,met%pmb)
            dq = ssoil%qstss - met%qvair
            ssoil%potev =  Humidity_deficit_method(dq,ssoil%qstss ) 
            !else
            !   canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux + ssoil%isflag*canopy%sghflux
            !   ssoil%potev =  Penman_Monteith(canopy%ghflux) 
            !endif

         
            ! Soil latent heat:
            call latent_heat_flux()

            ! Soil sensible heat:
            canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
            canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
            
            ! Set total latent heat:
            canopy%fe = canopy%fev + canopy%fes
            
            ! Set total sensible heat:
            canopy%fh = canopy%fhv + canopy%fhs
            !jhan:{ prev. in UM only
            !---diagnostic purposes
            do j=1,mp
               if (ssoil%potev(j) .ge. 0.) then
                  ssoil%potev(j) = max(0.00001,ssoil%potev(j))
               else
                  ssoil%potev(j) = min(-0.0002,ssoil%potev(j))
               endif
               !jhan:Eva has added this block
               if (canopy%fevw_pot(j) .ge. 0.) then
                  canopy%fevw_pot(j) = max(0.000001,canopy%fevw_pot(j))
               else
                  canopy%fevw_pot(j) = min(-0.002,canopy%fevw_pot(j))
               endif
            enddo 


            canopy%rnet = canopy%fnv + canopy%fns  
            canopy%epot = ((1.-rad%transd)*canopy%fevw_pot + rad%transd*ssoil%potev) * dels/air%rlam  
            ! convert to mm/day
            !rlower_limit = canopy%fevw_pot + ssoil%potev
            rlower_limit = canopy%epot * air%rlam / dels
            !do j=1,mp
            !   if( abs(rlower_limit(j) ) < 1.e-30) & 
            !      rlower_limit = rlower_limit + 1.e-6 
            !      !rlower_limit = rlower_limit + 1.e-29 
            !enddo
            canopy%wetfac_cs = max(0., min(1.0,canopy%fe / rlower_limit ))
            
!            do j=1,mp
!               if ( canopy%wetfac_cs(j) .le. 0. )  &
!                  canopy%wetfac_cs(j) = max(0.,min(1., &
!                     max(canopy%fev(j)/canopy%fevw_pot(j),canopy%fes(j)/ssoil%potev(j))))
!            enddo 

            call update_zetar()

      END DO           ! do iter = 1, niter

      canopy%cduv = canopy%us * canopy%us / (max(met%ua,umin))**2

      !jhan:prev. in UM only
      !---diagnostic purposes
      canopy%gswx_T = rad%fvlai(:,1)/max(LAI_THRESH,canopy%vlaiw(:))*canopy%gswx(:,1) &
             + rad%fvlai(:,2)/max(LAI_THRESH,canopy%vlaiw(:))*canopy%gswx(:,2)
      !jhan:Evachanged
      canopy%gswx_T = max(1.e-05,canopy%gswx_T )
             
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
       ftemp = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) + &
            psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk
   
       ! Calculate screen temperature:
       canopy%tscrn = met%tk - tfrz - tstar * ftemp
   
       ! Calculate radiative/skin temperature; at this stage old soil temperature is used
       !rad%trad = ( (1.-rad%transd)*canopy%tv**4 + rad%transd * ssoil%tss**4 )**0.25

      ! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman
      !jhan:Eva's mods in r2411_cableum_17 July 2011 reset rghlai
      !jhan:in merged version (canopy%)rghlai is init. in _roughness()
      !do j=1,mp
      !   if(ssoil%snowd(j) .lt. 0.001.and.veg%iveg(j) .ne. 1) &
      !      canopy%rghlai(j) = min(3.,canopy%vlaiw(j))
      !enddo
       
      term1=0.
      term2=0.
      term5=0.
      term3 = 0. ! Work around for Intel compiler problem with nested whres
      r_sc = 0.
      zscl = max(rough%z0soilsn,2.0)
 
      !---assume screen temp of bareground if all these conditions are not met
      do j=1,mp
         if ( canopy%vlaiw(j) > LAI_THRESH .and. rough%hruff(j) > 0.01) then
            if ( rough%disp(j)  > 0.0 ) then
               !jhan:Eva uses new var
               term1(j) = EXP(2*csw*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
               term2(j) = EXP(2*csw*canopy%rghlai(j)*(1-rough%disp(j)/rough%hruff(j)))
               term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)
            endif
            !jhan:Eva uses new var
            term3(j) = a33**2*ctl*2*csw*canopy%rghlai(j)
 
            if( zscl(j) < rough%disp(j) ) then
                r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) * &
                ( exp(2*csw*canopy%rghlai(j)) - term1(j) ) / term3(j)

            elseif ( rough%disp(j) <= zscl(j) .and. zscl(j) < rough%hruff(j) ) then
                r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) / term3(j)

            elseif ( rough%hruff(j) <= zscl(j) .and. zscl(j) <  rough%zruffs(j) ) then
                r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) * &
                ( zscl(j) - rough%hruff(j) ) /  (a33**2*ctl*rough%hruff(j) )

            elseif (zscl(j) >= rough%zruffs(j) ) then
                r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +  &
                  ( log( (zscl(j) - rough%disp(j))/MAX(rough%zruffs(j)-rough%disp(j), rough%z0soilsn(j)) ) &
                    - psis1( (zscl(j) - rough%disp(j)) / (rough%zref_tq(j)/canopy%zetar(j,iterplus)) )  &
                    + psis1( (rough%zruffs(j) - rough%disp(j)) / (rough%zref_tq(j)/canopy%zetar(j,iterplus)) )  &
                  ) / vonk
            endif
           !jhan:Eva uses 
           canopy%tscrn(j) = ssoil%tss(j) + (met%tk(j) - ssoil%tss(j)) * min(1.,r_sc(j) /            &
                  max(1.,rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j))) - tfrz 
         endif  
      enddo  
     
      call qsatfjh(rsts,canopy%tscrn,met%pmb)
        
      qtgnet = rsts * ssoil%wetfac - met%qv
      do j=1,mp
         if (qtgnet(j) .gt. 0. ) then
            qsurf(j) = rsts(j) * ssoil%wetfac(j)
         else
            qsurf(j) = 0.1*rsts(j)*ssoil%wetfac(j) + 0.9*met%qv(j)
         endif
         canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

         !jhan:Eva uses 
         if ( canopy%vlaiw(j) >LAI_THRESH .and. rough%hruff(j) > 0.01) &
            canopy%qscrn(j) =  qsurf(j) + (met%qv(j) - qsurf(j)) * min(1.,r_sc(j) / &
                            max(1.,rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j)))
      enddo 


      ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
      canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
           dels * 1.0e3 / (rhow*air%rlam)
      ! Add dewfall to canopy water storage:
      canopy%cansto = canopy%cansto + canopy%dewmm
      ! rml 28/2/11 moved from a few lines further down to improve water balance (BP)
      ! Modify canopy water storage for evaporation:
      canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw,r_1))*dels &
        *1.0e3/(rhow*air%rlam), 0.0)

      ! Calculate canopy water storage excess:
      canopy%spill=max(0.0, canopy%cansto-cansat)
  
      ! Move excess canopy water to throughfall:
      canopy%through = canopy%through + canopy%spill
      ! Initialise 'throughfall to soil' as 'throughfall from canopy'; snow may absorb
      canopy%precis = max(0.,canopy%through)
      ! Update canopy storage term:
      canopy%cansto=canopy%cansto - canopy%spill
      ! Modify canopy water storage for evaporation:
      !canopy%cansto = max(canopy%cansto-max(0.0,canopy%fevw)*dels*1.0e3/ &
      !(rhow*air%rlam), 0.0)
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
           +capp*rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009
      bal%wetbal=canopy%fevw+canopy%fhvw-SUM(rad%rniso,2)*canopy%fwet &
           +capp*rmair*(tlfy-met%tk)*SUM(rad%gradis,2)*canopy%fwet  ! YP nov2009
  
      DEALLOCATE(cansat,gbhu)
      DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
      DEALLOCATE(ecy, hcy, rny)
      DEALLOCATE(gbhf, csx)
      DEALLOCATE(ghwet)

      CONTAINS


      subroutine comp_friction_vel()
      real(r_1), dimension(mp)  :: lower_limit, rescale
            psim_1 = psim(canopy%zetar(:,iter)) 
            
               rescale = vonk * MAX(met%ua,umin)
               z_eff = rough%zref_uv / rough%z0m
            
            psim_arg = canopy%zetar(:,iter) / z_eff 
            !---fix for compiler limitation. bitwise reproducable whilst we  
            !---we know it to 11th decimal. psim_arg typically of a few 
            !psim_arg = nint(psim_arg * 1.e11)*1.e-11

            psim_2 = psim( psim_arg )
                   
            lower_limit = rescale / ( LOG(z_eff) - psim_1 + psim_2 )

            canopy%us = MAX(1.e-6, lower_limit )
      end subroutine comp_friction_vel


   function Penman_Monteith( ground_H_flux ) result(ssoilpotev)
      implicit none
      real(r_1), intent(in), dimension(mp)  :: ground_H_flux
      real(r_1), dimension(mp)  :: ssoilpotev 
      REAL(r_1), DIMENSION(mp) :: sss ! var for Penman-Monteith soil evap
      REAL(r_1), DIMENSION(mp) :: cc1 ! var for Penman-Monteith soil evap
      REAL(r_1), DIMENSION(mp) :: cc2 ! var for Penman-Monteith soil evap
      REAL(r_1), DIMENSION(mp) :: qsatfvar
      integer :: j
         ! Penman-Monteith formula
         sss=air%dsatdk
         cc1=sss/(sss+air%psyc )
         cc2=air%psyc /(sss+air%psyc )
         
         call qsatfjh(qsatfvar,met%tvair-tfrz,met%pmb)

         ssoilpotev = cc1 * (canopy%fns - ground_H_flux) + &
         cc2 * air%rho * air%rlam*(qsatfvar  - met%qvair)/ssoil%rtsoil
      return
   end function Penman_Monteith


   ! method alternative to P-M formula above
   function Humidity_deficit_method(dq,qstss ) result(ssoilpotev)
      implicit none
      real(r_1), dimension(mp)    :: ssoilpotev 
      REAL(r_1), DIMENSION(mp)    :: dq ! sat spec hum diff.
      REAL(r_1), DIMENSION(mp)    :: qstss !dummy var for compilation 
      integer :: j

       do j=1,mp
         !if(ssoil%snowd(j) > 1.0) dq(j) = max( -0.1e-3, dq(j))
         if(ssoil%snowd(j)>1.0 .or. ssoil%tgg(j,1).eq.tfrz)  dq(j) = max( -0.1e-3, dq(j))
       enddo 
      ssoilpotev =air%rho * air%rlam * dq /ssoil%rtsoil
      return
   end function Humidity_deficit_method

 
   subroutine latent_heat_flux() 
      use cable_common_module
      implicit none
      real(r_1), dimension(mp) ::  swilt_eff, frescale
      real(r_1), dimension(mp) ::  flower_limit, fupper_limit
      integer :: j

         ! Soil latent heat:
         !canopy%fes= ssoil%wetfac * ssoil%potev
         !jhan:Eva uses %fess
         canopy%fess= ssoil%wetfac * ssoil%potev
         where (ssoil%potev < 0. ) canopy%fess = ssoil%potev
         ! Reduce soil evap due to presence of puddle
         pwet = max(0.,min(0.2,ssoil%pudsto/max(1.,ssoil%pudsmx)))
         canopy%fess = canopy%fess * (1.-pwet)

         !jhan: is this necessary anymore       
         !--- decrease wilting point in UM 
!         if( cable_runtime%um) then
               swilt_eff = soil%swilt/2.0    
!         else
!               swilt_eff = soil%swilt    
!         endif
         !jhan:NB %evapfbl
         !frescale = soil%zse(1) * 1000.  - ssoil%evapfbl(:,1) * air%rlam / dels         
         frescale = soil%zse(1) * 1000. * air%rlam / dels         

         do j=1,mp
            !jhan:Eva uses %fess
            if(ssoil%snowd(j) < 0.1 .and. canopy%fess(j) .gt. 0. ) then

               !flower_limit(j) = REAL(ssoil%wb(j,1),r_1)-soil%swilt(j)
               flower_limit(j) = REAL(ssoil%wb(j,1),r_1)-swilt_eff(j)
               fupper_limit(j) = max(0._r_2,flower_limit(j)*frescale(j)-ssoil%evapfbl(j,1)*air%rlam(j)/dels)
               canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))
               
               fupper_limit(j) = REAL(ssoil%wb(j,1)-ssoil%wbice(j,1),r_1) * frescale(j)
               canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))
            END if
   
            ssoil%cls(j)=1.
            if (ssoil%snowd(j) >= 0.1 .and. ssoil%potev(j) > 0.) then
               ssoil%cls(j) = 1.1335
               canopy%fess(j) = min( (ssoil%wetfac(j)*ssoil%potev(j)),                           &
                     ssoil%snowd(j)/dels*air%rlam(j)*ssoil%cls(j))
            ENDif

         enddo 
         
         ! Evaporation form soil puddle
         canopy%fesp = min(ssoil%pudsto/dels*air%rlam,max(pwet*ssoil%potev,0.))
         canopy%fes = canopy%fess + canopy%fesp

      return
   end subroutine latent_heat_flux

subroutine within_canopy( gbhu, gbhf )
  implicit none
   real(r_2), intent(in), dimension(:,:) :: gbhu ! forcedConvectionBndryLayerCond
   real(r_2), intent(in), dimension(:,:):: gbhf ! freeConvectionBndryLayerCond
      REAL(r_1), DIMENSION(mp)            :: rrsw ! recipr. stomatal resistance for water
      REAL(r_1), DIMENSION(mp)            :: rrbw ! recipr. leaf boundary layer resistance for water
       REAL(r_1), DIMENSION(mp) :: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL(r_1), DIMENSION(mp) :: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL(r_1), DIMENSION(mp) :: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL(r_1), DIMENSION(mp) :: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL(r_1), DIMENSION(mp) :: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL(r_1), DIMENSION(mp) :: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      real(r_1)  :: lower_limit, upper_limit
      !real(r_1), dimension(mp)  :: lower_limit, upper_limit
      integer :: j
      
            !jhan:mk3l DID typecast again, and rrbw=1/rbw 
            rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT 
            ! leaf stomatal resistance for water
            rrsw = sum(canopy%gswx,2)/air%cmolar ! MJT
            do j=1,mp
            if(veg%meth(j) > 0 .and. canopy%vlaiw(j) > LAI_THRESH .and. &
                      rough%hruff(j) > rough%z0soilsn(j)) then
               !   use the dispersion matrix (DM) to find the air temperature 
               !   and specific humidity 
               !   (Raupach, Finkele and Zhang 1997, pp 17)
               ! leaf boundary layer resistance for water
               ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmah(j) = (rt0(j)+rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j)) &
                    + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
               ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmbh(j) = (-air%rlam(j)/capp)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
               ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)* &
                    (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*capp)
               ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmae(j) = (-air%epsi(j)*capp/air%rlam(j))*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
               ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmbe(j) = (rt0(j)+ssoil%wetfac(j)*rough%rt1(j))* &
                      ((1.+air%epsi(j))*rrsw(j) + rrbw(j))+(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
               ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
               dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)* &
                      (canopy%fev(j) + canopy%fes(j))/(air%rho(j)*air%rlam(j))
     
               ! Within canopy air temperature:
               met%tvair(j) = met%tk(j) + (dmbe(j)*dmch(j)-dmbh(j)*dmce(j))/(dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
                 
               !---set limits for comparisson
               lower_limit =  min( ssoil%tss(j), met%tk(j)) - 5.0
               upper_limit =  max( ssoil%tss(j), met%tk(j)) + 5.0
      
               !--- tvair within these limits
               met%tvair(j) = max(met%tvair(j) , lower_limit)
               met%tvair(j) = min(met%tvair(j) , upper_limit)
     
               ! recalculate using canopy within temperature
               !     whhere (veg%meth .eq. 0 )
               met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j))/(dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
               met%qvair(j) = max(0.0,met%qvair(j))

               !---set limits for comparisson
               lower_limit =  min( ssoil%qstss(j), met%qv(j)) 
               upper_limit =  max( ssoil%qstss(j), met%qv(j))
 
                  !--- qvair within these limits
               met%qvair(j) =  max(met%qvair(j),lower_limit)
               met%qvair(j) =  min(met%qvair(j), upper_limit)

               ! Saturated specific humidity in canopy:
               call qsatfjh2(qstvair(j),met%tvair(j)-tfrz,met%pmb(j))
               ! Saturated vapour pressure deficit in canopy:
               met%dva(j) = (qstvair(j) - met%qvair(j)) *  rmair/rmh2o * met%pmb(j) * 100.
            END if 
         enddo 
     
   return
end subroutine within_canopy

   subroutine update_zetar()
      implicit none
      integer :: j
            ! monin-obukhov stability parameter zetar=zref/l
            !        recompute zetar for the next iteration, except on last iteration
            IF (iter < niter) THEN ! dont compute zetar on the last iter
               iterplus = max(iter+1,2)
               canopy%zetar(:,iterplus) = -(vonk*grav*rough%zref_tq*(canopy%fh+0.07*canopy%fe))/ &
                     (air%rho*capp*met%tk*canopy%us**3)
               ! case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
               IF (niter == 2) THEN
                  canopy%zetar(:,2) = zetmul * canopy%zetar(:,2)
                  do j=1,mp
                     if ( (met%fsd(j,1)+met%fsd(j,2))  ==  0.0 ) &
                        canopy%zetar(j,2) = 0.5 * canopy%zetar(j,2)
                  enddo
               END IF
               !     constrain zeta to zetpos and zetneg (set in param0)
               canopy%zetar(:,iterplus) = min(zetpos,canopy%zetar(:,iterplus))        ! zetar too +
               canopy%zetar(:,iterplus) = max(zetneg,canopy%zetar(:,iterplus))        ! zetar too -
            END IF ! (iter < niter) 
   end subroutine update_zetar



!jhan:these are all NOT contained subrs in mk3l
    !--------------------------------------------------------------------------
    FUNCTION qsatf(j,tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1)           :: r    ! result; sat sp humidity
      integer, intent(in) :: j 
      r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END FUNCTION qsatf

    !---------------------------------------------------------
    subroutine qsatfjh(var,tair,pmb) 
      REAL(r_1), INTENT(IN), dimension(mp) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN), dimension(mp) :: pmb  ! pressure PMB (mb)
      REAL(r_1), intent(out), dimension(mp)           :: var    ! result; sat sp humidity
      integer :: j
      
      do j=1,mp
  
         var(j) = (rmh2o/rmair) * (tetena*EXP(tetenb*tair(j)/(tetenc+tair(j)))) / pmb(j)
      enddo   
    END subroutine qsatfjh

    !---------------------------------------------------------

    subroutine qsatfjh2(var,tair,pmb) 
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1), intent(out)           :: var    ! result; sat sp humidity
      
  
         var = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END subroutine qsatfjh2


     FUNCTION psim(zeta) result(r)
     ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
     ! computes integrated stability function psim(z/l) (z/l=zeta)
     ! for momentum, using the businger-dyer form for unstable cases
     ! and the Beljaars and Holtslag (1991) form for stable cases.
     USE math_constants
     REAL(r_1), INTENT(IN),dimension(mp)    :: zeta
     REAL(r_1), dimension(mp)                 :: r
     REAL(r_1), dimension(mp)                 :: x
     REAL(r_1), PARAMETER                     :: gu = 16.0
     REAL(r_1), PARAMETER                     :: gs = 5.0
     REAL(r_1), dimension(mp)                 :: z
     REAL(r_1), dimension(mp)                 :: stable
     REAL(r_1), dimension(mp)                 :: unstable

     REAL(r_1), PARAMETER      :: a = 1.0
     REAL(r_1), PARAMETER      :: b = 0.667
     REAL(r_1), PARAMETER      :: c = 5.0
     REAL(r_1), PARAMETER      :: d = 0.35
     !      x = (1.0 + gu*abs(zeta))**0.25
     !      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
     !           + pi_c*0.5, -gs*zeta, zeta < 0.0)
     z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
     stable = -gs*zeta
     ! replace with Beljaars and Holtslag (1991) for stable
     stable = -a*zeta - b*(zeta - c/d)*exp( -d*zeta) - b*c/d
     x      = (1.0 + gu*abs(zeta))**0.25
     unstable = alog((1.0+x*x)*(1.0+x)**2/8) - 2.0*atan(x) + pi_c*0.5
     r   = z*stable + (1.0-z)*unstable
    END FUNCTION psim


     ELEMENTAL FUNCTION psis(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      REAL(r_1), INTENT(IN)     :: zeta
      REAL(r_1)                 :: r
      REAL(r_1), PARAMETER      :: gu = 16.0
      REAL(r_1), PARAMETER      :: gs = 5.0

      REAL(r_1)                 :: stzeta
!      REAL(r_1)                 :: ustzeta
      REAL(r_1)                 :: z
      REAL(r_1)                 :: y
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable

      REAL(r_1), PARAMETER      :: a = 1.0
      REAL(r_1), PARAMETER      :: b = 0.667
      REAL(r_1), PARAMETER      :: c = 5.0
      REAL(r_1), PARAMETER      :: d = 0.35

      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
      stable = -gs*zeta
      ! replace with Beljaars and Holtslag (1991) for stable
      !stzeta = min(2.0,max(0.,zeta))
      stzeta = max(0.,zeta)
      stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
                b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
      y      = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis


    !---------------------------------------------------------
    FUNCTION psis1(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      REAL(r_1), INTENT(IN)     :: zeta
      REAL(r_1)                 :: r
      REAL(r_1), PARAMETER      :: gu = 16.0
      REAL(r_1), PARAMETER      :: gs = 5.0
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable
      REAL(r_1)                 :: stzeta

      REAL(r_1)                 :: z
      REAL(r_1)                 :: y
      !REAL(r_1)                 :: stable
      !REAL(r_1)                 :: unstable

      REAL(r_1), PARAMETER      :: a = 1.0
      REAL(r_1), PARAMETER      :: b = 0.667
      REAL(r_1), PARAMETER      :: c = 5.0
      REAL(r_1), PARAMETER      :: d = 0.35

      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
      stable = -gs*zeta
      ! replace with Beljaars and Holtslag (1991) for stable
      stzeta = max(0.,zeta)
      stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
                b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.

      y      = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis1

    !---------------------------------------------------------
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
      REAL(r_1), INTENT(IN)     :: rpconst
      REAL(r_1), INTENT(IN)     :: rpcoef
      REAL(r_1), INTENT(IN)     :: tair
      REAL(r_1)                 :: z
      z = rpconst * exp(rpcoef * tair)
    END FUNCTION rplant
    !---------------------------------------------------------

    !---------------------------------------------------------
  SUBROUTINE wetLeaf(dels,rad,rough,air,met,veg,canopy,cansat, tlfy, gbhu, gbhf,ghwet)

    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
   real(r_1),intent(in), dimension(:) :: tlfy ! leaf temp (K)
   real(r_2), intent(in), dimension(:,:) :: gbhu ! forcedConvectionBndryLayerCond
   real(r_2), intent(in), dimension(:,:) :: gbhf ! freeConvectionBndryLayerCond
      real(r_2), intent(out), dimension(:) :: ghwet  ! cond for heat for a wet canopy

      real(r_1),intent(in), dimension(:) :: cansat ! max canopy intercept. (mm)
   ! assuming the temperature of wet leaf is equal that of dry leaf ="tlfy"
    REAL(r_1), INTENT(IN)     :: dels ! integration time step (s)
    REAL(r_1), DIMENSION(mp)  :: ccfevw ! limitation term for
                                        ! wet canopy evaporation rate
    REAL(r_1), DIMENSION(mp)  :: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp)  :: ghrwet ! wet canopy cond: heat & thermal rad
    REAL(r_1), DIMENSION(mp)  :: sum_gbh, sum_rad_rniso, sum_rad_gradis 
    REAL(r_1), DIMENSION(mp)  :: xx1,xx2,xx3
    integer :: j

    ghwet = 1.0e-3
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    sum_gbh = SUM((gbhu+gbhf),2)
    sum_rad_rniso = SUM(rad%rniso,2)
    sum_rad_gradis = SUM(rad%gradis,2)

    do j=1,mp
      if(canopy%vlaiw(j) > LAI_THRESH) then
       ! VEG SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
       ! calculate total thermal resistance, rthv in s/m
          ghwet(j) = 2.0   * sum_gbh(j) 
          gwwet(j) = 1.075 * sum_gbh(j) 
          ghrwet(j) = sum_rad_gradis(j) + ghwet(j)
          ! Calculate fraction of canopy which is wet:
          canopy%fwet(j) = MAX(0.0,MIN(1.0,0.8*canopy%cansto(j)/MAX(cansat(j),0.01)))
          ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
          ! to avoid excessive evaporation:
          ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )
   
          canopy%fevw(j) = MIN(canopy%fwet(j) * (air%dsatdk(j)*(sum_rad_rniso(j)- &
                        capp*rmair*(met%tvair(j)-met%tk(j))*sum_rad_gradis(j)) &
                        + capp*rmair*met%dva(j)*ghrwet(j)) &
                            / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j)), ccfevw(j))
          !jhan:prev. in UM only
          canopy%fevw_pot(j) = (air%dsatdk(j)* (sum_rad_rniso(j)- &
                        capp*rmair*(met%tvair(j)-met%tk(j))*sum_rad_gradis(j)) &
                        + capp*rmair*met%dva(j)*ghrwet(j)) &
                            / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )
          ! alculate sens heat from wet canopy:
          !canopy%fhvw = (canopy%fwet*(SUM(rad%rniso,2)-capp*rmair *              &
          !(met%tvair(:)-met%tk(:))*sum(rad%gradis,2))-canopy%fevw)*ghwet/ghrwet
   
          canopy%fhvw(j) = canopy%fwet(j)*(sum_rad_rniso(j) -capp*rmair*(tlfy(j)-met%tk(j))* &
                        sum_rad_gradis(j))  - canopy%fevw(j)
          xx1(j) = canopy%fhvw(j)
          canopy%fhvw(j) = sum_rad_rniso(j) - canopy%fevc(j) - canopy%fevw(j)-  &
                           (1.0 - canopy%fwet(j)) *  REAL(hcy(j),r_1) 
          xx2(j) = canopy%fhvw(j)
          canopy%fhvw(j) =  xx1(j)

      ENDif 
    enddo        
  END SUBROUTINE wetLeaf

  END SUBROUTINE define_canopy


   subroutine surf_wetness_fact( cansat, canopy, ssoil,veg, met, soil, dels )
      use cable_common_module
      use define_dimensions
      use define_types
      use physical_constants, only : tfrz
      implicit none
      type (veg_parameter_type), INTENT(INOUT)    :: veg
      type (soil_snow_type), intent(inout):: ssoil
      type (soil_parameter_type), intent(inout)   :: soil
      type (canopy_type), INTENT(INOUT)   :: canopy
      type (met_type), INTENT(INOUT)   :: met
      real(r_1), intent(in)               :: dels ! integration time setp (s)
      real(r_1),intent(in), dimension(:) :: cansat ! max canopy intercept. (mm)
      real(r_1), dimension(mp)  :: lower_limit, upper_limit,ftemp
      real(r_1), dimension(mp) ::  swilt_eff
      integer :: j

         ! Rainfall variable is limited so canopy interception is limited,
         ! used to stabilise latent fluxes.
         ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
         upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 ) 
         ftemp =MIN(met%precip-met%precip_sn, upper_limit )
         ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
         lower_limit = cansat - canopy%cansto
         upper_limit = max(lower_limit, 0.0) 
         canopy%wcint = MERGE( MIN(upper_limit, ftemp), 0.0,     &
               ftemp > 0.0  .AND. met%tk > tfrz)  !EAK, 09/10

         ! Define canopy throughfall (100% of precip if temp < 0C, see above):
         canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn , &
              MAX(0.0, met%precip - met%precip_sn - canopy%wcint) )  ! EK nov2007
     
         ! Add canopy interception to canopy storage term:
         canopy%cansto = canopy%cansto + canopy%wcint
     
!         if( cable_runtime%um) then
               swilt_eff = soil%swilt/2.0
!         else
!               swilt_eff = soil%swilt
!         endif

         ! Calculate fraction of canopy which is wet:
         !jhan:rm hardwire
         canopy%fwet   = max(0.0,min(0.9,0.8*canopy%cansto/MAX(cansat,0.01)))
     
         ssoil%wetfac = MAX(1.e-6, MIN(1.0, &
              (REAL(ssoil%wb(:,1),r_1) - swilt_eff) / (soil%sfc - swilt_eff)))
        
            do j=1,mp
               if( ssoil%wbice(j,1) > 0. ) &
                   ssoil%wetfac(j) = ssoil%wetfac(j)*max(0.5,1.-min(0.2,(ssoil%wbice(j,1)/ssoil%wb(j,1))**2))
               if ( ssoil%snowd(j) > 0.1) ssoil%wetfac(j) = 0.9
               if ( veg%iveg(j) == 16 .and. met%tk(j) >= tfrz + 5.) ssoil%wetfac(j) = 1.0
               if ( veg%iveg(j) == 16 .and. met%tk(j) < tfrz + 5.) ssoil%wetfac(j) = 0.7
            enddo 
            
    
         ! owetfac introduced to reduce sharp changes in dry regions,
         ! especially in offline runs in which there may be discrepancies between
         ! timing of precip and temperature change (EAK apr2009)
         ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)

         !jhan:Eva has modified this
         ! Temporay fixer for accounting the reduction of soil evap due to freezing
         !WERE ( ssoil%wbice(:,1) > 0.0 ) ! Prevents divide by zero at glaciated
         !                                ! points whhere wb and wbice=0.
          ! ssoil%wetfac = ssoil%wetfac &
          !              * REAL(((1.0-ssoil%wbice(:,1)/ssoil%wb(:,1)))**2,r_1)
         !END WERE
    
      return
   end subroutine surf_wetness_fact


   SUBROUTINE dryLeaf(dels,rad,rough,air,met,veg,canopy, soil,ssoil, dsx,fwsoil, tlfx,  &
         tlfy, ecy, hcy, rny, gbhu, gbhf, csx, cansat, ghwet, iter )
      use define_dimensions
      use define_types
      use cable_common_module
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
      implicit none

      TYPE (radiation_type), INTENT(INOUT):: rad
      TYPE (roughness_type), INTENT(INOUT):: rough
      TYPE (air_type), INTENT(INOUT)      :: air
      TYPE (met_type), INTENT(INOUT)      :: met
      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
      TYPE (canopy_type), INTENT(INOUT)   :: canopy
      type (soil_parameter_type), intent(inout)   :: soil
      type (soil_snow_type), intent(inout):: ssoil
      real(r_1), intent(inout), dimension(:) :: dsx ! leaf surface vpd
      real(r_1), intent(inout), dimension(:):: fwsoil ! soil water modifier of stom. cond
      real(r_1),intent(inout), dimension(:) :: tlfx ! leaf temp prev. iter (K)
      real(r_1),intent(inout), dimension(:) :: tlfy ! leaf temp (K)
      real(r_2),intent(inout), dimension(:) :: ecy ! lat heat fl dry big leaf
      real(r_2),intent(inout), dimension(:) :: hcy ! veg. sens heat
      real(r_2),intent(inout), dimension(:) :: rny 
      real(r_2), intent(in), dimension(:,:) :: gbhu ! forcedConvectionBndryLayerCond
      real(r_2), intent(inout), dimension(:,:) :: gbhf ! freeConvectionBndryLayerCond
      real(r_2), intent(inout), dimension(:,:) :: csx ! leaf surface CO2 concentration
      real(r_1),intent(in), dimension(:) :: cansat
      real(r_2), intent(out), dimension(:) :: ghwet  ! cond for heat for a wet canopy

      integer,intent(in) :: iter
    
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
      REAL(r_2), DIMENSION(mp)  :: ecx_t ! lat. hflux big leaf
      REAL(r_2), DIMENSION(mp)  :: hcx ! sens heat fl big leaf prev iteration
      REAL(r_2), DIMENSION(mp)  :: rnx ! net rad prev timestep
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
      INTEGER(i_d)   :: k  ! iteration count
      INTEGER(i_d)   :: kk  ! iteration count
      REAL(r_1), DIMENSION(mp,ms)  :: oldevapfbl
      REAL(r_1), DIMENSION(mp)  :: sum_rad_rniso
      real(r_1), dimension(mp)  :: sum_rad_gradis

    REAL(r_1), DIMENSION(mp)  :: ccfevw ! limitation term for
                                        ! wet canopy evaporation rate
    REAL(r_1), DIMENSION(mp)  :: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp)  :: ghrwet ! wet canopy cond: heat & thermal rad
    REAL(r_1), DIMENSION(mp)  :: sum_gbh



      INTEGER(i_d)   :: i,j !loop counts - kdcorbin, 09/10

      !jhan:lower_limit = IS this doable outside stability loop
      real(r_1), dimension(:,:), pointer :: gswmin ! min stomatal conductance
      real(r_1), dimension(mp,2)  ::  gsw_term, lower_limit2  ! local temp var 
      REAL(r_1), DIMENSION(mp,mf) :: frac42 ! 2D frac4


      allocate( gswmin(mp,mf ))

      ! Soil water limitation on stomatal conductance:
      if( iter ==1) then
         if(cable_user%FWSOIL_SWITCH == 'standard') then
            call fwsoil_calc_std( fwsoil, soil, ssoil, veg) 
         elseif (cable_user%FWSOIL_SWITCH == 'non-linear extrapolation') then
            !EAK, 09/10 - replace linear approx by polynomial fitting
            call fwsoil_calc_non_linear(fwsoil, soil, ssoil, veg) 
         elseif(cable_user%FWSOIL_SWITCH == 'Lai and Ktaul 2000') then
            call fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssoil, veg) 
         else
            STOP 'fwsoil_switch failed.'
         endif
      endif

      ! weight min stomatal conductance by C3 an C4 plant fractions
      frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants

      gsw_term = gsw03 * (1. - frac42) + gsw04 * frac42
      lower_limit2 = rad%scalex * (gsw03 * (1. - frac42) + gsw04 * frac42)
      gswmin = max(1.e-6,lower_limit2)
            

      gw = 1.0e-3 ! default values of conductance
      gh = 1.0e-3
      ghr= 1.0e-3
      rdx = 0.0
      anx = 0.0
      rnx = SUM(rad%rniso,2)
      abs_deltlf = 999.0
      !jhan:Eva uses
      !jh:canopy%fevc = 0.0
  
  
      gras = 1.0e-6
      an_y= 0.0
      hcx = 0.0              ! init sens heat iteration memory variable
      hcy = 0.0
      rdy = 0.0
      ecx = SUM(rad%rniso,2) ! init lat heat iteration memory variable
      tlfxx = tlfx
      psycst(:,:) = SPREAD(air%psyc,2,mf)
      canopy%fevc = 0.0
      ssoil%evapfbl = 0.0

    ghwet = 1.0e-3
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    sum_gbh = SUM((gbhu+gbhf),2)
    sum_rad_rniso = SUM(rad%rniso,2)
    sum_rad_gradis = SUM(rad%gradis,2)

      DO kk=1,mp
         IF(canopy%vlaiw(kk) <= LAI_THRESH) THEN
            rnx(kk) = 0.0 ! intialise
            ecx(kk) = 0.0 ! intialise
            ecy(kk) = ecx(kk) ! store initial values
            abs_deltlf(kk)=0.0
            rny(kk) = rnx(kk) ! store initial values
       ! calculate total thermal resistance, rthv in s/m

         END IF
      ENDDO
      deltlfy = abs_deltlf
      k = 0

      !kdcorbin, 08/10 - doing all points all the time
      DO WHILE (k < maxiter)
         k = k + 1

         DO i=1,mp
            IF (canopy%vlaiw(i) > LAI_THRESH .AND. abs_deltlf(i) > 0.1) Then
          ghwet(i) = 2.0   * sum_gbh(i)
          gwwet(i) = 1.075 * sum_gbh(i)
          ghrwet(i) = sum_rad_gradis(i) + ghwet(i)
          ! Calculate fraction of canopy which is wet:
          canopy%fwet(i) = MAX(0.0,MIN(1.0,0.8*canopy%cansto(i)/MAX(cansat(i),0.01)))
          ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
          ! to avoid excessive evaporation:
          ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )
   
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
                  !jhan: here mk3l usess (i,:) {
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
               !jhan:}  
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
    
                !jhan: here mk3l usess (i,:) {
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
               !jhan:sort out braces here: x/() + y
               xleuning(i,1) = (fwsoil(i) / (csx(i,1)-co2cp3))  &
                         * ((1.0-veg%frac4(i)) * a1c3 / (1.0+dsx(i)/d0c3) &
                             + veg%frac4(i)    * a1c4 / (1.0+dsx(i)/d0c4))
               xleuning(i,2) = (fwsoil(i) / (csx(i,2)-co2cp3))  &
                         * ((1.0-veg%frac4(i)) * a1c3 / (1.0+dsx(i)/d0c3) &
                             + veg%frac4(i)    * a1c4 / (1.0+dsx(i)/d0c4))
    
               !jhan:}  
            ENDIF
         ENDDO !i=1,mp
   
         call photosynthesis(csx(:,:),SPREAD(cx1(:),2,mf), &
                SPREAD(cx2(:),2,mf),gswmin(:,:),rdx(:,:), &
                vcmxt3(:,:),vcmxt4(:,:),vx3(:,:),vx4(:,:),xleuning(:,:), &
                rad%fvlai(:,:),SPREAD(abs_deltlf,2,mf),anx(:,:))

         DO i=1,mp
            IF (canopy%vlaiw(i) > LAI_THRESH .AND. abs_deltlf(i) > 0.1) Then
               do kk=1,mf
                  if(rad%fvlai(i,kk)>LAI_THRESH) then
                     csx(i,kk) = met%ca(i) - rgbwc*anx(i,kk) / (gbhu(i,kk) + gbhf(i,kk))
                     csx(i,kk) = MAX(1.0e-4,csx(i,kk))

                     canopy%gswx(i,kk) = MAX(1.e-3, gswmin(i,kk) +                         &
                       MAX(0.0,rgswc*xleuning(i,kk)*anx(i,kk)))

                     !Recalculate conductance for water:
                     gw(i,kk) = 1.0/(1.0/canopy%gswx(i,kk)   +  1.0/(1.075*(gbhu(i,kk)+gbhf(i,kk))))
                     gw(i,kk) = MAX(gw(i,kk),0.00001)

                     !Modified psychrometric constant (Monteith and Unsworth, 1990)
                     psycst(i,kk) = air%psyc(i)*REAL(ghr(i,kk)/gw(i,kk),r_1)
                  endif
               enddo

               !ecx_t(i) = (air%dsatdk(i)*(rad%rniso(i,1) &
               ecx(i) = (air%dsatdk(i)*(rad%rniso(i,1) &
                     - capp*rmair*(met%tvair(i)-met%tk(i)) &
                     * rad%gradis(i,1)) + capp*rmair*met%dva(i)*ghr(i,1)) &
                     / (air%dsatdk(i)+psycst(i,1)) &
                     + (air%dsatdk(i)*(rad%rniso(i,2) &
                     -  capp*rmair*(met%tvair(i)-met%tk(i))*rad%gradis(i,2)) &
                     + capp*rmair*met%dva(i)*ghr(i,2)) &
                     / (air%dsatdk(i)+psycst(i,2))
               !ecx(i) = (1.0-canopy%fwet(i))* ecx_t(i)  + &
               !         MIN(canopy%fwet(i) * (air%dsatdk(i)*(sum_rad_rniso(i)- &
               !         capp*rmair*(met%tvair(i)-met%tk(i))*sum_rad_gradis(i)) &
               !         + capp*rmair*met%dva(i)*ghrwet(i)) &
               !             / (air%dsatdk(i)+air%psyc(i)*ghrwet(i)/gwwet(i)), ccfevw(i))

              !canopy%fevw(i) = MIN(canopy%fwet(i) * (air%dsatdk(i)*(sum_rad_rniso(i)- &
              !          capp*rmair*(met%tvair(i)-met%tk(i))*sum_rad_gradis(i)) &
              !          + capp*rmair*met%dva(i)*ghrwet(i)) &
              !              / (air%dsatdk(i)+air%psyc(i)*ghrwet(i)/gwwet(i)), ccfevw(i))

              ! canopy%fhvw(i) = 0.

               !IF (ecx_t(i) > 0.0 .AND. canopy%fwet(i) < 1.0) Then
               IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) Then
                  !evapfb(i) = (1.0-canopy%fwet(i))*REAL(ecx_t(i),r_1)*dels/air%rlam(i)
                  evapfb(i) = (1.0-canopy%fwet(i))*REAL(ecx(i),r_1)*dels/air%rlam(i)
                  DO kk = 1,ms
                     !ssoil%evapfbl(i,kk) = MIN(evapfb(i)*veg%froot(i,kk), &
                     !    MAX(0.0,MIN(REAL(ssoil%wb(i,kk),r_1) &
                     !      -soil%swilt(i),REAL(ssoil%wb(i,kk)  &
                     !      -1.05*ssoil%wbice(i,kk),r_1))) &
                     !      * soil%zse(kk)*1000.0)
                     !jhan:Eva uses
                     ssoil%evapfbl(i,kk) = MIN(evapfb(i)*veg%froot(i,kk), &
                         MAX(0.0,REAL(ssoil%wb(i,kk),r_1) - 1.1*soil%swilt(i)) &
                           * soil%zse(kk)*1000.0)
                  ENDDO

                  canopy%fevc(i) = SUM(ssoil%evapfbl(i,:))*air%rlam(i)/dels
       
                  ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                  !ecx_t(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))

               ENDIF
               !ecx(i) = (1.0-canopy%fwet(i))* ecx_t(i) + canopy%fevw(i)

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

         ! Whhere leaf temp change b/w iterations is significant, and
         ! difference is smaller than the previous iteration, store results:
         !jhan:Eva uses instead
         DO i=1,mp
            if (abs_deltlf(i) < ABS(deltlfy(i)) ) then
               !jh:of WhHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
               deltlfy(i) = deltlf(i)
               tlfy(i) = tlfx(i)
               rny(i) = rnx(i)
               hcy(i) = hcx(i)
               ecy(i) = ecx(i)
               rdy(i,1) = rdx(i,1)
               rdy(i,2) = rdx(i,2)
               an_y(i,1) = anx(i,1)
               an_y(i,2) = anx(i,2)
               ! save last values calculated for ssoil%evapfbl
               oldevapfbl(i,1) = ssoil%evapfbl(i,1)
               oldevapfbl(i,2) = ssoil%evapfbl(i,2)
               oldevapfbl(i,3) = ssoil%evapfbl(i,3)
               oldevapfbl(i,4) = ssoil%evapfbl(i,4)
               oldevapfbl(i,5) = ssoil%evapfbl(i,5)
               oldevapfbl(i,6) = ssoil%evapfbl(i,6)
            endif 
            if(abs_deltlf(i) > 0.1) &
               ! after 4 iterations, take mean value of current & previous estimates
               ! as the next estimate of leaf temperature, to avoid oscillation
               tlfx(i) = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx(i) + &
                     (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx(i)
        
            IF(k==1) THEN
               ! take the first iterated estimates as the defaults
               tlfy(i) = tlfx(i)
               rny(i) = rnx(i)
               hcy(i) = hcx(i)
               ecy(i) = ecx(i)
               rdy(i,:) = rdx(i,:)
               an_y(i,:) = anx(i,:)
               ! save last values calculated for ssoil%evapfbl
               oldevapfbl(i,:) = ssoil%evapfbl(i,:)
            END IF
         END DO !over mp 


      END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < maxiter)

      ! dry canopy flux
      canopy%fevc = (1.0-canopy%fwet) * ecy
      !canopy%fevc = (1.0-canopy%fwet) * ecx_t
      !canopy%fevc = ecx_t

      ! Recalculate ssoil%evapfbl as ecy may not be updated with the ecx
      ! calculated in the last iteration.
      ! DO NOT use simple scaling as there are times that ssoil%evapfbl is zero.
      ! ** ssoil%evapfbl(i,:) = ssoil%evapfbl(i,:) * ecy(i) / ecx(i) **
      DO i = 1, mp
         IF (ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
            IF (ABS(ecy(i)-ecx(i)) > 1.0e-6 ) THEN
               IF (ABS(canopy%fevc(i) - (SUM(oldevapfbl(i,:))*air%rlam(i)/dels)) &
                     > 1.0e-4 ) THEN
                  PRINT *, 'Error! oldevapfbl not right.', ktau_gl, i
                  PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
                  PRINT *, 'or in mm = ',ecx(i)*(1.0-canopy%fwet(i))/air%rlam(i)*dels, &
                                         ecy(i)*(1.0-canopy%fwet(i))/air%rlam(i)*dels
                  PRINT *,'fevc = ',canopy%fevc(i),SUM(oldevapfbl(i,:))*air%rlam(i)/dels
                  PRINT *, 'fwet = ', canopy%fwet(i)
                  PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
                  PRINT *, 'ssoil%evapfbl before rescaling: ', ssoil%evapfbl(i,:)
                  STOP
               ELSE
                  ssoil%evapfbl(i,:) = oldevapfbl(i,:)
               END IF
            END IF
         END IF
      END DO

      canopy%frday = 12.0 * SUM(rdy, 2)
      canopy%fpn = -12.0 * SUM(an_y, 2)
      canopy%evapfbl = ssoil%evapfbl


      return
   END SUBROUTINE dryLeaf

  SUBROUTINE photosynthesis(csxz,cx1z,cx2z,gswminz,rdxz,vcmxt3z,vcmxt4z, &
                    vx3z,vx4z,xleuningz,vlaiz,deltlfz,anxz)
      use define_dimensions
      use other_constants, only : LAI_THRESH
      !use math_constants
      !use physical_constants
      use photosynthetic_constants, only : rgswc
      implicit none
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

    ! rgswc - inherited from canopy_module's USE photosynthetic_constants
    ! mp - inherited from canopy_module
    ! mf - inherited from canopy_module
   
   DO i=1,mp
      IF (sum(vlaiz(i,:)) .gt. LAI_THRESH) Then
      DO j=1,mf
         IF (vlaiz(i,j) .gt. LAI_THRESH .AND. deltlfz(i,j) .gt. 0.1) Then

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
    ENDIF
   ENDDO
     
  END SUBROUTINE photosynthesis
    !---------------------------------------------------------



    FUNCTION ej3x(parx,x) result(z)
      use define_dimensions
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
      REAL(r_1), INTENT(IN)     :: parx
      REAL(r_1), INTENT(IN)     :: x
      REAL(r_1)                 :: z
      z = max(0.0, &
           0.25*((alpha3*parx+x-sqrt((alpha3*parx+x)**2 - &
           4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
    END FUNCTION ej3x
    !---------------------------------------------------------
    FUNCTION ej4x(parx,x) result(z)
      use define_dimensions
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
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
      use define_dimensions
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
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
      use define_dimensions
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
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
      use define_dimensions
      use other_constants
      use math_constants
      use physical_constants
      use photosynthetic_constants!, only : gsw03
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants

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

      subroutine fwsoil_calc_std(fwsoil, soil, ssoil, veg) 
      use define_types
      use define_dimensions
         implicit none
      type (soil_snow_type), intent(inout):: ssoil
      type (soil_parameter_type), intent(inout)   :: soil
      type (veg_parameter_type), intent(inout)    :: veg
         real(r_1), intent(out), dimension(:):: fwsoil ! soil water modifier of stom. cond
         REAL(r_1), DIMENSION(mp) :: rwater ! soil water availability

         rwater = MAX(1.0e-4_r_2, &
              SUM(veg%froot * MAX(0.024,MIN(1.0_r_2,ssoil%wb - &
                   SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))
        
            fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
            
         return   
      end subroutine fwsoil_calc_std 


      subroutine fwsoil_calc_non_linear(fwsoil, soil, ssoil, veg) 
      use define_types
      use define_dimensions
         implicit none
      type (soil_snow_type), intent(inout):: ssoil
      type (soil_parameter_type), intent(inout)   :: soil
      type (veg_parameter_type), intent(inout)    :: veg
         real(r_1), intent(out), dimension(:):: fwsoil ! soil water modifier of stom. cond
         REAL(r_1), DIMENSION(mp) :: rwater ! soil water availability
         REAL(r_1), DIMENSION(mp,3)          :: xi, ti, si
      integer :: j
 
         rwater = MAX(1.0e-4_r_2, &
              SUM(veg%froot * MAX(0.024,MIN(1.0_r_2,ssoil%wb - &
                   SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))
         fwsoil = 1.

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
         do j=1,mp
            if (rwater(j) < soil%sfc(j) - 0.02) &
               fwsoil(j) = max(0.,min(1., ti(j,1)*si(j,1) + &
                             ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))
         enddo

         return   
      end subroutine fwsoil_calc_non_linear 


      ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
      subroutine fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssoil, veg) 
      use define_types
      use define_dimensions
         implicit none
      type (soil_snow_type), intent(inout):: ssoil
      type (soil_parameter_type), intent(inout)   :: soil
      type (veg_parameter_type), intent(inout)    :: veg
         real(r_1), intent(out), dimension(:):: fwsoil ! soil water modifier of stom. cond
         REAL(r_1), DIMENSION(mp) :: rwater ! soil water availability
         INTEGER(i_d)   :: ns
         REAL(r_1), parameter ::rootgamma = 0.01   ! (19may2010)
         REAL(r_1), DIMENSION(mp)  :: dummy, normFac
         !--- local level dependent rwater 
         REAL(r_1), DIMENSION(mp,ms)  :: frwater
    
            fwsoil(:) = 1.0e-4
            normFac(:) = 0.0
            do ns=1,ms
               dummy(:) = rootgamma/max(1.0e-3,ssoil%wb(:,ns)-soil%swilt(:))
               frwater(:,ns) =MAX(1.0e-4,((ssoil%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
                     ** dummy)
               fwsoil(:) = min(1.0,max(fwsoil(:),frwater(:,ns)))
               normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)
             enddo
         return
      end subroutine fwsoil_calc_Lai_Ktaul








END MODULE canopy_module
