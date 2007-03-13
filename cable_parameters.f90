!!$ cable_parameters.f90
!!$
!!$ Parameter input module for CABLE land surface scheme offline driver; 
!!$
!!$ Gab Abramowitz 2006 CSIRO Marine and Atmospheric
!!$ Research/ Macquarie University; gabsun@gmail.com
!!$
!!$ This subroutine reads default parameter sets for CABLE from a global 
!!$ grid, using 13 vegetation types and 13 soil types.

MODULE abort_module
  IMPLICIT NONE
CONTAINS
  SUBROUTINE abort(message)
    CHARACTER(LEN=*), INTENT(IN)	:: message
    WRITE(*, *) message
    STOP 1
  END SUBROUTINE abort
END MODULE abort_module

MODULE parameter_module
  USE define_types
  USE abort_module
  IMPLICIT NONE
  REAL(r_1), POINTER,DIMENSION(:) :: latitude, longitude
  INTEGER(i_d),POINTER,DIMENSION(:) :: gdpt ! gridpoint number (default params)

CONTAINS
  !---------------------------------------------------------------------------
  SUBROUTINE default_params(met,air,ssoil,veg,bgc, &
       soil,canopy,rough,rad,sum_flux,bal,logn)     
    ! gets some initial conditions as well as parameters
    IMPLICIT NONE
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    INTEGER(i_d),INTENT(IN) :: logn     ! log file unit number
    INTEGER(i_d) :: numGdpt=13824 ! total number of default parameter grid gridpoints
    INTEGER(i_d) :: v_type ! vegetation type for current gridcell
    INTEGER(i_d) :: s_type ! soil type for current gridcell
    REAL(r_1),POINTER :: lat_temp(:),lon_temp(:) ! temporary lat lon values
    REAL(r_1)    :: lon2 ! temporary lat lon values
    REAL(r_1)    :: grid_lat, grid_lon ! gridcell centre for site
    LOGICAL      :: gridFlag ! for reporting issues with default parameters
    INTEGER(i_d) :: ioerror ! input error integer
    INTEGER(i_d),POINTER, DIMENSION(:) :: oldgdpt ! in case gridcell is sea
    INTEGER(i_d) :: max_it = 81 ! maximum # iterations to find gridqsuare default params
    INTEGER(i_d) :: nvegt ! Number of vegetation types
    REAL(r_1),POINTER :: tgg_temp(:,:),wb_temp(:,:) ! temp soil temp + moist
    REAL(r_1) :: distance ! used to find correct grid cell
    CHARACTER(LEN=2),POINTER :: land_temp(:) ! land/sea point read in (logical)
    CHARACTER(LEN=2) :: land ! land/sea point (logical)
    INTEGER(i_d),POINTER :: isoil(:), iveg(:) ! temporary soil and veg type
    CHARACTER(LEN=70), DIMENSION(13) :: veg_desc
    CHARACTER(LEN=70), DIMENSION(13) :: soil_desc
    TYPE vegin_type ! dimension will be # of vegetation types
       REAL(r_1), DIMENSION(13) :: canst1
       REAL(r_1), DIMENSION(13) :: dleaf
       REAL(r_1), DIMENSION(13) :: vcmax
       REAL(r_1), DIMENSION(13) :: ejmax
       REAL(r_1), DIMENSION(13) :: hc
       REAL(r_1), DIMENSION(13) :: xfang
       REAL(r_1), DIMENSION(13) :: rp20
       REAL(r_1), DIMENSION(13) :: rpcoef
       REAL(r_1), DIMENSION(13) :: rs20 
       REAL(r_1), DIMENSION(13) :: shelrb
       REAL(r_1), DIMENSION(13) :: frac4
       REAL(r_1), DIMENSION(13) :: tminvj
       REAL(r_1), DIMENSION(13) :: tmaxvj
       REAL(r_1), DIMENSION(13) :: vbeta
       REAL(r_1), DIMENSION(3,13) :: cplant
       REAL(r_1), DIMENSION(2,13) :: csoil
       REAL(r_1), DIMENSION(3,13) :: ratecp
       REAL(r_1), DIMENSION(2,13) :: ratecs
    END TYPE vegin_type
    TYPE soilin_type
       REAL(r_1), DIMENSION(6,13) :: froot
       REAL(r_1), DIMENSION(13) :: silt
       REAL(r_1), DIMENSION(13) :: clay
       REAL(r_1), DIMENSION(13) :: sand
       REAL(r_1), DIMENSION(13) :: swilt
       REAL(r_1), DIMENSION(13) :: sfc
       REAL(r_1), DIMENSION(13) :: ssat
       REAL(r_1), DIMENSION(13) :: bch
       REAL(r_1), DIMENSION(13) :: hyds
       REAL(r_1), DIMENSION(13) :: sucs
       REAL(r_1), DIMENSION(13) :: rhosoil
       REAL(r_1), DIMENSION(13) :: css
       REAL(r_1), DIMENSION(13) :: c3
    END TYPE soilin_type
    TYPE(soilin_type) :: soilin
    TYPE(vegin_type)  :: vegin
    INTEGER(i_d) :: a,b,c,d,e,j ! do loop counter

    ! Allocate default gridpoint ref number variable:
    ALLOCATE(gdpt(mp))

    ! Allocate CABLE's main variables:
    CALL alloc_cbm_var(air, mp)
    CALL alloc_cbm_var(bgc, mp)
    CALL alloc_cbm_var(canopy, mp)
    CALL alloc_cbm_var(met, mp)
    CALL alloc_cbm_var(bal, mp)
    CALL alloc_cbm_var(rad, mp)
    CALL alloc_cbm_var(rough, mp)
    CALL alloc_cbm_var(soil, mp)
    CALL alloc_cbm_var(ssoil, mp)
    CALL alloc_cbm_var(sum_flux, mp)
    CALL alloc_cbm_var(veg, mp)

    WRITE(logn,*) 'Using DEFAULT parameter set based on lat/lon:'

    !================= Read in vegetation type specifications: ============
    OPEN(40,FILE='surface_data/veg_parm.txt', &
         STATUS='old',ACTION='READ',IOSTAT=ioerror)
    IF(ioerror/=0) CALL abort ('Cannot open vegetation type definitions.')
    READ(40,*)
    READ(40,*)
    READ(40,*) nvegt ! read # vegetation types
    READ(40,*)
    READ(40,*)
    DO a = 1,nvegt 
       READ(40,40) veg_desc(a) ! Read description of each vegetation type
    END DO
40  FORMAT(8X,A70)
    READ(40,*) 
    READ(40,*) 
    READ(40,*) vegin%canst1
    READ(40,*) 
    READ(40,*) vegin%dleaf
    READ(40,*) vegin%vcmax
    READ(40,*) vegin%hc
    READ(40,*) vegin%xfang
    READ(40,*) vegin%rp20
    READ(40,*) vegin%rpcoef
    READ(40,*) vegin%rs20
    READ(40,*) vegin%shelrb
    READ(40,*) vegin%frac4
    READ(40,*) vegin%tminvj
    READ(40,*) vegin%tmaxvj
    READ(40,*) vegin%vbeta
    READ(40,*) vegin%cplant(1,:)
    READ(40,*) vegin%cplant(2,:)
    READ(40,*) vegin%cplant(3,:)
    READ(40,*) vegin%csoil(1,:)
    READ(40,*) vegin%csoil(2,:)
    READ(40,*) 
    READ(40,*) vegin%ratecp(:,1)
    ! Set ratecp to be the same for all veg types:
    vegin%ratecp(1,:)=vegin%ratecp(1,1)
    vegin%ratecp(2,:)=vegin%ratecp(2,1)
    vegin%ratecp(3,:)=vegin%ratecp(3,1)
    READ(40,*) 
    READ(40,*) vegin%ratecs(:,1)
    vegin%ratecs(1,:)=vegin%ratecs(1,1)
    vegin%ratecs(2,:)=vegin%ratecs(2,1)
    CLOSE(40)

    ! Description of soil types:
    soil_desc(1)="Coarse sand/Loamy sand"
    soil_desc(2)="Medium clay loam/silty clay loam/silt loam"
    soil_desc(3)="Fine clay"
    soil_desc(4)="Coarse-medium sandy loam/loam"
    soil_desc(5)="Coarse-fine sandy clay"
    soil_desc(6)="Medium-fine silty clay"
    soil_desc(7)="Coarse-medium-fine sandy clay loam"
    soil_desc(8)="Organic peat"
    soil_desc(9)="Permanent ice"

    ! Details from Eva's file soil.txt:
    soilin%froot(1,:) = (/ &
         .02,.04,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05/)
    soilin%froot(2,:) = (/ &
         .06,.11,.11,.11,.11,.15,.15,.10,.10,.10,.10,.15,.15/)
    soilin%froot(3,:) = (/ &
         .14,.20,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35/)
    soilin%froot(4,:) = (/ &
         .28,.26,.26,.26,.26,.39,.39,.35,.35,.35,.35,.38,.40/)
    soilin%froot(5,:) = (/ &
         .35,.24,.24,.24,.24,.05,.05,.10,.10,.10,.10,.06,.04/)
    soilin%froot(6,:) = (/ &
         .15,.15,.15,.15,.15,.01,.01,.05,.05,.05,.05,.02,.01/)
    soilin%silt = (/.08, .33, .17, .2, .06, .25, .15, .70, .33 &
         , .2, .33, .33, .17/)                        ! with mxst=13
    soilin%clay = (/.09, .3, .67, .2, .42, .48, .27, .17, .30 &
         , .2, .3, .3, .67/)                          ! with mxst=13
    soilin%sand = (/.83, .37, .16, .6, .52, .27, .58, .13, .37 &
         , .6, .37, .37, .16/)                        ! with mxst=13
    soilin%swilt = (/ .072, .216, .286, .135, .219, .283, .175, .395 &
         , .216, .1142, .1547, .2864, .2498/)
    soilin%sfc = (/  .143, .301, .367, .218, .31 , .37 , .255, .45 &
         , .301, .22 , .25 , .367, .294/)
    soilin%ssat = (/ .398, .479, .482, .443, .426, .482, .420, .451 &
         , .479, .435, .451, .482, .476/)
    soilin%bch = (/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9 &
         , 5.39, 11.4, 8.52/)      ! bch for gravity term
    soilin%hyds = (/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6, &
         800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/) ! ksat
    soilin%sucs = (/-.106, -.591, -.405, -.348, -.153, -.49, -.299, &
         -.356, -.153, -.218, -.478, -.405, -.63/)   ! phisat (m)
    soilin%rhosoil = (/1600.,1600.,1600.,1600.,1600.,1600.,1600., 1300., &
         910., 1600.,1600.,1600.,1600./) ! soil density
    soilin%css = (/850.,850.,850.,850.,850.,850.,850.,1920., 2100., &
         850.,850.,850.,850./)       ! soil specific heat capacity
    soilin%c3 = (/1.255,.334,.138,.521,.231,.199,.375,.623,.334,.199,.375, &
         .623,.334/)

    ! Set those parameters not dependent on vegetation/soil type in Eva's file:
    soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! layer thickness 20/11/03 
    soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints:
    soil%zshh(ms+1)=0.5*soil%zse(ms)
    soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
    rough%za = 40.0 ! lowest atm. model layer/reference height
    soil%albsoil = 0.1 ! soil albedo
    veg%meth = 1 ! canopy turbulence parameterisation method: 0 or 1

    ! Other initialisations (all gridcells):
    canopy%cansto = 0.0   ! canopy water storage (mm or kg/m2)
    canopy%sghflux = 0.0
    canopy%ghflux = 0.0
    ssoil%ssdnn  = 140.0 ! overall snow density (kg/m3)
    ssoil%snowd  = 0.0   ! snow liquid water equivalent depth (mm or kg/m2)
    ssoil%osnowd = 0.0   ! snow depth prev timestep (mm or kg/m2)
    ssoil%snage  = 0.0   ! snow age
    ssoil%isflag = 0     ! snow layer scheme flag (0 = no or little snow, 1=snow)
    ! snmin  = 0.11  ! switch b/w 1 and 3 layer snow at this depth (m)
    ssoil%wbice  = 0.0   ! soil ice 
    ssoil%tggsn  = 273.1 ! snow temperature per layer (K)
    ssoil%ssdn   = 140.0 ! snow density per layer (kg/m3)
    ssoil%smass  = 0.0   ! snow mass per layer (kg/m^2)
    ssoil%runoff = 0.0   ! runoff total = subsurface + surface runoff
    ssoil%rnof1  = 0.0   ! surface runoff (mm/timestepsize)
    ssoil%rnof2  = 0.0   ! deep drainage (mm/timestepsize)
    ssoil%rtsoil = 100.0 ! turbulent resistance for soil
    canopy%ga     = 0.0   ! ground heat flux (W/m2)
    canopy%dgdtg  = 0.0   ! derivative of ground heat flux wrt soil temp
    canopy%fev    = 0.0   ! latent heat flux from vegetation (W/m2)
    canopy%fes    = 0.0   ! latent heat flux from soil (W/m2)
    canopy%fhs    = 0.0   ! sensible heat flux from soil (W/m2)
    ssoil%albsoilsn(:,1) = 0.1 ! soil+snow albedo
    ssoil%albsoilsn(:,2) = 0.3
    ssoil%albsoilsn(:,3) = 0.05
    ! Initialise sum flux variables:
    sum_flux%sumpn = 0.0
    sum_flux%sumrp = 0.0
    sum_flux%sumrpw = 0.0
    sum_flux%sumrpr = 0.0
    sum_flux%sumrs = 0.0
    sum_flux%sumrd = 0.0
    sum_flux%dsumpn = 0.0
    sum_flux%dsumrp = 0.0
    sum_flux%dsumrd = 0.0
    ! Initialise conservation variables:
    bal%precip_tot = 0.0
    bal%rnoff_tot = 0.0
    bal%evap_tot = 0.0
    bal%wbal_tot = 0.0
    bal%ebal_tot = 0.0
    bal%drybal = 0.0
    bal%wetbal = 0.0

    !============ Decide vegetation types based on lon/lat ==============
    OPEN(41,FILE='surface_data/fort.22', &
         STATUS='old',ACTION='READ',IOSTAT=ioerror)
    IF(ioerror/=0) CALL abort('Cannot open veg/soil type file.')
    gridFlag=.FALSE. ! initialise
    ! Allocations for parameter read in:
    ALLOCATE(lon_temp(numGdpt),lat_temp(numGdpt),land_temp(numGdpt))
    ALLOCATE(isoil(numGdpt),iveg(numGdpt),tgg_temp(numGdpt,2),wb_temp(numGdpt,2))
    READ(41,*) ! skip header lines
    READ(41,*)
    DO e=1,mp ! over all land grid points
       distance = 10.0 ! initialise, units are degrees
       ALLOCATE(oldgdpt(max_it))
       oldgdpt = 0
       ! convert longitude:
       IF(longitude(e)<0.0) THEN
          lon2 = longitude(e) + 360.0
       ELSE
          lon2 = longitude(e)
       END IF
       
       DO d=1,max_it ! max # additional grid searches if site is sea
          DO b=1,numGdpt
             IF(e==1.AND.d==1) THEN ! only read file once:
                READ(41,41) lon_temp(b),lat_temp(b),land_temp(b),isoil(b),iveg(b), &
                     tgg_temp(b,1),tgg_temp(b,2),wb_temp(b,1),wb_temp(b,2) 
             END IF
             ! If current gridcell is closer, set stored to current
             IF(SQRT((latitude(e)-lat_temp(b))**2 +(lon2-lon_temp(b))**2) &
                  < distance.AND..NOT.ANY(oldgdpt==b)) THEN
                ! Reset distance:
                distance = SQRT((latitude(e)-lat_temp(b))**2 +(lon2-lon_temp(b))**2)
                ! Note grid centre:
                grid_lat = lat_temp(b)
                grid_lon = lon_temp(b)
                
                ! Set gridpoint number:
                gdpt(e) = b
                ! Record soil type and veg type:
                v_type  = iveg(b)
                s_type = isoil(b)
                land = land_temp(b) ! land sea flag
                ! Set initial soil temperature and moisture:
                ssoil%tgg(e,1) = tgg_temp(b,1) ! soil temperature, 6 layers (K)
                ssoil%tgg(e,2) = tgg_temp(b,1) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,3) = ssoil%tgg(e,2) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,4) = ssoil%tgg(e,3) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,5) = ssoil%tgg(e,4) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,6) = tgg_temp(b,2)
                ssoil%wb(e,1) = wb_temp(b,1) ! volumetric soil moisture,6 layers
                ssoil%wb(e,2) = wb_temp(b,1) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,3) = ssoil%wb(e,2) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,4) = ssoil%wb(e,3) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,5) = ssoil%wb(e,4) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,6) = wb_temp(b,2)
             END IF
          END DO ! over all default par grid points 

          IF(land/=' T') THEN
             IF(d==1) THEN
                ! If this is the first problem gdpt, write nature of problem to screen:
                IF(.NOT.gridFlag) WRITE(*,'(A39)') 'Problem with default parameter loading:'
                ! Write to screen details about current problem grid point:
                WRITE(*,'(A31,F6.2,A5,F6.2,A29,I7,A4,I7,A1)') '  Nearest default gridpoint to ', &
                     latitude(e),' lat ',lon2, &
                     ' lon is SEA - see log file. (', e, ' of ', mp, ')'
                ! Write again to log file:
                WRITE(logn,2) '     Site is close to gridcell centred:', &
                     grid_lat,'lat', grid_lon,'lon'
                WRITE(logn,'(A63,I7,A4,I7,A1)') &
                     '     this is SEA! ...looking for nearest land point... (land pt', e, ' of ', mp,')'
                gridFlag = .TRUE. ! ie there have been issues finding default pars
             ELSE
                WRITE(logn,*) 'trying ',grid_lat,'lat', grid_lon,'lon'
             END IF
             REWIND(41) ! rewind open file
           !  READ(41,*)
           !  READ(41,*)
             distance = distance + 3.0 ! reset distance
             oldgdpt(d) = gdpt(e)
          ELSE
             ! ie we've found an appropriate default par gridpt for this lat/lon
             EXIT
          END IF
       END DO
       IF(land/=' T') THEN
          WRITE(*,'(A16,I3,A59)') '    The nearest ',max_it, &
               ' default parameter grid locations to this site are all SEA:' 
          WRITE(*,'(A17,I6,A16,F6.2,A12,F6.2)') '    land point # ',e, &
               '      latitude: ',latitude(e),' longitude: ',lon2
          CALL abort('    *** This location is sea - ABORTING *** ')
       END IF
       IF(v_type<0.OR.v_type>13) CALL abort('Unknown veg type!')
       IF(s_type<0.OR.s_type>13) CALL abort('Unknown soil type!')

       ! Report veg type, soil type and gridcell number to log file:
       WRITE(logn,*)
       WRITE(logn,*)
       WRITE(logn,21) ' DETAILS FOR LAND POINT ', e
21     FORMAT(A24,I8)
       WRITE(logn,4) '     Latitude ',latitude(e),'Longitude',longitude(e)
       WRITE(logn,3) '     Site is closest to gridcell #',gdpt(e), &
            ' centred:',grid_lat,'lat', grid_lon,'lon'
       WRITE(logn,*) '    Vegetation type: ',v_type,TRIM(veg_desc(v_type))
       WRITE(logn,*) '    Soil type:       ',s_type,TRIM(soil_desc(s_type))

       ! Prescribe parameters for current gridcell based on veg/soil type:
       veg%iveg(e)   = v_type
       veg%canst1(e) = vegin%canst1(v_type)
       veg%dleaf(e)  = vegin%dleaf(v_type)
       veg%vcmax(e)  = vegin%vcmax(v_type)
       veg%ejmax(e)  = vegin%ejmax(v_type)
       veg%hc(e)     = vegin%hc(v_type)
       veg%xfang(e)  = vegin%xfang(v_type)
       veg%vbeta(e)  = vegin%vbeta(v_type)
       veg%rp20(e)   = vegin%rp20(v_type)
       veg%rpcoef(e) = vegin%rpcoef(v_type)
       soil%rs20(e)   = vegin%rs20(v_type)
       veg%shelrb(e) = vegin%shelrb(v_type)
       veg%frac4(e)  = vegin%frac4(v_type)
       veg%tminvj(e) = vegin%tminvj(v_type)
       veg%tmaxvj(e) = vegin%tmaxvj(v_type)
       bgc%cplant(e,:) = vegin%cplant(:,v_type)
       bgc%csoil(e,:)  = vegin%csoil(:,v_type)
       bgc%ratecp(:) = vegin%ratecp(:,v_type)
       bgc%ratecs(:) = vegin%ratecs(:,v_type)
       soil%froot(e,:)  = soilin%froot(:,s_type)
       soil%silt(e)   =  soilin%silt(s_type)
       soil%clay(e)   =  soilin%clay(s_type)
       soil%sand(e)  =  soilin%sand(s_type)
       soil%swilt(e)  =  soilin%swilt(s_type)
       soil%sfc(e)    =  soilin%sfc(s_type)
       soil%ssat(e)   =  soilin%ssat(s_type)
       soil%bch(e)    =  soilin%bch(s_type)
       soil%hyds(e)   =  soilin%hyds(s_type)
       soil%sucs(e)   =  soilin%sucs(s_type)
       soil%rhosoil(e)=  soilin%rhosoil(s_type)
       soil%css(e)    =  soilin%css(s_type)
       soil%c3(e)     =  soilin%c3(s_type)
       rad%latitude(e) = latitude(e)
       soil%isoilm(e)  = s_type
       veg%ejmax(e) = 2.0*veg%vcmax(e) 
       soil%cnsd(e)  = soil%sand(e)*0.3 + soil%clay(e)*0.25 &
            + soil%silt(e)*0.265 ! set dry soil thermal conductivity [W/m/K]
       soil%hsbh(e)  = soil%hyds(e)*ABS(soil%sucs(e))*soil%bch(e)        !difsat*etasat
       soil%ibp2(e)  = NINT(soil%bch(e))+2
       soil%i2bp3(e) = 2*NINT(soil%bch(e))+3
       ! Conservation variable initialisations:
       bal%wbtot0(e) = 0.0
       DO j=1,ms
          bal%wbtot0(e) = bal%wbtot0(e) + ssoil%wb(e,j) * soil%zse(j) * 1000.0
       END DO
       bal%osnowd0(e) = ssoil%osnowd(e)
       ! Write parameter set details to log file:
       WRITE(logn,*) '============================================================'
       WRITE(logn,*) '    Vegetation parameters: '
       WRITE(logn,11) 'Vegetation height (m): ',veg%hc(e)
       WRITE(logn,11) 'Fraction of plants which are C4 (-): ',veg%frac4(e)
       WRITE(logn,11) 'Maximum canopy water storage (mm/LAI): ',veg%canst1(e)
       WRITE(logn,15) 'Max pot elec transport rate top leaf (mol/m2/s): ', &
            veg%ejmax(e)
       WRITE(logn,15) 'Max RuBP carboxylation rate top leaf (mol/m^2/s): ', &
            veg%vcmax(e)
       WRITE(logn,11) 'Plant respiration coeff @ 20 C (mol/m^2/s): ', &
            veg%rp20(e)
       WRITE(logn,11) 'Temperature coef nonleaf plant respiration (1/C): ', &
            veg%rpcoef(e)
       WRITE(logn,11) 'Sheltering factor (-): ', veg%shelrb(e)
       WRITE(logn,11) 'Chararacteristic legnth of leaf (m): ', veg%dleaf(e)
       WRITE(logn,11) 'Leaf angle parameter (-): ', veg%xfang(e)
       WRITE(logn,11) 'Min temperature for start of photosynthesis (C): ', &
            veg%tminvj(e)
       WRITE(logn,11) 'Max temperature for start of photosynthesis (C): ', &
            veg%tmaxvj(e)
       WRITE(logn,11) ' vbeta: ', &
            veg%vbeta(e)
       WRITE(logn,14) ' Plant carbon rate constant (1/year): ',bgc%ratecp
       WRITE(logn,*) '============================================================'
       WRITE(logn,*) '    Soil parameters: '
       WRITE(logn,11) 'Fraction of soil which is sand (-): ',soil%sand(e)
       WRITE(logn,11) 'Fraction of soil which is silt (-): ',soil%silt(e)
       WRITE(logn,11) 'Fraction of soil which is clay (-): ',soil%clay(e)
       WRITE(logn,13) 'Fraction of roots in each soil layer (-): ', &
            soil%froot(e,:)
       WRITE(logn,11) 'Volumetric soil moisture at saturation (m^3/m^3): ', &
            soil%ssat(e)
       WRITE(logn,11) 'Volumetric soil moisture at field capacity (m^3/m^3): ', &
            soil%sfc(e)
       WRITE(logn,11) 'Volumetric soil moisture at wilting point (m^3/m^3): ', &
            soil%swilt(e)
       WRITE(logn,11) 'Soil respiration coeff @ 20C (mol/m^2/s): ', &
            soil%rs20(e)
       WRITE(logn,11) 'Suction at saturation (m): ', soil%sucs(e)
       WRITE(logn,11) 'Soil density (kg/m^3): ',soil%rhosoil(e)
       WRITE(logn,11) 'Soil specific heat capacity (kJ/kg/K): ', soil%css(e)
       WRITE(logn,11) 'parameter b in Campbell equation: ', soil%bch(e)
       WRITE(logn,'(4X,A50,E10.4)') 'hydraulic conductivity @ saturation (m/s): ',soil%hyds(e)
       WRITE(logn,11) 'c3 drainage coeff (-): ',soil%c3(e)
       WRITE(logn,14) 'Soil carbon rate constant (1/year): ',bgc%ratecs
       WRITE(logn,11) 'Bare soil albedo (-): ', soil%albsoil(e)
       WRITE(logn,*) '============================================================'
       WRITE(logn,*) '    Site dependent initialisations: '
       WRITE(logn,14) ' Plant carbon pool size (g C/m2)): ', bgc%cplant(e,:)
       WRITE(logn,14) ' Soil carbon pool size (g C/m2)): ', bgc%csoil(e,:)
       WRITE(logn,13) 'Soil moisture, by layer: ', ssoil%wb(e,:)
       IF(ANY(ssoil%wb(e,:)<soil%swilt(e))) &
            WRITE(logn,*) '    SOIL MOISTURE INITIALISED BELOW WILTING POINT!'
       WRITE(logn,13) 'Soil temperature, by layer (K): ', ssoil%tgg(e,:)
       WRITE(logn,14) 'Soil+snow albedo (-): ', ssoil%albsoilsn(e,:)
       WRITE(logn,*) '============================================================'
       WRITE(logn,*) '    Site independent initialisations: '
       WRITE(logn,11) 'Canopy surface water storage (mm): ',canopy%cansto(e)
       WRITE(logn,11) 'Snow liq water equiv depth (mm): ', ssoil%snowd(e)
       WRITE(logn,11) 'Overall snow density (kg/m^3): ',ssoil%ssdnn(e)
       WRITE(logn,11) 'Snow age (-): ',ssoil%snage(e)
       WRITE(logn,13) 'Volumetric soil ice, by layer (?): ', &
            ssoil%wbice(e,:)
       WRITE(logn,13) 'Snow temperature, by layer (K): ', &
            ssoil%tggsn(e,:)
       WRITE(logn,13) 'Snow density, by layer (kg/m^3): ', &
            ssoil%ssdn(e,:)
       WRITE(logn,13) 'Snow mass, by layer (kg/m^2): ', &
            ssoil%smass(e,:)
       WRITE(logn,11) 'Runoff total (mm/time step): ', ssoil%runoff(e)
       WRITE(logn,11) 'Surface runoff (mm/time step): ', ssoil%rnof1(e)
       WRITE(logn,11) 'Deep drainage runoff (mm/time step): ', ssoil%rnof2(e)
       WRITE(logn,*) '============================================================'
    END DO
    CLOSE(41)
41  FORMAT(14X,2F8.3,A2,14X,2I3,14X,2F7.1,2F6.2)
11  FORMAT(4X,A50,F10.4)
15  FORMAT(4X,A50,E10.4)
12  FORMAT(4X,A50,I8)
13  FORMAT(4X,A50,6F10.4)
14  FORMAT(4X,A50,3F11.4)
2   FORMAT(A42,2(1X,F6.2,1X,A3))
3   FORMAT(A34,I7,A9,2(1X,F8.3,1X,A3))
4   FORMAT(2(A14,1X,F9.4,1X))

    DEALLOCATE(lon_temp,lat_temp,land_temp,isoil,iveg,tgg_temp,wb_temp,oldgdpt) 

  END SUBROUTINE default_params
   
END MODULE parameter_module


