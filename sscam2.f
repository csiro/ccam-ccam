      SUBROUTINE SSCAM (timi,npij,kpoff,
     &                  fsd,fld,ta,qa,ca,ua,pmb,precip,
     &                  init,zeta0,t0v0,z0m)
c-----------------------------------------------------------------------
c Snow-Soil-Canopy-Atmosphere Model for energy and CO2 exchange between a
c vegetated land surface and the atmosphere
c Created by merging SCAM22b and the snow-soil-ice-temperature model by
c Ewa Kowalczyk
c-----------------------------------------------------------------------
c SSCAM1 (KF, 24-8-97)
c   (a) streamlining scam22b to accomodate soil water and temperature
c       soiltypes of global model, rootdistribution, albedo, similar
c       lai calculation with dependence on deep soil temperature
c       standalone version is not perseved needs work to reinstate
c       common blocks rationalised, hold now 2dim array holding patch 
c       information (lai varing due to snow and seasonality based on TS)
c sscam2
c   (a) cansat = 0 when ta<0, no interception when snow
c   (b) store now not used anymore, eta and ts direcly accesses darlam arrays
c   (c) soilc replaced by insoil (Ewa's)
c
c-----------------------------------------------------------------------
c INPUTS:
c   timi   ! time steps from start of simulation (for FC infiltration)
c   npij   ! number of patches for current gridcell (1=single column)
c   kpoff  ! offset of current gridcell in patch array (0=single column)
c   fsd    ! downward global solar irradiance              [W/m2]
c   fld    ! downward global longwave irradiance           [W/m2]
c   ta     ! air temp at height za                         [deg C]
c   qa     ! specific humidity at za                       [kg/kg]
c   ca     ! CO2 concentration at za                       [mumol/mol]
c   ua     ! wind speed at height za                       [m/s]
c   pmb    ! atmospheric pressure                          [mb]
c   precip ! precipitation during step                     [mm]
c   init   ! initialisation flag: init<0 for run, init>0 for section
c   param0 ! global (patch-independent) parameters
c   param1 ! default local (patch-dependent) parameters
c   param2 ! patch-dependent local params for 2D version
c   zeta0  ! start-of-step Monin-Obukhov parameter za/L    [-]
c   t0v0   ! start-of-step veg surface temperature         [C]
c OUTPUTS:
c   hold   ! array of parameters computed on first call and saved
c   hname  ! array of 8-char names assigned to variables in hold
c   outvar ! array of output variables returned at each time step
c   outij  ! array of patch-averaged outputs for current gridcell
c   Elements of hold:
c         1:  rlai   ! total LAI                           [-]
c         2:  hc     ! canopy height                       [m]
c         3:  disp
c         8:  rt0us  ! us*rt0 (z=0 to z=disp)              [-]
c         9:  rt1usa ! us*rt1a (z=disp to z=hc)            [-]
c   Elements of outvar:
c         1:  fn     ! ACTUAL net irradiance (+ down)      [W/m2]
c         2:  gfluxkf! soil heat flux (+ down)             [W/m2]
c         3:  fh     ! sensible heat flux                  [W/m2]
c         4:  fe     ! latent heat flux                    [W/m2]
c         5:  fc     ! mass CO2 flux                       [kgCO2/m2/s]
c         6:  us2    ! (friction velocity)**2              [m2/s2]
c         7:  zetar  ! Monin-Obukhov parameter za/L        [-]
c         8:  rbv    ! veg-bulked bound-layer resistance   [s/m]
c         9:  rsv    ! veg-bulked stomatal resistance      [s/m]
c        10:  rt0    ! turb resistance: (0,d)              [s/m]
c        11:  rt1    ! turb resistance: (d,za)             [s/m]
c        12:  runoff ! runoff in this step                 [mm]
c        13:  rinfil ! surface infiltration this step      [mm]
c        14:  rexfil ! exfiltration (deep drainage)        [mm]
c        15:  evap   ! total evaporation this step         [mm]
c        16:  fnv,fhv,fev    ! veg  parts of fn,fh,fe      [W/m2]
c        19:  fns,fhs,fes    ! soil parts of fn,fh,fe      [W/m2]
c        22:  fesatm,fessoi  ! atm-lim, soil-lim soil fe   [W/m2]
c        24:  fevc,fevw      ! dry-veg, wet-veg fe         [W/m2]
c        26:  fsu,flu        ! total up short, long rads   [W/m2]
c        28:  t0v,q0v,d0v    ! veg surface t,q,deficit     [C, kg/kg]
c        31:  tv, qv, dv     ! veg air layer t,q,deficit   [C, kg/kg]
c        34:  zetar(10)      ! iterations toward zeta      [-]
c        44:  rnofdb         ! dunne-black runoff in step  [mm]
c-----------------------------------------------------------------------
c Begin control statements
      include 'newmpar.h'   ! parameter statement darlam npatch
      include 'parm.h'      ! ktau,idjd for diags
      include 'scamdim.h'   ! parameter statement for all dimensions
      include 'scampar.h'   ! dimension statement for scam part
      include 'soilsnow.h'  ! tgg,wb soil temperature and moisture
      include 'soilv.h'
c     common/soilpr/
c    & rhos(9),css(9),cnsd(9),bch(9),ssat(9),swilt(9),hsbh(9),i2bp3(9),
c    & ibp2(9),sfc(9),sucs(9),hyds(9),silt(9),clay(9),sand(9),sormax(9),
c    & tgmax(9),scond1s(3),cap1s(3),ssatcur(ms+1),coef(ms+4),cnsw(ms),
c    & gamm(ms),zse(ms),zsh(npmax,ms+1)

c     vegetation properties set for each veg type at initialisation
c     (ivtype=44) is ST version, (ivtype=0,31) is 2D version:
      common /canopy/
     &    rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax(0:44),xjmax0(0:44)
c     link to subroutines needing constant soil properties:
c     link to subroutine swsto*
      common /swblk/ istype,zsw1,zse1,zsm,qes,qevc,qrunf,qinf,qexf
c     variables needed for infiltration calculation
      common /swblkb/ ksw,tsrt(npmax),tend(npmax),
     &     sorin(npmax),suminf(npmax)
c     root distribution passed to soilsnow
c     real froot(5)
      save nprint
      data nprint/0/

c Physical and mathematical constants:
      parameter(sboltz = 5.67e-8,  ! Stefan-Boltzmann constant (W/m2/K4)
     &          rgas   = 8.3143,   ! universal gas constant    (J/mol/K)
     &          rmair  = 0.02897,  ! molecular weights: dry air (kg/mol)
     &          rmh2o  = 0.018016, !                    WATER   (kg/mol)
     &          rmco2  = 0.044022, !                    CO2     (kg/mol)
     &          grav   = 9.806,    ! gravity acceleration (m/s2)
     &          pi     = 3.141592654)
c Aerodynamic parameters, diffusivities, water density:
      parameter(vonk   = 0.40,     ! von Karman constant
     &          a33    = 1.25,     ! inertial sublayer sw/us
     &          csw    = 0.50,     ! canopy sw decay (Weil theory)
     &          ctl    = 0.40,     ! Wagga wheat (RDD 1992, Challenges)
     &          apol   = 0.70,     ! Polhausen coeff: single-sided plate
     &          prandt = 0.71,     ! Prandtl number: visc/diffh
     &          schmid = 0.60,     ! Schmidt number: visc/diffw
     &          diffwc = 1.60,     ! diffw/diffc = H2O/CO2 diffusivity
     &          rhow   = 1000.0)   ! liquid water density   [kg/m3] 
c-----------------------------------------------------------------------
c Begin executable code

c Set GLOBAL input parameters, independent of gridcell or patch, passed 
c via PARAM0(20):
c   ref height, time step, switches:
      za     = param0(1)     ! ref height (z=0 at disp)    [m]
      dels   = param0(2)     ! time step                   [s]
c     istemp = param0(3)     ! soil temp:   1,2,3,4 = FR,kf,mrr,mrrkf
c     ismois = param0(4)     ! soil moist:  1,2,3   = MP84,NP89,Richards
c     isinf  = param0(5)     ! soil infilt: 1,2     = MP84, FC96
c     isevap = param0(6)     ! soil evap: 1,2,3 = alfa,beta,threshold
c     itherm = param0(7)     ! VW or KGK algorithm for hconds,rkapps
c     irktem = param0(8)     ! RK steps in soil temp schemes
c     irkmoi = param0(9)     ! RK steps in soil moisture schemes
c not supported anymore     istory = param0(10)    ! under/overstory are used
c   turbulence parameters:
c     niter  = param0(11)    ! number of iterations for za/L
c     zetmul = param0(12)    ! if niter=2, final zeta=zetmul*zetar(2)
c     zetneg = param0(13)    ! negative limit on za/L when niter>=3
c     zetpos = param0(14)    ! positive limit on za/L when niter>=3
c     umult  = param0(15)    ! wind speed multiplier
c     zdlin  = param0(16)    ! height frac of d below which TL linear
c     umin   = param0(17)    ! minimum wind speed          [m/s]
c   soil water parameters:
c     etarct = param0(18)    ! rel soil moisture for finding zst1,zst2
c     dbde   = param0(19)    ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
c     fdrain = param0(20)    ! exfilt param: qexf=fdrain*cond2
c     istom  = param0(21)    ! stomata: 1,2,3 = Jarvis, AssimF, AssimH
c     istsw  = param0(22)    !

      
c Set DEFAULT LOCAL input parameters, passed via param1(30):
c NB: (2D) parameters are overwritten by patch-dependent veg type or
c          soil type assignments from ivtype in 2D version
c     (ST) parameters are overwritten by section assignments in
c          standalone version
c   geometric, radiative and aerodynamic properties:
      hc     = param1(1)     ! canopy height       (2D,ST) [m]
      rlai   = param1(2)     ! leaf area index     (2D,ST) [-]
      canstl = param1(3)     ! canopy water store / LAI    [mm]
      z0soil = param1(4)     ! soil z0 (used when rlai=0)  [m]
      albedv = param1(5)     ! vegetation albedo   (2D)    [-]
      albeds = param1(6)     ! soil albedo                 [-]
      emissv = param1(7)     ! vegetation emissivity       [-]
      emisss = param1(8)     ! soil emissivity             [-]
      cextin = param1(9)     ! radiation extinction coeff  [-]
      dleaf  = param1(10)    ! leaf dimension              [m]
      shelrb = param1(11)    ! rbv multiplier for shelter  [-]
c   biological properties:
      resp20 = param1(14)    ! soil CO2 flux at 20 C       [kgCO2/m2/s]
      respdt = param1(15)    ! soil CO2 flux doubling temp [C]
      if (abs(resp20).gt.1.0) resp20 = 0.0                 ! cut 999
      gsmax  = param1(16)    ! maximum (unconstrained) gs  [m/s]
      fcfsd  = param1(17)    ! gs: solar constraint param  [W/m2]
      fcdef  = param1(18)    ! gs: satdef constraint param [kg/kg]
      fceta1 = param1(19)    ! gs: max eta/ssat          [-]
      fceta2 = param1(20)    ! gs: min eta/ssat          [-]
      fctemh = param1(21)    ! gs: max temp                [C]
      fctemL = param1(22)    ! gs: min temp                [C]
      fctemm = (fctemL+fctemh)/2.0 ! gs: midpoint temp     [C]
      cint   = param1(23)    ! intercellular [CO2]         [mol/mol]
      rjmax0 = param1(32)    ! potl elect transport at T0  [mol/m2/s]
      alfaqy = param1(33)    ! quantum yield               [-]
      g0     = param1(34)    ! intercept param in gs model [mol/m2/s]
      a1     = param1(35)    ! slope param in gs model     [-]


c   soil properties: assign in patch loop
c INITIALISATIONS: initialise run if init<0, section if init>0
      if (init.ne.0) then
        call isscam (init,npij,kpoff,zeta0,z0m)
      end if
c INPUT FOR soilscam for grid average input
      theta = ta + 273.16
      t1    = ta + 273.16

c LOOP through patches
      DO 10 kp=1,npij
        kij = kp + kpoff               ! offset of current gridcell 

c Assign ivtype, istype, zsm:
        ivtype = param2(kij,1)
c        if (ivtype.eq.44) then         ! standalone version
c          zsm    = param1(24)          ! soil depth [m]
c          istype = param1(25)          ! soil type
c        else                           ! 2D version
          istype = param2(kij,2)    ! depth-uniform soil type
          zsm    = param2(kij,3)
c        end if
c Assign rooting fraction within layer based on rooting depth zsm       
c uniform distribution
c       do i=1,nx
c          bottom   = min(zs(i),zsm)
c          froot(i) = max(0,(bottom-zs(i-1))/zsm)       
c       enddo
c exponential distribution where rho(zsm) = 0.1 => rho=exp(-2/zsm*z)
c and the remainder of the roots is in the layer below, in case
c where zsm is in last layer the remainder of the roots in last layer
c        sum    = 0.      
c        do i=1,nx
c           froot(i) = exp(-3./zsm*zs(i-1)) - exp(-3./zsm*zs(i))
c           if (zsm.gt.zs(i-1).and.zsm.le.zs(i)) iroot=i ! layer# of zsm
c        enddo    
c        iroot = min(iroot,nx-1)      ! in case zsm is in last layer
c        do i=1,iroot
c           sum   = sum + froot(i)
c        enddo
c        froot(iroot+1)= max(1.-sum,0)! rest of the roots in next layer
c        do i=iroot+2,nx
c           froot(i) = 0.             ! no roots below
c        enddo  
c        if (init.ne.0) print*,'95% froot =',(froot(i),i=1,nx)

ckf        sum    = 0.      
ckf        do i=1,4
ckf           froot(i) = exp(-3./zsm*zs(i-1)) - exp(-3./zsm*zs(i))
ckf           if (zsm.gt.zs(i-1).and.zsm.le.zs(i)) iroot=i ! layer# of zsm
ckf        enddo  
ckf        irootmax = min(iroot-1,3)  
ckf        do i=1,irootmax
ckf           sum   = sum + froot(i)
ckf        enddo
ckf        froot(irootmax+1)= max(1.-sum,0.0)  ! rest of roots in layer where zsm is
ckf        do i=irootmax+2,4 
ckf           froot(i) = 0.             ! no roots below
ckf       enddo
      froot(1) = 0.20
      froot(2) = 0.45
      froot(3) = 0.20
      froot(4) = 0.10
      froot(5) = 0.05
 
c      if (init.ne.0) write(*,160)  (froot(i),i=1,4)
c160       format('  froot:',6f10.4)
          
c Assign canopy properties:  variables for current patch
        gsmax  = xgsmax(ivtype)
        rjmax0 = xjmax0(ivtype)

        rlai   = hold(kij,1)           ! xrlai(ivtype)
        hc     = hold(kij,2)           ! variables for current patch
        disp   = hold(kij,3)           ! xdisp(ivtype)
        usuh   = hold(kij,4)           ! xusuh(ivtype)
        coexp  = hold(kij,5)           ! xcoexp(ivtype)
        zruffs = hold(kij,6)           ! xzrufs(ivtype)
        trans  = hold(kij,7)           ! xtrans(ivtype)
        rt0us  = hold(kij,8)           ! xrt0us(ivtype)
        rt1usa = hold(kij,9)           ! xt1usa(ivtype)
        rt1usb = hold(kij,10)          ! xt1usb(ivtype)
        cansat = canstl * rlai 

c ******************************************************Darlam test
        albedv = param1(12)
        albeds = param1(12)
c ******************************************************

c Assign input-then-modified parameters from array store:
        ts1    = tgg(kij,1)-273.16  ! layer 1 soil temp (surface) [C]
        eta1   = wb(kij,1)          ! layer 1 soil moisture       [-]
        eta2   = wb(kij,2)          ! layer 2 soil moisture       [-]
        cansto = store(kij,1)     ! canopy water store          [mm]
c        totsto = store(kij,2)     ! total water store           [mm]
c        totsti = totsto            ! hold totsto at start of this step 
         
c Assign special variables for common blocks:
        ksw    = kij               ! for common /swblkb/

c Calculate relative soil moisture
        denom = ssat(istype)-swilt(istype)
        eta1r = (eta1-swilt(istype))/denom
        eta2r = (eta2-swilt(istype))/denom
        eta1r = min(max(eta1r,0.0),1.0)
        eta2r = min(max(eta2r,0.0),1.0)

c MAIN COMPUTATIONS:                                                   
        da   = max(qsatf(ta,pmb)-qa,0.0)  
                                   ! air sat def [kg/kg]
c       t0v0=.5*(ts1+ta)
        t0v0=ts1
        t0v  = t0v0
        q0v  = qa                  ! first ests of surface concs
        d0v  = da
        c0v  = ca
        tv   = ta                  ! first ests of veg layer air concs
        qv   = qa
        dv   = da
        cv   = ca             
        zetar(1) = zeta0           ! initialise loop with input zeta0
        do iter=2,niter            ! reset zetar(iter)
          zetar(iter) = zetpos+1   ! = za/L used on iteration iter
        end do

c changed with new multilayer soil model 
c        if (ismois.le.2) then 
c           etamr = (zsw1*eta1r+(zsm-zsw1)*eta2r)/zsm 
c        elseif (ismois.eq.3) then
           etamr = 0.0
           rwater = 0.0 ! water available to roots
           twater = 0.0 ! maximal available water to roots
            do i=1,4
              rwater = rwater+ froot(i)*zse(i)*(wb(kij,i)-swilt(istype))
              twater = twater+ froot(i)*zse(i)*
     &                 (ssat(istype)-swilt(istype)) 
           enddo                                          
           etamr = rwater/twater  
c        endif

c        etamr=(etamr-(1-fceta2)*swilt(istype))
c     &         /(ssat(istype)-(1-fceta2)*swilt(istype))
        
c LOOP FOR AIR FLUXES and zetar=za/L: starting from previous values,
c iterate for t0v,q0v,d0v,c0v,tv,qv,dv,cv,zetar
        do 1 iter=1,niter
          t0v=max(t0v,min(ts1,ta)-3.)  ! fix to keep t0v OK for radiation
          t0v=min(t0v,max(ts1,ta)+3.)  ! fix to keep t0v OK for radiation
c AIR PROPERTIES rho, capp, rlam, epsi=(rlam/capp)*dqsdt, at midpoint
c temp tm = (t0m+ta)/2
          t0m = (1.0-trans)*t0v + trans*ts1
          tm  = (t0m + ta)*0.5
          call air(tm,pmb,rho,volm,capp,rlam,qsatm,epsi,visc)
c VEGETATION RADIATION FLUXES (factor 2 in fluv for up and down):
          fsdv = (1.0-trans)*fsd
          fsuv = albedv*fsdv
          fldv = (1.0-trans)*fld 
     &           + (1.0-trans)*emisss*sboltz*(ts1+273.16)**4
          fluv = 2.0*(1.0-trans)*emissv*sboltz*(t0v+273.16)**4
          fnv  = fsdv-fsuv+fldv-fluv   ! veg net irradiance
c SOIL RADIATION FLUXES:
          fsds = trans*fsd
          fsus = albeds*fsds
          flds = trans*fld 
     &           +(1.0-trans)*emissv*sboltz*(t0v+273.16)**4
          flus = emisss*sboltz*(ts1+273.16)**4
          fns  = fsds-fsus+flds-flus
c AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent 
c resistances rt0, rt1 (elements of dispersion matrix):
          uam    = max(ua,umin)        ! lower limit for wind speed
          us     = vonk*uam / ( alog(za/z0m) - psim(zetar(iter)) 
     &                          + psim(zetar(iter)*z0m/za) )
          us     = max(us,0.0001)      ! ensure no very small us
          term1  = alog(za / max(zruffs-disp,z0soil) )
          zetruf = zetar(iter)*(zruffs-disp)/za
c x=1 if za+disp>zruffs, 0 otherwise: thus rt1usc = 0 if za+disp<zruffs
          x      = 0.5 + sign(0.5,za+disp-zruffs)
          rt1usc = x*(term1 - psis(zetar(iter)) + psis(zetruf)) / vonk 
          rt1us  = rt1usa + rt1usb + rt1usc
          rt0    = rt0us/us            ! rt from z0 to z1=disp
          rt1    = rt1us/us            ! rt from z1=disp to za
      
c Vegetation boundary-layer resistance rbv:
          diffh  = visc/prandt
          uleaf  = us / usuh                     ! for leaf at z=hc
          releaf = uleaf*dleaf/visc
          gb     = (apol*diffh/dleaf) * releaf**0.5 * prandt**0.33333333
          gb     = gb / shelrb                   ! inter-element shelter
          gbv    = gb * (2.0*rlai/coexp)
          rbv    = 1.0/max(gbv,1e-6)
      
c SURFACE RESISTANCE at current t0v,d0v,c0v [calculate faceta (water
c stress factor between 0 and 1) here because it is needed elsewhere]:
          faceta = (etamr - fceta2) / (fceta1 - fceta2)
          faceta = etamr / fceta1
          faceta = max(faceta,0.0)     ! water constraint
          faceta = min(faceta,1.0)
          if (istom.le.1) then         ! Jarvis scheme
            call gsjarv (fsd,t0v,d0v,c0v,gsmax,cextin,trans,fcfsd,
     &               fcdef,fctemh,fctemm,faceta,cint,volm,gbv,gsv,fcv)
          else                         ! assimilation scheme
            parad  = 2e-6*max(fsd,1.0) ! full-sun PAR [mol/m2/s]
            call PHINT1 (rlai,cextin,parad,t0v,d0v,c0v,faceta,fcdef,
     &                   rjmax0,alfaqy,g0,a1,volm,istom-2,
     &                   an,anc,anq,rd,ci,gsc,gsv,fcv)
          end if
          rsv    = 1.0/max(gsv,1e-6) ! canopy rs

c CARBON FLUXES: vegetative flux fcv already found
c soil-respiration CO2 flux fcs [kgCO2/m2/s]:
          fcs = resp20 * faceta * (2.0**((ts1-20.0)/respdt))
c total CO2 flux [kgCO2/m2/s]:
          fc  = fcv + fcs

c deficit resistances and wet canopy fraction:             
          gdv    = gsv*gbv/((epsi+1.0)*gsv + gbv + 1e-12)
                                       ! veg deficit conductance
          rdv    = (epsi+1.0)*rbv+rsv  ! veg deficit resistance
          gdvw   = gbv/(epsi+1.0)      ! wet veg deficit conductance
          rdvw   = (epsi+1.0)*rbv      ! wet veg deficit resistance
          fwet   = cansto/max(cansat,1e-6)
                                       ! fraction of wet canopy
c canopy temp, humidity tv, qv: 
c   first find feva,fhva using ref level ta,qa, not tv,qv
          term1  = (epsi*rbv*fnv + rho*rlam*da)
          fevca  = (1.0-fwet) * term1*gdv        ! dry canopy fev
          fevwa  = fwet * term1*gdvw             ! wet canopy fev
          feva   = fevca+fevwa                   ! total fev
          fhva   = fnv - feva                    ! total fhv
c   likewise find fesa,fhsa using ref level ta,qa, not tv,qv
c   rtsoil = rt0+rt1 (bare), rt0 (veg)
          rtsoil = rt0 + rt1*(0.5+sign(0.5,0.001-rlai))

c   soil-air heat flux
          fhsa   = rho*capp*(ts1 - ta)/rtsoil
c   soil evaporation flux
          call sevap (qa,ts1,eta1r,rho,rlam,rtsoil,pmb,isevap,dbde,
     &                fesa,qesa,beta)
        
c   find coefficients of linear equations for tv, qv, then solve
          aaa    = rt0*rt1*gbv*gsv
          bbb    = (1.0 + epsi)*gsv + gbv
          aah    = (rt0 + rt1)*bbb + epsi*aaa
          bbh    = -(rlam/capp)*aaa
          cch    = bbb*rt0*rt1*(fhsa+fhva)/(rho*capp)
          aae    = -epsi*(capp/rlam)*aaa
          bbe    = (rt0 + beta*rt1)*bbb + aaa
          cce    = bbb*rt0*rt1*(fesa+feva)/(rho*rlam)
          denom  = aah*bbe-aae*bbh
          tv     = ta + (bbe*cch - bbh*cce) / (denom + 1e-12)
          qv     = qa + (aah*cce - aae*cch) / (denom + 1e-12)
          dv     = max(qsatf(tv,pmb) - qv,0.0)  
          cv     = ca

c VEGETATION SENSIBLE AND LATENT HEAT FLUXES fev, fhv (W/m2), with fev
c partitioned to wet, dry canopies using area fraction cansto/cansat:
          term1  = (epsi*rbv*fnv + rho*rlam*dv)
          fevc   = (1.0-fwet) * term1*gdv        ! dry canopy fev
          fevw   = fwet * term1*gdvw             ! wet canopy fev
          fevwmx = 0.001*cansat*rhow*rlam/dels
          fevw   = min(fevw,fevwmx)              ! limit to max wet fev
          fev    = fevc+fevw                     ! total fev
          fhv    = fnv - fev                     ! total fhv
c mean aerodynamic vegetation surface temperature, humidity, sat def 
c (for next step):
          t0v    = tv + fhv*rbv/(rho*capp)
          q0v    = qv + fev*rbv/(rho*rlam)
          d0v    = max(qsatf(t0v,pmb) - q0v,0.0)
          c0v    = ca
c         c0v    = cv + fcv*rbv/(rmco2/volm)

           if ((t0v.lt.ta-10..and.t0v.lt.ts1-10..and.nprint.lt.1000).or.
     .       (kij.eq.idjd.and.ktau.gt.300.and.ktau.lt.310)) then
c    .       (kij.eq.idjd                                )) then
              nprint=nprint+1
              print *
              print *,'small t0v'
              print *,'ktau,iter,iq: ',ktau,iter,kpoff+npij
              print*,'t0v,ta,ts1,tm      ',t0v,ta,ts1,tm
              print*,'tv,rbv,rsv,denom: ',tv,rbv,rsv,denom
              print*,'zetar: ',(zetar(nn),nn=1,iter)
              print*,'istype,ivtype,za,dels: ',istype,ivtype,za,dels
              print*,'bbe,cch,bbh,cce: ',bbe,cch,bbh,cce
              print*,'rt0,rt1,fhsa,fhva: ',rt0,rt1,fhsa,fhva
              print*,'rt0us,rt1us,bbb,us: ',rt0us,rt1us,bbb,us
              print*,'q0v,qa,eta1: ',q0v,qa,eta1
              print*,'fnv,fhv,fev: ',fnv,fhv,fev
              print*,'fsds,fsus,flds,flus',fsds,fsus,flds,flus
              print*,'hold: ',(hold(kij,ll),ll=1,10)
           endif

c   soil-air heat flux
          fhs   = rho*capp*(ts1 - tv)/rtsoil
c   soil evaporation flux
          call sevap (qv,ts1,eta1r,rho,rlam,rtsoil,pmb,isevap,dbde,
     &                fes,qes,beta)
c TOTAL ENERGY FLUXES:
          fe   = fev+fes
          fh   = fhv+fhs
          fn   = fnv+fns
          flu  = fluv/2.0 + trans*flus
          fsu  = fsuv + fsus                     ! NB: no trans in fsus
          evap = 1000.0*dels*fe/(rhow*rlam)      ! evap in step (mm)
          
c Monin-Obukhov stability parameter zetar=za/L
c   recompute zetar for the next iteration, except on last iteration
c          if (iter.eq.niter) goto 2
          if (iter.lt.niter) then ! dont compute zetar on the last iter
             iterplus = max(iter+1,niter)
             zetar(iterplus) =-(vonk*grav*za*(fh+0.07*fe))
     &                     /(rho*capp*(273.16+ta)*us**3)
c   case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
             if (niter.eq.2) zetar(2) = zetmul*zetar(2)
c   constrain zeta to zetpos and zetneg (set in param0)
             zetar(iterplus) = min(zetpos,zetar(iterplus))    ! zetar too +
             zetar(iterplus) = max(zetneg,zetar(iterplus))    ! zetar too -
           endif
1       continue                       ! loop termination
2       continue                       ! loop exit
c FINISHED LOOP FOR AIR FLUXES and zetar=za/L

c modify CANOPY WATER STORAGE:
        qevc   = fevc/(rhow*rlam)      ! water flux from dry veg (m/s)
        qevw   = fevw/(rhow*rlam)      ! water flux from wet veg (m/s)
        evapvw = qevw*dels*1000.0      ! wet veg evap in this step (mm)
        cansto = cansto - evapvw       ! withdraw evapvw from cansto
        cansto = max(cansto,0.0)       ! do not overdraw
        cansto = cansto + precip       ! add precip to cansto
        precis = 0.0                   ! precip reaching soil [mm]
        water  = 0.5 + sign(0.5,ta)    ! if ta< 0, no water in canopy store,
                                       ! all freezes and falls as snow
        cansat = cansat*water
        full   = 0.5 + sign(0.5,cansto-cansat)
                                       ! full=0,1 for cansto <,> cansat
        precis = full*(cansto-cansat)  ! precip to soil [mm]
        cansto = cansto-precis

c GROUND HEAT FLUXES (at current ts1)
        gfluxkf = fns-fhs-fes
        if ((abs(gfluxkf).gt.500.and.ktau.gt.20).or.
     .       (kij.eq.idjd.and.ktau.gt.300.and.ktau.lt.310)) then
c    .       (kij.eq.idjd                                )) then
          print*,'in sscam; gfluxkf,fns,fhs,fes: ',gfluxkf,fns,fhs,fes
          print *,'iq,init,istype,fh: ',kij,init,istype,fh
          print*,'rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,
     .    rt1usb: ',(hold(kpoff+npij,ll),ll=1,10)
          print*,'ta,tv,ts1,rho,capp,rtsoil: ',ta,tv,ts1,rho,capp,rtsoil
          print*,'tgg: ',(tgg(kij,nlev),nlev=1,ms)
          print*,'fsds,fsus,flds,flus: ',fsds,fsus,flds,flus
          print*,'trans,fld,trans,emissv,sboltz,t0v: '
     &     ,trans,fld,trans,emissv,sboltz,t0v
          print*,'qa,qv,eta1,eta1r: ',qa,qv,eta1,eta1r
        endif
c Find SOIL WATER FLUXES: 
C campbell soil model uses (eta/ssat) and not etarel=(eta-swilt)/(ssat-etadry)
c Infiltration and deep drainage parameterisations for bulk soil models        
        if (isinf.eq.2) then           ! event based infiltration
           etain = 0.5*(eta1+eta2)	 ! 
           if(ismois.eq.3) etain = eta1
           call swflxb (precis,dels,fdrain,timi,etain,eta2, 
     &                    runoff,rinfil,rexfil)
        else                           ! threshold infiltration
           call swflxa (precis,dels,fdrain,eta1,eta2,
     &                    runoff,rinfil,rexfil)
        end if

c Modify SOIL WATER STORAGES and THERMAL STORAGES:
c        qinf, qexf water flux in mm/s
!     call soilsnow(kij,istype,dels,timi,froot,gfluxkf,precis,fev,fes,
!    &               runoff, xgx )
      stop 'call soilsnowv not ready for scam'
c Constrain eta1, eta2 =< ssat (saturation excess runoff) 
c for Mahrt Pan (zsw1=0.1), NP (zsw1=0) water lost to Dunn Black runoff
c for Richards equation, water loss should be = 0, but cacluated using 
c eta(i) (real8) above
ckf        if (ismois.le.2) then
ckf           rnofdb = 1000.*zsw1*(eta1 - min(eta1,ssat(istype)))
ckf           rnofdb = rnofdb + 
ckf     *           1000.*(zsm-zsw1)*(eta2 - min(eta2,ssat(istype)))
ckf        endif
ckf        runoff = runoff + rnofdb

c Assign returned storage values:
        store(kij,1) = cansto
ckf        dummy = cansto 
ckf        do i=1,ms
ckf           dummy  = dummy + dz(i)*wb(kij,i)*1000.
ckf        enddo 
ckf        store(kij,2) = dummy
ckf        store(kij,3) = totsti                   ! totsto on entry       

c       print *,'set up scrn stuff'

c Screen Temperature, windspeed and relative humidity  at 1.8 m
        tstar = -fh/(rho*capp*us)
        qstar = -fe/(rho*rlam*us)
        zscrn = max(z0m,1.8-disp)
        denom = (alog(za/zscrn) - psis(zetar(niter)) 
     &           + psis(zetar(niter)*zscrn/za))/vonk
        tscrn = ta - tstar*denom
        qscrn = (qa - qstar*denom)/qsatf(tscrn,pmb)
        uscrn = max(0.,uam - us*denom)

c Values at 3.0 m
        z3m   = max(z0m,3.0-disp)
        denom = (alog(za/z3m) - psis(zetar(niter)) 
     &           + psis(zetar(niter)*z3m/za))/vonk
        u3m   =  max(0.,uam - us*denom)

c Values at 10 m
        z10m  = max(z0m,10.-disp)
        denom = (alog(za/z10m) - psis(zetar(niter)) 
     &           + psis(zetar(niter)*z10m/za))/vonk
        u10m  =  max(0.,uam - us*denom)

c Radiative temperature very similar to average         
        trad  = (1.-trans)*t0v+trans*ts1
        
c Individual-patch outputs:
c       print *,'set up outvar'
        outvar(kp,1)  = fn         ! ACTUAL net irrad (+ down)   [W/m2]
        outvar(kp,2)  = gfluxkf   ! soil heat flux (+ down)     [W/m2]
        outvar(kp,3)  = fh         ! sensible heat flux          [W/m2]
        outvar(kp,4)  = fe         ! latent heat flux            [W/m2]
        outvar(kp,5)  = fc         ! mass CO2 flux         [kgCO2/m2/s]
        outvar(kp,6)  = us*us      ! (friction velocity)**2      [m2/s2]
        outvar(kp,7)  = zetar(niter)   ! Monin-Obukhov za/L      [-]
        outvar(kp,8)  = rbv        ! veg-bulk bound-layer resist [s/m]
        outvar(kp,9)  = rsv        ! veg-bulk stomatal resist    [s/m]
        outvar(kp,10) = rt0        ! turb resistance: (0,d)      [s/m]
        outvar(kp,11) = rt1        ! turb resistance: (d,za)     [s/m]
        outvar(kp,12) = runoff     ! runoff in this step         [mm]
        outvar(kp,13) = rinfil     ! surface infilt this step    [mm]
        outvar(kp,14) = rexfil     ! exfilt (deep drainage)      [mm]
        outvar(kp,15) = evap       ! total evaporation this step [mm]
        outvar(kp,16) = fnv        ! veg  parts of fn,fh,fe      [W/m2]
        outvar(kp,17) = fhv
        outvar(kp,18) = fev 
        outvar(kp,19) = fns        ! soil parts of fn,fh,fe      [W/m2]
        outvar(kp,20) = fhs
        outvar(kp,21) = fes
        outvar(kp,22) = 0. ! 1000.0*dels*qesatm  ! atm-lim soil fe
        outvar(kp,23) = 0. ! 1000.0*dels*qessoi  ! soil-lim soil fe
        outvar(kp,24) = fevc
        outvar(kp,25) = fevw
        outvar(kp,26) = fsu
        outvar(kp,27) = flu        ! flu uses old t0v and ts1
        outvar(kp,28) = t0v        ! veg surface t[C], q,deficit [kg/kg]
        outvar(kp,29) = q0v
        outvar(kp,30) = d0v
        outvar(kp,31) = tv         ! veg air layer t,q,deficit
        outvar(kp,32) = qv
        outvar(kp,33) = dv     
        outvar(kp,34) = zetar(1)   ! iterations toward zeta
        outvar(kp,35) = zetar(2)
        outvar(kp,36) = zetar(3)
        outvar(kp,37) = zetar(4)
        outvar(kp,38) = zetar(5)
        outvar(kp,39) = zetar(6)
        outvar(kp,40) = zetar(7)
        outvar(kp,41) = zetar(8)
        outvar(kp,42) = zetar(9)
        outvar(kp,43) = zetar(10)
c temporary output fields
        outvar(kp,44) = 0.0 !rnofdb not defined    ! dunn black runoff 0.!zl      
        outvar(kp,45) = fcv        ! veg carbon flux
        outvar(kp,46) = fcs        ! soil carbon flux
        outvar(kp,47) = faceta     ! soil water constraint param
c DARLAM necessary output
        outvar(kp,48) = us*us/uam  ! transfer coeffient for darlam
        outvar(kp,49) = trad       ! effective radiative temperature
        outvar(kp,50) = tscrn      ! screen temperature 
        outvar(kp,51) = qscrn      ! screen relative humidity
        outvar(kp,52) = uscrn      ! screen windspeed
        outvar(kp,53) = u3m        ! 3m windspeed
        outvar(kp,54) = u10m       ! 10m windspeed
10    continue                     ! end of loop through patches

c Calculate patch averages of all output variables and stores, using
c fractional area weighting. This may not be consistent for all 
c variables but troublesome averages (such as z/L) should be ignored
      do j=1,noutva
        outij(j) = 0.0
        do kp=1,npij
          kij = kp + kpoff         ! offset of current gridcell 
          pfrac  = param2(kij,4)   ! patch area fraction of patch kp
          outij(j) = outij(j) + outvar(kp,j)*pfrac
        end do
      end do

c use this if grid averaged values are needed for darlam ie soil temp 
ckf      do j=1,nstore
ckf        storij(j) = 0.0
ckf        do kp=1,npij
ckf          kij = kp + kpoff         ! offset of current gridcell 
ckf          pfrac  = param2(kij,4)   ! patch area fraction of patch kp
ckf          storij(j) = storij(j) + store(kij,j)*pfrac
ckf        end do
ckf      end do
    
      return
      end

c#######################################################################

      subroutine isscam (init,npij,kpoff,zeta0,z0m)
c-----------------------------------------------------------------------
c INITIALISATION routine
c INPUTS:
c   init   ! initialisation flag: init<0 for run, init>0 for section
c   npij   ! number of patches for current gridcell (1=single column)
c   kpoff  ! offset of current gridcell in patch array (0=single column)
c   param0 ! global (patch-independent) parameters
c   param1 ! default local (patch-dependent) parameters
c   param2 ! patch-dependent local params for 2D version
c VARIABLES INPUT AND MODIFIED ON OUTPUT:
c   store  ! array of heat, moisture storages 
c   sname  ! names of variables in store
c OUTPUTS:
c   hold:  ! array of parameters computed on first call and saved
c   hname  ! names of variables in hold
c   zeta0  ! initial z/L
c OUTPUTS via common blocks:
c   all properties in /canopy/ computed and saved throughout
c   all properties in /soiblk/ computed and saved throughout
c   properties zsw1,zse1 in /swblk/ computed and saved throughout
c-----------------------------------------------------------------------
      include 'newmpar.h'   ! parameter statement darlam npatch
      include 'parm.h'      ! ktau,idjd for diags
      include 'scamdim.h'    ! parameter statement for all dimensions
      include 'scampar.h'   ! dimension statement for blow commented out var
      include 'soilsnow.h'  ! tgg,wb soil temperature and moisture

      common /canopy/
     &    rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax(0:44),xjmax0(0:44)
      common /swblk/ istype,zsw1,zse1,zsm,qes,qevc,qrunf,qinf,qexf
c-----------------------------------------------------------------------
c Section initialisations: carried out at start of each section
c return initialised global parameters for all patches, via array hold:
      do kp=1,npij 
        kij    = kp + kpoff        ! offset of current gridcell 
        ivtype = param2(kij,1)     ! 44 for ST, (0,31) D versions
        tsoil  = 0.5*(tgg(kij,ms-2)+tgg(kij,ms-1))
c       print*,'in isscam tsoil',tsoil
c       print*,'in isscam smass,ssdn',kij,
c    &  smass(kij,1),smass(kij,2),smass(kij,3),
c    &  ssdn(kij,1),ssdn(kij,2),ssdn(kij,3)
        sdep   = max(0., smass(kij,1)/ssdn(kij,1) + 
     &                   smass(kij,2)/ssdn(kij,2) + 
     &                   smass(kij,3)/ssdn(kij,3) )  ! snow depth 
        call vegc (ivtype,tsoil,sdep,z0m)

        hold(kij,1) = rlai        ! xrlai(ivtype)
        hold(kij,2) = hc          ! variables for current patch
        hold(kij,3) = disp        ! xdisp(ivtype)
        hold(kij,4) = usuh        ! xusuh(ivtype)
        hold(kij,5) = coexp       ! xcoexp(ivtype)
        hold(kij,6) = zruffs      ! xzrufs(ivtype)
        hold(kij,7) = trans       ! xtrans(ivtype)
        hold(kij,8) = rt0us       ! xrt0us(ivtype)
        hold(kij,9) = rt1usa      ! xt1usa(ivtype)
        hold(kij,10)= rt1usb      ! xt1usb(ivtype)
      end do
      
c Global initialisations: carried out only once at start of run
      if (init.lt.0) then
                       
c constant soil properties, returned via common /soiblk/
c     call insoil !soilc   call now in indata
c     print *,'in sscam after call to insoil',npij,kpoff
      do 3 kp=1,npij
        kij = kp + kpoff               ! offset of current gridcell 

c initial storages: convert etarel, canrel to eta, cansto
        canstl = param1(3)
        cansat = canstl * rlai
c       print*,canstl,cansat,store(kij,1),kij,kp,kpoff
        canrel = store(kij,1)       ! rel canopy storage        [-]
        cansto = canrel*cansat
        cansto = max(cansto,0.0)
        cansto = min(cansto,cansat)
ckf        totsto = cansto 
ckf        do i=1,ms
ckf           totsto = totsto + dz(i)*wb(kij,i)*1000.
ckf        enddo 
        store(kij,1) = cansto
ckf        store(kij,2) = totsto       ! to be value at end of step
ckf        store(kij,3) = totsto       ! to be value at start of step
3     continue

      end if
c End of global initialisations

c Write initial parameters for the first 5 patches in this cell
c to check that the correct ones are being passed
c      kpmax = min(npij,5)
c      kpmax = 1 
c      write(*, 104) 'at the end of iscam'
c      write(*, 100) (kp+kpoff,kp=1,kpmax)                  ! patch
c      write(*, 101) (nint(param2(kp+kpoff,1)),kp=1,kpmax)  ! ivtype
c      write(*, 102) (nint(param2(kp+kpoff,2)),kp=1,kpmax)  ! istype

c      do i=1,nhold
c        if (i.eq.1.or.i.eq.2.or.i.eq.7.or.i.eq.8)
c     &    write(*, 103) hname(i),(hold(kp,i),kp=1,kpmax)
c      end do
c      do i=1,nstore
c        if (i.eq.1.or.i.eq.2)
c     &    write(*, 103) sname(i),(store(kp+kpoff,i),kp=1,kpmax)
c      end do
c      write(99,104)
c      write(99,100) (kp+kpoff,kp=1,kpmax)                  ! patch
c      write(99,101) (nint(param2(kp+kpoff,1)),kp=1,kpmax)  ! ivtype
c      write(99,102) (nint(param2(kp+kpoff,2)),kp=1,kpmax)  ! istype
c      do i=1,nhold
c        write(99,103) hname(i),(hold(kp,i),kp=1,kpmax)
c      end do
c      do i=1,nstore
c        write(99,103) sname(i),(store(kp,i),kp=1,kpmax)
c      end do
104   format(' SCAM section initialisation for first 5 patches:')
100   format('  patch     ',5i11)
101   format('  ivtype    ',5i11)
102   format('  istype    ',5i11)
103   format((2x,a8,2x,5f11.5))

      return
      end

c#######################################################################

      subroutine gsjarv (fsd,ts,ds,cs,gsmax,cextin,trans,fcfsd,fcdef,
     &                   fctemh,fctemm,faceta,cint,volm,gbv,  gsv,fcv)
c-----------------------------------------------------------------------
c Jarvis formulation for canopy gs and assilimation
c mrr, 22-jan-97 (in SEBxx)
c 27-apr-97: modified for SCAM; 
c KF, 26-05-97
c-----------------------------------------------------------------------
c Inputs:
c   fsd    = downward solar irradiance on canopy [W/m2]
c   ts     = leaf temperature [deg C]
c   ds     = leaf surface saturation deficit [units as fcdef]
c   cs     = leaf surface [CO2] [mol/mol]
c   gsmax  = max LEAF conductance for H2O [m/s]
c   cextin = PAR extinction coefficient in canopy
c   trans  = canopy transmission = exp(-cextin*rlai)
c   fcfsd  = PAR constraint parameter [W/m2]
c   fcdef  = saturation deficit constraint parameter [kg/kg]
c   fctemh = high-temp limit for gs [deg C]
c   fctemm = median (optimal) temp for gs [deg C]
c   faceta = soil water constraint parameter (0 =< faceta =< 1)
c   cint   = intercellular CO2 concentration [mol/mol]
c   volm   = molar volume [m3/mol]
c   gbv    = boundary layer conductance for H2O [m/s]  KF
c Outputs:
c   gsv    = canopy conductance for H2O [m/s]
c   fcv    = mass photosynthetic flux of CO2 [kgCO2/m2/s, positive up]
c-----------------------------------------------------------------------
      parameter(rgas   = 8.3143,   ! universal gas constant    (J/mol/K)
     &          rmair  = 0.02897,  ! molecular weights: dry air (kg/mol)
     &          rmh2o  = 0.018016, !                    water   (kg/mol)
     &          rmco2  = 0.044022, !                    CO2     (kg/mol)
     &          diffr  = 1.60,     ! H2O/CO2 diffusivity ratio
     &          oxy    = 0.209)    ! ambient oxygen conc       (mol/mol)
      zz     = (cextin*fsd+fcfsd)/(cextin*trans*fsd+fcfsd)
      facfsd = alog(zz)/cextin
      facdef = 1.0/(1.0+ds/fcdef)  ! def constraint
      factem = 1.0 - ((ts-fctemm)/(fctemh-fctemm))**2
      factem = max(factem,0.0)     ! temp contraint
      factor = facfsd*facdef*factem*faceta
c      factor = max(factor,0.01)    ! constrain gsv >= 0.01*gsmax
      gsv    = gsmax*factor
      fcv    = (rmco2/volm)*(cint-cs)*gsv*gbv/(gbv*diffr+gsv+1e-12)
c      fcv    = (rmco2/volm)*(cint-cs)*gsv/diffr
      return
      end

c#######################################################################

      subroutine PHOTOL2 (par,ts,ds,cs,swpar,ds0,
     &                    rjmax0,alfaqy,g0,a1,volm,ifh,
     &                    an,anc,anq,rd,ci,gsc,gsq,fc)
c-----------------------------------------------------------------------
c mrr, 15-jan-97
c Operational version of FCBBL leaf photosynthesis model
c See notes, 15-jan-97
c Method: Solve quadratic equation for net assimilation rate an, under
c both Rubisco and light limitation. Correct quadratic root is the
c larger one. Choose smaller an (from the Rubisco and light limited
c cases, anc and anq), then find gsc, ci.
c-----------------------------------------------------------------------
c Inputs:
c   par    = incoming photon PAR irradiance on leaf [mol/m2/s]
c   ts     = leaf temperature [deg C]
c   ds     = leaf surface saturation deficit [units as ds0]
c   cs     = leaf surface [CO2] [mol/mol]
c   swpar  = soil water stress parameter (0 < swpar < 1)
c   ds0    = sat deficit function parameter [units as ds]
c   rjmax0 = potential electron transport at T0 [mol/m2/s]
c   alfaqy = quantum yield
c   g0,a1  = intercept, slope params in Leuning (1995) model for gsc
c   volm   = molar volume (if volm.le.0, use default value 0.025 m3/mol)
c   ifh    = 0,1 for (Farquhar etal 1980),(Harley etal 1992) temp coeffs
c Outputs:
c   an     = net C assimilation rate [mol/m2/s] = min(anc,anq)
c   anc    = Rubisco-catalytic-activity limited an [mol/m2/s]
c   anq    = electron transport (light) limited an [mol/m2/s]
c   rd     = repiration [mol/m2/s]
c   ci     = internal [CO2] [mol/mol]
c   gsc    = stomatal conductance for CO2 [mol/m2/s]
c   gsq    = stomatal conductance for H2O [m/s]
c   fc     = CO2 flux in mass units [kgCO2/m2/s], positive outgoing
c Parameters: Leuning (1995, 1997), Farquhar et al (1980),
c             Harley et al (1992)
c Units: all molar for main calcs, with an positive incoming to leaf.
c        Exception: returned quantities fc and gsq are in met units
c        (fc in kgCO2/m2/s, positive outgoing, and gsq in m/s).
c-----------------------------------------------------------------------
c General parameters:
      parameter(rgas   = 8.3143,   ! universal gas constant    (J/mol/K)
     &          rmair  = 0.02897,  ! molecular weights: dry air (kg/mol)
     &          rmh2o  = 0.018016, !                    water   (kg/mol)
     &          rmco2  = 0.044022, !                    CO2     (kg/mol)
     &          diffr  = 1.60,     ! H2O/CO2 diffusivity ratio
     &          oxy    = 0.209)    ! ambient oxygen conc       (mol/mol)
c Begin code:
c Parameters for PAR (photosynthetic electron transport) limited A:
c      rjmax0 = 100e-6        ! [mol/m2/s]  = potl elect transport at T0
c      alfaqy = 0.20          ! [-]         = quantum yield
      theta  = 0.95          ! [-]         = shape factor for J
      absorb = 0.86          ! [-]         = leaf absorptance to PAR
      if (ifh.eq.0) then     ! Farquhar et al (1980) temp dependences
        hvj    = 37000       ! [J/mol]     = jmax activation energy
        hdj    = 220000      ! [J/mol]     = jmax deactivation energy
        svj    = 710         ! [J/mol]     = jmax entropy term
      elseif (ifh.eq.1) then ! Harley et al (1992) temp dependences
        hvj    = 79500       ! [J/mol]     = jmax activation energy
        hdj    = 201000      ! [J/mol]     = jmax deactivation energy
        svj    = 650         ! [J/mol/K]   = jmax entropy term
      else
        stop 'photol2: illegal ifh'
      end if
c Parameters for Rubisco-limited assimilation anc:
c To set vcmax0, use result of Leuning (1997), Table 2
      vcmax0 = rjmax0/2.44   ! [mol/m2/s]  = max Rubisco activity at T0
      rkc0   = 300e-6        ! [mol/mol]   = Michaelis constant for CO2
      rko0   = 256e-3        ! [mol/mol]   = Michaelis constant for O2
      hkc    = 59430         ! [J/mol]     = activation energy for CO2
      hko    = 36000         ! [J/mol]     = activation energy for O2
      if (ifh.eq.0) then     ! Farquhar et al (1980) temp dependences
        hvc    = 58520       ! [J/mol]     = vcmax activation energy
        hdc    = 0           ! [J/mol]     = (undefined)
        svc    = 0           ! [J/mol/K]   = (undefined)
      elseif (ifh.eq.1) then ! Harley et al (1992) temp dependences
        hvc    = 116300      ! [J/mol]     = vcmax activation energy
        hdc    = 202900      ! [J/mol]     = vcmax deactivation energy
        svc    = 650         ! [J/mol/K]   = vcmax entropy term
      else
        stop 'photol2: illegal ifh'
      end if
c Parameters for respiration:
        rd0    = 0.32e-6     ! [mol/m2/s]  = respiration at T0
        hrd    = 53000       ! [J/mol]     = activation energy for Rd
        trd2   = 40 + 273.16 ! [K]         = cutoff high temp for Rd
        hrd2   = 200000      ! [J/mol]     = param for high-T Rd cutoff
        srd2   = hrd2/trd2   ! [J/mol/K]   = param for high-T Rd cutoff
c Parameters for temp dependence of no-resp compensation point ggamms
        gamma0 = 34.6e-6     ! [mol/mol]
        gamma1 = 0.0451      ! [1/degK]
        gamma2 = 0.000347    ! [1/degK2]
c Parameters for Leuning (1995) stomatal conductance model,
c gsc = g0 + a1*(1/(1+ds/ds0))*an/(cs-ggamma):
c        a1     = 5.0         ! [-]         = slope parameter
c        g0     = 0.02        ! [mol/m2/s]  = intercept parameter
c Calculate all temperature-dependent quantities at leaf temp ts
c Use ifh to set denom term2c=1 for Farquhar vcmax temp dependence
      t0     = 20            ! reference temperature
      tk0    = t0 + 273.16
      tks    = ts + 273.16
      term1c = exp(hvc*(1-tk0/tks)/(rgas*tk0))
      term2c = 1 + ifh * exp((svc*tks-hdc)/(rgas*tks))
      vcmax  = vcmax0 * term1c/term2c
      term1j = exp(hvj*(1-tk0/tks)/(rgas*tk0))
      term2j = 1 + exp((svj*tks-hdj)/(rgas*tks))
      rjmax  = rjmax0 * term1j/term2j
      rkc    = rkc0 * exp(hkc*(1-tk0/tks)/(rgas*tk0))
      rko    = rko0 * exp(hko*(1-tk0/tks)/(rgas*tk0))
      term1d = exp(hrd*(1-tk0/tks)/(rgas*tk0))
      term2d = 1 + exp((srd2*tks-hrd2)/(rgas*tks))
      rd     = rd0 * term1d/term2d
      ggamms = gamma0 * (1 + gamma1*(ts-t0) + gamma2*(ts-t0)**2)
      ggamma = 1.1 * ggamms
c apply water stress factor to vcmax, rjmax, rd
      vcmax  = vcmax * swpar
      rjmax  = rjmax * swpar
      rd     = rd * swpar
c Calculate electron transport rate rj (smaller root of quadratic)
      apar   = absorb*par
      discr  = (alfaqy*apar+rjmax)**2 - 4*theta*alfaqy*apar*rjmax
      rj     = ((alfaqy*apar+rjmax) - sqrt(discr)) / (2*theta*theta)
c Calculate function ff (0 <= ff <= 1) accounting for sat deficit
      ff     = 1/(1+ds/ds0)
c Solve quadratic equation for net assimilation rate an, under both 
c Rubisco and light limitation. Correct quadratic root is larger one. 
      do i=1,2
        if (i.eq.1) then
          beta = vcmax       ! Rubisco-limited photosynthesis
          gama = rkc * (1 + oxy/rko)
        else
          beta = rj/4        ! light-limited photosynthesis
          gama = 2*ggamms
        end if
        c1     = cs + gama
        c2     = a1*ff/(cs-ggamma)
        c3     = cs - ggamms
        aa     = c1*c2 - 1
        bb     = c1*g0 + rd*(c1*c2-1) -beta*(c3*c2-1)
        cc     = rd*c1*g0 - beta*c3*g0
        discr  = bb**2 - 4*aa*cc
        if (i.eq.1) then     ! Rubisco-limited photosynthesis
          anc  = (-bb + sqrt(abs(discr)))/(2*aa)
        else                 ! light-limited photosynthesis
          anq  = (-bb + sqrt(abs(discr)))/(2*aa)
        end if
      end do
c Choose the smaller an from the Rubisco and light limited cases
      an  = min(anc,anq)
c Find gs (without allowing an<0) and ci
      gsc = g0 + a1*ff*max(an,0.0)/(cs-ggamma)
      ci  = cs - an/gsc
c Outputs in meteorological units
      if (volm.le.0) volm = 0.025
      gsq = volm*diffr*gsc   ! gs for H2O [m/s]
      fc  = -rmco2*an        ! CO2 flux [kg/m2/s], positive out
      return
      end

c#######################################################################

      subroutine PHINT1 (rlai,cextin,par,ts,ds,cs,swpar,ds0,
     &                   rjmax0,alfaqy,g0,a1,volm,ifh,
     &                   aan,aanc,aanq,rrd,cci,ggsc,ggsq,ffc)
c-----------------------------------------------------------------------
c mrr, 25-jan-97, ammended KF 21-Aug-97
c Integration of leaf photosynthesis model through canopy with
c exponential PAR distribution
c-----------------------------------------------------------------------
      parameter(rgas   = 8.3143,   ! universal gas constant    (J/mol/K)
     &          rmair  = 0.02897,  ! molecular weights: dry air (kg/mol)
     &          rmh2o  = 0.018016, !                    water   (kg/mol)
     &          rmco2  = 0.044022, !                    CO2     (kg/mol)
     &          diffr  = 1.60,     ! H2O/CO2 diffusivity ratio
     &          oxy    = 0.209)    ! ambient oxygen conc       (mol/mol)
      real*4 q(3),an(3),anc(3),anq(3),rd(3),ci(3),gsc(3),y(3)
c find 3 leaf PAR values
      q(3) = cextin*par                ! brightest leaf
      q(1) = q(3)*exp(-cextin*rlai)    ! dimmest leaf
      q(2) = (q(1)+q(3))/2             ! median leaf
      q3   = q(3)
      q1   = q(1)
c      write(*,100) cextin,rlai,par,q
c100   format(' cextin,rlai,par,q:',6f8.3)
c find photosynthesis properties for these 3 leaves
      do i=1,3
        call PHOTOL2 (q(i),ts,ds,cs,swpar,ds0,
     &                rjmax0,alfaqy,g0,a1,volm,ifh,
     &                an(i),anc(i),anq(i),rd(i),ci(i),gsc(i),gsq,fc)
      end do
c calculate integrals using parabolic fits
      do i=1,3               ! an
        y(i) = an(i)/max(q(i),1.0e-6)
      end do
      call parab (q(1),q(2),q(3),y(1),y(2),y(3),pa,pb,pc,err)
      aan  = (pa*(q3**3-q1**3)/3 + pb*(q3**2-q1**2)/2 
     &       + pc*(q3-q1)) / cextin
      do i=1,3               ! anc
        y(i) = anc(i)/max(q(i),1.0e-6)
      end do
      call parab (q(1),q(2),q(3),y(1),y(2),y(3),pa,pb,pc,err)
      aanc = (pa*(q3**3-q1**3)/3 + pb*(q3**2-q1**2)/2
     &       + pc*(q3-q1)) / cextin
      do i=1,3               ! anq
        y(i) = anq(i)/max(q(i),1.0e-6)
      end do
      call parab (q(1),q(2),q(3),y(1),y(2),y(3),pa,pb,pc,err)
      aanq = (pa*(q3**3-q1**3)/3 + pb*(q3**2-q1**2)/2
     &       + pc*(q3-q1)) / cextin
      do i=1,3               ! rd
        y(i) = rd(i)/max(q(i),1.0e-6)
      end do
      call parab (q(1),q(2),q(3),y(1),y(2),y(3),pa,pb,pc,err)
      rrd  = (pa*(q3**3-q1**3)/3 + pb*(q3**2-q1**2)/2
     &       + pc*(q3-q1)) / cextin
      do i=1,3               ! gsc
        y(i) = gsc(i)/max(q(i),1.0e-6)
      end do
      call parab (q(1),q(2),q(3),y(1),y(2),y(3),pa,pb,pc,err)
      ggsc = (pa*(q3**3-q1**3)/3 + pb*(q3**2-q1**2)/2
     &       + pc*(q3-q1)) / cextin
      cci  = cs - aan/max(ggsc,1.0e-6)
c Outputs in meteorological units
      if (volm.le.0) volm = 0.025
      ggsq  = volm*diffr*ggsc      ! gs for H2O [m/s]
      ffc   = -rmco2*aan           ! CO2 flux [kgCO2/m2/s, positive out]
      return
      end

c#######################################################################

      subroutine ruff(h,rlai,usuh,z0,d,coexp)
c-----------------------------------------------------------------------
c m.r. raupach, 24-oct-92
c see: Raupach, 1992, BLM 60 375-395
c      MRR notes "Simplified wind model for canopy", 23-oct-92
c      MRR draft paper "Simplified expressions...", dec-92
c-----------------------------------------------------------------------
c inputs:
c   h     = roughness height
c   rlai  = leaf area index (assume rl = frontal area index = rlai/2)
c output:
c   usuh  = us/uh (us=friction velocity, uh = mean velocity at z=h)
c   z0    = roughness length
c   d     = zero-plane displacement
c   coexp = coefficient in exponential in-canopy wind profile
c           U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
c           canopy and roughness-sublayer U(z) at z=h
c-----------------------------------------------------------------------
c preset parameters:
      parameter (cr    = 0.3,          ! element drag coefficient
     &           cs    = 0.003,        ! substrate drag coefficient
     &           beta  = cr/cs,        ! ratio cr/cs
     &           ccd   = 15.0,         ! constant in d/h equation
     &           ccw   = 2.0,          ! ccw=(zw-d)/(h-d)
     &           usuhm = 0.3,          ! (max of us/uh)
     &           vonk  = 0.4)          ! von Karman constant
      rl = rlai*0.5
c find uh/us
      usuhl  = sqrt(cs+cr*rl)
      USUH   = min(usuhl,usuhm)
c find d/h and d 
      xx     = sqrt(ccd*max(rl,0.0005))
      dh     = 1.0 - (1.0 - exp(-xx))/xx
      d      = dh*h
c find z0h and z0:
      psih   = alog(ccw) - 1.0 + 1.0/ccw
      z0h    = (1.0 - dh) * exp(psih - vonk/usuh)
      z0     = z0h*h
c find coexp: see notes "simplified wind model ..." eq 34a
      coexp  = usuh / (vonk*ccw*(1.0 - dh))
      return
      end

c#######################################################################

      subroutine vegc (iv,tsoil,sdep,z0m)
c KF, 1997
c For each vegetation type (= iv), assign veg height, total LAI, albedo,
c and computed aerodynamic, radiative and interception properties.
c Jmax0 assinged due to table by Ray Leuning and estimates  21-08-97 
c Apply seasonal variations in LAI and height. Return via /canopy/
c NB: total LAI = xrlai, vegLAI = xvlai, veg cover fraction = xpfc,
c     with xrlai = xvlai*xpfc
c Type  0 to 31: 2D version with Graetz veg types
c Type 32 to 43: 2D version with GCM veg types
c Type 44:       standalone version
c-----------------------------------------------------------------------
c   Name                             Symbol  Code hc:cm PFC:%  vegLAI
c   Ocean                                Oc     0     0     0  0.0
c   Tall dense forest                    T4     1  4200   100  4.8
c   Tall mid-dense forest                T3     2  3650    85  6.3
c   Dense forest                         M4     3  2500    85  5.0  (?)
c   Mid-dense forest                     M3     4  1700    50  3.75
c   Sparse forest (woodland)             M2     5  1200    20  2.78
c   Very sparse forest (woodland)        M1     6  1000     5  2.5
c   Low dense forest                     L4     7   900    85  3.9
c   Low mid-dense forest                 L3     8   700    50  2.77
c   Low sparse forest (woodland)         L2     9   550    20  2.04
c   Tall mid-dense shrubland (scrub)     S3    10   300    50  2.6
c   Tall sparse shrubland                S2    11   250    20  1.69
c   Tall very sparse shrubland           S1    12   200     5  1.9
c   Low mid-dense shrubland              Z3    13   100    50  1.37
c   Low sparse shrubland                 Z2    14    60    20  1.5
c   Low very sparse shrubland            Z1    15    50     5  1.21
c   Sparse hummock grassland             H2    16    50    20  1.58
c   Very sparse hummock grassland        H1    17    45     5  1.41
c   Dense tussock grassland              G4    18    75    85  2.3
c   Mid-dense tussock grassland          G3    19    60    50  1.2
c   Sparse tussock grassland             G2    20    45    20  1.71
c   Very sparse tussock grassland        G1    21    40     5  1.21
c   Dense pasture/herbfield (perennial)  F4    22    60    85  2.3
c   Dense pasture/herbfield (seasonal)  F4S    23    60    85  2.3
c   Mid-dense pasture/herb (perennial)   F3    24    45    50  1.2
c   Mid-dense pasture/herb  (seasonal)  F3S    25    45    50  1.2
c   Sparse herbfield*                    F2    26    35    20  1.87
c   Very sparse herbfield                F1    27    30     5  1.0
c   Littoral                             LL    28   250    50  3.0
c   Permanent lake                       PL    29     0     0  0
c   Ephemeral Lake (salt)                SL    30     0     0  0
c   Urban                                 U    31     0     0  0
c   STAND ALONE: hc,rlai from param1      -    44     -   100  -
c-----------------------------------------------------------------------
      include 'newmpar.h'   ! parameter statement darlam npatch
      include 'scamdim.h' 
      include 'scampar.h'
      real xhc(0:44),xpfc(0:44),xvlai(0:44),xslveg(0:44)
      common /canopy/
     &    rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax(0:44),xjmax0(0:44)
c Aerodynamic parameters, diffusivities, water density:
      parameter(vonk   = 0.40,     ! von Karman constant
     &          a33    = 1.25,     ! inertial sublayer sw/us
     &          csw    = 0.50,     ! canopy sw decay (Weil theory)
     &          ctl    = 0.40)     ! Wagga wheat (RDD 1992, Challenges)
c Vegetation height
      data xhc    / 0.0,
     &             30.0,28.0,25.0,17.0,12.0,10.0, 9.0, 7.0, 5.5, 3.0,
     &              2.5, 2.0, 1.0, 0.6, 0.5, 0.5,0.45,0.75, 0.6,0.45,
     &              0.4, 0.6, 0.6,0.24,0.25,0.35, 0.3, 2.5, 0.0, 0.0,
     &              0.0, 
     & 32.0,20.0,20.0,17.0,17.0, 1.0, 1.0, 1.0, 0.5, 0.6, 0.0, 1.0,0./ !Sellers 1996 J.Climate 

c Vegetation fractional cover
      data xpfc   /0.00,
     &             1.00,0.85,0.85,0.50,0.20,0.05,0.85,0.50,0.20,0.50,
     &             0.20,0.05,0.50,0.20,0.05,0.20,0.05,0.85,0.50,0.20,
     &             0.05,0.85,0.85,0.50,0.50,0.20,0.05,0.50,0.00,0.00,
     &             0.00,
     &   .98,.75,.75,.75,.50,.86,.65,.79,.30,.42,.02,.54,1.0/

c Veg LAI from Graetz table of 283 veg types (iv=0 to 31), and maximum 
c Veg LAI for GCM veg types (iv=32 to 43)
      data xvlai  / 0.0,
     &              4.80,6.30,5.00,3.75,2.78,2.50,3.90,2.77,2.04,2.60,
     &              1.69,1.90,1.37,1.50,1.21,1.58,1.41,2.30,1.20,1.71,
     &              1.21,2.30,2.30,1.20,1.20,1.87,1.00,3.00,0.00,0.00,
     &              0.00,
     &   6.0,5.0,4.0,4.0,4.0,3.0,3.0,3.0,1.0,4.0,0.5,3.0,0.0/

c For seasonally varying LAI, amplitude of veg LAI seasonal change
      data xslveg  /0.00,
     &              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &              0.00,
     &   2.0,2.0,2.0,2.0,2.0,1.5,1.5,1.5,1.0,0.5,0.5,0.5,0.0/
! Leaf gsmax for forest (0.006), grass (0.008) and crop (0.012)
! littoral is regarded as forest, dense pasture between grass and crop
!     data xgsmax  / 0.0,
!    &     0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,
!    &     0.006,0.006,0.006,0.006,0.006,0.008,0.008,0.008,0.008,0.008,
!    &     0.008,0.010,0.010,0.012,0.012,0.008,0.008,0.006,0.000,0.000,
!    &     0.000,
!    &  .006,.006,.006,.006,.006,.006,.008,.008,.006,.006,.0,0.010,0.0/
! littoral is regarded as forest, dense pasture between grass and crop
!     data xjmax0 / 0.0,
!    &     5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,
!    &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,10e-5,10e-5,
!    &     10e-5,15e-5,15e-5,15e-5,15e-5,10e-5,10e-5,5e-5,1e-5,1e-5,
!    &     1e-5,
!    &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,5e-5,10e-5,
!    &     1e-5,15e-5,1e-5/
!     &     25e-5,25e-5,20e-5,15e-5,25e-5,7e-5,7e-5,7e-5,15e-5,15e-5,
!     &     7e-5,25e-5,1e-5/ !Sellers 1996 J.Climate, I think they are too high
c-----------------------------------------------------------------------
c Global aerodynamic parameters:
      za         = param0(1)           ! ref ht (z=0 at disp)    [m]
c     zdlin      = param0(16)          ! z/disp below which TL linear
      z0soil     = param1(4)           ! soil z0                 [m]
c Assign canopy properties for standalone runs:
ckf needs workd      if (iv.eq.44) then 
c         hc    = param1(1)           ! standalone canopy ht    [m]
c         lai   = param1(2)           ! standalone LAI          [-]
c         xgsmax(44) = param1(16)          ! standalone gsmax        [m/s]      
c         xjmax0(44) = param1(32)          ! standalone rjmax0             
c      endif

c Assign aerodynamic, radiative, stomatal, interception properties
c Assign total LAI (xrlai) from veg LAI and PFC, and assign seasonal 
c   variation in LAI and veg height where necessary. This is controlled
c   by the factor season (0 =< season =< 1).
c      tsoil =.5*(ts4+ts5) ! deep soil temp
      ftsoil=max(0.,1.-.0016*(298.-tsoil)**2)
      if( tsoil .ge. 298. ) ftsoil=1.
      vrlai = max(0.0,(xvlai(iv)-xslveg(iv)*(1.-ftsoil))*xpfc(iv))
      hc    = max(0.0,xhc(iv) - sdep)
      rlai  = vrlai*hc/max(0.01,xhc(iv))
c   Find roughness length z0m from hc and rlai:
      call ruff(hc,rlai,usuh,z0m,disp,coexp)
c   Set aerodynamic variables for bare soil and vegetated cases:
      if (rlai.lt.0.001 .or. hc.lt.z0soil) then
          z0m    = z0soil      ! BARE SOIL SURFACE
          hc     = 0.0  
          rlai   = 0.0
          rt0us  = 0.0                                             
          disp   = 0.0
          zruffs = 0.0
          rt1usa = 0.0
          rt1usb = 0.0
c these quantities cant be zero since they are in the denominator later on
c          usuh   = 0.0
c          coexp  = 0.0 
      else                       ! VEGETATED SURFACE
c   rt0 = turbulent resistance from soil (z0 = 0) to canopy
c   (z1 = zero-plane displacement), normalized as rt0us=rt0*us
c   nb: term4 added 13-sep-95 to make TL proportional to z near ground.
c   nb: term5 added 03-oct-96 to account for sparse canopies. Constant
c       length scale ctl*hc replaced by ctl*(3/2)*disp, taking effect
c       when disp<(2/3)*hc, or total LAI < 1.11. Otherwise, term5=1. 
          term1  = exp(2*csw*rlai)
          term2  = exp(2*csw*rlai*(1-disp/hc))
          term3  = a33**2*ctl*2*csw*rlai
          term4  = zdlin*alog(zdlin*disp/z0soil) + (1-zdlin)
          term5  = max((2.0/3.0)*hc/disp, 1.0)
          rt0us = term5*term4*(term1-term2)/term3
c   rt1 = turbulent resistance from canopy (z1 = disp) to
c   reference level za (from disp as origin). Normalisation:
c     rt1us = us*rt1 = rt1usa + rt1usb + rt1usc
c   with term a = resistance from disp to hc
c     term b = resistance from hc to zruffs (or za if za<zruffs)
c     term c = resistance from zruffs to za (= 0 if za<zruffs)
c   where zruffs = SCALAR roughness sublayer depth (ground=origin)
c          xzrufs(iv) = xdisp(iv) + xhc(iv)*a33**2*ctl/vonk
          zruffs = disp + hc*a33**2*ctl/vonk/term5
          rt1usa = term5*(term2 - 1.0)/term3
          rt1usb = term5*(min(za+disp,zruffs)-hc)
     &                  /(a33**2*ctl*hc)
          rt1usb = max(rt1usb,0.0)       ! in case zrufs < hc
      end if
c   Canopy transmission = shortwave radiation fraction to ground
      cextin = param1(9)
      trans  = exp(-cextin*rlai)
c   Canopy water storage capacity, from canstl = canopy water store 
c   per unit LAI [mm]
ckf      canstl = param1(3)
ckf      cansat = canstl * rlai
      return
      end

c#######################################################################

      SUBROUTINE AIR(TC,PMB,RHO,VOLM,CAPP,RLAM,qsat,epsi,visc)
C MRR, 09-NOV-88
c 27-jan-94: transfer qsat and epsi to their own function subroutines
c 04-mar-94: add kinematic viscosity of air
c-----------------------------------------------------------------------
c Find properties of air at given temperature and pressure:
c see Garratt (1992), "The Atmospheric Boundary Layer"
c Inputs:
c   TC   = Air temperature (DEG C)
c   PMB  = Air pressure (mb)
c Outputs:
c   RHO  = DRY AIR DENSITY (KG M-3)
c   VOLM = MOLAR VOLUME (M3 MOL-1)
c   CAPP = SPECIFIC HEAT OF AIR AT CONSTANT PRESSURE (J KG-1 K-1)
c   RLAM = LATENT HEAT FOR WATER (J KG-1)
c   qsat = saturation specific humidity (kg/kg)
c   epsi = d(qsat)/dT ((kg/kg)/K)
c   visc = kinematic viscosity of air (m2/s)
c-----------------------------------------------------------------------
c Physical and mathematical constants 
      parameter(rgas   = 8.3143,       ! universal gas const  (J/mol/K)
     &          rmair  = 0.02897)      ! molecular wt: dry air (kg/mol)
      CAPP = 1004.0                    ! J/kg/K
      RLAM = (2501.0-2.38*TC)*1000.0   ! kg/m3
      TK   = 273.16+TC
      PPA  = 100.0*PMB                 ! pressure in Pascals
      RHO  = RMAIR*PPA/(RGAS*TK)       ! kg/m3
      VOLM = RGAS*TK/PPA               ! m3/mol
      qsat = qsatf(tc,PMB)
      epsi = epsif(tc,PMB)
      visc = (1e-5)*(1.35+0.0092*tc)   ! m2/s
      visc = max(visc,1e-5)            ! avoid negative visc
      RETURN
      END

c#######################################################################

      FUNCTION qsatf(TC,PMB)
C MRR, 1987
C AT TEMP TC (DEG C) AND PRESSURE PMB (MB), GET SATURATION SPECIFIC
C HUMIDITY (KG/KG) FROM TETEN FORMULA
      parameter (rmair  = 0.02897,         ! molecular wt: dry air (kg/mol)
     &           rmh2o  = 0.018016)        ! molecular wt: water   (kg/mol)
      parameter (A=6.106,B=17.27,C=237.3)  ! Teten coefficients
      ES    = A*EXP(B*TC/(C+TC))           ! sat vapour pressure (mb)
      qsatf = (rmh2o/rmair)*ES/PMB         ! sat specific humidity (kg/kg)
      RETURN
      END

c#######################################################################

      FUNCTION epsif(TC,PMB)
C MRR, 1987, 27-jan-94
C AT TEMP TC (DEG C) AND PRESSURE PMB (MB), GET epsi = 
c d(sat spec humidity)/dT ((kg/kg)/K) FROM TETEN FORMULA
      parameter(rmair  = 0.02897,      ! molecular wt: dry air (kg/mol)
     &          rmh2o  = 0.018016)     ! molecular wt: water   (kg/mol)
      parameter                        ! Teten coefficients
     & (A=6.106,B=17.27,C=237.3)  
      CAPP  = 1004.0                   ! J/kg/K
      RLAM  = (2501.-2.38*TC)*1000.0   ! kg/m3
      ES    = A*EXP(B*TC/(C+TC))       ! sat vapour pressure (mb)
      DESDT = ES*B*C/(C+TC)**2         ! d(sat VP)/dT: (mb/K)
      DQSDT = (RMH2O/RMAIR)*DESDT/PMB  ! d(sat spec hum)/dT: (kg/kg)/K
      epsif = (RLAM/CAPP)*DQSDT        ! dimensionless
      return
      end

c#######################################################################

      function psim(zeta)
c mrr, 16-sep-92 (from function PSI: MRR, EDINBURGH 1977)
C COMPUTES INTEGRATED STABILITY FUNCTION PSIm(Z/L) (Z/L=zeta)
C FOR MOMENTUM, USING THE BUSINGER-DYER FORM FOR UNSTABLE CASES 
C AND THE WEBB FORM FOR STABLE CASES. SEE PAULSON (1970).
      GU=16.0
      GS=5.0
      PI=3.14159265358979
      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable          
      stable = -gs*zeta
      x      = (1.0 + gu*abs(zeta))**0.25
      unstab = alog((1.0+x*x)*(1.0+x)**2/8) - 2.0*atan(x) + pi*0.5     
      psim   = z*stable + (1.0-z)*unstab
      RETURN
      END

c#######################################################################

      function psis(zeta)
c mrr, 16-sep-92 (from function PSI: MRR, EDINBURGH 1977)
C COMPUTES INTEGRATED STABILITY FUNCTION PSIs(Z/L) (Z/L=zeta)
C FOR SCALARS, USING THE BUSINGER-DYER FORM FOR UNSTABLE CASES 
C AND THE WEBB FORM FOR STABLE CASES. SEE PAULSON (1970).
      GU=16.0
      GS=5.0
      PI=3.14159265358979
      z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable          
      stable = -gs*zeta
      y      = (1.0 + gu*abs(zeta))**0.5
      unstab = 2.0 * alog((1+y)*0.5)
      psis   = z*stable + (1.0-z)*unstab
      RETURN
      END

c#######################################################################

      subroutine sevap (qv,ts1,eta1r,rho,rlam,rtsoil,pmbh,isevap,dbde,
     &                  fes,qes,beta)
c-----------------------------------------------------------------------
c MRR/LZ, 06-feb-95, 04-mar-95
c Includes calculation of soil evaporation by several schemes (details:
c Mahfouf and Noilhan (1991) JAM 30, 1354-64).
c Inputs: 
c   qv
c   ts1
c   eta1r
c   rho
c   rlam 
c   rtsoil
c   pmbh
c   isevap
c   dbde
c Outputs:
c   fes    = soil evap (W/m2)
c   qes    = soil evap (water m/s)
c   beta   = 1 for alfa scheme
c-----------------------------------------------------------------------
      parameter(rhow   = 1000.0,   ! liquid water density (kg/m3)
     &          pi     = 3.141592654)
c soil evaporation
c      if (isevap.lt.3) then !stop 'isevap=3 not yet enabled'
c   x=1 if isevap=1 (alfa), and x=0 if isevap=2,3 (beta, threshold)
c   note: use beta scheme for atm-lim part of threshold scheme
      x = 0.5 + sign(0.5,1.5-isevap)
      etacr  = 0.75                    ! rel eta at field cap, etacap
      eta1rc = eta1r/etacr             ! rel eta wrt etacap
      eta1rc = min(eta1rc,1.0)
      eta1rc = max(eta1rc,0.0)
      rhs    = 0.5*(1.0-cos(pi*eta1rc))
      betaa  = 1.0                     ! for alfa scheme
      betab  = dbde * eta1r            ! dbde = 4/3, 8 for D78, KGK91
      betab  = min(betab,1.0)
      betab  = max(betab,0.0)
      beta   = x*betaa + (1.0-x)*betab
      dq     = (x*rhs + (1.0-x)*betab) * qsatf(ts1,pmbh) - beta*qv
      dq     = max(dq,0.0)        
      fes    = rho*rlam*dq/rtsoil      ! soil evap: [W/m2]
      qes    = fes/(rhow*rlam)         ! soil evap: [water m/s]
      
      return
      end

c#######################################################################

      subroutine swflxa (precis,dels,fdrain,eta1,eta2,
     &                   runoff,rinfil,rexfil)
c ----------------------------------------------------------------------
c Find SOIL WATER FLUXES: use soil moistures for current step.
c Surface infiltration expression from Mahrt and Pan (1984).
c Note soil evap flux qes [m/s] already found in SEB calculation.
c INPUTS IN PARAMETER LIST:  
c     precis    precipitation reaching the ground
c     dels      time step in sec                            
c     fdrain    exfiltation parameter
c     eta1      soil moisture of the top layer
c     eta2      soil moisture of the bulk layer
c OUTPUTS IN PARAMETER LIST:
c     runoff    runoff in this step [mm]
c     rinfil    surface infiltration in this step [mm]
c     rexfil    deep drainage in this step [mm]
c INPUTS IN COMMON BLOCK swblk:
c     istype:   soil type
c     zsm:      accessible soil moisture depth (m)
c     qes:      actual soil evap (m/s)
c     qevc:     actual dry veg evap through roots (m/s)
c     qrunf:    runoff rate (m/s)
c     qinf:     surface infiltration (m/s)
c     qexf:     deep exfiltration below soil layers (m/s)
c ----------------------------------------------------------------------
      include 'newmpar.h'   ! parameter statement darlam npatch
      include 'scamdim.h'   ! parameter statement for all dimensions
      include 'scampar.h'   ! dimension statement for scam part
      include 'soilsnow.h'  ! tgg,wb soil temperature and moisture
      include 'soilv.h'
c     common/soilpr/
c    & rhos(9),css(9),cnsd(9),bch(9),ssat(9),swilt(9),hsbh(9),i2bp3(9),
c    & ibp2(9),sfc(9),sucs(9),hyds(9),silt(9),clay(9),sand(9),sormax(9),
c    & tgmax(9),scond1s(3),cap1s(3),ssatcur(ms+1),coef(ms+4),cnsw(ms),
c    & gamm(ms),zse(ms),zsh(npmax,ms+1)
      common /swblk/ istype,zsw1,zse1,zsm,qes,qevc,qrunf,qinf,qexf

      qpres  = 0.001*precis/dels     ! soil precip flux [m/s]
      deleta = max(0.,ssat(istype)-eta1)
      difsat = hsbh(istype)/ssat(istype)
      qpond  = difsat*deleta/(zse1/2) + hyds(istype)
c     &         (ssat(istype(1))-eta1)/(zse1/2) + hyds(istype(1))
      qinf   = min(qpres,qpond)      ! MP84 ponded infilt flux [m/s]
c ----------------------------------------------------------------------
c Campbell Soil model (eta/ssat)^(2b+3) and not etarel^(2b+3)
      eta2r  = eta2/ssat(istype)
      cond2  = hyds(istype) * eta2r**i2bp3(istype)

      qexf   = fdrain*cond2          ! deep exfilt flux [m/s]
      qrunf  = qpres - qinf          ! surface runoff flux [m/s]
      runoff = qrunf*dels*1000.0     ! surface runoff [mm]
      rinfil = qinf*dels*1000.0      ! surface infiltration [mm]
      rexfil = qexf*dels*1000.0      ! exfilt (drainage) [mm]
      return
      end
      
c#######################################################################

      subroutine swflxb (precis,dels,fdrain,timi,etain,eta2,
     &                   runoff,rinfil,rexfil)
c ----------------------------------------------------------------------
c Find SOIL WATER FLUXES using soil moistures for current step.
c Infiltration calculation due to Freeman Cook, with duration of 
c rain event based on tgmax, soil type dependent
c Note: soil evap flux qes [m/s] already found in SEB calculation.   
c       dry veg evap flux [m/s] as previously (swflxa) 
c       deep exfiltration as previously (swflxa)
c INPUTS IN PARAMETER LIST:  
c     precis    precipitation reaching the ground
c     dels      time step in sec                            
c     fdrain    parameter for deep drainage
c     timi      time from start of simulation
c     eta1      upper soil moisture
c     eta2      bulk soil moisture
c OUTPUTS IN PARAMETER LIST:
c     runoff    runoff in this step                        [mm]
c     rinfil    surface infiltration in this step          [mm]
c     rexfil    deep drainage in this step                 [mm]
c INPUTS IN COMMON BLOCK swblk:
c     istype:   soil type
c     zsm:      accessible soil moisture depth             (m)
c     qes:      actual soil evap                           (m/s)
c     qevc:     actual dry veg evap through roots          (m/s)
c     qrunf:    runoff rate                                (m/s)
c     qinf:     surface infiltration                       (m/s)
c     qexf:     deep exfiltration below soil layers        (m/s)
c INPUTS IN COMMON BLOCK swblkb:
c     ntsrt     start of rain event                        [s]        
c     ntend     end of rain event                          [s]
c     sorin     sorptivity at beginning of rain event      [m/s^0.5]
c     suminf    accumulated infiltration                   [m]
c ----------------------------------------------------------------------
      include 'newmpar.h'   ! parameter statement darlam npatch
      include 'scamdim.h'   ! parameter statement for all dimensions
      include 'soilv.h'
c     common/soilpr/
c    & rhos(9),css(9),cnsd(9),bch(9),ssat(9),swilt(9),hsbh(9),i2bp3(9),
c    & ibp2(9),sfc(9),sucs(9),hyds(9),silt(9),clay(9),sand(9),sormax(9),
c    & tgmax(9),scond1s(3),cap1s(3),ssatcur(ms+1),coef(ms+4),cnsw(ms),
c    & gamm(ms),zse(ms),zsh(npmax,ms+1)

      common /swblk/ istype,zsw1,zse1,zsm,qes,qevc,qrunf,qinf,qexf
      common /swblkb/	ksw,tsrt(npmax),tend(npmax),
     &     sorin(npmax),suminf(npmax)
c
c infiltration calculation
      if (timi.le.1.001) tend(ksw)=0.0     ! initialise end timer
      qpres= 0.001*precis/dels         ! rainfall rate [m/s]
      if (qpres.gt.0.0) then           ! it is raining
        if (timi.gt.tend(ksw)) then    ! new event: reset event marker
          tsrt(ksw)  = timi
          sorin(ksw) = sormax(istype)*sqrt((ssat(istype)-etain)/
     &                 (ssat(istype)-swilt(istype)))
                                       ! inital sorptivity [m/s^0.5]
          suminf(ksw) = 0.d0           ! set infilt accumulator
        end if
        
        tend(ksw) = timi+tgmax(istype)/dels ! end marker for this rain event
        time0 = (timi-tsrt(ksw))*dels     ! rain time to start of time step
        time1 = time0 + dels
        tc    = (sorin(ksw) / (2.*qpres-hyds(istype)))**2  
        tjoin = (sorin(ksw) / hyds(istype))**2              
        tp1   = 0.
        if (time1.lt.tc) then             ! no compression needed yet
          qinf   = qpres
          rinfil = qinf*dels*1000.        ! surface infiltration [mm] 
          qrunf  = 0.0                    ! surface runoff flux [m/s]
          runoff = 0.0                    ! surface runoff [mm]
        else    ! tc < time1              ! compression analysis
          if (tc.ge.tjoin) then			! passed joining time i=Ks
             qinf   = min(qpres,hyds(istype))
          else  ! tc < tjoin  			! compression analysis set up
             potinf = sorin(ksw)*sqrt(tc) + 0.5*hyds(istype)*tc
                                          ! potential infiltration [m]
             qinfpot = qpres                
             tp1    =  time0 + (potinf-suminf(ksw))/qpres
             if (tp1.ge.time0 .and. tp1.le.time1) ! ponding starts in this step
     &          qinfpot = (qpres*(tp1-time0)
     &                   +sorin(ksw)*sqrt(time1-tp1+tc)
     &                   -sorin(ksw)*sqrt(tc)
     &                   + 0.5*hyds(istype)*(time1-tp1))/dels
             if (tp1.lt.time0) 	! ponding started in earlier time step
     &          qinfpot = (sorin(ksw)*sqrt(time1-tp1+tc)
     &                   - sorin(ksw)*sqrt(time0-tp1+tc)
     &                   + 0.5*hyds(istype)*dels )/dels
             qinf   = max(hyds(istype),qinfpot)                                                  
             qinf   = min(qinf,qpres)
          endif
          rinfil = qinf*dels*1000.0        ! surface infiltration [mm]
          qrunf  = max(0.,qpres - qinf)    ! surface runoff flux [m/s]
          runoff = qrunf*dels*1000.0       ! surface runoff [mm]
        end if             
        suminf(ksw) = suminf(ksw) + qpres*dels
                                       ! accumulated infiltration [m]
      else                             ! it is not raining
        qinf   = 0.0
        rinfil = 0.0
        qrunf  = 0.0                   ! surface runoff flux [m/s]
        runoff = 0.0                   ! surface runoff [mm]
      end if

c exfiltration
c ----------------------------------------------------------------------
c Campbell Soil model (eta/ssat)^(2b+3) and not etarel^(2b+3)
      eta2r  = eta2/ssat(istype)
      cond2  = hyds(istype) * eta2r**i2bp3(istype)

      qexf   = fdrain*cond2                ! deep exfilt flux [m/s]
      rexfil = qexf*dels*1000.0            ! exfilt (drainage) [mm]
      return
      end

c#######################################################################

      subroutine parab(x1,x2,x3,y1,y2,y3,a,b,c,err)
c mrr, 24-nov-96
c find the parabola (y = a*x*x + b*x + c) passing through the three
c points (x1,y1), (x2,y2), (x3,y3)
      if (x1.eq.x2.or.x1.eq.x3.or.x2.eq.x3) then
        err = 1
        a   = 0
        b   = 0
        c   = 0
        return
      end if
      err = 0
      t31 = (y3-y1)/(x3-x1)
      t21 = (y2-y1)/(x2-x1)
      a   = (t31-t21)/(x3-x2)
      b   = t21 - a*(x1+x2)
      c   = y1 - b*x1 - a*x1*x1
      return
      end

c***********************************************************************

