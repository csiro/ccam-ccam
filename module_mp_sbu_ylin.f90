!WRF:MODEL_LAYER:PHYSICS

!--- The code is based on Lin and Colle (A New Bulk Microphysical Scheme 
!             that Includes Riming Intensity and Temperature Dependent Ice Characteristics, 2011, MWR)
!             and Lin et al. (Parameterization of riming intensity and its impact on ice fall speed using ARM data, 2011, MWR)
!--- NOTE: 1) Prognose variables are: qi,PI(precipitating ice, qs, which includes snow, partially rimed snow and graupel),qw,qr
!---       2) Sedimentation flux is based on Prudue Lin scheme 
!---       2) PI has varying properties depending on riming intensity (Ri, diagnosed currently following Lin et al. (2011, MWR) and T 
!---       3) Autoconverion is based on Liu and Daum (2004)         
!---       4) PI size distribution assuming Gamma distribution, but mu_s=0 (Exponential) currently
!---       5) No density dependent fall speed since the V-D is derived using Best number approach, which already includes density effect 
!---       6) Future work will include radar equivalent reflectivity using the new PI property (A-D, M-D, N(D)). If you use RIP for reflectivity 
!---          computation, please note that snow is (1-Ri)*qs and graupel is Ri*qs. Otherwise, reflectivity will be underestimated.      
!---       7) The Liu and Daum autoconverion is quite sensitive on Nt_c. For mixed-phase cloud and marine environment, Nt_c of 10 or 20 is suggested.
!---          default value is 10E.6. Change accordingly for your use.


MODULE module_mp_sbu_ylin
!USE    module_wrf_error
!
!..Parameters user might change based on their need
  real, parameter, private :: RH = 1.0
  real, parameter, private :: xnor = 8.0e6
  real, parameter, private :: Nt_c = 100.E6 !250.E6 
!..Water vapor and air gas constants at constant pressure
  real, parameter, private :: Rvapor = 461.5
  real, parameter, private :: oRv    = 1./Rvapor
  real, parameter, private :: Rair   = 287.04
  real, parameter, private :: Cp     = 1004.0
  real, parameter, private :: grav   = 9.81
  real, parameter, private :: rhowater = 1000.0
  real, parameter, private :: rhosnow  = 100.0
    
  real, parameter, private :: SVP1=0.6112
  real, parameter, private :: SVP2=17.67
  real, parameter, private :: SVP3=29.65
  real, parameter, private :: SVPT0=273.15
  real, parameter, private :: EP1=Rvapor/Rair-1.
  real, parameter, private :: EP2=Rair/Rvapor
!..Enthalpy of sublimation, vaporization, and fusion at 0C.
  real, parameter, private :: XLS = 2.834E6
  real, parameter, private :: XLV = 2.5E6
  real, parameter, private :: XLF = XLS - XLV
    
  !REAL, SAVE, PUBLIC :: qi0 = 1.0e-3   
  real, parameter, private ::                               &
             qi0 = 1.0e-3,                                  &   !--- ice aggregation to snow threshold
             xmi50 = 4.8e-10, xmi40 = 2.46e-10,             &
             xni0 = 1.0e-2, xmnin = 1.05e-18, bni = 0.5,    &
             di50 = 1.0e-4, xmi = 4.19e-13,                 &   !--- parameters used in BF process
             bv_r = 0.8, bv_i = 0.25,                       &
             o6 = 1./6.,  cdrag = 0.6,                      &
             avisc = 1.49628e-6, adiffwv = 8.7602e-5,       &
             axka = 1.4132e3, cw = 4.187e3,  ci = 2.093e3

interface parama1
  module procedure parama1_v, parama1_s
end interface parama1

CONTAINS

SUBROUTINE clphy1d_ylin(dt, imax,                           &
                      qvz, qlz, qrz, qiz, qsz,              &
                      thz, tothz, rho, orho, sqrho,         &
                      prez, zz, dzw, zsfc,                  &
                      precrz, preciz, precsz,               & !zdc20220116
                      EFFC1D, EFFI1D, EFFS1D, EFFR1D,       & !zdc 20220208
                      pptrain, pptsnow,pptice,              &
                      kts, kte, riz,                        &
                      ncz, nrz, niz, nsz,                   &
                      fluxr, fluxi, fluxs, fluxg, fluxm,    &
                      fluxf, fevap, fsubl, fauto, fcoll,    &
                      faccr, vi, vs, vg,                    &
                      zpsnow,zpsaut,zpsfw,zpsfi,zpraci,     & !process rate 
                      zpiacr,zpsaci,zpsacw,zpsdep,          &
                      zpssub,zpracs,zpsacr,zpsmlt,          &
                      zpsmltevp,zprain,zpraut,zpracw,       &
                      zprevp,zpgfr,zpvapor,zpclw,           &
                      zpladj,zpcli,zpimlt,zpihom,           &
                      zpidw,zpiadj,zpmidep,                  &
                      zdrop,lin_aerosolmode)                  !aerosol feedback

!-----------------------------------------------------------------------
  use cc_mpi

  implicit none
!-----------------------------------------------------------------------
!  This program handles the vertical 1-D cloud micphysics
!-----------------------------------------------------------------------
! avisc: constant in empirical formular for dynamic viscosity of air
!         =1.49628e-6 [kg/m/s] = 1.49628e-5 [g/cm/s]
! adiffwv: constant in empirical formular for diffusivity of water
!          vapor in air
!          = 8.7602e-5 [kgm/s3] = 8.7602 [gcm/s3]
! axka: constant in empirical formular for thermal conductivity of air
!       = 1.4132e3 [m2/s2/K] = 1.4132e7 [cm2/s2/K]
! qi0: mixing ratio threshold for cloud ice aggregation [kg/kg]
! xmi50: mass of a 50 micron ice crystal
!        = 4.8e-10 [kg] =4.8e-7 [g]
! xmi40: mass of a 40 micron ice crystal
!        = 2.46e-10 [kg] = 2.46e-7 [g]
! di50: diameter of a 50 micro (radius) ice crystal
!       =1.0e-4 [m]
! xmi: mass of one cloud ice crystal
!      =4.19e-13 [kg] = 4.19e-10 [g]
! oxmi=1.0/xmi
!
! xni0=1.0e-2 [m-3] The value given in Lin et al. is wrong.(see
!                   Hsie et al.(1980) and Rutledge and Hobbs(1983) )
! bni=0.5 [K-1]
! xmnin: mass of a natural ice nucleus
!    = 1.05e-18 [kg] = 1.05e-15 [g] This values is suggested by
!                    Hsie et al. (1980)
!    = 1.0e-12 [kg] suggested by Rutlegde and Hobbs (1983)

! av_r: av_r in empirical formular for terminal
!         velocity of raindrop
!         =2115.0 [cm**(1-b)/s] = 2115.0*0.01**(1-b) [m**(1-b)/s]
! bv_r: bv_r in empirical formular for terminal
!         velocity of raindrop
!         =0.8
! av_i: av_i in empirical formular for terminal
!         velocity of snow
!         =152.93 [cm**(1-d)/s] = 152.93*0.01**(1-d) [m**(1-d)/s]
! bv_i: bv_i in empirical formular for terminal
!         velocity of snow
!         =0.25
! vf1r: ventilation factors for rain =0.78
! vf2r: ventilation factors for rain =0.31
! vf1s: ventilation factors for snow =0.65
! vf2s: ventilation factors for snow =0.44
!
!----------------------------------------------------------------------
! new 2D declaration
  integer,                         intent(in)    :: kts, kte
  integer,                         intent(in)    :: imax
  integer,                         intent(in)    :: lin_aerosolmode
  real,                            intent(in)    :: dt
  real, dimension(1:imax),         intent(in)    :: zsfc
  real, dimension(1:imax,kts:kte), intent(in)    :: zdrop, riz
  real, dimension(1:imax,kts:kte), intent(in)    :: tothz,rho,orho,sqrho,              &
                                                    prez,zz,dzw
  real, dimension(1:imax,kts:kte), intent(out)   :: precrz,preciz,precsz
  real, dimension(1:imax,kts:kte), intent(out)   :: EFFC1D,EFFI1D,                     &
                                                    EFFS1D,EFFR1D
  real, dimension(1:imax,kts:kte), intent(out)   :: fluxr,fluxi,fluxs,fluxg,           &
                                                    fluxm,fluxf,fevap,fsubl,           &
                                                    fauto,fcoll,faccr
  real, dimension(1:imax,kts:kte), intent(out)   :: vi,vs,vg
  real, dimension(1:imax,kts:kte), intent(out)   :: zpsnow,zpsaut,zpsfw,               &
                                                    zpsfi,zpraci,zpiacr,               &
                                                    zpsaci,zpsacw,zpsdep,              &
                                                    zpssub,zpracs,zpsacr,              &
                                                    zpsmlt,zpsmltevp,zprain,           &
                                                    zpraut,zpracw,zprevp,              &
                                                    zpgfr,zpvapor,zpclw,               &
                                                    zpladj,zpcli,zpimlt,               &
                                                    zpihom,zpidw,zpiadj,               &
                                                    zpmidep
  real, dimension(1:imax),         intent(inout) :: pptrain, pptsnow, pptice
  real, dimension(1:imax,kts:kte), intent(inout) :: qvz,qlz,qrz,qiz,qsz,thz
  real, dimension(1:imax,kts:kte), intent(inout) :: ncz,niz,nrz,nsz
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer                                        :: i, j, k, is, ie, iq, tile 
  ! local vars
  real                                           :: obp4, bp3, bp5, bp6, odp4,         &
                                                    dp3, dp5, dp5o2
  ! temperary vars
  real                                           :: tmp, tmp0, tmp1, tmp2,tmp3,        &
                                                    tmp4, tmpa,tmpb,tmpc,tmpd,alpha1,  &
                                                    qic, abi,abr, abg, odtberg,        &
                                                    vti50,eiw,eri,esi,esr, esw,        &
                                                    erw,delrs,term0,term1,             &
                                                    Ap, Bp,                            &
                                                    factor, tmp_r, tmp_s,tmp_g,        &
                                                    qlpqi, rsat, a1, a2, xnin
  real, dimension(1:imax)                        :: es1d
  real, dimension(1:imax,kts:kte)                :: nczodt, nizodt, nrzodt, nszodt
  real, dimension(1:imax,kts:kte)                :: oprez, tem, temcc, theiz, qswz,    &
                                                    qsiz, qvoqswz, qvoqsiz, qvzodt,    &
                                                    qlzodt, qizodt, qszodt, qrzodt
  real, dimension(1:imax,kts:kte)                :: tmp2d
  real, dimension(1:imax,kts:kte)                :: qvsbar,rs0,viscmu,visc,diffwv,     &
                                                    schmidt,xka
!--- microphysical processes
  real, dimension(1:imax,kts:kte)                :: psnow, psaut, psfw,  psfi,  praci, &
                                                    piacr, psaci, psacw, psdep, pssub, &
                                                    pracs, psacr, psmlt, psmltevp,     &
                                                    prain, praut, pracw, prevp, pvapor,&
                                                    pclw,  pladj, pcli,  pimlt, pihom, &
                                                    pidw,  piadj, pgfr,                &
                                                    qschg, pracis
!---- new snow parameters
  real                                          :: vf1s = 0.65,vf2s = 0.44,            &
                                                   vf1r =0.78,vf2r = 0.31 
  real                                          :: am_c1=0.004,am_c2= 6e-5,  am_c3=0.15
  real                                          :: bm_c1=1.85, bm_c2= 0.003, bm_c3=1.25
  real                                          :: aa_c1=1.28, aa_c2= -0.012,aa_c3=-0.6
  real                                          :: ba_c1=1.5,  ba_c2= 0.0075,ba_c3=0.5
  real                                          :: best_a=1.08 ,  best_b = 0.499
  real                                          :: disp, Dc_liu, eta, mu_c, R6c        !--- for Liu's autoconversion
  real, dimension(1:imax)                       :: tc0
  real, dimension(kts:kte)                      :: ab_s,ab_r,ab_riming 
  real, dimension(1:imax,kts:kte)               :: cap_s    !---- capacitance of snow
  real, dimension(1:imax,kts:kte)               :: am_s,bm_s,av_s,bv_s,Ri,tmp_ss,lams 
  real, dimension(1:imax,kts:kte)               :: aa_s,ba_s,tmp_sa 
  real                                          :: mu_s=0.,mu_i=0.,mu_r=0.
 
  ! Adding variable Riz, which will duplicate Ri but be a copy passed upward
  real                                          :: episp0k, dtb, odtb, pi, pio4,       &
                                                   pio6, oxLf, xLvocp, xLfocp, av_r,   &
                                                   av_i, ocdrag, gambp4, gamdp4,       &
                                                   gam4pt5, Cpor, oxmi, gambp3, gamdp3,&
                                                   gambp6, gam3pt5, gam2pt75, gambp5o2,&
                                                   gamdp5o2, cwoxlf, ocp, xni50, es
  real                                          :: qvmin=1.e-20
  real                                          :: temc1,save1,save2,xni50mx
  real, dimension(1:imax,kts:kte)               :: riz2d_temp
  real, dimension(1:imax,kts:kte)               :: vtr, vts,                           &
                                                   vtrold, vtsold, vtiold,             &
                                                   xlambdar, xlambdas,                 &
                                                   olambdar, olambdas
  ! for terminal velocity flux
  real                                          :: xmr,xms,xmc,dcs,xmr_i
  real                                          :: lamminr, lammaxr,lammins,           &
                                                   lammaxs,lammini, lammaxi,lammin,lammax
  real                                          :: gambvr1
  real                                          :: lvap
  real                                          :: mi0
  real, dimension(1:imax)                       :: max_ri
  real, dimension(1:imax)                       :: t_del_tv,del_tv,flux,               &
                                                   fluxin,fluxout,tmpqrz
  real, dimension(1:imax)                       :: nflux,nfluxin,nfluxout
  integer, dimension(1:imax)                    :: min_q, max_q, max_ri_k
  logical, dimension(1:imax)                    :: notlast, notlast_work
  real, dimension(1:imax,kts:kte)               :: npsaut, npraci, npiacr, npsaci,     &
                                                   npsacw, npssub, npsdep, npsacr,     &
                                                   npgfr,  npsmlt, npsmltevp,npraut,   &   
                                                   npracw, nprevp, nihom,  nimlt,      &
                                                   nsagg,  npraut_r
  real, dimension(1:imax,kts:kte)               :: nvtr,   nvts
  real, dimension(1:imax,kts:kte)               :: qisten, qrsten, qssten
  real, dimension(1:imax,kts:kte)               :: nisten, nrsten, nssten
  real, dimension(1:imax,kts:kte)               :: nidep, midep
  real, dimension(1:imax,kts:kte)               :: n0_r,n0_i,n0_c,n0_s                  
  real, dimension(1:imax,kts:kte)               :: lami,lamc
  !------------------------------------------------------------------------------------
  !call START_LOG(p1_begin)
  vtrold=0.
  vtsold=0.
  vtiold=0.

  mu_c    = AMIN1(15., (1000.E6/Nt_c + 2.))
  R6c     = 10.0E-6      !---- 10 micron, threshold radius of cloud droplet
  dtb     = dt                                                                         !sny
  odtb    =1./dtb
  pi      =acos(-1.)
  pio4    =acos(-1.)/4.
  pio6    =acos(-1.)/6.
  ocp     =1./cp
  oxLf    =1./xLf
  xLvocp  =xLv/cp
  xLfocp  =xLf/cp
  Cpor    =cp/Rair
  oxmi    =1.0/xmi
  cwoxlf  =cw/xlf 
  av_r    =2115.0*0.01**(1-bv_r)
  av_i    =152.93*0.01**(1-bv_i)
  ocdrag  =1./Cdrag
  episp0k =RH*ep2*1000.*svp1

  gambp4  =ggamma(bv_r+4.)
  gamdp4  =ggamma(bv_i+4.)
  gambp3  =ggamma(bv_r+3.)
  gambp6  =ggamma(bv_r+6)
  gambp5o2=ggamma((bv_r+5.)/2.)
  gamdp5o2=ggamma((bv_i+5.)/2.)
  gambvr1=ggamma(bv_r+1.)
  !
  !     oprez       1./prez ( prez : pressure)
  !     qsw         saturated mixing ratio on water surface
  !     qsi         saturated mixing ratio on ice surface
  !     episp0k     RH*e*saturated pressure at 273.15 K  = 611.2 hPa (Rogers and Yau 1989)
  !     qvoqsw      qv/qsw
  !     qvoqsi      qv/qsi
  !     qvzodt      qv/dt
  !     qlzodt      ql/dt
  !     qizodt      qi/dt
  !     qszodt      qs/dt

  obp4    =1.0/(bv_r+4.0)
  bp3     =bv_r+3.0
  bp5     =bv_r+5.0
  bp6     =bv_r+6.0
  odp4    =1.0/(bv_i+4.0)
  dp3     =bv_i+3.0
  dp5     =bv_i+5.0
  dp5o2   =0.5*(bv_i+5.0)
    
  dcs     = 125.E-6  ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
  xms     =pi*500.*(dcs)**3/6.   !morr =PI*RHOI*DCS**3/6.=5.11*10e-10
  xmr     =4./3.*pi*rhowater*(500.E-6)**3
  xmr_i     =4./3.*pi*rhowater*(25.E-6)**3
  mi0     = 4./3.*3.14*500.*(10.e-6)**3
  !    xmc     =4.17*10e-14 !4./3.*pi*(0.00001)**3*1000.
  lammaxr = 1./20.E-6
  !lamminr = 1./500.E-6
  lamminr = 1./2800.E-6
  lammaxs = 1./10.E-6
  lammins = 1./2000.E-6
  lammaxi = 1./1.E-6
  lammini = 1./(2.*dcs+100.E-6)

  fluxs = 0.                                        ! sny
  fluxi = 0.
  fluxr = 0.
  fluxs = 0.

  select case( lin_aerosolmode )
    case(0)
      do k=kts,kte
        !ncz(k) = 250.*1.E6/rho(k)
        ncz(1:imax,k) = Nt_c/rho(1:imax,k)
      end do
    case(1)
      do k=kts,kte
        ncz(1:imax,k) = zdrop(1:imax,k)/rho(1:imax,k)
      end do
    case default
      write(6,*) "ERROR: Unknown option aerosolmode"
      stop
  end select

  do k=kts,kte
    niz(1:imax,k) = min(niz(1:imax,k),0.3E6/rho(1:imax,k))
    nrz(1:imax,k) = amax1( 0.0,nrz(1:imax,k) )
    nsz(1:imax,k) = amax1( 0.0,nsz(1:imax,k) )
    nczodt(1:imax,k)=amax1( 0.0,odtb*ncz(1:imax,k) )
    nizodt(1:imax,k)=amax1( 0.0,odtb*niz(1:imax,k) )
    nrzodt(1:imax,k)=amax1( 0.0,odtb*nrz(1:imax,k) )
    nszodt(1:imax,k)=amax1( 0.0,odtb*nsz(1:imax,k) )
  end do

  do k=kts,kte
    oprez(1:imax,k)=1./prez(1:imax,k) 
    qlz(1:imax,k)  =amax1( 0.0,qlz(1:imax,k) )
    qiz(1:imax,k)  =amax1( 0.0,qiz(1:imax,k) )
    qvz(1:imax,k)  =amax1( qvmin,qvz(1:imax,k) )
    qsz(1:imax,k)  =amax1( 0.0,qsz(1:imax,k) )
    qrz(1:imax,k)  =amax1( 0.0,qrz(1:imax,k) )
    tem(1:imax,k)  =thz(1:imax,k)*tothz(1:imax,k)
    temcc(1:imax,k)=tem(1:imax,k)-273.15
    es1d(1:imax)     =1000.*svp1*exp( svp2*temcc(1:imax,k)/(tem(1:imax,k)-svp3) )  !--- RY89 Eq(2.17)
    qswz(1:imax,k) =ep2*es1d(1:imax)/(prez(1:imax,k)-es1d(1:imax))
 
    do iq=1,imax
      if (tem(iq,k) .lt. 273.15 ) then
        es1d(iq)=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
        qsiz(iq,k)=ep2*es1d(iq)/(prez(iq,k)-es1d(iq))
        if (temcc(iq,k) .lt. -40.0) qswz(iq,k)=qsiz(iq,k)
      else
        qsiz(iq,k)=qswz(iq,k)
      endif
    enddo

    qvoqswz(1:imax,k)  =qvz(1:imax,k)/qswz(1:imax,k)
    qvoqsiz(1:imax,k)  =qvz(1:imax,k)/qsiz(1:imax,k)
    qvzodt(1:imax,k)   =amax1( 0.0,odtb*qvz(1:imax,k) )
    qlzodt(1:imax,k)   =amax1( 0.0,odtb*qlz(1:imax,k) )
    qizodt(1:imax,k)   =amax1( 0.0,odtb*qiz(1:imax,k) )
    qszodt(1:imax,k)   =amax1( 0.0,odtb*qsz(1:imax,k) )
    qrzodt(1:imax,k)   =amax1( 0.0,odtb*qrz(1:imax,k) )
    theiz(1:imax,k)=thz(1:imax,k)+(xlvocp*qvz(1:imax,k)-xlfocp*qiz(1:imax,k))/tothz(1:imax,k)
  enddo

  do k=kts,kte
    psnow(1:imax,k)   =0.                  ! sum all process for snow
    psaut(1:imax,k)   =0.                  ! ice crystal aggregation to snow
    psfw(1:imax,k)    =0.                  ! BERGERON process to transfer cloud water to snow
    psfi(1:imax,k)    =0.                  ! BERGERON process to transfer cloud ice to snow
    praci(1:imax,k)   =0.                  ! cloud ice accretion by rain
    piacr(1:imax,k)   =0.                  ! rain accretion by cloud ice
    psaci(1:imax,k)   =0.                  ! ice crystal accretion by snow
    psacw(1:imax,k)   =0.                  ! accretion of cloud water by snow
    psdep(1:imax,k)   =0.                  ! deposition of snow
    pssub(1:imax,k)   =0.                  ! sublimation of snow (T<0)
    pracs(1:imax,k)   =0.                  ! accretion of snow by rain
    psacr(1:imax,k)   =0.                  ! accretion of rain by snow
    psmlt(1:imax,k)   =0.                  ! melting of snow
    psmltevp(1:imax,k)=0.                  ! evaporation of melting snow (T>0)
    prain(1:imax,k)   =0.                  ! sum all process for rain
    praut(1:imax,k)   =0.                  ! autoconversion of rain
    pracw(1:imax,k)   =0.                  ! accretion of cloud water by rain
    prevp(1:imax,k)   =0.                  ! evaporation of rain
    pgfr(1:imax,k)    =0.                  ! feezing of rain to form graupel (added to PI)
    pvapor(1:imax,k)  =0.                  ! sum all process for water vapor to determine qvz
    pclw(1:imax,k)    =0.                  ! sum all process for cloud liquid to determine qlz
    pladj(1:imax,k)   =0.                  ! saturation adjustment for ql
    pcli(1:imax,k)    =0.                  ! sum all process for cloud ice to determine qiz
    pimlt(1:imax,k)   =0.                  ! melting of ice crystal >0.
    pihom(1:imax,k)   =0.                  ! homogeneous nucleation <-40
    pidw(1:imax,k)    =0.                  ! production of cloud ice by BERGERON process
    piadj(1:imax,k)   =0.                  ! saturation adjustment for qi

    npsaut(1:imax,k)   =0.
    npraci(1:imax,k)   =0.
    npiacr(1:imax,k)   =0.
    npsaci(1:imax,k)   =0.
    npsacw(1:imax,k)   =0.
    npssub(1:imax,k)   =0.
    npsdep(1:imax,k)   =0.
    npsacr(1:imax,k)   =0.
    npgfr(1:imax,k)    =0.
    npsmlt(1:imax,k)   =0.
    npsmltevp(1:imax,k)=0.
    npraut(1:imax,k)   =0.
    npracw(1:imax,k)   =0.
    nprevp(1:imax,k)   =0.

    nimlt(1:imax,k)    =0.
    nihom(1:imax,k)    =0.
    nsagg(1:imax,k)    =0.
    npraut_r(1:imax,k) =0.

    n0_i(1:imax,k)     =0.
    n0_s(1:imax,k)     =0.
    n0_r(1:imax,k)     =0.
    n0_c(1:imax,k)     =0.
    lamc(1:imax,k)     =0.
    lami(1:imax,k)     =0.
    xlambdar(1:imax,k) =0.
    xlambdas(1:imax,k) =0.
    vtr(1:imax,k)      =0.
    vts(1:imax,k)      =0.
    vtiold(1:imax,k)   =0.
    nvtr(1:imax,k)     =0.
    nvts(1:imax,k)     =0.

    qisten(1:imax,k)   =0.
    qrsten(1:imax,k)   =0.
    qssten(1:imax,k)   =0.
    nisten(1:imax,k)   =0.
    nrsten(1:imax,k)   =0.
    nssten(1:imax,k)   =0.
    nidep(1:imax,k)    =0.
    midep(1:imax,k)    =0.
  end do

  !call END_LOG(p1_end)
  !call START_LOG(p2_begin)

  !***********************************************************************
  !*****  compute viscosity,difusivity,thermal conductivity, and    ******
  !*****  Schmidt number                                            ******
  !***********************************************************************
  !      viscmu: dynamic viscosity of air kg/m/s
  !      visc: kinematic viscosity of air = viscmu/rho (m2/s)
  !      avisc=1.49628e-6 kg/m/s=1.49628e-5 g/cm/s
  !      viscmu=1.718e-5 kg/m/s in RH
  !      diffwv: Diffusivity of water vapor in air
  !      adiffwv = 8.7602e-5 (8.794e-5 in MM5) kgm/s3
  !              = 8.7602 (8.794 in MM5) gcm/s3
  !      diffwv(k)=2.26e-5 m2/s
  !      schmidt: Schmidt number=visc/diffwv
  !      xka: thermal conductivity of air J/m/s/K (Kgm/s3/K)
  !      xka(k)=2.43e-2 J/m/s/K in RH
  !      axka=1.4132e3 (1.414e3 in MM5) m2/s2/k = 1.4132e7 cm2/s2/k
  !------------------------------------------------------------------
  do k=kts,kte
    viscmu(1:imax,k)=avisc*tem(1:imax,k)**1.5/(tem(1:imax,k)+120.0)
    visc(1:imax,k)=viscmu(1:imax,k)*orho(1:imax,k)
    diffwv(1:imax,k)=adiffwv*tem(1:imax,k)**1.81*oprez(1:imax,k)
    schmidt(1:imax,k)=visc(1:imax,k)/diffwv(1:imax,k)
    xka(1:imax,k)=axka*viscmu(1:imax,k)
    rs0(1:imax,k)=ep2*1000.*svp1/(prez(1:imax,k)-1000.*svp1)
  end do

! rewrite in "iq" by sny
!  do k=kts,kte
!    do iq=1,imax
!      viscmu(iq,k)=avisc*tem(iq,k)**1.5/(tem(iq,k)+120.0)
!      visc(iq,k)=viscmu(iq,k)*orho(iq,k)
!      diffwv(iq,k)=adiffwv*tem(iq,k)**1.81*oprez(iq,k)
!      schmidt(iq,k)=visc(iq,k)/diffwv(iq,k)
!      xka(iq,k)=axka*viscmu(iq,k)
!      rs0(iq,k)=ep2*1000.*svp1/(prez(iq,k)-1000.*svp1)
!    end do
!  end do

  ! ---- YLIN, set snow variables
  !
  !---- A+B in depositional growth, the first try just take from Rogers and Yau(1989)
  !         ab_s(k) = lsub*lsub*orv/(tcond(k)*temp(k))+&
  !                   rv*temp(k)/(diffu(k)*qvsi(k))

!  do k = kts, kte
!    tc0(1:imax)   = tem(1:imax,k)-273.15
!    do iq=1,imax
!      if (rho(iq,k)*qlz(iq,k) .gt. 1e-5 .AND. rho(iq,k)*qsz(iq,k) .gt. 1e-5) then
!        Ri(iq,k) = 1.0/(1.0+6e-5/(rho(iq,k)**1.170*qlz(iq,k)*qsz(iq,k)**0.170))
!      else
!        Ri(iq,k) = 0.
!      end if
!    end do
!  end do

! rewrite in "iq" by sny
  do k = kts, kte
    tc0(1:imax)   = tem(1:imax,k)-273.15
    where( rho(1:imax,k)*qlz(1:imax,k) .gt. 1e-5 .AND. rho(1:imax,k)*qsz(1:imax,k) .gt. 1e-5  )
      Ri(1:imax,k) = 1.0/(1.0+6e-5/(rho(1:imax,k)**1.170*qlz(1:imax,k)*qsz(1:imax,k)**0.170)) 
    elsewhere
      Ri(1:imax,k) = 0.
    end where
  end do


  !
  !--- make sure Ri does not decrease downward in a column
  !
  max_ri_k(:) = maxloc(Ri,dim=2)
  do iq = 1, imax
    max_ri(iq)   = Ri(iq,max_ri_k(iq))
  end do

  do k = kts,kte
    do iq = 1,imax
      if ( k<=max_ri_k(iq) ) then
        Ri(iq,k) = max_ri(iq)
      end if
    end do
  end do

  !--- YLIN, get PI properties
  do k = kts, kte
    Ri(1:imax,k) = AMAX1(0.,AMIN1(Ri(1:imax,k),1.0))
    ! Store the value of Ri(k) as Riz(k)
    riz2d_temp(1:imax,k) = Ri(1:imax,k)

    cap_s(1:imax,k)= 0.25*(1+Ri(1:imax,k))
    tc0(1:imax)    = AMIN1(-0.1, tem(1:imax,k)-273.15)
    n0_s(1:imax,k) = amin1(2.0E8, 2.0E6*exp(-0.12*tc0(1:imax)))
    am_s(1:imax,k) = am_c1+am_c2*tc0(1:imax)+am_c3*Ri(1:imax,k)*Ri(1:imax,k)   !--- Heymsfield 2007
    am_s(1:imax,k) = AMAX1(0.000023,am_s(1:imax,k))                                !--- use the a_min in table 1 of Heymsfield
    bm_s(1:imax,k) = bm_c1+bm_c2*tc0(1:imax)+bm_c3*Ri(1:imax,k)
    bm_s(1:imax,k) = AMIN1(bm_s(1:imax,k),3.0)                                     !---- capped by 3
    !--  converting from cgs to SI unit
    am_s(1:imax,k) =  10**(2*bm_s(1:imax,k)-3.0)*am_s(1:imax,k)
    aa_s(1:imax,k) = aa_c1 + aa_c2*tc0(1:imax) + aa_c3*Ri(1:imax,k)
    ba_s(1:imax,k) = ba_c1 + ba_c2*tc0(1:imax) + ba_c3*Ri(1:imax,k)
    !--  convert to SI unit as in paper
    aa_s(1:imax,k) = (1e-2)**(2.0-ba_s(1:imax,k))*aa_s(1:imax,k)
    !---- get v from Mitchell 1996
    av_s(1:imax,k) = best_a*viscmu(1:imax,k)*(2*grav*am_s(1:imax,k)/rho(1:imax,k)/ &
                       aa_s(1:imax,k)/(viscmu(1:imax,k)**2))**best_b
    bv_s(1:imax,k) = best_b*(bm_s(1:imax,k)-ba_s(1:imax,k)+2)-1

    tmp_ss(1:imax,k)= bm_s(1:imax,k)+mu_s+1
    tmp_sa(1:imax,k)= ba_s(1:imax,k)+mu_s+1
  end do

! rewrite in "iq" by sny
!  !--- YLIN, get PI properties
!  do k = kts, kte
!    do iq = 1, imax
!      Ri(iq,k) = AMAX1(0.,AMIN1(Ri(iq,k),1.0))
!      ! Store the value of Ri(k) as Riz(k)
!      riz2d_temp(iq,k) = Ri(iq,k)
!
!      cap_s(iq,k)= 0.25*(1+Ri(iq,k))
!      tc0(iq)    = AMIN1(-0.1, tem(iq,k)-273.15)
!      n0_s(iq,k) = amin1(2.0E8, 2.0E6*exp(-0.12*tc0(iq)))
!      am_s(iq,k) = am_c1+am_c2*tc0(iq)+am_c3*Ri(iq,k)*Ri(iq,k)   !--- Heymsfield 2007
!      am_s(iq,k) = AMAX1(0.000023,am_s(iq,k))                                !--- use the a_min in table 1 of Heymsfield
!      bm_s(iq,k) = bm_c1+bm_c2*tc0(iq)+bm_c3*Ri(iq,k)
!      bm_s(iq,k) = AMIN1(bm_s(iq,k),3.0)                                     !---- capped by 3
!      !--  converting from cgs to SI unit
!      am_s(iq,k) =  10**(2*bm_s(iq,k)-3.0)*am_s(iq,k)
!      aa_s(iq,k) = aa_c1 + aa_c2*tc0(iq) + aa_c3*Ri(iq,k)
!      ba_s(iq,k) = ba_c1 + ba_c2*tc0(iq) + ba_c3*Ri(iq,k)
!      !--  convert to SI unit as in paper
!      aa_s(iq,k) = (1e-2)**(2.0-ba_s(iq,k))*aa_s(iq,k)
!      !---- get v from Mitchell 1996
!      av_s(iq,k) = best_a*viscmu(iq,k)*(2*grav*am_s(iq,k)/rho(iq,k)/ &
!                       aa_s(iq,k)/(viscmu(iq,k)**2))**best_b
!      bv_s(iq,k) = best_b*(bm_s(iq,k)-ba_s(iq,k)+2)-1
!      tmp_ss(iq,k)= bm_s(iq,k)+mu_s+1
!      tmp_sa(iq,k)= ba_s(iq,k)+mu_s+1
!    end do
!  end do


  !call END_LOG(p2_end)
  !call START_LOG(p3_begin)

  !***********************************************************************
  ! Calculate precipitation fluxes due to terminal velocities.
  !***********************************************************************
  !
  !- Calculate termianl velocity (vt?)  of precipitation q?z
  !- Find maximum vt? to determine the small delta t
  !
  !-- rain
  !       CALL wrf_debug ( 100 , 'module_ylin, start precip fluxes' )

  t_del_tv(1:imax)=0.
  del_tv(1:imax)=dtb
  notlast(1:imax)=.true.
  DO while (any(notlast))

    notlast_work = notlast

    min_q(1:imax)=kte
    max_q(1:imax)=kts-1

    ! if no rain, --> minq>maxq --> notlast=False (only case minq>maxq)
    ! if rain --> minq<maxq (always), some vertical points norain--> no lamda, velocity
    do k=kts,kte-1
      do iq = 1,imax
        if (notlast(iq)) then
          if (qrz(iq,k) .gt. 1.0e-8) then
            min_q(iq)=min0(min_q(iq),k)
            max_q(iq)=max0(max_q(iq),k)
            ! tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrz(k))
            ! tmp1=sqrt(tmp1)
            ! vtrold(k)=o6*av_r*gambp4*sqrho(k)/tmp1**bv_r
            xlambdar(iq,k)=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
            n0_r(iq,k)=nrz(iq,k)*xlambdar(iq,k)
            if (xlambdar(iq,k).lt.lamminr) then
              xlambdar(iq,k) = lamminr
              n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k) 
            else if (xlambdar(iq,k).gt.lammaxr) then
              xlambdar(iq,k) = lammaxr
              n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k)
            end if
            olambdar(iq,k)=1.0/xlambdar(iq,k)
            vtrold(iq,k)=o6*av_r*gambp4*sqrho(iq,k)*olambdar(iq,k)**bv_r
            nvtr(iq,k)=av_r*gambvr1*sqrho(iq,k)*olambdar(iq,k)**bv_r
            
            if (k .eq. 1) then
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zsfc(iq))/vtrold(iq,k))
            else
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zz(iq,k-1))/vtrold(iq,k))
            endif
          else
            vtrold(iq,k)=0.
            nvtr(iq,k)=0.
            olambdar(iq,k)=0.
          endif
        endif     ! notlast
      enddo       ! iq
    enddo         ! k

    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)

    do iq = 1,imax
      if (notlast(iq)) then
        if (max_q(iq) .ge. min_q(iq)) then
          t_del_tv(iq)=t_del_tv(iq)+del_tv(iq)
          if ( t_del_tv(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv(iq)=dtb+del_tv(iq)-t_del_tv(iq)
          end if
          fluxin(iq)=0.
          nfluxin(iq)=0. ! sny
        end if ! maxq>minq
      end if   ! notlast
    end do     ! iq

    do k = kte,kts,-1
    !do k = maxval(max_q),minval(min_q),-1
      do iq = 1,imax
        if ( notlast(iq) ) then
          !if (max_q(iq) .ge. min_q(iq)) then      
          if ( k>=min_q(iq) .and. k<=max_q(iq) ) then
            fluxout(iq)=rho(iq,k)*vtrold(iq,k)*qrz(iq,k)
            flux(iq)=(fluxin(iq)-fluxout(iq))/rho(iq,k)/dzw(iq,k)
            tmpqrz(iq)=qrz(iq,k)
            qrz(iq,k)=qrz(iq,k)+del_tv(iq)*flux(iq)
            fluxin(iq)=fluxout(iq)

            nfluxout(iq)=rho(iq,k)*nvtr(iq,k)*nrz(iq,k)
            nflux(iq)=(nfluxin(iq)-nfluxout(iq))/rho(iq,k)/dzw(iq,k)
            nrz(iq,k)=nrz(iq,k)+del_tv(iq)*nflux(iq)
            nfluxin(iq)=nfluxout(iq)
            qrsten(iq,k)=flux(iq)
            nrsten(iq,k)=nflux(iq)
 
            fluxr(iq,k) = fluxout(iq)                     ! sny
          end if !maxq, minq
        end if   !notlast
      end do     !iq
    end do       !k

    do iq = 1, imax
      if ( notlast(iq) ) then
        if (max_q(iq) .ge. min_q(iq)) then      
          if (min_q(iq) .eq. 1) then
            pptrain(iq)=pptrain(iq)+fluxin(iq)*del_tv(iq)
          else
            qrz(iq,min_q(iq)-1)=qrz(iq,min_q(iq)-1)+del_tv(iq)*  &
                           fluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
            nrz(iq,min_q(iq)-1)=nrz(iq,min_q(iq)-1)+del_tv(iq)*  &
                           nfluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
          endif  !minq
        else
          notlast_work(iq)=.false.
        end if 
      end if   ! notlast
    end do     !iq


    notlast(:) = notlast_work(:)
  ENDDO      ! while(any(notlast))


  !
  !-- snow
  !
  t_del_tv(1:imax)=0.
  del_tv(1:imax)=dtb
  notlast(1:imax)=.true.
  DO while (any(notlast))

    notlast_work = notlast
    min_q(1:imax)=kte
    max_q(1:imax)=kts-1

    do k=kts,kte-1
      do iq = 1, imax
        if (notlast(iq)) then
          if (qsz(iq,k) .gt. 1.0e-8) then
            min_q(iq)=min0(min_q(iq),k)
            max_q(iq)=max0(max_q(iq),k)
            !   tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
            !   **(1./tmp_ss(k))
            !   vtsold(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
            !   ggamma(tmp_ss(k))/(tmp1**bv_s(k))
            ! Zhao 2022 - Row 2 Table 2 or Lin 2011 - Formula A3
            xlambdas(iq,k)=(am_s(iq,k)*ggamma(tmp_ss(iq,k))*nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
            ! Zhao 2022 - Row 1 Table 2
            n0_s(iq,k)=nsz(iq,k)*xlambdas(iq,k)
            if (xlambdas(iq,k).lt.lammins) then
              xlambdas(iq,k)= lammins
              n0_s(iq,k) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,k)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,k)/xlambdas(iq,k)
            else if (xlambdas(iq,k).gt.lammaxs) then
              xlambdas(iq,k) = lammaxs
              n0_s(iq,k) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,k)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,k)/xlambdas(iq,k)
            end if
            olambdas(iq,k)=1.0/xlambdas(iq,k)
            ! Zhao 2022 - Row 3 Table 2
            vtsold(iq,k)= sqrho(iq,k)*av_s(iq,k)*ggamma(bv_s(iq,k)+tmp_ss(iq,k))/ &
               ggamma(tmp_ss(iq,k))*(olambdas(iq,k)**bv_s(iq,k))
            ! Zhao 2022 - Row 4 Table 2
            nvts(iq,k)=sqrho(iq,k)*av_s(iq,k)*ggamma(bv_s(iq,k)+1)*(olambdas(iq,k)**bv_s(iq,k))
            if (k .eq. 1) then
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zsfc(iq))/vtsold(iq,k))
            else
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zz(iq,k-1))/vtsold(iq,k))
            endif
          else
            vtsold(iq,k)=0.
            nvts(iq,k)=0.
            olambdas(iq,k)=0.
          endif
        endif   ! notlast
      enddo     ! iq
    enddo       ! k

    !
    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)

    do iq = 1,imax
      if (notlast(iq)) then
        if (max_q(iq) .ge. min_q(iq)) then
          t_del_tv(iq)=t_del_tv(iq)+del_tv(iq)
          if ( t_del_tv(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv(iq)=dtb+del_tv(iq)-t_del_tv(iq)
          endif
          fluxin(iq) = 0.
          nfluxin(iq) = 0.
        end if ! maxq>minq
      end if   ! notlast
    end do     ! iq

    do k = kte,kts,-1
      do iq = 1,imax
        if ( notlast(iq) ) then
          if ( k>=min_q(iq) .and. k<=max_q(iq) ) then
            fluxout(iq)=rho(iq,k)*vtsold(iq,k)*qsz(iq,k)
            flux(iq)=(fluxin(iq)-fluxout(iq))/rho(iq,k)/dzw(iq,k)
            qsz(iq,k)=qsz(iq,k)+del_tv(iq)*flux(iq)
            qsz(iq,k)=amax1(0.,qsz(iq,k))
            fluxin(iq)=fluxout(iq)

            nfluxout(iq)=rho(iq,k)*nvts(iq,k)*nsz(iq,k)
            nflux(iq)   =(nfluxin(iq)-nfluxout(iq))/rho(iq,k)/dzw(iq,k)
            nsz(iq,k)  =nsz(iq,k)+del_tv(iq)*nflux(iq)
            nfluxin(iq) =nfluxout(iq)
            qssten(iq,k)=flux(iq)
            nssten(iq,k)=nflux(iq)

            fluxs(iq,k) = fluxout(iq)                     ! sny
          end if ! maxq, minq
        end if   ! notlast
      end do     ! iq
    end do       ! k

    do iq = 1, imax
      if ( notlast(iq) ) then
        if (max_q(iq) .ge. min_q(iq)) then
          if (min_q(iq) .eq. 1) then
            pptsnow(iq)=pptsnow(iq)+fluxin(iq)*del_tv(iq)
          else
            qsz(iq,min_q(iq)-1)=qsz(iq,min_q(iq)-1)+del_tv(iq)*  &
                     fluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
            nsz(iq,min_q(iq)-1)=nsz(iq,min_q(iq)-1)+del_tv(iq)*  &
                       nfluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
          endif ! minq
        else
          notlast_work(iq)=.false.
        endif
      endif     ! notlast
    enddo       ! iq

    notlast(:) = notlast_work(:)
  ENDDO       ! while(any(notlast))

  !
  !-- cloud ice  (03/21/02) using Heymsfield and Donner (1990) Vi=3.29*qi^0.16
  !
  t_del_tv(1:imax)=0.
  del_tv(1:imax)=dtb
  notlast(1:imax)=.true.
  DO while (any(notlast))

    notlast_work = notlast
    min_q(1:imax)=kte
    max_q(1:imax)=kts-1

    do k=kts,kte-1
      do iq = 1, imax
        if (notlast(iq)) then
          if (qiz(iq,k) .gt. 1.0e-8) then
            min_q(iq)=min0(min_q(iq),k)
            max_q(iq)=max0(max_q(iq),k)
            vtiold(iq,k)= 3.29 * (rho(iq,k)* qiz(iq,k))** 0.16  ! Heymsfield and Donner
            if (k .eq. 1) then
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zsfc(iq))/vtiold(iq,k))
            else
              del_tv(iq)=amin1(del_tv(iq),0.9*(zz(iq,k)-zz(iq,k-1))/vtiold(iq,k))
            endif
          else
            vtiold(iq,k)=0.
          endif
        endif   ! notlast
      enddo     ! iq
    enddo       ! k

    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)
    do iq = 1,imax
      if (notlast(iq)) then
        if (max_q(iq) .ge. min_q(iq)) then
          t_del_tv(iq)=t_del_tv(iq)+del_tv(iq)
          if ( t_del_tv(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv(iq)=dtb+del_tv(iq)-t_del_tv(iq)
          endif
          fluxin(iq) = 0.
          nfluxin(iq) = 0.
        end if   ! maxq>minq
      end if     ! notlast
    end do       ! iq

    do k = kte,kts,-1
      do iq = 1,imax
        if ( notlast(iq) ) then
          if ( k>=min_q(iq) .and. k<=max_q(iq) ) then
            fluxout(iq)=rho(iq,k)*vtiold(iq,k)*qiz(iq,k)
            flux(iq)=(fluxin(iq)-fluxout(iq))/rho(iq,k)/dzw(iq,k)
            qiz(iq,k)=qiz(iq,k)+del_tv(iq)*flux(iq)
            qiz(iq,k)=amax1(0.,qiz(iq,k))
            fluxin(iq)=fluxout(iq)

            nfluxout(iq)=rho(iq,k)*vtiold(iq,k)*niz(iq,k)
            nflux(iq)=(nfluxin(iq)-nfluxout(iq))/rho(iq,k)/dzw(iq,k)
            niz(iq,k)=niz(iq,k)+del_tv(iq)*nflux(iq)
            niz(iq,k)=amax1(0.,niz(iq,k))
            nfluxin(iq)=nfluxout(iq)
            qisten(iq,k)=flux(iq)
            nisten(iq,k)=nflux(iq)

            fluxi(iq,k) = fluxout(iq)                     ! sny
          end if ! maxq, minq
        end if   ! notlast
      end do     ! iq
    end do       ! k

    do iq = 1, imax
      if ( notlast(iq) ) then
        if (max_q(iq) .ge. min_q(iq)) then
          if (min_q(iq) .eq. 1) then
            pptice(iq)=pptice(iq)+fluxin(iq)*del_tv(iq)
          else
            qiz(iq,min_q(iq)-1)=qiz(iq,min_q(iq)-1)+del_tv(iq)*  &
                         fluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
            niz(iq,min_q(iq)-1)=niz(iq,min_q(iq)-1)+del_tv(iq)*  &
                         nfluxin(iq)/rho(iq,min_q(iq)-1)/dzw(iq,min_q(iq)-1)
          endif
        else
          notlast_work(iq)=.false.
        end if
      end if      ! notlast
    end do        ! iq
    
    notlast(:) = notlast_work(:)
  ENDDO       ! while(any(notlast))

  ! zdc 20220116
  do k=kts,kte-1                         !sg beg
    precrz(1:imax,k)=rho(1:imax,k)*vtrold(1:imax,k)*qrz(1:imax,k)
    preciz(1:imax,k)=rho(1:imax,k)*vtiold(1:imax,k)*qiz(1:imax,k)
    precsz(1:imax,k)=rho(1:imax,k)*vtsold(1:imax,k)*qsz(1:imax,k)
  enddo                                  !sg end
    precrz(1:imax,kte)=0. !wig - top level never set for vtXold vars
    preciz(1:imax,kte)=0. !wig
    precsz(1:imax,kte)=0. !wig
  !     CALL wrf_debug ( 100 , 'module_ylin: end precip flux' )

  !call END_LOG(p3_end)
  !call START_LOG(p4_begin)

  ! Microphpysics processes
  DO k=kts,kte
    qvzodt(1:imax,k)=amax1( 0.0,odtb*qvz(1:imax,k) )
    qlzodt(1:imax,k)=amax1( 0.0,odtb*qlz(1:imax,k) )
    qizodt(1:imax,k)=amax1( 0.0,odtb*qiz(1:imax,k) )
    qszodt(1:imax,k)=amax1( 0.0,odtb*qsz(1:imax,k) )
    qrzodt(1:imax,k)=amax1( 0.0,odtb*qrz(1:imax,k) )

    !***********************************************************************
    !*****   diagnose mixing ratios (qrz,qsz), terminal                *****
    !*****   velocities (vtr,vts), and slope parameters in size        *****
    !*****   distribution(xlambdar,xlambdas) of rain and snow          *****
    !*****   follows Nagata and Ogura, 1991, MWR, 1309-1337. Eq (A7)   *****
    !***********************************************************************
    !
    !**** assuming no cloud water can exist in the top two levels due to
    !**** radiation consideration
    !
    !!  if
    !!     unsaturated,
    !!     no cloud water, rain, ice, snow
    !!  then
    !!     skip these processes and jump to line 2000
    !
    !
    tmp2d(1:imax,k)=qiz(1:imax,k)+qlz(1:imax,k)+qsz(1:imax,k)+qrz(1:imax,k)
       
    do iq = 1, imax
      if( .not.(qvz(iq,k)+qlz(iq,k)+qiz(iq,k) .lt. qsiz(iq,k)  &
            .and. tmp2d(iq,k) .eq. 0.0) ) then !go to 2000
      !
      !! calculate terminal velocity of rain
      !
        if (qrz(iq,k) .gt. 1.0e-8) then
        !  tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
        !  xlambdar(k)=sqrt(tmp1)
        !  olambdar(k)=1.0/xlambdar(k)
        !  vtrold(k)=o6*av_r*gambp4*sqrho(k)*olambdar(k)**bv_r
          xlambdar(iq,k)=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
          n0_r(iq,k)=nrz(iq,k)*xlambdar(iq,k)
          if (xlambdar(iq,k).lt.lamminr) then
            xlambdar(iq,k) = lamminr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k)
          else if (xlambdar(iq,k).gt.lammaxr) then
            xlambdar(iq,k) = lammaxr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k)
          end if
          olambdar(iq,k)=1.0/xlambdar(iq,k)
          vtrold(iq,k)=o6*av_r*gambp4*sqrho(iq,k)*olambdar(iq,k)**bv_r
          nvtr(iq,k)=av_r*gambvr1*sqrho(iq,k)*olambdar(iq,k)**bv_r
        else
          vtrold(iq,k)=0.
          olambdar(iq,k)=0.
          nvtr(iq,k)=0.
        end if  ! qrz

        if (qrz(iq,k) .gt. 1.0e-8) then
        !  tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
        !  xlambdar(k)=sqrt(tmp1)
        !  olambdar(k)=1.0/xlambdar(k)
        !  vtr(k)=o6*av_r*gambp4*sqrho(k)*olambdar(k)**bv_r
        else
        !  vtr(k)=0.
        !  olambdar(k)=0.
        endif
        vtr(iq,k)=vtrold(iq,k)


        !!
        !!! calculate terminal velocity of snow
        !!
        if (qsz(iq,k) .gt. 1.0e-8) then
        !  tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
        !                   **(1./tmp_ss(k))
        !            xlambdas(k)=tmp1
        !            olambdas(k)=1.0/tmp1
        !            vtsold(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
        !                      ggamma(tmp_ss(k))/(tmp1**bv_s(k))
          xlambdas(iq,k)=(am_s(iq,k)*ggamma(tmp_ss(iq,k))*nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
          n0_s(iq,k)=nsz(iq,k)*xlambdas(iq,k)
          if (xlambdas(iq,k).lt.lammins) then
            xlambdas(iq,k)= lamminS
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,K)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          else if (xlambdas(iq,k).gt.lammaxs) then
            xlambdas(iq,k) = lammaxs
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,K)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          end if
          olambdas(iq,k)=1.0/xlambdas(iq,k)
          vtsold(iq,k)= sqrho(iq,k)*av_s(iq,k)*ggamma(bv_s(iq,k)+tmp_ss(iq,k))/ &
                 ggamma(tmp_ss(iq,k))*(olambdas(iq,k)**bv_s(iq,k))
          nvts(iq,k)=sqrho(iq,k)*av_s(iq,k)*ggamma(bv_s(iq,k)+1)*(olambdas(iq,k)**bv_s(iq,k))
        else
          vtsold(iq,k)=0.
          olambdas(iq,k)=0.
          xlambdas(iq,k)=0.
          nvts(iq,k)=0.
        endif

        if (qsz(iq,k) .gt. 1.0e-8) then
        ! tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
        !                   **(1./tmp_ss(k))
        !             olambdas(k)=1.0/tmp1
        !             vts(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
        !                      ggamma(tmp_ss(k))/(tmp1**bv_s(k))
        else
        !            vts(iq,k)=0.
        !            olambdas(iq,k)=0.
        endif
        vts(iq,k)=vtsold(iq,k)

        !---------- start of snow/ice processes below freezing

        if (tem(iq,k) .lt. 273.15) then

        !
        ! ice nucleation, cooper curve

          if ((qvoqswz(iq,k).ge.0.999.and.temcc(iq,k).le. -8.).or. &
            qvoqsiz(iq,k).ge.1.08) then
            nidep(iq,k) = 5.*exp(0.304*(273.15-tem(iq,k)))     ! m-3
            nidep(iq,k) = min(nidep(iq,k), 500.e3)             !5.e8) sny ! limit to 500 L-1
            nidep(iq,k) = max(nidep(iq,k)/rho(iq,k), 0.)       ! convert to kg-1
            nidep(iq,k) = (nidep(iq,k) - niz(iq,k))*odtb
            midep(iq,k) = nidep(iq,k)*mi0

            !**********PROCESS**********
            zpmidep(iq,k)= midep(iq,k)   ! wapor deposition to ice
            !---------------------------            
          end if
          !***********************************************************************
          !*********        snow production processes for T < 0 C       **********
          !***********************************************************************
          !
          ! (1) ICE CRYSTAL AGGREGATION TO SNOW (Psaut): Lin (21)
          !!    psaut=alpha1*(qi-qi0)
          !!    alpha1=1.0e-3*exp(0.025*(T-T0))
          !
          alpha1=1.0e-3*exp( 0.025*temcc(iq,k) )

          if(temcc(iq,k) .lt. -20.0) then
              tmp1=-7.6+4.0*exp( -0.2443e-3*MAX(-temcc(iq,k)-20,0.)**2.455 )
              qic=1.0e-3*exp(tmp1)*orho(iq,k)
          else
              qic=qi0
          end if

          tmp1=odtb*(qiz(iq,k)-qic)*(1.0-exp(-alpha1*dtb))
          psaut(iq,k)=amax1( 0.0,tmp1 )
          npsaut(iq,k)=amax1( 0.0,psaut(iq,k)/xms)

          !**********PROCESS**********
          zpsaut(iq,k) = psaut(iq,k)  ! ice crystal aggregation to snow
          !---------------------------

          !
          ! (2) BERGERON PROCESS TRANSFER OF CLOUD WATER TO SNOW (Psfw)
          !     this process only considered when -31 C < T < 0 C
          !     Lin (33) and Hsie (17)
          !!
          !!    parama1 and parama2 functions must be user supplied
          !!

          if( qlz(iq,k) .gt. 1.0e-10 ) then
            temc1=amax1(-30.99,temcc(iq,k))
            a1=parama1( temc1 )
            a2=parama2( temc1 )
            tmp1=1.0-a2
            !!   change unit from cgs to mks
            a1=a1*0.001**tmp1
            !!   dtberg is the time needed for a crystal to grow from 40 to 50 um
            !!   odtberg=1.0/dtberg
            odtberg=(a1*tmp1)/(xmi50**tmp1-xmi40**tmp1)
            !
            !!   compute terminal velocity of a 50 micron ice cystal
            !
            vti50=av_i*di50**bv_i*sqrho(iq,k)
            eiw=1.0
            save1=a1*xmi50**a2
            save2=0.25*pi*eiw*rho(iq,k)*di50*di50*vti50
            tmp2=( save1 + save2*qlz(iq,k) )
            !
            !!  maximum number of 50 micron crystals limited by the amount
            !!  of supercool water
            !
            xni50mx=qlzodt(iq,k)/tmp2
            !
            !!   number of 50 micron crystals produced
            !
            xni50=qiz(iq,k)*( 1.0-exp(-dtb*odtberg) )/xmi50
            xni50=amin1(xni50,xni50mx)
            !
            tmp3=odtb*tmp2/save2*( 1.0-exp(-save2*xni50*dtb) )
            psfw(iq,k)=amin1( tmp3,qlzodt(iq,k) )

            !**********PROCESS**********       
            zpsfw(iq,k) = psfw(iq,k)    ! BERGERON process to transfer cloud water to snow
            !---------------------------
            !
            ! (3) REDUCTION OF CLOUD ICE BY BERGERON PROCESS (Psfi): Lin (34)
            !     this process only considered when -31 C < T < 0 C
            !
            tmp1=xni50*xmi50-psfw(iq,k)
            psfi(iq,k)=amin1(tmp1,qizodt(iq,k))

            !**********PROCESS**********
            zpsfi(iq,k) = psfi(iq,k)    ! BERGERON process to transfer cloud ice to snow
            !---------------------------
          end if
        
          if(qrz(iq,k) .gt. 0.0) then  ! go to 1000

          !
          ! Processes (4) and (5) only need when qrz > 0.0
          !
          ! (4) CLOUD ICE ACCRETION BY RAIN (Praci): Lin (25)
          !     produce PI
          !
            eri=1.0
            save1=pio4*eri*n0_r(iq,k)*av_r*sqrho(iq,k)
            tmp1=save1*gambp3*olambdar(iq,k)**bp3
            praci(iq,k)=qizodt(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npraci(iq,k)=niz(iq,k)*tmp1

            !**********PROCESS**********
            zpraci(iq,k) = praci(iq,k)  ! cloud ice accretion by rain
            !---------------------------

          !
          ! (5) RAIN ACCRETION BY CLOUD ICE (Piacr): Lin (26)
          !
            tmp2=qiz(iq,k)*save1*rho(iq,k)*pio6*rhowater*gambp6*oxmi* &
                olambdar(iq,k)**bp6
            piacr(iq,k)=amin1( tmp2,qrzodt(iq,k) )
            npiacr(iq,k)=pio4*eri*nrz(iq,k)*av_r*niz(iq,k)*gambp3* &
                olambdar(iq,k)**bp3  !--wdm6

            !**********PROCESS**********
            zpiacr(iq,k) = piacr(iq,k)  ! rain accretion by cloud ice
            !---------------------------
          end if !1000    continue

          if(qsz(iq,k) .gt. 0.0) then !go to 1200
          !
          ! Compute the following processes only when qsz > 0.0
          !
          !
          ! (6) ICE CRYSTAL ACCRETION BY SNOW (Psaci): Lin (22)
          !
            esi=exp( 0.025*temcc(iq,k) )
            save1 = aa_s(iq,k)*sqrho(iq,k)*n0_s(iq,k)* &
            ggamma(bv_s(iq,k)+tmp_sa(iq,k))*         &
            olambdas(iq,k)**(bv_s(iq,k)+tmp_sa(iq,k))

            tmp1=esi*save1
            psaci(iq,k)=qizodt(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsaci(iq,k)=amin1( tmp1*niz(iq,k),nizodt(iq,k))

            !**********PROCESS**********
            zpsaci(iq,k) = psaci(iq,k)  ! ! ice crystal accretion by snow
            !---------------------------

          !
          ! (7) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            tmp1=esw*save1
            psacw(iq,k)=qlzodt(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsacw(iq,k)=amin1(tmp1*ncz(iq,k),ncz(iq,k))

            !**********PROCESS**********
            zpsacw(iq,k) = psacw(iq,k)  ! ! accretion of cloud water by snow
            !---------------------------

          ! recalculate the saturatuin temperature
          !
          ! (8) DEPOSITION/SUBLIMATION OF SNOW (Psdep/Pssub): Lin (31)
          !     includes consideration of ventilation effect
          !
            tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
            tmpb=xls*xls*rho(iq,k)*qsiz(iq,k)*diffwv(iq,k)
            tmpc=tmpa*qsiz(iq,k)*diffwv(iq,k)
            abi=4.0*pi*cap_s(iq,k)*(qvoqsiz(iq,k)-1.0)*tmpc/(tmpa+tmpb)
            tmp1=av_s(iq,k)*sqrho(iq,k)*        &
                olambdas(iq,k)**(5+bv_s(iq,k)+2*mu_s)/visc(iq,k)

          !---- YLIN, here there is some approximation assuming mu_s =1, so gamma(2)=1, etc.

            tmp2= abi*n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                 vf2s*schmidt(iq,k)**0.33334* &
                 ggamma(2.5+0.5*bv_s(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=odtb*( qvz(iq,k)-qsiz(iq,k) )
            tmp3=amin1(tmp3,0.)

            if( tmp2 .le. 0.0) then
              tmp2=amax1( tmp2,tmp3)
              pssub(iq,k)=amax1( tmp2,-qszodt(iq,k) )
              psdep(iq,k)=0.0
 
              !**********PROCESS**********
              zpssub(iq,k) = pssub(iq,k)  ! snow submimation to form water vapor
              zpsdep(iq,k) = psdep(iq,k)  ! water wapor depostion to from snow
              !---------------------------
            else
              psdep(iq,k)=amin1( tmp2,tmp3 )
              pssub(iq,k)=0.0
              !**********PROCESS**********
              zpssub(iq,k) = pssub(iq,k)  ! snow submimation to form water vapor
              zpsdep(iq,k) = psdep(iq,k)  ! water wapor depostion to from snow
              !---------------------------
            end if
            if(qsz(iq,k) .ge. 0.0) then
              npssub(iq,k)=pssub(iq,k)*nsz(iq,k)/qsz(iq,k)
              npsdep(iq,k)=npsdep(iq,k)*nsz(iq,k)/qsz(iq,k)
            else
              npssub(iq,k)=pssub(iq,k)/xms
              npsdep(iq,k)=npsdep(iq,k)/xms
            end if

            if(qrz(iq,k) .gt. 0.0) then !go to 1200
            !
            ! Compute processes (9) and (10) only when qsz > 0.0 and qrz > 0.0
            ! these two terms need to be refined in the future, they should be equal
            !
            ! (9) ACCRETION OF SNOW BY RAIN (Pracs): Lin (27)
            !
              esr=1.0
              tmpa=olambdar(iq,k)*olambdar(iq,k)
              tmpb=olambdas(iq,k)*olambdas(iq,k)
              tmpc=olambdar(iq,k)*olambdas(iq,k)
              tmp1=pi*pi*esr*n0_r(iq,k)*n0_s(iq,k)*    &
                        abs( vtr(iq,k)-vts(iq,k) )*orho(iq,k)
            ! tmp1=pi*pi*esr*n0_r(k)*N0_s(k)*            &
            ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5*orho(k)
              tmp2=tmpb*tmpb*olambdar(iq,k)*(5.0*tmpb+2.0*tmpc+0.5*tmpa)
              tmp3=tmp1*rhosnow*tmp2
              pracs(iq,k)=amin1( tmp3,qszodt(iq,k) )

              !**********PROCESS**********
              zpracs(iq,k) = pracs(iq,k)  ! accretion of snow by rain
              !---------------------------
            !
            ! (10) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
            !
              tmp3=tmpa*tmpa*olambdas(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
              tmp4=tmp1*rhowater*tmp3
              psacr(iq,k)=amin1( tmp4,qrzodt(iq,k) )

              !**********PROCESS**********
              zpracs(iq,k) = psacr(iq,k)  ! accretion of rain by snow
              !---------------------------

              tmp1=0.25*pi*esr*n0_r(iq,k)*n0_s(iq,k)*abs( vtr(iq,k)-vts(iq,k) )
              tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
              tmp3=tmp1*tmp2
              npsacr(iq,k)=amin1( tmp3,nrzodt(iq,k) )
            !
            !
            ! (2) FREEZING OF RAIN TO FORM GRAUPEL  (pgfr): Lin (45), added to PI
            !     positive value
            !     Constant in Bigg freezing Aplume=Ap=0.66 /k
            !     Constant in raindrop freezing equ. Bplume=Bp=100./m/m/m/s
            !

              if (qrz(iq,k) .gt. 1.e-8 ) then
                Bp=100.
                Ap=0.66
                tmp1=olambdar(iq,k)*olambdar(iq,k)*olambdar(iq,k)
                tmp2=20.*pi*pi*Bp*n0_r(iq,k)*rhowater*orho(iq,k)*  &
                  (exp(-Ap*temcc(iq,k))-1.0)*tmp1*tmp1*olambdar(iq,k)
                pgfr(iq,k)=amin1( tmp2,qrzodt(iq,k) )
                npgfr(iq,k)=pi*Bp*n0_r(iq,k)*tmpa*tmpa*(exp(-Ap*temcc(iq,k))-1.0)

                !**********PROCESS**********
                zpgfr(iq,k) = pgfr(iq,k)    ! feezing of rain to form graupel (added to PI)
                !---------------------------
              else
                pgfr(iq,k)=0
                npgfr(iq,k)=0.
                !**********PROCESS**********
                zpgfr(iq,k) = pgfr(iq,k)    ! feezing of rain to form graupel (added to PI)
                !---------------------------
              endif
            end if ! for the go to 1200
          end if   !1200    continue

        else  !(T>273.15)    

        !
        !***********************************************************************
        !*********        snow production processes for T > 0 C       **********
        !***********************************************************************
        !
          if (qsz(iq,k) .gt. 0.0)  then !go to 1400
          !
          ! (1) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            save1 =aa_s(iq,k)*sqrho(iq,k)*n0_s(iq,k)* &
                   ggamma(bv_s(iq,k)+tmp_sa(iq,k))*     &
                   olambdas(iq,k)**(bv_s(iq,k)+tmp_sa(iq,k))

            tmp1=esw*save1
            psacw(iq,k)=qlzodt(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsacw(iq,k)=tmp1*ncz(iq,k)
            !**********PROCESS**********
            zpsacw(iq,k) = psacw(iq,k)  ! ! accretion of cloud water by snow
            !---------------------------
          !
          ! (2) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
          !
            esr=1.0
            tmpa=olambdar(iq,k)*olambdar(iq,k)
            tmpb=olambdas(iq,k)*olambdas(iq,k)
            tmpc=olambdar(iq,k)*olambdas(iq,k)
            tmp1=pi*pi*esr*n0_r(iq,k)*n0_s(iq,k)*   &
                    abs( vtr(iq,k)-vts(iq,k) )*orho(iq,k)
          ! tmp1=pi*pi*esr*n0_r(k)*N0_s(k)*            &
          ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5*orho(k)
            tmp2=tmpa*tmpa*olambdas(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
            tmp3=tmp1*rhowater*tmp2
            psacr(iq,k)=amin1( tmp3,qrzodt(iq,k) )

            !**********PROCESS**********
            zpracs(iq,k) = psacr(iq,k)  ! accretion of rain by snow
            !---------------------------

            tmp1=0.25*pi*esr*n0_r(iq,k)*n0_s(iq,k)*abs( vtr(iq,k)-vts(iq,k) )
            tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
            tmp3=tmp1*tmp2
            npsacr(iq,k)=amin1( tmp3,nrzodt(iq,k) )
          !
          ! (3) MELTING OF SNOW (Psmlt): Lin (32)
          !     Psmlt is negative value
          !
            delrs=rs0(iq,k)-qvz(iq,k)
            term1=2.0*pi*orho(iq,k)*( xlv*diffwv(iq,k)*rho(iq,k)*delrs- &
                  xka(iq,k)*temcc(iq,k) )
            tmp1= av_s(iq,k)*sqrho(iq,k)*        &
                  olambdas(iq,k)**(5+bv_s(iq,k)+2*mu_s)/visc(iq,k)
            tmp2= n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                  vf2s*schmidt(iq,k)**0.33334* &
                  ggamma(2.5+0.5*bv_s(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=term1*oxlf*tmp2-cwoxlf*temcc(iq,k)*( psacw(iq,k)+psacr(iq,k) )
            tmp4=amin1(0.0,tmp3)
            psmlt(iq,k)=amax1( tmp4,-qszodt(iq,k) )

            !**********PROCESS**********
            zpsmlt(iq,k) = psmlt(iq,k)  ! melting of snow
            !---------------------------

            if(qsz(iq,k) .ge. 0.0) then
              npsmlt(iq,k)=psmlt(iq,k)*nsz(iq,k)/qsz(iq,k)
            else
              npsmlt(iq,k)=psmlt(iq,k)/xms
            end if
          !
          ! (4) EVAPORATION OF MELTING SNOW (Psmltevp): HR (A27)
          !     but use Lin et al. coefficience
          !     Psmltevp is a negative value
          !
            tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
            tmpb=xlv*xlv*rho(iq,k)*qswz(iq,k)*diffwv(iq,k)
            tmpc=tmpa*qswz(iq,k)*diffwv(iq,k)
            tmpd=amin1( 0.0,(qvoqswz(iq,k)-0.90)*qswz(iq,k)*odtb )

            abr=2.0*pi*(qvoqswz(iq,k)-0.90)*tmpc/(tmpa+tmpb)
          !
          !**** allow evaporation to occur when RH less than 90%
          !**** here not using 100% because the evaporation cooling
          !**** of temperature is not taking into account yet; hence,
          !**** the qsw value is a little bit larger. This will avoid
          !**** evaporation can generate cloud.
          
            tmp1=av_s(iq,k)*sqrho(iq,k)*    &
                 olambdas(iq,k)**(5+bv_s(iq,k)+2*mu_s)/visc(iq,k)
            tmp2=n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                 vf2s*schmidt(iq,k)**0.33334* &
                 ggamma(2.5+0.5*bv_s(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=amin1(0.0,tmp2)
            tmp3=amax1( tmp3,tmpd )
            psmltevp(iq,k)=amax1( tmp3,-qszodt(iq,k) )
 
            !**********PROCESS**************
            zpsmltevp(iq,k) = psmltevp(iq,k)! evaporation of melting snow (T>0)
            !-------------------------------

            if(qsz(iq,k) .ge. 0.0) then
              npsmltevp(iq,k)=psmltevp(iq,k)*nsz(iq,k)/qsz(iq,k)
            else
              npsmltevp(iq,k)=psmltevp(iq,k)/xmr
            end if
          end if !          1400     continue
        end if      !---- end of snow/ice processes   if (tem(iq,k) .lt. 273.15) then
        !---------- end of snow/ice processes below freezing

        !***********************************************************************
        !*********           rain production processes                **********
        !***********************************************************************
        !
        ! (1) AUTOCONVERSION OF RAIN (Praut): using Liu and Daum (2004)
        !
        !---- YLIN, autoconversion use Liu and Daum (2004), unit = g cm-3 s-1, in the scheme kg/kg s-1, so

        if (qlz(iq,k) .gt. 1e-6) then
          !mu_c    = AMIN1(15., (1000.E6/ncz(iq,k) + 2.))                              ! spectral shape parameter for droplets

          tmp = prez(iq,k)/(287.15*tem(iq,k))                      ! temperature in K
          mu_c=0.0005714*(ncz(iq,k)/1.E6*tmp)+0.2714
          mu_c=1./(mu_c**2)-1.
          mu_c=MAX(mu_c,2.)
          mu_c=MIN(mu_c,10.)

          lammin = (mu_c+1.)/60.E-6
          lammax = (mu_c+1.)/1.E-6                                  ! refered from Hugh Morrison
          lamc(iq,k) = (ncz(iq,k)*rhowater*pi*ggamma(4.+mu_c)/(6.*qlz(iq,k)*ggamma(1+mu_c)))**(1./3)

          lamc(iq,k) = max(lamc(iq,k),lammin)
          lamc(iq,k) = min(lamc(iq,k),lammax)

          Dc_liu  = (ggamma(6+1+mu_c)/ggamma(1+mu_c))**(1./6.)/lamc(iq,k)             !----- R6 in m A3 with p=6
          if (Dc_liu .gt. R6c) then
            disp = 1./(mu_c+1.)      !--- square of relative dispersion               ! A7
            eta  = (0.75/pi/(1e-3*rhowater))**2*1.9e11*((1+3*disp)*(1+4*disp)*&
                   (1+5*disp)/(1+disp)/(1+2*disp))                                    ! part of 28c LIU&DAM
            praut(iq,k) = eta*(1e-3*rho(iq,k)*qlz(iq,k))**3/(1e-6*ncz(iq,k))          !--- g cm-3 s-1
            praut(iq,k) = praut(iq,k)/(1e-3*rho(iq,k))                                !--- kg kg-1 s-1
            npraut_r(iq,k) = praut(iq,k)/xmr                                          !--- kg kg-1 s-1
            npraut(iq,k) = praut(iq,k)/qlz(iq,k)*ncz(iq,k)                            !--- kg kg-1 s-1
            npraut(iq,k) = praut(iq,k)/xmr                                            !--- kg kg-1 s-1
            !**********PROCESS**********
            zpraut(iq,k) = praut(iq,k)  ! autoconversion of rain
            !---------------------------
          else
            praut(iq,k) = 0.0
            npraut(iq,k) = 0.0
            npraut_r(iq,k) = 0.0
            !**********PROCESS**********
            zpraut(iq,k) = praut(iq,k)  ! autoconversion of rain
            !---------------------------
          endif
        else
          praut(iq,k) = 0.0
          npraut(iq,k) = 0.0
          npraut_r(iq,k) = 0.0
          !**********PROCESS**********
          zpraut(iq,k) = praut(iq,k)  ! autoconversion of rain
          !---------------------------
        endif
        ! if (qlz(k) .gt. 1e-6) then
        ! praut(k)=1350.*qlz(k)**2.47*  &
        ! (ncz(k)/1.e6*rho(k))**(-1.79)
        ! npraut_r(k) = praut(k)/xmr
        ! npraut(k) = praut(k)/(qlz(k)/ncz(k))
        ! npraut(K) = MIN(npraut(k),nczodt(k))
        ! npraut_r(K) = MIN(npraut_r(k),npraut(k))
        ! endif

        !
        ! (2) ACCRETION OF CLOUD WATER BY RAIN (Pracw): Lin (51)
        !
        erw=1.0
        tmp1=pio4*erw*n0_r(iq,k)*av_r*sqrho(iq,k)* &
             gambp3*olambdar(iq,k)**bp3 ! son
        pracw(iq,k)=qlzodt(iq,k)*( 1.0-exp(-tmp1*dtb) )
        npracw(iq,k)=tmp1*ncz(iq,k)

        !**********PROCESS**********
        zpracw(iq,k) = pracw(iq,k)  ! accretion of cloud water by rain
        !---------------------------
        
        !
        ! (3) EVAPORATION OF RAIN (Prevp): Lin (52)
        !     Prevp is negative value
        !
        !     Sw=qvoqsw : saturation ratio
        !
        tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
        tmpb=xlv*xlv*rho(iq,k)*qswz(iq,k)*diffwv(iq,k)
        tmpc=tmpa*qswz(iq,k)*diffwv(iq,k)
        tmpd=amin1(0.0,(qvoqswz(iq,k)-0.99)*qswz(iq,k)*odtb)

        abr=2.0*pi*(qvoqswz(iq,k)-0.99)*tmpc/(tmpa+tmpb)
        tmp1=av_r*sqrho(iq,k)*olambdar(iq,k)**bp5/visc(iq,k) !son
        tmp2=abr*n0_r(iq,k)*( vf1r*olambdar(iq,k)*olambdar(iq,k)+  &
             vf2r*schmidt(iq,k)**0.33334*gambp5o2*sqrt(tmp1) )
        tmp3=amin1( 0.0,tmp2 )
        tmp3=amax1( tmp3,tmpd )
        prevp(iq,k)=amax1( tmp3,-qrzodt(iq,k) )

        !**********PROCESS************
        zprevp(iq,k) = prevp(iq,k)! evaporation of rain
        !-----------------------------

        if (qrz(iq,k).gt.0.) then
          nprevp(iq,k)=prevp(iq,k)*nrz(iq,k)/qrz(iq,k)
        else
          nprevp(iq,k)=prevp(iq,k)*xmr
        end if

        ! CALL wrf_debug ( 100 , 'module_ylin: finish rain processes' )
        !
        !**********************************************************************
        !*****     combine all processes together and avoid negative      *****
        !*****     water substances
        !***********************************************************************
        !
        if ( temcc(iq,k) .lt. 0.0) then
        !
        !  combined water vapor depletions
        !
          tmp=psdep(iq,k) + midep(iq,k)
          if ( tmp .gt. qvzodt(iq,k) ) then
            factor=qvzodt(iq,k)/tmp
            psdep(iq,k)=psdep(iq,k)*factor
            midep(iq,k)=midep(iq,k)*factor
          end if
        !
        !  combined cloud water depletions
        !
          tmp=praut(iq,k)+psacw(iq,k)+psfw(iq,k)+pracw(iq,k)
          if ( tmp .gt. qlzodt(iq,k) ) then
            factor=qlzodt(iq,k)/tmp
            praut(iq,k)=praut(iq,k)*factor
            psacw(iq,k)=psacw(iq,k)*factor
            psfw(iq,k)=psfw(iq,k)*factor
            pracw(iq,k)=pracw(iq,k)*factor
          end if
        !
        !  combined cloud ice depletions
        !
          tmp=psaut(iq,k)+psaci(iq,k)+praci(iq,k)+psfi(iq,k)
          if (tmp .gt. qizodt(iq,k) ) then
            factor=qizodt(iq,k)/tmp
            psaut(iq,k)=psaut(iq,k)*factor
            psaci(iq,k)=psaci(iq,k)*factor
            praci(iq,k)=praci(iq,k)*factor
            psfi(iq,k)=psfi(iq,k)*factor
          endif

        !
        !  combined all rain processes
        !
          tmp_r=piacr(iq,k)+psacr(iq,k)-prevp(iq,k)-  & 
                praut(iq,k)-pracw(iq,k)+pgfr(iq,k)
          if (tmp_r .gt. qrzodt(iq,k) ) then
            factor=qrzodt(iq,k)/tmp_r
            piacr(iq,k)=piacr(iq,k)*factor
            psacr(iq,k)=psacr(iq,k)*factor
            prevp(iq,k)=prevp(iq,k)*factor
            pgfr(iq,k)=pgfr(iq,k)*factor
          endif
        !
        !   combined all snow processes
        !
          tmp_s=-pssub(iq,k)-(psaut(iq,k)+psaci(iq,k)+ &
                 psacw(iq,k)+psfw(iq,k)+pgfr(iq,k)+ &
                 psfi(iq,k)+praci(iq,k)+piacr(iq,k)+ &
                 psdep(iq,k)+psacr(iq,k)-pracs(iq,k))
          if ( tmp_s .gt. qszodt(iq,k) ) then
            factor=qszodt(iq,k)/tmp_s
            pssub(iq,k)=pssub(iq,k)*factor
            pracs(iq,k)=pracs(iq,k)*factor
          endif

        !
        ! calculate new water substances, thetae, tem, and qvsbar
        !

          pvapor(iq,k)=-pssub(iq,k)-psdep(iq,k)-prevp(iq,k)-midep(iq,k)
          qvz(iq,k)=amax1( qvmin,qvz(iq,k)+dtb*pvapor(iq,k) )
          pclw(iq,k)=-praut(iq,k)-pracw(iq,k)-psacw(iq,k)-psfw(iq,k)
          qlz(iq,k)=amax1( 0.0,qlz(iq,k)+dtb*pclw(iq,k) )
          pcli(iq,k)=-psaut(iq,k)-psfi(iq,k)-psaci(iq,k)- & 
                    praci(iq,k)+midep(iq,k)
          qiz(iq,k)=amax1( 0.0,qiz(iq,k)+dtb*pcli(iq,k) )
          tmp_r=piacr(iq,k)+psacr(iq,k)-prevp(iq,k)-praut(iq,k)- &
                    pracw(iq,k)+pgfr(iq,k)-pracs(iq,k)
          prain(iq,k)=-tmp_r
          qrz(iq,k)=amax1( 0.0,qrz(iq,k)+dtb*prain(iq,k) )
          tmp_s=-pssub(iq,k)-(psaut(iq,k)+psaci(iq,k)+ &
                 psacw(iq,k)+psfw(iq,k)+pgfr(iq,k)+  &
                 psfi(iq,k)+praci(iq,k)+piacr(iq,k)+  &
                 psdep(iq,k)+psacr(iq,k)-pracs(iq,k))
          psnow(iq,k)=-tmp_s
          qsz(iq,k)=amax1( 0.0,qsz(iq,k)+dtb*psnow(iq,k) )
          qschg(iq,k)=qschg(iq,k)+psnow(iq,k)
          qschg(iq,k)=psnow(iq,k)

          !**********PROCESS**********
          zpvapor(iq,k)=pvapor(iq,k)  ! sum all process for water vapor to determine qvz
          zpclw(iq,k)=  pclw(iq,k)    ! sum all process for cloud liquid to determine qlz
          zpcli(iq,k)=  pcli(iq,k)    ! sum all process for cloud ice to determine qiz
          zprain(iq,k)= prain(iq,k)   ! sum all process for rain
          zpsnow(iq,k)= psnow(iq,k)   ! sum all process for snow
          !---------------------------

          tmp=ocp/tothz(iq,k)*xLf*qschg(iq,k)
          theiz(iq,k)=theiz(iq,k)+dtb*tmp
          ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
          ! tem(k)=thz(k)*tothz(k)
          ! temcc(k)=tem(k)-273.15
          !==================update temperature=================================================
          temcc(iq,k)=tem(iq,k)-273.15
          lvap = xlv + (2106.0 - 4218.0)*temcc(iq,k)  !Enthalpy of vaporization
          tmp1=(pssub(iq,k)+psdep(iq,k))*xls*ocp + prevp(iq,k)*lvap*ocp+  &
               (psfw(iq,k)+pgfr(iq,k)+psacr(iq,k)-pracs(iq,k))*xlfocp
          !bug fixed 20191126
          tem(iq,k)=tem(iq,k)+tmp1*dtb
          temcc(iq,k)=tem(iq,k)-273.15
          thz(iq,k)=tem(iq,k)/tothz(iq,k)
          !===================================================================
          if( temcc(iq,k) .lt. -40.0 ) qswz(iq,k)=qsiz(iq,k)
            qlpqi=qlz(iq,k)+qiz(iq,k)
            if ( qlpqi .eq. 0.0 ) then
              qvsbar(iq,k)=qsiz(iq,k)
            else
              qvsbar(iq,k)=(qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k))/qlpqi
            endif
            tmp1=-npraut(iq,k)-npracw(iq,k)-npsacw(iq,k)
            ncz(iq,k)=amax1( 0.0,ncz(iq,k)+dtb*tmp1 )
            tmp1=-npsaut(iq,k)-npsaci(iq,k)-npraci(iq,k)+nidep(iq,k)
            niz(iq,k)=amax1( 0.0,niz(iq,k)+dtb*tmp1 )
            tmp1=npiacr(iq,k)+npsacr(iq,k)-nprevp(iq,k)-npraut_r(iq,k)+npgfr(iq,k)
            nrz(iq,k)=amax1( 0.0,nrz(iq,k)-dtb*tmp1 )
            tmp1=-(npsaut(iq,k)+npgfr(iq,k)+  &
                   npraci(iq,k)+npiacr(iq,k)+  &
                   npsdep(iq,k)+npsacr(iq,k))
            nsz(iq,k)=amax1( 0.0,nsz(iq,k)-dtb*tmp1 )

          !----------------------------------------
          else                  !>0 C
          !----------------------------------------
          !
          !  combined cloud water depletions
          !
            tmp=praut(iq,k)+psacw(iq,k)+pracw(iq,k)
            if ( tmp .gt. qlzodt(iq,k) ) then
              factor=qlzodt(iq,k)/tmp
              praut(iq,k)=praut(iq,k)*factor
              psacw(iq,k)=psacw(iq,k)*factor
              pracw(iq,k)=pracw(iq,k)*factor
            end if
          !
          !  combined all snow processes
          !
            tmp_s=-(psmlt(iq,k)+psmltevp(iq,k))
            if (tmp_s .gt. qszodt(iq,k) ) then
              factor=qszodt(iq,k)/tmp_s
              psmlt(iq,k)=psmlt(iq,k)*factor
              psmltevp(iq,k)=psmltevp(iq,k)*factor
            endif
          !
          !  combined all rain processes
          !
            tmp_r=-prevp(iq,k)-(praut(iq,k)+pracw(iq,k)+psacw(iq,k)-psmlt(iq,k))
            if (tmp_r .gt. qrzodt(iq,k) ) then
              factor=qrzodt(iq,k)/tmp_r
              prevp(iq,k)=prevp(iq,k)*factor
            endif
          !
          !  calculate new water substances and thetae
          !
            pvapor(iq,k)=-psmltevp(iq,k)-prevp(iq,k)
            qvz(iq,k)=amax1( qvmin,qvz(iq,k)+dtb*pvapor(iq,k))
            pclw(iq,k)=-praut(iq,k)-pracw(iq,k)-psacw(iq,k)
            qlz(iq,k)=amax1( 0.0,qlz(iq,k)+dtb*pclw(iq,k) )
            pcli(iq,k)=0.0
            qiz(iq,k)=amax1( 0.0,qiz(iq,k)+dtb*pcli(iq,k) )
            tmp_r=-prevp(iq,k)-(praut(iq,k)+pracw(iq,k)+psacw(iq,k)-psmlt(iq,k))
            prain(iq,k)=-tmp_r
            tmpqrz=qrz(iq,k)
            qrz(iq,k)=amax1( 0.0,qrz(iq,k)+dtb*prain(iq,k) )
            tmp_s=-(psmlt(iq,k)+psmltevp(iq,k))
            psnow(iq,k)=-tmp_s
            qsz(iq,k)=amax1( 0.0,qsz(iq,k)+dtb*psnow(iq,k) )
            qschg(iq,k)=psnow(iq,k)

            !**********PROCESS**********
            zpvapor(iq,k)=pvapor(iq,k)  ! sum all process for water vapor to determine qvz
            zpclw(iq,k)=  pclw(iq,k)    ! sum all process for cloud liquid to determine qlz
            zpcli(iq,k)=  pcli(iq,k)    ! sum all process for cloud ice to determine qiz
            zprain(iq,k)= prain(iq,k)   ! sum all process for rain
            zpsnow(iq,k)= psnow(iq,k)   ! sum all process for snow
            !---------------------------

            tmp=ocp/tothz(iq,k)*xLf*qschg(iq,k)
            theiz(iq,k)=theiz(iq,k)+dtb*tmp
            ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
            ! tem(k)=thz(k)*tothz(k)
            ! temcc(k)=tem(k)-273.15
            !==================update tmperature=================================================
            temcc(iq,k)=tem(iq,k)-273.15
            lvap = xlv + (2106.0 - 4218.0)*temcc(iq,k)  !Enthalpy of vaporization
            tmp1=psmltevp(iq,k)*xls*ocp + prevp(iq,k)*lvap*ocp+  &
                 psmlt(iq,k)*xlfocp
            !tmp1 =  ! 1. evaporation of rain formed by melting snow ??? (-)
                     ! 2. evaporation of rain (-)
                     ! 3. melting of snow to form rain (+)
            tem(iq,k)=tem(iq,k)+tmp1*dtb
            !bugfix 20191126

            !tem(k)=tem(k)+tmp1*dtb
            temcc(iq,k)=tem(iq,k)-273.15

            thz(iq,k)=tem(iq,k)/tothz(iq,k)

            !===================================================================
            es=1000.*svp1*exp( svp2*temcc(iq,k)/(tem(iq,k)-svp3) )
            qswz(iq,k)=ep2*es/(prez(iq,k)-es)
            qsiz(iq,k)=qswz(iq,k)
            qvsbar(iq,k)=qswz(iq,k)
            ! tmp1=-(npraut(iq,k)+npsacw(iq,k)+npracw(iq,k))
            ncz(iq,k)=amax1( 0.0,ncz(iq,k)+dtb*tmp1)
            tmp1=-nprevp(iq,k)-(npraut_r(iq,k)-npsmlt(iq,k))
            ! tmp1=-nprevp(k)-(nprautr(k)+npracwr(k)+npsacw(k)-npsmltr(k))
            nrz(iq,k)=amax1(0.0,nrz(iq,k)-dtb*tmp1)
            tmp1=-(npsmlt(iq,k)+npsmltevp(iq,k))
            nsz(iq,k)=amax1( 0.0,nsz(iq,k)-dtb*tmp1 )
          end if    !T seperate for source and sink terms
          ! CALL wrf_debug ( 100 , 'module_ylin: finish sum of all processes' )

        !rain
        if (qrz(iq,k) .gt. 1.0e-8) then
          xlambdar(iq,k)=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
          if (xlambdar(iq,k).lt.lamminr) then
            xlambdar(iq,k) = lamminr
            n0_r(iq,K) = xlambdar(iq,K)**4*qrz(iq,K)/(pi*rhowater)
            nrz(iq,K) = n0_r(iq,K)/xlambdar(iq,K)
          else if (xlambdar(iq,K).gt.lammaxr) then
            xlambdar(iq,K) = lammaxr
            n0_r(iq,K) = xlambdar(iq,K)**4*qrz(iq,K)/(pi*rhowater)
            nrz(iq,K) = n0_r(iq,K)/xlambdar(iq,K)
          end if
        end if

        !snow
        if (qsz(iq,k) .gt. 1.0e-8) then
          xlambdas(iq,k)=(am_s(iq,k)*ggamma(tmp_ss(iq,k))*     &
          nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
          if (xlambdas(iq,k).lt.lammins) then
            xlambdas(iq,k)= lamminS
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*      &
                           qsz(iq,K)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          else if (xlambdas(iq,k).gt.lammaxs) then
            xlambdas(iq,k) = lammaxs
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*      & 
                           qsz(iq,K)/ggamma(1+bm_s(iq,k))/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          end if
        end if

        !cloud ice
        if (qiz(iq,k).ge.1.0e-8) then
          lami(iq,k) = max((ggamma(1.+3.)*500.*pi/6.)*niz(iq,k)/qiz(iq,k),1.e-20)**(1./3) !fixed zdc
          if (lami(iq,k).lt.lammini) then
            lami(iq,k)= lammini
            n0_i(iq,K) = lami(iq,k)**4./ggamma(1.+3.)*500.*pi/6.
            niz(iq,K) = n0_i(iq,K)/lami(iq,k)
          else if (lami(iq,k).gt.lammaxi) then
            lami(iq,k) = lammaxi
            n0_i(iq,K) = lami(iq,k)**4./ggamma(1.+3.)*500.*pi/6.
            niz(iq,K) = n0_i(iq,K)/lami(iq,k)
          end if
        end if

        !cloud water zdc 20220208
        if (qlz(iq,k).ge.1.0e-8) then
          lamc(iq,k) = (ncz(iq,k)*rhowater*pi*ggamma(4.+mu_c)/(6.*qlz(iq,k)*ggamma(1+mu_c)))**(1./3)
          if (lamc(iq,k).lt.lammin) then
            lamc(iq,k)= lammin
            n0_c(iq,k)= lamc(iq,k)**(mu_c+4.)*6.*qlz(iq,k)/(pi*rhowater*ggamma(mu_c+4))
            ncz(iq,k) = n0_c(iq,k)/lamc(iq,k)
          else if (lamc(iq,k).gt.lammax) then
            lamc(iq,k)= lammax
            n0_c(iq,k)= lamc(iq,k)**(mu_c+4.)*6.*qlz(iq,k)/(pi*rhowater*ggamma(mu_c+4))
            ncz(iq,k) = n0_c(iq,k)/lamc(iq,k)
          end if
        end if

        !
        !***********************************************************************
        !**********              saturation adjustment                **********
        !***********************************************************************
        !
        !    allow supersaturation exits linearly from 0% at 500 mb to 50%
        !    above 300 mb
        !    5.0e-5=1.0/(500mb-300mb)
        !
        rsat=1.0
        if( qvz(iq,k)+qlz(iq,k)+qiz(iq,k) .lt. rsat*qvsbar(iq,k) ) then ! goto 1800

        !
        !   unsaturated
        !
          qvz(iq,k)=qvz(iq,k)+qlz(iq,k)+qiz(iq,k)
          qlz(iq,k)=0.0
          qiz(iq,k)=0.0

          thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)

          tem(iq,k)=thz(iq,k)*tothz(iq,k)
          temcc(iq,k)=tem(iq,k)-273.15

        else
        !
        !   saturated
        !
          pladj(iq,k)=qlz(iq,k)
          piadj(iq,k)=qiz(iq,k)
        !

          CALL satadj(qvz(iq,:), qlz(iq,:), qiz(iq,:), prez(iq,:), &
                      theiz(iq,:), thz(iq,:), tothz(iq,:), kts, kte, &
                      k, xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

          pladj(iq,k)=odtb*(qlz(iq,k)-pladj(iq,k))
          piadj(iq,k)=odtb*(qiz(iq,k)-piadj(iq,k))
          pclw(iq,k)=pclw(iq,k)+pladj(iq,k)
          pcli(iq,k)=pcli(iq,k)+piadj(iq,k)
          pvapor(iq,k)=pvapor(iq,k)-( pladj(iq,k)+piadj(iq,k) )
          thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
          tem(iq,k)=thz(iq,k)*tothz(iq,k)
          temcc(iq,k)=tem(iq,k)-273.15

          es=1000.*svp1*exp( svp2*temcc(iq,k)/(tem(iq,k)-svp3) )
          qswz(iq,k)=ep2*es/(prez(iq,k)-es)
          if (tem(iq,k) .lt. 273.15 ) then
            es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
            qsiz(iq,k)=ep2*es/(prez(iq,k)-es)
            if (temcc(iq,k) .lt. -40.0) qswz(iq,k)=qsiz(iq,k)
          else
            qsiz(iq,k)=qswz(iq,k)
          endif
          qlpqi=qlz(iq,k)+qiz(iq,k)
          if ( qlpqi .eq. 0.0 ) then
            qvsbar(iq,k)=qsiz(iq,k)
          else
            qvsbar(iq,k)=( qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k) )/qlpqi
          endif

        !
        !***********************************************************************
        !*****     melting and freezing of cloud ice and cloud water       *****
        !***********************************************************************
          qlpqi=qlz(iq,k)+qiz(iq,k)
          if( qlpqi .gt. 0.0 ) then !go to 1800
          !
          !
          ! (1)  HOMOGENEOUS NUCLEATION WHEN T< -40 C (Pihom)
          !
            if(temcc(iq,k) .lt. -40.0) then
              pihom(iq,k)=qlz(iq,k)*odtb
              nihom(iq,k)=ncz(iq,k)*odtb

              !**********PROCESS**********
              zpihom(iq,k) = pihom(iq,k)  ! homogeneous nucleation <-40
              !---------------------------
            end if
          !
          ! (2)  MELTING OF ICE CRYSTAL WHEN T> 0 C (Pimlt)
          !
            if(temcc(iq,k) .gt. 0.0) then
              pimlt(iq,k)=qiz(iq,k)*odtb
              nimlt(iq,k)=niz(iq,k)*odtb
              !**********PROCESS**********
              zpimlt(iq,k) = pimlt(iq,k)  ! melting of ice crystal >0.
              !---------------------------
            end if
          !
          ! (3) PRODUCTION OF CLOUD ICE BY BERGERON PROCESS (Pidw): Hsie (p957)
          !     this process only considered when -31 C < T < 0 C
          !
            if(temcc(iq,k) .lt. 0.0 .and. temcc(iq,k) .gt. -31.0) then
          !!
          !!   parama1 and parama2 functions must be user supplied
          !!
              a1=parama1( temcc(iq,k) )
              a2=parama2( temcc(iq,k) )
              !! change unit from cgs to mks
              a1=a1*0.001**(1.0-a2)
              xnin=xni0*exp(-bni*temcc(iq,k))
              pidw(iq,k)=xnin*orho(iq,k)*(a1*xmnin**a2)

              !**********PROCESS**********
              zpidw(iq,k)= pidw(iq,k)     ! production of cloud ice by BERGERON process
              !---------------------------
            end if

            pcli(iq,k)=pcli(iq,k)+pihom(iq,k)-pimlt(iq,k)+pidw(iq,k)
            pclw(iq,k)=pclw(iq,k)-pihom(iq,k)+pimlt(iq,k)-pidw(iq,k)

            !**********PROCESS*******
            zpiadj(iq,k)= pcli(iq,k) ! zpiadj - zpcli = +pihom(iq,k)-pimlt(iq,k)+pidw(iq,k) 
            zpladj(iq,k)= pclw(iq,k) ! zpladj - zpclw = -pihom(iq,k)+pimlt(iq,k)-pidw(iq,k) 
            !------------------------

            qlz(iq,k)=amax1( 0.0,qlz(iq,k)+dtb*(-pihom(iq,k)+pimlt(iq,k)-pidw(iq,k)) )
            qiz(iq,k)=amax1( 0.0,qiz(iq,k)+dtb*(pihom(iq,k)-pimlt(iq,k)+pidw(iq,k)) )

            ncz(iq,k)=amax1( 0.0,ncz(iq,k)+dtb*(-nihom(iq,k)+nimlt(iq,k)) )
            niz(iq,k)=amax1( 0.0,niz(iq,k)+dtb*( nihom(iq,k)-nimlt(iq,k)) )

            CALL satadj(qvz(iq,:), qlz(iq,:), qiz(iq,:), prez(iq,:), &
                    theiz(iq,:), thz(iq,:), tothz(iq,:), kts, kte, &
                    k, xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

            thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
            tem(iq,k)=thz(iq,k)*tothz(iq,k)
            temcc(iq,k)=tem(iq,k)-273.15
            es=1000.*svp1*exp( svp2*temcc(iq,k)/(tem(iq,k)-svp3) )
            qswz(iq,k)=ep2*es/(prez(iq,k)-es)

            if (tem(iq,k) .lt. 273.15 ) then
              es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
              qsiz(iq,k)=ep2*es/(prez(iq,k)-es)
              if (temcc(iq,k) .lt. -40.0) qswz(iq,k)=qsiz(iq,k)
            else
              qsiz(iq,k)=qswz(iq,k)
            endif
            qlpqi=qlz(iq,k)+qiz(iq,k)

            if ( qlpqi .eq. 0.0 ) then
              qvsbar(iq,k)=qsiz(iq,k)
            else
              qvsbar(iq,k)=( qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k) )/qlpqi
            endif
          end if ! 1800  continue
        end if ! 1800  continue
        !
        !***********************************************************************
        !**********    integrate the productions of rain and snow     **********
        !***********************************************************************
        !
      end if
    end do     ! iq
  ENDDO        ! k

  !call END_LOG(p4_end)
  !call START_LOG(p5_begin)

  do k=kts+1,kte
    where ( qvz(1:imax,k) .lt. qvmin ) 
      qlz(1:imax,k)=0.0
      qiz(1:imax,k)=0.0
      ncz(1:imax,k)=0.0
      niz(1:imax,k)=0.0
      qvz(1:imax,k)=amax1( qvmin,qvz(1:imax,k)+qlz(1:imax,k)+qiz(1:imax,k) )
    end where
      niz(1:imax,k) = min(niz(1:imax,k),0.3E6/rho(1:imax,k))
      ncz(1:imax,k) = min(ncz(1:imax,k),250000.E6/rho(1:imax,k))
      ncz(1:imax,k) = max(ncz(1:imax,k),0.01E6/rho(1:imax,k))
  end do


  ! CALCULATE EFFECTIVE RADIUS zdc 20220208
  do k=kts,kte
    where (qiz(1:imax,k) .gt. 1.0e-8 .and. lami(1:imax,k)>0. ) 
      EFFI1D(1:imax,k) = 3./LAMI(1:imax,k)/2.
    elsewhere
      EFFI1D(1:imax,k) = 25.E-6
    end where

    where (qsz(1:imax,k) .gt. 1.0e-8) 
      EFFS1D(1:imax,k) = 3./xlambdas(1:imax,k)/2.
    elsewhere
      EFFS1D(1:imax,k) = 25.E-6
    end where

    where (qrz(1:imax,k) .gt. 1.0e-8) 
      EFFR1D(1:imax,k) = 3./xlambdar(1:imax,k)/2.
    elsewhere
      EFFR1D(1:imax,k) = 25.E-6
    end where

    where (qlz(1:imax,k) .gt. 1.0e-8 .and. lamc(1:imax,k) >0.)
      EFFC1D(1:imax,k) = GAMMA(mu_c+4.)/GAMMA(mu_c+3.)/lamc(1:imax,k)/2.
    elsewhere
      EFFC1D(1:imax,k) = 25.E-6
    end where
  end do

  ! save all process rate for understanding cloud microphysics
!  do k=kts,kte
!    zpsnow(1:imax,k)   = psnow(1:imax,k)    ! sum all process for snow
!    zpsaut(1:imax,k)   = psaut(1:imax,k)    ! ice crystal aggregation to snow
!    zpsfw(1:imax,k)    = psfw(1:imax,k)     ! BERGERON process to transfer cloud water to snow
!    zpsfi(1:imax,k)    = psfi(1:imax,k)     ! BERGERON process to transfer cloud ice to snow
!    zpraci(1:imax,k)   = praci(1:imax,k)    ! cloud ice accretion by rain
!    zpiacr(1:imax,k)   = piacr(1:imax,k)    ! rain accretion by cloud ice
!    zpsaci(1:imax,k)   = psaci(1:imax,k)    ! ice crystal accretion by snow
!    zpsacw(1:imax,k)   = psacw(1:imax,k)    ! accretion of cloud water by snow
!    zpsdep(1:imax,k)   = psdep(1:imax,k)    ! deposition of snow
!    zpssub(1:imax,k)   = pssub(1:imax,k)    ! sublimation of snow (T<0)
!    zpracs(1:imax,k)   = pracs(1:imax,k)    ! accretion of snow by rain
!    zpsacr(1:imax,k)   = psacr(1:imax,k)    ! accretion of rain by snow
!    zpsmlt(1:imax,k)   = psmlt(1:imax,k)    ! melting of snow
!    zpsmltevp(1:imax,k)= psmltevp(1:imax,k) ! evaporation of melting snow (T>0)
!    zprain(1:imax,k)   = prain(1:imax,k)    ! sum all process for rain
!    zpraut(1:imax,k)   = praut(1:imax,k)    ! autoconversion of rain
!    zpracw(1:imax,k)   = pracw(1:imax,k)    ! accretion of cloud water by rain
!    zprevp(1:imax,k)   = prevp(1:imax,k)    ! evaporation of rain
!    zpgfr(1:imax,k)    = pgfr(1:imax,k)     ! feezing of rain to form graupel (added to PI)
!    zpvapor(1:imax,k)  = pvapor(1:imax,k)   ! sum all process for water vapor to determine qvz
!    zpclw(1:imax,k)    = pclw(1:imax,k)     ! sum all process for cloud liquid to determine qlz
!    zpladj(1:imax,k)   = pladj(1:imax,k)    ! saturation adjustment for ql
!    zpcli(1:imax,k)    = pcli(1:imax,k)     ! sum all process for cloud ice to determine qiz
!    zpimlt(1:imax,k)   = pimlt(1:imax,k)    ! melting of ice crystal >0.
!    zpihom(1:imax,k)   = pihom(1:imax,k)    ! homogeneous nucleation <-40
!    zpidw(1:imax,k)    = pidw(1:imax,k)     ! production of cloud ice by BERGERON process
!    zpiadj(1:imax,k)   = piadj(1:imax,k)    ! saturation adjustment for qi
!    zpmidep(1:imax,k)  = midep(1:imax,k)    ! 
!  enddo


  ! save process rate for aerisol scheme
  do k=kts,kte
    fluxi(1:imax,k) = fluxi(1:imax,k)                         ! - ice flux leaving layer k to k-1 (kg/m2/s)
    fluxs(1:imax,k) = fluxs(1:imax,k)                         ! - snow flux leaving layer k to k-1 (kg/m2/s)
    fluxr(1:imax,k) = fluxr(1:imax,k)                         ! - rain flux leaving layer k to k-1 (kg/m2/s)
    fluxg(1:imax,k) = 0.                                        ! - graupel flux leving layer k to k-1 (kg/m2/s)
    fluxm(1:imax,k) = -1.*(psmlt(1:imax,k)*dzw(1:imax,k)*  &!
                        rho(1:imax,k)+psmltevp(1:imax,k)*    &!
                        dzw(1:imax,k)*rho(1:imax,k))          ! - ice melting flux in layer k (kg/m2/s)
    fluxf(1:imax,k) = pgfr(1:imax,k)*dzw(1:imax,k)*        &!
                        rho(1:imax,k)                           ! - liquid freezing flux in layer k (kg/m2/s)
    fevap(1:imax,k) = -1.*prevp(1:imax,k)*dzw(1:imax,k)*   &!
                        rho(1:imax,k)                           ! - evaporation of rainfall flux (kg/m2/s)
    fsubl(1:imax,k) = -1.*pssub(1:imax,k)*dzw(1:imax,k)*   &!
                        rho(1:imax,k)                           ! - sublimation of snow, ice and graupel flux (kg/m2/s)
    fauto(1:imax,k) = praut(1:imax,k)*dzw(1:imax,k)*       &!
                        rho(1:imax,k)                           ! - autoconversion flux for rainfall (kg/m2/s)
    fcoll(1:imax,k) = pracw(1:imax,k)*dzw(1:imax,k)*       &!
                        rho(1:imax,k)                           ! - collection of cloud liquid water by rain (kg/m2/s)
    faccr(1:imax,k) = psacw(1:imax,k)*dzw(1:imax,k)*       &!
                        rho(1:imax,k) + psfw(1:imax,k)*      &! - accretion of cloud liq water by snow,ice and graupel (kg/m2/s)
                        dzw(1:imax,k)*rho(1:imax,k)
    vi(1:imax,k)    = vtiold(1:imax,k)
    vs(1:imax,k)    = vtsold(1:imax,k)
    vg(1:imax,k)    = 0.
  end do

  !call END_LOG(p5_end) 

END SUBROUTINE clphy1d_ylin



!---------------------------------------------------------------------
!                         SATURATED ADJUSTMENT
!---------------------------------------------------------------------
SUBROUTINE satadj(qvz, qlz, qiz, prez, theiz, thz, tothz,      &
                  kts, kte, k, xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)
!---------------------------------------------------------------------
      IMPLICIT NONE
!---------------------------------------------------------------------
!  This program use Newton's method for finding saturated temperature
!  and saturation mixing ratio.
!
! In this saturation adjustment scheme we assume
! (1)  the saturation mixing ratio is the mass weighted average of
!      saturation values over liquid water (qsw), and ice (qsi)
!      following Lord et al., 1984 and Tao, 1989
!
! (2) the percentage of cloud liquid and cloud ice will
!      be fixed during the saturation calculation
!---------------------------------------------------------------------


  INTEGER, INTENT(IN   )             :: kts, kte, k
  REAL,      DIMENSION( kts:kte ),                                   &
                       INTENT(INOUT) :: qvz, qlz, qiz
  REAL,      DIMENSION( kts:kte ),                                   &
                       INTENT(IN   ) :: prez, theiz, tothz
  REAL,     INTENT(IN   )            :: xLvocp, xLfocp, episp0k
  REAL,     INTENT(IN   )            :: EP2,SVP1,SVP2,SVP3,SVPT0

  ! LOCAL VARS

  INTEGER                            :: n
  REAL, DIMENSION( kts:kte )         :: thz, tem, temcc, qsiz,       &
                                        qswz, qvsbar
  REAL :: qsat, qlpqi, ratql, t0, t1, tmp1, ratqi, tsat, absft,    &
             denom1, denom2, dqvsbar, ftsat, dftsat, qpz,es             
!---------------------------------------------------------------------

  thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
  tem(k)=tothz(k)*thz(k)
  if (tem(k) .gt. 273.15) then
  ! qsat=episp0k/prez(k)*  &
  ! exp( svp2*(tem(k)-273.15)/(tem(k)-svp3) )
    es=1000.*svp1*exp( svp2*(tem(k)-svpt0)/(tem(k)-svp3) )
    qsat=ep2*es/(prez(k)-es)
  else
    qsat=episp0k/prez(k)*  &
    exp( 21.8745584*(tem(k)-273.15)/(tem(k)-7.66) )
  end if
  qpz=qvz(k)+qlz(k)+qiz(k)
  if (qpz .lt. qsat) then
    qvz(k)=qpz
    qiz(k)=0.0
    qlz(k)=0.0
    !  go to 400
    ! end if
  else     ! this else to remove the go to above
    qlpqi=qlz(k)+qiz(k)
    if( qlpqi .ge. 1.0e-5) then
      ratql=qlz(k)/qlpqi
      ratqi=qiz(k)/qlpqi
    else
      t0=273.15
      ! t1=233.15
      t1=248.15
      tmp1=( t0-tem(k) )/(t0-t1)
      tmp1=amin1(1.0,tmp1)
      tmp1=amax1(0.0,tmp1)
      ratqi=tmp1
      ratql=1.0-tmp1
    end if
    !
    !
    !--  saturation mixing ratios over water and ice
    !--  at the outset we will follow Bolton 1980 MWR for
    !--  the water and Murray JAS 1967 for the ice
    !
    !-- dqvsbar=d(qvsbar)/dT
    !-- ftsat=F(Tsat)
    !-- dftsat=d(F(T))/dT
    !
    !  First guess of tsat

    tsat=tem(k)
    absft=1.0
    !
    do n=1,20   !do 200 n=1,20
      denom1=1.0/(tsat-svp3)
      denom2=1.0/(tsat-7.66)
      ! qswz(k)=episp0k/prez(k)*  &
      ! exp( svp2*denom1*(tsat-273.15) )
      es=1000.*svp1*exp( svp2*denom1*(tsat-svpt0) )
      qswz(k)=ep2*es/(prez(k)-es)
      if (tem(k) .lt. 273.15) then
      ! qsiz(k)=episp0k/prez(k)*  &
      ! exp( 21.8745584*denom2*(tsat-273.15) )
        es=1000.*svp1*exp( 21.8745584*denom2*(tsat-273.15) )
        qsiz(k)=ep2*es/(prez(k)-es)
        if (tem(k) .lt. 233.15) qswz(k)=qsiz(k)
      else
        qsiz(k)=qswz(k)
      endif
      qvsbar(k)=ratql*qswz(k)+ratqi*qsiz(k)
      ! if( absft .lt. 0.01 .and. n .gt. 3 ) go to 300
      if( absft .ge. 0.01 ) then !go to 300
        dqvsbar=ratql*qswz(k)*svp2*243.5*denom1*denom1+  &
                ratqi*qsiz(k)*21.8745584*265.5*denom2*denom2
        ftsat=tsat+(xlvocp+ratqi*xlfocp)*qvsbar(k)-  &
              tothz(k)*theiz(k)-xlfocp*ratqi*(qvz(k)+qlz(k)+qiz(k))
        dftsat=1.0+(xlvocp+ratqi*xlfocp)*dqvsbar
        tsat=tsat-ftsat/dftsat
        absft=abs(ftsat)
      end if !300   continue
    end do !200   continue
    9020  format(1x,'point can not converge, absft,n=',e12.5,i5)
    !300   continue

    if( qpz .gt. qvsbar(k) ) then
      qvz(k)=qvsbar(k)
      qiz(k)=ratqi*( qpz-qvz(k) )
      qlz(k)=ratql*( qpz-qvz(k) )
    else
      qvz(k)=qpz
      qiz(k)=0.0
      qlz(k)=0.0
    end if
  end if !400  continue
END SUBROUTINE satadj


!----------------------------------------------------------------
FUNCTION parama1_s(temp) result(ans)
!----------------------------------------------------------------
  IMPLICIT NONE
  !----------------------------------------------------------------
  !  This program calculate the parameter for crystal growth rate
  !  in Bergeron process
  !----------------------------------------------------------------

  REAL, INTENT (IN   )   :: temp
  real                   :: ans
  REAL, DIMENSION(32)    :: a1
  INTEGER                :: i1, i1p1
  REAL                   :: ratio

  data a1/0.100e-10,0.7939e-7,0.7841e-6,0.3369e-5,0.4336e-5, &
              0.5285e-5,0.3728e-5,0.1852e-5,0.2991e-6,0.4248e-6, &
              0.7434e-6,0.1812e-5,0.4394e-5,0.9145e-5,0.1725e-4, &
              0.3348e-4,0.1725e-4,0.9175e-5,0.4412e-5,0.2252e-5, &
              0.9115e-6,0.4876e-6,0.3473e-6,0.4758e-6,0.6306e-6, &
              0.8573e-6,0.7868e-6,0.7192e-6,0.6513e-6,0.5956e-6, &
              0.5333e-6,0.4834e-6/

  i1=int(-temp)+1
  i1p1=i1+1
  ratio=-(temp)-float(i1-1)
  ans=a1(i1)+ratio*( a1(i1p1)-a1(i1) )
  END FUNCTION parama1_s

!----------------------------------------------------------------
FUNCTION parama1_v(temp) result(ans)
!----------------------------------------------------------------
  IMPLICIT NONE
  !----------------------------------------------------------------
  !  This program calculate the parameter for crystal growth rate
  !  in Bergeron process
  !----------------------------------------------------------------

  REAL, dimension(:), INTENT (IN   )   :: temp
  real, dimension(size(temp)) :: ans
  REAL, DIMENSION(32)    :: a1
  INTEGER                :: i1, i1p1, ilen, iq
  REAL                   :: ratio

  data a1/0.100e-10,0.7939e-7,0.7841e-6,0.3369e-5,0.4336e-5, &
              0.5285e-5,0.3728e-5,0.1852e-5,0.2991e-6,0.4248e-6, &
              0.7434e-6,0.1812e-5,0.4394e-5,0.9145e-5,0.1725e-4, &
              0.3348e-4,0.1725e-4,0.9175e-5,0.4412e-5,0.2252e-5, &
              0.9115e-6,0.4876e-6,0.3473e-6,0.4758e-6,0.6306e-6, &
              0.8573e-6,0.7868e-6,0.7192e-6,0.6513e-6,0.5956e-6, &
              0.5333e-6,0.4834e-6/

  do iq = 1,size(temp)
    i1=int(-temp(iq))+1
    i1p1=i1+1
    ratio=-(temp(iq))-float(i1-1)
    ans(iq)=a1(i1)+ratio*( a1(i1p1)-a1(i1) )
  end do  
END FUNCTION parama1_v

!----------------------------------------------------------------
REAL FUNCTION parama2(temp)
!----------------------------------------------------------------
  IMPLICIT NONE
  !----------------------------------------------------------------
  !  This program calculate the parameter for crystal growth rate
  !  in Bergeron process
  !----------------------------------------------------------------

  REAL, INTENT (IN   )   :: temp
  REAL, DIMENSION(32)    :: a2
  INTEGER                :: i1, i1p1
  REAL                   :: ratio

  data a2/0.0100,0.4006,0.4831,0.5320,0.5307,0.5319,0.5249, &
              0.4888,0.3849,0.4047,0.4318,0.4771,0.5183,0.5463, &
              0.5651,0.5813,0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506,0.4483,0.4460, &
              0.4433,0.4413,0.4382,0.4361/
  i1=int(-temp)+1
  i1p1=i1+1
  ratio=-(temp)-float(i1-1)
  parama2=a2(i1)+ratio*( a2(i1p1)-a2(i1) )

END FUNCTION parama2

!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
REAL FUNCTION RSLF(P,T)
  IMPLICIT NONE
  REAL, INTENT(IN):: P, T
  REAL:: ESL,X
  REAL, PARAMETER:: C0= .611583699E03
  REAL, PARAMETER:: C1= .444606896E02
  REAL, PARAMETER:: C2= .143177157E01
  REAL, PARAMETER:: C3= .264224321E-1
  REAL, PARAMETER:: C4= .299291081E-3
  REAL, PARAMETER:: C5= .203154182E-5
  REAL, PARAMETER:: C6= .702620698E-8
  REAL, PARAMETER:: C7= .379534310E-11
  REAL, PARAMETER:: C8=-.321582393E-13

  X=MAX(-80.,T-273.16)

  ! ESL=612.2*EXP(17.67*X/(T-29.65))
  ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
  RSLF=.622*ESL/(P-ESL)
END FUNCTION RSLF

!
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
REAL FUNCTION RSIF(P,T)

  IMPLICIT NONE
  REAL, INTENT(IN):: P, T
  REAL:: ESI,X
  REAL, PARAMETER:: C0= .609868993E03
  REAL, PARAMETER:: C1= .499320233E02
  REAL, PARAMETER:: C2= .184672631E01
  REAL, PARAMETER:: C3= .402737184E-1
  REAL, PARAMETER:: C4= .565392987E-3
  REAL, PARAMETER:: C5= .521693933E-5
  REAL, PARAMETER:: C6= .307839583E-7
  REAL, PARAMETER:: C7= .105785160E-9
  REAL, PARAMETER:: C8= .161444444E-12

  X=MAX(-80.,T-273.16)
  ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
  RSIF=.622*ESI/(P-ESI)

  END FUNCTION RSIF
!+---+-----------------------------------------------------------------+

!----------------------------------------------------------------
REAL FUNCTION ggamma(X)
!----------------------------------------------------------------
  IMPLICIT NONE
  !----------------------------------------------------------------
  REAL, INTENT(IN   ) :: x
  REAL, DIMENSION(8)  :: B
  INTEGER             ::j, K1
  REAL                ::PF, G1TO2 ,TEMP

  DATA B/-.577191652,.988205891,-.897056937,.918206857,  &
             -.756704078,.482199394,-.193527818,.035868343/

  PF=1.
  TEMP=X
  DO 10 J=1,200
    IF (TEMP .LE. 2) GO TO 20
      TEMP=TEMP-1.
   10 PF=PF*TEMP
   !  100 FORMAT(//,5X,'module_mp_lin: INPUT TO GAMMA FUNCTION TOO LARGE, X=',E12.5)
   !      WRITE(wrf_err_message,100)X
   !      CALL wrf_error_fatal(wrf_err_message)
   20 G1TO2=1.
      TEMP=TEMP - 1.
      DO 30 K1=1,8
   30 G1TO2=G1TO2 + B(K1)*TEMP**K1
      ggamma=PF*G1TO2

END FUNCTION ggamma

!----------------------------------------------------------------

END MODULE module_mp_sbu_ylin

!WRF:MODEL_LAYER:PHYSICS

!--- The code is based on Lin and Colle (A New Bulk Microphysical Scheme 
!             that Includes Riming Intensity and Temperature Dependent Ice Characteristics, 2011, MWR)
