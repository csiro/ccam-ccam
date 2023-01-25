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
  real, parameter, private :: Nt_c = 250.E6 
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

SUBROUTINE clphy1d_ylin(dt2D, imax2D,                             &
                      qvz2D, qlz2D, qrz2D, qiz2D, qsz2D,          &
                      thz2D, tothz2D, rho2D, orho2D, sqrho2D,     &
                      prez2D, zz2D, dzw2D, zsfc2D,                &
                      precrz2D, preciz2D, precsz2D,               & !zdc20220116
                      EFFC1D2D, EFFI1D2D, EFFS1D2D, EFFR1D2D,     & !zdc 20220208
                      pptrain1D, pptsnow1D,pptice1D,              &
                      kts2D, kte2D, i2D, j2D, riz2D,              &
                      ncz2D, nrz2D, niz2D, nsz2D,                 &
                      fluxr2D, fluxi2D, fluxs2D, fluxg2D, fluxm2D,&
                      fluxf2D, fevap2D, fsubl2D, fauto2D, fcoll2D,&
                      faccr2D, vi2D, vs2D, vg2D,                  &
                      zpsnow2D,zpsaut2D,zpsfw2D,zpsfi2D,zpraci2D, & !process rate 
                      zpiacr2D,zpsaci2D,zpsacw2D,zpsdep2D,        &
                      zpssub2D,zpracs2D,zpsacr2D,zpsmlt2D,        &
                      zpsmltevp2D,zprain2D,zpraut2D,zpracw2D,     &
                      zprevp2D,zpgfr2D,zpvapor2D,zpclw2D,         &
                      zpladj2D,zpcli2D,zpimlt2D,zpihom2D,         &
                      zpidw2D,zpiadj2D,zqschg2D,                  &
                      zdrop2D,lin_aerosolmode2D)                    !aerosol feedback

!-----------------------------------------------------------------------
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
  integer,                             intent(in)    :: kts2D, kte2D, i2D, j2D
  integer,                             intent(in)    :: imax2D
  integer,                             intent(in)    :: lin_aerosolmode2D
  real,                                intent(in)    :: dt2D
  real, dimension(1:imax2D),           intent(in)    :: zsfc2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(in)  :: zdrop2D, riz2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(in)  :: tothz2D,rho2D,orho2D,sqrho2D,   &
                                                          prez2D,zz2D,dzw2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(out) :: precrz2D,preciz2D,precsz2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(out) :: EFFC1D2D,EFFI1D2D,              &
                                                          EFFS1D2D,EFFR1D2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(out) :: fluxr2D,fluxi2D,fluxs2D,fluxg2D,&
                                                          fluxm2D,fluxf2D,fevap2D,fsubl2D,&
                                                          fauto2D,fcoll2D,faccr2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(out) :: vi2D,vs2D,vg2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(out) :: zpsnow2D,zpsaut2D,zpsfw2D,      &
                                                          zpsfi2D,zpraci2D,zpiacr2D,      &
                                                          zpsaci2D,zpsacw2D,zpsdep2D,     &
                                                          zpssub2D,zpracs2D,zpsacr2D,     &
                                                          zpsmlt2D,zpsmltevp2D,zprain2D,  &
                                                          zpraut2D,zpracw2D,zprevp2D,     &
                                                          zpgfr2D,zpvapor2D,zpclw2D,      &
                                                          zpladj2D,zpcli2D,zpimlt2D,      &
                                                          zpihom2D,zpidw2D,zpiadj2D,      &
                                                          zqschg2D
  real, dimension(1:imax2D),             intent(inout) :: pptrain1D, pptsnow1D, pptice1D
  real, dimension(1:imax2D,kts2D:kte2D), intent(inout) :: qvz2D,qlz2D,qrz2D,qiz2D,qsz2D,  &
                                                          thz2D
  real, dimension(1:imax2D,kts2D:kte2D), intent(inout) :: ncz2D,niz2D,nrz2D,nsz2D
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer                               :: kts, kte
  integer                               :: i, j, is, ie, iq, tile, imax
  real, dimension(kts2D:kte2D)          :: zdrop
  integer                               :: lin_aerosolmode
  integer                               :: cnt_sny 
  real, dimension(1:imax2D,kts2D:kte2D) :: nczodt2D,nizodt2D,nrzodt2D,nszodt2D
  real, dimension(1:imax2D,kts2D:kte2D) :: oprez2D,tem2D,temcc2D,qswz2D,qsiz2D,  &
                                             theiz2D,qvoqswz2D,qvoqsiz2D,qvzodt2D, &
                                             qlzodt2D,qizodt2D,qszodt2D,qrzodt2D
  real, dimension(1:imax2D)             :: es1D 
  real, dimension(kts2D:kte2D)          :: qvz,qlz,qrz,qiz,qsz,thz 
  real, dimension(kts2D:kte2D)          :: tothz,rho,orho,sqrho,prez,zz,dzw

  !zdc 20220116
  real, dimension( kts2D:kte2D )        :: precrz, preciz, precsz    
  real, dimension( kts2D:kte2D )        :: EFFC1D, EFFI1D, EFFS1D, EFFR1D
  real, dimension( kts2D:kte2D)         :: fluxr,fluxi,fluxs,fluxg,fluxm, &
                                           fluxf,fevap,fsubl,fauto,fcoll, &
                                           faccr
  real, dimension( kts2D:kte2D)         :: vi, vs, vg        
  real, dimension( kts2D:kte2D)         :: zpsnow,zpsaut,zpsfw,zpsfi,zpraci, & !process rate to understand cloud microphysic
                                           zpiacr,zpsaci,zpsacw,zpsdep,      &
                                           zpssub,zpracs,zpsacr,zpsmlt,      &
                                           zpsmltevp,zprain,zpraut,zpracw,   &
                                           zprevp,zpgfr,zpvapor,zpclw,       &
                                           zpladj,zpcli,zpimlt,zpihom,       &
                                           zpidw,zpiadj,zqschg
  real                                  :: pptrain, pptsnow, pptice
  real                                  :: dt, zsfc
  ! local vars
  real                                  :: obp4, bp3, bp5, bp6, odp4,  &
                                           dp3, dp5, dp5o2
  ! temperary vars
  real                                  :: tmp, tmp0, tmp1, tmp2,tmp3,  &
                                           tmp4, tmpa,tmpb,tmpc,tmpd,alpha1,  &
                                           qic, abi,abr, abg, odtberg,  &
                                           vti50,eiw,eri,esi,esr, esw,  &
                                           erw,delrs,term0,term1,       &
                                           Ap, Bp,                      &
                                           factor, tmp_r, tmp_s,tmp_g,  &
                                           qlpqi, rsat, a1, a2, xnin
  real, dimension( kts2D:kte2D )        :: tmp1D
  real, dimension(1:imax2D,kts2D:kte2D) :: tmp2D
  real, dimension( kts2D:kte2D )        :: oprez, tem, temcc, theiz, qswz,    &
                                           qsiz, qvoqswz, qvoqsiz, qvzodt,    &
                                           qlzodt, qizodt, qszodt, qrzodt
!--- microphysical processes
  real, dimension( kts2D:kte2D )        :: psnow, psaut, psfw,  psfi,  praci,  &
                                           piacr, psaci, psacw, psdep, pssub,  &
                                           pracs, psacr, psmlt, psmltevp,      &
                                           prain, praut, pracw, prevp, pvapor, &
                                           pclw,  pladj, pcli,  pimlt, pihom,  &
                                           pidw,  piadj, pgfr,                 &
                                           qschg, pracis
  real, dimension(1:imax2D,kts2D:kte2D) :: psnow2D, psaut2D, psfw2D,  psfi2D,  praci2D,  &
                                           piacr2D, psaci2D, psacw2D, psdep2D, pssub2D,  &
                                           pracs2D, psacr2D, psmlt2D, psmltevp2D,      &
                                           prain2D, praut2D, pracw2D, prevp2D, pvapor2D, &
                                           pclw2D,  pladj2D, pcli2D,  pimlt2D, pihom2D,  &
                                           pidw2D,  piadj2D, pgfr2D,                 &
                                           qschg2D, pracis2D
  real, dimension(kts2D:kte2D)          :: qvsbar,rs0,viscmu,visc,diffwv,schmidt,xka
  real, dimension(1:imax2D,kts2D:kte2D) :: qvsbar2D,rs02D,viscmu2D,visc2D,    &
                                           diffwv2D,schmidt2D,xka2D
!---- new snow parameters
  real, dimension( kts2D:kte2D )        :: ab_s,ab_r,ab_riming 
  real, dimension( kts2D:kte2D )        :: cap_s       !---- capacitance of snow
  real, dimension(1:imax2D,kts2D:kte2D) :: cap_s2D
  real                                  :: vf1s = 0.65, vf2s = 0.44, vf1r =0.78 , vf2r = 0.31 
  real                                  :: am_c1=0.004, am_c2= 6e-5,    am_c3=0.15
  real                                  :: bm_c1=1.85,  bm_c2= 0.003,   bm_c3=1.25
  real                                  :: aa_c1=1.28,  aa_c2= -0.012,  aa_c3=-0.6
  real                                  :: ba_c1=1.5,   ba_c2= 0.0075,  ba_c3=0.5
  real                                  :: best_a=1.08 ,  best_b = 0.499
  real, dimension(kts2D:kte2D)          :: am_s,bm_s,av_s,bv_s,Ri,tmp_ss,lams 
  real, dimension(1:imax2D,kts2D:kte2D) :: am_s2D,bm_s2D,av_s2D,bv_s2D,Ri2D,tmp_ss2D,lams2D 
  real, dimension(kts2D:kte2D)          :: aa_s,ba_s,tmp_sa 
  real, dimension(1:imax2D,kts2D:kte2D) :: aa_s2D,ba_s2D,tmp_sa2D 
  real                                  :: mu_s=0.,mu_i=0.,mu_r=0.
   
  real                                  :: tc0, disp, Dc_liu, eta, mu_c, R6c      !--- for Liu's autoconversion
  real, dimension(1:imax2D)             :: tc01D, disp1D, Dc_liu1D, eta1D, mu_c1D, R6c1D
  ! Adding variable Riz, which will duplicate Ri but be a copy passed upward
  real, dimension(kts2D:kte2D)          :: Riz
  real, dimension(1:imax2D,kts2D:kte2D) :: Riz2D_lc
  real, dimension( kts2D:kte2D )        :: vtr, vts,                       &
                                           vtrold, vtsold, vtiold,         &
                                           xlambdar, xlambdas,             &
                                           olambdar, olambdas
  real, dimension(1:imax2D,kts2D:kte2D) :: vtr2D, vts2D,                &
                                           vtrold2D, vtsold2D, vtiold2D,&
                                           xlambdar2D, xlambdas2D,      &
                                           olambdar2D, olambdas2D
  !real, dimension( kts2D:kte2D )       :: fluxrain,fluxs,fluxice
  real                                  :: episp0k, dtb, odtb, pi, pio4,       &
                                           pio6, oxLf, xLvocp, xLfocp, av_r,   &
                                           av_i, ocdrag, gambp4, gamdp4,       &
                                           gam4pt5, Cpor, oxmi, gambp3, gamdp3,&
                                           gambp6, gam3pt5, gam2pt75, gambp5o2,&
                                           gamdp5o2, cwoxlf, ocp, xni50, es
  real                                  :: qvmin=1.e-20
  real                                  :: temc1,save1,save2,xni50mx
  ! for terminal velocity flux
  integer                               :: min_q, max_q, max_ri_k, k
  integer, dimension(1:imax2D)          :: min_q1D, max_q1D, max_ri_k1D
  real                                  :: max_ri
  real, dimension(1:imax2D)             :: max_ri1D
  real                                  :: t_del_tv, del_tv, flux, fluxin, fluxout ,tmpqrz
  real, dimension(1:imax2D)             :: t_del_tv1D,del_tv1D,flux1D,fluxin1D,fluxout1D,tmpqrz1D
  logical                               :: notlast
  logical, dimension(1:imax2D)          :: notlast1D, notlast_work
  !
  !zx add
  real, dimension( kts2D:kte2D )        :: ncz,niz,nrz,nsz
  real, dimension( kts2D:kte2D )        :: nczodt, nizodt, nrzodt, nszodt
  real, dimension( kts2D:kte2D )        :: npsaut, npraci, npiacr, npsaci,    &
                                           npsacw, npssub, npsdep, npsacr,    &
                                           npgfr,  npsmlt, npsmltevp,npraut,  &   
                                           npracw, nprevp, nihom,  nimlt,     &
                                           nsagg,  npraut_r
  real, dimension(1:imax2D,kts2D:kte2D) :: npsaut2D, npraci2D, npiacr2D, npsaci2D,    &
                                           npsacw2D, npssub2D, npsdep2D, npsacr2D,    &
                                           npgfr2D,  npsmlt2D, npsmltevp2D,npraut2D,  &
                                           npracw2D, nprevp2D, nihom2D,  nimlt2D,     &
                                           nsagg2D,  npraut_r2D
  real, dimension( kts2D:kte2D )        :: nvtr,   nvts
  real, dimension(1:imax2D,kts2D:kte2D) :: nvtr2D, nvts2D 
  real, dimension( kts2D:kte2D )        :: qisten, qrsten, qssten
  real, dimension( kts2D:kte2D )        :: nisten, nrsten, nssten
  real, dimension(1:imax2D,kts2D:kte2D) :: qisten2D, qrsten2D, qssten2D
  real, dimension(1:imax2D,kts2D:kte2D) :: nisten2D, nrsten2D, nssten2D
  real                                  :: nflux,  nfluxin,nfluxout
  real, dimension(1:imax2D)             :: nflux1D,nfluxin1D,nfluxout1D
  real, dimension( kts2D:kte2D )        :: n0_r,n0_i,n0_c,n0_s                    
  real, dimension( kts2D:kte2D )        :: lami,lamc
  real, dimension(1:imax2D,kts2D:kte2D) :: n0_r2D,n0_i2D,n0_c2D,n0_s2D
  real, dimension(1:imax2D,kts2D:kte2D) :: lami2D,lamc2D
  real                                  :: xmr,xms,xmc,dcs,xmr_i
  real                                  :: lamminr, lammaxr,lammins, lammaxs,lammini, lammaxi
  real                                  :: gambvr1
  real                                  :: lvap
  real, dimension( kts2D:kte2D )        :: nidep, midep
  real, dimension(1:imax2D,kts2D:kte2D) :: nidep2D, midep2D
  real                                  :: mi0

  vtrold=0.
  vtsold=0.
  vtiold=0.

  mu_c    = AMIN1(15., (1000.E6/Nt_c + 2.))
  R6c     = 10.0E-6      !---- 10 micron, threshold radius of cloud droplet
  dtb     = dt2D                                                                         !sny
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
  lamminr = 1./500.E-6
  lamminr = 1./2800.E-6
  lammaxs = 1./10.E-6
  lammins = 1./2000.E-6
  lammaxi = 1./1.E-6
  lammini = 1./(2.*dcs+100.E-6)


  fluxs2D = 0.                                        ! sny
  fluxi2D = 0.
  fluxr2D = 0.
  fluxs2D = 0.

  select case( lin_aerosolmode2D )
    case(0)
      do k=kts2D,kte2D
        !ncz(k) = 250.*1.E6/rho(k)
        ncz2D(1:imax2D,k) = Nt_c/rho2D(1:imax2D,k)
      end do
    case(1)
      do k=kts2D,kte2D
        ncz2D(1:imax2D,k) = zdrop2D(1:imax2D,k)/rho2D(1:imax2D,k)
      end do
    case default
      write(6,*) "ERROR: Unknown option aerosolmode"
      stop
  end select

  kts   = kts2D
  kte   = kte2D
  imax  = imax2D

  do k=kts,kte
    niz2D(1:imax,k) = min(niz2D(1:imax,k),0.3E6/rho2D(1:imax,k))
    nrz2D(1:imax,k) = amax1( 0.0,nrz2D(1:imax,k) )
    nsz2D(1:imax,k) = amax1( 0.0,nsz2D(1:imax,k) )
    nczodt2D(1:imax,k)=amax1( 0.0,odtb*ncz2D(1:imax,k) )
    nizodt2D(1:imax,k)=amax1( 0.0,odtb*niz2D(1:imax,k) )
    nrzodt2D(1:imax,k)=amax1( 0.0,odtb*nrz2D(1:imax,k) )
    nszodt2D(1:imax,k)=amax1( 0.0,odtb*nsz2D(1:imax,k) )
  end do

  do k=kts,kte
    oprez2D(1:imax,k)=1./prez2D(1:imax,k) 
    qlz2D(1:imax,k)  =amax1( 0.0,qlz2D(1:imax,k) )
    qiz2D(1:imax,k)  =amax1( 0.0,qiz2D(1:imax,k) )
    qvz2D(1:imax,k)  =amax1( qvmin,qvz2D(1:imax,k) )
    qsz2D(1:imax,k)  =amax1( 0.0,qsz2D(1:imax,k) )
    qrz2D(1:imax,k)  =amax1( 0.0,qrz2D(1:imax,k) )
    tem2D(1:imax,k)  =thz2D(1:imax,k)*tothz2D(1:imax,k)
    temcc2D(1:imax,k)=tem2D(1:imax,k)-273.15
    es1D(1:imax)     =1000.*svp1*exp( svp2*temcc2D(1:imax,k)/(tem2D(1:imax,k)-svp3) )  !--- RY89 Eq(2.17)
    qswz2D(1:imax,k) =ep2*es1D(1:imax)/(prez2D(1:imax,k)-es1D(1:imax))
 
    do iq=1,imax
      if (tem2D(iq,k) .lt. 273.15 ) then
        es1D(iq)=1000.*svp1*exp( 21.8745584*(tem2D(iq,k)-273.16)/(tem2D(iq,k)-7.66) )
        qsiz2D(iq,k)=ep2*es1D(iq)/(prez2D(iq,k)-es1D(iq))
        if (temcc2D(iq,k) .lt. -40.0) qswz2D(iq,k)=qsiz2D(iq,k)
      else
        qsiz2D(iq,k)=qswz2D(iq,k)
      endif
    enddo

    qvoqswz2D(1:imax,k)  =qvz2D(1:imax,k)/qswz2D(1:imax,k)
    qvoqsiz2D(1:imax,k)  =qvz2D(1:imax,k)/qsiz2D(1:imax,k)
    qvzodt2D(1:imax,k)   =amax1( 0.0,odtb*qvz2D(1:imax,k) )
    qlzodt2D(1:imax,k)   =amax1( 0.0,odtb*qlz2D(1:imax,k) )
    qizodt2D(1:imax,k)   =amax1( 0.0,odtb*qiz2D(1:imax,k) )
    qszodt2D(1:imax,k)   =amax1( 0.0,odtb*qsz2D(1:imax,k) )
    qrzodt2D(1:imax,k)   =amax1( 0.0,odtb*qrz2D(1:imax,k) )
    theiz2D(1:imax,k)=thz2D(1:imax,k)+(xlvocp*qvz2D(1:imax,k)-xlfocp*qiz2D(1:imax,k))/tothz2D(1:imax,k)
  enddo

  do k=kts,kte
    psnow2D(1:imax,k)   =0.                  ! sum all process for snow
    psaut2D(1:imax,k)   =0.                  ! ice crystal aggregation to snow
    psfw2D(1:imax,k)    =0.                  ! BERGERON process to transfer cloud water to snow
    psfi2D(1:imax,k)    =0.                  ! BERGERON process to transfer cloud ice to snow
    praci2D(1:imax,k)   =0.                  ! cloud ice accretion by rain
    piacr2D(1:imax,k)   =0.                  ! rain accretion by cloud ice
    psaci2D(1:imax,k)   =0.                  ! ice crystal accretion by snow
    psacw2D(1:imax,k)   =0.                  ! accretion of cloud water by snow
    psdep2D(1:imax,k)   =0.                  ! deposition of snow
    pssub2D(1:imax,k)   =0.                  ! sublimation of snow (T<0)
    pracs2D(1:imax,k)   =0.                  ! accretion of snow by rain
    psacr2D(1:imax,k)   =0.                  ! accretion of rain by snow
    psmlt2D(1:imax,k)   =0.                  ! melting of snow
    psmltevp2D(1:imax,k)=0.                  ! evaporation of melting snow (T>0)
    prain2D(1:imax,k)   =0.                  ! sum all process for rain
    praut2D(1:imax,k)   =0.                  ! autoconversion of rain
    pracw2D(1:imax,k)   =0.                  ! accretion of cloud water by rain
    prevp2D(1:imax,k)   =0.                  ! evaporation of rain
    pgfr2D(1:imax,k)    =0.                  ! feezing of rain to form graupel (added to PI)
    pvapor2D(1:imax,k)  =0.                  ! sum all process for water vapor to determine qvz
    pclw2D(1:imax,k)    =0.                  ! sum all process for cloud liquid to determine qlz
    pladj2D(1:imax,k)   =0.                  ! saturation adjustment for ql
    pcli2D(1:imax,k)    =0.                  ! sum all process for cloud ice to determine qiz
    pimlt2D(1:imax,k)   =0.                  ! melting of ice crystal >0.
    pihom2D(1:imax,k)   =0.                  ! homogeneous nucleation <-40
    pidw2D(1:imax,k)    =0.                  ! production of cloud ice by BERGERON process
    piadj2D(1:imax,k)   =0.                  ! saturation adjustment for qi
    qschg2D(1:imax,k)   =0.                  ! = psnow / unsure

    npsaut2D(1:imax,k)   =0.
    npraci2D(1:imax,k)   =0.
    npiacr2D(1:imax,k)   =0.
    npsaci2D(1:imax,k)   =0.
    npsacw2D(1:imax,k)   =0.
    npssub2D(1:imax,k)   =0.
    npsdep2D(1:imax,k)   =0.
    npsacr2D(1:imax,k)   =0.
    npgfr2D(1:imax,k)    =0.
    npsmlt2D(1:imax,k)   =0.
    npsmltevp2D(1:imax,k)=0.
    npraut2D(1:imax,k)   =0.
    npracw2D(1:imax,k)   =0.
    nprevp2D(1:imax,k)   =0.

    nimlt2D(1:imax,k)    =0.
    nihom2D(1:imax,k)    =0.
    nsagg2D(1:imax,k)    =0.
    npraut_r2D(1:imax,k) =0.

    n0_i2D(1:imax,k)     =0.
    n0_s2D(1:imax,k)     =0.
    n0_r2D(1:imax,k)     =0.
    n0_c2D(1:imax,k)     =0.
    lamc2D(1:imax,k)     =0.
    lami2D(1:imax,k)     =0.
    xlambdar2D(1:imax,k) =0.
    xlambdas2D(1:imax,k) =0.
    vtr2D(1:imax,k)      =0.
    vts2D(1:imax,k)      =0.
    vtiold2D(1:imax,k)   =0.
    nvtr2D(1:imax,k)     =0.
    nvts2D(1:imax,k)     =0.

    qisten2D(1:imax,k)   =0.
    qrsten2D(1:imax,k)   =0.
    qssten2D(1:imax,k)   =0.
    nisten2D(1:imax,k)   =0.
    nrsten2D(1:imax,k)   =0.
    nssten2D(1:imax,k)   =0.
    nidep2D(1:imax,k)    =0.
    midep2D(1:imax,k)    =0.
  end do

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
    viscmu2D(1:imax,k)=avisc*tem2D(1:imax,k)**1.5/(tem2D(1:imax,k)+120.0)
    visc2D(1:imax,k)=viscmu2D(1:imax,k)*orho2D(1:imax,k)
    diffwv2D(1:imax,k)=adiffwv*tem2D(1:imax,k)**1.81*oprez2D(1:imax,k)
    schmidt2D(1:imax,k)=visc2D(1:imax,k)/diffwv2D(1:imax,k)
    xka2D(1:imax,k)=axka*viscmu2D(1:imax,k)
    rs02D(1:imax,k)=ep2*1000.*svp1/(prez2D(1:imax,k)-1000.*svp1)
  end do

  ! ---- YLIN, set snow variables
  !
  !---- A+B in depositional growth, the first try just take from Rogers and Yau(1989)
  !         ab_s(k) = lsub*lsub*orv/(tcond(k)*temp(k))+&
  !                   rv*temp(k)/(diffu(k)*qvsi(k))

  do k = kts, kte
    tc01D(1:imax)   = tem2D(1:imax,k)-273.15
    do iq=1,imax
      if (rho2D(iq,k)*qlz2D(iq,k) .gt. 1e-5 .AND. rho2D(iq,k)*qsz2D(iq,k) .gt. 1e-5) then
        Ri2D(iq,k) = 1.0/(1.0+6e-5/(rho2D(iq,k)**1.170*qlz2D(iq,k)*qsz2D(iq,k)**0.170))
      else
        Ri2D(iq,k) = 0
      end if
    end do
  end do

  !
  !--- make sure Ri does not decrease downward in a column
  !
  max_ri_k1D(:) = maxloc(Ri2D,dim=2)
  do iq = 1, imax
    max_ri1D(iq)   = Ri2D(iq,max_ri_k1D(iq))
  end do

  do k = kts,kte
    do iq = 1,imax
      if ( k<=max_ri_k1D(iq) ) then
        Ri2D(iq,k) = max_ri1D(iq)
      end if
    end do
  end do

  !--- YLIN, get PI properties
  do k = kts, kte
    Ri2D(1:imax,k) = AMAX1(0.,AMIN1(Ri2D(1:imax,k),1.0))
    ! Store the value of Ri(k) as Riz(k)
    Riz2D_lc(1:imax,k) = Ri2D(1:imax,k)

    cap_s2D(1:imax,k)= 0.25*(1+Ri2D(1:imax,k))
    tc01D(1:imax)    = AMIN1(-0.1, tem2D(1:imax,k)-273.15)
    n0_s2D(1:imax,k) = amin1(2.0E8, 2.0E6*exp(-0.12*tc01D(1:imax)))
    am_s2D(1:imax,k) = am_c1+am_c2*tc01D(1:imax)+am_c3*Ri2D(1:imax,k)*Ri2D(1:imax,k)   !--- Heymsfield 2007
    am_s2D(1:imax,k) = AMAX1(0.000023,am_s2D(1:imax,k))                                !--- use the a_min in table 1 of Heymsfield
    bm_s2D(1:imax,k) = bm_c1+bm_c2*tc01D(1:imax)+bm_c3*Ri2D(1:imax,k)
    bm_s2D(1:imax,k) = AMIN1(bm_s2D(1:imax,k),3.0)                                     !---- capped by 3
    !--  converting from cgs to SI unit
    am_s2D(1:imax,k) =  10**(2*bm_s2D(1:imax,k)-3.0)*am_s2D(1:imax,k)
    aa_s2D(1:imax,k) = aa_c1 + aa_c2*tc01D(1:imax) + aa_c3*Ri2D(1:imax,k)
    ba_s2D(1:imax,k) = ba_c1 + ba_c2*tc01D(1:imax) + ba_c3*Ri2D(1:imax,k)
    !--  convert to SI unit as in paper
    aa_s2D(1:imax,k) = (1e-2)**(2.0-ba_s2D(1:imax,k))*aa_s2D(1:imax,k)
    !---- get v from Mitchell 1996
    av_s2D(1:imax,k) = best_a*viscmu2D(1:imax,k)*(2*grav*am_s2D(1:imax,k)/rho2D(1:imax,k)/ &
                       aa_s2D(1:imax,k)/(viscmu2D(1:imax,k)**2))**best_b
    bv_s2D(1:imax,k) = best_b*(bm_s2D(1:imax,k)-ba_s2D(1:imax,k)+2)-1

    tmp_ss2D(1:imax,k)= bm_s2D(1:imax,k)+mu_s+1
    tmp_sa2D(1:imax,k)= ba_s2D(1:imax,k)+mu_s+1
  end do


  !***********************************************************************
  ! Calculate precipitation fluxes due to terminal velocities.
  !***********************************************************************
  !
  !- Calculate termianl velocity (vt?)  of precipitation q?z
  !- Find maximum vt? to determine the small delta t
  !
  !-- rain
  !       CALL wrf_debug ( 100 , 'module_ylin, start precip fluxes' )

  t_del_tv1D(1:imax)=0.
  del_tv1D(1:imax)=dtb
  notlast1D(1:imax)=.true.
  DO while (any(notlast1D))

    notlast_work = notlast1D

    min_q1D(1:imax)=kte
    max_q1D(1:imax)=kts-1

    ! if no rain, --> minq>maxq --> notlast=False (only case minq>maxq)
    ! if rain --> minq<maxq (always), some vertical points norain--> no lamda, velocity
    do k=kts,kte-1
      do iq = 1,imax
        if (notlast1D(iq)) then
          if (qrz2D(iq,k) .gt. 1.0e-8) then
            min_q1D(iq)=min0(min_q1D(iq),k)
            max_q1D(iq)=max0(max_q1D(iq),k)
            ! tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrz(k))
            ! tmp1=sqrt(tmp1)
            ! vtrold(k)=o6*av_r*gambp4*sqrho(k)/tmp1**bv_r
            xlambdar2D(iq,k)=(pi*rhowater*nrz2D(iq,k)/qrz2D(iq,k))**(1./3.)   !zx
            n0_r2D(iq,k)=nrz2D(iq,k)*xlambdar2D(iq,k)
            if (xlambdar2D(iq,k).lt.lamminr) then
              xlambdar2D(iq,k) = lamminr
              n0_r2D(iq,k) = xlambdar2D(iq,k)**4*qrz2D(iq,k)/(pi*rhowater)
              nrz2D(iq,k) = n0_r2D(iq,k)/xlambdar2D(iq,k) 
            else if (xlambdar2D(iq,k).gt.lammaxr) then
              xlambdar2D(iq,k) = lammaxr
              n0_r2D(iq,k) = xlambdar2D(iq,k)**4*qrz2D(iq,k)/(pi*rhowater)
              nrz2D(iq,k) = n0_r2D(iq,k)/xlambdar2D(iq,k)
            end if
            olambdar2D(iq,k)=1.0/xlambdar2D(iq,k)
            vtrold2D(iq,k)=o6*av_r*gambp4*sqrho2D(iq,k)*olambdar2D(iq,k)**bv_r
            nvtr2D(iq,k)=av_r*gambvr1*sqrho2D(iq,k)*olambdar2D(iq,k)**bv_r
            
            if (k .eq. 1) then
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zsfc2D(iq))/vtrold2D(iq,k))
            else
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zz2D(iq,k-1))/vtrold2D(iq,k))
            endif
          else
            vtrold2D(iq,k)=0.
            nvtr2D(iq,k)=0.
            olambdar2D(iq,k)=0.
          endif
        endif     ! notlast1D
      enddo       ! iq
    enddo         ! k

    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)

    do iq = 1,imax
      if (notlast1D(iq)) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then
          t_del_tv1D(iq)=t_del_tv1D(iq)+del_tv1D(iq)
          if ( t_del_tv1D(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv1D(iq)=dtb+del_tv1D(iq)-t_del_tv1D(iq)
          end if
          fluxin1D(iq)=0.
          nfluxin1D(iq)=0. ! sny
        end if ! maxq>minq
      end if   ! notlast
    end do     ! iq

    do k = kte,kts,-1
    !do k = maxval(max_q1D),minval(min_q1D),-1
      do iq = 1,imax
        if ( notlast1D(iq) ) then
          !if (max_q1D(iq) .ge. min_q1D(iq)) then      
          if ( k>=min_q1D(iq) .and. k<=max_q1D(iq) ) then
            fluxout1D(iq)=rho2D(iq,k)*vtrold2D(iq,k)*qrz2D(iq,k)
            flux1D(iq)=(fluxin1D(iq)-fluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            tmpqrz1D(iq)=qrz2D(iq,k)
            qrz2D(iq,k)=qrz2D(iq,k)+del_tv1D(iq)*flux1D(iq)
            fluxin1D(iq)=fluxout1D(iq)

            nfluxout1D(iq)=rho2D(iq,k)*nvtr2D(iq,k)*nrz2D(iq,k)
            nflux1D(iq)=(nfluxin1D(iq)-nfluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            nrz2D(iq,k)=nrz2D(iq,k)+del_tv1D(iq)*nflux1D(iq)
            nfluxin1D(iq)=nfluxout1D(iq)
            qrsten2D(iq,k)=flux1D(iq)
            nrsten2D(iq,k)=nflux1D(iq)
 
            fluxr2D(iq,k) = fluxout1D(iq)                     ! sny
          end if !maxq, minq
        end if   !notlast
      end do     !iq
    end do       !k

    do iq = 1, imax
      if ( notlast1D(iq) ) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then      
          if (min_q1D(iq) .eq. 1) then
            pptrain1D(iq)=pptrain1D(iq)+fluxin1D(iq)*del_tv1D(iq)
          else
            qrz2D(iq,min_q1D(iq)-1)=qrz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                           fluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
            nrz2D(iq,min_q1D(iq)-1)=nrz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                           nfluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
          endif  !minq
        else
          notlast_work(iq)=.false.
        end if 
      end if   ! notlast
    end do     !iq


    notlast1D(:) = notlast_work(:)
  ENDDO      ! while(any(notlast))


  !
  !-- snow
  !
  t_del_tv1D(1:imax)=0.
  del_tv1D(1:imax)=dtb
  notlast1D(1:imax)=.true.
  DO while (any(notlast1D))

    notlast_work = notlast1D
    min_q1D(1:imax)=kte
    max_q1D(1:imax)=kts-1

    do k=kts,kte-1
      do iq = 1, imax
        if (notlast1D(iq)) then
          if (qsz2D(iq,k) .gt. 1.0e-8) then
            min_q1D(iq)=min0(min_q1D(iq),k)
            max_q1D(iq)=max0(max_q1D(iq),k)
            !   tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
            !   **(1./tmp_ss(k))
            !   vtsold(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
            !   ggamma(tmp_ss(k))/(tmp1**bv_s(k))
            ! Zhao 2022 - Row 2 Table 2 or Lin 2011 - Formula A3
            xlambdas2D(iq,k)=(am_s2D(iq,k)*ggamma(tmp_ss2D(iq,k))*nsz2D(iq,k)/qsz2D(iq,k))**(1./bm_s2D(iq,k))
            ! Zhao 2022 - Row 1 Table 2
            n0_s2D(iq,k)=nsz2D(iq,k)*xlambdas2D(iq,k)
            if (xlambdas2D(iq,k).lt.lammins) then
              xlambdas2D(iq,k)= lammins
              n0_s2D(iq,k) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*qsz2D(iq,k)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
              nsz2D(iq,k) = n0_s2D(iq,k)/xlambdas2D(iq,k)
            else if (xlambdas2D(iq,k).gt.lammaxs) then
              xlambdas2D(iq,k) = lammaxs
              n0_s2D(iq,k) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*qsz2D(iq,k)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
              nsz2D(iq,k) = n0_s2D(iq,k)/xlambdas2D(iq,k)
            end if
            olambdas2D(iq,k)=1.0/xlambdas2D(iq,k)
            ! Zhao 2022 - Row 3 Table 2
            vtsold2D(iq,k)= sqrho2D(iq,k)*av_s2D(iq,k)*ggamma(bv_s2D(iq,k)+tmp_ss2D(iq,k))/ &
               ggamma(tmp_ss2D(iq,k))*(olambdas2D(iq,k)**bv_s2D(iq,k))
            ! Zhao 2022 - Row 4 Table 2
            nvts2D(iq,k)=sqrho2D(iq,k)*av_s2D(iq,k)*ggamma(bv_s2D(iq,k)+1)*(olambdas2D(iq,k)**bv_s2D(iq,k))
            if (k .eq. 1) then
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zsfc2D(iq))/vtsold2D(iq,k))
            else
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zz2D(iq,k-1))/vtsold2D(iq,k))
            endif
          else
            vtsold2D(iq,k)=0.
            nvts2D(iq,k)=0.
            olambdas2D(iq,k)=0.
          endif
        endif   ! notlast1D
      enddo     ! iq
    enddo       ! k

    !
    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)

    do iq = 1,imax
      if (notlast1D(iq)) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then
          t_del_tv1D(iq)=t_del_tv1D(iq)+del_tv1D(iq)
          if ( t_del_tv1D(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv1D(iq)=dtb+del_tv1D(iq)-t_del_tv1D(iq)
          endif
          fluxin1D(iq) = 0.
          nfluxin1D(iq) = 0.
        end if ! maxq>minq
      end if   ! notlast
    end do     ! iq

    do k = kte,kts,-1
      do iq = 1,imax
        if ( notlast1D(iq) ) then
          if ( k>=min_q1D(iq) .and. k<=max_q1D(iq) ) then
            fluxout1D(iq)=rho2D(iq,k)*vtsold2D(iq,k)*qsz2D(iq,k)
            flux1D(iq)=(fluxin1D(iq)-fluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            qsz2D(iq,k)=qsz2D(iq,k)+del_tv1D(iq)*flux1D(iq)
            qsz2D(iq,k)=amax1(0.,qsz2D(iq,k))
            fluxin1D(iq)=fluxout1D(iq)

            nfluxout1D(iq)=rho2D(iq,k)*nvts2D(iq,k)*nsz2D(iq,k)
            nflux1D(iq)   =(nfluxin1D(iq)-nfluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            nsz2D(iq,k)  =nsz2D(iq,k)+del_tv1D(iq)*nflux1D(iq)
            nfluxin1D(iq) =nfluxout1D(iq)
            qssten2D(iq,k)=flux1D(iq)
            nssten2D(iq,k)=nflux1D(iq)

            fluxs2D(iq,k) = fluxout1D(iq)                     ! sny
          end if ! maxq, minq
        end if   ! notlast
      end do     ! iq
    end do       ! k

    do iq = 1, imax
      if ( notlast1D(iq) ) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then
          if (min_q1D(iq) .eq. 1) then
            pptsnow1D(iq)=pptsnow1D(iq)+fluxin1D(iq)*del_tv1D(iq)
          else
            qsz2D(iq,min_q1D(iq)-1)=qsz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                     fluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
            nsz2D(iq,min_q1D(iq)-1)=nsz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                       nfluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
          endif ! minq
        else
          notlast_work(iq)=.false.
        endif
      endif     ! notlast
    enddo       ! iq

    notlast1D(:) = notlast_work(:)
  ENDDO       ! while(any(notlast))

  !
  !-- cloud ice  (03/21/02) using Heymsfield and Donner (1990) Vi=3.29*qi^0.16
  !
  t_del_tv1D(1:imax)=0.
  del_tv1D(1:imax)=dtb
  notlast1D(1:imax)=.true.
  DO while (any(notlast1D))

    notlast_work = notlast1D
    min_q1D(1:imax)=kte
    max_q1D(1:imax)=kts-1

    do k=kts,kte-1
      do iq = 1, imax
        if (notlast1D(iq)) then
          if (qiz2D(iq,k) .gt. 1.0e-8) then
            min_q1D(iq)=min0(min_q1D(iq),k)
            max_q1D(iq)=max0(max_q1D(iq),k)
            vtiold2D(iq,k)= 3.29 * (rho2D(iq,k)* qiz2D(iq,k))** 0.16  ! Heymsfield and Donner
            if (k .eq. 1) then
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zsfc2D(iq))/vtiold2D(iq,k))
            else
              del_tv1D(iq)=amin1(del_tv1D(iq),0.9*(zz2D(iq,k)-zz2D(iq,k-1))/vtiold2D(iq,k))
            endif
          else
            vtiold2D(iq,k)=0.
          endif
        endif   ! notlast1D
      enddo     ! iq
    enddo       ! k

    !
    !- Check if the summation of the small delta t >=  big delta t
    !             (t_del_tv)          (del_tv)             (dtb)
    do iq = 1,imax
      if (notlast1D(iq)) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then
          t_del_tv1D(iq)=t_del_tv1D(iq)+del_tv1D(iq)
          if ( t_del_tv1D(iq) .ge. dtb ) then
            notlast_work(iq)=.false.
            del_tv1D(iq)=dtb+del_tv1D(iq)-t_del_tv1D(iq)
          endif
          fluxin1D(iq) = 0.
          nfluxin1D(iq) = 0.
        end if   ! maxq>minq
      end if     ! notlast
    end do       ! iq

    do k = kte,kts,-1
      do iq = 1,imax
        if ( notlast1D(iq) ) then
          if ( k>=min_q1D(iq) .and. k<=max_q1D(iq) ) then
            fluxout1D(iq)=rho2D(iq,k)*vtiold2D(iq,k)*qiz2D(iq,k)
            flux1D(iq)=(fluxin1D(iq)-fluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            qiz2D(iq,k)=qiz2D(iq,k)+del_tv1D(iq)*flux1D(iq)
            qiz2D(iq,k)=amax1(0.,qiz2D(iq,k))
            fluxin1D(iq)=fluxout1D(iq)

            nfluxout1D(iq)=rho2D(iq,k)*vtiold2D(iq,k)*niz2D(iq,k)
            nflux1D(iq)=(nfluxin1D(iq)-nfluxout1D(iq))/rho2D(iq,k)/dzw2D(iq,k)
            niz2D(iq,k)=niz2D(iq,k)+del_tv1D(iq)*nflux1D(iq)
            niz2D(iq,k)=amax1(0.,niz2D(iq,k))
            nfluxin1D(iq)=nfluxout1D(iq)
            qisten2D(iq,k)=flux1D(iq)
            nisten2D(iq,k)=nflux1D(iq)

            fluxi2D(iq,k) = fluxout1D(iq)                     ! sny
          end if ! maxq, minq
        end if   ! notlast
      end do     ! iq
    end do       ! k

    do iq = 1, imax
      if ( notlast1D(iq) ) then
        if (max_q1D(iq) .ge. min_q1D(iq)) then
          if (min_q1D(iq) .eq. 1) then
            pptice1D(iq)=pptice1D(iq)+fluxin1D(iq)*del_tv1D(iq)
          else
            qiz2D(iq,min_q1D(iq)-1)=qiz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                         fluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
            niz2D(iq,min_q1D(iq)-1)=niz2D(iq,min_q1D(iq)-1)+del_tv1D(iq)*  &
                         nfluxin1D(iq)/rho2D(iq,min_q1D(iq)-1)/dzw2D(iq,min_q1D(iq)-1)
          endif
        else
          notlast_work(iq)=.false.
        end if
      end if      ! notlast
    end do        ! iq
    
    notlast1D(:) = notlast_work(:)
  ENDDO       ! while(any(notlast))

  ! zdc 20220116
  do k=kts,kte-1                         !sg beg
    precrz2D(1:imax,k)=rho2D(1:imax,k)*vtrold2D(1:imax,k)*qrz2D(1:imax,k)
    preciz2D(1:imax,k)=rho2D(1:imax,k)*vtiold2D(1:imax,k)*qiz2D(1:imax,k)
    precsz2D(1:imax,k)=rho2D(1:imax,k)*vtsold2D(1:imax,k)*qsz2D(1:imax,k)
  enddo                                  !sg end
    precrz2D(1:imax,kte)=0. !wig - top level never set for vtXold vars
    preciz2D(1:imax,kte)=0. !wig
    precsz2D(1:imax,kte)=0. !wig
  !     CALL wrf_debug ( 100 , 'module_ylin: end precip flux' )


  ! Microphpysics processes
  DO k=kts,kte
    qvzodt2D(1:imax,k)=amax1( 0.0,odtb*qvz2D(1:imax,k) )
    qlzodt2D(1:imax,k)=amax1( 0.0,odtb*qlz2D(1:imax,k) )
    qizodt2D(1:imax,k)=amax1( 0.0,odtb*qiz2D(1:imax,k) )
    qszodt2D(1:imax,k)=amax1( 0.0,odtb*qsz2D(1:imax,k) )
    qrzodt2D(1:imax,k)=amax1( 0.0,odtb*qrz2D(1:imax,k) )

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
    tmp2D(1:imax,k)=qiz2D(1:imax,k)+qlz2D(1:imax,k)+qsz2D(1:imax,k)+qrz2D(1:imax,k)
       
    do iq = 1, imax
      if( .not.(qvz2D(iq,k)+qlz2D(iq,k)+qiz2D(iq,k) .lt. qsiz2D(iq,k)  &
            .and. tmp2D(iq,k) .eq. 0.0) ) then !go to 2000
      !
      !! calculate terminal velocity of rain
      !
        if (qrz2D(iq,k) .gt. 1.0e-8) then
        !  tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
        !  xlambdar(k)=sqrt(tmp1)
        !  olambdar(k)=1.0/xlambdar(k)
        !  vtrold(k)=o6*av_r*gambp4*sqrho(k)*olambdar(k)**bv_r
          xlambdar2D(iq,k)=(pi*rhowater*nrz2D(iq,k)/qrz2D(iq,k))**(1./3.)   !zx
          n0_r2D(iq,k)=nrz2D(iq,k)*xlambdar2D(iq,k)
          if (xlambdar2D(iq,k).lt.lamminr) then
            xlambdar2D(iq,k) = lamminr
            n0_r2D(iq,k) = xlambdar2D(iq,k)**4*qrz2D(iq,k)/(pi*rhowater)
            nrz2D(iq,k) = n0_r2D(iq,k)/xlambdar2D(iq,k)
          else if (xlambdar2D(iq,k).gt.lammaxr) then
            xlambdar2D(iq,k) = lammaxr
            n0_r2D(iq,k) = xlambdar2D(iq,k)**4*qrz2D(iq,k)/(pi*rhowater)
            nrz2D(iq,k) = n0_r2D(iq,k)/xlambdar2D(iq,k)
          end if
          olambdar2D(iq,k)=1.0/xlambdar2D(iq,k)
          vtrold2D(iq,k)=o6*av_r*gambp4*sqrho2D(iq,k)*olambdar2D(iq,k)**bv_r
          nvtr2D(iq,k)=av_r*gambvr1*sqrho2D(iq,k)*olambdar2D(iq,k)**bv_r
        else
          vtrold2D(iq,k)=0.
          olambdar2D(iq,k)=0.
          nvtr2D(iq,k)=0.
        end if  ! qrz2D

        if (qrz2D(iq,k) .gt. 1.0e-8) then
        !  tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
        !  xlambdar(k)=sqrt(tmp1)
        !  olambdar(k)=1.0/xlambdar(k)
        !  vtr(k)=o6*av_r*gambp4*sqrho(k)*olambdar(k)**bv_r
        else
        !  vtr(k)=0.
        !  olambdar(k)=0.
        endif
        vtr2D(iq,k)=vtrold2D(iq,k)


        !!
        !!! calculate terminal velocity of snow
        !!
        if (qsz2D(iq,k) .gt. 1.0e-8) then
        !  tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
        !                   **(1./tmp_ss(k))
        !            xlambdas(k)=tmp1
        !            olambdas(k)=1.0/tmp1
        !            vtsold(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
        !                      ggamma(tmp_ss(k))/(tmp1**bv_s(k))
          xlambdas2D(iq,k)=(am_s2D(iq,k)*ggamma(tmp_ss2D(iq,k))*nsz2D(iq,k)/qsz2D(iq,k))**(1./bm_s2D(iq,k))
          n0_s2D(iq,k)=nsz2D(iq,k)*xlambdas2D(iq,k)
          if (xlambdas2D(iq,k).lt.lammins) then
            xlambdas2D(iq,k)= lamminS
            n0_s2D(iq,K) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*qsz2D(iq,K)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
            nsz2D(iq,K) = n0_s2D(iq,K)/xlambdas2D(iq,k)
          else if (xlambdas2D(iq,k).gt.lammaxs) then
            xlambdas2D(iq,k) = lammaxs
            n0_s2D(iq,K) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*qsz2D(iq,K)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
            nsz2D(iq,K) = n0_s2D(iq,K)/xlambdas2D(iq,k)
          end if
          olambdas2D(iq,k)=1.0/xlambdas2D(iq,k)
          vtsold2D(iq,k)= sqrho2D(iq,k)*av_s2D(iq,k)*ggamma(bv_s2D(iq,k)+tmp_ss2D(iq,k))/ &
                 ggamma(tmp_ss2D(iq,k))*(olambdas2D(iq,k)**bv_s2D(iq,k))
          nvts2D(iq,k)=sqrho2D(iq,k)*av_s2D(iq,k)*ggamma(bv_s2D(iq,k)+1)*(olambdas2D(iq,k)**bv_s2D(iq,k))
        else
          vtsold2D(iq,k)=0.
          olambdas2D(iq,k)=0.
          xlambdas2D(iq,k)=0.
          nvts2D(iq,k)=0.
        endif

        if (qsz2D(iq,k) .gt. 1.0e-8) then
        ! tmp1= (am_s(k)*N0_s(k)*ggamma(tmp_ss(k))*orho(k)/qsz(k))&
        !                   **(1./tmp_ss(k))
        !             olambdas(k)=1.0/tmp1
        !             vts(k)= sqrho(k)*av_s(k)*ggamma(bv_s(k)+tmp_ss(k))/ &
        !                      ggamma(tmp_ss(k))/(tmp1**bv_s(k))
        else
        !            vts2D(iq,k)=0.
        !            olambdas2D(iq,k)=0.
        endif
        vts2D(iq,k)=vtsold2D(iq,k)

        !---------- start of snow/ice processes below freezing

        if (tem2D(iq,k) .lt. 273.15) then

        !
        ! ice nucleation, cooper curve

          if ((qvoqswz2D(iq,k).ge.0.999.and.temcc2D(iq,k).le. -8.).or. &
            qvoqsiz2D(iq,k).ge.1.08) then
            nidep2D(iq,k) = 5.*exp(0.304*(273.15-temcc2D(iq,k)))     ! m-3
            nidep2D(iq,k) = min(nidep2D(iq,k), 500.e3)               !5.e8) sny ! limit to 500 L-1
            nidep2D(iq,k) = max(nidep2D(iq,k)/rho2D(iq,k), 0.)       ! convert to kg-1
            nidep2D(iq,k) = (nidep2D(iq,k) - niz2D(iq,k))*odtb
            midep2D(iq,k) = nidep2D(iq,k)*mi0
          end if
          !***********************************************************************
          !*********        snow production processes for T < 0 C       **********
          !***********************************************************************
          !
          ! (1) ICE CRYSTAL AGGREGATION TO SNOW (Psaut): Lin (21)
          !!    psaut=alpha1*(qi-qi0)
          !!    alpha1=1.0e-3*exp(0.025*(T-T0))
          !
          alpha1=1.0e-3*exp( 0.025*temcc2D(iq,k) )

          ! BELOW SECTION TURN OFF BY SONNY   sny: on temp
          ! ---------------------------------------------------------------
          !if(temcc(k) .lt. -20.0) then
          !    tmp1=-7.6+4.0*exp( -0.2443e-3*(abs(temcc(k))-20)**2.455 )
          !    qic=1.0e-3*exp(tmp1)*orho(k)
          !else
          !    qic=qi0
          !end if
          !----------------------------------------------------------------

          qic = qi0  ! sny: OFF temp

          tmp1=odtb*(qiz2D(iq,k)-qic)*(1.0-exp(-alpha1*dtb))
          psaut2D(iq,k)=amax1( 0.0,tmp1 )
          npsaut2D(iq,k)=amax1( 0.0,psaut2D(iq,k)/xms)
          !
          ! (2) BERGERON PROCESS TRANSFER OF CLOUD WATER TO SNOW (Psfw)
          !     this process only considered when -31 C < T < 0 C
          !     Lin (33) and Hsie (17)
          ! 
          !!
          !!    parama1 and parama2 functions must be user supplied
          !!

          if( qlz2D(iq,k) .gt. 1.0e-10 ) then
            temc1=amax1(-30.99,temcc2D(iq,k))
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
            vti50=av_i*di50**bv_i*sqrho2D(iq,k)
            eiw=1.0
            save1=a1*xmi50**a2
            save2=0.25*pi*eiw*rho2D(iq,k)*di50*di50*vti50
            tmp2=( save1 + save2*qlz2D(iq,k) )
            !
            !!  maximum number of 50 micron crystals limited by the amount
            !!  of supercool water
            !
            xni50mx=qlzodt2D(iq,k)/tmp2
            !
            !!   number of 50 micron crystals produced
            !
            xni50=qiz2D(iq,k)*( 1.0-exp(-dtb*odtberg) )/xmi50
            xni50=amin1(xni50,xni50mx)
            !
            tmp3=odtb*tmp2/save2*( 1.0-exp(-save2*xni50*dtb) )
            psfw2D(iq,k)=amin1( tmp3,qlzodt2D(iq,k) )
            !
            ! (3) REDUCTION OF CLOUD ICE BY BERGERON PROCESS (Psfi): Lin (34)
            !     this process only considered when -31 C < T < 0 C
            !
            tmp1=xni50*xmi50-psfw2D(iq,k)
            psfi2D(iq,k)=amin1(tmp1,qizodt2D(iq,k))
          end if
        
          if(qrz2D(iq,k) .gt. 0.0) then  ! go to 1000

          !
          ! Processes (4) and (5) only need when qrz > 0.0
          !
          ! (4) CLOUD ICE ACCRETION BY RAIN (Praci): Lin (25)
          !     produce PI
          !
            eri=1.0
            save1=pio4*eri*n0_r2D(iq,k)*av_r*sqrho2D(iq,k)
            tmp1=save1*gambp3*olambdar2D(iq,k)**bp3
            praci2D(iq,k)=qizodt2D(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npraci2D(iq,k)=niz2D(iq,k)*tmp1

          !
          ! (5) RAIN ACCRETION BY CLOUD ICE (Piacr): Lin (26)
          !
            tmp2=qiz2D(iq,k)*save1*rho2D(iq,k)*pio6*rhowater*gambp6*oxmi* &
                olambdar2D(iq,k)**bp6
            piacr2D(iq,k)=amin1( tmp2,qrzodt2D(iq,k) )
            npiacr2D(iq,k)=pio4*eri*nrz2D(iq,k)*av_r*niz2D(iq,k)*gambp3* &
                olambdar2D(iq,k)**bp3  !--wdm6
          end if !1000    continue

          if(qsz2D(iq,k) .gt. 0.0) then !go to 1200
          !
          ! Compute the following processes only when qsz > 0.0
          !
          !
          ! (6) ICE CRYSTAL ACCRETION BY SNOW (Psaci): Lin (22)
          !
            esi=exp( 0.025*temcc2D(iq,k) )
            save1 = aa_s2D(iq,k)*sqrho2D(iq,k)*N0_s2D(iq,k)* &
            ggamma(bv_s2D(iq,k)+tmp_sa2D(iq,k))*         &
            olambdas2D(iq,k)**(bv_s2D(iq,k)+tmp_sa2D(iq,k))

            tmp1=esi*save1
            psaci2D(iq,k)=qizodt2D(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsaci2D(iq,k)=amin1( tmp1*niz2D(iq,k),nizodt2D(iq,k))
          !
          ! (7) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            tmp1=esw*save1
            psacw2D(iq,k)=qlzodt2D(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsacw2D(iq,k)=amin1(tmp1*ncz2D(iq,k),ncz2D(iq,k))

          ! recalculate the saturatuin temperature
          !
          ! (8) DEPOSITION/SUBLIMATION OF SNOW (Psdep/Pssub): Lin (31)
          !     includes consideration of ventilation effect
          !
            tmpa=rvapor*xka2D(iq,k)*tem2D(iq,k)*tem2D(iq,k)
            tmpb=xls*xls*rho2D(iq,k)*qsiz2D(iq,k)*diffwv2D(iq,k)
            tmpc=tmpa*qsiz2D(iq,k)*diffwv2D(iq,k)
            abi=4.0*pi*cap_s2D(iq,k)*(qvoqsiz2D(iq,k)-1.0)*tmpc/(tmpa+tmpb)
            tmp1=av_s2D(iq,k)*sqrho2D(iq,k)*        &
                olambdas2D(iq,k)**(5+bv_s2D(iq,k)+2*mu_s)/visc2D(iq,k)

          !---- YLIN, here there is some approximation assuming mu_s =1, so gamma(2)=1, etc.

            tmp2= abi*N0_s2D(iq,k)*( vf1s*olambdas2D(iq,k)*olambdas2D(iq,k)+ &
                 vf2s*schmidt2D(iq,k)**0.33334* &
                 ggamma(2.5+0.5*bv_s2D(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=odtb*( qvz2D(iq,k)-qsiz2D(iq,k) )
            tmp3=amin1(tmp3,0.)

            if( tmp2 .le. 0.0) then
              tmp2=amax1( tmp2,tmp3)
              pssub2D(iq,k)=amax1( tmp2,-qszodt2D(iq,k) )
              psdep2D(iq,k)=0.0
            else
              psdep2D(iq,k)=amin1( tmp2,tmp3 )
              pssub2D(iq,k)=0.0
            end if
            if(qsz2D(iq,k) .ge. 0.0) then
              npssub2D(iq,k)=pssub2D(iq,k)*nsz2D(iq,k)/qsz2D(iq,k)
              npsdep2D(iq,k)=npsdep2D(iq,k)*nsz2D(iq,k)/qsz2D(iq,k)
            else
              npssub2D(iq,k)=pssub2D(iq,k)/xms
              npsdep2D(iq,k)=npsdep2D(iq,k)/xms
            end if

            if(qrz2D(iq,k) .gt. 0.0) then !go to 1200
            !
            ! Compute processes (9) and (10) only when qsz > 0.0 and qrz > 0.0
            ! these two terms need to be refined in the future, they should be equal
            !
            ! (9) ACCRETION OF SNOW BY RAIN (Pracs): Lin (27)
            !
              esr=1.0
              tmpa=olambdar2D(iq,k)*olambdar2D(iq,k)
              tmpb=olambdas2D(iq,k)*olambdas2D(iq,k)
              tmpc=olambdar2D(iq,k)*olambdas2D(iq,k)
              tmp1=pi*pi*esr*n0_r2D(iq,k)*N0_s2D(iq,k)*    &
                        abs( vtr2D(iq,k)-vts2D(iq,k) )*orho2D(iq,k)
            ! tmp1=pi*pi*esr*n0_r(k)*N0_s(k)*            &
            ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5*orho(k)
              tmp2=tmpb*tmpb*olambdar2D(iq,k)*(5.0*tmpb+2.0*tmpc+0.5*tmpa)
              tmp3=tmp1*rhosnow*tmp2
              pracs2D(iq,k)=amin1( tmp3,qszodt2D(iq,k) )
            !
            ! (10) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
            !
              tmp3=tmpa*tmpa*olambdas2D(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
              tmp4=tmp1*rhowater*tmp3
              psacr2D(iq,k)=amin1( tmp4,qrzodt2D(iq,k) )
              tmp1=0.25*pi*esr*n0_r2D(iq,k)*N0_s2D(iq,k)*abs( vtr2D(iq,k)-vts2D(iq,k) )
              tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
              tmp3=tmp1*tmp2
              npsacr2D(iq,k)=amin1( tmp3,nrzodt2D(iq,k) )
            !
            !
            ! (2) FREEZING OF RAIN TO FORM GRAUPEL  (pgfr): Lin (45), added to PI
            !     positive value
            !     Constant in Bigg freezing Aplume=Ap=0.66 /k
            !     Constant in raindrop freezing equ. Bplume=Bp=100./m/m/m/s
            !

              if (qrz2D(iq,k) .gt. 1.e-8 ) then
                Bp=100.
                Ap=0.66
                tmp1=olambdar2D(iq,k)*olambdar2D(iq,k)*olambdar2D(iq,k)
                tmp2=20.*pi*pi*Bp*n0_r2D(iq,k)*rhowater*orho2D(iq,k)*  &
                  (exp(-Ap*temcc2D(iq,k))-1.0)*tmp1*tmp1*olambdar2D(iq,k)
                pgfr2D(iq,k)=amin1( tmp2,qrzodt2D(iq,k) )
                npgfr2D(iq,k)=pi*Bp*n0_r2D(iq,k)*tmpa*tmpa*(exp(-Ap*temcc2D(iq,k))-1.0)
              else
                pgfr2D(iq,k)=0
                npgfr2D(iq,k)=0.
              endif
            end if ! for the go to 1200
          end if   !1200    continue

        else       ! if (tem2D(iq,k) .lt. 273.15) then

        !
        !***********************************************************************
        !*********        snow production processes for T > 0 C       **********
        !***********************************************************************
        !
          if (qsz2D(iq,k) .gt. 0.0)  then !go to 1400
          !
          ! (1) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            save1 =aa_s2D(iq,k)*sqrho2D(iq,k)*N0_s2D(iq,k)* &
                   ggamma(bv_s2D(iq,k)+tmp_sa2D(iq,k))*     &
                   olambdas2D(iq,k)**(bv_s2D(iq,k)+tmp_sa2D(iq,k))

            tmp1=esw*save1
            psacw2D(iq,k)=qlzodt2D(iq,k)*( 1.0-exp(-tmp1*dtb) )
            npsacw2D(iq,k)=tmp1*ncz2D(iq,k)
          !
          ! (2) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
          !
            esr=1.0
            tmpa=olambdar2D(iq,k)*olambdar2D(iq,k)
            tmpb=olambdas2D(iq,k)*olambdas2D(iq,k)
            tmpc=olambdar2D(iq,k)*olambdas2D(iq,k)
            tmp1=pi*pi*esr*n0_r2D(iq,k)*N0_s2D(iq,k)*   &
                    abs( vtr2D(iq,k)-vts2D(iq,k) )*orho2D(iq,k)
          ! tmp1=pi*pi*esr*n0_r(k)*N0_s(k)*            &
          ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5*orho(k)
            tmp2=tmpa*tmpa*olambdas2D(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
            tmp3=tmp1*rhowater*tmp2
            psacr2D(iq,k)=amin1( tmp3,qrzodt2D(iq,k) )

            tmp1=0.25*pi*esr*n0_r2D(iq,k)*N0_s2D(iq,k)*abs( vtr2D(iq,k)-vts2D(iq,k) )
            tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
            tmp3=tmp1*tmp2
            npsacr2D(iq,k)=amin1( tmp3,nrzodt2D(iq,k) )
          !
          ! (3) MELTING OF SNOW (Psmlt): Lin (32)
          !     Psmlt is negative value
          !
            delrs=rs02D(iq,k)-qvz2D(iq,k)
            term1=2.0*pi*orho2D(iq,k)*( xlv*diffwv2D(iq,k)*rho2D(iq,k)*delrs- &
                  xka2D(iq,k)*temcc2D(iq,k) )
            tmp1= av_s2D(iq,k)*sqrho2D(iq,k)*        &
                  olambdas2D(iq,k)**(5+bv_s2D(iq,k)+2*mu_s)/visc2D(iq,k)
            tmp2= N0_s2D(iq,k)*( vf1s*olambdas2D(iq,k)*olambdas2D(iq,k)+ &
                  vf2s*schmidt2D(iq,k)**0.33334* &
                  ggamma(2.5+0.5*bv_s2D(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=term1*oxlf*tmp2-cwoxlf*temcc2D(iq,k)*( psacw2D(iq,k)+psacr2D(iq,k) )
            tmp4=amin1(0.0,tmp3)
            psmlt2D(iq,k)=amax1( tmp4,-qszodt2D(iq,k) )

            if(qsz2D(iq,k) .ge. 0.0) then
              npsmlt2D(iq,k)=psmlt2D(iq,k)*nsz2D(iq,k)/qsz2D(iq,k)
            else
              npsmlt2D(iq,k)=psmlt2D(iq,k)/xms
            end if
          !
          ! (4) EVAPORATION OF MELTING SNOW (Psmltevp): HR (A27)
          !     but use Lin et al. coefficience
          !     Psmltevp is a negative value
          !
            tmpa=rvapor*xka2D(iq,k)*tem2D(iq,k)*tem2D(iq,k)
            tmpb=xlv*xlv*rho2D(iq,k)*qswz2D(iq,k)*diffwv2D(iq,k)
            tmpc=tmpa*qswz2D(iq,k)*diffwv2D(iq,k)
            tmpd=amin1( 0.0,(qvoqswz2D(iq,k)-0.90)*qswz2D(iq,k)*odtb )

            abr=2.0*pi*(qvoqswz2D(iq,k)-0.90)*tmpc/(tmpa+tmpb)
          !
          !**** allow evaporation to occur when RH less than 90%
          !**** here not using 100% because the evaporation cooling
          !**** of temperature is not taking into account yet; hence,
          !**** the qsw value is a little bit larger. This will avoid
          !**** evaporation can generate cloud.
          
            tmp1=av_s2D(iq,k)*sqrho2D(iq,k)*    &
                 olambdas2D(iq,k)**(5+bv_s2D(iq,k)+2*mu_s)/visc2D(iq,k)
            tmp2=N0_s2D(iq,k)*( vf1s*olambdas2D(iq,k)*olambdas2D(iq,k)+ &
                 vf2s*schmidt2D(iq,k)**0.33334* &
                 ggamma(2.5+0.5*bv_s2D(iq,k)+mu_s)*sqrt(tmp1) )
            tmp3=amin1(0.0,tmp2)
            tmp3=amax1( tmp3,tmpd )
            psmltevp2D(iq,k)=amax1( tmp3,-qszodt2D(iq,k) )
            if(qsz2D(iq,k) .ge. 0.0) then
              npsmltevp2D(iq,k)=psmltevp2D(iq,k)*nsz2D(iq,k)/qsz2D(iq,k)
            else
              npsmltevp2D(iq,k)=psmltevp2D(iq,k)/xmr
            end if
          end if !          1400     continue
        end if      !---- end of snow/ice processes   if (tem2D(iq,k) .lt. 273.15) then
        !---------- end of snow/ice processes below freezing

        !***********************************************************************
        !*********           rain production processes                **********
        !***********************************************************************

        !
        ! (1) AUTOCONVERSION OF RAIN (Praut): using Liu and Daum (2004)
        !

        !---- YLIN, autoconversion use Liu and Daum (2004), unit = g cm-3 s-1, in the scheme kg/kg s-1, so

        !if (AEROSOL .eq. .true.) then
        !  print*, ' AEROSOL option here'

        !else
        if (qlz2D(iq,k) .gt. 1e-6) then
          mu_c    = AMIN1(15., (1000.E6/ncz2D(iq,k) + 2.))
          lamc2D(iq,k) = (ncz2D(iq,k)*rhowater*pi*ggamma(4.+mu_c)/(6.*qlz2D(iq,k)*ggamma(1+mu_c)))**(1./3)
          Dc_liu  = (ggamma(6+1+mu_c)/ggamma(1+mu_c))**(1./6.)/lamc2D(iq,k)             !----- R6 in m
          if (Dc_liu .gt. R6c) then
            disp = 1./(mu_c+1.)      !--- square of relative dispersion
            eta  = (0.75/pi/(1e-3*rhowater))**2*1.9e11*((1+3*disp)*(1+4*disp)*&
                   (1+5*disp)/(1+disp)/(1+2*disp))
            praut2D(iq,k) = eta*(1e-3*rho2D(iq,k)*qlz2D(iq,k))**3/(1e-6*ncz2D(iq,k))  !--- g cm-3 s-1
            praut2D(iq,k) = praut2D(iq,k)/(1e-3*rho2D(iq,k))                          !--- kg kg-1 s-1
            npraut_r2D(iq,k) = praut2D(iq,k)/xmr                                                !--- kg kg-1 s-1
            npraut2D(iq,k) = praut2D(iq,k)/qlz2D(iq,k)*ncz2D(iq,k)                                        !--- kg kg-1 s-1
            npraut2D(iq,k) = praut2D(iq,k)/xmr                                                  !--- kg kg-1 s-1
          else
            praut2D(iq,k) = 0.0
            npraut2D(iq,k) = 0.0
            npraut_r2D(iq,k) = 0.0
          endif
        else
          praut2D(iq,k) = 0.0
          npraut2D(iq,k) = 0.0
          npraut_r2D(iq,k) = 0.0
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
        tmp1=pio4*erw*n0_r2D(iq,k)*av_r*sqrho2D(iq,k)* &
             gambp3*olambdar2D(iq,k)**bp3 ! son
        pracw2D(iq,k)=qlzodt2D(iq,k)*( 1.0-exp(-tmp1*dtb) )
        npracw2D(iq,k)=tmp1*ncz2D(iq,k)
        
        !
        ! (3) EVAPORATION OF RAIN (Prevp): Lin (52)
        !     Prevp is negative value
        !
        !     Sw=qvoqsw : saturation ratio
        !
        tmpa=rvapor*xka2D(iq,k)*tem2D(iq,k)*tem2D(iq,k)
        tmpb=xlv*xlv*rho2D(iq,k)*qswz2D(iq,k)*diffwv2D(iq,k)
        tmpc=tmpa*qswz2D(iq,k)*diffwv2D(iq,k)
        tmpd=amin1(0.0,(qvoqswz2D(iq,k)-0.99)*qswz2D(iq,k)*odtb)

        abr=2.0*pi*(qvoqswz2D(iq,k)-0.99)*tmpc/(tmpa+tmpb)
        tmp1=av_r*sqrho2D(iq,k)*olambdar2D(iq,k)**bp5/visc2D(iq,k) !son
        tmp2=abr*n0_r2D(iq,k)*( vf1r*olambdar2D(iq,k)*olambdar2D(iq,k)+  &
             vf2r*schmidt2D(iq,k)**0.33334*gambp5o2*sqrt(tmp1) )
        tmp3=amin1( 0.0,tmp2 )
        tmp3=amax1( tmp3,tmpd )
        prevp2D(iq,k)=amax1( tmp3,-qrzodt2D(iq,k) )
        if (qrz2D(iq,k).gt.0.) then
          nprevp2D(iq,k)=prevp2D(iq,k)*nrz2D(iq,k)/qrz2D(iq,k)
        else
          nprevp2D(iq,k)=prevp2D(iq,k)*xmr
        end if

        ! CALL wrf_debug ( 100 , 'module_ylin: finish rain processes' )
        !
        !**********************************************************************
        !*****     combine all processes together and avoid negative      *****
        !*****     water substances
        !***********************************************************************
        !
        if ( temcc2D(iq,k) .lt. 0.0) then
        !
        !  combined water vapor depletions
        !
          tmp=psdep2D(iq,k) + midep2D(iq,k)
          if ( tmp .gt. qvzodt2D(iq,k) ) then
            factor=qvzodt2D(iq,k)/tmp
            psdep2D(iq,k)=psdep2D(iq,k)*factor
            midep2D(iq,k)=midep2D(iq,k)*factor
          end if
        !
        !  combined cloud water depletions
        !
          tmp=praut2D(iq,k)+psacw2D(iq,k)+psfw2D(iq,k)+pracw2D(iq,k)
          if ( tmp .gt. qlzodt2D(iq,k) ) then
            factor=qlzodt2D(iq,k)/tmp
            praut2D(iq,k)=praut2D(iq,k)*factor
            psacw2D(iq,k)=psacw2D(iq,k)*factor
            psfw2D(iq,k)=psfw2D(iq,k)*factor
            pracw2D(iq,k)=pracw2D(iq,k)*factor
          end if
        !
        !  combined cloud ice depletions
        !
          tmp=psaut2D(iq,k)+psaci2D(iq,k)+praci2D(iq,k)+psfi2D(iq,k)
          if (tmp .gt. qizodt2D(iq,k) ) then
            factor=qizodt2D(iq,k)/tmp
            psaut2D(iq,k)=psaut2D(iq,k)*factor
            psaci2D(iq,k)=psaci2D(iq,k)*factor
            praci2D(iq,k)=praci2D(iq,k)*factor
            psfi2D(iq,k)=psfi2D(iq,k)*factor
          endif

        !
        !  combined all rain processes
        !
          tmp_r=piacr2D(iq,k)+psacr2D(iq,k)-prevp2D(iq,k)-  & 
                praut2D(iq,k)-pracw2D(iq,k)+pgfr2D(iq,k)
          if (tmp_r .gt. qrzodt2D(iq,k) ) then
            factor=qrzodt2D(iq,k)/tmp_r
            piacr2D(iq,k)=piacr2D(iq,k)*factor
            psacr2D(iq,k)=psacr2D(iq,k)*factor
            prevp2D(iq,k)=prevp2D(iq,k)*factor
            pgfr2D(iq,k)=pgfr2D(iq,k)*factor
          endif
        !
        !   combined all snow processes
        !
          tmp_s=-pssub2D(iq,k)-(psaut2D(iq,k)+psaci2D(iq,k)+ &
                 psacw2D(iq,k)+psfw2D(iq,k)+pgfr2D(iq,k)+ &
                 psfi2D(iq,k)+praci2D(iq,k)+piacr2D(iq,k)+ &
                 psdep2D(iq,k)+psacr2D(iq,k)-pracs2D(iq,k))
          if ( tmp_s .gt. qszodt2D(iq,k) ) then
            factor=qszodt2D(iq,k)/tmp_s
            pssub2D(iq,k)=pssub2D(iq,k)*factor
            pracs2D(iq,k)=Pracs2D(iq,k)*factor
          endif

        !
        ! calculate new water substances, thetae, tem, and qvsbar
        !

          pvapor2D(iq,k)=-pssub2D(iq,k)-psdep2D(iq,k)-prevp2D(iq,k)-midep2D(iq,k)
          qvz2D(iq,k)=amax1( qvmin,qvz2D(iq,k)+dtb*pvapor2D(iq,k) )
          pclw2D(iq,k)=-praut2D(iq,k)-pracw2D(iq,k)-psacw2D(iq,k)-psfw2D(iq,k)
          qlz2D(iq,k)=amax1( 0.0,qlz2D(iq,k)+dtb*pclw2D(iq,k) )
          pcli2D(iq,k)=-psaut2D(iq,k)-psfi2D(iq,k)-psaci2D(iq,k)- & 
                    praci2D(iq,k)+midep2D(iq,k)
          qiz2D(iq,k)=amax1( 0.0,qiz2D(iq,k)+dtb*pcli2D(iq,k) )
          tmp_r=piacr2D(iq,k)+psacr2D(iq,k)-prevp2D(iq,k)-praut2D(iq,k)- &
                    pracw2D(iq,k)+pgfr2D(iq,k)-pracs2D(iq,k)
          prain2D(iq,k)=-tmp_r
          qrz2D(iq,k)=amax1( 0.0,qrz2D(iq,k)+dtb*prain2D(iq,k) )
          tmp_s=-pssub2D(iq,k)-(psaut2D(iq,k)+psaci2D(iq,k)+ &
                 psacw2D(iq,k)+psfw2D(iq,k)+pgfr2D(iq,k)+  &
                 psfi2D(iq,k)+praci2D(iq,k)+piacr2D(iq,k)+  &
                 psdep2D(iq,k)+psacr2D(iq,k)-pracs2D(iq,k))
          psnow2D(iq,k)=-tmp_s
          qsz2D(iq,k)=amax1( 0.0,qsz2D(iq,k)+dtb*psnow2D(iq,k) )
          qschg2D(iq,k)=qschg2D(iq,k)+psnow2D(iq,k)
          qschg2D(iq,k)=psnow2D(iq,k)

          tmp=ocp/tothz2D(iq,k)*xLf*qschg2D(iq,k)
          theiz2D(iq,k)=theiz2D(iq,k)+dtb*tmp
          ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
          ! tem(k)=thz(k)*tothz(k)
          ! temcc(k)=tem(k)-273.15
          !==================update temperature=================================================
          temcc2D(iq,k)=tem2D(iq,k)-273.15
          lvap = xlv + (2106.0 - 4218.0)*temcc2D(iq,k)  !Enthalpy of vaporization
          tmp1=(pssub2D(iq,k)+psdep2D(iq,k))*xls*ocp + prevp2D(iq,k)*lvap*ocp+  &
               (psfw2D(iq,k)+pgfr2D(iq,k)+psacr2D(iq,k)-pracs2D(iq,k))*xlfocp
          !bug fixed 20191126
          tem2D(iq,k)=tem2D(iq,k)+tmp1*dtb
          temcc2D(iq,k)=tem2D(iq,k)-273.15
          thz2D(iq,k)=tem2D(iq,k)/tothz2D(iq,k)
          !===================================================================
          if( temcc2D(iq,k) .lt. -40.0 ) qswz2D(iq,k)=qsiz2D(iq,k)
            qlpqi=qlz2D(iq,k)+qiz2D(iq,k)
            if ( qlpqi .eq. 0.0 ) then
              qvsbar2D(iq,k)=qsiz2D(iq,k)
            else
              qvsbar2D(iq,k)=(qiz2D(iq,k)*qsiz2D(iq,k)+qlz2D(iq,k)*qswz2D(iq,k))/qlpqi
            endif
            tmp1=-npraut2D(iq,k)-npracw2D(iq,k)-npsacw2D(iq,k)
            ncz2D(iq,k)=amax1( 0.0,ncz2D(iq,k)+dtb*tmp1 )
            tmp1=-npsaut2D(iq,k)-npsaci2D(iq,k)-npraci2D(iq,k)+nidep2D(iq,k)
            niz2D(iq,k)=amax1( 0.0,niz2D(iq,k)+dtb*tmp1 )
            tmp1=npiacr2D(iq,k)+npsacr2D(iq,k)-nprevp2D(iq,k)-npraut_r2D(iq,k)+npgfr2D(iq,k)
            nrz2D(iq,k)=amax1( 0.0,nrz2D(iq,k)-dtb*tmp1 )
            tmp1=-(npsaut2D(iq,k)+npgfr2D(iq,k)+  &
                   npraci2D(iq,k)+npiacr2D(iq,k)+  &
                   npsdep2D(iq,k)+npsacr2D(iq,k))
            nsz2D(iq,k)=amax1( 0.0,nsz2D(iq,k)-dtb*tmp1 )
          else                  !>0 C
          !
          !  combined cloud water depletions
          !
            tmp=praut2D(iq,k)+psacw2D(iq,k)+pracw2D(iq,k)
            if ( tmp .gt. qlzodt2D(iq,k) ) then
              factor=qlzodt2D(iq,k)/tmp
              praut2D(iq,k)=praut2D(iq,k)*factor
              psacw2D(iq,k)=psacw2D(iq,k)*factor
              pracw2D(iq,k)=pracw2D(iq,k)*factor
            end if
          !
          !  combined all snow processes
          !
            tmp_s=-(psmlt2D(iq,k)+psmltevp2D(iq,k))
            if (tmp_s .gt. qszodt2D(iq,k) ) then
              factor=qszodt2D(iq,k)/tmp_s
              psmlt2D(iq,k)=psmlt2D(iq,k)*factor
              psmltevp2D(iq,k)=psmltevp2D(iq,k)*factor
            endif
          !
          !  combined all rain processes
          !
            tmp_r=-prevp2D(iq,k)-(praut2D(iq,k)+pracw2D(iq,k)+psacw2D(iq,k)-psmlt2D(iq,k))
            if (tmp_r .gt. qrzodt2D(iq,k) ) then
              factor=qrzodt2D(iq,k)/tmp_r
              prevp2D(iq,k)=prevp2D(iq,k)*factor
            endif
          !
          !  calculate new water substances and thetae
          !
            pvapor2D(iq,k)=-psmltevp2D(iq,k)-prevp2D(iq,k)
            qvz2D(iq,k)=amax1( qvmin,qvz2D(iq,k)+dtb*pvapor2D(iq,k))
            pclw2D(iq,k)=-praut2D(iq,k)-pracw2D(iq,k)-psacw2D(iq,k)
            qlz2D(iq,k)=amax1( 0.0,qlz2D(iq,k)+dtb*pclw2D(iq,k) )
            pcli2D(iq,k)=0.0
            qiz2D(iq,k)=amax1( 0.0,qiz2D(iq,k)+dtb*pcli2D(iq,k) )
            tmp_r=-prevp2D(iq,k)-(praut2D(iq,k)+pracw2D(iq,k)+psacw2D(iq,k)-psmlt2D(iq,k))
            prain2D(iq,k)=-tmp_r
            tmpqrz=qrz2D(iq,k)
            qrz2D(iq,k)=amax1( 0.0,qrz2D(iq,k)+dtb*prain2D(iq,k) )
            tmp_s=-(psmlt2D(iq,k)+psmltevp2D(iq,k))
            psnow2D(iq,k)=-tmp_s
            qsz2D(iq,k)=amax1( 0.0,qsz2D(iq,k)+dtb*psnow2D(iq,k) )
            qschg2D(iq,k)=psnow2D(iq,k)

            tmp=ocp/tothz2D(iq,k)*xLf*qschg2D(iq,k)
            theiz2D(iq,k)=theiz2D(iq,k)+dtb*tmp
            ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
            ! tem(k)=thz(k)*tothz(k)
            ! temcc(k)=tem(k)-273.15
            !==================update tmperature=================================================
            temcc2D(iq,k)=tem2D(iq,k)-273.15
            lvap = xlv + (2106.0 - 4218.0)*temcc2D(iq,k)  !Enthalpy of vaporization
            tmp1=psmltevp2D(iq,k)*xls*ocp + prevp2D(iq,k)*lvap*ocp+  &
                 psmlt2D(iq,k)*xlfocp
            !tmp1 =  ! 1. evaporation of rain formed by melting snow ??? (-)
                     ! 2. evaporation of rain (-)
                     ! 3. melting of snow to form rain (+)
            tem2D(iq,k)=tem2D(iq,k)+tmp1*dtb
            !bugfix 20191126

            !tem(k)=tem(k)+tmp1*dtb
            temcc2D(iq,k)=tem2D(iq,k)-273.15

            thz2D(iq,k)=tem2D(iq,k)/tothz2D(iq,k)

            !===================================================================
            es=1000.*svp1*exp( svp2*temcc2D(iq,k)/(tem2D(iq,k)-svp3) )
            qswz2D(iq,k)=ep2*es/(prez2D(iq,k)-es)
            qsiz2D(iq,k)=qswz2D(iq,k)
            qvsbar2D(iq,k)=qswz2D(iq,k)
            ! tmp1=-(npraut2D(iq,k)+npsacw2D(iq,k)+npracw2D(iq,k))
            ncz2D(iq,k)=amax1( 0.0,ncz2D(iq,k)+dtb*tmp1)
            tmp1=-nprevp2D(iq,k)-(npraut_r2D(iq,k)-npsmlt2D(iq,k))
            ! tmp1=-nprevp(k)-(nprautr(k)+npracwr(k)+npsacw(k)-npsmltr(k))
            nrz2D(iq,k)=amax1(0.0,nrz2D(iq,k)-dtb*tmp1)
            tmp1=-(npsmlt2D(iq,k)+npsmltevp2D(iq,k))
            nsz2D(iq,k)=amax1( 0.0,nsz2D(iq,k)-dtb*tmp1 )
          end if    !T seperate for source and sink terms
          ! CALL wrf_debug ( 100 , 'module_ylin: finish sum of all processes' )

        !rain
        if (qrz2D(iq,k) .gt. 1.0e-8) then
          xlambdar2D(iq,k)=(pi*rhowater*nrz2D(iq,k)/qrz2D(iq,k))**(1./3.)   !zx
          if (xlambdar2D(iq,k).lt.lamminr) then
            xlambdar2D(iq,k) = lamminr
            n0_r2D(iq,K) = xlambdar2D(iq,K)**4*qrz2D(iq,K)/(pi*rhowater)
            nrz2D(iq,K) = n0_r2D(iq,K)/xlambdar2D(iq,K)
          else if (xlambdar2D(iq,K).gt.lammaxr) then
            xlambdar2D(iq,K) = lammaxr
            n0_r2D(iq,K) = xlambdar2D(iq,K)**4*qrz2D(iq,K)/(pi*rhowater)
            nrz2D(iq,K) = n0_r2D(iq,K)/xlambdar2D(iq,K)
          end if
        end if

        !snow
        if (qsz2D(iq,k) .gt. 1.0e-8) then
          xlambdas2D(iq,k)=(am_s2D(iq,k)*ggamma(tmp_ss2D(iq,k))*     &
          nsz2D(iq,k)/qsz2D(iq,k))**(1./bm_s2D(iq,k))
          if (xlambdas2D(iq,k).lt.lammins) then
            xlambdas2D(iq,k)= lamminS
            n0_s2D(iq,K) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*      &
                           qsz2D(iq,K)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
            nsz2D(iq,K) = n0_s2D(iq,K)/xlambdas2D(iq,k)
          else if (xlambdas2D(iq,k).gt.lammaxs) then
            xlambdas2D(iq,k) = lammaxs
            n0_s2D(iq,K) = xlambdas2D(iq,k)**(bm_s2D(iq,k)+1)*      & 
                           qsz2D(iq,K)/ggamma(1+bm_s2D(iq,k))/am_s2D(iq,k)
            nsz2D(iq,K) = n0_s2D(iq,K)/xlambdas2D(iq,k)
          end if
        end if

        !cloud ice
        if (qiz2D(iq,k).ge.1.0e-8) then
          lami2D(iq,k) = max((ggamma(1.+3.)*500.*pi/6.)*niz2D(iq,k)/qiz2D(iq,k),1.e-20)**(1./3) !fixed zdc
          if (lami2D(iq,k).lt.lammini) then
            lami2D(iq,k)= lammini
            n0_i2D(iq,K) = lami2D(iq,k)**4./ggamma(1.+3.)*500.*pi/6.
            niz2D(iq,K) = n0_i2D(iq,K)/lami2D(iq,k)
          else if (lami2D(iq,k).gt.lammaxi) then
            lami2D(iq,k) = lammaxi
            n0_i2D(iq,K) = lami2D(iq,k)**4./ggamma(1.+3.)*500.*pi/6.
            niz2D(iq,K) = n0_i2D(iq,K)/lami2D(iq,k)
          end if
        end if

        !cloud water zdc 20220208
        if (qlz2D(iq,k).ge.1.0e-8) then
          lamc2D(iq,k) = (ncz2D(iq,k)*rhowater*pi*ggamma(4.+mu_c)/(6.*qlz2D(iq,k)*ggamma(1+mu_c)))**(1./3)
          if (lamc2D(iq,k).lt.lammini) then
            lamc2D(iq,k)= lammini
            n0_c2D(iq,k)= lamc2D(iq,k)**(mu_c+4.)*6.*qlz2D(iq,k)/(pi*rhowater*ggamma(mu_c+4))
            ncz2D(iq,k) = n0_c2D(iq,k)/lamc2D(iq,k)
          else if (lamc2D(iq,k).gt.lammaxi) then
            lamc2D(iq,k)= lammaxi
            n0_c2D(iq,k)= lamc2D(iq,k)**(mu_c+4.)*6.*qlz2D(iq,k)/(pi*rhowater*ggamma(mu_c+4))
            ncz2D(iq,k) = n0_c2D(iq,k)/lamc2D(iq,k)
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
        if( qvz2D(iq,k)+qlz2D(iq,k)+qiz2D(iq,k) .lt. rsat*qvsbar2D(iq,k) ) then ! goto 1800

        !
        !   unsaturated
        !
          qvz2D(iq,k)=qvz2D(iq,k)+qlz2D(iq,k)+qiz2D(iq,k)
          qlz2D(iq,k)=0.0
          qiz2D(iq,k)=0.0

          thz2D(iq,k)=theiz2D(iq,k)-(xLvocp*qvz2D(iq,k)-xLfocp*qiz2D(iq,k))/tothz2D(iq,k)

          tem2D(iq,k)=thz2D(iq,k)*tothz2D(iq,k)
          temcc2D(iq,k)=tem2D(iq,k)-273.15

        else
        !
        !   saturated
        !
          pladj2D(iq,k)=qlz2D(iq,k)
          piadj2D(iq,k)=qiz2D(iq,k)
        !

          CALL satadj(qvz2D(iq,:), qlz2D(iq,:), qiz2D(iq,:), prez2D(iq,:), &
                      theiz2D(iq,:), thz2D(iq,:), tothz2D(iq,:), kts, kte, &
                      k, xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

          pladj2D(iq,k)=odtb*(qlz2D(iq,k)-pladj2D(iq,k))
          piadj2D(iq,k)=odtb*(qiz2D(iq,k)-piadj2D(iq,k))
          pclw2D(iq,k)=pclw2D(iq,k)+pladj2D(iq,k)
          pcli2D(iq,k)=pcli2D(iq,k)+piadj2D(iq,k)
          pvapor2D(iq,k)=pvapor2D(iq,k)-( pladj2D(iq,k)+piadj2D(iq,k) )
          thz2D(iq,k)=theiz2D(iq,k)-(xLvocp*qvz2D(iq,k)-xLfocp*qiz2D(iq,k))/tothz2D(iq,k)
          tem2D(iq,k)=thz2D(iq,k)*tothz2D(iq,k)
          temcc2D(iq,k)=tem2D(iq,k)-273.15

          es=1000.*svp1*exp( svp2*temcc2D(iq,k)/(tem2D(iq,k)-svp3) )
          qswz2D(iq,k)=ep2*es/(prez2D(iq,k)-es)
          if (tem2D(iq,k) .lt. 273.15 ) then
            es=1000.*svp1*exp( 21.8745584*(tem2D(iq,k)-273.16)/(tem2D(iq,k)-7.66) )
            qsiz2D(iq,k)=ep2*es/(prez2D(iq,k)-es)
            if (temcc2D(iq,k) .lt. -40.0) qswz2D(iq,k)=qsiz2D(iq,k)
          else
            qsiz2D(iq,k)=qswz2D(iq,k)
          endif
          qlpqi=qlz2D(iq,k)+qiz2D(iq,k)
          if ( qlpqi .eq. 0.0 ) then
            qvsbar2D(iq,k)=qsiz2D(iq,k)
          else
            qvsbar2D(iq,k)=( qiz2D(iq,k)*qsiz2D(iq,k)+qlz2D(iq,k)*qswz2D(iq,k) )/qlpqi
          endif

        !
        !***********************************************************************
        !*****     melting and freezing of cloud ice and cloud water       *****
        !***********************************************************************
          qlpqi=qlz2D(iq,k)+qiz2D(iq,k)
          if( qlpqi .gt. 0.0 ) then !go to 1800
          !
          !
          ! (1)  HOMOGENEOUS NUCLEATION WHEN T< -40 C (Pihom)
          !
            if(temcc2D(iq,k) .lt. -40.0) then
              pihom2D(iq,k)=qlz2D(iq,k)*odtb
              nihom2D(iq,k)=ncz2D(iq,k)*odtb
            end if
          !
          ! (2)  MELTING OF ICE CRYSTAL WHEN T> 0 C (Pimlt)
          !
            if(temcc2D(iq,k) .gt. 0.0) then
              pimlt2D(iq,k)=qiz2D(iq,k)*odtb
              nimlt2D(iq,k)=niz2D(iq,k)*odtb
            end if
          !
          ! (3) PRODUCTION OF CLOUD ICE BY BERGERON PROCESS (Pidw): Hsie (p957)
          !     this process only considered when -31 C < T < 0 C
          !
            if(temcc2D(iq,k) .lt. 0.0 .and. temcc2D(iq,k) .gt. -31.0) then
          !!
          !!   parama1 and parama2 functions must be user supplied
          !!
              a1=parama1( temcc2D(iq,k) )
              a2=parama2( temcc2D(iq,k) )
              !! change unit from cgs to mks
              a1=a1*0.001**(1.0-a2)
              xnin=xni0*exp(-bni*temcc2D(iq,k))
              pidw2D(iq,k)=xnin*orho2D(iq,k)*(a1*xmnin**a2)
            end if

            pcli2D(iq,k)=pcli2D(iq,k)+pihom2D(iq,k)-pimlt2D(iq,k)+pidw2D(iq,k)
            pclw2D(iq,k)=pclw2D(iq,k)-pihom2D(iq,k)+pimlt2D(iq,k)-pidw2D(iq,k)
            qlz2D(iq,k)=amax1( 0.0,qlz2D(iq,k)+dtb*(-pihom2D(iq,k)+pimlt2D(iq,k)-pidw2D(iq,k)) )
            qiz2D(iq,k)=amax1( 0.0,qiz2D(iq,k)+dtb*(pihom2D(iq,k)-pimlt2D(iq,k)+pidw2D(iq,k)) )

            ncz2D(iq,k)=amax1( 0.0,ncz2D(iq,k)+dtb*(-nihom2D(iq,k)+nimlt2D(iq,k)) )
            niz2D(iq,k)=amax1( 0.0,niz2D(iq,k)+dtb*( nihom2D(iq,k)-nimlt2D(iq,k)) )

            CALL satadj(qvz2D(iq,:), qlz2D(iq,:), qiz2D(iq,:), prez2D(iq,:), &
                    theiz2D(iq,:), thz2D(iq,:), tothz2D(iq,:), kts, kte, &
                    k, xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

            thz2D(iq,k)=theiz2D(iq,k)-(xLvocp*qvz2D(iq,k)-xLfocp*qiz2D(iq,k))/tothz2D(iq,k)
            tem2D(iq,k)=thz2D(iq,k)*tothz2D(iq,k)
            temcc2D(iq,k)=tem2D(iq,k)-273.15
            es=1000.*svp1*exp( svp2*temcc2D(iq,k)/(tem2D(iq,k)-svp3) )
            qswz2D(iq,k)=ep2*es/(prez2D(iq,k)-es)

            if (tem2D(iq,k) .lt. 273.15 ) then
              es=1000.*svp1*exp( 21.8745584*(tem2D(iq,k)-273.16)/(tem2D(iq,k)-7.66) )
              qsiz2D(iq,k)=ep2*es/(prez2D(iq,k)-es)
              if (temcc2D(iq,k) .lt. -40.0) qswz2D(iq,k)=qsiz2D(iq,k)
            else
              qsiz2D(iq,k)=qswz2D(iq,k)
            endif
            qlpqi=qlz2D(iq,k)+qiz2D(iq,k)

            if ( qlpqi .eq. 0.0 ) then
              qvsbar2D(iq,k)=qsiz2D(iq,k)
            else
              qvsbar2D(iq,k)=( qiz2D(iq,k)*qsiz2D(iq,k)+qlz2D(iq,k)*qswz2D(iq,k) )/qlpqi
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


  do k=kts+1,kte
    where ( qvz2D(1:imax,k) .lt. qvmin ) 
      qlz2D(1:imax,k)=0.0
      qiz2D(1:imax,k)=0.0
      ncz2D(1:imax,k)=0.0
      niz2D(1:imax,k)=0.0
      qvz2D(1:imax,k)=amax1( qvmin,qvz2D(1:imax,k)+qlz2D(1:imax,k)+qiz2D(1:imax,k) )
    end where
      niz2D(1:imax,k) = min(niz2D(1:imax,k),0.3E6/rho2D(1:imax,k))
      ncz2D(1:imax,k) = min(ncz2D(1:imax,k),250000.E6/rho2D(1:imax,k))
      ncz2D(1:imax,k) = max(ncz2D(1:imax,k),0.01E6/rho2D(1:imax,k))
  end do


  ! CALCULATE EFFECTIVE RADIUS zdc 20220208
  do k=kts,kte
    where (qiz2D(1:imax,k) .gt. 1.0e-8 .and. lami2D(1:imax,k)>0. ) 
      EFFI1D2D(1:imax,k) = 3./LAMI2D(1:imax,k)/2.
    elsewhere
      EFFI1D2D(1:imax,k) = 25.E-6
    end where

    where (qsz2D(1:imax,k) .gt. 1.0e-8) 
      EFFS1D2D(1:imax,k) = 3./xlambdas2D(1:imax,k)/2.
    elsewhere
      EFFS1D2D(1:imax,k) = 25.E-6
    end where

    where (qrz2D(1:imax,k) .gt. 1.0e-8) 
      EFFR1D2D(1:imax,k) = 3./xlambdar2D(1:imax,k)/2.
    elsewhere
      EFFR1D2D(1:imax,k) = 25.E-6
    end where

    where (qlz2D(1:imax,k) .gt. 1.0e-8 .and. lamc2D(1:imax,k) >0.)
      EFFC1D2D(1:imax,k) = GAMMA(mu_c+4.)/GAMMA(mu_c+3.)/LAMC2D(1:imax,k)/2.
    elsewhere
      EFFC1D2D(1:imax,k) = 25.E-6
    end where
  end do

  ! save all process rate for understanding cloud microphysics
  do k=kts,kte
    zpsnow2D(1:imax,k)   = psnow2D(1:imax,k)    ! sum all process for snow
    zpsaut2D(1:imax,k)   = psaut2D(1:imax,k)    ! ice crystal aggregation to snow
    zpsfw2D(1:imax,k)    = psfw2D(1:imax,k)     ! BERGERON process to transfer cloud water to snow
    zpsfi2D(1:imax,k)    = psfi2D(1:imax,k)     ! BERGERON process to transfer cloud ice to snow
    zpraci2D(1:imax,k)   = praci2D(1:imax,k)    ! cloud ice accretion by rain
    zpiacr2D(1:imax,k)   = piacr2D(1:imax,k)    ! rain accretion by cloud ice
    zpsaci2D(1:imax,k)   = psaci2D(1:imax,k)    ! ice crystal accretion by snow
    zpsacw2D(1:imax,k)   = psacw2D(1:imax,k)    ! accretion of cloud water by snow
    zpsdep2D(1:imax,k)   = psdep2D(1:imax,k)    ! deposition of snow
    zpssub2D(1:imax,k)   = pssub2D(1:imax,k)    ! sublimation of snow (T<0)
    zpracs2D(1:imax,k)   = pracs2D(1:imax,k)    ! accretion of snow by rain
    zpsacr2D(1:imax,k)   = psacr2D(1:imax,k)    ! accretion of rain by snow
    zpsmlt2D(1:imax,k)   = psmlt2D(1:imax,k)    ! melting of snow
    zpsmltevp2D(1:imax,k)= psmltevp2D(1:imax,k) ! evaporation of melting snow (T>0)
    zprain2D(1:imax,k)   = prain2D(1:imax,k)    ! sum all process for rain
    zpraut2D(1:imax,k)   = praut2D(1:imax,k)    ! autoconversion of rain
    zpracw2D(1:imax,k)   = pracw2D(1:imax,k)    ! accretion of cloud water by rain
    zprevp2D(1:imax,k)   = prevp2D(1:imax,k)    ! evaporation of rain
    zpgfr2D(1:imax,k)    = pgfr2D(1:imax,k)     ! feezing of rain to form graupel (added to PI)
    zpvapor2D(1:imax,k)  = pvapor2D(1:imax,k)   ! sum all process for water vapor to determine qvz
    zpclw2D(1:imax,k)    = pclw2D(1:imax,k)     ! sum all process for cloud liquid to determine qlz
    zpladj2D(1:imax,k)   = pladj2D(1:imax,k)    ! saturation adjustment for ql
    zpcli2D(1:imax,k)    = pcli2D(1:imax,k)     ! sum all process for cloud ice to determine qiz
    zpimlt2D(1:imax,k)   = pimlt2D(1:imax,k)    ! melting of ice crystal >0.
    zpihom2D(1:imax,k)   = pihom2D(1:imax,k)    ! homogeneous nucleation <-40
    zpidw2D(1:imax,k)    = pidw2D(1:imax,k)     ! production of cloud ice by BERGERON process
    zpiadj2D(1:imax,k)   = piadj2D(1:imax,k)    ! saturation adjustment for qi
    zqschg2D(1:imax,k)   = qschg2D(1:imax,k)    ! = psnow / unsure
  enddo


  ! save process rate for aerisol scheme
  do k=kts,kte
    fluxi2D(1:imax,k) = fluxi2D(1:imax,k)                         ! - ice flux leaving layer k to k-1 (kg/m2/s)
    fluxs2D(1:imax,k) = fluxs2D(1:imax,k)                         ! - snow flux leaving layer k to k-1 (kg/m2/s)
    fluxr2D(1:imax,k) = fluxr2D(1:imax,k)                         ! - rain flux leaving layer k to k-1 (kg/m2/s)
    fluxg2D(1:imax,k) = 0.                                        ! - graupel flux leving layer k to k-1 (kg/m2/s)
    fluxm2D(1:imax,k) = -1.*(psmlt2D(1:imax,k)*dzw2D(1:imax,k)*  &!
                        rho2D(1:imax,k)+psmltevp2D(1:imax,k)*    &!
                        dzw2D(1:imax,k)*rho2D(1:imax,k))          ! - ice melting flux in layer k (kg/m2/s)
    fluxf2D(1:imax,k) = pgfr2D(1:imax,k)*dzw2D(1:imax,k)*        &!
                        rho2D(1:imax,k)                           ! - liquid freezing flux in layer k (kg/m2/s)
    fevap2D(1:imax,k) = -1.*prevp2D(1:imax,k)*dzw2D(1:imax,k)*   &!
                        rho2D(1:imax,k)                           ! - evaporation of rainfall flux (kg/m2/s)
    fsubl2D(1:imax,k) = -1.*pssub2D(1:imax,k)*dzw2D(1:imax,k)*   &!
                        rho2D(1:imax,k)                           ! - sublimation of snow, ice and graupel flux (kg/m2/s)
    fauto2D(1:imax,k) = praut2D(1:imax,k)*dzw2D(1:imax,k)*       &!
                        rho2D(1:imax,k)                           ! - autoconversion flux for rainfall (kg/m2/s)
    fcoll2D(1:imax,k) = pracw2D(1:imax,k)*dzw2D(1:imax,k)*       &!
                        rho2D(1:imax,k)                           ! - collection of cloud liquid water by rain (kg/m2/s)
    faccr2D(1:imax,k) = psacw2D(1:imax,k)*dzw2D(1:imax,k)*       &!
                        rho2D(1:imax,k) + psfw2D(1:imax,k)*      &! - accretion of cloud liq water by snow,ice and graupel (kg/m2/s)
                        dzw2D(1:imax,k)*rho2D(1:imax,k)
    vi2D(1:imax,k)    = vtiold2D(1:imax,k)
    vs2D(1:imax,k)    = vtsold2D(1:imax,k)
    vg2D(1:imax,k)    = 0.
  end do

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
