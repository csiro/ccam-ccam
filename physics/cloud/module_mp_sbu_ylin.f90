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

! Zhao, X., Lin, Y., Luo, Y., Qian, Q., Liu, X., Liu, X., & Colle, B. A. (2021).
! A double-moment SBU-YLIN cloud microphysics scheme and its impact
! on a squall line simulation. Journal of Advances in Modeling Earth Systems,
! 13, e2021MS002545. https://doi.org/10.1029/2021MS002545

! Modified code to improve vectorization with Sonny Truong and Marcus Thatcher
! in CCAM
    
MODULE module_mp_sbu_ylin

  implicit none

  private
  public clphy1d_ylin

!..Parameters user might change based on their need
  real, parameter :: RH = 1.0
  real, parameter :: xnor = 8.0e6
  real, parameter :: Nt_c = 250.E6 
!..Water vapor and air gas constants at constant pressure
  real, parameter :: Rvapor = 461.5
  real, parameter :: oRv    = 1./Rvapor
  real, parameter :: Rair   = 287.04
  real, parameter :: Cp     = 1004.0
  real, parameter :: grav   = 9.81
  real, parameter :: rhowater = 1000.0
  real, parameter :: rhosnow  = 100.0
    
  real, parameter :: SVP1 = 0.6112
  real, parameter :: SVP2 = 17.67
  real, parameter :: SVP3 = 29.65
  real, parameter :: SVPT0 = 273.15
  real, parameter :: EP1 = Rvapor/Rair-1.
  real, parameter :: EP2 = Rair/Rvapor
!..Enthalpy of sublimation, vaporization, and fusion at 0C.
  real, parameter :: XLS = 2.834E6
  real, parameter :: XLV = 2.5E6
  real, parameter :: XLF = XLS - XLV
    
  !REAL, SAVE, PUBLIC :: qi0 = 1.0e-3   
  real, parameter ::                                        &
             qi0 = 1.0e-3,                                  &   !--- ice aggregation to snow threshold
             xmi50 = 4.8e-10, xmi40 = 2.46e-10,             &
             xni0 = 1.0e-2, xmnin = 1.05e-18, bni = 0.5,    &
             di50 = 1.0e-4, xmi = 4.19e-13,                 &   !--- parameters used in BF process
             bv_r = 0.8, bv_i = 0.25,                       &
             o6 = 1./6.,  cdrag = 0.6,                      &
             avisc = 1.49628e-6, adiffwv = 8.7602e-5,       &
             axka = 1.4132e3, cw = 4.187e3,  ci = 2.093e3

  REAL, DIMENSION(32), parameter :: a1 = (/                      &
              0.100e-10,0.7939e-7,0.7841e-6,0.3369e-5,0.4336e-5, &
              0.5285e-5,0.3728e-5,0.1852e-5,0.2991e-6,0.4248e-6, &
              0.7434e-6,0.1812e-5,0.4394e-5,0.9145e-5,0.1725e-4, &
              0.3348e-4,0.1725e-4,0.9175e-5,0.4412e-5,0.2252e-5, &
              0.9115e-6,0.4876e-6,0.3473e-6,0.4758e-6,0.6306e-6, &
              0.8573e-6,0.7868e-6,0.7192e-6,0.6513e-6,0.5956e-6, &
              0.5333e-6,0.4834e-6 /)

  REAL, DIMENSION(32), parameter :: a2 = (/                     &
              0.0100,0.4006,0.4831,0.5320,0.5307,0.5319,0.5249, &
              0.4888,0.3849,0.4047,0.4318,0.4771,0.5183,0.5463, &
              0.5651,0.5813,0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506,0.4483,0.4460, &
              0.4433,0.4413,0.4382,0.4361 /)
  
  REAL, DIMENSION(8), parameter :: B = (/                    &
             -.577191652,.988205891,-.897056937,.918206857,  &
             -.756704078,.482199394,-.193527818,.035868343 /)
  
  
contains

subroutine clphy1d_ylin(dt, imax,                           &
                      qvz, qlz, qrz, qiz, qsz,              &
                      thz, tothz, rho,                      &
                      prez, zz, dzw,                        &
!                      precrz, preciz, precsz,               & !zdc20220116
                      EFFC1D, EFFI1D, EFFS1D, EFFR1D,       & !zdc 20220208
                      pptrain, pptsnow,pptice,              &
                      kts, kte, Ri,                         &
                      ncz, nrz, niz, nsz,                   &
                      fluxr, fluxi, fluxs, fluxm,           &
                      fluxf, fevap, fsubl, fauto, fcoll,    &
                      faccr, vi,                            &
#ifdef diagnostic
                      zpsnow,zpsaut,zpsfw,zpsfi,zpraci,     & !process rate 
                      zpiacr,zpsaci,zpsacw,zpsdep,          &
                      zpssub,zpracs,zpsacr,zpsmlt,          &
                      zpsmltevp,zprain,zpraut,zpracw,       &
                      zprevp,zpgfr,zpvapor,zpclw,           &
                      zpladj,zpcli,zpimlt,zpihom,           &
                      zpidw,zpiadj,zqschg,                  &
#endif
                      zdrop,lin_aerosolmode)                  !aerosol feedback

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

  real, parameter :: sqrho=1.
! new 2D declaration
  integer,                         intent(in)    :: kts, kte
  integer,                         intent(in)    :: imax
  integer,                         intent(in)    :: lin_aerosolmode
  real,                            intent(in)    :: dt
  real, dimension(1:imax,kts:kte), intent(in)    :: zdrop
  real, dimension(1:imax,kts:kte), intent(inout) :: Ri
  real, dimension(1:imax,kts:kte), intent(in)    :: tothz,rho,              &
                                                    prez,dzw
  real, dimension(1:imax,0:kte), intent(in)      :: zz
  !real, dimension(1:imax,kts:kte), intent(out)   :: precrz,preciz,precsz
  real, dimension(1:imax,kts:kte), intent(out)   :: EFFC1D,EFFI1D,                     &
                                                    EFFS1D,EFFR1D
  real, dimension(1:imax,kts:kte), intent(out)   :: fluxr,fluxi,fluxs,                 &
                                                    fluxm,fluxf,fevap,fsubl,           &
                                                    fauto,fcoll,faccr
  real, dimension(1:imax,kts:kte), intent(out)   :: vi
#ifdef diagnostic
  real, dimension(1:imax,kts:kte), intent(out)   :: zpsnow,zpsaut,zpsfw,               &
                                                    zpsfi,zpraci,zpiacr,               &
                                                    zpsaci,zpsacw,zpsdep,              &
                                                    zpssub,zpracs,zpsacr,              &
                                                    zpsmlt,zpsmltevp,zprain,           &
                                                    zpraut,zpracw,zprevp,              &
                                                    zpgfr,zpvapor,zpclw,               &
                                                    zpladj,zpcli,zpimlt,               &
                                                    zpihom,zpidw,zpiadj,               &
                                                    zqschg
#endif
  real, dimension(1:imax),         intent(inout) :: pptrain, pptsnow, pptice
  real, dimension(1:imax,kts:kte), intent(inout) :: qvz,qlz,qrz,qiz,qsz,thz
  real, dimension(1:imax,kts:kte), intent(inout) :: ncz,niz,nrz,nsz
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer                                        :: k, iq 
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
                                                    qlpqi, rsat, aa1, aa2, xnin
  real, dimension(1:imax)                        :: es1d
  real, dimension(1:imax,kts:kte)                :: nczodt, nizodt, nrzodt, nszodt
  real, dimension(1:imax,kts:kte)                :: tem, temcc, theiz, qswz,    &
                                                    qsiz 
  real, dimension(1:imax)                        :: qvoqswz, qvoqsiz
  real, dimension(1:imax)                        :: qvzodt, qlzodt, qizodt, qszodt, qrzodt
  real                                           :: tmp2d
  real, dimension(1:imax,kts:kte)                :: rs0,viscmu,visc,diffwv,     &
                                                    schmidt,xka
  real, dimension(1:imax)                        :: qvsbar
!--- microphysical processes
  real, dimension(1:imax)                        :: psacw, psaut, psfw, psfi, praci
  real, dimension(1:imax)                        :: piacr, psaci, psdep, pssub, pracs
  real, dimension(1:imax)                        :: psacr, psmlt, psmltevp, prain, praut
  real, dimension(1:imax)                        :: pracw, prevp, pvapor, pclw, pladj
  real, dimension(1:imax)                        :: pcli, pimlt, pihom, pidw, piadj
  real, dimension(1:imax)                        :: pgfr, psnow, qschg
!---- new snow parameters
  real, parameter                               :: vf1s = 0.65,vf2s = 0.44,            &
                                                   vf1r =0.78,vf2r = 0.31 
  real, parameter                               :: am_c1=0.004,am_c2= 6e-5,  am_c3=0.15
  real, parameter                               :: bm_c1=1.85, bm_c2= 0.003, bm_c3=1.25
  real, parameter                               :: aa_c1=1.28, aa_c2= -0.012,aa_c3=-0.6
  real, parameter                               :: ba_c1=1.5,  ba_c2= 0.0075,ba_c3=0.5
  real, parameter                               :: best_a=1.08 ,  best_b = 0.499
  real                                          :: disp, Dc_liu, eta, R6c        !--- for Liu's autoconversion
  real, dimension(1:imax)                       :: tc0
  real, dimension(1:imax,kts:kte)               :: mu_c
  !real, dimension(kts:kte)                      :: ab_s,ab_r,ab_riming 
  real, dimension(1:imax,kts:kte)               :: cap_s    !---- capacitance of snow
  real, dimension(1:imax,kts:kte)               :: am_s,bm_s,av_s,bv_s,tmp_ss
  real, dimension(1:imax,kts:kte)               :: aa_s,ba_s,tmp_sa 
  real                                          :: mu_s=0.,mu_i=0.,mu_r=0.
 
  ! Adding variable Riz, which will duplicate Ri but be a copy passed upward
  real                                          :: episp0k, dtb, odtb, pi, pio4,       &
                                                   pio6, oxLf, xLvocp, xLfocp, av_r,   &
                                                   av_i, ocdrag, gambp4, gamdp4,       &
                                                   gam4pt5, Cpor, oxmi, gambp3, gamdp3,&
                                                   gambp6, gam3pt5, gam2pt75, gambp5o2,&
                                                   gamdp5o2, cwoxlf, ocp, xni50, es
  real                                          :: gam13
  real, parameter                               :: qvmin=1.e-20
  real                                          :: temc1,save1,save2,xni50mx
  real, dimension(1:imax,kts:kte)               :: vtr, vts,                           &
                                                   vtrold, vtsold, vtiold,             &
                                                   xlambdar, xlambdas,                 &
                                                   olambdar, olambdas
  ! for terminal velocity flux
  real                                          :: xmr,xms,xmc,dcs,xmr_i
  real                                          :: lamminr, lammaxr,lammins,           &
                                                   lammaxs,lammini, lammaxi
  real                                          :: gambvr1
  real                                          :: lvap
  real                                          :: mi0
  real t_del_tv,del_tv
  real flux, fluxin, fluxout
  real nflux, nfluxin, nfluxout
  integer                                       :: min_q, max_q
  logical                                       :: notlast
  logical, dimension(1:imax)                    :: mask
  real                                          :: nimlt, nihom
  real, dimension(1:imax)                       :: npraut_r, npracw, nprevp, npraut, npsmltevp
  real, dimension(1:imax)                       :: npsmlt, npgfr, npsacr, npsdep, npsacw, npsaci
  real, dimension(1:imax)                       :: npiacr, npraci, npsaut 
  real, dimension(1:imax)                       :: nidep, midep
  real, dimension(1:imax,kts:kte)               :: nvtr, nvts
  real, dimension(1:imax,kts:kte)               :: n0_s, n0_r
  real, dimension(1:imax)                       :: n0_i, n0_c                  
  real, dimension(1:imax,kts:kte)               :: lami, lamc
  real, dimension(1:imax)                       :: gg21, gg22
  real, dimension(1:imax)                       :: gg31, gg32, gg33
  real, dimension(1:imax,kts:kte)               :: gam_ss, gam_bm_s, gam_bv_ss
  real, dimension(1:imax,kts:kte)               :: gam_bv_s
  real ratio
  integer i1, i1p1  

  real, dimension(kts:kte)                      :: nrz_k, qrz_k, nsz_k, qsz_k, niz_k, qiz_k
  real, dimension(kts:kte)                      :: vtrold_k, nvtr_k, olambdar_k
  real, dimension(kts:kte)                      :: xlambdar_k, n0_r_k 
  real, dimension(kts:kte)                      :: vtsold_k, nvts_k, olambdas_k
  real, dimension(kts:kte)                      :: xlambdas_k, n0_s_k
  real, dimension(kts:kte)                      :: vtiold_k
  real, dimension(kts:kte)                      :: rho_k, dzw_k
  real, dimension(0:kte)                        :: zz_k
  real, dimension(kts:kte)                      :: am_s_k, av_s_k, bm_s_k, bv_s_k, gam_bm_s_k, gam_bv_s_k, gam_bv_ss_k, gam_ss_k 
  real                                          :: pptrain_k, pptsnow_k, pptice_k

  !------------------------------------------------------------------------------------

  vtrold=0.
  vtsold=0.
  vtiold=0.

  mu_c    = MIN(15., (1000.E6/Nt_c + 2.))
  R6c     = 10.E-6      !---- 10 micron, threshold radius of cloud droplet
  dtb     = dt                                                                         !sny
  odtb    =1./dtb
  pi      =acos(-1.)
  pio4    =pi/4.
  pio6    =pi/6.
  ocp     =1./cp
  oxLf    =1./xLf
  xLvocp  =xLv/cp
  xLfocp  =xLf/cp
  Cpor    =cp/Rair
  oxmi    =1./xmi
  cwoxlf  =cw/xlf 
  av_r    =2115.0*0.01**(1.-bv_r)
  av_i    =152.93*0.01**(1.-bv_i)
  ocdrag  =1./Cdrag
  episp0k =RH*ep2*1000.*svp1

  gambp4  =ggamma(bv_r+4.)
  gamdp4  =ggamma(bv_i+4.)
  gambp3  =ggamma(bv_r+3.)
  gambp6  =ggamma(bv_r+6.)
  gambp5o2=ggamma((bv_r+5.)/2.)
  gamdp5o2=ggamma((bv_i+5.)/2.)
  gambvr1 =ggamma(bv_r+1.)
  gam13   =ggamma(1.+3.)
  !
  !     qsw         saturated mixing ratio on water surface
  !     qsi         saturated mixing ratio on ice surface
  !     episp0k     RH*e*saturated pressure at 273.15 K  = 611.2 hPa (Rogers and Yau 1989)
  !     qvoqsw      qv/qsw
  !     qvoqsi      qv/qsi
  !     qvzodt      qv/dt
  !     qlzodt      ql/dt
  !     qizodt      qi/dt
  !     qszodt      qs/dt

  obp4    =1./(bv_r+4.)
  bp3     =bv_r+3.
  bp5     =bv_r+5.
  bp6     =bv_r+6.
  odp4    =1./(bv_i+4.)
  dp3     =bv_i+3.
  dp5     =bv_i+5.
  dp5o2   =0.5*(bv_i+5.)
    
  dcs     = 125.E-6  ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
  xms     =pi*500.*(dcs)**3/6.   !morr =PI*RHOI*DCS**3/6.=5.11*10e-10
  xmr     =4./3.*pi*rhowater*(500.E-6)**3
  xmr_i   =4./3.*pi*rhowater*(25.E-6)**3
  mi0     = 4./3.*3.14*500.*(10.e-6)**3
  !    xmc     =4.17*10e-14 !4./3.*pi*(0.00001)**3*1000.
  lammaxr = 1./20.E-6
  !lamminr = 1./500.E-6
  lamminr = 1./2800.E-6
  lammaxs = 1./10.E-6
  lammins = 1./2000.E-6
  lammaxi = 1./1.E-6
  lammini = 1./(2.*dcs+100.E-6)

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
    nrz(1:imax,k) = max( 0.0,nrz(1:imax,k) )
    nsz(1:imax,k) = max( 0.0,nsz(1:imax,k) )
    nczodt(1:imax,k)=max( 0.0,odtb*ncz(1:imax,k) )
    nizodt(1:imax,k)=max( 0.0,odtb*niz(1:imax,k) )
    nrzodt(1:imax,k)=max( 0.0,odtb*nrz(1:imax,k) )
    nszodt(1:imax,k)=max( 0.0,odtb*nsz(1:imax,k) )
  end do

  do k=kts,kte
    qlz(1:imax,k)  =max( 0.0,qlz(1:imax,k) )
    qiz(1:imax,k)  =max( 0.0,qiz(1:imax,k) )
    qvz(1:imax,k)  =max( qvmin,qvz(1:imax,k) )
    qsz(1:imax,k)  =max( 0.0,qsz(1:imax,k) )
    qrz(1:imax,k)  =max( 0.0,qrz(1:imax,k) )
    tem(1:imax,k)  =thz(1:imax,k)*tothz(1:imax,k)
    temcc(1:imax,k)=tem(1:imax,k)-273.15
    es1d(1:imax)   =1000.*svp1*exp( svp2*temcc(1:imax,k)/(tem(1:imax,k)-svp3) )  !--- RY89 Eq(2.17)
    qswz(1:imax,k) =ep2*es1d(1:imax)/max(prez(1:imax,k)-es1d(1:imax),0.1)
 
    do iq=1,imax
      if (tem(iq,k) .lt. 233.15 ) then
        es1d(iq)=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
        qsiz(iq,k)=ep2*es1d(iq)/max(prez(iq,k)-es1d(iq),0.1)
        qswz(iq,k)=qsiz(iq,k)
      else if (tem(iq,k) .lt. 273.15 ) then
        es1d(iq)=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
        qsiz(iq,k)=ep2*es1d(iq)/max(prez(iq,k)-es1d(iq),0.1)
      else
        qsiz(iq,k)=qswz(iq,k)
      endif
    enddo

    theiz(1:imax,k)=thz(1:imax,k)+(xlvocp*qvz(1:imax,k)-xlfocp*qiz(1:imax,k))/tothz(1:imax,k)
  enddo

  do k=kts,kte

    n0_s(1:imax,k)     =0.
    lamc(1:imax,k)     =0.
    lami(1:imax,k)     =0.
    xlambdar(1:imax,k) =0.
    xlambdas(1:imax,k) =0.
    vtr(1:imax,k)      =0.
    vts(1:imax,k)      =0.
    vtiold(1:imax,k)   =0.

    !qisten(1:imax,k)   =0.
    !qrsten(1:imax,k)   =0.
    !qssten(1:imax,k)   =0.
    !nisten(1:imax,k)   =0.
    !nrsten(1:imax,k)   =0.
    !nssten(1:imax,k)   =0.
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
    viscmu(1:imax,k)=avisc*tem(1:imax,k)**1.5/(tem(1:imax,k)+120.0)
    visc(1:imax,k)=viscmu(1:imax,k)/rho(1:imax,k)
    diffwv(1:imax,k)=adiffwv*tem(1:imax,k)**1.81/prez(1:imax,k)
    schmidt(1:imax,k)=visc(1:imax,k)/diffwv(1:imax,k)
    xka(1:imax,k)=axka*viscmu(1:imax,k)
    rs0(1:imax,k)=ep2*1000.*svp1/(prez(1:imax,k)-1000.*svp1)
  end do

  ! ---- YLIN, set snow variables
  !
  !---- A+B in depositional growth, the first try just take from Rogers and Yau(1989)
  !         ab_s(k) = lsub*lsub*orv/(tcond(k)*temp(k))+&
  !                   rv*temp(k)/(diffu(k)*qvsi(k))

  do k = kts, kte
    !tc0(1:imax)   = tem(1:imax,k)-273.15
    where( rho(1:imax,k)*qlz(1:imax,k) .gt. 1e-5 .AND. rho(1:imax,k)*qsz(1:imax,k) .gt. 1e-5  )
      Ri(1:imax,k) = 1.0/(1.0+6e-5/(rho(1:imax,k)**1.170*qlz(1:imax,k)*qsz(1:imax,k)**0.170)) 
    elsewhere
      Ri(1:imax,k) = 0.
    end where
  end do


  !
  !--- make sure Ri does not decrease downward in a column
  !
  do k = kte-1,kts,-1
    do iq = 1,imax
      Ri(iq,k) = max( Ri(iq,k), Ri(iq,k+1) )
    end do
  end do

  !--- YLIN, get PI properties
  do k = kts, kte
    Ri(1:imax,k) = max(0.,min(Ri(1:imax,k),1.0))

    cap_s(1:imax,k)= 0.25*(1.+Ri(1:imax,k))
    tc0(1:imax)    = min(-0.1, tem(1:imax,k)-273.15)
    n0_s(1:imax,k) = min(2.0E8, 2.0E6*exp(-0.12*tc0(1:imax)))
    am_s(1:imax,k) = am_c1+am_c2*tc0(1:imax)+am_c3*Ri(1:imax,k)*Ri(1:imax,k)   !--- Heymsfield 2007
    am_s(1:imax,k) = max(0.000023,am_s(1:imax,k))                                !--- use the a_min in table 1 of Heymsfield
    bm_s(1:imax,k) = bm_c1+bm_c2*tc0(1:imax)+bm_c3*Ri(1:imax,k)
    bm_s(1:imax,k) = min(bm_s(1:imax,k),3.0)                                     !---- capped by 3
    !--  converting from cgs to SI unit
    am_s(1:imax,k) =  10.**(2.*bm_s(1:imax,k)-3.)*am_s(1:imax,k)
    aa_s(1:imax,k) = aa_c1 + aa_c2*tc0(1:imax) + aa_c3*Ri(1:imax,k)
    ba_s(1:imax,k) = ba_c1 + ba_c2*tc0(1:imax) + ba_c3*Ri(1:imax,k)
    !--  convert to SI unit as in paper
    aa_s(1:imax,k) = (1.e-2)**(2.-ba_s(1:imax,k))*aa_s(1:imax,k)
    !---- get v from Mitchell 1996
    av_s(1:imax,k) = best_a*viscmu(1:imax,k)*(2*grav*am_s(1:imax,k)/rho(1:imax,k)/ &
                       aa_s(1:imax,k)/(viscmu(1:imax,k)**2))**best_b
    bv_s(1:imax,k) = best_b*(bm_s(1:imax,k)-ba_s(1:imax,k)+2.)-1.

    tmp_ss(1:imax,k)= bm_s(1:imax,k)+mu_s+1.
    tmp_sa(1:imax,k)= ba_s(1:imax,k)+mu_s+1.
  end do
  
  do k = kts,kte
    do iq = 1,imax
      gam_ss(iq,k) = ggamma(tmp_ss(iq,k)) 
      gam_bm_s(iq,k) = ggamma(1.+bm_s(iq,k))
      gam_bv_ss(iq,k) = ggamma(bv_s(iq,k)+tmp_ss(iq,k))
      gam_bv_s(iq,k) = ggamma(bv_s(iq,k)+1.)
    end do
  end do

  !***********************************************************************
  ! Calculate precipitation fluxes due to terminal velocities.
  !***********************************************************************
  !
  !- Calculate termianl velocity (vt?)  of precipitation q?z
  !- Find maximum vt? to determine the small delta t
  !
  !-- rain

  do iq = 1,imax
    qrz_k(kts:kte) = qrz(iq,kts:kte)    
    notlast=any( qrz_k(kts:kte)>1.e-8 )
    if ( notlast ) then
      vtrold_k(kts:kte) = 0.
      n0_r_k(kts:kte) = 0.
      nrz_k(kts:kte) = nrz(iq,kts:kte)
      pptrain_k = pptrain(iq)
      zz_k(0:kte) = zz(iq,0:kte)
      rho_k(kts:kte) = rho(iq,kts:kte)
      dzw_k(kts:kte) = dzw(iq,kts:kte)
      t_del_tv=0.
      del_tv=dtb
      do while (notlast)
        
        min_q=kte
        max_q=kts-1

        ! if no rain, --> minq>maxq --> notlast=False (only case minq>maxq)
        ! if rain --> minq<maxq (always), some vertical points norain--> no lamda, velocity
        do k = kts,kte-1
          if ( qrz_k(k) > 1.e-8 ) then
            min_q = min(min_q,k)
            max_q = max(max_q,k)
            xlambdar_k(k) = (pi*rhowater*nrz_k(k)/qrz_k(k))**(1./3.)   !zx
            n0_r_k(k) = nrz_k(k)*xlambdar_k(k)
            if ( xlambdar_k(k) < lamminr ) then
              xlambdar_k(k) = lamminr
              n0_r_k(k) = xlambdar_k(k)**4*qrz_k(k)/(pi*rhowater)
              nrz_k(k) = n0_r_k(k)/xlambdar_k(k) 
            else if ( xlambdar_k(k) > lammaxr ) then
              xlambdar_k(k) = lammaxr
              n0_r_k(k) = xlambdar_k(k)**4*qrz_k(k)/(pi*rhowater)
              nrz_k(k) = n0_r_k(k)/xlambdar_k(k)
            end if
            olambdar_k(k) = 1./xlambdar_k(k)
            tmp1 = olambdar_k(k)**bv_r
            vtrold_k(k) = o6*av_r*gambp4*sqrho*tmp1
            nvtr_k(k) = av_r*gambvr1*sqrho*tmp1
            del_tv = min(del_tv,0.9*(zz_k(k)-zz_k(k-1))/vtrold_k(k))
          else
            vtrold_k(k)=0.
            nvtr_k(k)=0.
            olambdar_k(k)=0.
          end if
        end do         ! k
      
        !
        !- Check if the summation of the small delta t >=  big delta t
        !             (t_del_tv)          (del_tv)             (dtb)

        if (max_q >= min_q) then

          fluxin=0.
          nfluxin=0. ! sny
          t_del_tv=t_del_tv+del_tv
          if ( t_del_tv >= dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
          end if

          do k = max_q,min_q,-1
            fluxout=rho_k(k)*vtrold_k(k)*qrz_k(k)
            flux=(fluxin-fluxout)/rho_k(k)/dzw_k(k)
            !tmpqrz(iq)=qrz(iq,k)
            qrz_k(k)=qrz_k(k)+del_tv*flux
            fluxin=fluxout

            nfluxout=rho_k(k)*nvtr_k(k)*nrz_k(k)
            nflux=(nfluxin-nfluxout)/rho_k(k)/dzw_k(k)
            nrz_k(k)=nrz_k(k)+del_tv*nflux
            nfluxin=nfluxout
            !qrsten(iq,k)=flux
            !nrsten(iq,k)=nflux
          end do       !k

          if ( min_q == 1 ) then
            pptrain_k = pptrain_k + fluxin*del_tv  
          else      
            qrz_k(min_q-1)=qrz_k(min_q-1)+del_tv*  &
                           fluxin/rho_k(min_q-1)/dzw_k(min_q-1)
            nrz_k(min_q-1)=nrz_k(min_q-1)+del_tv*  &
                           nfluxin/rho_k(min_q-1)/dzw_k(min_q-1)
          end if

        else
          notlast=.false.
        end if ! maxq>minq

      END DO      ! while(notlast)
      vtrold(iq,kts:kte)   = vtrold_k(kts:kte)
      n0_r(iq,kts:kte)     = n0_r_k(kts:kte)
      nrz(iq,kts:kte)      = nrz_k(kts:kte)
      qrz(iq,kts:kte)      = qrz_k(kts:kte)
      pptrain(iq)          = pptrain_k
    end if
  end do ! iq

  !
  !-- snow
  !
  do iq = 1,imax
    qsz_k(kts:kte) = qsz(iq,kts:kte)    
    notlast=any( qsz_k(kts:kte)>1.e-8 )
    if ( notlast ) then
      vtsold_k(kts:kte) = 0.
      n0_s_k(kts:kte) = n0_s(iq,kts:kte)
      nsz_k(kts:kte) = nsz(iq,kts:kte)
      pptsnow_k = pptsnow(iq)
      zz_k(0:kte) = zz(iq,0:kte)
      rho_k(kts:kte) = rho(iq,kts:kte)
      dzw_k(kts:kte) = dzw(iq,kts:kte)
      am_s_k(kts:kte) = am_s(iq,kts:kte)
      av_s_k(kts:kte) = av_s(iq,kts:kte)
      bm_s_k(kts:kte) = bm_s(iq,kts:kte)
      bv_s_k(kts:kte) = bv_s(iq,kts:kte)
      gam_bm_s_k(kts:kte) = gam_bm_s(iq,kts:kte)
      gam_bv_s_k(kts:kte) = gam_bv_s(iq,kts:kte)
      gam_bv_ss_k(kts:kte) = gam_bv_ss(iq,kts:kte)
      gam_ss_k(kts:kte) = gam_ss(iq,kts:kte)
      t_del_tv=0.
      del_tv=dtb
      DO while (notlast)

        min_q=kte
        max_q=kts-1

        do k=kts,kte-1
          if (qsz_k(k) > 1.e-8) then
            min_q = min(min_q,k)
            max_q = max(max_q,k)
            ! Zhao 2022 - Row 2 Table 2 or Lin 2011 - Formula A3
            xlambdas_k(k)=(am_s_k(k)*gam_ss_k(k)*nsz_k(k)/qsz_k(k))**(1./bm_s_k(k))
            ! Zhao 2022 - Row 1 Table 2
            n0_s_k(k)=nsz_k(k)*xlambdas_k(k)
            if (xlambdas_k(k)<lammins) then
              xlambdas_k(k)= lammins
              n0_s_k(k) = xlambdas_k(k)**(bm_s_k(k)+1.)*qsz_k(k)/gam_bm_s_k(k)/am_s_k(k)
              nsz_k(k) = n0_s_k(k)/xlambdas_k(k)
            else if (xlambdas_k(k)>lammaxs) then
              xlambdas_k(k) = lammaxs
              n0_s_k(k) = xlambdas_k(k)**(bm_s_k(k)+1.)*qsz_k(k)/gam_bm_s_k(k)/am_s_k(k)
              nsz_k(k) = n0_s_k(k)/xlambdas_k(k)
            end if
            olambdas_k(k)=1./xlambdas_k(k)
            tmp1 = olambdas_k(k)**bv_s_k(k)
            ! Zhao 2022 - Row 3 Table 2
            vtsold_k(k)= sqrho*av_s_k(k)*gam_bv_ss_k(k)/ &
               gam_ss_k(k)*tmp1
            ! Zhao 2022 - Row 4 Table 2
            nvts_k(k)=sqrho*av_s_k(k)*gam_bv_s_k(k)*tmp1
            del_tv=min(del_tv,0.9*(zz_k(k)-zz_k(k-1))/vtsold_k(k))
          else
            vtsold_k(k)=0.
            nvts_k(k)=0.
            olambdas_k(k)=0.
          endif
        end do       ! k

        !
        !- Check if the summation of the small delta t >=  big delta t
        !             (t_del_tv)          (del_tv)             (dtb)

        if (max_q >= min_q) then

          fluxin = 0.
          nfluxin = 0.
          t_del_tv=t_del_tv+del_tv
          if ( t_del_tv >= dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
          endif

          do k = max_q,min_q,-1
            fluxout=rho_k(k)*vtsold_k(k)*qsz_k(k)
            flux=(fluxin-fluxout)/rho_k(k)/dzw_k(k)
            qsz_k(k)=qsz_k(k)+del_tv*flux
            qsz_k(k)=max(0.,qsz_k(k))
            fluxin=fluxout

            nfluxout=rho_k(k)*nvts_k(k)*nsz_k(k)
            nflux   =(nfluxin-nfluxout)/rho_k(k)/dzw_k(k)
            nsz_k(k)  =nsz_k(k)+del_tv*nflux
            nfluxin =nfluxout
            !qssten(iq,k)=flux
            !nssten(iq,k)=nflux
          end do       ! k

          if ( min_q == 1 ) then
            pptsnow_k = pptsnow_k + fluxin*del_tv  
          else
            qsz_k(min_q-1)=qsz_k(min_q-1)+del_tv*  &
                     fluxin/rho_k(min_q-1)/dzw_k(min_q-1)
            nsz_k(min_q-1)=nsz_k(min_q-1)+del_tv*  &
                       nfluxin/rho_k(min_q-1)/dzw_k(min_q-1)
          end if

        else
          notlast=.false.
        end if ! maxq>minq

      END DO       ! while(notlast)
      vtsold(iq,kts:kte)   = vtsold_k(kts:kte)
      n0_s(iq,kts:kte)     = n0_s_k(kts:kte)
      nsz(iq,kts:kte)      = nsz_k(kts:kte)
      qsz(iq,kts:kte)      = qsz_k(kts:kte)
      pptsnow(iq)          = pptsnow_k
    end if
  end do        ! iq = 1,imax
 
  !
  !-- cloud ice  (03/21/02) using Heymsfield and Donner (1990) Vi=3.29*qi^0.16
  !
  do iq = 1,imax
    qiz_k(kts:kte) = qiz(iq,kts:kte)
    notlast=any( qiz_k(kts:kte)>1.e-8 )
    if ( notlast ) then
      vtiold_k(kts:kte) = 0.
      niz_k(kts:kte) = niz(iq,kts:kte)
      pptice_k = pptice(iq)
      zz_k(0:kte) = zz(iq,0:kte)
      rho_k(kts:kte) = rho(iq,kts:kte)
      dzw_k(kts:kte) = dzw(iq,kts:kte)
      t_del_tv=0.
      del_tv=dtb
      DO while (notlast)

        min_q=kte
        max_q=kts-1

        do k=kts,kte-1
          if (qiz_k(k) > 1.e-8) then
            min_q = min(min_q,k)
            max_q = max(max_q,k)
            vtiold_k(k) = 3.29 * (rho_k(k)* qiz_k(k))** 0.16  ! Heymsfield and Donner
            del_tv=min(del_tv,0.9*(zz_k(k)-zz_k(k-1))/vtiold_k(k))
          else
            vtiold_k(k)=0.
          endif
        enddo       ! k
      
        !
        !- Check if the summation of the small delta t >=  big delta t
        !             (t_del_tv)          (del_tv)             (dtb)
        if (max_q >= min_q) then

          fluxin = 0.
          nfluxin = 0.
          t_del_tv=t_del_tv+del_tv
          if ( t_del_tv >= dtb ) then
            notlast=.false.
            del_tv=dtb+del_tv-t_del_tv
          endif

          do k = max_q,min_q,-1
            fluxout=rho_k(k)*vtiold_k(k)*qiz_k(k)
            flux=(fluxin-fluxout)/rho_k(k)/dzw_k(k)
            qiz_k(k)=qiz_k(k)+del_tv*flux
            qiz_k(k)=max(0.,qiz_k(k))
            fluxin=fluxout

            nfluxout=rho_k(k)*vtiold_k(k)*niz_k(k)
            nflux=(nfluxin-nfluxout)/rho_k(k)/dzw_k(k)
            niz_k(k)=niz_k(k)+del_tv*nflux
            niz_k(k)=max(0.,niz_k(k))
            nfluxin=nfluxout
            !qisten(iq,k)=flux
            !nisten(iq,k)=nflux
          end do       ! k

          if ( min_q == 1 ) then
            pptice_k = pptice_k + fluxin*del_tv  
          else
            qiz_k(min_q-1)=qiz_k(min_q-1)+del_tv*  &
                         fluxin/rho_k(min_q-1)/dzw_k(min_q-1)
            niz_k(min_q-1)=niz_k(min_q-1)+del_tv*  &
                         nfluxin/rho_k(min_q-1)/dzw_k(min_q-1)
          end if

        else
          notlast=.false.
        end if   ! maxq>minq

    
      END DO       ! while(notlast)
      vtiold(iq,kts:kte)   = vtiold_k(kts:kte)
      niz(iq,kts:kte)      = niz_k(kts:kte)
      qiz(iq,kts:kte)      = qiz_k(kts:kte)
      pptice(iq)           = pptice_k
    end if
  end do        ! iq = 1,imax
    
  ! Microphpysics processes
  DO k=kts,kte

    qvoqswz(1:imax)  =qvz(1:imax,k)/qswz(1:imax,k)
    qvoqsiz(1:imax)  =qvz(1:imax,k)/qsiz(1:imax,k)
    qvzodt(1:imax)=max( 0.,odtb*qvz(1:imax,k) )
    qlzodt(1:imax)=max( 0.,odtb*qlz(1:imax,k) )
    qizodt(1:imax)=max( 0.,odtb*qiz(1:imax,k) )
    qszodt(1:imax)=max( 0.,odtb*qsz(1:imax,k) )
    qrzodt(1:imax)=max( 0.,odtb*qrz(1:imax,k) )

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

    n0_r(1:imax,k)     = 0.
    n0_i(1:imax)       = 0.
    n0_c(1:imax)       = 0.

    nvtr(1:imax,k)     = 0.
    nvts(1:imax,k)     = 0.

    xlambdar(1:imax,k) = 0.
    xlambdas(1:imax,k) = 0.

    olambdar(1:imax,k) = 0.
    olambdas(1:imax,k) = 0.

    psacw(1:imax)   =0.                  ! accretion of cloud water by snow
    psaut(1:imax)   =0.                  ! ice crystal aggregation to snow
    psfw(1:imax)    =0.                  ! BERGERON process to transfer cloud water to snow
    psfi(1:imax)    =0.                  ! BERGERON process to transfer cloud ice to snow
    praci(1:imax)   =0.                  ! cloud ice accretion by rain
    piacr(1:imax)   =0.                  ! rain accretion by cloud ice
    psaci(1:imax)   =0.                  ! ice crystal accretion by snow
    psdep(1:imax)   =0.                  ! deposition of snow
    pssub(1:imax)   =0.                  ! sublimation of snow (T<0)
    psacr(1:imax)   =0.                  ! accretion of rain by snow
    psmlt(1:imax)   =0.                  ! melting of snow
    psmltevp(1:imax)=0.                  ! evaporation of melting snow (T>0)
    prain(1:imax)   =0.                  ! sum all process for rain
    praut(1:imax)   =0.                  ! autoconversion of rain
    pracw(1:imax)   =0.                  ! accretion of cloud water by rain
    prevp(1:imax)   =0.                  ! evaporation of rain
    pvapor(1:imax)  =0.                  ! sum all process for water vapor to determine qvz
    pclw(1:imax)    =0.                  ! sum all process for cloud liquid to determine qlz
    pladj(1:imax)   =0.                  ! saturation adjustment for ql
    pcli(1:imax)    =0.                  ! sum all process for cloud ice to determine qiz
    pimlt(1:imax)   =0.                  ! melting of ice crystal >0.
    pihom(1:imax)   =0.                  ! homogeneous nucleation <-40
    pidw(1:imax)    =0.                  ! production of cloud ice by BERGERON process
    piadj(1:imax)   =0.                  ! saturation adjustment for qi
    pgfr(1:imax)    =0.                  ! feezing of rain to form graupel (added to PI)
    psnow(1:imax)   =0.                  ! sum all process for snow
    qschg(1:imax)   =0.                  ! = psnow / unsure
    pracs(1:imax)   =0.

    nidep(:) = 0.
    midep(:) = 0.
    npraut_r(:) = 0.
    nprevp(:) = 0.
    npracw(:) = 0.
    npraut(:) = 0.
    npsmltevp(:) = 0.
    npsmlt(:) = 0.
    npgfr(:) = 0.
    npsacr(:) = 0.
    npsdep(:) = 0.
    npsacw(:) = 0.
    npsaci(:) = 0.
    npiacr(:) = 0.
    npraci(:) = 0.
    npsaut(:) = 0.
    !npssub(1:imax,k)   =0.
    
    fluxr(1:imax,k) = qrzodt(1:imax)*dzw(1:imax,k)*rho(1:imax,k) ! MJT suggestion
    fluxi(1:imax,k) = qizodt(1:imax)*dzw(1:imax,k)*rho(1:imax,k) ! MJT suggestion
    fluxs(1:imax,k) = qszodt(1:imax)*dzw(1:imax,k)*rho(1:imax,k) ! MJT suggestion
    
    do iq = 1,imax
      tmp2d=qiz(iq,k)+qlz(iq,k)+qsz(iq,k)+qrz(iq,k)
      mask(iq) = .not.(qvz(iq,k)+qlz(iq,k)+qiz(iq,k) < qsiz(iq,k) .and. &
                       tmp2d == 0.)
      
      if ( mask(iq) ) then
        !gg11(iq) = ggamma(tmp_ss(iq,k))
        !gg12(iq) = ggamma(1.+bm_s(iq,k))
        !gg13(iq) = ggamma(bv_s(iq,k)+tmp_ss(iq,k))
        !gg14(iq) = ggamma(bv_s(iq,k)+1.)
        gg21(iq) = ggamma(bv_s(iq,k)+tmp_sa(iq,k))
        gg22(iq) = ggamma(2.5+0.5*bv_s(iq,k)+mu_s)
        gg31(iq) = ggamma(4.+mu_c(iq,k))
        gg32(iq) = ggamma(1.+mu_c(iq,k))
        gg33(iq) = ggamma(6.+1.+mu_c(iq,k))
        !gg41(iq) = ggamma(tmp_ss(iq,k))
        !gg42(iq) = ggamma(1.+bm_s(iq,k))
        !gg51(iq) = ggamma(4.+mu_c(iq,k))
        !gg52(iq) = ggamma(1.+mu_c(iq,k))      
      end if  
    end do  

    !
    !! calculate terminal velocity of rain
    !    
    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000

        if (qrz(iq,k) .gt. 1.e-8) then
        !  tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrz(k))
        !  xlambdar(k)=sqrt(tmp1)
        !  olambdar(k)=1.0/xlambdar(k)
        !  vtrold(k)=o6*av_r*gambp4*sqrho*olambdar(k)**bv_r
          xlambdar(iq,k)=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
          n0_r(iq,k)=nrz(iq,k)*xlambdar(iq,k)
          if (xlambdar(iq,k)<lamminr) then
            xlambdar(iq,k) = lamminr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k)
          else if (xlambdar(iq,k)>lammaxr) then
            xlambdar(iq,k) = lammaxr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,k)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,k)
          end if
          olambdar(iq,k)=1.0/xlambdar(iq,k)
          tmp1 = olambdar(iq,k)**bv_r
          vtrold(iq,k)=o6*av_r*gambp4*sqrho*tmp1
          nvtr(iq,k)=av_r*gambvr1*sqrho*tmp1
        else
          vtrold(iq,k)=0.
          olambdar(iq,k)=0.
          nvtr(iq,k)=0.
        end if  ! qrz

        vtr(iq,k)=vtrold(iq,k)
        
      end if
    end do
    
    !!
    !!! calculate terminal velocity of snow
    !!
   
    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000

        if (qsz(iq,k) > 1.e-8) then
          xlambdas(iq,k)=(am_s(iq,k)*gam_ss(iq,k)*nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
          n0_s(iq,k)=nsz(iq,k)*xlambdas(iq,k)
          if (xlambdas(iq,k).lt.lammins) then
            xlambdas(iq,k)= lamminS
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          else if (xlambdas(iq,k).gt.lammaxs) then
            xlambdas(iq,k) = lammaxs
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1)*qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          end if
          olambdas(iq,k)=1.0/xlambdas(iq,k)
          tmp1 = olambdas(iq,k)**bv_s(iq,k)
          vtsold(iq,k)= sqrho*av_s(iq,k)*gam_bv_ss(iq,k)/ &
                 gam_ss(iq,k)*tmp1
          nvts(iq,k)=sqrho*av_s(iq,k)*gam_bv_s(iq,k)*tmp1
        else
          vtsold(iq,k)=0.
          olambdas(iq,k)=0.
          xlambdas(iq,k)=0.
          nvts(iq,k)=0.
        endif

        vts(iq,k)=vtsold(iq,k)
        
      end if
    end do
    
    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
        
        !---------- start of snow/ice processes below freezing

        if (tem(iq,k) < 273.15) then

        !
        ! ice nucleation, cooper curve

          if ((qvoqswz(iq)>=0.999.and.temcc(iq,k)<=-8.).or. &
            qvoqsiz(iq)>=1.08) then
            nidep(iq) = 5.*exp(0.304*(273.15-tem(iq,k))) ! m-3
            nidep(iq) = min(nidep(iq), 500.e3)               !5.e8) sny ! limit to 500 L-1
            nidep(iq) = max(nidep(iq)/rho(iq,k), 0.)       ! convert to kg-1
            nidep(iq) = (nidep(iq) - niz(iq,k))*odtb
            midep(iq) = nidep(iq)*mi0
          end if
          !***********************************************************************
          !*********        snow production processes for T < 0 C       **********
          !***********************************************************************
          !
          ! (1) ICE CRYSTAL AGGREGATION TO SNOW (Psaut): Lin (21)
          !!    psaut=alpha1*(qi-qi0)
          !!    alpha1=1.0e-3*exp(0.025*(T-T0))
          !
          alpha1=1.e-3*exp( 0.025*temcc(iq,k) )

          ! ---------------------------------------------------------------
          if(temcc(iq,k) .lt. -20.0) then
              tmp1=-7.6+4.*exp( -0.2443e-3*(abs(temcc(iq,k))-20.)**2.455 )
              qic=1.e-3*exp(tmp1)/rho(iq,k)
          else
              qic=qi0
          end if
          !----------------------------------------------------------------

          !qic = qi0  ! sny: OFF temp

          tmp1=odtb*(qiz(iq,k)-qic)*(1.0-exp(-alpha1*dtb))
          psaut(iq)=max( 0., tmp1 )
          npsaut(iq)=max( 0., psaut(iq)/xms)

          !
          ! (2) BERGERON PROCESS TRANSFER OF CLOUD WATER TO SNOW (Psfw)
          !     this process only considered when -31 C < T < 0 C
          !     Lin (33) and Hsie (17)
          ! 
          !!
          !!    parama1 and parama2 functions must be user supplied
          !!

          if( qlz(iq,k) > 1.e-10 ) then
            temc1=max(-30.99,temcc(iq,k))
            ! aa1=parama1(temc1)
            ! aa2=parama2(temc1)
            i1=int(-temc1)+1
            i1p1=i1+1
            ratio=-(temc1)-real(i1-1)
            aa1=a1(i1)+ratio*( a1(i1p1)-a1(i1) )
            aa2=a2(i1)+ratio*( a2(i1p1)-a2(i1) )
            
            tmp1=1.-aa2
            !!   change unit from cgs to mks
            aa1=aa1*0.001**tmp1
            !!   dtberg is the time needed for a crystal to grow from 40 to 50 um
            !!   odtberg=1.0/dtberg
            odtberg=(aa1*tmp1)/(xmi50**tmp1-xmi40**tmp1)
            !
            !!   compute terminal velocity of a 50 micron ice cystal
            !
            vti50=av_i*di50**bv_i*sqrho
            eiw=1.0
            save1=aa1*xmi50**aa2
            save2=0.25*pi*eiw*rho(iq,k)*di50*di50*vti50
            tmp2=( save1 + save2*qlz(iq,k) )
            !
            !!  maximum number of 50 micron crystals limited by the amount
            !!  of supercool water
            !
            xni50mx=qlzodt(iq)/tmp2
            !
            !!   number of 50 micron crystals produced
            !
            xni50=qiz(iq,k)*( 1.0-exp(-dtb*odtberg) )/xmi50
            xni50=min(xni50,xni50mx)
            !
            tmp3=odtb*tmp2/save2*( 1.0-exp(-save2*xni50*dtb) )
            psfw(iq)=min( tmp3,qlzodt(iq) )
            !
            ! (3) REDUCTION OF CLOUD ICE BY BERGERON PROCESS (Psfi): Lin (34)
            !     this process only considered when -31 C < T < 0 C
            !
            tmp1=xni50*xmi50-psfw(iq)
            psfi(iq)=min(tmp1,qizodt(iq))
          end if
        
          if(qrz(iq,k) > 0.) then  ! go to 1000

          !
          ! Processes (4) and (5) only need when qrz > 0.0
          !
          ! (4) CLOUD ICE ACCRETION BY RAIN (Praci): Lin (25)
          !     produce PI
          !
            eri=1.0
            save1=pio4*eri*n0_r(iq,k)*av_r*sqrho
            tmp1= save1*gambp3*olambdar(iq,k)**bp3
            praci(iq)=qizodt(iq)*( 1.0-exp(-tmp1*dtb) )
            npraci(iq)=niz(iq,k)*tmp1

          !
          ! (5) RAIN ACCRETION BY CLOUD ICE (Piacr): Lin (26)
          !
            tmp2=qiz(iq,k)*save1*rho(iq,k)*pio6*rhowater*gambp6*oxmi* &
                olambdar(iq,k)**bp6
            piacr(iq)=min( tmp2,qrzodt(iq) )
            tmp1 = olambdar(iq,k)**bp3 
            npiacr(iq)=pio4*eri*nrz(iq,k)*av_r*niz(iq,k)*gambp3* &
                  tmp1 !--wdm6
          end if !1000    continue

          if(qsz(iq,k) > 0.) then !go to 1200
          !
          ! Compute the following processes only when qsz > 0.0
          !
          !
          ! (6) ICE CRYSTAL ACCRETION BY SNOW (Psaci): Lin (22)
          !
            esi=exp( 0.025*temcc(iq,k) )
            tmp1 = olambdas(iq,k)**(bv_s(iq,k)+tmp_sa(iq,k))
            save1 = aa_s(iq,k)*sqrho*n0_s(iq,k)* &
              gg21(iq)*tmp1

            tmp1=esi*save1
            psaci(iq)=qizodt(iq)*( 1.-exp(-tmp1*dtb) )
            npsaci(iq)=min( tmp1*niz(iq,k),nizodt(iq,k))
          !
          ! (7) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            tmp1=esw*save1
            psacw(iq)=qlzodt(iq)*( 1.-exp(-tmp1*dtb) )
            npsacw(iq)=min(tmp1*ncz(iq,k),ncz(iq,k))

          ! recalculate the saturatuin temperature
          !
          ! (8) DEPOSITION/SUBLIMATION OF SNOW (Psdep/Pssub): Lin (31)
          !     includes consideration of ventilation effect
          !
            tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
            tmpb=xls*xls*rho(iq,k)*qsiz(iq,k)*diffwv(iq,k)
            tmpc=tmpa*qsiz(iq,k)*diffwv(iq,k)
            abi=4.0*pi*cap_s(iq,k)*(qvoqsiz(iq)-1.0)*tmpc/(tmpa+tmpb)
            tmp1=av_s(iq,k)*sqrho*olambdas(iq,k)**(5.+bv_s(iq,k)+2.*mu_s)/visc(iq,k)

          !---- YLIN, here there is some approximation assuming mu_s =1, so gamma(2)=1, etc.

            tmp2= abi*n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                 vf2s*schmidt(iq,k)**0.33334* &
                 gg22(iq)*sqrt(tmp1) )
            tmp3=odtb*( qvz(iq,k)-qsiz(iq,k) )
            tmp3=min(tmp3,0.)

            if( tmp2 <= 0.) then
              tmp2=max( tmp2,tmp3)
              pssub(iq)=max( tmp2,-qszodt(iq) )
              psdep(iq)=0.0
            else
              psdep(iq)=min( tmp2,tmp3 )
              pssub(iq)=0.0
            end if
            if(qsz(iq,k) >= 0.) then
              !npssub(iq)=pssub(iq)*nsz(iq,k)/qsz(iq,k)
              npsdep(iq)=npsdep(iq)*nsz(iq,k)/qsz(iq,k)
            else
              !npssub(iq)=pssub(iq)/xms
              npsdep(iq)=npsdep(iq)/xms
            end if

            if (qrz(iq,k) > 0.) then !go to 1200
            !
            ! Compute processes (9) and (10) only when qsz > 0.0 and qrz > 0.0
            ! these two terms need to be refined in the future, they should be equal
            !
            ! (9) ACCRETION OF SNOW BY RAIN (Pracs): Lin (27)
            !
              esr=1.0
              tmpa=olambdar(iq,k)**2
              tmpb=olambdas(iq,k)**2
              tmpc=olambdar(iq,k)*olambdas(iq,k)
              tmp1=pi**2*esr*n0_r(iq,k)*n0_s(iq,k)*    &
                        abs( vtr(iq,k)-vts(iq,k) )/rho(iq,k)
            ! tmp1=pi*pi*esr*n0_r(iq,k)*N0_s(k)*            &
            ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5/rho(k)
              tmp2=tmpb**2*olambdar(iq,k)*(5.0*tmpb+2.0*tmpc+0.5*tmpa)
              tmp3=tmp1*rhosnow*tmp2
              pracs(iq)=min( tmp3,qszodt(iq) )
            !
            ! (10) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
            !
              tmp3=tmpa**2*olambdas(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
              tmp4=tmp1*rhowater*tmp3
              psacr(iq)=min( tmp4,qrzodt(iq) )
              tmp1=0.25*pi*esr*n0_r(iq,k)*n0_s(iq,k)*abs( vtr(iq,k)-vts(iq,k) )
              tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2.*tmpb)
              tmp3=tmp1*tmp2
              npsacr(iq)=min( tmp3,nrzodt(iq,k) )
            !
            !
            ! (2) FREEZING OF RAIN TO FORM GRAUPEL  (pgfr): Lin (45), added to PI
            !     positive value
            !     Constant in Bigg freezing Aplume=Ap=0.66 /k
            !     Constant in raindrop freezing equ. Bplume=Bp=100./m/m/m/s
            !

              if (qrz(iq,k) > 1.e-8 ) then
                Bp=100.
                Ap=0.66
                tmp1=olambdar(iq,k)**3
                tmp2=20.*pi**2*Bp*n0_r(iq,k)*rhowater/rho(iq,k)*  &
                  (exp(-Ap*temcc(iq,k))-1.)*tmp1**2*olambdar(iq,k)
                pgfr(iq)=min( tmp2,qrzodt(iq) )
                npgfr(iq)=pi*Bp*n0_r(iq,k)*tmpa**2*(exp(-Ap*temcc(iq,k))-1.0)
              else
                pgfr(iq)=0
                npgfr(iq)=0.
              endif
            end if ! for the go to 1200
          end if   !1200    continue

        else ! if ( tem(iq,k) .lt. 273.15) ..else..

        !
        !***********************************************************************
        !*********        snow production processes for T > 0 C       **********
        !***********************************************************************
        !
          if (qsz(iq,k) > 0.)  then !go to 1400
          !
          ! (1) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
          !
            esw=1.0
            tmp1 = olambdas(iq,k)**(bv_s(iq,k)+tmp_sa(iq,k))
            save1 =aa_s(iq,k)*sqrho*n0_s(iq,k)* &
                   gg21(iq)*     &
                   tmp1

            tmp1=esw*save1
            psacw(iq)=qlzodt(iq)*( 1.0-exp(-tmp1*dtb) )
            npsacw(iq)=tmp1*ncz(iq,k)
          !
          ! (2) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
          !
            esr=1.
            tmpa=olambdar(iq,k)*olambdar(iq,k)
            tmpb=olambdas(iq,k)*olambdas(iq,k)
            tmpc=olambdar(iq,k)*olambdas(iq,k)
            tmp1=pi*pi*esr*n0_r(iq,k)*n0_s(iq,k)*   &
                    abs( vtr(iq,k)-vts(iq,k) )/rho(iq,k)
          ! tmp1=pi*pi*esr*n0_r(iq,k)*N0_s(k)*            &
          ! ( (1.2*vtr(k)-0.95*vts(k))**2+0.08*vtr(k)*vts(k))**0.5/rho(k)
            tmp2=tmpa*tmpa*olambdas(iq,k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
            tmp3=tmp1*rhowater*tmp2
            psacr(iq)=min( tmp3,qrzodt(iq) )

            tmp1=0.25*pi*esr*n0_r(iq,k)*n0_s(iq,k)*abs( vtr(iq,k)-vts(iq,k) )
            tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
            tmp3=tmp1*tmp2
            npsacr(iq)=min( tmp3,nrzodt(iq,k) )
          !
          ! (3) MELTING OF SNOW (Psmlt): Lin (32)
          !     Psmlt is negative value
          !
            delrs=rs0(iq,k)-qvz(iq,k)
            term1=2.0*pi/rho(iq,k)*( xlv*diffwv(iq,k)*rho(iq,k)*delrs- &
                  xka(iq,k)*temcc(iq,k) )
            tmp1= av_s(iq,k)*sqrho*olambdas(iq,k)**(5.+bv_s(iq,k)+2.*mu_s)/visc(iq,k)
            tmp2= n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                  vf2s*schmidt(iq,k)**0.33334* &
                  gg22(iq)*sqrt(tmp1) )
            tmp3=term1*oxlf*tmp2-cwoxlf*temcc(iq,k)*( psacw(iq)+psacr(iq) )
            tmp4=min(0.0,tmp3)
            psmlt(iq)=max( tmp4,-qszodt(iq) )

            if(qsz(iq,k) >= 0.) then
              npsmlt(iq)=psmlt(iq)*nsz(iq,k)/qsz(iq,k)
            else
              npsmlt(iq)=psmlt(iq)/xms
            end if
          !
          ! (4) EVAPORATION OF MELTING SNOW (Psmltevp): HR (A27)
          !     but use Lin et al. coefficience
          !     Psmltevp is a negative value
          !
            tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
            tmpb=xlv*xlv*rho(iq,k)*qswz(iq,k)*diffwv(iq,k)
            tmpc=tmpa*qswz(iq,k)*diffwv(iq,k)
            tmpd=min( 0.0,(qvoqswz(iq)-0.90)*qswz(iq,k)*odtb )

            abr=2.0*pi*(qvoqswz(iq)-0.90)*tmpc/(tmpa+tmpb)
          !
          !**** allow evaporation to occur when RH less than 90%
          !**** here not using 100% because the evaporation cooling
          !**** of temperature is not taking into account yet; hence,
          !**** the qsw value is a little bit larger. This will avoid
          !**** evaporation can generate cloud.
          
            tmp1=av_s(iq,k)*sqrho*olambdas(iq,k)**(5.+bv_s(iq,k)+2.*mu_s)/visc(iq,k)
            tmp2=n0_s(iq,k)*( vf1s*olambdas(iq,k)*olambdas(iq,k)+ &
                 vf2s*schmidt(iq,k)**0.33334* &
                 gg22(iq)*sqrt(tmp1) )
            tmp3=min(0.0,tmp2)
            tmp3=max( tmp3,tmpd )
            psmltevp(iq)=max( tmp3,-qszodt(iq) )
            if(qsz(iq,k) >= 0.) then
              npsmltevp(iq)=psmltevp(iq)*nsz(iq,k)/qsz(iq,k)
            else
              npsmltevp(iq)=psmltevp(iq)/xmr
            end if
          end if !          1400     continue
        end if      !---- end of snow/ice processes   if (tem(iq,k) .lt. 273.15) then
        !---------- end of snow/ice processes below freezing

      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
        
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
        if (qlz(iq,k) > 1.e-6) then
          mu_c(iq,k) = min(15., (1000.E6/ncz(iq,k) + 2.))
          lamc(iq,k) = (ncz(iq,k)*rhowater*pi*gg31(iq)/(6.*qlz(iq,k)*gg32(iq)))**(1./3)
          Dc_liu = (gg33(iq)/gg32(iq))**(1./6.)/lamc(iq,k) !----- R6 in m
          if (Dc_liu > R6c) then
            disp = 1./(mu_c(iq,k)+1.)                      !--- square of relative dispersion
            eta  = (0.75/pi/(1.e-3*rhowater))**2*1.9e11*((1.+3.*disp)*(1.+4.*disp)*&
                   (1.+5.*disp)/(1.+disp)/(1.+2.*disp))
            praut(iq) = eta*(1.e-3*rho(iq,k)*qlz(iq,k))**3/(1.e-6*ncz(iq,k))  !--- g cm-3 s-1
            praut(iq) = praut(iq)/(1.e-3*rho(iq,k))                           !--- kg kg-1 s-1
            npraut_r(iq) = praut(iq)/xmr                                                !--- kg kg-1 s-1
            npraut(iq) = praut(iq)/qlz(iq,k)*ncz(iq,k)                                        !--- kg kg-1 s-1
            npraut(iq) = praut(iq)/xmr                                                  !--- kg kg-1 s-1
          else
            praut(iq) = 0.0
            npraut(iq) = 0.0
            npraut_r(iq) = 0.0
          endif
        else
          praut(iq) = 0.0
          npraut(iq) = 0.0
          npraut_r(iq) = 0.0
        endif
        ! if (qlz(k) .gt. 1e-6) then
        ! praut(k)=1350.*qlz(k)**2.47*  &
        ! (ncz(k)/1.e6*rho(k))**(-1.79)
        ! npraut_r(k) = praut(k)/xmr
        ! npraut(k) = praut(k)/(qlz(k)/ncz(k))
        ! npraut(K) = MIN(npraut(k),nczodt(k))
        ! npraut_r(K) = MIN(npraut_r(k),npraut(k))
        ! endif

      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000

        !
        ! (2) ACCRETION OF CLOUD WATER BY RAIN (Pracw): Lin (51)
        !
        erw=1.
        tmp1=pio4*erw*n0_r(iq,k)*av_r*sqrho* &
             gambp3*olambdar(iq,k)**bp3 ! son
        pracw(iq)=qlzodt(iq)*( 1.-exp(-tmp1*dtb) )
        npracw(iq)=tmp1*ncz(iq,k)

      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
        
        !
        ! (3) EVAPORATION OF RAIN (Prevp): Lin (52)
        !     Prevp is negative value
        !
        !     Sw=qvoqsw : saturation ratio
        !
        tmpa=rvapor*xka(iq,k)*tem(iq,k)*tem(iq,k)
        tmpb=xlv*xlv*rho(iq,k)*qswz(iq,k)*diffwv(iq,k)
        tmpc=tmpa*qswz(iq,k)*diffwv(iq,k)
        tmpd=min(0.0,(qvoqswz(iq)-0.99)*qswz(iq,k)*odtb)

        abr=2.0*pi*(qvoqswz(iq)-0.99)*tmpc/(tmpa+tmpb)
        tmp1=av_r*sqrho*olambdar(iq,k)**bp5/visc(iq,k) !son
        tmp2=abr*n0_r(iq,k)*( vf1r*olambdar(iq,k)*olambdar(iq,k)+  &
             vf2r*schmidt(iq,k)**0.33334*gambp5o2*sqrt(tmp1) )
        tmp3=min( 0.0,tmp2 )
        tmp3=max( tmp3,tmpd )
        prevp(iq)=max( tmp3,-qrzodt(iq) )
        if (qrz(iq,k).gt.0.) then
          nprevp(iq)=prevp(iq)*nrz(iq,k)/qrz(iq,k)
        else
          nprevp(iq)=prevp(iq)*xmr
        end if

      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
    
        !
        !**********************************************************************
        !*****     combine all processes together and avoid negative      *****
        !*****     water substances
        !***********************************************************************
        !
        if ( temcc(iq,k) < 0.) then
        !
        !  combined water vapor depletions
        !
          tmp=psdep(iq) + midep(iq)
          if ( tmp > qvzodt(iq) ) then
            factor=qvzodt(iq)/tmp
            psdep(iq)=psdep(iq)*factor
            midep(iq)=midep(iq)*factor
          end if
        !
        !  combined cloud water depletions
        !
          tmp=praut(iq)+psacw(iq)+psfw(iq)+pracw(iq)
          if ( tmp > qlzodt(iq) ) then
            factor=qlzodt(iq)/tmp
            praut(iq)=praut(iq)*factor
            psacw(iq)=psacw(iq)*factor
            psfw(iq)=psfw(iq)*factor
            pracw(iq)=pracw(iq)*factor
          end if
        !
        !  combined cloud ice depletions
        !
          tmp=psaut(iq)+psaci(iq)+praci(iq)+psfi(iq)
          if (tmp > qizodt(iq) ) then
            factor=qizodt(iq)/tmp
            psaut(iq)=psaut(iq)*factor
            psaci(iq)=psaci(iq)*factor
            praci(iq)=praci(iq)*factor
            psfi(iq)=psfi(iq)*factor
          endif

        !
        !  combined all rain processes
        !
          tmp_r=piacr(iq)+psacr(iq)-prevp(iq)-  & 
                praut(iq)-pracw(iq)+pgfr(iq)
          if (tmp_r > qrzodt(iq) ) then
            factor=qrzodt(iq)/tmp_r
            piacr(iq)=piacr(iq)*factor
            psacr(iq)=psacr(iq)*factor
            prevp(iq)=prevp(iq)*factor
            pgfr(iq)=pgfr(iq)*factor
          endif
        !
        !   combined all snow processes
        !
          tmp_s=-pssub(iq)-(psaut(iq)+psaci(iq)+ &
                 psacw(iq)+psfw(iq)+pgfr(iq)+ &
                 psfi(iq)+praci(iq)+piacr(iq)+ &
                 psdep(iq)+psacr(iq)-pracs(iq))
          if ( tmp_s > qszodt(iq) ) then
            factor=qszodt(iq)/tmp_s
            pssub(iq)=pssub(iq)*factor
            pracs(iq)=pracs(iq)*factor
          endif

        !
        ! calculate new water substances, thetae, tem, and qvsbar
        !

          pvapor(iq)=-pssub(iq)-psdep(iq)-prevp(iq)-midep(iq)
          qvz(iq,k)=max( qvmin,qvz(iq,k)+dtb*pvapor(iq) )
          pclw(iq)=-praut(iq)-pracw(iq)-psacw(iq)-psfw(iq)
          qlz(iq,k)=max( 0.0,qlz(iq,k)+dtb*pclw(iq) )
          pcli(iq)=-psaut(iq)-psfi(iq)-psaci(iq)- & 
                    praci(iq)+midep(iq)
          qiz(iq,k)=max( 0.0,qiz(iq,k)+dtb*pcli(iq) )
          tmp_r=piacr(iq)+psacr(iq)-prevp(iq)-praut(iq)- &
                    pracw(iq)+pgfr(iq)-pracs(iq)
          prain(iq)=-tmp_r
          qrz(iq,k)=max( 0.0,qrz(iq,k)+dtb*prain(iq) )
          tmp_s=-pssub(iq)-(psaut(iq)+psaci(iq)+ &
                 psacw(iq)+psfw(iq)+pgfr(iq)+  &
                 psfi(iq)+praci(iq)+piacr(iq)+  &
                 psdep(iq)+psacr(iq)-pracs(iq))
          psnow(iq)=-tmp_s
          qsz(iq,k)=max( 0.0,qsz(iq,k)+dtb*psnow(iq) )
          !qschg=qschg+psnow(iq)
          qschg(iq)=psnow(iq)

          tmp=ocp/tothz(iq,k)*xLf*qschg(iq)
          theiz(iq,k)=theiz(iq,k)+dtb*tmp
          ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
          ! tem(k)=thz(k)*tothz(k)
          ! temcc(k)=tem(k)-273.15
          !==================update temperature=================================================
          temcc(iq,k)=tem(iq,k)-273.15
          lvap = xlv + (2106.0 - 4218.0)*temcc(iq,k)  !Enthalpy of vaporization
          tmp1=(pssub(iq)+psdep(iq))*xls*ocp + prevp(iq)*lvap*ocp+  &
               (psfw(iq)+pgfr(iq)+psacr(iq)-pracs(iq))*xlfocp
          !bug fixed 20191126
          tem(iq,k)=tem(iq,k)+tmp1*dtb
          temcc(iq,k)=tem(iq,k)-273.15
          thz(iq,k)=tem(iq,k)/tothz(iq,k)
          !===================================================================
          if ( temcc(iq,k) < -40. ) qswz(iq,k)=qsiz(iq,k)
          qlpqi=qlz(iq,k)+qiz(iq,k)
          if ( qlpqi == 0. ) then
            qvsbar(iq)=qsiz(iq,k)
          else
            qvsbar(iq)=(qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k))/qlpqi
          endif
          tmp1=-npraut(iq)-npracw(iq)-npsacw(iq)
          ncz(iq,k)=max( 0.0,ncz(iq,k)+dtb*tmp1 )
          tmp1=-npsaut(iq)-npsaci(iq)-npraci(iq)+nidep(iq)
          niz(iq,k)=max( 0.0,niz(iq,k)+dtb*tmp1 )
          tmp1=npiacr(iq)+npsacr(iq)-nprevp(iq)-npraut_r(iq)+npgfr(iq)
          nrz(iq,k)=max( 0.0,nrz(iq,k)-dtb*tmp1 )
          tmp1=-(npsaut(iq)+npgfr(iq)+  &
                 npraci(iq)+npiacr(iq)+  &
                 npsdep(iq)+npsacr(iq))
          nsz(iq,k)=max( 0.0,nsz(iq,k)-dtb*tmp1 )
            
        else ! if ( temcc(iq,k) .lt. 0.0) ..else..
          !
          !  combined cloud water depletions
          !
          tmp=praut(iq)+psacw(iq)+pracw(iq)
          if ( tmp > qlzodt(iq) ) then
            factor=qlzodt(iq)/tmp
            praut(iq)=praut(iq)*factor
            psacw(iq)=psacw(iq)*factor
            pracw(iq)=pracw(iq)*factor
          end if
          !
          !  combined all snow processes
          !
          tmp_s=-(psmlt(iq)+psmltevp(iq))
          if (tmp_s > qszodt(iq) ) then
            factor=qszodt(iq)/tmp_s
            psmlt(iq)=psmlt(iq)*factor
            psmltevp(iq)=psmltevp(iq)*factor
          endif
          !
          !  combined all rain processes
          !
          tmp_r=-prevp(iq)-(praut(iq)+pracw(iq)+psacw(iq)-psmlt(iq))
          if (tmp_r > qrzodt(iq) ) then
            factor=qrzodt(iq)/tmp_r
            prevp(iq)=prevp(iq)*factor
          endif
          !
          !  calculate new water substances and thetae
          !
          pvapor(iq)=-psmltevp(iq)-prevp(iq)
          qvz(iq,k)=max( qvmin,qvz(iq,k)+dtb*pvapor(iq))
          pclw(iq)=-praut(iq)-pracw(iq)-psacw(iq)
          qlz(iq,k)=max( 0.0,qlz(iq,k)+dtb*pclw(iq) )
          pcli(iq)=0.0
          qiz(iq,k)=max( 0.0,qiz(iq,k)+dtb*pcli(iq) )
          tmp_r=-prevp(iq)-(praut(iq)+pracw(iq)+psacw(iq)-psmlt(iq))
          prain(iq)=-tmp_r
          !tmpqrz=qrz(iq,k)
          qrz(iq,k)=max( 0.0,qrz(iq,k)+dtb*prain(iq) )
          tmp_s=-(psmlt(iq)+psmltevp(iq))
          psnow(iq)=-tmp_s
          qsz(iq,k)=max( 0.0,qsz(iq,k)+dtb*psnow(iq) )
          qschg(iq)=psnow(iq)

          tmp=ocp/tothz(iq,k)*xLf*qschg(iq)
          theiz(iq,k)=theiz(iq,k)+dtb*tmp
          ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
          ! tem(k)=thz(k)*tothz(k)
          ! temcc(k)=tem(k)-273.15
          !==================update tmperature=================================================
          temcc(iq,k)=tem(iq,k)-273.15
          lvap = xlv + (2106.0 - 4218.0)*temcc(iq,k)  !Enthalpy of vaporization
          tmp1=psmltevp(iq)*xls*ocp + prevp(iq)*lvap*ocp+  &
               psmlt(iq)*xlfocp
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
          qswz(iq,k)=ep2*es/max(prez(iq,k)-es,0.1)
          qsiz(iq,k)=qswz(iq,k)
          qvsbar(iq)=qswz(iq,k)
          ! tmp1=-(npraut(iq)+npsacw(iq)+npracw(iq))
          ncz(iq,k)=max( 0.0,ncz(iq,k)+dtb*tmp1)
          tmp1=-nprevp(iq)-(npraut_r(iq)-npsmlt(iq))
          ! tmp1=-nprevp(k)-(nprautr(k)+npracwr(k)+npsacw(k)-npsmltr(k))
          nrz(iq,k)=max(0.0,nrz(iq,k)-dtb*tmp1)
          tmp1=-(npsmlt(iq)+npsmltevp(iq))
          nsz(iq,k)=max( 0.0,nsz(iq,k)-dtb*tmp1 )
        end if    !T seperate for source and sink terms

      end if
    end do
    
    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000        

        !rain
        if (qrz(iq,k) > 1.e-8) then
          xlambdar(iq,k)=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
          if (xlambdar(iq,k).lt.lamminr) then
            xlambdar(iq,k) = lamminr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,K)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,K)
          else if (xlambdar(iq,k).gt.lammaxr) then
            xlambdar(iq,k) = lammaxr
            n0_r(iq,k) = xlambdar(iq,k)**4*qrz(iq,K)/(pi*rhowater)
            nrz(iq,k) = n0_r(iq,k)/xlambdar(iq,K)
          end if
        end if

      end if
    end do
    
    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
        
        !snow
        if (qsz(iq,k) .gt. 1.0e-8) then
          xlambdas(iq,k)=(am_s(iq,k)*gam_ss(iq,k)*     &
          nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
          if (xlambdas(iq,k).lt.lammins) then
            xlambdas(iq,k)= lamminS
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1.)*      &
                           qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          else if (xlambdas(iq,k).gt.lammaxs) then
            xlambdas(iq,k) = lammaxs
            n0_s(iq,K) = xlambdas(iq,k)**(bm_s(iq,k)+1.)*      & 
                           qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
            nsz(iq,K) = n0_s(iq,K)/xlambdas(iq,k)
          end if
        end if

      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000
        
        !cloud ice
        if (qiz(iq,k) >= 1.e-8) then
          lami(iq,k) = max((gam13*500.*pi/6.)*niz(iq,k)/qiz(iq,k),1.e-20)**(1./3) !fixed zdc
          if (lami(iq,k).lt.lammini) then
            lami(iq,k)= lammini
            n0_i(iq) = lami(iq,k)**4./gam13*500.*pi/6.
            niz(iq,K) = n0_i(iq)/lami(iq,k)
          else if (lami(iq,k).gt.lammaxi) then
            lami(iq,k) = lammaxi
            n0_i(iq) = lami(iq,k)**4./gam13*500.*pi/6.
            niz(iq,K) = n0_i(iq)/lami(iq,k)
          end if
        end if
        
      end if
    end do

    do iq = 1,imax 
      if( mask(iq) ) then !go to 2000

        !cloud water zdc 20220208
        if (qlz(iq,k) >= 1.e-8) then
          lamc(iq,k) = (ncz(iq,k)*rhowater*pi*gg31(iq)/(6.*qlz(iq,k)*gg32(iq)))**(1./3)
          if (lamc(iq,k).lt.lammini) then
            lamc(iq,k)= lammini
            tmp1 = lamc(iq,k)**(mu_c(iq,k)+4.)
            tmp2 = 6.*qlz(iq,k)/(pi*rhowater*gg31(iq))
            n0_c(iq)= tmp1*tmp2
            ncz(iq,k) = n0_c(iq)/lamc(iq,k)
          else if (lamc(iq,k).gt.lammaxi) then
            lamc(iq,k)= lammaxi
            tmp1 = lamc(iq,k)**(mu_c(iq,k)+4.)
            tmp2 = 6.*qlz(iq,k)/(pi*rhowater*gg31(iq))
            n0_c(iq)= tmp1*tmp2
            ncz(iq,k) = n0_c(iq)/lamc(iq,k)
          end if
        end if
        
      end if
    end do

    do iq = 1,imax
      if( mask(iq) ) then !go to 2000    
        
        nimlt = 0.
        nihom = 0.
    
        !
        !***********************************************************************
        !**********              saturation adjustment                **********
        !***********************************************************************
        !
        !    allow supersaturation exits linearly from 0% at 500 mb to 50%
        !    above 300 mb
        !    5.0e-5=1.0/(500mb-300mb)
        !
        rsat=1.
        if( qvz(iq,k)+qlz(iq,k)+qiz(iq,k) < rsat*qvsbar(iq) ) then ! goto 1800

        !
        !   unsaturated
        !
          qvz(iq,k)=qvz(iq,k)+qlz(iq,k)+qiz(iq,k)
          qlz(iq,k)=0.
          qiz(iq,k)=0.

          thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)

          tem(iq,k)=thz(iq,k)*tothz(iq,k)
          temcc(iq,k)=tem(iq,k)-273.15

        else
        !
        !   saturated
        !
          pladj(iq)=qlz(iq,k)
          piadj(iq)=qiz(iq,k)
        !

          CALL satadj(qvz(iq,k), qlz(iq,k), qiz(iq,k), prez(iq,k), &
                      theiz(iq,k), thz(iq,k), tothz(iq,k), &
                      xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

          pladj(iq)=odtb*(qlz(iq,k)-pladj(iq))
          piadj(iq)=odtb*(qiz(iq,k)-piadj(iq))
          pclw(iq)=pclw(iq)+pladj(iq)
          pcli(iq)=pcli(iq)+piadj(iq)
          pvapor(iq)=pvapor(iq)-( pladj(iq)+piadj(iq) )
          thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
          tem(iq,k)=thz(iq,k)*tothz(iq,k)
          temcc(iq,k)=tem(iq,k)-273.15

          es=1000.*svp1*exp( svp2*temcc(iq,k)/(tem(iq,k)-svp3) )
          qswz(iq,k)=ep2*es/max(prez(iq,k)-es,0.1)
          if (tem(iq,k) < 273.15 ) then
            es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
            qsiz(iq,k)=ep2*es/max(prez(iq,k)-es,0.1)
            if (temcc(iq,k) < -40.) qswz(iq,k)=qsiz(iq,k)
          else
            qsiz(iq,k)=qswz(iq,k)
          endif
          qlpqi=qlz(iq,k)+qiz(iq,k)
          if ( qlpqi .eq. 0. ) then
            qvsbar(iq)=qsiz(iq,k)
          else
            qvsbar(iq)=( qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k) )/qlpqi
          endif

        !
        !***********************************************************************
        !*****     melting and freezing of cloud ice and cloud water       *****
        !***********************************************************************
          qlpqi=qlz(iq,k)+qiz(iq,k)
          if( qlpqi .gt. 0. ) then !go to 1800
          !
          !
          ! (1)  HOMOGENEOUS NUCLEATION WHEN T< -40 C (Pihom)
          !
            if(temcc(iq,k) .lt. -40.) then
              pihom(iq)=qlz(iq,k)*odtb
              nihom=ncz(iq,k)*odtb
            end if
          !
          ! (2)  MELTING OF ICE CRYSTAL WHEN T> 0 C (Pimlt)
          !
            if(temcc(iq,k) .gt. 0.) then
              pimlt(iq)=qiz(iq,k)*odtb
              nimlt=niz(iq,k)*odtb
            end if
          !
          ! (3) PRODUCTION OF CLOUD ICE BY BERGERON PROCESS (Pidw): Hsie (p957)
          !     this process only considered when -31 C < T < 0 C
          !
            if(temcc(iq,k) < 0. .and. temcc(iq,k) > -31.) then
          !!
          !!   parama1 and parama2 functions must be user supplied
          !!
              ! aa1=parama1(temcc(iq,k))
              ! aa2=parama2(temcc(iq,k))
              i1=int(-temcc(iq,k))+1
              i1p1=i1+1
              ratio=-(temcc(iq,k))-real(i1-1)
              aa1=a1(i1)+ratio*( a1(i1p1)-a1(i1) )
              aa2=a2(i1)+ratio*( a2(i1p1)-a2(i1) )
              
              !! change unit from cgs to mks
              aa1=aa1*0.001**(1.0-aa2)
              xnin=xni0*exp(-bni*temcc(iq,k))
              pidw(iq)=xnin/rho(iq,k)*(aa1*xmnin**aa2)
            end if

            pcli(iq)=pcli(iq)+pihom(iq)-pimlt(iq)+pidw(iq)
            pclw(iq)=pclw(iq)-pihom(iq)+pimlt(iq)-pidw(iq)
            qlz(iq,k)=max( 0.,qlz(iq,k)+dtb*(-pihom(iq)+pimlt(iq)-pidw(iq)) )
            qiz(iq,k)=max( 0.,qiz(iq,k)+dtb*(pihom(iq)-pimlt(iq)+pidw(iq)) )

            ncz(iq,k)=max( 0.,ncz(iq,k)+dtb*(-nihom+nimlt) )
            niz(iq,k)=max( 0.,niz(iq,k)+dtb*( nihom-nimlt) )

            CALL satadj(qvz(iq,k), qlz(iq,k), qiz(iq,k), prez(iq,k), &
                    theiz(iq,k), thz(iq,k), tothz(iq,k), &
                    xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

            thz(iq,k)=theiz(iq,k)-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
            tem(iq,k)=thz(iq,k)*tothz(iq,k)
            temcc(iq,k)=tem(iq,k)-273.15
            es=1000.*svp1*exp( svp2*temcc(iq,k)/(tem(iq,k)-svp3) )
            qswz(iq,k)=ep2*es/max(prez(iq,k)-es,0.1)

            if (tem(iq,k) < 273.15 ) then
              es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
              qsiz(iq,k)=ep2*es/max(prez(iq,k)-es,0.1)
              if (temcc(iq,k) < -40.0) qswz(iq,k)=qsiz(iq,k)
            else
              qsiz(iq,k)=qswz(iq,k)
            endif
            qlpqi=qlz(iq,k)+qiz(iq,k)

            if ( qlpqi == 0. ) then
              qvsbar(iq)=qsiz(iq,k)
            else
              qvsbar(iq)=( qiz(iq,k)*qsiz(iq,k)+qlz(iq,k)*qswz(iq,k) )/qlpqi
            endif
          end if ! 1800  continue
        end if ! 1800  continue

      end if
    end do     ! iq

    !
    !***********************************************************************
    !**********    integrate the productions of rain and snow     **********
    !***********************************************************************
    !    
    
    where ( qvz(1:imax,k) < qvmin ) 
      qvz(1:imax,k)=max( qvmin,qvz(1:imax,k)+qlz(1:imax,k)+qiz(1:imax,k) )
      qlz(1:imax,k)=0.
      qiz(1:imax,k)=0.
      ncz(1:imax,k)=0.
      niz(1:imax,k)=0.
    end where
    niz(1:imax,k) = min( niz(1:imax,k), 0.3E6/rho(1:imax,k) )
    ncz(1:imax,k) = min( ncz(1:imax,k), 250000.E6/rho(1:imax,k) )
    ncz(1:imax,k) = max( ncz(1:imax,k), 0.01E6/rho(1:imax,k) )


  ! CALCULATE EFFECTIVE RADIUS zdc 20220208
    where (qiz(1:imax,k) > 1.e-8 .and. lami(1:imax,k)>0. ) 
      EFFI1D(1:imax,k) = 3./LAMI(1:imax,k)/2.
    elsewhere
      EFFI1D(1:imax,k) = 25.E-6
    end where

    where (qsz(1:imax,k) > 1.0e-8) 
      EFFS1D(1:imax,k) = 3./xlambdas(1:imax,k)/2.
    elsewhere
      EFFS1D(1:imax,k) = 25.E-6
    end where

    where (qrz(1:imax,k) > 1.e-8) 
      EFFR1D(1:imax,k) = 3./xlambdar(1:imax,k)/2.
    elsewhere
      EFFR1D(1:imax,k) = 25.E-6
    end where

    where (qlz(1:imax,k)>1.e-8 .and. lamc(1:imax,k)>0.)
      EFFC1D(1:imax,k) = GAMMA(mu_c(1:imax,k)+4.)/GAMMA(mu_c(1:imax,k)+3.)/LAMC(1:imax,k)/2.
    elsewhere
      EFFC1D(1:imax,k) = 25.E-6
    end where


#ifdef diagnostic
  ! save all process rate for understanding cloud microphysics
    zpsaut(1:imax,k)   = psaut(1:imax)    ! ice crystal aggregation to snow
    zpsfw(1:imax,k)    = psfw(1:imax)     ! BERGERON process to transfer cloud water to snow
    zpsfi(1:imax,k)    = psfi(1:imax)     ! BERGERON process to transfer cloud ice to snow
    zpraci(1:imax,k)   = praci(1:imax)    ! cloud ice accretion by rain
    zpiacr(1:imax,k)   = piacr(1:imax)    ! rain accretion by cloud ice
    zpsaci(1:imax,k)   = psaci(1:imax)    ! ice crystal accretion by snow
    zpsacw(1:imax,k)   = psacw(1:imax)    ! accretion of cloud water by snow
    zpsdep(1:imax,k)   = psdep(1:imax)    ! deposition of snow
    zpssub(1:imax,k)   = pssub(1:imax)    ! sublimation of snow (T<0)
    zpracs(1:imax,k)   = pracs(1:imax)    ! accretion of snow by rain
    zpsacr(1:imax,k)   = psacr(1:imax)    ! accretion of rain by snow
    zpsmlt(1:imax,k)   = psmlt(1:imax)    ! melting of snow
    zpsmltevp(1:imax,k)= psmltevp(1:imax) ! evaporation of melting snow (T>0)
    zprain(1:imax,k)   = prain(1:imax)    ! sum all process for rain
    zpraut(1:imax,k)   = praut(1:imax)    ! autoconversion of rain
    zpracw(1:imax,k)   = pracw(1:imax)    ! accretion of cloud water by rain
    zprevp(1:imax,k)   = prevp(1:imax)    ! evaporation of rain
    zpgfr(1:imax,k)    = pgfr(1:imax)     ! feezing of rain to form graupel (added to PI)
    zpvapor(1:imax,k)  = pvapor(1:imax)   ! sum all process for water vapor to determine qvz
    zpclw(1:imax,k)    = pclw(1:imax)     ! sum all process for cloud liquid to determine qlz
    zpladj(1:imax,k)   = pladj(1:imax)    ! saturation adjustment for ql
    zpcli(1:imax,k)    = pcli(1:imax)     ! sum all process for cloud ice to determine qiz
    zpimlt(1:imax,k)   = pimlt(1:imax)    ! melting of ice crystal >0.
    zpihom(1:imax,k)   = pihom(1:imax)    ! homogeneous nucleation <-40
    zpidw(1:imax,k)    = pidw(1:imax)     ! production of cloud ice by BERGERON process
    zpiadj(1:imax,k)   = piadj(1:imax)    ! saturation adjustment for qi
    zpsnow(1:imax,k)   = psnow(1:imax)
    zqschg(1:imax,k)   = qschg(1:imax)
#endif


  ! save process rate for aerisol scheme
    fluxm(1:imax,k) = -1.*psmlt(1:imax)*dzw(1:imax,k)*     &!
                        rho(1:imax,k)                       ! - ice melting flux in layer k (kg/m2/s)
    fluxf(1:imax,k) = pgfr(1:imax)*dzw(1:imax,k)*          &!
                        rho(1:imax,k)                       ! - liquid freezing flux in layer k (kg/m2/s)
    fevap(1:imax,k) = -1.*prevp(1:imax)*dzw(1:imax,k)*     &!
                        rho(1:imax,k)                       ! - evaporation of rainfall flux (kg/m2/s)
    fsubl(1:imax,k) = -1.*pssub(1:imax)*dzw(1:imax,k)*     &!
                        rho(1:imax,k)                       ! - sublimation of snow, ice and graupel flux (kg/m2/s)
    fauto(1:imax,k) = praut(1:imax)*dzw(1:imax,k)*         &!
                        rho(1:imax,k)                       ! - autoconversion flux for rainfall (kg/m2/s)
    fcoll(1:imax,k) = pracw(1:imax)*dzw(1:imax,k)*         &!
                        rho(1:imax,k)                       ! - collection of cloud liquid water by rain (kg/m2/s)
    faccr(1:imax,k) = psacw(1:imax)*dzw(1:imax,k)*         &!
                        rho(1:imax,k) + psfw(1:imax)*      &! - accretion of cloud liq water by snow,ice and graupel (kg/m2/s)
                        dzw(1:imax,k)*rho(1:imax,k)
    vi(1:imax,k)    = vtiold(1:imax,k)
    !vs(1:imax,k)    = vtsold(1:imax,k)

  ENDDO        ! k

END SUBROUTINE clphy1d_ylin



!---------------------------------------------------------------------
!                         SATURATED ADJUSTMENT
!---------------------------------------------------------------------
pure SUBROUTINE satadj(qvz, qlz, qiz, prez, theiz, thz, tothz,      &
                  xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

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


  REAL,     INTENT(INOUT) :: qvz, qlz, qiz, thz
  REAL,     INTENT(IN   ) :: prez, theiz, tothz
  REAL,     INTENT(IN   ) :: xLvocp, xLfocp, episp0k
  REAL,     INTENT(IN   ) :: EP2, SVP1, SVP2, SVP3, SVPT0

  ! LOCAL VARS

  INTEGER                            :: n
  REAL                               :: tem, temcc, qsiz,       &
                                        qswz, qvsbar
  REAL :: qsat, qlpqi, ratql, t0, t1, tmp1, ratqi, tsat, absft,    &
             denom1, denom2, dqvsbar, ftsat, dftsat, qpz,es             
!---------------------------------------------------------------------
  
  thz=theiz-(xLvocp*qvz-xLfocp*qiz)/tothz
  tem=tothz*thz
  if (tem .gt. 273.15) then
  ! qsat=episp0k/prez*  &
  ! exp( svp2*(tem-273.15)/(tem-svp3) )
    es=1000.*svp1*exp( svp2*(tem-svpt0)/(tem-svp3) )
    qsat=ep2*es/max(prez-es,0.1)
  else
    qsat=episp0k/prez*  &
    exp( 21.8745584*(tem-273.15)/(tem-7.66) )
  end if
  qpz=qvz+qlz+qiz
  if (qpz < qsat) then
    qvz=qpz
    qiz=0.0
    qlz=0.0
    !  go to 400
    ! end if
  else     ! this else to remove the go to above
    qlpqi=qlz+qiz
    if( qlpqi .ge. 1.0e-5) then
      ratql=qlz/qlpqi
      ratqi=qiz/qlpqi
    else
      t0=273.15
      ! t1=233.15
      t1=248.15
      tmp1=( t0-tem )/(t0-t1)
      tmp1=min(1.0,tmp1)
      tmp1=max(0.0,tmp1)
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

    tsat=tem
    absft=1.
    !
    do n = 1,20
      denom1=1.0/(tsat-svp3)
      denom2=1.0/(tsat-7.66)
      ! qswz(k)=episp0k/prez(k)*  &
      ! exp( svp2*denom1*(tsat-273.15) )
      es=1000.*svp1*exp( svp2*denom1*(tsat-svpt0) )
      qswz=ep2*es/max(prez-es,0.1)
      if (tem < 233.15 ) then
        es=1000.*svp1*exp( 21.8745584*denom2*(tsat-273.15) )
        qsiz=ep2*es/max(prez-es,0.1)
        qswz=qsiz
      else if (tem < 273.15) then
        es=1000.*svp1*exp( 21.8745584*denom2*(tsat-273.15) )
        qsiz=ep2*es/max(prez-es,0.1)
      else
        qsiz=qswz
      endif
      qvsbar = ratql*qswz + ratqi*qsiz
      ! if( absft .lt. 0.01 .and. n .gt. 3 ) go to 300
      if( absft >= 0.01 ) then !go to 300
        dqvsbar=ratql*qswz*svp2*243.5*denom1*denom1+  &
                ratqi*qsiz*21.8745584*265.5*denom2*denom2
        ftsat=tsat+(xlvocp+ratqi*xlfocp)*qvsbar-  &
              tothz*theiz-xlfocp*ratqi*(qvz+qlz+qiz)
        dftsat=1.+(xlvocp+ratqi*xlfocp)*dqvsbar
        tsat=tsat-ftsat/dftsat
        absft=abs(ftsat)
      else   ! original
        exit ! original
      end if !300   continue
    end do !200   continue
    !9020  format(1x,'point can not converge, absft,n=',e12.5,i5)
    !300   continue

    if( qpz > qvsbar ) then
      qvz=qvsbar
      qiz=ratqi*( qpz-qvz )
      qlz=ratql*( qpz-qvz )
    else
      qvz=qpz
      qiz=0.
      qlz=0.
    end if
  end if !400  continue

END SUBROUTINE satadj

!----------------------------------------------------------------
PURE FUNCTION ggamma(X) result(ans)
!----------------------------------------------------------------
  IMPLICIT NONE
  !----------------------------------------------------------------
  REAL, INTENT(IN   ) :: x
  INTEGER             ::diff, j
  REAL                ::PF, G1TO2 ,TEMP
  real                :: ans
  
  TEMP=X
  diff = max(int(temp-2.), 0)
  if ( temp-real(diff) > 2. ) diff = diff + 1

  ! Original method
  PF=1.
  do J=1,diff
    TEMP=TEMP-1.
    PF=PF*TEMP
  end do

  ! Alternative method
  !temp = temp - real(diff)
  !pf = gamma( x ) / gamma( temp )
  
  TEMP=TEMP - 1.
  G1TO2=1. + B(1)*TEMP + B(2)*TEMP**2 + B(3)*TEMP**3 + B(4)*TEMP**4 &
           + B(5)*TEMP**5 + B(6)*TEMP**6 + B(7)*TEMP**7 + B(8)*TEMP**8
  ans=PF*G1TO2
  
END FUNCTION ggamma

!----------------------------------------------------------------

END MODULE module_mp_sbu_ylin
