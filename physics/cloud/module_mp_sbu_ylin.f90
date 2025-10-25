!--- The code is based on Lin and Colle (A New Bulk Microphysical Scheme 
!             that Includes Riming Intensity and Temperature Dependent Ice Characteristics, 2011, MWR)
!             and Lin et al. (Parameterization of riming intensity and its impact on ice fall speed using ARM data, 2011, MWR)
!--- NOTE: 1) Prognose variables are: qi,PI(precipitating ice, qs, which includes snow, partially rimed snow and graupel),qw,qr
!---       2) Sedimentation flux is based on Prudue Lin scheme 
!---       2) PI has varying properties depending on riming intensity (Ri, diagnosed currently following Lin et al. (2011, MWR)
!--           and T 
!---       3) Autoconverion is based on Liu and Daum (2004)         
!---       4) PI size distribution assuming Gamma distribution, but mu_s=0 (Exponential) currently
!---       5) No density dependent fall speed since the V-D is derived using Best number approach, which already includes density
!---          effect 
!---       6) Future work will include radar equivalent reflectivity using the new PI property (A-D, M-D, N(D)). If you use RIP
!---          for reflectivity computation, please note that snow is (1-Ri)*qs and graupel is Ri*qs. Otherwise, reflectivity will
!---          be underestimated.
!---       7) The Liu and Daum autoconverion is quite sensitive on Nt_c. For mixed-phase cloud and marine environment, Nt_c of 10
!---          or 20 is suggested. default value is 10E.6. Change accordingly for your use.

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
  
  real, parameter :: qvmin=1.e-20
  real, parameter :: sqrho=1.

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

subroutine clphy1d_ylin(dt_in, imax,                        &
                      qvz, qlz, qrz, qiz, qsz,              &
                      thz, tothz, rho,                      &
                      prez, zz, dzw,                        &
                      EFFC1D, EFFI1D, EFFS1D, EFFR1D,       & !zdc 20220208
                      pptrain, pptsnow,pptice,              &
                      kts, kte, Ri,                         &
                      ncz, nrz, niz, nsz,                   &
                      fluxr, fluxi, fluxs, fluxm,           &
                      fluxf, fevap, fsubl, fauto, fcoll,    &
                      faccr, vi,                            &
#ifndef GPU    
                      zpsnow,zpsaut,zpsfw,zpsfi,zpraci,     & !process rate 
                      zpiacr,zpsaci,zpsacw,zpsdep,          &
                      zpssub,zpracs,zpsacr,zpsmlt,          &
                      zpsmltevp,zprain,zpraut,zpracw,       &
                      zprevp,zpgfr,zpvapor,zpclw,           &
                      zpladj,zpcli,zpimlt,zpihom,           &
                      zpidw,zpiadj,zqschg,                  &
#endif                          
                      zdrop,lin_aerosolmode,lin_adv,        &
                      njumps)

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
  integer,                          intent(in)    :: kts, kte
  integer,                          intent(in)    :: imax, njumps
  integer,                          intent(in)    :: lin_aerosolmode, lin_adv
  real,                             intent(in)    :: dt_in
  real, dimension(1:imax,kts:kte), intent(in)     :: zdrop
  real, dimension(1:imax,kts:kte), intent(inout)  :: Ri
  real, dimension(1:imax,kts:kte), intent(in)     :: tothz,rho,prez,dzw
  real, dimension(1:imax,0:kte),   intent(in)     :: zz
  real, dimension(1:imax,kts:kte), intent(out)    :: EFFC1D,EFFI1D,                     &
                                                     EFFS1D,EFFR1D
  real, dimension(1:imax,kts:kte), intent(out)    :: fluxr,fluxi,fluxs,                 &
                                                     fluxm,fluxf,fevap,fsubl,           &
                                                     fauto,fcoll,faccr
  real, dimension(1:imax,kts:kte), intent(out)    :: vi
#ifndef GPU  
  real, dimension(1:imax,kts:kte), intent(out)    :: zpsnow,zpsaut,zpsfw,               &
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
  real, dimension(1:imax),         intent(inout)  :: pptrain, pptsnow, pptice
  real, dimension(1:imax,kts:kte), intent(inout)  :: qvz,qlz,qrz,qiz,qsz,thz
  real, dimension(1:imax,kts:kte), intent(inout)  :: ncz,niz,nrz,nsz
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
  real                                           :: es1d
  real                                           :: nczodt, nizodt, nrzodt, nszodt
  real, dimension(1:imax,kts:kte)                :: tem
  real                                           :: qswz, qsiz
  real                                           :: theiz, temcc
  real                                           :: qvoqswz, qvoqsiz
  real                                           :: qvzodt, qlzodt, qizodt, qszodt, qrzodt
  real                                           :: tmp2d
  real                                           :: rs0, visc, xka, diffwv, schmidt
  real, dimension(1:imax,kts:kte)                :: viscmu
  real                                           :: qvsbar
!--- microphysical processes
  real                                          :: psacw, psaut, psfw, psfi, praci
  real                                          :: piacr, psaci, psdep, pssub
  real                                          :: psacr, psmlt, psmltevp, prain, praut
  real                                          :: pracw, prevp, pvapor, pclw, pladj
  real                                          :: pcli, pimlt, pihom, pidw, piadj
  real                                          :: pgfr, psnow, qschg, pracs
!---- new snow parameters
  real, parameter                               :: vf1s = 0.65,vf2s = 0.44,            &
                                                   vf1r =0.78,vf2r = 0.31 
  real, parameter                               :: am_c1=0.004,am_c2= 6e-5,  am_c3=0.15
  real, parameter                               :: bm_c1=1.85, bm_c2= 0.003, bm_c3=1.25
  real, parameter                               :: aa_c1=1.28, aa_c2= -0.012,aa_c3=-0.6
  real, parameter                               :: ba_c1=1.5,  ba_c2= 0.0075,ba_c3=0.5
  real, parameter                               :: best_a=1.08 ,  best_b = 0.499
  real                                          :: disp, Dc_liu, eta, R6c        !--- for Liu's autoconversion
  real                                          :: tc0
  real                                          :: mu_c
  !real, dimension(kts:kte)                     :: ab_s,ab_r,ab_riming 
  real                                          :: cap_s    !---- capacitance of snow
  real, dimension(1:imax,kts:kte)               :: am_s,bm_s,av_s,bv_s
  real                                          :: tmp_ss
  real, dimension(1:imax,kts:kte)               :: aa_s,tmp_sa
  real                                          :: ba_s
  real                                          :: mu_s=0.,mu_i=0.,mu_r=0.
 
  ! Adding variable Riz, which will duplicate Ri but be a copy passed upward
  real                                          :: episp0k, dtb, odtb, pi, pio4,       &
                                                   pio6, oxLf, xLvocp, xLfocp, av_r,   &
                                                   av_i, ocdrag, gambp4, gamdp4,       &
                                                   gam4pt5, Cpor, oxmi, gambp3, gamdp3,&
                                                   gambp6, gam3pt5, gam2pt75, gambp5o2,&
                                                   gamdp5o2, cwoxlf, ocp, xni50, es
  real                                          :: gam13
  real                                          :: temc1,save1,save2,xni50mx
  real                                          :: vtr, vts
  real, dimension(1:imax,kts:kte)               :: vtrold, vtsold, vtiold
  real                                          :: xlambdar, xlambdas
  real                                          :: olambdar, olambdas
  ! for terminal velocity flux
  real                                          :: xmr,xms,xmc,dcs,xmr_i
  real                                          :: lamminr, lammaxr,lammins,           &
                                                   lammaxs,lammini, lammaxi
  real                                          :: gambvr1
  real                                          :: lvap
  real                                          :: mi0
  real t_del_tv,del_tv
  real flux, fluxout
  real nflux, nfluxout
  real, dimension(1:imax)                       :: fluxin, nfluxin
  integer                                       :: min_q, max_q
  logical                                       :: notlast
  logical                                       :: mask
  real                                          :: nimlt, nihom, npgfr
  real                                          :: npsacw, npsaut, npraci, npiacr, npsaci
  real                                          :: npsdep, npsacr, npsmlt, npsmltevp
  real                                          :: npraut, npraut_r, npracw, nprevp
  real                                          :: nidep, midep
  real, dimension(1:imax,kts:kte)               :: nvtr, nvts
  real, dimension(1:imax,kts:kte)               :: n0_s 
  real                                          :: n0_r, n0_i, n0_c                  
  real                                          :: lami, lamc
  real, dimension(1:imax,kts:kte)               :: gam_ss, gam_bm_s, gam_bv_ss
  real, dimension(1:imax,kts:kte)               :: gam_bv_s
  real gg21, gg22, gg31, gg32, gg33
  real ratio, gg31c, gg32c, gg33c, mu_c_s
  integer i1, i1p1, n
  real dt

  real fout, fthru, nfout, nfthru, alph
  real qnew
  
  !------------------------------------------------------------------------------------
  
  dt = dt_in/real(njumps)  
  
  mu_c_s = MIN(15., (1000.E6/Nt_c + 2.))
  gg31c = ggamma(4.+mu_c_s)
  gg32c = ggamma(1.+mu_c_s)
  gg33c = ggamma(7.+mu_c_s)

  R6c     = 10.E-6      !---- 10 micron, threshold radius of cloud droplet
  dtb     = dt          !sny
  odtb    = 1./dtb
  pi      = acos(-1.)
  pio4    = pi/4.
  pio6    = pi/6.
  ocp     = 1./cp
  oxLf    = 1./xLf
  xLvocp  = xLv/cp
  xLfocp  = xLf/cp
  Cpor    = cp/Rair
  oxmi    = 1./xmi
  cwoxlf  = cw/xlf 
  av_r    = 2115.0*0.01**(1.-bv_r)
  av_i    = 152.93*0.01**(1.-bv_i)
  ocdrag  = 1./Cdrag
  episp0k = RH*ep2*1000.*svp1

  gambp4  =ggamma(bv_r+4.)
  gamdp4  =ggamma(bv_i+4.)
  gambp3  =ggamma(bv_r+3.)
  gambp6  =ggamma(bv_r+6.)
  gambp5o2=ggamma((bv_r+5.)/2.)
  gamdp5o2=ggamma((bv_i+5.)/2.)
  gambvr1 =ggamma(bv_r+1.)
  gam13   =ggamma(1.+3.)
  
  obp4    =1./(bv_r+4.)
  bp3     =bv_r+3.
  bp5     =bv_r+5.
  bp6     =bv_r+6.
  odp4    =1./(bv_i+4.)
  dp3     =bv_i+3.
  dp5     =bv_i+5.
  dp5o2   =0.5*(bv_i+5.)
    
  dcs     = 125.E-6  ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
  xms     = pi*500.*(dcs)**3/6.   !morr =PI*RHOI*DCS**3/6.=5.11*10e-10
  xmr     = 4./3.*pi*rhowater*(500.E-6)**3
  xmr_i   = 4./3.*pi*rhowater*(25.E-6)**3
  mi0     = 4./3.*3.14*500.*(10.e-6)**3
  !    xmc     =4.17*10e-14 !4./3.*pi*(0.00001)**3*1000.
  lammaxr = 1./20.E-6
  !lamminr = 1./500.E-6
  lamminr = 1./2800.E-6
  lammaxs = 1./10.E-6
  lammins = 1./2000.E-6
  lammaxi = 1./1.E-6
  lammini = 1./(2.*dcs+100.E-6)

  !$acc data create(tem,prez,tothz,rho,dzw,Ri,viscmu,am_s,bm_s,aa_s,av_s,bv_s,tmp_sa, &
  !$acc             gam_ss,gam_bm_s,gam_bv_ss,                                        &
  !$acc             fluxr,fluxi,fluxs,fluxm,fluxf,fevap,fsubl,fauto,fcoll,faccr,      &
  !$acc             effc1d,effr1d,effs1d,effi1d)
  !$acc update device(prez,tothz,rho,dzw)

  do n = 1,njumps
  
    vtrold(:,:)=0.
    vtsold(:,:)=0.
    vtiold(:,:)=0.

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

    select case( lin_aerosolmode )
      case(0)
        do k=kts,kte
          !ncz(k) = 250.*1.E6/rho(k)
          ncz(:,k) = Nt_c/rho(:,k)
        end do
      case(1)
        do k=kts,kte
          ncz(:,k) = zdrop(:,k)/rho(:,k)
        end do
      !case default
      !  write(6,*) "ERROR: Unknown option aerosolmode"
      !  stop
    end select

    do k = kts,kte
      niz(:,k) = min(niz(:,k),0.3E6/rho(:,k))
      nrz(:,k) = max( 0.0,nrz(:,k) )
      nsz(:,k) = max( 0.0,nsz(:,k) )

      qlz(:,k)  =max( 0.0,qlz(:,k) )
      qiz(:,k)  =max( 0.0,qiz(:,k) )
      qvz(:,k)  =max( qvmin,qvz(:,k) )
      qsz(:,k)  =max( 0.0,qsz(:,k) )
      qrz(:,k)  =max( 0.0,qrz(:,k) )
      tem(:,k)  =thz(:,k)*tothz(:,k)

      n0_s(:,k)     =0.
      xlambdar =0.
      xlambdas =0.
      vtr      =0.
      vts      =0.
      vtiold(:,k)   =0.

      !qisten(:,k)   =0.
      !qrsten(:,k)   =0.
      !qssten(:,k)   =0.
      !nisten(:,k)   =0.
      !nrsten(:,k)   =0.
      !nssten(:,k)   =0.

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
      viscmu(:,k)=avisc*tem(:,k)**1.5/(tem(:,k)+120.0)

      ! ---- YLIN, set snow variables
      !
      !---- A+B in depositional growth, the first try just take from Rogers and Yau(1989)
      !         ab_s(k) = lsub*lsub*orv/(tcond(k)*temp(k))+&
      !                   rv*temp(k)/(diffu(k)*qvsi(k))

      !tc0(:)   = tem(:,k)-273.15
      where( rho(:,k)*qlz(:,k) .gt. 1e-5 .AND. rho(:,k)*qsz(:,k) .gt. 1e-5  )
        Ri(:,k) = 1.0/(1.0+6e-5/(rho(:,k)**1.170*qlz(:,k)*qsz(:,k)**0.170)) 
      elsewhere
        Ri(:,k) = 0.
      end where

    end do ! k

    !
    !--- make sure Ri does not decrease downward in a column
    !
    do k = kte-1,kts,-1
      Ri(:,k) = max( Ri(:,k), Ri(:,k+1) )
    end do

    !$acc update device(Ri,tem,viscmu)
    !$acc parallel loop collapse(2) copyout(n0_s,gam_bv_s)                    &
    !$acc   present(tem,rho,Ri,viscmu,am_s,bm_s,aa_s,av_s,bv_s,tmp_sa,gam_ss) &
    !$acc   present(gam_bm_s,gam_bv_ss)
    do k = kts, kte
      do iq = 1,imax

        !--- YLIN, get PI properties
        Ri(iq,k) = max(0.,min(Ri(iq,k),1.0))

        tc0    = min(-0.1, tem(iq,k)-273.15)
        n0_s(iq,k) = min(2.0E8, 2.0E6*exp(-0.12*tc0))
        am_s(iq,k) = am_c1+am_c2*tc0+am_c3*Ri(iq,k)*Ri(iq,k)    !--- Heymsfield 2007
        am_s(iq,k) = max(0.000023,am_s(iq,k))                   !--- use the a_min in table 1 of Heymsfield
        bm_s(iq,k) = bm_c1+bm_c2*tc0+bm_c3*Ri(iq,k)
        bm_s(iq,k) = min(bm_s(iq,k),3.0)                        !---- capped by 3
        !--  converting from cgs to SI unit
        am_s(iq,k) =  10.**(2.*bm_s(iq,k)-3.)*am_s(iq,k)
        aa_s(iq,k) = aa_c1 + aa_c2*tc0 + aa_c3*Ri(iq,k)
        ba_s       = ba_c1 + ba_c2*tc0 + ba_c3*Ri(iq,k)
        !--  convert to SI unit as in paper
        aa_s(iq,k) = (1.e-2)**(2.-ba_s)*aa_s(iq,k)
        !---- get v from Mitchell 1996
        av_s(iq,k) = best_a*viscmu(iq,k)*(2*grav*am_s(iq,k)/rho(iq,k)/ &
                         aa_s(iq,k)/(viscmu(iq,k)**2))**best_b
        bv_s(iq,k) = best_b*(bm_s(iq,k)-ba_s+2.)-1.
        tmp_ss      = bm_s(iq,k)+mu_s+1.
        tmp_sa(iq,k)= ba_s+mu_s+1.
      
        gam_ss(iq,k) = ggamma(tmp_ss)
        gam_bm_s(iq,k) = ggamma(1.+bm_s(iq,k))
        gam_bv_ss(iq,k) = ggamma(bv_s(iq,k)+tmp_ss)
        gam_bv_s(iq,k) = ggamma(bv_s(iq,k)+1.)

      end do  ! iq
    end do ! k
    !$acc end parallel loop

    ! The following variables are needed for the calculation of precipitation fluxes
    !$acc update self(am_s,bm_s,av_s,bv_s,gam_ss,gam_bm_s,gam_bv_ss)

    !***********************************************************************
    ! Calculate precipitation fluxes due to terminal velocities.
    !***********************************************************************
  
    if ( lin_adv==0 ) then
      
      !
      !- Calculate termianl velocity (vt?)  of precipitation q?z
      !- Find maximum vt? to determine the small delta t

      !
      !-- rain
      !
      do iq = 1,imax
        notlast = any( qrz(iq,kts:kte)>1.e-8 )
        if ( notlast ) then
          vtrold(iq,kts:kte) = 0.
          n0_r = 0.
          t_del_tv = 0.
          del_tv = dtb
          do while ( notlast )
        
            min_q = kte
            max_q = kts-1

            ! if no rain, --> minq>maxq --> notlast=False (only case minq>maxq)
            ! if rain --> minq<maxq (always), some vertical points norain--> no lamda, velocity
            do k = kts,kte-1
              if ( qrz(iq,k) > 1.e-8 ) then
                min_q = min(min_q,k)
                max_q = max(max_q,k)
                xlambdar = (pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
                n0_r = nrz(iq,k)*xlambdar
                if ( xlambdar < lamminr ) then
                  xlambdar = lamminr
                  n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
                  nrz(iq,k) = n0_r/xlambdar
                else if ( xlambdar > lammaxr ) then
                  xlambdar = lammaxr
                  n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
                  nrz(iq,k) = n0_r/xlambdar
                end if
                olambdar = 1./xlambdar
                tmp1 = olambdar**bv_r
                vtrold(iq,k) = o6*av_r*gambp4*sqrho*tmp1
                nvtr(iq,k) = av_r*gambvr1*sqrho*tmp1
                del_tv = min(del_tv,0.9*(zz(iq,k)-zz(iq,k-1))/vtrold(iq,k))
              else
                vtrold(iq,k)=0.
                nvtr(iq,k)=0.
                olambdar=0.
              end if
            end do         ! k
      
            !
            !- Check if the summation of the small delta t >=  big delta t
            !             (t_del_tv)          (del_tv)             (dtb)

            if (max_q >= min_q) then

              fluxin(iq)=0.
              nfluxin(iq)=0. ! sny
              t_del_tv=t_del_tv+del_tv
              if ( t_del_tv >= dtb ) then
                notlast=.false.
                del_tv=dtb+del_tv-t_del_tv
              end if

              do k = max_q,min_q,-1
                fluxout=rho(iq,k)*vtrold(iq,k)*qrz(iq,k)
                flux=(fluxin(iq)-fluxout)/rho(iq,k)/dzw(iq,k)
                !tmpqrz(iq)=qrz(iq,k)
                qrz(iq,k)=qrz(iq,k)+del_tv*flux
                fluxin(iq)=fluxout

                nfluxout=rho(iq,k)*nvtr(iq,k)*nrz(iq,k)
                nflux=(nfluxin(iq)-nfluxout)/rho(iq,k)/dzw(iq,k)
                nrz(iq,k)=nrz(iq,k)+del_tv*nflux
                nfluxin(iq)=nfluxout
                !qrsten(iq,k)=flux
                !nrsten(iq,k)=nflux
              end do       !k

              if ( min_q == 1 ) then
                pptrain(iq) = pptrain(iq) + fluxin(iq)*del_tv  
              else      
                qrz(iq,min_q-1)=qrz(iq,min_q-1)+del_tv*  &
                               fluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
                nrz(iq,min_q-1)=nrz(iq,min_q-1)+del_tv*  &
                                nfluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
              end if

            else
              notlast=.false.
            end if ! maxq>minq

          END DO      ! while(notlast)
        end if
      end do ! iq
 
      !
      !-- snow
      !
      do iq = 1,imax
        notlast=any( qsz(iq,kts:kte)>1.e-8 )
        if ( notlast ) then
          vtsold(iq,kts:kte) = 0.
          t_del_tv=0.
          del_tv=dtb
          do while (notlast)

            min_q=kte
            max_q=kts-1

            do k=kts,kte-1
              if (qsz(iq,k) > 1.e-8) then
                min_q = min(min_q,k)
                max_q = max(max_q,k)
                ! Zhao 2022 - Row 2 Table 2 or Lin 2011 - Formula A3
                xlambdas=(am_s(iq,k)*gam_ss(iq,k)*nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
                ! Zhao 2022 - Row 1 Table 2
                n0_s(iq,k)=nsz(iq,k)*xlambdas
                if (xlambdas<lammins) then
                  xlambdas= lammins
                  n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1.)*qsz(iq,k)/gam_bm_s(iq,k)/am_s(iq,k)
                  nsz(iq,k) = n0_s(iq,k)/xlambdas
                else if (xlambdas>lammaxs) then
                  xlambdas = lammaxs
                  n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1.)*qsz(iq,k)/gam_bm_s(iq,k)/am_s(iq,k)
                  nsz(iq,k) = n0_s(iq,k)/xlambdas
                end if
                olambdas=1./xlambdas
                tmp1 = olambdas**bv_s(iq,k)
                ! Zhao 2022 - Row 3 Table 2
                vtsold(iq,k)= sqrho*av_s(iq,k)*gam_bv_ss(iq,k)/ &
                   gam_ss(iq,k)*tmp1
                ! Zhao 2022 - Row 4 Table 2
                nvts(iq,k)=sqrho*av_s(iq,k)*gam_bv_s(iq,k)*tmp1
                del_tv=min(del_tv,0.9*(zz(iq,k)-zz(iq,k-1))/vtsold(iq,k))
              else
                vtsold(iq,k)=0.
                nvts(iq,k)=0.
                olambdas=0.
              endif
            end do       ! k

            !
            !- Check if the summation of the small delta t >=  big delta t
            !             (t_del_tv)          (del_tv)             (dtb)

            if (max_q >= min_q) then

              fluxin(iq) = 0.
              nfluxin(iq) = 0.
              t_del_tv=t_del_tv+del_tv
              if ( t_del_tv >= dtb ) then
                notlast=.false.
                del_tv=dtb+del_tv-t_del_tv
              endif

              do k = max_q,min_q,-1
                fluxout=rho(iq,k)*vtsold(iq,k)*qsz(iq,k)
                flux=(fluxin(iq)-fluxout)/rho(iq,k)/dzw(iq,k)
                qsz(iq,k)=qsz(iq,k)+del_tv*flux
                qsz(iq,k)=max(0.,qsz(iq,k))
                fluxin(iq)=fluxout

                nfluxout=rho(iq,k)*nvts(iq,k)*nsz(iq,k)
                nflux   =(nfluxin(iq)-nfluxout)/rho(iq,k)/dzw(iq,k)
                nsz(iq,k)  =nsz(iq,k)+del_tv*nflux
                nfluxin(iq) =nfluxout
                !qssten(iq,k)=flux
                !nssten(iq,k)=nflux
              end do       ! k

              if ( min_q == 1 ) then
                pptsnow(iq) = pptsnow(iq) + fluxin(iq)*del_tv  
              else
                qsz(iq,min_q-1)=qsz(iq,min_q-1)+del_tv*  &
                         fluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
                nsz(iq,min_q-1)=nsz(iq,min_q-1)+del_tv*  &
                         nfluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
              end if

            else
              notlast=.false.
            end if ! maxq>minq

          END DO       ! while(notlast)
        end if
      end do        ! iq = 1,imax
 
      !
      !-- cloud ice  (03/21/02) using Heymsfield and Donner (1990) Vi=3.29*qi^0.16
      !
      do iq = 1,imax
        qiz(iq,kts:kte) = qiz(iq,kts:kte)
        notlast=any( qiz(iq,kts:kte)>1.e-8 )
        if ( notlast ) then
          vtiold(iq,kts:kte) = 0.
          t_del_tv=0.
          del_tv=dtb
          DO while (notlast)

            min_q=kte
            max_q=kts-1

            do k=kts,kte-1
              if (qiz(iq,k) > 1.e-8) then
                min_q = min(min_q,k)
                max_q = max(max_q,k)
                vtiold(iq,k) = 3.29 * (rho(iq,k)* qiz(iq,k))** 0.16  ! Heymsfield and Donner
                del_tv=min(del_tv,0.9*(zz(iq,k)-zz(iq,k-1))/vtiold(iq,k))
              else
                vtiold(iq,k)=0.
              endif
            enddo       ! k
      
            !
            !- Check if the summation of the small delta t >=  big delta t
            !             (t_del_tv)          (del_tv)             (dtb)
            if (max_q >= min_q) then

              fluxin(iq) = 0.
              nfluxin(iq) = 0.
              t_del_tv=t_del_tv+del_tv
              if ( t_del_tv >= dtb ) then
                notlast=.false.
                del_tv=dtb+del_tv-t_del_tv
              endif

              do k = max_q,min_q,-1
                fluxout=rho(iq,k)*vtiold(iq,k)*qiz(iq,k)
                flux=(fluxin(iq)-fluxout)/rho(iq,k)/dzw(iq,k)
                qiz(iq,k)=qiz(iq,k)+del_tv*flux
                qiz(iq,k)=max(0.,qiz(iq,k))
                fluxin(iq)=fluxout

                nfluxout=rho(iq,k)*vtiold(iq,k)*niz(iq,k)
                nflux=(nfluxin(iq)-nfluxout)/rho(iq,k)/dzw(iq,k)
                niz(iq,k)=niz(iq,k)+del_tv*nflux
                niz(iq,k)=max(0.,niz(iq,k))
                !niz(iq,k) = min(niz(iq,k),0.3E6/rho(iq,k))
                nfluxin(iq)=nfluxout
                !qisten(iq,k)=flux
                !nisten(iq,k)=nflux
              end do       ! k

              if ( min_q == 1 ) then
                pptice(iq) = pptice(iq) + fluxin(iq)*del_tv  
              else
                qiz(iq,min_q-1)=qiz(iq,min_q-1)+del_tv*  &
                             fluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
                niz(iq,min_q-1)=niz(iq,min_q-1)+del_tv*  &
                             nfluxin(iq)/rho(iq,min_q-1)/dzw(iq,min_q-1)
              end if

            else
              notlast=.false.
            end if   ! maxq>minq

          END DO       ! while(notlast)
        end if
      end do        ! iq = 1,imax
    
    else ! lin_adv==1
      ! Flux method based on Rotstayn 1997
      
      !
      !-- rain
      !
      vtrold(:,:) = 0.
      n0_r = 0.  
      fluxin(:) = 0.
      nfluxin(:) = 0.
      do k = kte-1,kts,-1
        do iq = 1,imax
          ! if no rain, --> minq>maxq --> notlast=False (only case minq>maxq)
          ! if rain --> minq<maxq (always), some vertical points norain--> no lamda, velocity
          qnew = qrz(iq,k) + fluxin(iq)/rho(iq,k)  
          if ( qnew > 1.e-8 ) then  
            xlambdar = (pi*rhowater*nrz(iq,k)/qnew)**(1./3.)   !zx
            n0_r = nrz(iq,k)*xlambdar
            if ( xlambdar < lamminr ) then
              xlambdar = lamminr
              n0_r = xlambdar**4*qnew/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            else if ( xlambdar > lammaxr ) then
              xlambdar = lammaxr
              n0_r = xlambdar**4*qnew/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            end if
            olambdar = 1./xlambdar
            tmp1 = olambdar**bv_r
            vtrold(iq,k) = o6*av_r*gambp4*sqrho*tmp1
            nvtr(iq,k) = av_r*gambvr1*sqrho*tmp1
            alph = dtb*vtrold(iq,k)/dzw(iq,k)
            fout = 1. - exp(-alph)
            fthru = 1. - fout/alph
            alph = dtb*nvtr(iq,k)/dzw(iq,k)
            nfout = 1. - exp(-alph)
            nfthru = 1. - nfout/alph           
          else
            olambdar = 0.
            vtrold(iq,k) = 0.
            nvtr(iq,k) = 0.
            fout = 0.
            fthru = 0.
            nfout = 0.
            nfthru = 0.
          end if
          fluxout = rho(iq,k)*qrz(iq,k)*fout*dzw(iq,k)
          flux = fluxin(iq)*(1.-fthru) - fluxout
          qrz(iq,k) = qrz(iq,k) + flux/(rho(iq,k)*dzw(iq,k))
          fluxin(iq) = fluxout + fluxin(iq)*fthru
          nfluxout = rho(iq,k)*nrz(iq,k)*nfout*dzw(iq,k)
          nflux = nfluxin(iq)*(1.-nfthru) - nfluxout
          nrz(iq,k) = nrz(iq,k) + nflux/(rho(iq,k)*dzw(iq,k))
          nfluxin(iq) = nfluxout + nfluxin(iq)*nfthru
        end do ! iq  
      end do   ! k
      do iq = 1,imax
        pptrain(iq) = pptrain(iq) + fluxin(iq)
      end do ! iq
    
      !
      !-- snow
      !
      vtsold(:,:) = 0.
      fluxin(:) = 0.
      nfluxin(:) = 0. ! sny
      do k = kte-1,kts,-1
        do iq = 1,imax
          qnew = qsz(iq,k) + fluxin(iq)/rho(iq,k)  
          if ( qnew > 1.e-8 ) then
            ! Zhao 2022 - Row 2 Table 2 or Lin 2011 - Formula A3
            xlambdas = (am_s(iq,k)*gam_ss(iq,k)*nsz(iq,k)/qnew)**(1./bm_s(iq,k))
            ! Zhao 2022 - Row 1 Table 2
            n0_s(iq,k) = nsz(iq,k)*xlambdas
            if ( xlambdas<lammins ) then
              xlambdas= lammins
              n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1.)*qnew/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,k)/xlambdas
            else if (xlambdas>lammaxs) then
              xlambdas = lammaxs
              n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1.)*qnew/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,k)/xlambdas
            end if
            olambdas = 1./xlambdas
            tmp1 = olambdas**bv_s(iq,k)
            ! Zhao 2022 - Row 3 Table 2
            vtsold(iq,k)= sqrho*av_s(iq,k)*gam_bv_ss(iq,k)/gam_ss(iq,k)*tmp1
            ! Zhao 2022 - Row 4 Table 2
            nvts(iq,k) = sqrho*av_s(iq,k)*gam_bv_s(iq,k)*tmp1
            alph = dtb*vtsold(iq,k)/dzw(iq,k)
            fout = 1. - exp(-alph)
            fthru = 1. - fout/alph
            alph = dtb*nvts(iq,k)/dzw(iq,k)
            nfout = 1. - exp(-alph)
            nfthru = 1. - nfout/alph  
          else
            olambdas = 0.
            vtsold(iq,k) = 0.
            nvts(iq,k) = 0.
            fout = 0.
            fthru = 0.
            nfout = 0.
            nfthru = 0.  
          endif
          fluxout = rho(iq,k)*qsz(iq,k)*fout*dzw(iq,k)
          flux = fluxin(iq)*(1.-fthru)-fluxout
          qsz(iq,k) = qsz(iq,k) + flux/(rho(iq,k)*dzw(iq,k))
          qsz(iq,k) = max(0.,qsz(iq,k))
          fluxin(iq) = fluxout + fluxin(iq)*fthru
          nfluxout = rho(iq,k)*nsz(iq,k)*nfout*dzw(iq,k)
          nflux = nfluxin(iq)*(1.-nfthru)-nfluxout
          nsz(iq,k) = nsz(iq,k) + nflux/(rho(iq,k)*dzw(iq,k))
          nfluxin(iq) = nfluxout + nfluxin(iq)*nfthru
        end do ! iq  
      end do   ! k
      do iq = 1,imax
        pptsnow(iq) = pptsnow(iq) + fluxin(iq)
      end do
 
      !
      !-- cloud ice  (03/21/02) using Heymsfield and Donner (1990) Vi=3.29*qi^0.16
      !
      vtiold(:,:) = 0.
      fluxin(:) = 0.
      nfluxin(:) = 0.
      do k = kte-1,kts,-1
        do iq = 1,imax
          qnew = qiz(iq,k) + fluxin(iq)/rho(iq,k)  
          if ( qnew > 1.e-8 ) then
            vtiold(iq,k) = 3.29*(rho(iq,k)*qnew)**0.16  ! Heymsfield and Donner
            alph = dtb*vtiold(iq,k)/dzw(iq,k)
            fout = 1. - exp(-alph)
            fthru = 1. - fout/alph
            nfout = fout
            nfthru = fthru
          else
            vtiold(iq,k) = 0.
            fout = 0.
            fthru = 0.
            nfout = 0.
            nfthru = 0.
          end if
          fluxout = rho(iq,k)*qiz(iq,k)*fout*dzw(iq,k)
          flux = fluxin(iq)*(1.-fthru)-fluxout
          qiz(iq,k) = qiz(iq,k) + flux/(rho(iq,k)*dzw(iq,k))
          qiz(iq,k) = max(0.,qiz(iq,k))
          fluxin(iq) = fluxout + fluxin(iq)*fthru
          nfluxout = rho(iq,k)*niz(iq,k)*nfout*dzw(iq,k)
          nflux = nfluxin(iq)*(1.-nfthru)-nfluxout
          niz(iq,k) = niz(iq,k) + nflux/(rho(iq,k)*dzw(iq,k))
          niz(iq,k) = max(0.,niz(iq,k))
          !niz(iq,k) = min(niz(iq,k),0.3E6/rho(iq,k))
          nfluxin(iq) = nfluxout + nfluxin(iq)*nfthru
        end do ! iq  
      end do   ! k
      do iq = 1,imax
        pptice(iq) = pptice(iq) + fluxin(iq)
      end do
    
    end if  ! if lin_adv==0 ..else..  

    ! save for output
    vi(:,:)    = vtiold(:,:)
    !vs(:,:)    = vtsold(:,:)

    !$acc parallel loop collapse(2) copy(ncz,niz,nrz,nsz,thz,qvz,qlz,qiz,qsz,qrz)  &
    !$acc   copyin(n0_s)                                                           &
    !$acc   present(fluxr,fluxi,fluxs,fluxm,fluxf,fevap,fsubl,fauto,fcoll,faccr)   &
    !$acc   present(effc1d,effr1d,effs1d,effi1d)                                   &
    !$acc   present(tem,prez,tothz,rho,dzw,Ri,viscmu,am_s,bm_s,aa_s,av_s,bv_s)     &
    !$acc   present(tmp_sa,gam_ss,gam_bm_s,gam_bv_ss)
    do k = kts,kte
      do iq = 1,imax

        ! Microphpysics processes
          
        temcc = tem(iq,k)-273.15  
        
        es1d =1000.*svp1*exp( svp2*temcc/(tem(iq,k)-svp3) )  !--- RY89 Eq(2.17)
        qswz =ep2*es1d/max(prez(iq,k)-es1d,0.1)
 
        if (tem(iq,k) .lt. 233.15 ) then
          es1d=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
          qsiz=ep2*es1d/max(prez(iq,k)-es1d,0.1)
          qswz=qsiz
        else if (tem(iq,k) .lt. 273.15 ) then
          es1d=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
          qsiz=ep2*es1d/max(prez(iq,k)-es1d,0.1)
        else
          qsiz=qswz
        endif
          
        nczodt = max( 0.0,odtb*ncz(iq,k) )
        nizodt = max( 0.0,odtb*niz(iq,k) )
        nrzodt = max( 0.0,odtb*nrz(iq,k) )
        nszodt = max( 0.0,odtb*nsz(iq,k) )

        qvoqswz  =qvz(iq,k)/qswz
        qvoqsiz  =qvz(iq,k)/qsiz
        qvzodt=max( 0.,odtb*qvz(iq,k) )
        qlzodt=max( 0.,odtb*qlz(iq,k) )
        qizodt=max( 0.,odtb*qiz(iq,k) )
        qszodt=max( 0.,odtb*qsz(iq,k) )
        qrzodt=max( 0.,odtb*qrz(iq,k) )
        
        theiz = thz(iq,k)+(xlvocp*qvz(iq,k)-xlfocp*qiz(iq,k))/tothz(iq,k)
        
        rs0 = ep2*1000.*svp1/(prez(iq,k)-1000.*svp1)
        visc = viscmu(iq,k)/rho(iq,k)
        xka = axka*viscmu(iq,k)
        diffwv = adiffwv*tem(iq,k)**1.81/prez(iq,k)
        schmidt = visc/diffwv
        
        cap_s = 0.25*(1.+Ri(iq,k))

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

        n0_r     = 0.
        n0_i     = 0.
        n0_c     = 0.
        
        mu_c = mu_c_s
        
        lamc =0.
        lami =0.

        xlambdar = 0.
        xlambdas = 0.

        olambdar = 0.
        olambdas = 0.

        psacw   =0.                  ! accretion of cloud water by snow
        psaut   =0.                  ! ice crystal aggregation to snow
        psfw    =0.                  ! BERGERON process to transfer cloud water to snow
        psfi    =0.                  ! BERGERON process to transfer cloud ice to snow
        praci   =0.                  ! cloud ice accretion by rain
        piacr   =0.                  ! rain accretion by cloud ice
        psaci   =0.                  ! ice crystal accretion by snow
        psdep   =0.                  ! deposition of snow
        pssub   =0.                  ! sublimation of snow (T<0)
        psacr   =0.                  ! accretion of rain by snow
        psmlt   =0.                  ! melting of snow
        psmltevp=0.                  ! evaporation of melting snow (T>0)
        prain   =0.                  ! sum all process for rain
        praut   =0.                  ! autoconversion of rain
        pracw   =0.                  ! accretion of cloud water by rain
        prevp   =0.                  ! evaporation of rain
        pvapor  =0.                  ! sum all process for water vapor to determine qvz
        pclw    =0.                  ! sum all process for cloud liquid to determine qlz
        pladj   =0.                  ! saturation adjustment for ql
        pcli    =0.                  ! sum all process for cloud ice to determine qiz
        pimlt   =0.                  ! melting of ice crystal >0.
        pihom   =0.                  ! homogeneous nucleation <-40
        pidw    =0.                  ! production of cloud ice by BERGERON process
        piadj   =0.                  ! saturation adjustment for qi
        pgfr    =0.                  ! feezing of rain to form graupel (added to PI)
        psnow   =0.                  ! sum all process for snow
        qschg   =0.                  ! = psnow / unsure
        pracs   =0.

        nidep = 0.
        midep = 0.
        npraut_r = 0.
        nprevp = 0.
        npracw = 0.
        npraut = 0.
        npsmltevp = 0.
        npsmlt = 0.
        npgfr = 0.
        npsacr = 0.
        npsdep = 0.
        npsacw = 0.
        npsaci = 0.
        npiacr = 0.
        npraci = 0.
        npsaut = 0.
        !npssub = 0.
    
        fluxr(iq,k) = qrzodt*dzw(iq,k)*rho(iq,k)
        fluxi(iq,k) = qizodt*dzw(iq,k)*rho(iq,k)
        fluxs(iq,k) = qszodt*dzw(iq,k)*rho(iq,k)
          
        tmp2d=qiz(iq,k)+qlz(iq,k)+qsz(iq,k)+qrz(iq,k)
        mask = .not.(qvz(iq,k)+qlz(iq,k)+qiz(iq,k) < qsiz .and. &
                     tmp2d==0. )      
    
        if ( mask ) then
          !gg11(iq) = ggamma(tmp_ss(iq,k))
          !gg12(iq) = ggamma(1.+bm_s(iq,k))
          !gg13(iq) = ggamma(bv_s(iq,k)+tmp_ss(iq,k))
          !gg14(iq) = ggamma(bv_s(iq,k)+1.)
          gg21 = ggamma(bv_s(iq,k)+tmp_sa(iq,k))
          gg22 = ggamma(2.5+0.5*bv_s(iq,k)+mu_s)
          gg31 = gg31c
          gg32 = gg32c
          gg33 = gg33c
          !gg41(iq) = ggamma(tmp_ss(iq,k))
          !gg42(iq) = ggamma(1.+bm_s(iq,k))
          !gg51(iq) = ggamma(4.+mu_c)
          !gg52(iq) = ggamma(1.+mu_c)      

          !
          !! calculate terminal velocity of rain
          !    
        
          if (qrz(iq,k) .gt. 1.e-8) then
            !  tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrz(k))
            !  xlambdar=sqrt(tmp1)
            !  olambdar=1.0/xlambdar
            !  vtr=o6*av_r*gambp4*sqrho*olambdar**bv_r
            xlambdar=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
            n0_r=nrz(iq,k)*xlambdar
            if (xlambdar<lamminr) then
              xlambdar = lamminr
              n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            else if (xlambdar>lammaxr) then
              xlambdar = lammaxr
              n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            end if
            olambdar=1.0/xlambdar
            tmp1 = olambdar**bv_r
            vtr=o6*av_r*gambp4*sqrho*tmp1
            !nvtr(iq,k)=av_r*gambvr1*sqrho*tmp1
          else
            vtr=0.
            olambdar=0.
            xlambdar=0.
            !nvtr(iq,k)=0.
          end if  ! qrz

          !!
          !!! calculate terminal velocity of snow
          !!

          if (qsz(iq,k) > 1.e-8) then
            xlambdas=(am_s(iq,k)*gam_ss(iq,k)*nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
            n0_s(iq,k)=nsz(iq,k)*xlambdas
            if (xlambdas.lt.lammins) then
              xlambdas= lamminS
              n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1)*qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,K)/xlambdas
            else if (xlambdas.gt.lammaxs) then
              xlambdas = lammaxs
              n0_s(iq,k) = xlambdas**(bm_s(iq,k)+1)*qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,k) = n0_s(iq,K)/xlambdas
            end if
            olambdas=1.0/xlambdas
            tmp1 = olambdas**bv_s(iq,k)
            vts= sqrho*av_s(iq,k)*gam_bv_ss(iq,k)/ &
                   gam_ss(iq,k)*tmp1
            !nvts(iq,k)=sqrho*av_s(iq,k)*gam_bv_s(iq,k)*tmp1
          else
            vts=0.
            olambdas=0.
            xlambdas=0.
            !nvts(iq,k)=0.
          endif
        
          !---------- start of snow/ice processes below freezing

          if (tem(iq,k) < 273.15) then

            !
            ! ice nucleation, cooper curve

            if ((qvoqswz>=0.999.and.temcc<=-8.).or. &
              qvoqsiz>=1.08) then
              nidep = 5.*exp(0.304*(273.15-tem(iq,k))) ! m-3
              nidep = min(nidep, 500.e8)           !5.e8) sny ! limit to 500 L-1
              nidep = max(nidep/rho(iq,k), 0.)     ! convert to kg-1
              nidep = (nidep - niz(iq,k))*odtb
              midep = nidep*mi0
            end if
            !***********************************************************************
            !*********        snow production processes for T < 0 C       **********
            !***********************************************************************
            !
            ! (1) ICE CRYSTAL AGGREGATION TO SNOW (Psaut): Lin (21)
            !!    psaut=alpha1*(qi-qi0)
            !!    alpha1=1.0e-3*exp(0.025*(T-T0))
            !
            alpha1=1.e-3*exp( 0.025*temcc )

            ! ---------------------------------------------------------------
            if(temcc .lt. -20.0) then
              tmp1=-7.6+4.*exp( -0.2443e-3*(abs(temcc)-20.)**2.455 )
              qic=1.e-3*exp(tmp1)/rho(iq,k)
            else
              qic=qi0
            end if
            !----------------------------------------------------------------

            !qic = qi0  ! sny: OFF temp

            tmp1=odtb*(qiz(iq,k)-qic)*(1.0-exp(-alpha1*dtb))
            psaut=max( 0., tmp1 )
            npsaut=max( 0., psaut/xms)

            !
            ! (2) BERGERON PROCESS TRANSFER OF CLOUD WATER TO SNOW (Psfw)
            !     this process only considered when -31 C < T < 0 C
            !     Lin (33) and Hsie (17)
            ! 
            !!
            !!    parama1 and parama2 functions must be user supplied
            !!

            if( qlz(iq,k) > 1.e-8 ) then
              temc1=max(-30.99,temcc)
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
              xni50mx=qlzodt/tmp2
              !
              !!   number of 50 micron crystals produced
              !
              xni50=qiz(iq,k)*( 1.0-exp(-dtb*odtberg) )/xmi50
              xni50=min(xni50,xni50mx)
              !
              tmp3=odtb*tmp2/save2*( 1.0-exp(-save2*xni50*dtb) )
              psfw=min( tmp3,qlzodt )
              !
              ! (3) REDUCTION OF CLOUD ICE BY BERGERON PROCESS (Psfi): Lin (34)
              !     this process only considered when -31 C < T < 0 C
              !
              tmp1=xni50*xmi50-psfw
              psfi=min(tmp1,qizodt)
            end if
        
            if(qrz(iq,k) > 0.) then  ! go to 1000

              !
              ! Processes (4) and (5) only need when qrz > 0.0
              !
              ! (4) CLOUD ICE ACCRETION BY RAIN (Praci): Lin (25)
              !     produce PI
              !
              eri=1.0
              save1=pio4*eri*n0_r*av_r*sqrho
              tmp1= save1*gambp3*olambdar**bp3
              praci=qizodt*( 1.0-exp(-tmp1*dtb) )
              npraci=niz(iq,k)*tmp1

              !
              ! (5) RAIN ACCRETION BY CLOUD ICE (Piacr): Lin (26)
              !
              tmp2=qiz(iq,k)*save1*rho(iq,k)*pio6*rhowater*gambp6*oxmi* &
                   olambdar**bp6
              piacr=min( tmp2,qrzodt )
              tmp1 = olambdar**bp3 
              npiacr=pio4*eri*nrz(iq,k)*av_r*niz(iq,k)*gambp3*tmp1 !--wdm6
            end if !1000    continue

            if(qsz(iq,k) > 0.) then !go to 1200
              !
              ! Compute the following processes only when qsz > 0.0
              !
              !
              ! (6) ICE CRYSTAL ACCRETION BY SNOW (Psaci): Lin (22)
              !
              esi=exp( 0.025*temcc )
              tmp1 = olambdas**(bv_s(iq,k)+tmp_sa(iq,k))
              save1 = aa_s(iq,k)*sqrho*n0_s(iq,k)*gg21*tmp1

              tmp1=esi*save1
              psaci=qizodt*( 1.-exp(-tmp1*dtb) )
              npsaci=min( tmp1*niz(iq,k),nizodt)
              !
              ! (7) CLOUD WATER ACCRETION BY SNOW (Psacw): Lin (24)
              !
              esw=1.0
              tmp1=esw*save1
              psacw=qlzodt*( 1.-exp(-tmp1*dtb) )
              npsacw=min(tmp1*ncz(iq,k),ncz(iq,k))

              ! recalculate the saturatuin temperature
              !
              ! (8) DEPOSITION/SUBLIMATION OF SNOW (Psdep/Pssub): Lin (31)
              !     includes consideration of ventilation effect
              !
              tmpa=rvapor*xka*tem(iq,k)*tem(iq,k)
              tmpb=xls*xls*rho(iq,k)*qsiz*diffwv
              tmpc=tmpa*qsiz*diffwv
              abi=4.0*pi*cap_s*(qvoqsiz-1.0)*tmpc/(tmpa+tmpb)
              tmp1=av_s(iq,k)*sqrho*olambdas**(5.+bv_s(iq,k)+2.*mu_s)/visc

              !---- YLIN, here there is some approximation assuming mu_s =1, so gamma(2)=1, etc.

              tmp2= abi*n0_s(iq,k)*( vf1s*olambdas*olambdas+ &
                   vf2s*schmidt**0.33334* &
                   gg22*sqrt(tmp1) )
              tmp3=odtb*( qvz(iq,k)-qsiz )
              tmp3=min(tmp3,0.)

              if( tmp2 <= 0.) then
                tmp2=max( tmp2,tmp3)
                pssub=max( tmp2,-qszodt )
                psdep=0.0
              else
                psdep=min( tmp2,tmp3 )
                pssub=0.0
              end if
              if(qsz(iq,k) >= 0.) then
                !npssub=pssub*nsz(iq,k)/qsz(iq,k)
                npsdep=npsdep*nsz(iq,k)/qsz(iq,k)
              else
                !npssub=pssub/xms
                npsdep=npsdep/xms
              end if

              if (qrz(iq,k) > 0.) then !go to 1200
                !
                ! Compute processes (9) and (10) only when qsz > 0.0 and qrz > 0.0
                ! these two terms need to be refined in the future, they should be equal
                !
                ! (9) ACCRETION OF SNOW BY RAIN (Pracs): Lin (27)
                !
                esr=1.0
                tmpa=olambdar**2
                tmpb=olambdas**2
                tmpc=olambdar*olambdas
                tmp1=pi**2*esr*n0_r*n0_s(iq,k)*    &
                          abs( vtr-vts )/rho(iq,k)
                ! tmp1=pi*pi*esr*n0_r*N0_s(k)*            &
                ! ( (1.2*vtr-0.95*vts)**2+0.08*vtr*vts)**0.5/rho(k)
                tmp2=tmpb**2*olambdar*(5.0*tmpb+2.0*tmpc+0.5*tmpa)
                tmp3=tmp1*rhosnow*tmp2
                pracs=min( tmp3,qszodt )
                !
                ! (10) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
                !
                tmp3=tmpa**2*olambdas*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
                tmp4=tmp1*rhowater*tmp3
                psacr=min( tmp4,qrzodt )
                tmp1=0.25*pi*esr*n0_r*n0_s(iq,k)*abs( vtr-vts )
                tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2.*tmpb)
                tmp3=tmp1*tmp2
                npsacr=min( tmp3,nrzodt )
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
                  tmp1=olambdar**3
                  tmp2=20.*pi**2*Bp*n0_r*rhowater/rho(iq,k)*  &
                     (exp(-Ap*temcc)-1.)*tmp1**2*olambdar
                  pgfr=min( tmp2,qrzodt )
                  npgfr=pi*Bp*n0_r*tmpa**2*(exp(-Ap*temcc)-1.0)
                else
                  pgfr=0
                  npgfr=0.
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
              tmp1 = olambdas**(bv_s(iq,k)+tmp_sa(iq,k))
              save1 =aa_s(iq,k)*sqrho*n0_s(iq,k)*gg21*tmp1

              tmp1=esw*save1
              psacw=qlzodt*( 1.0-exp(-tmp1*dtb) )
              npsacw=tmp1*ncz(iq,k)
              !
              ! (2) ACCRETION OF RAIN BY SNOW (Psacr): Lin (28)
              !
              esr=1.
              tmpa=olambdar*olambdar
              tmpb=olambdas*olambdas
              tmpc=olambdar*olambdas
              tmp1=pi*pi*esr*n0_r*n0_s(iq,k)*   &
                      abs( vtr-vts )/rho(iq,k)
              ! tmp1=pi*pi*esr*n0_r*N0_s(k)*            &
              ! ( (1.2*vtr-0.95*vts)**2+0.08*vtr*vts)**0.5/rho(k)
              tmp2=tmpa*tmpa*olambdas*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
              tmp3=tmp1*rhowater*tmp2
              psacr=min( tmp3,qrzodt )

              tmp1=0.25*pi*esr*n0_r*n0_s(iq,k)*abs( vtr-vts )
              tmp2=tmpc*(2.0*tmpa+1.0*tmpc+2*tmpb)
              tmp3=tmp1*tmp2
              npsacr=min( tmp3,nrzodt )
              !
              ! (3) MELTING OF SNOW (Psmlt): Lin (32)
              !     Psmlt is negative value
              !
              delrs=rs0 - qvz(iq,k)
              term1=2.0*pi/rho(iq,k)*( xlv*diffwv*rho(iq,k)*delrs- &
                    xka*temcc )
              tmp1= av_s(iq,k)*sqrho*olambdas**(5.+bv_s(iq,k)+2.*mu_s)/visc
              tmp2= n0_s(iq,k)*( vf1s*olambdas*olambdas+ &
                    vf2s*schmidt**0.33334* &
                    gg22*sqrt(tmp1) )
              tmp3=term1*oxlf*tmp2-cwoxlf*temcc*( psacw+psacr )
              tmp4=min(0.0,tmp3)
              psmlt=max( tmp4,-qszodt )

              if(qsz(iq,k) >= 0.) then
                npsmlt=psmlt*nsz(iq,k)/qsz(iq,k)
              else
                npsmlt=psmlt/xms
              end if
              !
              ! (4) EVAPORATION OF MELTING SNOW (Psmltevp): HR (A27)
              !     but use Lin et al. coefficience
              !     Psmltevp is a negative value
              !
              tmpa=rvapor*xka*tem(iq,k)*tem(iq,k)
              tmpb=xlv*xlv*rho(iq,k)*qswz*diffwv
              tmpc=tmpa*qswz*diffwv
              tmpd=min( 0.0,(qvoqswz-0.90)*qswz*odtb )

              abr=2.0*pi*(qvoqswz-0.90)*tmpc/(tmpa+tmpb)
              !
              !**** allow evaporation to occur when RH less than 90%
              !**** here not using 100% because the evaporation cooling
              !**** of temperature is not taking into account yet; hence,
              !**** the qsw value is a little bit larger. This will avoid
              !**** evaporation can generate cloud.
          
              tmp1=av_s(iq,k)*sqrho*olambdas**(5.+bv_s(iq,k)+2.*mu_s)/visc
              tmp2=n0_s(iq,k)*( vf1s*olambdas*olambdas+ &
                   vf2s*schmidt**0.33334* &
                   gg22*sqrt(tmp1) )
              tmp3=min(0.0,tmp2)
              tmp3=max( tmp3,tmpd )
              psmltevp=max( tmp3,-qszodt )
              if(qsz(iq,k) >= 0.) then
                npsmltevp=psmltevp*nsz(iq,k)/qsz(iq,k)
              else
                npsmltevp=psmltevp/xmr
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

          praut = 0.0
          npraut = 0.0
          npraut_r = 0.0
        
          if (qlz(iq,k) > 1.e-8) then
            
            mu_c = min(15., (1000.E6/ncz(iq,k) + 2.))
            gg31 = ggamma(4.+mu_c)
            gg32 = ggamma(1.+mu_c)
            gg33 = ggamma(7.+mu_c)
            
            lamc = (ncz(iq,k)*rhowater*pi*gg31/(6.*qlz(iq,k)*gg32))**(1./3)
            Dc_liu = (gg33/gg32)**(1./6.)/lamc !----- R6 in m
            if (Dc_liu > R6c) then
              disp = 1./(mu_c+1.)              !--- square of relative dispersion
              eta  = (0.75/pi/(1.e-3*rhowater))**2*1.9e11*((1.+3.*disp)*(1.+4.*disp)*&
                     (1.+5.*disp)/(1.+disp)/(1.+2.*disp))
              praut = eta*(1.e-3*rho(iq,k)*qlz(iq,k))**3/(1.e-6*ncz(iq,k))  !--- g cm-3 s-1
              praut = praut/(1.e-3*rho(iq,k))                               !--- kg kg-1 s-1
              npraut_r = praut/xmr                                          !--- kg kg-1 s-1
              npraut = praut/qlz(iq,k)*ncz(iq,k)                            !--- kg kg-1 s-1
              npraut = praut/xmr                                            !--- kg kg-1 s-1
            end if

          end if
          ! if (qlz(k) .gt. 1e-6) then
          ! praut=1350.*qlz(k)**2.47*  &
          ! (ncz(k)/1.e6*rho(k))**(-1.79)
          ! npraut_r = praut/xmr
          ! npraut = praut/(qlz(k)/ncz(k))
          ! npraut = MIN(npraut,nczodt)
          ! npraut_r = MIN(npraut_r,npraut)
          ! endif

          !
          ! (2) ACCRETION OF CLOUD WATER BY RAIN (Pracw): Lin (51)
          !
          erw=1.
          tmp1=pio4*erw*n0_r*av_r*sqrho* &
               gambp3*olambdar**bp3 ! son
          pracw=qlzodt*( 1.-exp(-tmp1*dtb) )
          npracw=tmp1*ncz(iq,k)

          !
          ! (3) EVAPORATION OF RAIN (Prevp): Lin (52)
          !     Prevp is negative value
          !
          !     Sw=qvoqsw : saturation ratio
          !
          tmpa=rvapor*xka*tem(iq,k)*tem(iq,k)
          tmpb=xlv*xlv*rho(iq,k)*qswz*diffwv
          tmpc=tmpa*qswz*diffwv
          tmpd=min(0.0,(qvoqswz-0.99)*qswz*odtb)

          abr=2.0*pi*(qvoqswz-0.99)*tmpc/(tmpa+tmpb)
          tmp1=av_r*sqrho*olambdar**bp5/visc !son
          tmp2=abr*n0_r*( vf1r*olambdar*olambdar+  &
               vf2r*schmidt**0.33334*gambp5o2*sqrt(tmp1) )
          tmp3=min( 0.0,tmp2 )
          tmp3=max( tmp3,tmpd )
          prevp=max( tmp3,-qrzodt )
          if (qrz(iq,k).gt.0.) then
            nprevp=prevp*nrz(iq,k)/qrz(iq,k)
          else
            nprevp=prevp*xmr
          end if

          !
          !**********************************************************************
          !*****     combine all processes together and avoid negative      *****
          !*****     water substances
          !***********************************************************************
          !
          if ( temcc < 0.) then
            !
            !  combined water vapor depletions
            !
            tmp=psdep + midep
            if ( tmp > qvzodt ) then
              factor=qvzodt/tmp
              psdep=psdep*factor
              midep=midep*factor
            end if
            !
            !  combined cloud water depletions
            !
            tmp=praut+psacw+psfw+pracw
            if ( tmp > qlzodt ) then
              factor=qlzodt/tmp
              praut=praut*factor
              psacw=psacw*factor
              psfw=psfw*factor
              pracw=pracw*factor
            end if
            !
            !  combined cloud ice depletions
            !
            tmp=psaut+psaci+praci+psfi
            if (tmp > qizodt ) then
              factor=qizodt/tmp
              psaut=psaut*factor
              psaci=psaci*factor
              praci=praci*factor
              psfi=psfi*factor
            endif

            !
            !  combined all rain processes
            !
            tmp_r=piacr+psacr-prevp-  & 
                  praut-pracw+pgfr
            if (tmp_r > qrzodt ) then
              factor=qrzodt/tmp_r
              piacr=piacr*factor
              psacr=psacr*factor
              prevp=prevp*factor
              pgfr=pgfr*factor
            endif
            !
            !   combined all snow processes
            !
            tmp_s=-pssub-(psaut+psaci+   &
                   psacw+psfw+pgfr+      &
                   psfi+praci+piacr+     &
                   psdep+psacr-pracs)
            if ( tmp_s > qszodt ) then
              factor=qszodt/tmp_s
              pssub=pssub*factor
              pracs=pracs*factor
            endif

            !
            ! calculate new water substances, thetae, tem, and qvsbar
            !

            pvapor=-pssub-psdep-prevp-midep
            qvz(iq,k)=max( qvmin,qvz(iq,k)+dtb*pvapor )
            pclw=-praut-pracw-psacw-psfw
            qlz(iq,k)=max( 0.0,qlz(iq,k)+dtb*pclw )
            pcli=-psaut-psfi-psaci- & 
                      praci+midep
            qiz(iq,k)=max( 0.0,qiz(iq,k)+dtb*pcli )
            tmp_r=piacr+psacr-prevp-praut- &
                      pracw+pgfr-pracs
            prain=-tmp_r
            qrz(iq,k)=max( 0.0,qrz(iq,k)+dtb*prain )
            tmp_s=-pssub-(psaut+psaci+ &
                   psacw+psfw+pgfr+          &
                   psfi+praci+piacr+         &
                   psdep+psacr-pracs)
            psnow=-tmp_s
            qsz(iq,k)=max( 0.0,qsz(iq,k)+dtb*psnow )
            !qschg=qschg+psnow
            qschg=psnow

            tmp=ocp/tothz(iq,k)*xLf*qschg
            theiz = theiz +dtb*tmp
            ! thz(k)=theiz-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
            ! tem(k)=thz(k)*tothz(k)
            ! temcc=tem(k)-273.15
            !==================update temperature=================================================
            temcc=tem(iq,k)-273.15
            lvap = xlv + (2106.0 - 4218.0)*temcc  !Enthalpy of vaporization
            tmp1=(pssub+psdep)*xls*ocp + prevp*lvap*ocp+  &
                 (psfw+pgfr+psacr-pracs)*xlfocp
            !bug fixed 20191126
            tem(iq,k)=tem(iq,k)+tmp1*dtb
            temcc=tem(iq,k)-273.15
            thz(iq,k)=tem(iq,k)/tothz(iq,k)
            !===================================================================
            if ( temcc < -40. ) qswz=qsiz
            qlpqi=qlz(iq,k)+qiz(iq,k)
            if ( qlpqi == 0. ) then
              qvsbar=qsiz
            else
              qvsbar=(qiz(iq,k)*qsiz+qlz(iq,k)*qswz)/qlpqi
            endif
            tmp1=-npraut-npracw-npsacw
            ncz(iq,k)=max( 0.0,ncz(iq,k)+dtb*tmp1 )
            tmp1=-npsaut-npsaci-npraci+nidep
            niz(iq,k)=max( 0.0,niz(iq,k)+dtb*tmp1 )
            tmp1=npiacr+npsacr-nprevp-npraut_r+npgfr
            nrz(iq,k)=max( 0.0,nrz(iq,k)-dtb*tmp1 )
            tmp1=-(npsaut+npgfr+   &
                   npraci+npiacr+  &
                   npsdep+npsacr)
            nsz(iq,k)=max( 0.0,nsz(iq,k)-dtb*tmp1 )

          else ! if ( temcc .lt. 0.0) ..else..
            !
            !  combined cloud water depletions
            !
            tmp=praut+psacw+pracw
            if ( tmp > qlzodt ) then
              factor=qlzodt/tmp
              praut=praut*factor
              psacw=psacw*factor
              pracw=pracw*factor
            end if
            !
            !  combined all snow processes
            !
            tmp_s=-(psmlt+psmltevp)
            if (tmp_s > qszodt ) then
              factor=qszodt/tmp_s
              psmlt=psmlt*factor
              psmltevp=psmltevp*factor
            endif
            !
            !  combined all rain processes
            !
            tmp_r=-prevp-(praut+pracw+psacw-psmlt)
            if (tmp_r > qrzodt ) then
              factor=qrzodt/tmp_r
              prevp=prevp*factor
            endif
            !
            !  calculate new water substances and thetae
            !
            pvapor=-psmltevp-prevp
            qvz(iq,k)=max( qvmin,qvz(iq,k)+dtb*pvapor)
            pclw=-praut-pracw-psacw
            qlz(iq,k)=max( 0.0,qlz(iq,k)+dtb*pclw )
            pcli=0.0
            qiz(iq,k)=max( 0.0,qiz(iq,k)+dtb*pcli )
            tmp_r=-prevp-(praut+pracw+psacw-psmlt)
            prain=-tmp_r
            !tmpqrz=qrz(iq,k)
            qrz(iq,k)=max( 0.0,qrz(iq,k)+dtb*prain )
            tmp_s=-(psmlt+psmltevp)
            psnow=-tmp_s
            qsz(iq,k)=max( 0.0,qsz(iq,k)+dtb*psnow )
            qschg=psnow

            tmp=ocp/tothz(iq,k)*xLf*qschg
            theiz = theiz+dtb*tmp
            ! thz(k)=theiz(k)-(xLvocp*qvz(k)-xLfocp*qiz(k))/tothz(k)
            ! tem(k)=thz(k)*tothz(k)
            ! temcc=tem(k)-273.15
            !==================update tmperature=================================================
            temcc=tem(iq,k)-273.15
            lvap = xlv + (2106.0 - 4218.0)*temcc  !Enthalpy of vaporization
            tmp1=psmltevp*xls*ocp + prevp*lvap*ocp+  &
                 psmlt*xlfocp
            !tmp1 =  ! 1. evaporation of rain formed by melting snow ??? (-)
                     ! 2. evaporation of rain (-)
                     ! 3. melting of snow to form rain (+)
            tem(iq,k)=tem(iq,k)+tmp1*dtb
            !bugfix 20191126

            !tem(k)=tem(k)+tmp1*dtb
            temcc=tem(iq,k)-273.15
            thz(iq,k)=tem(iq,k)/tothz(iq,k)

            !===================================================================
            es=1000.*svp1*exp( svp2*temcc/(tem(iq,k)-svp3) )
            qswz=ep2*es/max(prez(iq,k)-es,0.1)
            qsiz=qswz
            qvsbar=qswz
            ! tmp1=-(npraut+npsacw+npracw)
            ncz(iq,k)=max( 0.0,ncz(iq,k)+dtb*tmp1)
            tmp1=-nprevp-(npraut_r-npsmlt)
            ! tmp1=-nprevp-(nprautr+npracwr+npsacw-npsmltr(k))
            nrz(iq,k)=max(0.0,nrz(iq,k)-dtb*tmp1)
            tmp1=-(npsmlt+npsmltevp)
            nsz(iq,k)=max( 0.0,nsz(iq,k)-dtb*tmp1 )
          end if    !T seperate for source and sink terms

          !rain
          if (qrz(iq,k) > 1.e-8) then
            xlambdar=(pi*rhowater*nrz(iq,k)/qrz(iq,k))**(1./3.)   !zx
            if (xlambdar.lt.lamminr) then
              xlambdar = lamminr
              n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            else if (xlambdar.gt.lammaxr) then
              xlambdar = lammaxr
              n0_r = xlambdar**4*qrz(iq,k)/(pi*rhowater)
              nrz(iq,k) = n0_r/xlambdar
            end if
          end if

          !snow
          if (qsz(iq,k) .gt. 1.e-8) then
            xlambdas=(am_s(iq,k)*gam_ss(iq,k)*     &
            nsz(iq,k)/qsz(iq,k))**(1./bm_s(iq,k))
            if (xlambdas.lt.lammins) then
              xlambdas= lamminS
              n0_s(iq,K) = xlambdas**(bm_s(iq,k)+1.)*      &
                             qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,K) = n0_s(iq,K)/xlambdas
            else if (xlambdas.gt.lammaxs) then
              xlambdas = lammaxs
              n0_s(iq,K) = xlambdas**(bm_s(iq,k)+1.)*      & 
                             qsz(iq,K)/gam_bm_s(iq,k)/am_s(iq,k)
              nsz(iq,K) = n0_s(iq,K)/xlambdas
            end if
          end if

          !cloud ice
          if (qiz(iq,k) >= 1.e-8) then
            lami = max((gam13*500.*pi/6.)*niz(iq,k)/qiz(iq,k),1.e-20)**(1./3) !fixed zdc
            if (lami.lt.lammini) then
              lami= lammini
              n0_i = lami**4./gam13*500.*pi/6.
              niz(iq,K) = n0_i/lami
            else if (lami.gt.lammaxi) then
              lami = lammaxi
              n0_i = lami**4./gam13*500.*pi/6.
              niz(iq,K) = n0_i/lami
            end if
          end if
        
          !cloud water zdc 20220208
          if (qlz(iq,k) >= 1.e-8) then
            lamc = (ncz(iq,k)*rhowater*pi*gg31/(6.*qlz(iq,k)*gg32))**(1./3)
            if (lamc.lt.lammini) then
              lamc= lammini
              tmp1 = lamc**(mu_c+4.)
              tmp2 = 6.*qlz(iq,k)/(pi*rhowater*gg31)
              n0_c= tmp1*tmp2
              ncz(iq,k) = n0_c/lamc
            else if (lamc.gt.lammaxi) then
              lamc= lammaxi
              tmp1 = lamc**(mu_c+4.)
              tmp2 = 6.*qlz(iq,k)/(pi*rhowater*gg31)
              n0_c= tmp1*tmp2
              ncz(iq,k) = n0_c/lamc
            end if
          end if

        
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
          if( qvz(iq,k)+qlz(iq,k)+qiz(iq,k) < rsat*qvsbar ) then ! goto 1800

            !
            !   unsaturated
            !
            qvz(iq,k)=qvz(iq,k)+qlz(iq,k)+qiz(iq,k)
            qlz(iq,k)=0.
            qiz(iq,k)=0.

            thz(iq,k)=theiz-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)

            tem(iq,k)=thz(iq,k)*tothz(iq,k)
            temcc=tem(iq,k)-273.15

          else
            !
            !   saturated
            !
            pladj=qlz(iq,k)
            piadj=qiz(iq,k)
            !

            CALL satadj(qvz(iq,k), qlz(iq,k), qiz(iq,k), prez(iq,k), &
                        theiz, thz(iq,k), tothz(iq,k),               &
                        xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

            pladj=odtb*(qlz(iq,k)-pladj)
            piadj=odtb*(qiz(iq,k)-piadj)
            pclw=pclw+pladj
            pcli=pcli+piadj
            pvapor=pvapor-( pladj+piadj )
            thz(iq,k)=theiz-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
            tem(iq,k)=thz(iq,k)*tothz(iq,k)
            temcc=tem(iq,k)-273.15

            es=1000.*svp1*exp( svp2*temcc/(tem(iq,k)-svp3) )
            qswz=ep2*es/max(prez(iq,k)-es,0.1)
            if (tem(iq,k) < 273.15 ) then
              es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
              qsiz=ep2*es/max(prez(iq,k)-es,0.1)
              if (temcc < -40.) qswz=qsiz
            else
              qsiz=qswz
            endif
            qlpqi=qlz(iq,k)+qiz(iq,k)
            if ( qlpqi .eq. 0. ) then
              qvsbar=qsiz
            else
              qvsbar=( qiz(iq,k)*qsiz+qlz(iq,k)*qswz )/qlpqi
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
              if(temcc .lt. -40.) then
                pihom=qlz(iq,k)*odtb
                nihom=ncz(iq,k)*odtb
              end if
              !
              ! (2)  MELTING OF ICE CRYSTAL WHEN T> 0 C (Pimlt)
              !
              if(temcc .gt. 0.) then
                pimlt=qiz(iq,k)*odtb
                nimlt=niz(iq,k)*odtb
              end if
              !
              ! (3) PRODUCTION OF CLOUD ICE BY BERGERON PROCESS (Pidw): Hsie (p957)
              !     this process only considered when -31 C < T < 0 C
              !
              if(temcc < 0. .and. temcc > -31.) then
                !!
                !!   parama1 and parama2 functions must be user supplied
                !!
                ! aa1=parama1(temcc)
                ! aa2=parama2(temcc)
                i1=int(-temcc)+1
                i1p1=i1+1
                ratio=-(temcc)-real(i1-1)
                aa1=a1(i1)+ratio*( a1(i1p1)-a1(i1) )
                aa2=a2(i1)+ratio*( a2(i1p1)-a2(i1) )
              
                !! change unit from cgs to mks
                aa1=aa1*0.001**(1.0-aa2)
                xnin=xni0*exp(-bni*temcc)
                pidw=xnin/rho(iq,k)*(aa1*xmnin**aa2)
              end if

              pcli=pcli+pihom-pimlt+pidw
              pclw=pclw-pihom+pimlt-pidw
              qlz(iq,k)=max( 0.,qlz(iq,k)+dtb*(-pihom+pimlt-pidw) )
              qiz(iq,k)=max( 0.,qiz(iq,k)+dtb*(pihom-pimlt+pidw) )

              ncz(iq,k)=max( 0.,ncz(iq,k)+dtb*(-nihom+nimlt) )
              niz(iq,k)=max( 0.,niz(iq,k)+dtb*( nihom-nimlt) )

              CALL satadj(qvz(iq,k), qlz(iq,k), qiz(iq,k), prez(iq,k), &
                      theiz, thz(iq,k), tothz(iq,k),                   &
                      xLvocp, xLfocp, episp0k, EP2,SVP1,SVP2,SVP3,SVPT0)

              thz(iq,k)=theiz-(xLvocp*qvz(iq,k)-xLfocp*qiz(iq,k))/tothz(iq,k)
              tem(iq,k)=thz(iq,k)*tothz(iq,k)
              temcc=tem(iq,k)-273.15
              es=1000.*svp1*exp( svp2*temcc/(tem(iq,k)-svp3) )
              qswz=ep2*es/max(prez(iq,k)-es,0.1)

              if (tem(iq,k) < 273.15 ) then
                es=1000.*svp1*exp( 21.8745584*(tem(iq,k)-273.16)/(tem(iq,k)-7.66) )
                qsiz=ep2*es/max(prez(iq,k)-es,0.1)
                if (temcc < -40.0) qswz=qsiz
              else
                qsiz=qswz
              endif
              qlpqi=qlz(iq,k)+qiz(iq,k)

              if ( qlpqi == 0. ) then
                qvsbar=qsiz
              else
                qvsbar=( qiz(iq,k)*qsiz+qlz(iq,k)*qswz )/qlpqi
              endif
            end if ! 1800  continue
          end if ! 1800  continue

        end if   ! mask
    
#ifndef GPU      
        ! save all process rate for understanding cloud microphysics
        zpsaut(iq,k)   = psaut    ! ice crystal aggregation to snow
        zpsfw(iq,k)    = psfw     ! BERGERON process to transfer cloud water to snow
        zpsfi(iq,k)    = psfi     ! BERGERON process to transfer cloud ice to snow
        zpraci(iq,k)   = praci    ! cloud ice accretion by rain
        zpiacr(iq,k)   = piacr    ! rain accretion by cloud ice
        zpsaci(iq,k)   = psaci    ! ice crystal accretion by snow
        zpsacw(iq,k)   = psacw    ! accretion of cloud water by snow
        zpsdep(iq,k)   = psdep    ! deposition of snow
        zpssub(iq,k)   = pssub    ! sublimation of snow (T<0)
        zpracs(iq,k)   = pracs    ! accretion of snow by rain
        zpsacr(iq,k)   = psacr    ! accretion of rain by snow
        zpsmlt(iq,k)   = psmlt    ! melting of snow
        zpsmltevp(iq,k)= psmltevp ! evaporation of melting snow (T>0)
        zprain(iq,k)   = prain    ! sum all process for rain
        zpraut(iq,k)   = praut    ! autoconversion of rain
        zpracw(iq,k)   = pracw    ! accretion of cloud water by rain
        zprevp(iq,k)   = prevp    ! evaporation of rain
        zpgfr(iq,k)    = pgfr     ! feezing of rain to form graupel (added to PI)
        zpvapor(iq,k)  = pvapor   ! sum all process for water vapor to determine qvz
        zpclw(iq,k)    = pclw     ! sum all process for cloud liquid to determine qlz
        zpladj(iq,k)   = pladj    ! saturation adjustment for ql
        zpcli(iq,k)    = pcli     ! sum all process for cloud ice to determine qiz
        zpimlt(iq,k)   = pimlt    ! melting of ice crystal >0.
        zpihom(iq,k)   = pihom    ! homogeneous nucleation <-40
        zpidw(iq,k)    = pidw     ! production of cloud ice by BERGERON process
        zpiadj(iq,k)   = piadj    ! saturation adjustment for qi
        zpsnow(iq,k)   = psnow
        zqschg(iq,k)   = qschg
#endif      

        ! save process rate for aerisol scheme
        fluxm(iq,k) = -1.*psmlt*dzw(iq,k)*rho(iq,k)         ! - ice melting flux in layer k (kg/m2/s)
        fluxf(iq,k) = pgfr*dzw(iq,k)*rho(iq,k)              ! - liquid freezing flux in layer k (kg/m2/s)
        fevap(iq,k) = -1.*prevp*dzw(iq,k)*rho(iq,k)         ! - evaporation of rainfall flux (kg/m2/s)
        fsubl(iq,k) = -1.*pssub*dzw(iq,k)*rho(iq,k)         ! - sublimation of snow, ice and graupel flux (kg/m2/s)
        fauto(iq,k) = praut*dzw(iq,k)*rho(iq,k)             ! - autoconversion flux for rainfall (kg/m2/s)
        fcoll(iq,k) = pracw*dzw(iq,k)*rho(iq,k)             ! - collection of cloud liquid water by rain (kg/m2/s)
        faccr(iq,k) = psacw*dzw(iq,k)*rho(iq,k)          &  !
                    + psfw*dzw(iq,k)*rho(iq,k)              ! - accretion of cloud liq water by snow,ice and graupel (kg/m2/s)

        !
        !***********************************************************************
        !**********    integrate the productions of rain and snow     **********
        !***********************************************************************
        !    

        if ( qvz(iq,k) < qvmin ) then 
          qvz(iq,k)=max( qvmin,qvz(iq,k)+qlz(iq,k)+qiz(iq,k) )
          qlz(iq,k)=0.
          qiz(iq,k)=0.
          ncz(iq,k)=0.
          niz(iq,k)=0.
        end if
        niz(iq,k) = min( niz(iq,k), 0.3E6/rho(iq,k) )
        ncz(iq,k) = min( ncz(iq,k), 250000.E6/rho(iq,k) )
        ncz(iq,k) = max( ncz(iq,k), 0.01E6/rho(iq,k) )

        ! CALCULATE EFFECTIVE RADIUS zdc 20220208
        if (qiz(iq,k) > 1.e-8 .and. lami>0. ) then
          EFFI1D(iq,k) = 3./LAMI/2.
        else
          EFFI1D(iq,k) = 25.E-6
        end if

        if (qsz(iq,k) > 1.e-8) then
          EFFS1D(iq,k) = 3./xlambdas/2.
        else
          EFFS1D(iq,k) = 25.E-6
        end if

        if (qrz(iq,k) > 1.e-8) then 
           EFFR1D(iq,k) = 3./xlambdar/2.
        else
           EFFR1D(iq,k) = 25.E-6
        end if
    
        if (qlz(iq,k)>1.e-8 .and. lamc>0.) then
          !EFFC1D(iq,k) = GAMMA(mu_c+4.)/GAMMA(mu_c+3.)/LAMC/2.
          EFFC1D(iq,k) = 0.5*(mu_c+3.)/LAMC
        else
          EFFC1D(iq,k) = 25.E-6
        end if

      end do     ! iq  
    END DO       ! k
    !$acc end parallel loop
  
  end do ! n = 1,njumps
  
  !$acc update self(Ri)
  !$acc update self(fluxr,fluxi,fluxs,fluxm,fluxf,fevap,fsubl,fauto,fcoll,faccr)
  !$acc update self(effc1d,effr1d,effs1d,effi1d)
  !$acc end data

END SUBROUTINE clphy1d_ylin



!---------------------------------------------------------------------
!                         SATURATED ADJUSTMENT
!---------------------------------------------------------------------
PURE SUBROUTINE satadj(qvz, qlz, qiz, prez, theiz, thz, tothz,      &
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
    do n = 1,20 ! MJT notes, usually 2-3 iterations are required
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
  real                :: temp2, temp3, temp4
  real                :: ans
    
  TEMP = X
  diff = max(int(temp-2.), 0)
  if ( temp-real(diff) > 2. ) diff = diff + 1

  ! Original method
  PF = 1.
  do J = 1,diff
    TEMP = TEMP - 1.
    PF = PF*TEMP
  end do

  ! Alternative method
  !temp = temp - real(diff)
  !pf = gamma( x ) / gamma( temp )

  TEMP = TEMP - 1.
  
  ! method suggested by Mark Dwyer and Lindsay Brebber
  TEMP2 = TEMP*TEMP
  TEMP3 = TEMP2*TEMP
  TEMP4 = TEMP2*TEMP2
  
  G1TO2=1. + B(1)*TEMP + B(2)*TEMP2 + B(3)*TEMP3 + B(4)*TEMP4 &
           + B(5)*TEMP4*TEMP + B(6)*TEMP3*TEMP3 + B(7)*TEMP4*TEMP3 + B(8)*TEMP4*TEMP4
  ans=PF*G1TO2
  
END FUNCTION ggamma

!----------------------------------------------------------------

END MODULE module_mp_sbu_ylin
