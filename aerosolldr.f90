
! This is a prognostic aerosol model for CCAM based on the LDR scheme used in Mk3.6 (Rotstayn and Lohmann 2002)
    
! Plan to replace diagnosed sea-salt with prognostic version based on CTM.  Eventually employ GLOMAP for modes.

module aerosolldr

implicit none

private
public aldrcalc,aldrinit,aldrend,aldrloademiss,aldrloaderod,aldrloadoxidant,cldrop,convscav
public xtg,xtgsav,xtosav,naero,ssn,ticeu
public itracdu,ndust,dustdd,dustwd,duste,Ch_dust
public itracbc,bce,bcdd,bcwd
public itracoc,oce,ocdd,ocwd
public itracso2,dmse,dmsso2o,so2e,so2so4o,so2dd,so2wd,so4e,so4dd,so4wd

integer, save :: ifull,kl
integer, save :: jk2,jk3,jk4,jk5,jk6,jk8,jk9           ! levels for injection
real, dimension(:,:,:), allocatable, save :: xtg       ! prognostic aerosols (see indexing below)
real, dimension(:,:,:), allocatable, save :: xtgsav    ! save for mass conservation in semi-Lagrangian models
real, dimension(:,:,:), allocatable, save :: xtosav    ! aerosol mixing ratio outside convective cloud
real, dimension(:,:,:), allocatable, save :: ssn       ! diagnostic sea salt concentration
real, dimension(:,:), allocatable, save :: erod        ! sand, clay and silt fraction that can erode
real, dimension(:,:), allocatable, save :: emissfield  ! non-volcanic emissions
real, dimension(:,:), allocatable, save :: zoxidant    ! oxidant fields
real, dimension(:), allocatable, save :: vso2          ! volcanic emissions
real, dimension(:), allocatable, save :: duste         ! Diagnostic - dust emissions
real, dimension(:), allocatable, save :: dustdd        ! Diagnostic - dust dry deposition
real, dimension(:), allocatable, save :: dustwd        ! Diagnostic - dust wet deposition
real, dimension(:), allocatable, save :: bce           ! Diagnostic - black carbon emissions
real, dimension(:), allocatable, save :: bcdd          ! Diagnostic - black carbon dry deposition
real, dimension(:), allocatable, save :: bcwd          ! Diagnostic - black carbon wet deposition
real, dimension(:), allocatable, save :: oce           ! Diagnostic - organic carbon emissions
real, dimension(:), allocatable, save :: ocdd          ! Diagnostic - organic carbon dry deposition
real, dimension(:), allocatable, save :: ocwd          ! Diagnostic - organic carbon wet deposition
real, dimension(:), allocatable, save :: dmse          ! Diagnostic - DMS emissions
real, dimension(:), allocatable, save :: dmsso2o       ! Diagnostic - DMS->so2 oxidation
real, dimension(:), allocatable, save :: so2e          ! Diagnostic - so2 emissions
real, dimension(:), allocatable, save :: so2so4o       ! Diagnostic - so2->so4 oxidation
real, dimension(:), allocatable, save :: so2dd         ! Diagnostic - so2 dry deposition
real, dimension(:), allocatable, save :: so2wd         ! Diagnostic - so2 wet deposition
real, dimension(:), allocatable, save :: so4e          ! Diagnostic - so4 emissions
real, dimension(:), allocatable, save :: so4dd         ! Diagnostic - so4 dry deposition
real, dimension(:), allocatable, save :: so4wd         ! Diagnostic - so4 wet deposition

! parameters
integer, parameter :: nsulf = 3
integer, parameter :: ncarb = 4
integer, parameter :: ndust = 4
integer, parameter :: naero = nsulf+ncarb+ndust ! Tracers: DMS, SO2, SO4, BCO, BCI, OCO, OCI, DUST(4)
integer, parameter :: itracso2 = 2              ! Index for SO2 tracer
integer, parameter :: itracbc = nsulf+1         ! Index for BC    "    (hydrophobicm, hydrophillic)
integer, parameter :: itracoc = nsulf+3         ! Index for OC    "    (hydrophobicm, hydrophillic)
integer, parameter :: itracdu = nsulf+ncarb+1   ! Index for dust  "
integer, parameter :: ndcls = 3                 ! No. of dust emission classes (sand, silt, clay)

integer, parameter :: enhanceu10 = 0            ! Modify 10m wind speed (0=none, 1=quadrature, 2=linear)

! physical constants
real, parameter :: grav      = 9.80616          ! Gravitation constant
real, parameter :: rdry      = 287.04           ! Specific gas const for dry air
real, parameter :: cp        = 1004.64          ! Heat capacity of air
real, parameter :: hl        = 2.5104e6         ! Latent heat of vaporisation
real, parameter :: vkar      = 0.4              ! von Karman constant
real, parameter :: rhos      = 100.             ! Assumed density of snow in kg/m^3

! emission constants
real, save :: Ch_dust        = 1.e-9            ! Transfer coeff for type natural source (kg*s2/m5)

! scavenging constants
real, parameter :: ticeu     = 263.16           ! Temperature for freezing in convective updraft

! Following array determines for which tracers the XTWETDEP routine is called.
! We now set it to .true. for BCO and OCO, since they do experience below-cloud scavenging.
! The efficiency for in-cloud scavenging is still set to zero for these.
logical, dimension(naero), parameter :: lwetdep = (/.false.,.true.,.true.,       & ! DMS, SO2, SO4
                                                 .true.,.true.,.true.,.true.,    & ! BCO, BCI, OCO, OCI
                                                 .true.,.true.,.true.,.true. /)    ! Dust
real, dimension(ndust), parameter :: dustden = (/ 2500., 2650., 2650., 2650. /)    ! Density of dust (kg/m3)
                                                                                   ! (Clay, small silt, small slit, small silt)
real, dimension(ndust), parameter :: dustreff = (/ 0.73e-6,1.4e-6,2.4e-6,4.5e-6 /) ! Main effective radius (m)
                                                                                   ! (Clay, small silt, small slit, small silt)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation

subroutine aldrinit(ifin,iextra,klin,sig)

implicit none

integer, intent(in) :: ifin,iextra,klin
integer pos(1)
real, dimension(klin), intent(in) :: sig

ifull=ifin
kl=klin
allocate(xtg(ifull+iextra,kl,naero),xtgsav(ifull,kl,naero))
allocate(xtosav(ifull,kl,naero),vso2(ifull))
allocate(emissfield(ifull,15),ssn(ifull,kl,2))
allocate(zoxidant(ifull,4*kl),erod(ifull,ndcls))
allocate(duste(ifull),dustdd(ifull),dustwd(ifull))
allocate(bce(ifull),bcdd(ifull),bcwd(ifull))
allocate(oce(ifull),ocdd(ifull),ocwd(ifull))
allocate(dmse(ifull),dmsso2o(ifull))
allocate(so2e(ifull),so2so4o(ifull),so2dd(ifull),so2wd(ifull))
allocate(so4e(ifull),so4dd(ifull),so4wd(ifull))

xtg=0.
xtgsav=0.
xtosav=0.
vso2=0.
emissfield=0.
ssn=0.
zoxidant=0.
erod=0.
duste=0.
dustdd=0.
dustwd=0.
bce=0.
bcdd=0.
bcwd=0.
oce=0.
ocdd=0.
ocwd=0.
dmse=0.
dmsso2o=0.
so2e=0.
so2so4o=0.
so2dd=0.
so2wd=0.
so4e=0.
so4dd=0.
so4wd=0.

! MJT - define injection levels

! scale to CSIRO9 18 levels
! jk2 is top of level=1, bottom of level=2
jk2=2
! jk3 is top of level=2, bottom of level=3
pos=maxloc(sig,sig<=0.965) ! 300m
jk3=pos(1)
! jk4 is top of level=3, bottom of level=4
pos=maxloc(sig,sig<=0.925) ! 600m
jk4=pos(1)
! jk5 is top of level=4, bottom of level=5
pos=maxloc(sig,sig<=0.870) ! 1,150m
jk5=pos(1)
! jk6 is top of level=5, bottom of level=6
pos=maxloc(sig,sig<=0.800) ! 1,800m
jk6=pos(1)
! jk8 is top of level=7, bottom of level=8
pos=maxloc(sig,sig<=0.650) ! 3,500m
jk8=pos(1)
! jk9 is top of level=8, bottom of level=9
pos=maxloc(sig,sig<=0.500) ! 5,500m
jk9=pos(1)

return
end subroutine aldrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End

subroutine aldrend

implicit none

deallocate(xtg,xtgsav,xtosav)
deallocate(vso2)
deallocate(emissfield)
deallocate(ssn)
deallocate(zoxidant,erod)
deallocate(duste,dustdd,dustwd)
deallocate(bce,bcdd,bcwd)
deallocate(oce,ocdd,ocwd)
deallocate(dmse,dmsso2o)
deallocate(so2e,so2so4o,so2dd,so2wd)
deallocate(so4e,so4dd,so4wd)

return
end subroutine aldrend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load emission arrays

subroutine aldrloademiss(index,aa)

implicit none

integer, intent(in) :: index
real, dimension(ifull), intent(in) :: aa

if (index<16) then
  emissfield(:,index)=aa ! Then follow SO2, BC and OC from anthro (a) and biomass-burning (b) levels 1 and 2
elseif (index==16) then
  vso2(:)=aa        ! volcanic
else
  write(6,*) "ERROR: index out-of-range for aldrloademiss"
  stop
end if

return
end subroutine aldrloademiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load oxidant arrays

subroutine aldrloadoxidant(index,aa)

implicit none

integer, intent(in) :: index
real, dimension(ifull), intent(in) :: aa

! First four are 3d oxidant fields (oh, h2o2, o3, no2)
zoxidant(:,index)=aa

return
end subroutine aldrloadoxidant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load soil data

subroutine aldrloaderod(inda,aa)

implicit none

integer, intent(in) :: inda
real, dimension(ifull), intent(in) :: aa

! EROD is the soil fraction of Sand (inda=1), Silt (inda=2) and Clay (inda=3) that can erode
erod(:,inda)=aa

return
end subroutine aldrloaderod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine

subroutine aldrcalc(dt,sig,sigh,dsig,zz,dz,fwet,wg,pblh,prf,ts,ttg,condc,snowd,sg,fg,eg,v10m,        &
                    ustar,zo,land,sicef,tsigmf,qvg,qlg,qfg,cfrac,clcon,cldcon,pccw,rhoa,vt,ppfprec,  &
                    ppfmelt,ppfsnow,ppfconv,ppfevap,ppfsubl,pplambs,ppmrate,ppmaccr,ppfstay,ppqfsed, &
                    pprscav,zdayfac,kbsav)

implicit none

integer, dimension(ifull), intent(in) :: kbsav ! Bottom of convective cloud
real, intent(in) :: dt                         ! Time step
real, dimension(kl), intent(in) :: sig         ! Sigma levels
real, dimension(kl), intent(in) :: dsig        ! Sigma level width
real, dimension(kl+1), intent(in) :: sigh      ! Sigma half levels
real, dimension(ifull), intent(in) :: fwet     ! Fraction of water on leaf
real, dimension(ifull), intent(in) :: wg       ! Soil moisture fraction of field capacity
real, dimension(ifull), intent(in) :: prf      ! Surface pressure
real, dimension(ifull), intent(in) :: ts       ! Surface temperture
real, dimension(ifull), intent(in) :: pblh     ! Boundary layer height
real, dimension(ifull), intent(in) :: v10m     ! 10m wind speed
real, dimension(ifull), intent(in) :: condc    ! Convective rainfall
real, dimension(ifull), intent(in) :: snowd    ! Snow depth
real, dimension(ifull), intent(in) :: sg       ! Downwelling shortwave radiation
real, dimension(ifull), intent(in) :: fg       ! Sensible heat flux
real, dimension(ifull), intent(in) :: eg       ! Latent heat flux
real, dimension(ifull), intent(in) :: ustar    ! Friction velocity
real, dimension(ifull), intent(in) :: zo       ! Roughness length
real, dimension(ifull), intent(in) :: sicef    ! Sea-ice fraction
real, dimension(ifull), intent(in) :: tsigmf   ! Vegetation fraction
real, dimension(ifull), intent(in) :: vt       ! transfer velocity
real, dimension(ifull), intent(in) :: zdayfac  ! scale factor for day length
real, dimension(ifull,kl), intent(in) :: zz    ! Height of vertical level (meters)
real, dimension(ifull,kl), intent(in) :: dz
real, dimension(:,:), intent(in) :: ttg        ! Air temperature
real, dimension(:,:), intent(in) :: qvg        ! liquid water mixing ratio
real, dimension(:,:), intent(in) :: qlg        ! liquid water mixing ratio
real, dimension(:,:), intent(in) :: qfg        ! frozen water mixing ratio
real, dimension(ifull,kl), intent(in) :: cfrac ! cloud fraction
real, dimension(ifull,kl), intent(in) :: clcon ! convective cloud fraction
real, dimension(ifull), intent(in) :: cldcon   ! convective rain fraction
real, dimension(ifull,kl), intent(in) :: pccw
real, dimension(ifull,kl), intent(in) :: rhoa  ! density of air
real, dimension(ifull,kl), intent(inout) :: ppfconv                      ! from LDR prog cloud
real, dimension(ifull,kl), intent(in) :: ppfprec,ppfmelt,ppfsnow         ! from LDR prog cloud
real, dimension(ifull,kl), intent(in) :: ppfevap,ppfsubl,pplambs,ppmrate ! from LDR prog cloud
real, dimension(ifull,kl), intent(in) :: ppmaccr,ppfstay,ppqfsed,pprscav ! from LDR prog cloud
logical, dimension(ifull), intent(in) :: land  ! land/sea mask (t=land)
real, dimension(ifull,naero) :: conwd          ! Diagnostic only: Convective wet deposition
real, dimension(ifull,naero) :: xtem
real, dimension(ifull,kl,naero) :: xte,xtu,xtm1
real, dimension(ifull,kl+1) :: aphp1
real, dimension(ifull,kl) :: pclcon
real, dimension(ifull,kl) :: prhop1,ptp1,pfevap,pfsubl,plambs
real, dimension(ifull,kl) :: pclcover,pcfcover,pmlwc,pmiwc
real, dimension(ifull) :: bbem
real, dimension(ifull) :: so2oh,so2h2,so2o3,dmsoh,dmsn3
real, dimension(ifull) :: cgssnowd
real, dimension(ifull) :: veff,vefn
real, dimension(ifull) :: cstrat,qtot
real, dimension(ifull) :: rrate,Wstar3,Vgust_free,Vgust_deep
real, dimension(ifull) :: v10n,thetav
real, parameter :: beta = 0.65
integer nt,k

conwd=0.
cgssnowd=1.E-3*snowd

! Calculate sub-grid Vgust
v10n=ustar*log(10./zo)/vkar ! neutral wind speed
! Mesoscale enhancement follows Redelsperger et al. (2000), J. Climate 13, 402-421.
! Equation numbers follow Fairall et al. 1996, JGR 101, 3747-3764.

! Calculate convective scaling velocity (Eq.17) and gustiness velocity (Eq.16)
thetav = ttg(1:ifull,1)*(1.+0.61*qvg(1:ifull,1))
Wstar3 = max(0.,(grav*pblh/thetav)*(fg/cp+0.61*ttg(1:ifull,1)*eg/hl)/rhoa(:,1))
Vgust_free = beta*Wstar3**(1./3.)
! Calculate the Redelsperger-based Vgust_deep if deep convection is present.
! Note that Redelspreger gives two other parameterizations, based on
! the updraft or downdraught mass fluxes respectively.
rrate = 8640.*condc/dt   !Rainfall rate in cm/day
Vgust_deep = (19.8*rrate*rrate/(1.5+rrate+rrate*rrate))**0.4
! Calculate effective 10m wind (Eq. 15)
! These can plausibly be added in quadrature, or linearly, the latter giving a much larger
! effect (Lunt & Valdes, JGR, 2002).
select case(enhanceu10)
  case(0)
    veff = v10m
    vefn = v10n
  case(1)
    veff = sqrt( v10m*v10m + Vgust_free*Vgust_free + Vgust_deep*Vgust_deep )
    vefn = sqrt( v10n*v10n + Vgust_free*Vgust_free + Vgust_deep*Vgust_deep )
  case(2)
    veff = v10m + Vgust_free + Vgust_deep
    vefn = v10n + Vgust_free + Vgust_deep
  case DEFAULT
    write(6,*) "Unknown option for enhanceu10 ",enhanceu10
    stop
end select

! Emission and dry deposition (sulfur cycle and carbonaceous aerosols)
do k=1,kl+1
  aphp1(:,k)=prf(:)*sigh(k)
enddo
! MJT notes - replace vefn with v10n for regional model
call xtemiss(dt, rhoa(:,1), ts, sicef, vefn, aphp1,                       & !Inputs
             land, tsigmf, cgssnowd, fwet, wg,                            & !Inputs
             xte, xtem, bbem)                                               !Outputs
xtg(1:ifull,:,:)=max(xtg(1:ifull,:,:)+xte(:,:,:)*dt,0.)

! Emission and dry deposition of dust
do k=1,kl
  aphp1(:,k)=prf(:)*sig(k)*0.01 ! hPa
end do
! Calculate the settling of large dust particles
call dsettling(dt,rhoa,ttg,dz,aphp1(:,1:kl))
! Calculate dust emission and turbulent dry deposition at the surface
call dustem(dt,rhoa(:,1),wg,veff,dz(:,1),vt,snowd,land)

! Decay of hydrophobic black and organic carbon into hydrophilic forms
call xtsink(dt,xte)
xtg(1:ifull,:,:)=max(xtg(1:ifull,:,:)+xte(:,:,:)*dt,0.)

! Compute diagnostic sea salt aerosol
call seasalt(land,sicef,zz,pblh,veff)

! Aerosol chemistry and wet deposition
! Need to invert vertical levels for ECHAM code... Don't you hate that?
do nt=1,naero
  do k=1,kl
    xtm1(:,kl+1-k,nt)=xtg(1:ifull,k,nt)
    ! Convert from aerosol concentration outside convective cloud (used by CCAM)
    ! to aerosol concentration inside convective cloud
    xtu(:,kl+1-k,nt)=max(xtg(1:ifull,k,nt)-(1.-clcon(:,k))*xtosav(:,k,nt),0.)/max(clcon(:,k),1.E-8)
  end do
end do
do k=1,kl
  aphp1(:,kl+1-k) =-prf(:)*dsig(k)                            ! delta pressure
  prhop1(:,kl+1-k)=rhoa(:,k)                                  ! air density
  ptp1(:,kl+1-k)  =ttg(1:ifull,k)                             ! air temperature
  pclcon(:,kl+1-k)=clcon(:,k)                                 ! convective cloud fraction
  qtot=qlg(1:ifull,k)+qfg(1:ifull,k)                          ! total liquid and ice mixing ratio
  cstrat(:)=max(min(cfrac(:,k)-clcon(:,k),1.),0.)             ! strat cloud fraction (i.e., ccov from leoncld.f)
  pclcover(:,kl+1-k)=cstrat(:)*qlg(1:ifull,k)/max(qtot,1.E-8) ! Liquid-cloud fraction
  pcfcover(:,kl+1-k)=cstrat(:)*qfg(1:ifull,k)/max(qtot,1.E-8) ! Ice-cloud fraction
  pmlwc(:,kl+1-k)=qlg(1:ifull,k)
  pmiwc(:,kl+1-k)=qfg(1:ifull,k)
  where (k<=kbsav)
    ppfconv(:,kl+1-k)=condc(:)/dt
  elsewhere
    ppfconv(:,kl+1-k)=0.
  end where
end do
call xtchemie (1, dt, zdayfac, aphp1(:,1:kl), ppmrate, ppfprec,                 & !Inputs
               pclcover, pmlwc, prhop1, ptp1, sg, xtm1, ppfevap,                & !Inputs
               ppfsnow,ppfsubl,pcfcover,pmiwc,ppmaccr,ppfmelt,ppfstay,ppqfsed,  & !Inputs
               pplambs,pprscav,pclcon,cldcon,pccw,ppfconv,xtu,                  & !Input
               conwd,                                                           & !In and Out
               xte, so2oh, so2h2, so2o3, dmsoh, dmsn3)                            !Output
do nt=1,naero
  do k=1,kl
    xtg(1:ifull,k,nt)=max(xtg(1:ifull,k,nt)+xte(:,kl+1-k,nt)*dt,0.)
  end do
enddo
dmsso2o=dmsoh+dmsn3        ! oxidation of DMS to SO2
so2so4o=so2oh+so2h2+so2o3  ! oxidation of SO2 to SO4

return
end subroutine aldrcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt emiss

SUBROUTINE XTEMISS(ztmst, P1MXTM1, TSM1M, SEAICEM, G3X01, APHP1,                 & !Inputs
                   LOLAND, PFOREST, PSNOW, fwet, WSM1M,                          & !Inputs
                   XTE, PXTEMS, bbem)                                              !Outputs
!
!    THIS ROUTINE CALCULATES THE LOWER BOUNDARY CONDITIONS
!    FOR VDIFF DEPENDING ON THE SURFACE EMISSION AND THE
!    DRY DEPOSITION FLUX.
!
!    JOHANN FEICHTER          UNI/HAMBURG         08/91
!    MODIFIED  U. SCHLESE    DKRZ-HAMBURG        JAN-95
!    Adapted for CSIRO GCM by Leon Rotstayn, 12/99
!    Adapted for CCAM by Marcus Thatcher 2012
!
!    PURPOSE
!   ---------
!    THE LOWER BOUNDARY CONDITION FOR CALCULATING THE
!    TURBULENT EXCHANGE IN THE BOUNDARY LAYER IS
!    DETERMINED BY THE EMISSION AND THE DRY DEPOSITION FLUX.
!

implicit none

! Argument list
REAL, intent(in) :: ztmst                           !Timestep [s]
REAL, dimension(ifull), intent(in) :: P1MXTM1       !Density of air in surface layer
real, dimension(ifull), intent(in) :: TSM1M         !Surface temp
real, dimension(ifull), intent(in) :: SEAICEM       !Sea-ice fraction
real, dimension(ifull), intent(in) :: G3X01         !10m wind (corrected to neutral for Nightingale scheme)
real, dimension(ifull,kl+1), intent(in) :: APHP1    !P at half levels at current timestep
LOGICAL, dimension(ifull), intent(in) :: LOLAND     !Land flag
REAL, dimension(ifull), intent(in) :: PFOREST       !Fractional vegetation cover
REAL, dimension(ifull), intent(in) :: PSNOW         !Snow depth [m]
! Land-surface details needed to specify dry deposition velocity
REAL, dimension(ifull), intent(in) :: fwet          !skin reservoir content of plant [fraction of maximum]
real, dimension(ifull), intent(in) :: WSM1M         !surface wetness [vol fraction for CSIRO GCM, not m]
real, dimension(ifull,kl,naero), intent(out) :: XTE !Tracer tendencies (kg/kg/s)
REAL, dimension(ifull,naero), intent(out) :: PXTEMS !Sfc. flux of tracer passed to vertical mixing [kg/m2/s]
! Some diagnostics
real, dimension(ifull), intent(out) :: bbem

real zdmscon,ZSST,ZZSPEED,VpCO2,VpCO2liss
real wtliss,ScDMS,zVdms,ddt
integer jl,jk,jt,nstep,ii

REAL, dimension(ifull,2) :: ZVDRD
REAL, dimension(ifull) :: ZMAXVDRY,gdp,zdmsemiss,pxtm1new
real, dimension(ifull) :: zhilbco,zhilbcy,zhiloco,zhilocy
real, dimension(ifull) :: dmsdd
real zvolcemi1,zvolcemi2,zvolcemi3
real zvd2ice,zvd4ice,zvd2nof,zvd4nof

!     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT (FROM *SURF*)
real, parameter :: ZSNCRI = 0.025
!     COEFFICIENTS FOR ZVDRD = FUNCTION OF SOIL MOISTURE
real, parameter :: ZVWC2 = (0.8E-2 - 0.2E-2)/(1. - 0.9)
real, parameter :: ZVW02 = ZVWC2-0.8E-2
real, parameter :: ZVWC4 = (0.2E-2 - 0.025E-2)/(1. - 0.9)
real, parameter :: ZVW04 = ZVWC4-0.2E-2
real, parameter :: tmelt = 273.05
real, parameter :: zvolcemi  = 8.             ! ZVOLCEMI  TOTAL EMISSION FROM VOLCANOES IN TG/YR
real, parameter :: ZVDPHOBIC = 0.025E-2       ! dry deposition
real, parameter :: ScCO2     = 600.

! Then follow SO2, BC and OC from anthro (a) and biomass-burning (b) levels 1 and 2
integer, parameter :: iso2a1=1
integer, parameter :: iso2a2=2
integer, parameter :: ibca1 =3
integer, parameter :: ibca2 =4
integer, parameter :: ioca1 =5
integer, parameter :: ioca2 =6
integer, parameter :: iso2b1=7
integer, parameter :: iso2b2=8
integer, parameter :: ibcb1 =9
integer, parameter :: ibcb2 =10
integer, parameter :: iocb1 =11
integer, parameter :: iocb2 =12
integer, parameter :: idmso =13     ! DMS ocean
integer, parameter :: idmst =14     ! DMS terr
integer, parameter :: iocna =15     ! Nat org

! Start code : ----------------------------------------------------------

pxtems(:,:)=0.
xte(:,:,:)=0.

! --------------------------------------------------------------
!
!*     1.   SURFACE EMISSION.
!           ------- --------
!
!   CALCULATE DMS EMISSIONS FOLLOWING LISS+MERLIVAT
!   DMS SEAWATER CONC. FROM KETTLE ET AL.
DO JL=1,ifull
  IF (LOLAND(JL)) THEN
    zdmsemiss(jl)=emissfield(jl,idmst) !kg/m2/s
  ELSE
    ! Reduce the effective DMS concentration more strongly at seaice points, since the flux should
    ! probably depend on wave breaking (e.g., Erickson, JGR, 1993; Woolf, Tellus B, 2005),
    ! which will be much reduced over leads.
    ZDMSCON=EMISSFIELD(JL,idmso)*(1.-SEAICEM(JL))**2
    ZSST=TSM1M(JL)-273.15 ! DegC
    ZSST=min(ZSST, 45.)   ! Even Saltzman Sc formula has trouble over 45 deg C
    !  G3X01:  10-M WINDS
    ZZSPEED=G3X01(JL)
    ! Nightingale (2000) scheme (J. Biogeochem. Cycles, 14, 373-387)
    ! For this scheme, zzspeed is the 10m wind adjusted to neutral stability.
    ! The formula for ScDMS from Saltzman et al (1993) is given by Kettle & Andreae (ref below)
    VpCO2 = 0.222*ZZSPEED**2 + 0.333*ZZSPEED !Nightingale et al
    ! Phase in Liss & Merlivat from 13 to 18 m/s, since Nightingale is doubtful for high windspeeds,
    ! due to limited data.
    VpCO2liss=5.9*ZZSPEED-49.3
    wtliss=dim(min(18.,ZZSPEED),13.)/5. ! dim(x,y) is the difference between x and y if the difference is positive
    VpCO2=wtliss*VpCO2liss+(1.-wtliss)*VpCO2
    ScDMS = 2674. - 147.12*ZSST + 3.726*ZSST**2 - 0.038*ZSST**3 !Sc for DMS (Saltzman et al.)
    if ( ZZSPEED<3.6 ) then
      zVdms = VpCO2 * (ScCO2/ScDMS)**(2./3.)
    else
      zVdms = VpCO2 * sqrt(ScCO2/ScDMS)
    end if
    ZDMSEMISS(jl)=ZDMSCON*ZVDMS*32.064E-11/3600.
    ! NANOMOL/LTR*CM/HOUR --> KG/M**2/SEC
  END IF
end do
gdp=grav/(aphp1(:,1)-aphp1(:,2))
xte(:,1,itracso2-1)=xte(:,1,itracso2-1)+zdmsemiss*gdp

! Other biomass emissions of SO2 are done below (with the non-surface S emissions)
PXTEMS(:,ITRACSO2)  =(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.97
PXTEMS(:,ITRACSO2+1)=(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.03
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
xte(:,1,itracso2)  =xte(:,1,itracso2)  +pxtems(:,itracso2)*gdp
xte(:,1,itracso2+1)=xte(:,1,itracso2+1)+pxtems(:,itracso2+1)*gdp

!Do carbonaceous aerosols
! Inject the low-level fossil-fuel and natural SOA emissions into layer 1
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACOC)  =0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
PXTEMS(:,ITRACOC+1)=0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
xte(:,1,itracbc)  =xte(:,1,itracbc)  +pxtems(:,itracbc)*gdp
xte(:,1,itracbc+1)=xte(:,1,itracbc+1)+pxtems(:,itracbc+1)*gdp
xte(:,1,itracoc)  =xte(:,1,itracoc)  +pxtems(:,itracoc)*gdp
xte(:,1,itracoc+1)=xte(:,1,itracoc+1)+pxtems(:,itracoc+1)*gdp
! Inject the upper-level fossil-fuel emissions into layer 2
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,ioca2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,ioca2)
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
do jk=jk2,jk3-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk3-jk2)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+pxtems(:,itracoc+1)*gdp
end do
! Inject the lower-level biomass emissions into layer 2 (NB: Doesn't include biofuel any more)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb1)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb1)
! Apply these here as a tendency (XTE)
do jk=jk2,jk3-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk3-jk2)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+pxtems(:,itracoc+1)*gdp
end do
! Inject the upper-level biomass emissions into layers 3-5 (30%, 40%, 30%)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb2)
! Apply these here as a tendency (XTE)
do jk=jk3,jk4-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk4-jk3)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.3*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.3*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.3*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.3*pxtems(:,itracoc+1)*gdp
end do
do jk=jk4,jk5-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk5-jk4)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.4*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.4*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.4*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.4*pxtems(:,itracoc+1)*gdp
end do
do jk=jk5,jk6-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk6-jk5)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.3*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.3*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.3*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.3*pxtems(:,itracoc+1)*gdp
end do

!  EMISSION OF ANTHROPOGENIC SO2 IN THE NEXT HIGHER LEVEL PLUS BIOMASS BURNING
do jk=jk2,jk3-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk3-jk2)
  XTE(:,JK,ITRACSO2)  =XTE(:,JK,ITRACSO2)  +0.97*EMISSFIELD(:,iso2a2)*gdp !100% of the "above 100m" SO2 emission
  XTE(:,JK,ITRACSO2+1)=XTE(:,JK,ITRACSO2+1)+0.03*EMISSFIELD(:,iso2a2)*gdp !100% of the "above 100m" SO4 emission
end do
do jk=jk3,jk4-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk4-jk3)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.3*emissfield(:,iso2b2)*gdp
end do
do jk=jk4,jk5-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk5-jk4)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.4*emissfield(:,iso2b2)*gdp
end do
do jk=jk5,jk6-1
  gdp=grav/(aphp1(:,jk)-aphp1(:,jk+1))/real(jk6-jk5)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.3*emissfield(:,iso2b2)*gdp
end do
  
!    VOLCANIC BACKGROUND EMISSIONS 
!
!   3 EMISSION LEVELS: 
!    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO)
!    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)
!    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)
ZVOLCEMI1=ZVOLCEMI*0.36
ZVOLCEMI2=ZVOLCEMI*0.36
ZVOLCEMI3=ZVOLCEMI*0.28
gdp=Grav/(APHP1(:,1)-APHP1(:,2))
XTE(:,1,ITRACSO2)=XTE(:,1,ITRACSO2)+ZVOLCEMI1*vso2*gdp
do jk=jk3,jk4-1
  gdp=Grav/(APHP1(:,jk)-APHP1(:,jk+1))/real(jk4-jk3)
  XTE(:,jk,ITRACSO2)=XTE(:,jk,ITRACSO2)+ZVOLCEMI2*vso2*gdp
end do
do jk=jk8,jk9-1
  gdp=Grav/(APHP1(:,jk)-APHP1(:,jk+1))/real(jk9-jk8)
  XTE(:,jk,ITRACSO2)=XTE(:,jk,ITRACSO2)+ZVOLCEMI3*vso2*gdp
end do
  
!   --------------------------------------------------------------
!
!*      2.    DRY DEPOSITION.
!             --- ----------
ZMAXVDRY(:)=(APHP1(:,1)-APHP1(:,2))/(Grav*P1MXTM1*ZTMST)

!      DRY DEPOSITION OF SO2, SO4
do jl=1,ifull
!     -  SEA -
  IF(.NOT.LOLAND(JL)) THEN
!         - SEA ICE -
!           - MELTING/NOT MELTING SEAICE-
    IF(TSM1M(JL)>=(TMELT-0.1)) THEN
      ZVD2ICE=0.8E-2
      ZVD4ICE=0.2E-2
    ELSE
      ZVD2ICE=0.1E-2
      ZVD4ICE=0.025E-2
    ENDIF
    ZVDRD(JL,1)=(1.-SEAICEM(JL))*0.8E-2+SEAICEM(JL)*ZVD2ICE !So leads agree with ocean
    ZVDRD(JL,2)=(1.-SEAICEM(JL))*0.2E-2+SEAICEM(JL)*ZVD4ICE
  ELSE
!      - LAND -
!        - NON-FOREST AREAS -
!         -  SNOW/NO SNOW -
    IF(PSNOW(JL)>ZSNCRI) THEN
!            - MELTING/NOT MELTING SNOW -
      if(tsm1m(jl)>=tmelt) then !This is a simplification of above line
        ZVD2NOF=0.8E-2
        ZVD4NOF=0.2E-2
      ELSE
        ZVD2NOF=0.1E-2
        ZVD4NOF=0.025E-2
      ENDIF
    ELSE
      IF(TSM1M(JL)<=TMELT) THEN
!           -  FROZEN SOIL -
        ZVD2NOF=0.2E-2
        ZVD4NOF=0.025E-2
      ELSE IF(fwet(JL)>=0.01.OR.WSM1M(JL)==1.) THEN
!               - COMPLETELY WET -
        ZVD2NOF=0.8E-2
        ZVD4NOF=0.2E-2
      ELSE IF(WSM1M(JL)<0.9) THEN
!                  - DRY -          
        ZVD2NOF=0.2E-2
        ZVD4NOF=0.025E-2
      ELSE
!                  - PARTLY WET -
        ZVD2NOF=ZVWC2*WSM1M(JL)-ZVW02
        ZVD4NOF=ZVWC4*WSM1M(JL)-ZVW04
      ENDIF
    ENDIF
    ZVDRD(JL,1)=PFOREST(JL)*0.8E-2+(1.-PFOREST(JL))*ZVD2NOF
    ZVDRD(JL,2)=PFOREST(JL)*0.2E-2+(1.-PFOREST(JL))*ZVD4NOF
  ENDIF
  ! Apply lower and upper bounds.
  ZVDRD(JL,1)=AMIN1(ZVDRD(JL,1),ZMAXVDRY(JL)) !SO2
  ZVDRD(JL,2)=AMIN1(ZVDRD(JL,2),ZMAXVDRY(JL)) !aerosols
end do

! Sulfur emission diagnostic (hard-coded for 3 sulfur variables)
do jk=1,kl
  dmse=dmse+xte(:,jk,1)*(aphp1(:,jk)-aphp1(:,jk+1))/grav            !Above surface
enddo
do jk=1,kl
  so2e=so2e+xte(:,jk,ITRACSO2)*(aphp1(:,jk)-aphp1(:,jk+1))/grav     !Above surface
enddo
do jk=1,kl
  so4e=so4e+xte(:,jk,ITRACSO2+1)*(aphp1(:,jk)-aphp1(:,jk+1))/grav   !Above surface
enddo

! Assume that BC and OC emissions are passed in through xte()
do jt=ITRACBC,ITRACBC+1
  do jk=1,kl
    bce=bce+xte(:,jk,jt)*(aphp1(:,jk)-aphp1(:,jk+1))/grav
  enddo
enddo
do jt=ITRACOC,ITRACOC+1
  do jk=1,kl
    oce=oce+xte(:,jk,jt)*(aphp1(:,jk)-aphp1(:,jk+1))/grav
  enddo
enddo

! Total biomass burning primary emissions
bbem=emissfield(:,ibcb1)+emissfield(:,ibcb2)+1.3*(emissfield(:,iocb1)+emissfield(:,iocb2))

! ZVDRD   DRY DEPOSITION VELOCITY IN M/S
! ZVDRD(JL,1)  FOR SO2 GAS
! ZVDRD(JL,2)  FOR AEROSOLS
gdp=grav/(aphp1(:,1)-aphp1(:,2))
nstep=int(ztmst/120.01)+1
ddt=ztmst/real(nstep)

pxtm1new=xtg(1:ifull,1,itracso2)
do ii=1,nstep
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*zvdrd(:,1)*gdp)+ddt*xte(:,1,ITRACSO2))        &
      /(1.+0.5*ddt*p1mxtm1*zvdrd(:,1)*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
so2dd=(xtg(1:ifull,1,itracso2)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACSO2)/gdp
xte(:,1,ITRACSO2)  =xte(:,1,ITRACSO2)  -so2dd*gdp
  
pxtm1new=xtg(1:ifull,1,itracso2+1)
do ii=1,nstep  
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*zvdrd(:,2)*gdp)+ddt*xte(:,1,ITRACSO2+1))      &
        /(1.+0.5*ddt*p1mxtm1*zvdrd(:,2)*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
so4dd=(xtg(1:ifull,1,itracso2+1)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACSO2+1)/gdp
xte(:,1,ITRACSO2+1)=xte(:,1,ITRACSO2+1)-so4dd*gdp

pxtm1new=xtg(1:ifull,1,ITRACBC)
do ii=1,nstep  
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*ZVDPHOBIC*gdp)+ddt*xte(:,1,ITRACBC))          &
          /(1.+0.5*ddt*p1mxtm1*ZVDPHOBIC*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
ZHILBCO=(xtg(1:ifull,1,ITRACBC)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACBC)/gdp
xte(:,1,itracbc)  =xte(:,1,itracbc)  -zhilbco*gdp

pxtm1new=xtg(1:ifull,1,ITRACBC+1)
do ii=1,nstep
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*ZVDRD(:,2)*gdp)+ddt*xte(:,1,ITRACBC+1))       &
          /(1.+0.5*ddt*p1mxtm1*ZVDRD(:,2)*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
ZHILBCY=(xtg(1:ifull,1,ITRACBC+1)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACBC+1)/gdp
xte(:,1,itracbc+1)=xte(:,1,itracbc+1)-zhilbcy*gdp

pxtm1new=xtg(1:ifull,1,ITRACOC)
do ii=1,nstep
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*ZVDPHOBIC*gdp)+ddt*xte(:,1,ITRACOC))          &
          /(1.+0.5*ddt*p1mxtm1*ZVDPHOBIC*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
ZHILOCO=(xtg(1:ifull,1,ITRACOC)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACOC)/gdp
xte(:,1,itracoc)  =xte(:,1,itracoc)  -zhiloco*gdp

pxtm1new=xtg(1:ifull,1,ITRACOC+1)
do ii=1,nstep
  pxtm1new=(pxtm1new*(1.-0.5*ddt*p1mxtm1*ZVDRD(:,2)*gdp)+ddt*xte(:,1,ITRACOC+1))       &
          /(1.+0.5*ddt*p1mxtm1*ZVDRD(:,2)*gdp)
  pxtm1new=max(0.,pxtm1new)
end do
ZHILOCY=(xtg(1:ifull,1,ITRACOC+1)-pxtm1new)/(ztmst*gdp)+xte(:,1,ITRACOC+1)/gdp
xte(:,1,itracoc+1)=xte(:,1,itracoc+1)-zhilocy*gdp

bcdd=bcdd+zhilbco+zhilbcy
ocdd=ocdd+zhiloco+zhilocy

return
end subroutine xtemiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt sink

SUBROUTINE XTSINK (PTMST,PXTE)
!
!   *XTSINK*  CALCULATES THE DECREASE OF TRACER CONCENTRATION
!             FOR  A GIVEN HALF-LIFE-TIME
!
!   JOHANN FEICHTER               UNI-HAMBURG    08-91
!
!   PURPOSE
!  ---------------
!   THE MASS MIXING-RATIO OF TRACERS IS MULTIPLIED WITH
!   EXP(ALOG(0.5)*TIME-STEP/HALF-LIFE-TIME).
!   THIS ROUTINE COULD ALSO BE USED FOR EMISSION OR SINK
!   ABOVE THE SURFACE
!

implicit none

! Argument list
REAL, intent(in) :: PTMST
REAL, dimension(ifull,kl,naero), intent(out) :: PXTE

! Local data, functions etc
real, dimension(ifull) :: zxtp1,zdxtdt
real pqtmst,zfac,zdecay

integer jk

! Start code : ----------------------------------------------------------

pxte(:,:,:)=0. !Very important!

PQTMST=1./PTMST
ZFAC=ALOG(0.5)*PTMST

ZDECAY=EXP(ZFAC/86400.) ! 1 day
DO JK=1,kl
  ZXTP1=xtg(1:ifull,JK,ITRACBC)+PXTE(1:ifull,JK,ITRACBC)*PTMST
  ZXTP1=ZXTP1*ZDECAY
  ZDXTDT=(ZXTP1-xtg(1:ifull,JK,ITRACBC))*PQTMST-PXTE(1:ifull,JK,ITRACBC)
  PXTE(1:ifull,JK,ITRACBC)  =PXTE(1:ifull,JK,ITRACBC)  +ZDXTDT
  PXTE(1:ifull,JK,ITRACBC+1)=PXTE(1:ifull,JK,ITRACBC+1)-ZDXTDT 
end do

ZDECAY=EXP(ZFAC/86400.) ! 1 day
DO JK=1,kl
  ZXTP1=xtg(1:ifull,JK,ITRACOC)+PXTE(1:ifull,JK,ITRACOC)*PTMST
  ZXTP1=ZXTP1*ZDECAY
  ZDXTDT=(ZXTP1-xtg(1:ifull,JK,ITRACOC))*PQTMST-PXTE(1:ifull,JK,ITRACOC)
  PXTE(1:ifull,JK,ITRACOC)  =PXTE(1:ifull,JK,ITRACOC)  +ZDXTDT
  PXTE(1:ifull,JK,ITRACOC+1)=PXTE(1:ifull,JK,ITRACOC+1)-ZDXTDT 
end do

RETURN
END subroutine xtsink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt chemie

SUBROUTINE XTCHEMIE(KTOP, PTMST,zdayfac,PDPP1, PMRATEP, PFPREC,                      & !Inputs
                    PCLCOVER, PMLWC, PRHOP1, PTP1, sg, xtm1, pfevap,                 & !Inputs
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs, & !Inputs
                    prscav,pclcon,fracc,pccw,pfconv,xtu,                             & !Inputs
                    conwd,                                                           & !In and Out
                    xte,so2oh,so2h2,so2o3,dmsoh,dmsn3)                                 !Outputs

! Inputs
! ktop: top level for aerosol processes (set to 1, counting downwards from top)
! ptmst: timestep (seconds; tdt in main program)
! pdpp1: delta p (si units)
! pmratep: precip formation rate (kg/kg/s)
! pfprec: rainfall flux (entering from above) (kg/m2/s)
! pclcover: liquid-water cloud fraction (input; don't pass in cfrac though)
! pmlwc: liquid-water mixing ratio (kg/kg)
! prhop1: density of air (kg/m3)
! ptp1: temperature (k)
! sg: net solar radiation at ground (W/m2; used to determine if daytime)
! xtm1: tracer mixing ratio (kg/kg)
! pfevap: rainfall flux evaporating in layer k (kg/m2/s)
! pfsnow: snowfall flux (entering from above) (kg/m2/s)
! pfsubl: snowfall flux evaporating in layer k (kg/m2/s)
! pcfcover: ice cloud fraction
! pmiwc: ice mixing ratio (kg/kg)
! pmaccr: accretion rate (kg/kg/s)
! pfmelt: snowfall flux melting in layer k (kg/m2/s)
! pfstay: snowfall flux staying in layer k (kg/m2/s)
! pqfsed: fractional ice sedimentation in timestep
! plambs: slope (lambda) for snow crystal size distribution (m**-1)
! prscav: fractional rain scavenging rate in time step (needs to be mult. by coll. eff.)
! pclcon: convective cloud fraction
! fracc: Convective rain fraction
! pccw: convective cloud water mixing ratio (kg/kg)
! pfconv: convective rainfall flux (kg/m2/s)
! xtu: tracer mixing ratio in convective updraught (kg/kg)

! In & Out
! conwd: convective wet scavenging (diagnostic: kg/m2/s)

! Outputs
! xte: tracer tendency (kg/kg/s)
! so2oh: oxidation of SO2 by OH (diagnostic)
! so2h2: oxidation of SO2 by H2O2 (diagnostic)
! so2o3: oxidation of SO2 by O3 (diagnostic)
! dmsoh: oxidation of DMS by OH (diagnostic)
! dmsn3: oxidation of DMS by NO3 (diagnostic)

!**** *XTCHEMIE*  CALCULATES DRY AND WET CHEMISTRY
!
!      J. FEICHTER             UNI HAMBURG    30/06/92
!
!      PURPOSE
!      ---------
!      THIS ROUTINE COMPUTES THE OXIDATION AND THE WET SCAVENGING
!      OF CHEMICAL SPECIES.
!
!**    INTERFACE
!      -----------
!      *XTCHEMIE*   IS CALLED FROM  progcld in CSIRO GCM
!
!      EXTERNALS
!      ------------
!          *XTWETDEP*  CALCULATES THE WET DEPOSITION

implicit none

! Argument list
INTEGER KTOP
REAL PTMST
REAL PDPP1(ifull,kl)
REAL PMRATEP(ifull,kl)
REAL PFPREC(ifull,kl)
REAL PFEVAP(ifull,kl)
REAL PCLCOVER(ifull,kl)
REAL PMLWC(ifull,kl)
REAL PRHOP1(ifull,kl)
REAL PTP1(ifull,kl)
real sg(ifull)
REAL XTM1(ifull,kl,naero)
real xtu(ifull,kl,naero)
REAL XTE(ifull,kl,naero)
real pfsnow(ifull,kl)
real pfconv(ifull,kl)
real pfsubl(ifull,kl)
real pcfcover(ifull,kl)
real pmiwc(ifull,kl)
real pmaccr(ifull,kl)
real pfmelt(ifull,kl)
real pfstay(ifull,kl)
real pqfsed(ifull,kl)
real plambs(ifull,kl)
real prscav(ifull,kl)
real pclcon(ifull,kl)
real fracc(ifull)
real pccw(ifull,kl)
real conwd(ifull,naero)

real dmsoh(ifull) !Diagnostic output
real dmsn3(ifull) !Diagnostic output
real so2oh(ifull) !Diagnostic output
real so2h2(ifull) !Diagnostic output
real so2o3(ifull) !Diagnostic output

! Local work arrays and variables
real so2oh3d(ifull,kl),dmsoh3d(ifull,kl),dmsn33d(ifull,kl)
!
REAL ZXTP10(ifull,kl),          ZXTP1C(ifull,kl),   &
     ZHENRY(ifull,kl),          ZKII(ifull,kl),     &
     ZSO4(ifull,kl),            ZRKH2O2(ifull,kl),  &
     ZSO4i(ifull,kl),           ZSO4C(ifull,kl),    &
     ZHENRYC(ifull,kl),         ZXTP1CON(ifull,kl), &
     zsolub(ifull,kl)
REAL ZZOH(ifull,kl),            ZZH2O2(ifull,kl),   &
     ZZO3(ifull,kl),            ZDUMMY(ifull,kl),   &
     ZZNO2(ifull,kl),           ZDXTE(ifull,kl,naero)
REAL ZDEP3D(ifull,kl),                              &
     ZAMU0(ifull),zlwcic(ifull,kl),ziwcic(ifull,kl)
real zrevap(ifull,kl),zso2ev(ifull,kl)
real xto(ifull,kl,naero),zx(ifull)
integer ZRDAYL(ifull)
real, dimension(ifull), intent(in) :: zdayfac
integer, parameter :: nfastox=0 !1 for "fast" in-cloud oxidation; 0 for "slow"
real, parameter :: zmin=1.e-20
real x,pcons2,pqtmst
real zlwcl,zlwcv,zhp,zqtp1,zrk,zrke
real zh_so2,zpfac,zp_so2,zf_so2,zp_h2o2
real zf_h2o2,zxtp1,ze1,ze2,ze3,zfac1,zrkfac
real zza,za21,za22,zph_o3,zf_o3,zdt
real zh2o2m,zso2m,zso4m,zsumh2o2,zsumo3
real zh_h2o2,zq,zso2mh,zdso2h,zso2l,zso4l
real zzb,zzp,zzq,zzp2,zqhp,za2,zheneff
real zrko3,zso2mo,zdso2o,zdso2tot,zfac
real zfac2,zxtp1dms,zso2,ztk23b
real zhil,zexp,zm,zdms,t,ztk1,tk3,zqt,zqt3
real zrhoair,zkno2o3,zkn2o5aq,zrx1,zrx2
real zkno2no3,zeqn2o5,ztk3,ztk2,zkn2o5
real zno3,zxtp1so2
integer jt,jk,jl,js1,js2,js3,js4,jn,niter

!    REACTION RATE SO2-OH
real, parameter :: ZK2I=2.0E-12
real, parameter :: ZK2=4.0E-31
real, parameter :: ZK2F=0.45

!   REACTION RATE DMS-NO3
real, parameter :: ZK3=1.9E-13

!   MOLECULAR WEIGHTS IN G
real, parameter :: ZMOLGS=32.064
real, parameter :: ZMOLGH2O2=34.01474
real, parameter :: ZMOLGAIR=28.84
real, parameter :: ZMOLGW=18.015

real, parameter :: ZHPBASE=2.5E-06
real, parameter :: ZE1K=1.1E-02
real, parameter :: ZE1H=2300.
real, parameter :: ZE2K=1.23
real, parameter :: ZE2H=3020.
real, parameter :: ZE3K=1.2E-02
real, parameter :: ZE3H=2010.
real, parameter :: ZQ298=1./298.
real, parameter :: ZRGAS=8.2E-02

real, parameter :: ZAVO=6.022E+23
real, parameter :: ZNAMAIR=1.E-03*ZAVO/ZMOLGAIR

real, parameter :: ZLWCMIN=1.E-07


! Start code : ----------------------------------------------------------
dmsoh(:)=0.
dmsn3(:)=0.
so2oh(:)=0.
so2h2(:)=0.
so2o3(:)=0.
so2oh3d(:,:)=0.
dmsoh3d(:,:)=0.
dmsn33d(:,:)=0.
xte(:,:,:)=0.
pcons2=1./(ptmst*grav)
where ( sg(:)>0. )
  zrdayl(:)=1
elsewhere
  zrdayl(:)=0  
end where

! Calculate xto, tracer mixing ratio outside convective updraughts
! Assumes pclcon < 1, but this shouldn't be a problem.
do jt=1,naero
  xto(:,:,jt)=(xtm1(:,:,jt)-pclcon(:,:)*xtu(:,:,jt))/(1.-pclcon(:,:))
enddo
xto=max(0.,xto)

!   CALCULATE THE ZRDAYL (=0 --> NIGHT; =1 --> DAY) AND
!                 ZAMUO  =  ZENITH ANGLE

!    CONSTANTS
PQTMST=1./PTMST
if(nfastox==0)then
   NITER=5  !Slow in-cloud oxidation
else 
   NITER=1  !Fast
endif

! Calculate in-cloud ql
where (pclcover(:,:)>zmin)
  zlwcic(:,:)=pmlwc(:,:)/pclcover(:,:)
elsewhere
  zlwcic(:,:)=0.
end where
where (pcfcover(:,:)>zmin)
  ziwcic(:,:)=pmiwc(:,:)/pcfcover(:,:)
elsewhere
  ziwcic(:,:)=0.
end where

!  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
DO JK=1,kl
  JS1=JK
  JS2=kl+JK
  JS3=2*kl+JK
  JS4=3*kl+JK

  ZX(:)=PRHOP1(:,JK)*1.E-03
  ZZOH(:,JK)=ZOXIDANT(:,JS1)
  ZZH2O2(:,JK)=ZOXIDANT(:,JS2)*ZX(:)
  ZZO3(:,JK)=ZOXIDANT(:,JS3)*ZX(:)
  ZZNO2(:,JK)=ZOXIDANT(:,JS4)*ZX(:)
end do

zhenry=0.
zhenryc=0.
zdxte=0.

!   PROCESSES WHICH ARE DIFERENT INSIDE AND OUTSIDE OF CLOUDS
ZSO4(:,ktop:kl)=XTO(:,ktop:kl,ITRACSO2+1)
ZSO4(:,ktop:kl)=AMAX1(0.,ZSO4(:,ktop:kl))
!
!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
DO JK=KTOP,kl
  DO JL=1,ifull
    IF(ZLWCIC(JL,JK)>ZMIN) THEN
      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4(JL,JK)*1000./(ZLWCIC(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
    ELSE
      ZRKH2O2(JL,JK)=0.
    ENDIF
  end do
end do

!   HETEROGENEOUS CHEMISTRY
DO JK=KTOP,kl
  DO JL=1,ifull
    ZXTP1=XTO(JL,JK,ITRACSO2)
    ZXTP10(JL,JK)=XTO(JL,JK,ITRACSO2)
    ZXTP1C(JL,JK)=XTO(JL,JK,ITRACSO2)
    IF(ZXTP1>ZMIN.AND.ZLWCIC(JL,JK)>ZMIN) THEN
      X=PRHOP1(JL,JK)

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
      ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)

      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
!    ZLWCL = LWC IN L/CM**3
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
!   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/FLOAT(NITER)

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
      ZSO4M=ZSO4(JL,JK)*XTOC(X,ZMOLGS)

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,NITER
        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
        ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT) & ! = zero if nfastox==1
               +nfastox * max (0., zso2m - zh2o2m )

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2

        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      end do  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1C(JL,JK)=ZXTP1-ZDSO2TOT
      ZSO4(JL,JK)=ZSO4(JL,JK)+ZDSO2TOT

      ZHENRY(JL,JK)=ZF_SO2
! Diagnostic only...
      ZFAC=PQTMST*CTOX(X,ZMOLGS)*PCLCOVER(JL,JK)
      ZFAC1=ZFAC*PDPP1(JL,JK)/grav
      ZFAC2=ZFAC*PRHOP1(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    ENDIF
  end do
end do

! Repeat the aqueous oxidation calculation for ice clouds.
ZSO4i(:,ktop:kl)=XTO(:,ktop:kl,ITRACSO2+1)
ZSO4i(:,ktop:kl)=AMAX1(0.,ZSO4i(:,ktop:kl))

! Repeat the aqueous oxidation calculation for convective clouds.
ZXTP1CON(:,ktop:kl)=XTU(:,ktop:kl,ITRACSO2)
ZSO4C(:,ktop:kl)=XTU(:,ktop:kl,ITRACSO2+1)
ZSO4C(:,ktop:kl)=AMAX1(0.,ZSO4C(:,ktop:kl))

! Comment from here when not using convective oxidation...

!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
DO JK=KTOP,kl
  DO JL=1,ifull
    IF(PCCW(JL,JK)>ZMIN) THEN
      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4C(JL,JK)*1000./(PCCW(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
    ELSE
      ZRKH2O2(JL,JK)=0.
    ENDIF
  ENDDO
ENDDO

!   HETEROGENEOUS CHEMISTRY
DO JK=KTOP,kl
  DO JL=1,ifull
    ZXTP1=XTU(JL,JK,ITRACSO2)
    IF(ZXTP1>ZMIN.AND.PCCW(JL,JK)>ZMIN) THEN
      X=PRHOP1(JL,JK)

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
      ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)

      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
!    ZLWCL = LWC IN L/CM**3
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
!   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/FLOAT(NITER)

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
      ZSO4M=ZSO4C(JL,JK)*XTOC(X,ZMOLGS)

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,NITER
        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
        ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT)+nfastox*max(0.,zso2m-zh2o2m) ! = zero if nfastox==1

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2
!
        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      ENDDO  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1CON(JL,JK)=ZXTP1CON(JL,JK)-ZDSO2TOT
      ZSO4C(JL,JK)=ZSO4C(JL,JK)+ZDSO2TOT
      ZHENRYC(JL,JK)=ZF_SO2
      ! Diagnostic only...
      ZFAC=PQTMST*CTOX(X,ZMOLGS)*pclcon(jl,jk)
      ZFAC1=ZFAC*PDPP1(JL,JK)/grav
      ZFAC2=ZFAC*PRHOP1(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    ENDIF
  ENDDO
ENDDO

!*******************************************************************************
!
!    CALCULATE THE WET DEPOSITION
!
DO JT=ITRACSO2,naero
  zdep3d=0.

  IF (LWETDEP(JT)) THEN          !True for all except DMS

    if(jt==itracso2) then        !SO2
      zsolub(:,:)=zhenry(:,:)

    elseif(jt==itracso2+1) then  !sulfate
      zxtp1c(:,:)=zso4(:,:)
      zxtp10(:,:)=zso4i(:,:)
      zxtp1con(:,:)=zso4c(:,:)
      zsolub (:,:)=0.6

    else !Carbonaceous aerosol and mineral dust
      zxtp10(:,:)=xto(:,:,jt)
      zxtp1c(:,:)=xto(:,:,jt)
      zxtp1con(:,:)=xtu(:,:,jt)
        
      if(jt==itracbc.or.jt==itracoc)then  !hydrophobic BC and OC
        zsolub(:,:)=0.
      elseif(jt==itracbc+1.or.jt==itracoc+1)then !hydrophilic BC and OC
        zsolub(:,:)=0.2
      elseif(jt>=itracdu.and.jt<itracdu+ndust)then !hydrophobic dust (first 4 dust vars)
        zsolub(:,:)=0.05
!      elseif(jt>=itracdu+ndust)then !hydrophilic dust !hydrophilic dust (last 4 dust vars)
!        zsolub(:,:)=1.
      endif

    endif

    CALL XTWETDEP( JT,                                         &
                   PTMST, PCONS2,                              &
                   PDPP1,                                      &
                   PMRATEP, PFPREC, PFEVAP,                    &
                   PCLCOVER, PRHOP1, zsolub, pmlwc, ptp1,      &
                   pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt, &
                   pfstay,pqfsed,plambs,prscav,pfconv,pclcon,  & 
                   fracc,                                      & !Inputs
                   ZXTP10, ZXTP1C,ZDEP3D,conwd,zxtp1con,       & !In and Out
                   zrevap )

!   CALCULATE NEW CHEMISTRY AND SCAVENGING TENDENCIES
    DO JK=KTOP,kl
      DO JL=1,ifull
        ZXTP1=(1.-pclcover(jl,jk)-pclcon(jl,jk))*ZXTP10(JL,JK)+ &
               PCLCOVER(JL,JK)*ZXTP1C(JL,JK)+                   &
               pclcon(jl,jk)*zxtp1con(jl,jk)
        zxtp1=max(zxtp1,0.)
        ZDXTE(JL,JK,JT)=(ZXTP1-XTM1(JL,JK,JT))*PQTMST  !Total tendency (Dep + chem)
      end do
    end do

! Note that wd as coded here includes the below-cloud convective scavenging/evaporation
    if(jt==itracso2)then
      do jk=1,kl
        so2wd(:)=so2wd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
    elseif(jt==itracso2+1)then
      do jk=1,kl
        so4wd(:)=so4wd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
    elseif(jt==itracbc.or.jt==itracbc+1) then
      do jk=1,kl
        bcwd(:)=bcwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      end do
    elseif(jt==itracoc.or.jt==itracoc+1) then
      do jk=1,kl
        ocwd(:)=ocwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      end do
    elseif(jt>=itracdu.and.jt<itracdu+ndust)then
      do jk=1,kl
        dustwd(:)=dustwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
    endif

  ENDIF  !lwetdep
end do

!    CHANGE THE TOTAL TENDENCIES

!  ZDXTE(ITRACSO2) = TENDENCY OF SO2
!  ZDXTE(ITRACSO2+1) = CHEMISTRY AND SCAVENGING TENDENCY OF TOTAL
!       SULFATE
xte(:,ktop:kl,ITRACSO2)=xte(:,ktop:kl,ITRACSO2)+zdxte(:,ktop:kl,ITRACSO2)
XTE(:,ktop:kl,ITRACSO2+1)=XTE(:,ktop:kl,ITRACSO2+1)+ZDXTE(:,ktop:kl,ITRACSO2+1)
      
! Update wet-scavenging tendencies for non-sulfate aerosols (JT >= ITRACBC)
! This covers BC, OC and dust in the current version of the model.
do jt = itracbc, naero
  if(lwetdep(jt)) then
    xte(:,:,jt) = xte(:,:,jt) + zdxte(:,:,jt)
  endif
enddo

!   CALCULATE THE DAY-LENGTH
! Need to hack this because of irritating CSIRO coding! (NH+SH latitudes concatenated!)
!      ZDAYL=0.
!      DO 402 JL=1,lon
!        IF(ZRDAYL(JL)==1) THEN
!          ZDAYL=ZDAYL+1.
!        ENDIF
!  402 CONTINUE
!      ZDAYFAC(:)=0.
!      ZNLON=FLOAT(lon)
!      IF(ZDAYL/=0.) ZDAYFAC(1)=ZNLON/ZDAYL !NH
!      IF(ZDAYL/=znlon) ZDAYFAC(2)=ZNLON/(znlon-ZDAYL) !SH

!   DAY-TIME CHEMISTRY
DO JK=1,kl
  DO JL=1,ifull
    IF(ZRDAYL(JL)==1) THEN
      !ins=(jl-1)/lon + 1 !hemisphere index
      X=PRHOP1(JL,JK)
      ZXTP1SO2=XTM1(JL,JK,ITRACSO2)+XTE(JL,JK,ITRACSO2)*PTMST
      ZTK2=ZK2*(PTP1(JL,JK)/300.)**(-3.3)
      ZM=X*ZNAMAIR
      ZHIL=ZTK2*ZM/ZK2I
      ZEXP=ALOG10(ZHIL)
      ZEXP=1./(1.+ZEXP*ZEXP)
      ZTK23B=ZTK2*ZM/(1.+ZHIL)*ZK2F**ZEXP
      !ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(ins)
      ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(jl)
      ZSO2=AMIN1(ZSO2,ZXTP1SO2*PQTMST)
      ZSO2=AMAX1(ZSO2,0.)
      XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)-ZSO2
      XTE(JL,JK,ITRACSO2+1)=XTE(JL,JK,ITRACSO2+1)+ZSO2
      so2oh3d(jl,jk)=zso2

      ZXTP1DMS=XTM1(JL,JK,ITRACSO2-1)+XTE(JL,JK,ITRACSO2-1)*PTMST
      IF(ZXTP1DMS<=ZMIN) THEN
        ZDMS=0.
      ELSE
        T=PTP1(JL,JK)
        ztk1=1.646e-10-1.850e-12*t+8.151e-15*t**2-1.253e-17*t**3 !Cubic fit good enough
        ztk1=max(ztk1,5.e-12) !Because cubic falls away for T > 300 K
        ztk1=1.5*ztk1         !This is the fudge factor to account for other oxidants
        !ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(ins)
        ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(jl)
        ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
        XTE(JL,JK,ITRACSO2-1)=XTE(JL,JK,ITRACSO2-1)-ZDMS
        XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)+ZDMS
        dmsoh3d(jl,jk)=zdms
      ENDIF
    ENDIF
  end do
end do

!   NIGHT-TIME CHEMISTRY
DO JK=1,kl
  DO JL=1,ifull
    IF(ZRDAYL(JL)/=1) THEN
      X=PRHOP1(JL,JK)
      ZXTP1DMS=XTM1(JL,JK,ITRACSO2-1)+XTE(JL,JK,ITRACSO2-1)*PTMST
      IF(ZXTP1DMS<=ZMIN) THEN
        ZDMS=0.
      ELSE
        ZTK3=ZK3*EXP(520./PTP1(JL,JK))
!    CALCULATE THE STEADY STATE CONCENTRATION OF NO3
        ZQT=1./PTP1(JL,JK)
        ZQT3=300.*ZQT
        ZRHOAIR=PRHOP1(JL,JK)*ZNAMAIR
        ZKNO2O3=1.2E-13*EXP(-2450.*ZQT)
        ZKN2O5AQ=0.1E-04
        ZRX1=2.2E-30*ZQT3**3.9*ZRHOAIR
        ZRX2=1.5E-12*ZQT3**0.7
        ZKNO2NO3=ZRX1/(1.+ZRX1/ZRX2)*0.6**(1./(1.+(ALOG10(ZRX1/ZRX2))**2))
        ZEQN2O5=4.E-27*EXP(10930.*ZQT)
        ZKN2O5=ZKNO2NO3/ZEQN2O5

        ZNO3=ZKNO2O3*(ZKN2O5+ZKN2O5AQ)*ZZNO2(JL,JK)*ZZO3(JL,JK)
        ZZQ=ZKNO2NO3*ZKN2O5AQ*ZZNO2(JL,JK)+(ZKN2O5+ZKN2O5AQ)*ZTK3*ZXTP1DMS*XTOC(X,ZMOLGS)
        IF(ZZQ>0.) THEN
          ZNO3=ZNO3/ZZQ
        ELSE
          ZNO3=0.
        ENDIF
        ZDMS=ZXTP1DMS*ZNO3*ZTK3
        ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
        XTE(JL,JK,ITRACSO2-1)=XTE(JL,JK,ITRACSO2-1)-ZDMS
        XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)+ZDMS
        dmsn33d(jl,jk)=zdms
      ENDIF
    ENDIF
  end do
end do

! Calculate tendency of SO2 due to oxidation by OH (diagnostic) and ox. tendencies of DMS
do jk=1,kl
  so2oh(:)=so2oh(:)+so2oh3d(:,jk)*pdpp1(:,jk)/grav
  dmsoh(:)=dmsoh(:)+dmsoh3d(:,jk)*pdpp1(:,jk)/grav
  dmsn3(:)=dmsn3(:)+dmsn33d(:,jk)*pdpp1(:,jk)/grav
enddo

RETURN
END subroutine xtchemie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt wetdep

SUBROUTINE XTWETDEP(KTRAC,                                                           &
                    PTMST, PCONS2,                                                   &
                    PDPP1,                                                           &
                    PMRATEP, PFPREC, PFEVAP,                                         &
                    PCLCOVER, PRHOP1, PSOLUB, pmlwc, ptp1,                           &
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs, &
                    prscav,pfconv,pclcon,fracc,                                      & !Inputs
                    PXTP10, PXTP1C, PDEP3D, conwd,pxtp1con,                          & !In & Out
                    prevap)                                                            !Outputs

!
!   *XTWETDEP* CALCULATES THE WET DEPOSITION OF TRACE GASES OR AEROSOLS
!
!   JOHANN FEICHTER              UNI HAMBURG            08-91
!
!   PURPOSE
!  ---------
!   TO CALCULATE THE WET SCAVENGING OF GASES OR AEROSOLS IN CLOUDS
!
!   INTERFACE
!  -------------
!   THIS ROUTINE IS CALLED FROM *XTCHEM*
!
!  METHOD
!  -------
!
!   NO EXTERNALS
!---------------
!

implicit none

! Argument list
INTEGER KTRAC
REAL PTMST
REAL PCONS2
REAL PXTP10(ifull,kl)   !Tracer m.r. outside liquid-water cloud (clear air/ice cloud)
REAL PXTP1C(ifull,kl)   !Tracer m.r.  inside liquid-water cloud
real pxtp1con(ifull,kl) !Tracer m.r.  inside convective cloud
REAL PDPP1(ifull,kl)
REAL PMRATEP(ifull,kl)
REAL PFPREC(ifull,kl)
REAL PFEVAP(ifull,kl)
REAL PDEP3D(ifull,kl)
REAL PCLCOVER(ifull,kl)
REAL PRHOP1(ifull,kl)
REAL PSOLUB(ifull,kl)
real pmlwc(ifull,kl)
real ptp1(ifull,kl)  !temperature
real pfsnow(ifull,kl)
real pfconv(ifull,kl)
real pclcon(ifull,kl)
real fracc(ifull)       !Convective rain fraction (originially set to 0.1)
real pfsubl(ifull,kl)
real pcfcover(ifull,kl)
real pmiwc(ifull,kl)
real pmaccr(ifull,kl)
real pfmelt(ifull,kl)
real pfstay(ifull,kl)
real pqfsed(ifull,kl)
real plambs(ifull,kl)
real prscav(ifull,kl)
real prevap(ifull,kl)
real conwd(ifull,naero)

! Local work arrays and variables
REAL ZDEP(ifull) !Only needed for old code
REAL ZDEPS(ifull),ZDEPR(ifull),    &
     ZMTOF(ifull),   ZFTOM(ifull), &
     ZCLEAR(ifull), ZCLR0(ifull)

integer, parameter :: ktop = 2    !Top level for wet deposition (counting from top)
! Allow in-cloud scavenging in ice clouds for hydrophobic BC and OC, and dust
real, dimension(naero), parameter :: Ecols = (/ 0.00,0.00,0.00,0.05,0.00,0.05,0.00,0.05,0.05,0.05,0.05/)
!Below-cloud collection eff. for rain
real, dimension(naero), parameter :: zcollefr = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.10,0.20,0.40/)
!Below-cloud collection eff. for snow
real, dimension(naero), parameter :: zcollefs = (/0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.04,0.08/)
!Retention coeff. on riming
real, dimension(naero), parameter :: Rcoeff = (/1.00,0.62,1.00,0.00,1.00,0.00,1.00,1.00,1.00,1.00,1.00/)
!Relative reevaporation rate
real, dimension(naero), parameter :: Evfac = (/0.25,1.00,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25/)

real, dimension(ifull,kl) :: zcollefc  !Collection efficiency for convection

integer, dimension(ifull) :: kbase
real pqtmst,ziicscav,xdep,zicscav,xicscav,zevap
real zilcscav,plambda,zbcscav,xbcscav,zstay,xstay
real xevap,frc,pcevap

integer jk,jl

real, parameter :: zmin = 1.e-20

! Start code : ----------------------------------------------------------

PQTMST=1./PTMST

zdepr(:)=0.
zdeps(:)=0.
prevap(:,:)=0.

! Search for convective cloud base
kbase(:)=kl+1
do jk=ktop,kl
  where (pclcon(:,jk)>zmin)
    kbase(:)=jk
  end where
enddo

!     BEGIN OF VERTICAL LOOP
do JK=KTOP,kl
  ZCLEAR(:)=1.-PCLCOVER(:,JK)-pcfcover(:,jk)-pclcon(:,jk)
  ZCLR0(:)=1.-PCLCOVER(:,JK)-pclcon(:,jk) !Clear air or ice cloud (applies to pxtp10)
  ZMTOF(:)=PDPP1(:,JK)*PCONS2
  ZFTOM(:)=1./ZMTOF(:)
  PXTP1C(:,JK)=AMAX1(0.,PXTP1C(:,JK))
  PXTP10(:,JK)=AMAX1(0.,PXTP10(:,JK))

! In-cloud ice scavenging (including vertical redistribution when snow
! evaporates or falls into a layer). Include accretion of ql by snow.
  if ( Ecols(ktrac)>zmin ) then
    do jl=1,ifull
      if ( pmiwc(jl,jk)>zmin ) then
        ziicscav=Ecols(ktrac)*pqfsed(jl,jk) !qfsed is the fractional sedimentation in dt
        ziicscav=min(ziicscav, 1.)        
        xdep=pxtp10(jl,jk)*ziicscav*pcfcover(jl,jk)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
        !pxtp10(jl,jk)=pxtp10(jl,jk)*(zclear(jl)+(1.-ziicscav)*pcfcover(jl,jk))/(1.-pclcover(jl,jk))
        pxtp10(jl,jk)=pxtp10(jl,jk)-xdep/zclr0(jl) ! MJT suggestion
        zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
      endif
    enddo
  endif

! This loop does riming (accretion of liquid water by falling snow)
  if ( Rcoeff(ktrac)>zmin ) then
    do jl=1,ifull
      if ( pmlwc(jl,jk)>zmin ) then
        zilcscav=Rcoeff(ktrac)*psolub(jl,jk)*(pmaccr(jl,jk)*ptmst/pmlwc(jl,jk))
        zilcscav=min(zilcscav, 1.)        
        xdep=pxtp1c(jl,jk)*zilcscav*pclcover(jl,jk)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
        pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zilcscav)
        zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
      endif
    enddo
  endif

! Below-cloud scavenging by snow
  if(zcollefs(ktrac)>zmin)then
    do jl=1,ifull
      if(pfsnow(jl,jk)>zmin.and.pclcover(jl,jk)<1.-zmin)then
        plambda=min(plambs(jl,jk),8.e3) !Cut it off at about -30 deg. C
        zbcscav=zcollefs(ktrac)*plambda*pfsnow(jl,jk)*ptmst/(2.*rhos)
        zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
        xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
        pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
        zdeps(jl)=zdeps(jl)+xbcscav*zmtof(jl)
      end if
    enddo
  endif

! Redistribution by snow that evaporates or stays in layer
  do jl=1,ifull
    if (pfsnow(jl,jk)>zmin) then
      zstay=(pfsubl(jl,jk)+pfstay(jl,jk))/pfsnow(jl,jk)
      zstay=min(1., zstay)
      xstay=zdeps(jl)*zstay*zftom(jl)
      zdeps(jl)=zdeps(jl)*(1.-zstay)
      zdeps(jl)=max(0.,zdeps(jl))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xstay
      if(zclr0(jl)>zmin)then
        pxtp10(jl,jk)=pxtp10(jl,jk)+xstay/zclr0(jl)
      else
        pxtp1c(jl,jk)=pxtp1c(jl,jk)+xstay/pclcover(jl,jk)
      endif
    end if
  enddo

  ! Melting of snow... 
  where (pfmelt(:,jk)>zmin)
    zdepr(:)=zdepr(:)+zdeps(:)
    zdeps(:)=0.
  end where

  !  In-cloud scavenging by warm-rain processes (autoconversion and collection)
  do jl=1,ifull
    !if(pmratep(jl,jk)>zmin) then
    if(pmratep(jl,jk)>zmin.and.pmlwc(jl,jk)>zmin) then ! MJT suggestion
      zicscav=psolub(jl,jk)*(pmratep(jl,jk)*ptmst/pmlwc(jl,jk))
      zicscav=min(zicscav,1.)
      xicscav=pxtp1c(jl,jk)*zicscav*pclcover(jl,jk) !gridbox mean
      pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zicscav)
      pdep3d(jl,jk)=pdep3d(jl,jk)+xicscav
      zdepr(jl)=zdepr(jl)+xicscav*zmtof(jl)
    end if
  enddo

  ! Below-cloud scavenging by stratiform rain (conv done below)
  if(zcollefr(ktrac)>zmin)then
    do jl=1,ifull
      if(pfprec(jl,jk)>zmin.and.zclr0(jl)>zmin)then
        zbcscav=zcollefr(ktrac)*prscav(jl,jk)
        zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
        xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
        pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
        zdepr(jl)=zdepr(jl)+xbcscav*zmtof(jl)
      end if
    enddo
  endif

  ! Reevaporation of rain
  do jl=1,ifull
    if (pfprec(jl,jk)>zmin.and.zclear(jl)>zmin) then
      zevap=pfevap(jl,jk)/pfprec(jl,jk)
      zevap=max(0.,min(1., zevap))
      if(zevap<1.)zevap=Evfac(ktrac)*zevap
      xevap=zdepr(jl)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
      zdepr(jl)=max(0.,zdepr(jl)*(1.-zevap))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
      prevap(jl,jk)=xevap
      pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
    end if
  enddo

end do !   END OF VERTICAL LOOP

! Use collection efficiencies for rain below melting level, snow above

where (ptp1(1:ifull,ktop:kl)>273.15)
  zcollefc(1:ifull,ktop:kl) = zcollefr(ktrac)
elsewhere
  zcollefc(1:ifull,ktop:kl) = zcollefs(ktrac)  
end where

! Now do the convective below-cloud bit...
! In-cloud convective bit was done in convjlm.
do jk=ktop,kl
  do jl=1,ifull
    zmtof(jl)=pdpp1(jl,jk)*pcons2
    zftom(jl)=1./zmtof(jl)
    zclr0(jl)=1.-pclcover(jl,jk)-pclcon(jl,jk)
          
! Below-cloud scavenging by convective precipitation (assumed to be rain)
    if(pfconv(jl,jk-1)>zmin.and.zclr0(jl)>zmin)then
      Frc=max(0.,pfconv(jl,jk-1)/fracc(jl))
      zbcscav=zcollefc(jl,jk)*fracc(jl)*0.24*ptmst*sqrt(Frc*sqrt(Frc))
      zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
      xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
      conwd(jl,ktrac)=conwd(jl,ktrac)+xbcscav*zmtof(jl)
      pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
      pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
    endif

! Below-cloud reevaporation of convective rain
    if(jk>kbase(jl).and.pfconv(jl,jk-1)>zmin.and.zclr0(jl)>zmin)then
      pcevap=pfconv(jl,jk-1)-pfconv(jl,jk)
      zevap=pcevap/pfconv(jl,jk-1)
      zevap=max(0.,min(1.,zevap))
      if(zevap<1.)zevap=Evfac(ktrac)*zevap
      xevap=conwd(jl,ktrac)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
      conwd(jl,ktrac)=max(0.,conwd(jl,ktrac)*(1.-zevap))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
      prevap(jl,jk)=prevap(jl,jk)+xevap
      pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
    end if
  enddo
enddo

RETURN
END subroutine xtwetdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust settling

subroutine dsettling(tdt,rhoa,tmp,delz,prf)

implicit none

!     Inputs:
real, intent(in) :: tdt                         !timestep (s)
real, dimension(1:ifull,kl), intent(in) :: rhoa !air density (kg/m3)
real, dimension(:,:), intent(in) :: tmp         !temperature (K)
real, dimension(ifull,kl), intent(in) :: delz   !Lowest layer thickness (m)
real, dimension(ifull,kl), intent(in) :: prf    !Pressure (hPa)

! Local work arrays and variables
real, dimension(ifull) :: dcol1, dcol2
real, dimension(kl) :: vd_cor
real dzmin,dtmax,vsettl,dt_settl
real pres,c_stokes,corr
real c_cun
integer n,k,l,iq,ndt_settl

real, parameter :: dyn_visc = 1.5E-5

! Start code : ----------------------------------------------------------

! Calculate integrated column dust before settling
dcol1 = 0.
do n=itracdu,itracdu+ndust-1
  do k=1,kl
    dcol1 = dcol1 + rhoa(:,k) * xtg(1:ifull,k,n)* delz (:,k)
  enddo
enddo

! MJT notes - need to calculate dzmin for each iq gridpoint as
! otherwise the answer changes with different numbers of processors

do iq=1,ifull

  dzmin = delz(iq,1) ! MJT suggestion, as delz(:,1) < delz(:,2:kl)    

  do k = 1, NDUST
    ! Settling velocity (m/s) for each soil classes (Stokes Law)
    ! DUSTDEN     soil class density             (kg/m3)
    ! DUSTREFF    effective radius of soil class (m)
    ! dyn_visc    dynamic viscosity              (kg/m2/s)
    ! grav        gravity                        (m/s2)
    ! 0.5         upper limit with temp correction
    vsettl = 2./9. * grav * DUSTDEN(k) * DUSTREFF(k)**2 / (0.5*dyn_visc)

    ! Determine the maximum time-step satisying the CFL condition:
    ! dt <= (dz)_min / v_settl (use 0.3 as MJT suggestion)
    dtmax = 0.3 * dzmin / vsettl
    ndt_settl = max(1, int( tdt /dtmax ) )


    ! Solve the bidiagonal matrix (l,l)
    dt_settl = tdt / real(ndt_settl)
    do n = 1, ndt_settl

      ! Solve at the model top
      ! Dynamic viscosity
      C_Stokes = 1.458E-6 * TMP(iq,kl)**1.5/(TMP(iq,kl)+110.4) 
      ! Cuningham correction
      Corr = 6.6E-8*prf(iq,kl)/1013.*TMP(iq,kl)/293.15
      C_Cun = 1. + 1.249*corr/dustreff(k)
      ! Settling velocity
      Vd_cor(kl) =2./9.*grav*dustden(k)*dustreff(k)**2/C_Stokes*C_Cun
      ! Update mixing ratio
      xtg(iq,kl,k+itracdu-1) = xtg(iq,kl,k+itracdu-1) / (1. + dt_settl*VD_cor(kl)/DELZ(iq,kl))

      ! Solve each vertical layer successively (layer l)
      do l = kl-1,1,-1
        ! Dynamic viscosity
        C_Stokes = 1.458E-6*TMP(iq,l)**1.5/(TMP(iq,l)+110.4) 
        ! Cuningham correction
        Corr = 6.6E-8*prf(iq,l)/1013.*TMP(iq,l)/293.15
        C_Cun = 1. + 1.249*corr/dustreff(k)
        ! Settling velocity
        Vd_cor(l) = 2./9.*grav*dustden(k)*dustreff(k)*dustreff(k)/C_Stokes*C_Cun
        ! Update mixing ratio
        xtg(iq,l,k+itracdu-1) = 1./(1. + dt_settl*Vd_cor(l)/DELZ(iq,l))                          &
            *(xtg(iq,l,k+itracdu-1) + dt_settl*Vd_cor(l+1)/DELZ(iq,l)*xtg(iq,l+1,k+itracdu-1)    &
             *rhoa(iq,l+1)/rhoa(iq,l))  ! MJT suggestion
      end do
    end do
  end do
end do
  
! Calculate integrated column dust after settling
dcol2 = 0.
do n=itracdu,itracdu+ndust-1
  do k=1,kl
    dcol2 = dcol2 + rhoa(:,k) * xtg(1:ifull,k,n) * delz (:,k)
  enddo
enddo

! Calculate deposition flux to surface
dustdd = dustdd + (dcol1-dcol2)/tdt

return
end subroutine dsettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust emissions

subroutine dustem(tdt,rhoa,wg,w10m,dz1,vt,snowd,land)

implicit none

!     Inputs:
real, intent(in) :: tdt                         !Leapfrog timestep (s)
real, dimension(ifull), intent(in) :: rhoa      !air density (kg/m3)
real, dimension(ifull), intent(in) :: wg        !ground wetness (fraction of field capacity)
real, dimension(ifull), intent(in) :: w10m      !10m windspeed (m/s)
real, dimension(ifull), intent(in) :: dz1       !Lowest layer thickness (m)
real, dimension(ifull), intent(in) :: vt        !Transfer velocity at surface for dry deposition (m/s)
real, dimension(ifull), intent(in) :: snowd     !Snow depth (mm equivalent water)
logical, dimension(ifull), intent(in) :: land

! Local work arrays and variables
integer, dimension(ndust), parameter :: ipoint = (/ 3, 2, 2, 2 /) !Pointer used for dust classes (sand, silt, clay)
! This array gives fraction of source in each size bin.
! Doesn't quite add to 1, because larger sizes (omitted) account for some too.
! All source is in first four bins, even when using eight bins, since next four are hydrophilic.
real, dimension(ndust), parameter :: frac_s = (/ 0.1, 0.25, 0.25, 0.25 /)
real, dimension(ifull) :: snowa     !Estimated snow areal coverage
real, dimension(ifull) :: u_ts0,u_ts
real, dimension(ifull) :: srce,dsrc,airmas
real, dimension(ifull) :: a,b,xold,xtendd,veff
real, dimension(ifull) :: airden

real g,den,diam,ddt
integer n,m,ii,nstep

! Start code : ----------------------------------------------------------

g = grav*1.e2
airden = rhoa*1.e-3

! Convert snow depth to estimated areal coverage (see Zender et al 2003, JGR 108, D14, 4416)
! Must convert from mm to m, then adjust by rho_l/rho_s=10.
! 0.05 m is the geometrical snow thickness for 100% areal coverage.
!hsnow = snowd*0.01 !Geometrical snow thickness in metres
snowa = min( 1., snowd/5. )

do n = 1, ndust
  ! Threshold velocity as a function of the dust density and the diameter
  ! from Bagnold (1941)
  den = dustden(n)*1.e-3
  diam = 2.*dustreff(n)*1.e2
  ! Pointer to the 3 classes considered in the source data files
  m = ipoint(n)
  
  ! Following is from Ginoux et al (2004) Env. Modelling & Software.
  u_ts0 = 0.13*1.e-2*sqrt(den*g*diam/airden)*sqrt(1.+0.006/den/g/diam**2.5)/ &
          sqrt(1.928*(1331.*diam**1.56+0.38)**0.092-1.)
  
  ! Case of surface dry enough to erode
  where ( wg<0.1 )
    !Tuning suggested for Asian source by P. Ginoux
    u_ts = max( 0., u_ts0*(1.2+0.2*alog10(max( 1.e-3, wg ))) )
  elsewhere
    ! Case of wet surface, no erosion
    u_ts = 100.
  end where
  
  ! MJT notes - erod should be zero for ocean points
    
  !srce = frac_s(n)*erod(i,m)*dxy(i) ! (m2)
  srce = frac_s(n)*erod(:,m) ! (fraction) - MJT suggestion
  dsrc = (1.-snowa)*Ch_dust*srce*W10m*W10m*(W10m-u_ts) ! (kg/s/m2)
  dsrc = max( 0., dsrc )

! Calculate dust mixing ratio tendency at first model level.
  airmas = dz1 * rhoa ! kg/m2 - MJT suggestion
  a = dsrc / airmas
  duste = duste + dsrc ! MJT suggestion

! Calculate turbulent dry deposition at surface
! Use the tau-1 value of dust m.r. for now, but may modify this...

! Use full layer thickness for CSIRO model (should be correct if Vt is relative to mid-layer)
  veff = Vt*(wg+(1.-wg)*exp(-max( 0., w10m-u_ts0 )))
  b = Veff / dz1

! Update mixing ratio
! Write in form dx/dt = a - bx (a = source term, b = drydep term)
  xold = xtg(1:ifull,1,n+itracdu-1)
  
  nstep=int(tdt/120.01)+1
  ddt=tdt/real(nstep)
  do ii=1,nstep
    xtg(1:ifull,1,n+itracdu-1) = (xtg(1:ifull,1,n+itracdu-1)*(1.-0.5*b*ddt)+a*ddt)/(1.+0.5*b*ddt)
    xtg(1:ifull,1,n+itracdu-1) = max( 0., xtg(1:ifull,1,n+itracdu-1) )
  end do

  xtendd = (xtg(1:ifull,1,n+itracdu-1)-xold)/tdt - a
  dustdd = dustdd - xtendd*airmas ! Diagnostic
end do

return
end subroutine dustem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dustage

!subroutine dustage(tdt,rhg) !Inputs
!
!implicit none
!
!!     Inputs:
!real tdt              !Leapfrog timestep
!real rhg(ifull,kl)    !RH (fraction 0 to 1)
!
!! Local work arrays and variables
!real rrate(ifull,kl)
!real vvso2,rk,xd,dx
!integer k,mg,nt
!
!! Start code : ----------------------------------------------------------
!
!! Reference for this scheme is
!! Fan et al. (2004): Impact of air pollution on wet deposition of mineral dust aerosols, GRL 31,
!! L02104, doi:10.1029/2003GL018501.
!
!! Loop over 4 dust size bins
!! Reaction rate is a linear function of SO2 concentration
!
!do k=1,kl
!  do mg=1,ifull
!    vvso2 = xtg(mg,k,itracso2)*29./64. !Vol. mixing ratio of SO2
!    rk = 0.01 * max (0.1, dim(rhg(mg,k),0.5))
!    xd = 0.
!    do nt = itracdu, itracdu+ndust-1
!      xd = xd + xtg(mg,k,nt)
!    enddo
!    rrate(mg,k) = rk * vvso2 / max (1.e-12, xd)
!  enddo
!enddo
!
!! Reduce hydrophobic dust and increase hydrophilic dust
!do nt = itracdu, itracdu+ndust-1
!  do k=1,kl
!    do mg=1,ifull
!      dx = min (xtg(mg,k,nt), rrate(mg,k)*tdt*xtg(mg,k,nt))
!      xtg(mg,k,nt) = xtg(mg,k,nt) - dx
!      xtg(mg,k,nt+ndust) = xtg(mg,k,nt+ndust) + dx
!    enddo
!  enddo
!enddo
!
!return
!end subroutine dustage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A simple diagnostic treatment of seasalt aerosol (LDR 3/02)

subroutine seasalt(land,sicef,zmid,pblh,v10m) !Inputs

implicit none

! Argument list
logical, dimension(ifull), intent(in) :: land  !True for land points
real, dimension(ifull), intent(in) :: sicef    !Sea-ice fraction
real, dimension(ifull,kl), intent(in) :: zmid  !Height of full level (m)
real, dimension(ifull), intent(in) :: pblh     !PBL height (m)
real, dimension(ifull), intent(in) :: v10m     !10m windpseed, including effect of sub-grid gustiness (m/s)

real Veff

integer k,mg

! Calculate number and mass concentration of seasalt within the marine BL.
! Set seasalt conc. to zero elsewhere.
! Height of BL taken from ncarpbl scheme, so needs this turned on (although we 
! set it to 2000m in hvertmx if ncarpbl=F)
! The first mode is the "film-drop" mode, and the second is the "jet-drop" mode.
! Smaller number mode radii are given by Nilsson et al. (2001) cf. O'Dowd's.
!
! References: O'Dowd et al. (1997) Atmos. Environ. 31, 73-80
!             Jones et al. (2001) JGR 106, 20293-20310.
!             Nilsson et al. (2001) JGR 106, 32139-32154.



! Number-to-mass conversion factors are derived from the parameters of the two
! lognormal modes given by Nilsson, together with rhosalt=2.0e3 kg/m3.
do mg=1,ifull
  if (.not.land(mg)) then
    ! Jones et al. give different windspeed relations for v10m < 2, 2 <= v10m <= 17.5,
    ! and v10m > 17.5, but let's apply the middle one everywhere, since a min. windspeed 
    ! of 2 m/s seems reasonable, and the model gives few points with v10m > 17.5.
    Veff=max(2.,v10m(mg))

    ssn(mg,1,1)=10.**(0.0950*Veff+6.2830) 
    ssn(mg,1,2)=10.**(0.0422*Veff+5.7122)  
    do k=2,kl
      if (zmid(mg,k)<pblh(mg)) then
        ssn(mg,k,1)=ssn(mg,1,1)
        ssn(mg,k,2)=ssn(mg,1,2)
      else
        ssn(mg,k,1)=1.e7*exp(-zmid(mg,k)/3000.)
        ssn(mg,k,2)=1.e6*exp(-zmid(mg,k)/1500.)
      end if
    end do
  else
    ssn(mg,:,1)=0.
    ssn(mg,:,2)=0.
  end if
end do

! Reduce over sea ice...
do k=1,kl
  ssn(:,k,1)=(1.-sicef(:))*ssn(:,k,1)
  ssn(:,k,2)=(1.-sicef(:))*ssn(:,k,2)
enddo

! These relations give ssm in kg/m3 based on ssn in m^{-3}...
! Using the size distributions from Nillson et al.
!ssm(:,:,1)=5.3e-17*ssn(:,:,1) !number mode radius = 0.1 um, sd=2
!ssm(:,:,2)=9.1e-15*ssn(:,:,2) !number mode radius = 0.5 um, sd=2

! Using the size distributions as assumed by Herzog in the radiation scheme (dry sea salt)
!ssm(:,:,1)=2.64e-18*ssn(:,:,1) !number mode radius = 0.035 um, sd=1.92, rho=2.165 g/cm3
!ssm(:,:,2)=1.38e-15*ssn(:,:,2) !number mode radius = 0.35 um, sd=1.7, rho=2.165

! Using the size distributions as assumed by Herzog in the radiation scheme (dry sea salt)
! Use ssn(3) to hold diagnostic of mass conc. for now in kg/m3
!ssn(:,:,3)=2.64e-18*ssn(:,:,1)+1.38e-15*ssn(:,:,2)
      
return
end subroutine seasalt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cloud droplet concentration

subroutine cldrop(istart,imax,cdn,rhoa,convmode)

implicit none

integer, intent(in) :: istart,imax
integer k,is,ie
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(out) :: cdn
real, dimension(imax,kl) :: xtgso4,xtgbc,xtgoc
real, dimension(imax) :: so4_n,cphil_n,Atot
logical, intent(in) :: convmode

is=istart
ie=istart+imax-1

if (convmode) then
  ! total grid-box
  xtgso4 = xtg(is:ie,:,itracso2+1)
  xtgbc  = xtg(is:ie,:,itracbc+1)
  xtgoc  = xtg(is:ie,:,itracoc+1)
else
  ! outside convective fraction of grid-box
  xtgso4 = xtosav(is:ie,:,itracso2+1)
  xtgbc  = xtosav(is:ie,:,itracbc+1)
  xtgoc  = xtosav(is:ie,:,itracoc+1)
end if

do k=1,kl
  ! Factor of 132.14/32.06 converts from sulfur to ammmonium sulfate
  ! 1.24e17 converts from mass (kg/m3) to number concentration (/m3) for dist'n 
  ! from Penner et al (1998).
  ! 1.69e17 converts from mass (kg/m3) to number concentration (/m3) for dist'n 
  ! from IPCC (2001), Table 5.1, as used by Minghuai Wang for lookup optical properties.
  so4_n = 1.24e17 * (132.14/32.06) * rhoa(:,k) * xtgso4(:,k)

  ! Factor of 1.3 converts from OC to organic matter (OM) 
  ! 1.25e17 converts from hydrophilic mass (kg/m3) to number concentration (/m3) for
  ! Hardiman lognormal distribution for carbonaceous aerosols (Penner et al, 1998).
  ! 1.21e17 converts from hydrophilic mass (kg/m3) to number concentration (/m3) for
  ! biomass regional haze distribution from IPCC (2001), Table 5.1. Using rho_a=1250 kg/m3.

  ! Following line counts Aitken mode as well as accumulation mode carb aerosols
  cphil_n = 2.30e17 * rhoa(:,k) * (xtgbc(:,k)+1.3*xtgoc(:,k))

  ! The dust particles are the accumulation mode only (80.2% of the hydrophilic 
  ! "small dust" particles)
  !dust_n(:) = 0.
  !aero_n(:) = max (10.e6, so4_n(:) + cphil_n(:) + ssn(is:ie,k,1) + ssn(is:ie,k,2) + dust_n(:))

  ! Jones et al., modified to account for hydrophilic carb aerosols as well
  Atot = so4_n + cphil_n + ssn(is:ie,k,1) + ssn(is:ie,k,2)
  cdn(:,k)=max(10.e6, 375.e6*(1.-exp(-2.5e-9*Atot)))

enddo

return
end subroutine cldrop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aerosol scavenging fraction for convective clouds

subroutine convscav(fscav,xpkp1,xpold,tt,xs,rho)

implicit none

real, dimension(ifull,naero), intent(out) :: fscav ! scavenging fraction
real, dimension(ifull), intent(in) :: xpkp1 ! cloud liquid water after precipitation
real, dimension(ifull), intent(in) :: xpold ! cloud liquid water before precipitation
real, dimension(ifull), intent(in) :: tt    ! parcel temperature
real, dimension(ifull), intent(in) :: xs    ! xtg(:,k,3) = so4
real, dimension(ifull), intent(in) :: rho   ! air density
real, dimension(ifull) :: f_so2
logical, dimension(ifull) :: bwkp1 
! In-cloud scavenging efficiency for liquid and frozen convective clouds follows.
! Hard-coded for 3 sulfur variables, 4 carbonaceous, 4 mineral dust.
! Note that value for SO2 (index 2) is overwritten by Henry coefficient f_so2 below.
! These ones are for 3 SULF, 4 CARB and 4 or 8 DUST (and include dummy variables at end)
real, parameter, dimension(naero) :: scav_effl = (/0.00,1.00,0.90,0.00,0.30,0.00,0.30,0.05,0.05,0.05,0.05/) ! liquid
real, parameter, dimension(naero) :: scav_effi = (/0.00,0.00,0.00,0.05,0.00,0.05,0.00,0.05,0.05,0.05,0.05/) ! ice
real, dimension(ifull) :: scav_eff
real, dimension(ifull) :: zqtp1,ze2,ze3,zfac,zso4l,zso2l,zqhp
real, dimension(ifull) :: zza,zzb,zzp,zzq,zzp2,zhp,zqhr,zheneff,p_so2
integer nt

bwkp1=tt>=ticeu ! condensate in parcel is liquid (true) or ice (false)

! CALCULATE THE SOLUBILITY OF SO2
! TOTAL SULFATE  IS ONLY USED TO CALCULATE THE PH OF CLOUD WATER
where ( xpold>1.e-20 .and. bwkp1 )
  ZQTP1=1./tt-1./298.
  ZE2=1.23*EXP(3020.*ZQTP1)
  ZE3=1.2E-02*EXP(2010.*ZQTP1)
  ZFAC=1000./(xpold*32.064)
  ZSO4L=xs*ZFAC
  ZSO4L=AMAX1(ZSO4L,0.)
  ZSO2L=xs*ZFAC
  ZSO2L=AMAX1(ZSO2L,0.)
  ZZA=ZE2*8.2E-02*tt*xpold*rho*1.E-03
  ZZB=2.5E-06+ZSO4L
  ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
  ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
  ZZP=0.5*ZZP
  ZZP2=ZZP*ZZP
  ZHP=-ZZP+SQRT(max(ZZP2-ZZQ,0.))
  ZQHP=1./ZHP
  ZHENEFF=1.+ZE3*ZQHP
  P_SO2=ZZA*ZHENEFF
  F_SO2=P_SO2/(1.+P_SO2)
  F_SO2=min(max(0.,F_SO2),1.)
elsewhere
  f_so2=0.
end where

do nt=1,naero
  where ( bwkp1 .and. nt==ITRACSO2 )
    scav_eff=f_so2(:)
  elsewhere ( bwkp1 )
    scav_eff=scav_effl(nt)
  elsewhere
    scav_eff=scav_effi(nt)
  end where
  ! Wet deposition scavenging fraction
  fscav(:,nt)=scav_eff*min(max(xpold-xpkp1,0.)/max(xpold,1.E-20),1.)
end do

return
end subroutine convscav

!     DEFINE FUNCTION FOR CHANGING THE UNITS
!     FROM MASS-MIXING RATIO TO MOLECULES PER CM**3 AND VICE VERSA

function xtoc(x,y) result(ans)
implicit none
real, intent(in) :: x, y
real ans
ans=X*6.022E+20/Y
end function xtoc

function ctox(x,y) result(ans)
implicit none
real, intent(in) :: x, y
real ans
ans=Y/(6.022E+20*X)
end function ctox

function zfarr(zk,zh,ztpq) result(ans)
implicit none
real, intent(in) :: zk, zh, ztpq
real ans
!   X = DENSITY OF AIR, Y = MOL WEIGHT IN GRAMM
ans=ZK*EXP(ZH*ZTPQ)
end function zfarr

end module aerosolldr
