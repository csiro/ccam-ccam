! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! This module calculates the turblent kinetic energy and mixing for the boundary layer based on Hurley 2007
! (eddy dissipation) and Angevine et al 2010 (mass flux).  Specifically, this version is modified for
! clouds and saturated air following Durran and Klemp JAS 1982.

! Usual procedure

! call tkeinit
! ...
! do t=1,tmax
!   ...
!   (host horizontal advection routines for tke and eps)
!   ...
!   shear=...       ! Calculate horizontal shear for tkemix
!   call tkemix     ! Updates TKE and eps source terms, updates theta and qg non-local terms and outputs kdiff
!   ...
!   (host vertical advection for TKE, eps, theta and mixing ratio)
!   ...
! end do
! ...
    
module tkeeps

implicit none

private
public tkeinit,tkemix,tkeend,tke,eps,shear
public cm0,ce0,ce1,ce2,ce3,cq,be,ent0,ent1,entc0,ezmin,dtrc0
public m0,b1,b2,qcmf
public buoymeth,maxdts,mintke,mineps,minl,maxl,stabmeth
public tke_umin,tkemeth
#ifdef offline
public wthl,wqv,wql,wqf
public mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
public ents,dtrs
#endif

integer, save :: ifull,iextra,kl
integer, save, public :: nb
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
#ifdef offline
real, dimension(:,:), allocatable, save :: wthl,wqv,wql,wqf
real, dimension(:,:), allocatable, save :: mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
real, dimension(:,:), allocatable, save :: u,v,ents,dtrs
#endif

! model ED constants
real, save :: cm0     = 0.09   ! Hurley (2007) 0.09, Duynkerke (1988) 0.03, Duynkerke (1987) 0.09
real, save :: ce0     = 0.69   ! Hurley (2007) 0.69, Duynkerke (1988) 0.42, Duynkerke (1987) 0.77
real, save :: ce1     = 1.46
real, save :: ce2     = 1.83
real, save :: ce3     = 0.45   ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, save :: cq      = 2.5    ! Adjustment to ED in absence of MF
! model MF constants
real, save :: be      = 0.1    ! Surface boundary condition (Hurley (2007) 1., Soares et al (2004) 0.3)
real, save :: ent0    = 0.25   ! Entrainment constant (Controls height of boundary layer)
real, save :: ent1    = 0.25
real, save :: ezmin   = 100.   ! Limits entrainment at plume top
real, save :: entc0   = 2.e-3  ! Saturated entrainment constant for mass flux
real, save :: dtrc0   = 3.e-3  ! Saturated detrainment constant for mass flux
real, save :: m0      = 0.1    ! Mass flux amplitude constant
real, save :: b1      = 2.     ! Updraft entrainment coeff (Soares et al (2004) 1., Siebesma et al (2003) 2.)
real, save :: b2      = 1./3.  ! Updraft buoyancy coeff (Soares et al (2004) 2., Siebesma et al (2003) 1./3.)
real, save :: qcmf    = 1.e-4  ! Critical mixing ratio of liquid water before autoconversion
! numerical constants
integer, save :: buoymeth = 1        ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0        ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: tkemeth  = 1        ! Method for TKE calculation (0=D&K84, 1=Hurley)
real, save :: maxdts      = 120.     ! max timestep for split
real, save :: mintke      = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps      = 1.E-11   ! min value for eps (1.0e-6 in TAPM)
real, save :: minl        = 1.       ! min value for L   (5. in TAPM)
real, save :: maxl        = 1000.    ! max value for L   (500. in TAPM)
real, save :: tke_umin    = 0.1      ! minimum wind speed (m/s) for drag calculation

! physical constants
real, parameter :: grav  = 9.80616    ! (m s^-2)
real, parameter :: lv    = 2.5104e6   ! (J kg^-1)
real, parameter :: lf    = 3.36e5     ! (J kg^-1)
real, parameter :: ls    = lv+lf      ! (J kg^-1)
real, parameter :: rd    = 287.04
real, parameter :: rv    = 461.5
real, parameter :: cp    = 1004.64    ! (J kg^-1 K^-1)
real, parameter :: vkar  = 0.4
real, parameter :: pi    = 3.14159265

! MOST constants
real, parameter :: a_1   = 1.
real, parameter :: b_1   = 2./3.
real, parameter :: c_1   = 5.
real, parameter :: d_1   = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

real, dimension(0:220), save :: table

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifullin,iextrain,klin,diag,nbin)

implicit none

integer, intent(in) :: ifullin,iextrain,klin,diag,nbin
real cm34

if (diag>0) write(6,*) "Initialise TKE-eps scheme"

nb=nbin
ifull=ifullin
iextra=iextrain
kl=klin

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(shear(ifull,kl))

cm34=cm0**0.75
tke=mintke
eps=mineps
shear=0.

#ifdef offline
allocate(wthl(ifull,kl),wqv(ifull,kl),wql(ifull,kl),wqf(ifull,kl))
allocate(mf(ifull,kl),w_up(ifull,kl),tl_up(ifull,kl),qv_up(ifull,kl))
allocate(ql_up(ifull,kl),qf_up(ifull,kl),cf_up(ifull,kl))
allocate(ents(ifull,kl),dtrs(ifull,kl))
wthl=0.
wqv=0.
wql=0.
wqf=0.
mf=0.
w_up=0.
tl_up=0.
qv_up=0.
ql_up=0.
qf_up=0.
cf_up=0.
ents=0.
dtrs=0.
#endif
ifull=ifullin/nb

table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)                                                 !70C

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

#ifdef scm
subroutine tkemix(kmo,theta,qvg,qlg,qfg,cfrac,uo,vo,zi,fg,eg,ps,zom,zz,zzh,sig,rhos,  &
                  dt,qgmin,mode,diag,naero,aero,cgmap,wthflux,wqvflux,uwflux,vwflux,mfout)
#else
subroutine tkemix(kmo,theta,qvg,qlg,qfg,cfrac,uo,vo,zi,fg,eg,ps,zom,zz,zzh,sig,rhos, &
                  dt,qgmin,mode,diag,naero,aero,cgmap,tkel,epsl,shearl)
#endif

implicit none

integer, intent(in) :: diag,mode,naero
integer k, i, j, ktopmax
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(:,:,:), intent(inout) :: aero
real, dimension(:,:), intent(inout) :: theta,cfrac,uo,vo,tkel,epsl,shearl
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl), intent(in) :: zz,zzh
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: fg,eg,ps,zom,rhos,cgmap
real, dimension(kl), intent(in) :: sig
real, dimension(ifull,kl,naero) :: arup
real, dimension(ifull,kl) :: km,thetav,thetal,qsat
real, dimension(ifull,kl) :: qsatc,qgnc,ff
real, dimension(ifull,kl) :: thetalhl,thetavhl
real, dimension(ifull,kl) :: quhl,qshl,qlhl,qfhl
real, dimension(ifull,kl) :: bb,cc,dd,rr
real, dimension(ifull,kl) :: rhoa,rhoahl
real, dimension(ifull,kl) :: qtot,qthl
real, dimension(ifull,kl) :: tlup,qvup,qlup,qfup
real, dimension(ifull,kl) :: cfup,mflx
real, dimension(ifull,2:kl) :: idzm
real, dimension(ifull,1:kl-1) :: idzp
real, dimension(ifull,2:kl) :: aa,qq,pps,ppt,ppb
real, dimension(ifull,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,kl-1) :: fzzh
real, dimension(ifull) :: wt0,wq0,wtv0
real, dimension(ifull) :: wstar,z_on_l,phim
real, dimension(ifull) :: tff,tgg,tempc,thetac,pres,temp
real, dimension(ifull) :: cdrag,umag,ustar
real, dimension(ifull) :: tempv,rvar,bvf,dc,mc,fc
real, dimension(ifull) :: tbb,tcc,tqq
real, dimension(ifull) :: avearray
real, dimension(kl) :: sigkap
real, dimension(kl) :: w2up,nn,dqdash,qupsat
real, dimension(kl) :: qtup,tvup,thup
real, dimension(1) :: templ
real xp, as, bs, cs, cm12, cm34, qcup
real dzht, ent
real ddts, upf
real lx, tempd, fice, qxup, dqsdt, al
real sigqtup, rng
logical, dimension(ifull,kl) :: lta

#ifdef offline
real dtr
#endif

#ifdef scm
real, dimension(ifull,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(ifull,kl-1), intent(out) :: mfout
real, dimension(ifull,kl) :: wthlflux, wqlflux
real, dimension(ifull,kl) :: wqfflux
#endif

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))

if ( diag>0 ) write(6,*) "Update PBL mixing with TKE-eps + MF turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

do k = 1,kl
  ! Impose limits on tke and eps after being advected by the host model
  tkel(1:ifull,k) = max(tkel(1:ifull,k), mintke)
  tff(1:ifull) = cm34*tkel(1:ifull,k)*sqrt(tkel(1:ifull,k))
  epsl(1:ifull,k) = min(epsl(1:ifull,k), tff(1:ifull)/minl)
  epsl(1:ifull,k) = max(epsl(1:ifull,k), tff(1:ifull)/maxl, mineps)

  ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
  sigkap(k) = sig(k)**(-rd/cp)
  pres(1:ifull) = ps(1:ifull)*sig(k) ! pressure
  ! Density must be updated when dz is updated so that rho*dz is conserved
  thetav(1:ifull,k) = theta(1:ifull,k)*(1.+0.61*qvg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k))
  rhoa(1:ifull,k) = sigkap(k)*pres(1:ifull)/(rd*thetav(1:ifull,k))

  ! Transform to thetal as it is the conserved variable
  thetal(1:ifull,k) = theta(1:ifull,k) - sigkap(k)*(lv*qlg(1:ifull,k)+ls*qfg(1:ifull,k))/cp
  
  ! Calculate first approximation to diffusion coeffs
  km(1:ifull,k) = cm0*tkel(1:ifull,k)*tkel(1:ifull,k)/epsl(1:ifull,k)
end do

! Calculate surface fluxes
wt0(1:ifull) = fg(1:ifull)/(rhos(1:ifull)*cp)  ! theta flux
wq0(1:ifull) = eg(1:ifull)/(rhos(1:ifull)*lv)  ! qtot flux (=qv flux)

do k = 1,kl-1
  ! Fraction for interpolation
  fzzh(1:ifull,k) = (zzh(1:ifull,k)-zz(1:ifull,k))/(zz(1:ifull,k+1)-zz(1:ifull,k))
  ! Calculate dz at half levels
  dz_hl(1:ifull,k) = zz(1:ifull,k+1) - zz(1:ifull,k)
end do

! Calculate dz at full levels
dz_fl(1:ifull,1)    = zzh(1:ifull,1)
dz_fl(1:ifull,2:kl) = zzh(1:ifull,2:kl) - zzh(1:ifull,1:kl-1)

! Calculate shear term on full levels (see hordifg.f for calculation of horizontal shear)
pps(1:ifull,2:kl-1) = km(1:ifull,2:kl-1)*shearl(1:ifull,2:kl-1)

! set top BC for TKE-eps source terms
pps(1:ifull,kl) = 0.
ppb(1:ifull,kl) = 0.
ppt(1:ifull,kl) = 0.

! interpolate diffusion coeff and air density to half levels
call updatekmo(kmo,   km,  fzzh)
call updatekmo(rhoahl,rhoa,fzzh)
idzm(1:ifull,2:kl)   = rhoahl(1:ifull,1:kl-1)/(rhoa(1:ifull,2:kl)*dz_fl(1:ifull,2:kl))
idzp(1:ifull,1:kl-1) = rhoahl(1:ifull,1:kl-1)/(rhoa(1:ifull,1:kl-1)*dz_fl(1:ifull,1:kl-1))

! Main loop to prevent time splitting errors
mcount = int(dt/(maxdts+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  ! Set-up thermodynamic variables temp, theta_v and surface fluxes
  do k = 1,kl
    temp(1:ifull) = theta(1:ifull,k)/sigkap(k)
    ! calculate saturated air mixing ratio
    pres(1:ifull) = ps(1:ifull)*sig(k)
    call getqsat(qsat(:,k),temp(:),pres(:))
    thetav(1:ifull,k) = theta(1:ifull,k)*(1.+0.61*qvg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k))
    qtot(1:ifull,k) = qvg(1:ifull,k) + qlg(1:ifull,k) + qfg(1:ifull,k)
  end do

  ! Update momentum flux
  wtv0(1:ifull) = wt0(1:ifull) + theta(1:ifull,1)*0.61*wq0(1:ifull) ! thetav flux
  umag(1:ifull) = sqrt(max( uo(1:ifull,1)*uo(1:ifull,1)+vo(1:ifull,1)*vo(1:ifull,1), tke_umin**2 ))
  call dyerhicks(cdrag(1:ifull),wtv0(1:ifull),zom(1:ifull),umag(1:ifull),thetav(1:ifull,1),zz(1:ifull,1))
  ustar(1:ifull) = sqrt(cdrag(1:ifull))*umag(1:ifull)
  
  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (1:ifull.e., volume conserving)
  mflx(1:ifull,:) = 0.
  tlup(1:ifull,:) = thetal(1:ifull,:)
  qvup(1:ifull,:) = qvg(1:ifull,:)
  qlup(1:ifull,:) = qlg(1:ifull,:)
  qfup(1:ifull,:) = qfg(1:ifull,:)
  cfup(1:ifull,:) = cfrac(1:ifull,:)
  if ( naero>0 ) then
    arup(1:ifull,:,:) = aero(1:ifull,:,:)
  end if

#ifdef scm
  mfout(1:ifull,:)=0.
#endif
#ifdef offline
  mf(1:ifull,:)=0.
  w_up(1:ifull,:)=0.
  tl_up(1:ifull,:)=thetal(1:ifull,:)
  qv_up(1:ifull,:)=qvg(1:ifull,:)
  ql_up(1:ifull,:)=qlg(1:ifull,:)
  qf_up(1:ifull,:)=qfg(1:ifull,:)
  cf_up(1:ifull,:)=0.
  ents(1:ifull,:)=0.
  dtrs(1:ifull,:)=0.
#endif

  if ( mode/=1 ) then ! mass flux

    do i = 1,ifull
      wstar(i) = (grav*zi(i)*max(wtv0(i),0.)/thetav(i,1))**(1./3.)
      if ( wtv0(i)>0. ) then ! unstable
        ! Initialise updraft
        tkel(i,1)=cm12*ustar(i)*ustar(i)+ce3*wstar(i)*wstar(i)
        tkel(i,1)=max(tkel(i,1),mintke)
        ktopmax=0
        w2up=0.
        dzht=zz(i,1)
        ! Entrainment and detrainment rates
        ent=entfn(zz(i,1),zi(i))
          
        ! first level -----------------
        ! initial thermodynamic state
        ! split qtot into components (conservation of thetal and qtot is maintained)
        tlup(i,1)=thetal(i,1)+be*wt0(i)/sqrt(tkel(i,1))       ! Hurley 2007
        qvup(i,1)=qvg(i,1)   +be*wq0(i)/sqrt(tkel(i,1))       ! Hurley 2007
        qlup(i,1)=qlg(i,1)
        qfup(i,1)=qfg(i,1)
        ! diagnose thermodynamic variables assuming no condensation
        qtup(1)=qvup(i,1)+qlup(i,1)+qfup(i,1)                    ! qtot,up
        ! state of plume after evaporation
        qvup(i,1)=qtup(1)
        qlup(i,1)=0.
        qfup(i,1)=0.
        thup(1)=tlup(i,1) !+sigkap(1)*(lv*qlup(1)+ls*qfup(1))/cp ! theta,up
        tvup(1)=thup(1)+theta(i,1)*0.61*qtup(1)                  ! thetav,up
        templ(1)=tlup(i,1)/sigkap(1)                             ! templ,up
        ! update updraft velocity and mass flux
        nn(1)  =grav*be*wtv0(i)/(thetav(i,1)*sqrt(tkel(i,1))) ! Hurley 2007
        w2up(1)=2.*dzht*b2*nn(1)/(1.+2.*dzht*b1*ent)         ! Hurley 2007
        ! estimate variance of qtup in updraft
        pres(i) = ps(i)*sig(1)
        call getqsat(qupsat(1:1),templ(1:1),pres(i:i))
        sigqtup=1.E-5
        rng=sqrt(6.)*sigqtup               ! variance of triangle distribution
        dqdash(1)=(qtup(1)-qupsat(1))/rng  ! scaled variance
        dqdash(1)=min(dqdash(1),-1.)
        cfup(i,1) = 0.
        
        ! updraft with condensation
        do k=2,kl
          dzht=dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          ent=entfn(zz(i,k),zi(i))
          ! update thermodynamics of plume
          ! split qtot into components (conservation is maintained)
          tlup(i,k)=(tlup(i,k-1)+dzht*ent*thetal(i,k))/(1.+dzht*ent)
          qvup(i,k)=(qvup(i,k-1)+dzht*ent*qvg(i,k)   )/(1.+dzht*ent)
          qlup(i,k)=(qlup(i,k-1)+dzht*ent*qlg(i,k)   )/(1.+dzht*ent)
          qfup(i,k)=(qfup(i,k-1)+dzht*ent*qfg(i,k)   )/(1.+dzht*ent)
          ! calculate conserved variables
          qtup(k)=qvup(i,k)+qlup(i,k)+qfup(i,k)                ! qtot,up
          ! estimate air temperature
          templ(1)=tlup(i,k)/sigkap(k)                         ! templ,up
          pres(i) = ps(i)*sig(k)
          call getqsat(qupsat(k:k),templ(1:1),pres(i:i))
          ! estimate variance of qtup in updraft (following Hurley and TAPM)
          sigqtup=sqrt(max(1.E-10, 1.6*tkel(i,k)/epsl(i,k)*cq*km(i,k)*((qtup(k)-qtup(k-1))/dzht)**2))
          ! MJT condensation scheme -  follow Smith 1990 and assume
          ! triangle distribution for qtup.  The average qvup is qxup
          ! after accounting for saturation
          rng=sqrt(6.)*sigqtup               ! variance of triangle distribution
          dqdash(k)=(qtup(k)-qupsat(k))/rng  ! scaled variance
          if (dqdash(k)<-1.) then
            ! gridbox all unsaturated
            qxup=qtup(k)
            cfup(i,k)=0.
          else if (dqdash(k)<0.) then
            ! gridbox minority saturated
            qxup=qtup(k)+0.5*rng*(-1./3.-dqdash(k)-dqdash(k)**2-1./3.*dqdash(k)**3)
            cfup(i,k)=0.5*(dqdash(k)+1.)**2
          else if (dqdash(k)<1.) then
            ! gridbox majority saturated
            qxup=qtup(k)+0.5*rng*(-1./3.-dqdash(k)-dqdash(k)**2+1./3.*dqdash(k)**3)
            cfup(i,k)=1.-0.5*(dqdash(k)-1.)**2
          else
            ! gridbox all saturated
            qxup=qupsat(k)
            cfup(i,k)=1.
          end if
          thup(k)=tlup(i,k)+sigkap(k)*(lv*qlup(i,k)+ls*qfup(i,k))/cp ! theta,up before redistribution
          tempd  =thup(k)/sigkap(k)                                  ! temp,up before redistribution
          fice=min(max(273.16-tempd,0.),40.)/40. ! approximate ice fraction based on temperature (not templ)
          lx=lv+lf*fice
          dqsdt=qupsat(k)*lx/(rv*templ(1)*templ(1))
          al=cp/(cp+lx*dqsdt)
          qcup=max(al*(qtup(k)-qxup), 0.)                        ! qcondensate,up after redistribution
          qcup=min(qcup, qcmf)                                   ! limit condensation with autoconversion
          qxup=qtup(k)-qcup                                      ! qv,up after redistribution
          thup(k)=tlup(i,k)+sigkap(k)*qcup*lx/cp                 ! theta,up after redistribution
          tvup(k)=thup(k)+theta(i,k)*(1.61*qxup-qtup(k))         ! thetav,up after redistribution
          ! state of plume after redistribution
          qvup(i,k)=qxup                                         ! qv,up after redistribution
          qlup(i,k)=qcup*(1.-fice)                               ! ql,up after redistribution
          qfup(i,k)=qcup*fice                                    ! qf,up after redistribution
          ! calculate buoyancy
          nn(k)  =grav*(tvup(k)-thetav(i,k))/thetav(i,k)
          ! update updraft velocity
          w2up(k)=(w2up(k-1)+2.*dzht*b2*nn(k))/(1.+2.*dzht*b1*ent)
          ! test if maximum plume height is reached
          if (w2up(k)<=0.) then
            as=min(2.*b2*(nn(k)-nn(k-1))/dzht,-1.e-20)
            bs=2.*b2*nn(k-1)
            cs=w2up(k-1)
            xp=0.5*(-bs-sqrt(max(bs*bs-4.*as*cs,0.)))/as
            xp=min(max(xp,0.),dzht)
            zi(i)=xp+zz(i,k-1)
            ktopmax=max(ktopmax,k-1)
            exit
          end if
        end do
          
        wstar(i)=(grav*zi(i)*max(wtv0(i),0.)/thetav(i,1))**(1./3.)
          
        ! update mass flux
        mflx(i,1) = m0*sqrt(max(w2up(1), 0.))
        do k = 2,ktopmax
          dzht = dz_hl(i,k-1)
          upf = mflx(i,k-1)/sqrt(max(w2up(k-1), 1.e-8))
          mflx(i,k) = (1.-cfup(i,k))*m0*sqrt(max(w2up(k), 0.))         &
                    + cfup(i,k)*mflx(i,k-1)/(1.+dzht*(dtrc0-entc0))
          mflx(i,k) = min( mflx(i,k), upf*sqrt(max(w2up(k), 0.)) )
        end do
        
#ifdef offline
        do k=1,ktopmax
          mf(i,k)=mflx(i,k)
          w_up(i,k)=sqrt(w2up(k))
          tl_up(i,k)=tlup(i,k)
          qv_up(i,k)=qvup(i,k)
          ql_up(i,k)=qlup(i,k)
          qf_up(i,k)=qfup(i,k)
          cf_up(i,k)=cfup(i,k)*min(mflx(i,k)/sqrt(w2up(k)),1.)
        end do

        ! Entrainment and detrainment rates
        dzht = zz(i,1)
        ent = entfn(zz(i,1),zi(i))
        dtr = -1./dzht + ent
        dtr = max( dtr, 0. )
        ents(i,1) = ent
        dtrs(i,1) = dtr
        do k = 2,ktopmax
          dzht = dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          ent = entfn(zz(i,k),zi(i))
          dtr = mflx(i,k-1)/(mflx(i,k)*dzht) - 1./dzht + ent
          dtr = max( dtr, 0. )
          ents(i,k) = ent
          dtrs(i,k) = dtr
        end do
#endif

        ! update reamining scalars which are not used in the iterative loop
        do j = 1,naero
          arup(i,1,j) = aero(i,1,j)
          do k = 2,ktopmax
            dzht = dz_hl(i,k-1)
            ent = entfn(zz(i,k),zi(i))
            arup(i,k,j) = (arup(i,k-1,j)+dzht*ent*aero(i,k,j))/(1.+dzht*ent)
          end do
        end do
        
      else                   ! stable
        !wpv_flux is calculated at half levels
        !wpv_flux(1)=-kmo(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1)
        !do k=2,kl-1
        !  wpv_flux(k)=-kmo(i,k)*(thetav(i,k+1)-thetav(i,k))/dz_hl(i,k)
        !  if (wpv_flux(k)*wpv_flux(1)<0.) then
        !    xp=(0.05*wpv_flux(1)-wpv_flux(k-1))/(wpv_flux(k)-wpv_flux(k-1))
        !    xp=min(max(xp,0.),1.)
        !    zi(i)=zzh(i,k-1)+xp*(zzh(i,k)-zzh(i,k-1))
        !    exit
        !  else if (abs(wpv_flux(k))<0.05*abs(wpv_flux(1))) then
        !    xp=(0.05*abs(wpv_flux(1))-abs(wpv_flux(k-1)))/(abs(wpv_flux(k))-abs(wpv_flux(k-1)))
        !    xp=min(max(xp,0.),1.)
        !    zi(i)=zzh(i,k-1)+xp*(zzh(i,k)-zzh(i,k-1))
        !    exit
        !  end if
        !end do
        zi(i)    = zz(i,1) ! MJT suggestion
      end if
    end do
  
  else
    wstar(1:ifull) = (grav*zi(1:ifull)*max(wtv0(1:ifull),0.)/thetav(1:ifull,1))**(1./3.)   
  end if

  
  ! turn off MF term if small grid spacing
  do k = 1,kl
    mflx(1:ifull,k) = mflx(1:ifull,k)*cgmap(1:ifull)
  end do
#ifdef scm  
  mfout(1:ifull,1:kl-1) = mflx(1:ifull,1:kl-1)*(1.-fzzh(1:ifull,1:kl-1)) &
                        + mflx(1:ifull,2:kl)*fzzh(1:ifull,1:kl-1)
#endif
  
  
  ! calculate tke and eps at 1st level
  z_on_l(1:ifull)=-vkar*zz(1:ifull,1)*grav*wtv0(1:ifull)/(thetav(1:ifull,1)*max(ustar(1:ifull)*ustar(1:ifull)*ustar(1:ifull),1.E-10))
  z_on_l(1:ifull)=min(z_on_l(1:ifull),10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  select case(stabmeth)
    case(0)
      where (z_on_l(1:ifull)<0.)
        phim(1:ifull)=(1.-16.*z_on_l(1:ifull))**(-0.25)
      elsewhere
        phim(1:ifull)=1.+z_on_l(1:ifull)*(a_1+b_1*exp(-d_1*z_on_l(1:ifull))*(1.+c_1-d_1*z_on_l(1:ifull))) ! Beljarrs and Holtslag (1991)
      end where
    case(1)
      where (z_on_l(1:ifull)<0.)
        phim(1:ifull)=(1.-16.*z_on_l(1:ifull))**(-0.25)
      elsewhere (z_on_l(1:ifull)<=0.4)
        phim(1:ifull)=1.+z_on_l(1:ifull)*(a_1+b_1*exp(-d_1*z_on_l(1:ifull))*(1.+c_1-d_1*z_on_l(1:ifull))) ! Beljarrs and Holtslag (1991)
      elsewhere
        phim(1:ifull)=aa1*bb1*(z_on_l(1:ifull)**bb1)*(1.+cc1/bb1*z_on_l(1:ifull)**(1.-bb1)) ! Luhar
      end where
    case default
      write(6,*) "ERROR: Invalid option for stabmeth in tkeeps ",stabmeth
      stop
  end select
  tkel(1:ifull,1)=cm12*ustar(1:ifull)*ustar(1:ifull)+ce3*wstar(1:ifull)*wstar(1:ifull)
  epsl(1:ifull,1)=ustar(1:ifull)*ustar(1:ifull)*ustar(1:ifull)*phim(1:ifull)/(vkar*zz(1:ifull,1))+grav*wtv0(1:ifull)/thetav(1:ifull,1)
  tkel(1:ifull,1)=max(tkel(1:ifull,1),mintke)
  tff(1:ifull)=cm34*tkel(1:ifull,1)*sqrt(tkel(1:ifull,1))
  epsl(1:ifull,1)=min(epsl(1:ifull,1),tff(1:ifull)/minl)
  epsl(1:ifull,1)=max(epsl(1:ifull,1),tff(1:ifull)/maxl,mineps)


  ! Update TKE and eps terms

  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tkel(1:ifull,kl)=mintke
  epsl(1:ifull,kl)=mineps
  
  ! Calculate buoyancy term
  select case(buoymeth)
    case(0) ! saturated from Durran and Klemp JAS 1982 (see also WRF)
      qsatc(1:ifull,:)=max(qsat(1:ifull,:),qvg(1:ifull,:))                                             ! assume qvg is saturated inside cloud
      ff(1:ifull,:)=qfg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                              ! inside cloud value  assuming max overlap
      dd(1:ifull,:)=qlg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                              ! inside cloud value assuming max overlap
      do k=1,kl
        tbb(1:ifull)=max(1.-cfrac(1:ifull,k),1.E-8)
        qgnc(1:ifull,k)=(qvg(1:ifull,k)-(1.-tbb(1:ifull))*qsatc(1:ifull,k))/tbb(1:ifull)                       ! outside cloud value
        qgnc(1:ifull,k)=min(max(qgnc(1:ifull,k),qgmin),qsatc(1:ifull,k))
      end do
      call updatekmo(thetalhl,thetal,fzzh)                                       ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh)                                             ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh)                                            ! inside cloud value
      call updatekmo(qlhl,dd,fzzh)                                               ! inside cloud value
      call updatekmo(qfhl,ff,fzzh)                                               ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(1:ifull,2:kl)=cfrac(1:ifull,2:kl)<=1.E-6
      do k=2,kl-1
        where(lta(1:ifull,k).and..not.lta(1:ifull,k+1))
          qlhl(1:ifull,k)=dd(1:ifull,k+1)
          qfhl(1:ifull,k)=ff(1:ifull,k+1)
        elsewhere (.not.lta(1:ifull,k).and.lta(1:ifull,k+1))
          qlhl(1:ifull,k)=dd(1:ifull,k)
          qfhl(1:ifull,k)=ff(1:ifull,k)
        end where
      end do
      do k=2,kl-1
        ! saturated
        thetac(1:ifull)=thetal(1:ifull,k)+sigkap(k)*(lv*dd(1:ifull,k)+ls*ff(1:ifull,k))/cp              ! inside cloud value
        tempc(1:ifull)=thetac(1:ifull)/sigkap(k)                                            ! inside cloud value          
        tqq(1:ifull)=(1.+lv*qsatc(1:ifull,k)/(rd*tempc(1:ifull)))/(1.+lv*lv*qsatc(1:ifull,k)/(cp*rv*tempc(1:ifull)*tempc(1:ifull)))
        tbb(1:ifull)=-grav*km(1:ifull,k)*(tqq(1:ifull)*((thetalhl(1:ifull,k)-thetalhl(1:ifull,k-1)+sigkap(k)/cp*(lv*(qlhl(1:ifull,k)-qlhl(1:ifull,k-1))  &
            +ls*(qfhl(1:ifull,k)-qfhl(1:ifull,k-1))))/thetac(1:ifull)+lv/cp*(qshl(1:ifull,k)-qshl(1:ifull,k-1))/tempc(1:ifull))              &
            -qshl(1:ifull,k)-qlhl(1:ifull,k)-qfhl(1:ifull,k)+qshl(1:ifull,k-1)+qlhl(1:ifull,k-1)+qfhl(1:ifull,k-1))/dz_fl(1:ifull,k)
        ! unsaturated
        tcc(1:ifull)=-grav*km(1:ifull,k)*(thetalhl(1:ifull,k)-thetalhl(1:ifull,k-1)+thetal(1:ifull,k)*0.61*(quhl(1:ifull,k)-quhl(1:ifull,k-1))) &
                         /(thetal(1:ifull,k)*dz_fl(1:ifull,k))
        ppb(1:ifull,k)=(1.-cfrac(1:ifull,k))*tcc(1:ifull)+cfrac(1:ifull,k)*tbb(1:ifull) ! cloud fraction weighted (e.g., Smith 1990)
      end do
      
    case(1) ! follow Marquet and Geleyn QJRMS (2012)
      call updatekmo(thetalhl,thetal,fzzh)
      call updatekmo(qthl,qtot,fzzh)
      do k=2,kl-1
        temp(1:ifull)  = theta(1:ifull,k)/sigkap(k)
        tempv(1:ifull) = thetav(1:ifull,k)/sigkap(k)
        rvar(1:ifull)=rd*tempv(1:ifull)/temp(1:ifull) ! rvar = qd*rd+qv*rv
        fc(1:ifull)=(1.-cfrac(1:ifull,k))+cfrac(1:ifull,k)*(lv*rvar(1:ifull)/(cp*rv*temp(1:ifull)))
        dc(1:ifull)=(1.+0.61*qvg(1:ifull,k))*lv*qvg(1:ifull,k)/(rd*tempv(1:ifull))
        mc(1:ifull)=(1.+dc(1:ifull))/(1.+lv*qlg(1:ifull,k)/(cp*temp(1:ifull))+dc(1:ifull)*fc(1:ifull))
        bvf(1:ifull)=grav*mc(1:ifull)*(thetalhl(1:ifull,k)-thetalhl(1:ifull,k-1))/(thetal(1:ifull,k)*dz_fl(1:ifull,k))           &
           +grav*(mc(1:ifull)*fc(1:ifull)*1.61-1.)*(temp(1:ifull)/tempv(1:ifull))*(qthl(1:ifull,k)-qthl(1:ifull,k-1))/dz_fl(1:ifull,k)
        ppb(1:ifull,k)=-km(1:ifull,k)*bvf(1:ifull)
      end do
    
    case(2) ! dry convection
      call updatekmo(thetavhl,thetav,fzzh)
      do k=2,kl-1
        tcc(1:ifull)=-grav*km(1:ifull,k)*(thetavhl(1:ifull,k)-thetavhl(1:ifull,k-1))/(thetav(1:ifull,k)*dz_fl(1:ifull,k))
        ppb(1:ifull,k)=tcc(1:ifull)
      end do
    
   case default
     write(6,*) "ERROR: Unknown buoymeth option ",buoymeth
     stop
  end select

  ! Calculate transport source term on full levels
  do k=2,kl-1
    ppt(1:ifull,k)= kmo(1:ifull,k)*idzp(1:ifull,k)*(tkel(1:ifull,k+1)-tkel(1:ifull,k))/dz_hl(1:ifull,k)                &
             -kmo(1:ifull,k-1)*idzm(1:ifull,k)*(tkel(1:ifull,k)-tkel(1:ifull,k-1))/dz_hl(1:ifull,k-1)
  end do
  
  qq(1:ifull,2:kl-1)=-ddts*idzm(1:ifull,2:kl-1)/dz_hl(1:ifull,1:kl-2)
  rr(1:ifull,2:kl-1)=-ddts*idzp(1:ifull,2:kl-1)/dz_hl(1:ifull,2:kl-1)
  
  ! eps vertical mixing
  aa(1:ifull,2:kl-1)=ce0*kmo(1:ifull,1:kl-2)*qq(1:ifull,2:kl-1)
  cc(1:ifull,2:kl-1)=ce0*kmo(1:ifull,2:kl-1)*rr(1:ifull,2:kl-1)
  ! follow PH to make scheme more numerically stable
  bb(1:ifull,2:kl-1)=1.-aa(1:ifull,2:kl-1)-cc(1:ifull,2:kl-1)+ddts*ce2*epsl(1:ifull,2:kl-1)/tkel(1:ifull,2:kl-1)
  dd(1:ifull,2:kl-1)=epsl(1:ifull,2:kl-1)+ddts*epsl(1:ifull,2:kl-1)/tkel(1:ifull,2:kl-1)                    &
              *ce1*(pps(1:ifull,2:kl-1)+max(ppb(1:ifull,2:kl-1),0.)+max(ppt(1:ifull,2:kl-1),0.))
  dd(1:ifull,2)     =dd(1:ifull,2)   -aa(1:ifull,2)*epsl(1:ifull,1)
  dd(1:ifull,kl-1)  =dd(1:ifull,kl-1)-cc(1:ifull,kl-1)*mineps
  call thomas(epsl(1:ifull,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1))

  ! TKE vertical mixing
  aa(1:ifull,2:kl-1)=kmo(1:ifull,1:kl-2)*qq(1:ifull,2:kl-1)
  cc(1:ifull,2:kl-1)=kmo(1:ifull,2:kl-1)*rr(1:ifull,2:kl-1)
  bb(1:ifull,2:kl-1)=1.-aa(1:ifull,2:kl-1)-cc(1:ifull,2:kl-1)
  dd(1:ifull,2:kl-1)=tkel(1:ifull,2:kl-1)+ddts*(pps(1:ifull,2:kl-1)+ppb(1:ifull,2:kl-1)-epsl(1:ifull,2:kl-1))
  dd(1:ifull,2)     =dd(1:ifull,2)   -aa(1:ifull,2)*tkel(1:ifull,1)
  dd(1:ifull,kl-1)  =dd(1:ifull,kl-1)-cc(1:ifull,kl-1)*mintke
  call thomas(tkel(1:ifull,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1))

  ! limit decay of TKE and EPS with coupling to MF term
  if ( tkemeth==1 ) then
    do k = 2,kl-1
      tbb(1:ifull) = max(1.-0.05*dz_hl(1:ifull,k-1)/250.,0.)
      where ( wstar(1:ifull)>0.5 .and. zz(1:ifull,k)>0.5*zi(1:ifull) .and. zz(1:ifull,k)<0.95*zi(1:ifull)   )
        tkel(1:ifull,k) = max( tkel(1:ifull,k), tbb(1:ifull)*tkel(1:ifull,k-1) )
        epsl(1:ifull,k) = max( epsl(1:ifull,k), tbb(1:ifull)*epsl(1:ifull,k-1) )
      end where
    end do
  end if

  do k=2,kl-1
    tkel(1:ifull,k)=max(tkel(1:ifull,k),mintke)
    tff(1:ifull)=cm34*tkel(1:ifull,k)*sqrt(tkel(1:ifull,k))
    epsl(1:ifull,k)=min(epsl(1:ifull,k),tff(1:ifull)/minl)
    epsl(1:ifull,k)=max(epsl(1:ifull,k),tff(1:ifull)/maxl,mineps)
  end do
    
  km(1:ifull,:)=cm0*tkel(1:ifull,:)*tkel(1:ifull,:)/epsl(1:ifull,:)
  call updatekmo(kmo,km,fzzh) ! interpolate diffusion coeffs to half levels
  
  ! update scalars
  qq(1:ifull,2:kl)  =-ddts*kmo(1:ifull,1:kl-1)*idzm(1:ifull,2:kl)/dz_hl(1:ifull,1:kl-1)
  rr(1:ifull,1:kl-1)=-ddts*kmo(1:ifull,1:kl-1)*idzp(1:ifull,1:kl-1)/dz_hl(1:ifull,1:kl-1)

  ! updating diffusion and non-local terms for qtot and thetal
  ! Note that vertical interpolation is linear so that qtot can be
  ! decomposed into qv, ql, qf and qr.

  cc(1:ifull,1)=rr(1:ifull,1)-ddts*mflx(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1)
  bb(1:ifull,1)=1.-rr(1:ifull,1)-ddts*mflx(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)
  aa(1:ifull,2:kl-1)=qq(1:ifull,2:kl-1)+ddts*mflx(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)
  cc(1:ifull,2:kl-1)=rr(1:ifull,2:kl-1)-ddts*mflx(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1)
  bb(1:ifull,2:kl-1)=1.-qq(1:ifull,2:kl-1)-rr(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)                   &
                                               -mflx(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1))
  aa(1:ifull,kl)=qq(1:ifull,kl)+ddts*mflx(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)
  bb(1:ifull,kl)=1.-qq(1:ifull,kl)+ddts*mflx(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl)

  
  avearray(1:ifull)=0.5*(maxval(thetal(1:ifull,:),dim=2)+minval(thetal(1:ifull,:),dim=2))
  do k=1,kl
    thetal(1:ifull,k)=thetal(1:ifull,k)-avearray(1:ifull)
    tlup(1:ifull,k)=tlup(1:ifull,k)-avearray(1:ifull)
  end do
  dd(1:ifull,1)=thetal(1:ifull,1)-ddts*(mflx(1:ifull,1)*tlup(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                                   &
                                       +mflx(1:ifull,2)*tlup(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1))                                       &
                                 +ddts*rhos(1:ifull)*wt0(1:ifull)/(rhoa(1:ifull,1)*dz_fl(1:ifull,1))
  dd(1:ifull,2:kl-1)=thetal(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,1:kl-2)*tlup(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)     &
                                                 +mflx(1:ifull,2:kl-1)*tlup(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)          &
                                                 -mflx(1:ifull,2:kl-1)*tlup(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1)     &
                                                 -mflx(1:ifull,3:kl)*tlup(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
  dd(1:ifull,kl)=thetal(1:ifull,kl)+ddts*(mflx(1:ifull,kl-1)*tlup(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                       &
                                         +mflx(1:ifull,kl)*tlup(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
  call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  do k=1,kl
    thetal(1:ifull,k)=thetal(1:ifull,k)+avearray(1:ifull)
    tlup(1:ifull,k)=tlup(1:ifull,k)+avearray(1:ifull)
  end do
#ifdef scm  
  wthlflux(1:ifull,1)=wt0(1:ifull)
  wthlflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(thetal(1:ifull,2:kl)-thetal(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                &
                         +mflx(1:ifull,1:kl-1)*(tlup(1:ifull,1:kl-1)-thetal(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))           &
                         +mflx(1:ifull,2:kl)*(tlup(1:ifull,2:kl)-thetal(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif
#ifdef offline
  wthl(1:ifull,1:kl-1)=-kmo(1:ifull,1:kl-1)*(thetal(1:ifull,2:kl)-thetal(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                    &
                       +mflx(1:ifull,1:kl-1)*(tlup(1:ifull,1:kl-1)-thetal(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))               &
                       +mflx(1:ifull,2:kl)*(tlup(1:ifull,2:kl)-thetal(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif


  avearray(1:ifull)=0.5*(maxval(qvg(1:ifull,:),dim=2)+minval(qvg(1:ifull,:),dim=2))
  do k=1,kl
    qvg(1:ifull,k)=qvg(1:ifull,k)-avearray(1:ifull)
    qvup(1:ifull,k)=qvup(1:ifull,k)-avearray(1:ifull)
  end do
  dd(1:ifull,1)=qvg(1:ifull,1)-ddts*(mflx(1:ifull,1)*qvup(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                                      &
                                    +mflx(1:ifull,2)*qvup(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1))                                          &
                              +ddts*rhos(1:ifull)*wq0(1:ifull)/(rhoa(1:ifull,1)*dz_fl(1:ifull,1))
  dd(1:ifull,2:kl-1)=qvg(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,1:kl-2)*qvup(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)        &
                                              +mflx(1:ifull,2:kl-1)*qvup(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)             &
                                              -mflx(1:ifull,2:kl-1)*qvup(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1)        &
                                              -mflx(1:ifull,3:kl)*qvup(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
  dd(1:ifull,kl)=qvg(1:ifull,kl)+ddts*(mflx(1:ifull,kl-1)*qvup(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                          &
                                      +mflx(1:ifull,kl)*qvup(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
  call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  do k=1,kl
    qvg(1:ifull,k)=qvg(1:ifull,k)+avearray(1:ifull)
    qvup(1:ifull,k)=qvup(1:ifull,k)+avearray(1:ifull)
  end do  
#ifdef scm
  wqvflux(1:ifull,1)=wq0(1:ifull)
  wqvflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(qvg(1:ifull,2:kl)-qvg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                        &
                        +mflx(1:ifull,1:kl-1)*(qvup(1:ifull,1:kl-1)-qvg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                &
                        +mflx(1:ifull,2:kl)*(qvup(1:ifull,2:kl)-qvg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif
#ifdef offline
  wqv(1:ifull,1:kl-1)=-kmo(1:ifull,1:kl-1)*(qvg(1:ifull,2:kl)-qvg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                           &
                      +mflx(1:ifull,1:kl-1)*(qvup(1:ifull,1:kl-1)-qvg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                   &
                      +mflx(1:ifull,2:kl)*(qvup(1:ifull,2:kl)-qvg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif


  dd(1:ifull,1)=qlg(1:ifull,1)-ddts*(mflx(1:ifull,1)*qlup(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                                      &
                                    +mflx(1:ifull,2)*qlup(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1))
  dd(1:ifull,2:kl-1)=qlg(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,1:kl-2)*qlup(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)        &
                                              +mflx(1:ifull,2:kl-1)*qlup(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)             &
                                              -mflx(1:ifull,2:kl-1)*qlup(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1)        &
                                              -mflx(1:ifull,3:kl)*qlup(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
  dd(1:ifull,kl)=qlg(1:ifull,kl)+ddts*(mflx(1:ifull,kl-1)*qlup(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                          &
                                      +mflx(1:ifull,kl)*qlup(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef scm
  wqlflux(1:ifull,1)=0.
  wqlflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(qlg(1:ifull,2:kl)-qlg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                        &
                        +mflx(1:ifull,1:kl-1)*(qlup(1:ifull,1:kl-1)-qlg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                &
                        +mflx(1:ifull,2:kl)*(qlup(1:ifull,2:kl)-qlg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif
#ifdef offline
  wql(1:ifull,1:kl-1)=-kmo(1:ifull,1:kl-1)*(qlg(1:ifull,2:kl)-qlg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                           &
                      +mflx(1:ifull,1:kl-1)*(qlup(1:ifull,1:kl-1)-qlg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                   &
                      +mflx(1:ifull,2:kl)*(qlup(1:ifull,2:kl)-qlg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif


  dd(1:ifull,1)=qfg(1:ifull,1)-ddts*(mflx(1:ifull,1)*qfup(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                                      &
                                    +mflx(1:ifull,2)*qfup(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1))
  dd(1:ifull,2:kl-1)=qfg(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,1:kl-2)*qfup(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)        &
                                              +mflx(1:ifull,2:kl-1)*qfup(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)             &
                                              -mflx(1:ifull,2:kl-1)*qfup(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1)        &
                                              -mflx(1:ifull,3:kl)*qfup(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
  dd(1:ifull,kl)=qfg(1:ifull,kl)+ddts*(mflx(1:ifull,kl-1)*qfup(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                          &
                                      +mflx(1:ifull,kl)*qfup(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef scm
  wqfflux(1:ifull,1)=0.
  wqfflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(qfg(1:ifull,2:kl)-qfg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                        &
                        +mflx(1:ifull,1:kl-1)*(qfup(1:ifull,1:kl-1)-qfg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                &
                        +mflx(1:ifull,2:kl)*(qfup(1:ifull,2:kl)-qfg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif
#ifdef offline
  wqf(1:ifull,1:kl-1)=-kmo(1:ifull,1:kl-1)*(qfg(1:ifull,2:kl)-qfg(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)                           &
                      +mflx(1:ifull,1:kl-1)*(qfup(1:ifull,1:kl-1)-qfg(1:ifull,1:kl-1))*(1.-fzzh(1:ifull,1:kl-1))                   &
                      +mflx(1:ifull,2:kl)*(qfup(1:ifull,2:kl)-qfg(1:ifull,2:kl))*fzzh(1:ifull,1:kl-1)
#endif


  ! update cloud fraction terms
  dd(1:ifull,1)=cfrac(1:ifull,1)-ddts*(mflx(1:ifull,1)*cfup(1:ifull,1)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                                    &
                                      +mflx(1:ifull,2)*cfup(1:ifull,2)*fzzh(1:ifull,1)*idzp(1:ifull,1))
  dd(1:ifull,2:kl-1)=cfrac(1:ifull,2:kl-1)+ddts*(mflx(1:ifull,1:kl-2)*cfup(1:ifull,1:kl-2)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1)      &
                                                +mflx(1:ifull,2:kl-1)*cfup(1:ifull,2:kl-1)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)           &
                                                -mflx(1:ifull,2:kl-1)*cfup(1:ifull,2:kl-1)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1)      &
                                                -mflx(1:ifull,3:kl)*cfup(1:ifull,3:kl)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
  dd(1:ifull,kl)=cfrac(1:ifull,kl)+ddts*(mflx(1:ifull,kl-1)*cfup(1:ifull,kl-1)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                        &
                                        +mflx(1:ifull,kl)*cfup(1:ifull,kl)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
  call thomas(cfrac,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfrac(1:ifull,:)=min(max(cfrac(1:ifull,:),0.),1.)


  ! Aerosols
  do j=1,naero
    dd(1:ifull,1)=aero(1:ifull,1,j)-ddts*(mflx(1:ifull,1)*arup(1:ifull,1,j)*(1.-fzzh(1:ifull,1))*idzp(1:ifull,1)                               &
                                         +mflx(1:ifull,2)*arup(1:ifull,2,j)*fzzh(1:ifull,1)*idzp(1:ifull,1))
    dd(1:ifull,2:kl-1)=aero(1:ifull,2:kl-1,j)+ddts*(mflx(1:ifull,1:kl-2)*arup(1:ifull,1:kl-2,j)*(1.-fzzh(1:ifull,1:kl-2))*idzm(1:ifull,2:kl-1) &
                                                   +mflx(1:ifull,2:kl-1)*arup(1:ifull,2:kl-1,j)*fzzh(1:ifull,1:kl-2)*idzm(1:ifull,2:kl-1)      &
                                                   -mflx(1:ifull,2:kl-1)*arup(1:ifull,2:kl-1,j)*(1.-fzzh(1:ifull,2:kl-1))*idzp(1:ifull,2:kl-1) &
                                                   -mflx(1:ifull,3:kl)*arup(1:ifull,3:kl,j)*fzzh(1:ifull,2:kl-1)*idzp(1:ifull,2:kl-1))
    dd(1:ifull,kl)=aero(1:ifull,kl,j)+ddts*(mflx(1:ifull,kl-1)*arup(1:ifull,kl-1,j)*(1.-fzzh(1:ifull,kl-1))*idzm(1:ifull,kl)                   &
                                           +mflx(1:ifull,kl)*arup(1:ifull,kl,j)*fzzh(1:ifull,kl-1)*idzm(1:ifull,kl))
    call thomas(aero(:,:,j),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
    aero(1:ifull,:,j) = max( aero(1:ifull,:,j), 0. )    
  end do


  ! Winds
  aa(1:ifull,2:kl)  =qq(1:ifull,2:kl)
  cc(1:ifull,1:kl-1)=rr(1:ifull,1:kl-1)
  bb(1:ifull,1)=1.-cc(1:ifull,1)+ddts*rhos(1:ifull)*cdrag(1:ifull)*umag(1:ifull)/(rhoa(1:ifull,1)*dz_fl(1:ifull,1)) ! implicit
  bb(1:ifull,2:kl-1)=1.-aa(1:ifull,2:kl-1)-cc(1:ifull,2:kl-1)
  bb(1:ifull,kl)=1.-aa(1:ifull,kl)
  dd(1:ifull,1:kl)=uo(1:ifull,1:kl)
  ! bb(1:ifull,1)=1.-cc(1:ifull,1)                                           ! explicit
  ! dd(1:ifull,1:kl)=uo(1:ifull,1:kl)-ddts*taux(1:ifull)/(rhoa(1:ifull,1)*dz_fl(:i1)) ! explicit
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  dd(1:ifull,1:kl)=vo(1:ifull,1:kl)
  ! dd(1:ifull,1:kl)=vo(1:ifull,1:kl)-ddts*tauy(1:ifull)/(rhoa(1:ifull,1)*dz_fl(1:ifull,1)) ! explicit
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  
  ! umag=sqrt(max(uo(1:ifull,1)*uo(1:ifull,1)+vo(1:ifull,1)*vo(1:ifull,1),1.e-4)) ! explicit
  ! call dyerhicks(cdrag,wtv0,zom,umag,thetav(:,1),zz(:,1))                       ! explicit
  ! taux=rhos*cdrag*umag*uo(1:ifull,1)                                            ! explicit
  ! tauy=rhos*cdrag*umag*vo(1:ifull,1)                                            ! explicit
  ! ustar=sqrt(sqrt(taux*taux+tauy*tauy)/rhos)                                    ! explicit

  ! account for phase transitions
  do k=1,kl
    tgg(1:ifull)=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k),qgmin) ! qtot before phase transition
    qvg(1:ifull,k)=max(qvg(1:ifull,k),0.)    
    qlg(1:ifull,k)=max(qlg(1:ifull,k),0.)
    qfg(1:ifull,k)=max(qfg(1:ifull,k),0.)
    tff(1:ifull)=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k),qgmin)
    tgg(1:ifull)=tgg(1:ifull)/tff(1:ifull) ! scale factor for conservation
    qvg(1:ifull,k)=qvg(1:ifull,k)*tgg(1:ifull)
    qlg(1:ifull,k)=qlg(1:ifull,k)*tgg(1:ifull)
    qfg(1:ifull,k)=qfg(1:ifull,k)*tgg(1:ifull)
    ! update theta for output or next time step
    theta(1:ifull,k)=thetal(1:ifull,k)+sigkap(k)*(lv*qlg(1:ifull,k)+ls*qfg(1:ifull,k))/cp
    where (qlg(1:ifull,k)+qfg(1:ifull,k)>1.E-12)
      cfrac(1:ifull,k)=max(cfrac(1:ifull,k),1.E-8)
    end where
  end do
  
#ifdef scm
  uwflux(1:ifull,1)=0.
  uwflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(uo(1:ifull,2:kl)-uo(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)
  vwflux(1:ifull,1)=0.
  vwflux(1:ifull,2:kl)=-kmo(1:ifull,1:kl-1)*(vo(1:ifull,2:kl)-vo(1:ifull,1:kl-1))/dz_hl(1:ifull,1:kl-1)
  do k=1,kl
    wthflux(1:ifull,k) = wthlflux(1:ifull,k) - (sigkap(k-1)*(1.-fzzh(1:ifull,k)+sigkap(k)*fzzh(1:ifull,k)) &
                                             *(lv*wqlflux(1:ifull,k)+ls*wqfflux(1:ifull,k)))
  end do
#endif
  
end do


return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

subroutine thomas(outdat,aai,bbi,cci,ddi)

implicit none

real, dimension(:,2:), intent(in) :: aai
real, dimension(:,:), intent(in) :: bbi,ddi
real, dimension(:,:), intent(in) :: cci
real, dimension(:,:), intent(out) :: outdat
real, dimension(ifull,size(outdat,2)) :: cc,dd
real, dimension(ifull) :: n
integer k,klin

klin=size(outdat,2)
cc(1:ifull,1)=cci(1:ifull,1)/bbi(1:ifull,1)
dd(1:ifull,1)=ddi(1:ifull,1)/bbi(1:ifull,1)

do k=2,klin-1
  n(1:ifull)=bbi(1:ifull,k)-cc(1:ifull,k-1)*aai(1:ifull,k)
  cc(1:ifull,k)=cci(1:ifull,k)/n(1:ifull)
  dd(1:ifull,k)=(ddi(1:ifull,k)-dd(1:ifull,k-1)*aai(1:ifull,k))/n(1:ifull)
end do
n(1:ifull)=bbi(1:ifull,klin)-cc(1:ifull,klin-1)*aai(1:ifull,klin)
dd(1:ifull,klin)=(ddi(1:ifull,klin)-dd(1:ifull,klin-1)*aai(1:ifull,klin))/n(1:ifull)
outdat(1:ifull,klin)=dd(1:ifull,klin)
do k=klin-1,1,-1
  outdat(1:ifull,k)=dd(1:ifull,k)-cc(1:ifull,k)*outdat(1:ifull,k+1)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(size(temp)), intent(in) :: ps
real, dimension(size(temp)), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix
logical, save :: first=.true.


tdiff=min(max( temp(:)-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat(:)=0.622*esatf/max(ps(:)-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

subroutine updatekmo(kmo,km,fzhl)

implicit none

real, dimension(:,:), intent(in) :: km
real, dimension(ifull,kl-1), intent(in) :: fzhl
real, dimension(ifull,kl), intent(out) :: kmo

kmo(1:ifull,1:kl-1)=km(1:ifull,1:kl-1)+fzhl(1:ifull,1:kl-1)*(km(1:ifull,2:kl)-km(1:ifull,1:kl-1))
! These terms are never used
kmo(1:ifull,kl)=0.

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

real function entfn(zht,zi)

implicit none

real, intent(in) :: zht,zi

!entfn=0.002                                               ! Angevine (2005)
!entfn=2./max(100.,zi)                                     ! Angevine et al (2010)
!entfn=1./zht                                              ! Siebesma et al (2003)
!entfn=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))       ! Soares et al (2004)
entfn = ent0/max( zht, 1. ) + ent1/max( zi-zht, ezmin )

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral detrainment

real function dtrfn(zht,zi,rat)

implicit none

real, intent(in) :: zht,zi,rat

!dtrfn=ent+0.05/max(zi-zht,zmin)   ! Angevine et al (2010)
dtrfn=rat/max(zi-zht,1.)+ent1/max(zi-zht,ezmin)

! results in analytic solution
!mflx(k)=A*(zht**ent0)*((zi-zht)**rat)


return
end function dtrfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate drag coeff

subroutine dyerhicks(cd,wtv0,zom,umag,thetav,zmin)

implicit none

integer ic
integer, parameter :: icmax = 10
real, dimension(ifull), intent(in) :: umag,thetav,zom,wtv0,zmin
real, dimension(ifull), intent(out) :: cd
real, dimension(ifull) :: ustar,thetavstar,ilzom
real, dimension(ifull) :: z_on_l,z0_on_l
real, dimension(ifull) :: pm0,pm1,integralm
real, dimension(ifull) :: dumzmin

dumzmin    = max(zmin,zom+0.2)
ilzom      = log(dumzmin/zom)
ustar      = vkar*umag/ilzom ! first guess

select case(stabmeth)
  case(0)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*dumzmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/dumzmin
      where ( z_on_l<0. )
        pm0     = (1.-16.*z0_on_l)**(-0.25)
        pm1     = (1.-16.*z_on_l )**(-0.25)
        integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                   +2.*(atan(1./pm1)-atan(1./pm0))
      elsewhere
        !--------------Beljaars and Holtslag (1991) momentum & heat            
        pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
        pm1 = -(a_1*z_on_l +b_1*(z_on_l -(c_1/d_1))*exp(-d_1*z_on_l )+b_1*c_1/d_1)
        integralm = ilzom-(pm1-pm0)
      end where
      ustar = vkar*umag/integralm
    end do
    
  case(1)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*dumzmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/dumzmin
      where ( z_on_l<0. )
        pm0     = (1.-16.*z0_on_l)**(-0.25)
        pm1     = (1.-16.*z_on_l )**(-0.25)
        integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                   +2.*(atan(1./pm1)-atan(1./pm0))
      elsewhere
        !--------------Beljaars and Holtslag (1991) momentum & heat            
        pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
        pm1 = -(a_1*z_on_l +b_1*(z_on_l -(c_1/d_1))*exp(-d_1*z_on_l )+b_1*c_1/d_1)
        integralm = ilzom-(pm1-pm0)
      end where
      where ( z_on_l<=0.4 )
        ustar = vkar*umag/integralm
      elsewhere ! Luhar
        ustar = vkar*umag/(aa1*(( z_on_l**bb1)*(1.+ cc1*z_on_l**(1.-bb1)) &
                             -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1))))
      end where
    end do
    
  case default
    write(6,*) "ERROR: Invalid option for stabmeth in tkeeps ",stabmeth
    stop
    
end select

cd = (vkar/integralm)**2

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE-eps

subroutine tkeend(diag)

implicit none

integer, intent(in) :: diag

if (diag>0) write(6,*) "Terminate TKE-eps scheme"

deallocate(tke,eps)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps

