! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! clouds and saturated air following Marquet and Geleyn 2012.

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
public m0,b1,b2,qcmf,ent_min,mfbeta
public buoymeth,maxdts,mintke,mineps,minl,maxl,stabmeth
public tke_umin,tkemeth,zimax
#ifdef offline
public wthl,wqv,wql,wqf
public mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
public ents,dtrs
#endif

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
#ifdef offline
real, dimension(:,:), allocatable, save :: wthl,wqv,wql,wqf
real, dimension(:,:), allocatable, save :: mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
real, dimension(:,:), allocatable, save :: u,v,ents,dtrs
#endif

! model ED constants
real, save :: cm0      = 0.09     ! Hurley (2007) 0.09, Duynkerke (1988) 0.03, Duynkerke (1987) 0.09
real, save :: ce0      = 0.69     ! Hurley (2007) 0.69, Duynkerke (1988) 0.42, Duynkerke (1987) 0.77
real, save :: ce1      = 1.46
real, save :: ce2      = 1.83
real, save :: ce3      = 0.45     ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, save :: cq       = 2.5      ! Adjustment to ED in absence of MF
real, save :: mintke   = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps   = 1.E-11   ! min value for eps (1.0e-6 in TAPM)
real, save :: minl     = 1.       ! min value for L   (5. in TAPM)
real, save :: maxl     = 1000.    ! max value for L   (500. in TAPM)
real, save :: tke_umin = 0.1      ! minimum wind speed (m/s) for drag calculation
real, save :: zimax    = 10000.   ! maximum PBL height feeding back into model calculations
! model MF constants
real, save :: be       = 0.1      ! Surface boundary condition (Hurley (2007) 1., Soares et al (2004) 0.3)
real, save :: ent0     = 0.25     ! Entrainment constant (Controls height of boundary layer) (Hurley (2007) 0.5)
real, save :: ent1     = 0.25
real, save :: ent_min  = 0.       ! Minimum entrainment
real, save :: ezmin    = 100.     ! Limits entrainment at plume top
real, save :: entc0    = 2.e-3    ! Saturated entrainment constant for mass flux
real, save :: dtrc0    = 3.e-3    ! Saturated detrainment constant for mass flux
real, save :: m0       = 0.1      ! Mass flux area constant (Hurley (2007) 0.1)
real, save :: b1       = 2.       ! Updraft entrainment coeff (Soares et al (2004) 1., Siebesma et al (2003) 2.)
real, save :: b2       = 1./3.    ! Updraft buoyancy coeff (Soares et al (2004) 2., Siebesma et al (2003) 1./3.)
real, save :: qcmf     = 1.e-4    ! Critical mixing ratio of liquid water before autoconversion
real, save :: mfbeta   = 0.15     ! Horizontal scale factor
! generic constants
integer, save :: buoymeth = 1     ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0     ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: tkemeth  = 1     ! Method for TKE calculation (0=D&K84, 1=Hurley)
real, save :: maxdts      = 120.  ! max timestep for split

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
real, parameter :: tfrz  = 273.16     ! (K)
real, parameter :: tice  = 233.1      ! (K)

! MOST constants
real, parameter :: a_1   = 1.
real, parameter :: b_1   = 2./3.
real, parameter :: c_1   = 5.
real, parameter :: d_1   = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

real, dimension(0:220), save :: tablei
real, dimension(-40:2), save :: esdiff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifullin,iextrain,klin,diag)

implicit none

integer, intent(in) :: ifullin,iextrain,klin,diag
real cm34

if (diag>0) write(6,*) "Initialise TKE-eps scheme"

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

tablei(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                               !-146C
tablei(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                            !-141C
tablei(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                         !-136C
tablei(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /) !-131C
tablei(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /) !-126C
tablei(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)  !-121C
tablei(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)      !-116C
tablei(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)        !-111C
tablei(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)           !-106C
tablei(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)            !-101C
tablei(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)            !-95C
tablei(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /)!-87C
tablei(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /)!-78C
tablei(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)  !-69C
tablei(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)   !-60C
tablei(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)   !-51C
tablei(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)           !-43C
tablei(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)   !-34C
tablei(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /)!-24C 
tablei(127:134)=(/ 77.09, 85.02, 93.70, 103.06, 113.40, 124.68, 136.98, 150.39 /)     !-16C
tablei(135:142)=(/ 164.99, 180.88, 198.16, 216.94, 237.34, 259.47, 283.49, 309.51 /)  !-8C
tablei(143:150)=(/ 337.71, 368.23, 401.25, 436.96, 475.54, 517.21, 562.19, 610.70 /)  !0C
tablei(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)  !8C
tablei(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)  !16C
tablei(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)  !24C
tablei(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)  !32C
tablei(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)  !40C
tablei(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)        !47C
tablei(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)   !54C
tablei(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)   !61C
tablei(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)   !68C
tablei(219:220)=(/ 29845.0, 31169.0 /)                                                !70C 

esdiff(-40:-31)=(/ 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27 /)
esdiff(-30:-21)=(/ 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61 /)
esdiff(-20:-11)=(/ 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27 /)
esdiff(-10:-1)= (/ 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65 /)
esdiff(0:2)=    (/ 0.08, 0., 0. /)

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,stratcloud,uo,vo,zi,fg,eg,cduv,ps,zz,zzh,sig,rhos, &
                  ustar_ave,dt,qgmin,mode,diag,tke,eps,shear,dx,                           &
#ifdef scm
                  wthflux,wqvflux,uwflux,vwflux,mfout,buoyproduction,                      &
                  shearproduction,totaltransport,                                          &
#endif
#ifdef offline
                  mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up,ents_up,dtrs_up,wthl,wqv,wql,wqf,  &
#endif
                  imax)

implicit none

integer, intent(in) :: imax
integer, intent(in) :: diag,mode
integer k, j, imax_p
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(:,:), intent(inout) :: theta,stratcloud,uo,vo
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg
real, dimension(imax,kl), intent(out) :: kmo
real, dimension(imax,kl), intent(in) :: zz,zzh
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(inout) :: eps
real, dimension(imax), intent(inout) :: zi
real, dimension(imax), intent(in) :: fg,eg,cduv,ps,rhos,dx
real, dimension(imax), intent(out) :: ustar_ave
real, dimension(kl), intent(in) :: sig
real, dimension(imax,kl) :: km,thetav,thetal,qsat
real, dimension(imax,kl) :: qsatc,qgnc,ff
real, dimension(imax,kl) :: thetalhl,thetavhl
real, dimension(imax,kl) :: quhl,qshl,qlhl,qfhl
real, dimension(imax,kl) :: bb,cc,dd,rr
real, dimension(imax,kl) :: rhoa,rhoahl
real, dimension(imax,kl) :: qtot,qthl
real, dimension(imax,kl) :: tlup,qvup,qlup,qfup
real, dimension(imax,kl) :: cfup,mflx
real, dimension(imax,kl) :: pps,ppt,ppb
real, dimension(imax,2:kl) :: idzm
real, dimension(imax,1:kl-1) :: idzp
real, dimension(imax,2:kl) :: aa,qq
real, dimension(imax,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(imax,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(imax,kl-1) :: fzzh
real, dimension(imax) :: wt0,wq0,wtv0
real, dimension(imax) :: wstar,z_on_l,phim
real, dimension(imax) :: tff,tgg,tempc,thetac,pres,temp
real, dimension(imax) :: cdrag,umag,ustar
real, dimension(imax) :: tempv,rvar,bvf,dc,mc,fc
real, dimension(imax) :: tbb,tcc,tqq
real, dimension(imax) :: zi_save, zturb, cgmap
real, dimension(imax) :: templ, dqsdt, al, fice, qc, qt, lx
real, dimension(kl) :: sigkap
real cm12, cm34
real ddts
logical, dimension(imax,kl) :: lta
logical, dimension(imax) :: lmask

#ifdef offline
real, dimension(imax,kl), intent(inout) :: mf
real, dimension(imax,kl), intent(inout) :: w_up
real, dimension(imax,kl), intent(inout) :: tl_up
real, dimension(imax,kl), intent(inout) :: qv_up
real, dimension(imax,kl), intent(inout) :: ql_up
real, dimension(imax,kl), intent(inout) :: qf_up
real, dimension(imax,kl), intent(inout) :: cf_up
real, dimension(imax,kl), intent(inout) :: ents
real, dimension(imax,kl), intent(inout) :: dtrs
real, dimension(imax,kl), intent(inout) :: wthl
real, dimension(imax,kl), intent(inout) :: wqv
real, dimension(imax,kl), intent(inout) :: wql
real, dimension(imax,kl), intent(inout) :: wqf
#endif

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl), intent(out) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(out) :: totaltransport
real, dimension(imax,kl-1), intent(out) :: mfout
real, dimension(imax,kl) :: wthlflux, wqlflux
real, dimension(imax,kl) :: wqfflux
#endif

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))

if ( diag>0 ) write(6,*) "Update PBL mixing with TKE-eps + MF turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

sigkap(1:kl) = sig(1:kl)**(-rd/cp)

do concurrent (k = 1:kl)
  ! Impose limits on tke and eps after being advected by the host model
  tke(:,k) = max(tke(:,k), mintke)
  tff(:)   = cm34*tke(:,k)*sqrt(tke(:,k))
  eps(:,k) = min(eps(:,k), tff/minl)
  eps(:,k) = max(eps(:,k), tff/maxl, mineps)
  
  ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
  pres(:) = ps(:)*sig(k) ! pressure
  ! Density must be updated when dz is updated so that rho*dz is conserved
  thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
  rhoa(:,k) = sigkap(k)*pres(:)/(rd*thetav(:,k))

  ! Transform to thetal as it is the conserved variable
  thetal(:,k) = theta(:,k) - sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
end do

! Calculate first approximation to diffusion coeffs
km = cm0*tke*tke/eps
  
! Calculate surface fluxes
wt0 = fg/(rhos*cp)  ! theta flux
wq0 = eg/(rhos*lv)  ! qtot flux

! Fraction for interpolation from full levels to half levels
fzzh(:,1:kl-1) = (zzh(:,1:kl-1)-zz(:,1:kl-1))/(zz(:,2:kl)-zz(:,1:kl-1))

! Calculate dz at half levels
dz_hl(:,1:kl-1) = max( zz(:,2:kl) - zz(:,1:kl-1), 1. )

! Calculate dz at full levels
dz_fl(:,1)    = zzh(:,1)
dz_fl(:,2:kl) = zzh(:,2:kl) - zzh(:,1:kl-1)
dz_fl = max( dz_fl, 1. )

! Calculate shear term on full levels
pps(:,1:kl-1) = km(:,1:kl-1)*shear(:,1:kl-1)

! set top boundary condition for TKE-eps source terms
pps(:,kl) = 0.
ppb(:,kl) = 0.
ppt(:,1)  = 0.
ppt(:,kl) = 0.

! interpolate diffusion coeff and air density to half levels
call updatekmo(kmo,   km,  fzzh,imax)
call updatekmo(rhoahl,rhoa,fzzh,imax)
! eddy diffusion terms to account for air density with level thickness
idzm(:,2:kl)   = rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1) = rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

ustar_ave(:) = 0.

! Main loop to prevent time splitting errors
mcount = int(dt/(maxdts+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  ! Set-up thermodynamic variables temp, theta_v, qtot and saturated mixing ratio
  do concurrent (k = 1:kl)
    templ(:) = thetal(:,k)/sigkap(k)
    pres(:) = ps(:)*sig(k)
    ! calculate saturated air mixing ratio
    where ( qfg(:,k)>1.e-12 )
      fice = min( qfg(:,k)/(qfg(:,k)+qlg(:,k)), 1. )
    elsewhere
      fice = 0.
    end where
    call getqsat(qsat(:,k),templ(:),pres(:),fice)
    thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
    qtot(:,k) = qvg(:,k) + qlg(:,k) + qfg(:,k)
  end do
  
  ! Update thetav flux
  wtv0 = wt0 + theta(:,1)*0.61*wq0 ! thetav flux
  
  ! Update momentum flux
  ustar = sqrt(cduv*sqrt(uo(:,1)**2+vo(:,1)**2))  
  
  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  mflx(:,:) = 0.
  tlup(:,:) = thetal(:,:)
  qvup(:,:) = qvg(:,:)
  qlup(:,:) = qlg(:,:)
  qfup(:,:) = qfg(:,:)
  cfup(:,:) = stratcloud(:,:)

#ifdef scm
  mfout(:,:)=0.
#endif
#ifdef offline
  mf=0.
  w_up=0.
  tl_up=thetal(:,:)
  qv_up=qvg(:,:)
  ql_up=qlg(:,:)
  qf_up=qfg(:,:)
  cf_up=0.
  ents=0.
  dtrs=0.
#endif

  wstar = (grav*min(zi,zimax)*max(wtv0,0.)/thetav(:,1))**(1./3.)   

  if ( mode/=1 ) then ! mass flux
      
    zi_save = zi  

    lmask = wtv0(:)>0. ! unstable
    imax_p = count(lmask)
    ! unstable boundary layer
    call plumerise(imax_p,lmask,cm12,                            &
                   zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,       &
                   zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,     &
                   stratcloud,km,ustar,wt0,wq0,wtv0,ps,          &
                   sig,sigkap,tke,eps,imax)

    ! stable boundary layer
    lmask = .not.lmask ! stable
#ifdef scm
    where ( lmask )
      zi = zz(:,1)
    end where  
#else
    call stablepbl(lmask,zi,zzh,dz_hl,thetav,kmo,imax)
#endif
    
    ! Turn off MF term if small grid spacing (beta=0. implies MF is always non-zero)
    ! Based on Boutle et al 2014
    zturb = min( 0.5*(zi_save + zi), zimax )
#ifdef scm
    cgmap(:) = 1.
#else
    cgmap(:) = 1. - tanh(mfbeta*zturb/dx)*max(0.,1.-0.25*dx/zturb)
#endif
    do concurrent (k = 1:kl)
      mflx(:,k) = mflx(:,k)*cgmap(:)
    end do
    
  end if

#ifdef scm  
  mfout(:,1:kl-1) = mflx(:,1:kl-1)*(1.-fzzh(:,1:kl-1)) &
                  + mflx(:,2:kl)*fzzh(:,1:kl-1)
#endif
  
  
  ! calculate tke and eps boundary condition at 1st vertical level
  z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-10))
  z_on_l=min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  call calc_phi(phim,z_on_l)
  tke(:,1) = cm12*ustar*ustar+ce3*wstar*wstar
  eps(:,1) = ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
  tke(:,1) = max( tke(:,1), mintke )
  tff = cm34*tke(:,1)*sqrt(tke(:,1))
  eps(:,1) = min( eps(:,1), tff/minl )
  eps(:,1) = max( eps(:,1), tff/maxl, mineps )


  ! Update TKE and eps terms

  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tke(:,kl) = mintke
  eps(:,kl) = mineps
  
  ! Calculate buoyancy term
  select case(buoymeth)
    case(0) ! Blend staturated and unsaturated terms - saturated method from Durran and Klemp JAS 1982 (see also WRF)
      qsatc=max(qsat,qvg(:,:))                          ! assume qvg is saturated inside cloud
      ff=qfg(:,:)/max(stratcloud(:,:),1.E-8)            ! inside cloud value  assuming max overlap
      dd=qlg(:,:)/max(stratcloud(:,:),1.E-8)            ! inside cloud value assuming max overlap
      do concurrent (k = 1:kl)
        tbb=max(1.-stratcloud(:,k),1.E-8)
        qgnc(:,k)=(qvg(:,k)-(1.-tbb)*qsatc(:,k))/tbb    ! outside cloud value
        qgnc(:,k)=min(max(qgnc(:,k),qgmin),qsatc(:,k))
      end do
      call updatekmo(thetalhl,thetal,fzzh,imax)              ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh,imax)                    ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh,imax)                   ! inside cloud value
      call updatekmo(qlhl,dd,fzzh,imax)                      ! inside cloud value
      call updatekmo(qfhl,ff,fzzh,imax)                      ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl)=stratcloud(:,2:kl)<=1.E-6
      do concurrent (k = 2:kl-1)
        where(lta(:,k).and..not.lta(:,k+1))
          qlhl(:,k)=dd(:,k+1)
          qfhl(:,k)=ff(:,k+1)
        elsewhere (.not.lta(:,k).and.lta(:,k+1))
          qlhl(:,k)=dd(:,k)
          qfhl(:,k)=ff(:,k)
        end where
      end do
      do concurrent (k = 2:kl-1)
        ! saturated
        thetac(:)=thetal(:,k)+sigkap(k)*(lv*dd(:,k)+ls*ff(:,k))/cp              ! inside cloud value
        tempc(:)=thetac(:)/sigkap(k)                                            ! inside cloud value          
        tqq=(1.+lv*qsatc(:,k)/(rd*tempc(:)))/(1.+lv*lv*qsatc(:,k)/(cp*rv*tempc(:)*tempc(:)))
        tbb=-grav*km(:,k)*(tqq*((thetalhl(:,k)-thetalhl(:,k-1)+sigkap(k)/cp*(lv*(qlhl(:,k)-qlhl(:,k-1))  &
            +ls*(qfhl(:,k)-qfhl(:,k-1))))/thetac(:)+lv/cp*(qshl(:,k)-qshl(:,k-1))/tempc(:))              &
            -qshl(:,k)-qlhl(:,k)-qfhl(:,k)+qshl(:,k-1)+qlhl(:,k-1)+qfhl(:,k-1))/dz_fl(:,k)
        ! unsaturated
        tcc=-grav*km(:,k)*(thetalhl(:,k)-thetalhl(:,k-1)+thetal(:,k)*0.61*(quhl(:,k)-quhl(:,k-1)))  &
                         /(thetal(:,k)*dz_fl(:,k))
        ppb(:,k)=(1.-stratcloud(:,k))*tcc+stratcloud(:,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)
      end do
      ! saturated
      thetac(:)=thetal(:,1)+sigkap(1)*(lv*dd(:,1)+ls*ff(:,1))/cp              ! inside cloud value
      tempc(:)=thetac(:)/sigkap(1)                                            ! inside cloud value          
      tqq=(1.+lv*qsatc(:,1)/(rd*tempc(:)))/(1.+lv*lv*qsatc(:,1)/(cp*rv*tempc(:)*tempc(:)))
      tbb=-grav*km(:,1)*(tqq*((thetalhl(:,1)-thetal(:,1)+sigkap(1)/cp*(lv*(qlhl(:,1)-qlg(:,1))         &
          +ls*(qfhl(:,1)-qfg(:,1))))/thetac(:)+lv/cp*(qshl(:,1)-qsatc(:,1))/tempc(:))                  &
          -qshl(:,1)-qlhl(:,1)-qfhl(:,1)+qsatc(:,1)+qlg(:,1)+qfg(:,1))/(zzh(:,1)-zz(:,1))
      ! unsaturated
      tcc=-grav*km(:,1)*(thetalhl(:,1)-thetal(:,1)+thetal(:,1)*0.61*(quhl(:,1)-qgnc(:,1)))        &
                       /(thetal(:,1)*(zzh(:,1)-zz(:,1)))
      ppb(:,1)=(1.-stratcloud(:,1))*tcc+stratcloud(:,1)*tbb ! cloud fraction weighted (e.g., Smith 1990)

      
    case(1) ! Marquet and Geleyn QJRMS (2012) for partially saturated
      call updatekmo(thetalhl,thetal,fzzh,imax)
      call updatekmo(qthl,qtot,fzzh,imax)
      do concurrent (k = 2:kl-1)
        temp(:)  = theta(:,k)/sigkap(k)
        tempv(:) = thetav(:,k)/sigkap(k)
        rvar=rd*tempv/temp ! rvar = qd*rd+qv*rv
        fc=(1.-stratcloud(:,k))+stratcloud(:,k)*(lv*rvar/(cp*rv*temp))
        dc=(1.+0.61*qvg(:,k))*lv*qvg(:,k)/(rd*tempv)
        mc=(1.+dc)/(1.+lv*qlg(:,k)/(cp*temp)+dc*fc)
        bvf=grav*mc*(thetalhl(:,k)-thetalhl(:,k-1))/(thetal(:,k)*dz_fl(:,k))           &
           +grav*(mc*fc*1.61-1.)*(temp/tempv)*(qthl(:,k)-qthl(:,k-1))/dz_fl(:,k)
        ppb(:,k)=-km(:,k)*bvf
      end do
      temp(:)  = theta(:,1)/sigkap(1)
      tempv(:) = thetav(:,1)/sigkap(1)
      rvar=rd*tempv/temp ! rvar = qd*rd+qv*rv
      fc=(1.-stratcloud(:,1))+stratcloud(:,1)*(lv*rvar/(cp*rv*temp))
      dc=(1.+0.61*qvg(:,1))*lv*qvg(:,1)/(rd*tempv)
      mc=(1.+dc)/(1.+lv*qlg(:,1)/(cp*temp)+dc*fc)
      bvf=grav*mc*(thetalhl(:,1)-thetal(:,1))/(thetal(:,1)*(zzh(:,1)-zz(:,1)))         &
         +grav*(mc*fc*1.61-1.)*(temp/tempv)*(qthl(:,1)-qtot(:,1))/(zzh(:,1)-zz(:,1))
      ppb(:,1) = -km(:,1)*bvf
      
    case(2) ! dry convection from Hurley 2007
      call updatekmo(thetavhl,thetav,fzzh,imax)
      do concurrent (k = 2:kl-1)
        tcc=-grav*km(:,k)*(thetavhl(:,k)-thetavhl(:,k-1))/(thetav(:,k)*dz_fl(:,k))
        ppb(:,k)=tcc
      end do
      tcc=-grav*km(:,1)*(thetavhl(:,1)-thetav(:,1))/(thetav(:,1)*(zzh(:,1)-zz(:,1)))
      ppb(:,1)=tcc 
      
   case default
     write(6,*) "ERROR: Unknown buoymeth option ",buoymeth
     stop
  end select

  ! Calculate transport source term on full levels
  ppt(:,2:kl-1)= kmo(:,2:kl-1)*idzp(:,2:kl-1)*(tke(:,3:kl)-tke(:,2:kl-1))/dz_hl(:,2:kl-1)  &
               -kmo(:,1:kl-2)*idzm(:,2:kl-1)*(tke(:,2:kl-1)-tke(:,1:kl-2))/dz_hl(:,1:kl-2)
  
  ! Pre-calculate eddy diffusivity mixing terms
  ! -ve because gradient is calculated at t+1
  qq(:,2:kl-1)=-ddts*idzm(:,2:kl-1)/dz_hl(:,1:kl-2)
  rr(:,2:kl-1)=-ddts*idzp(:,2:kl-1)/dz_hl(:,2:kl-1)
  
  ! eps vertical mixing
  aa(:,2:kl-1)=ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=ce0*kmo(:,2:kl-1)*rr(:,2:kl-1)
  ! follow Hurley 2007 to make scheme more numerically stable
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)+ddts*ce2*eps(:,2:kl-1)/tke(:,2:kl-1)
  dd(:,2:kl-1)=eps(:,2:kl-1)+ddts*eps(:,2:kl-1)/tke(:,2:kl-1)                    &
              *ce1*(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
  dd(:,2)     =dd(:,2)   -aa(:,2)*eps(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mineps
  call thomas(eps(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax)

  ! TKE vertical mixing
  aa(:,2:kl-1)=kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=kmo(:,2:kl-1)*rr(:,2:kl-1)
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  dd(:,2:kl-1)=tke(:,2:kl-1)+ddts*(pps(:,2:kl-1)+ppb(:,2:kl-1)-eps(:,2:kl-1))
  dd(:,2)     =dd(:,2)   -aa(:,2)*tke(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mintke
  call thomas(tke(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax)

  ! limit decay of TKE and EPS with coupling to mass flux term
  if ( tkemeth==1 ) then
    do concurrent (k = 2:kl-1)
      tbb(:) = max(1.-0.05*dz_hl(:,k-1)/250.,0.)
      where ( wstar(:)>0.5 .and. zz(:,k)>0.5*min(zi,zimax) .and. zz(:,k)<0.95*min(zi,zimax) )
        tke(:,k) = max( tke(:,k), tbb(:)*tke(:,k-1) )
        eps(:,k) = max( eps(:,k), tbb(:)*eps(:,k-1) )
      end where
    end do
  end if
  
  ! apply limits and corrections to tke and eps terms
  do concurrent (k = 2:kl-1)
    tke(:,k) = max(tke(:,k),mintke)
    tff = cm34*tke(:,k)*sqrt(tke(:,k))
    eps(:,k) = min(eps(:,k),tff/minl)
    eps(:,k) = max(eps(:,k),tff/maxl,mineps)
  end do
    
  ! estimate eddy diffusivity mixing coeff
  km = cm0*tke(:,:)*tke(:,:)/eps(:,:)
  call updatekmo(kmo,km,fzzh,imax) ! interpolate diffusion coeffs to half levels
  
  ! update scalars
  
  ! Pre-calculate eddy diffiusivity mixing terms using updated kmo values
  ! -ve because gradient is calculated at t+1
  qq(:,2:kl)  =-ddts*kmo(:,1:kl-1)*idzm(:,2:kl)/dz_hl(:,1:kl-1)
  rr(:,1:kl-1)=-ddts*kmo(:,1:kl-1)*idzp(:,1:kl-1)/dz_hl(:,1:kl-1)

  ! updating diffusion and non-local terms for qtot and thetal
  ! Note that vertical interpolation is linear so that qtot can be
  ! decomposed into qv, ql, qf and qr.
  cc(:,1)=rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)
  bb(:,1)=1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)
  aa(:,2:kl-1)=qq(:,2:kl-1)+ddts*mflx(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)
  cc(:,2:kl-1)=rr(:,2:kl-1)-ddts*mflx(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1)
  bb(:,2:kl-1)=1.-qq(:,2:kl-1)-rr(:,2:kl-1)+ddts*(mflx(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)                   &
                                                 -mflx(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1))
  aa(:,kl)=qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
  bb(:,kl)=1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)


  ! theta vertical mixing
  dd(:,1)=thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                    &
                                 +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                       &
                           +ddts*rhos*wt0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=thetal(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*tlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)      &
                                           +mflx(:,2:kl-1)*tlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)          &
                                           -mflx(:,2:kl-1)*tlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)     &
                                           -mflx(:,3:kl)*tlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=thetal(:,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                        &
                                   +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm  
  wthlflux(:,1)=wt0(:)
  wthlflux(:,2:kl)=-kmo(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)                    &
                   +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                         &
                   +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wthl(:,1:kl-1)=-kmo(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)                      &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                           &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! qv (part of qtot) vertical mixing
  dd(:,1)=qvg(:,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                                          &
                           +ddts*rhos*wq0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=qvg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qvup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qvup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qvup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qvup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qvg(:,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm
  wqvflux(:,1)=wq0(:)
  wqvflux(:,2:kl)=-kmo(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)                          &
                  +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                            &
                  +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wqv(:,1:kl-1)=-kmo(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)                            &
                 +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                 +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! ql (part of qtot) vertical mixing
  dd(:,1)=qlg(:,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qlg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qlg(:,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm
  wqlflux(:,1)=0.
  wqlflux(:,2:kl)=-kmo(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)                           &
                  +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                  +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wql(:,1:kl-1)=-kmo(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)                             &
                 +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! qf (part of qtot) vertical mixing
  dd(:,1)=qfg(:,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qfg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qfg(:,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm
  wqfflux(:,1)=0.
  wqfflux(:,2:kl)=-kmo(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)                           &
                  +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                  +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wqf(:,1:kl-1)=-kmo(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)                             &
                 +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! cloud fraction vertical mixing
  dd(:,1) = stratcloud(:,1) - ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                 &
                                   +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1) = stratcloud(:,2:kl-1) + ddts*(mflx(:,1:kl-2)*cfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)   &
                                             +mflx(:,2:kl-1)*cfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                             -mflx(:,2:kl-1)*cfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                             -mflx(:,3:kl)*cfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl) = stratcloud(:,kl) + ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                     &
                                     +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(stratcloud,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  stratcloud(:,:) = min( max( stratcloud(:,:), 0. ), 1. )

  
  ! momentum vertical mixing
  aa(:,2:kl)   = qq(:,2:kl)
  cc(:,1:kl-1) = rr(:,1:kl-1)
  bb(:,1) = 1. - cc(:,1) + ddts*rhos*cduv/(rhoa(:,1)*dz_fl(:,1)) ! implicit  
  bb(:,2:kl-1) = 1. - aa(:,2:kl-1) - cc(:,2:kl-1)
  bb(:,kl) = 1. - aa(:,kl)
  dd(:,1:kl) = uo(:,1:kl)
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  dd(:,1:kl) = vo(:,1:kl)
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  
  
  ! update surface momentum flux
  ustar = sqrt(cduv*sqrt(uo(:,1)**2+vo(:,1)**2))
  ustar_ave = ustar_ave + ustar/real(mcount)

  ! Account for phase transistions
  do concurrent (k = 1:kl)
    ! Check for -ve values  
    qt(:) = max( qfg(:,k) + qlg(:,k) + qvg(:,k), 0. )
    qc(:) = max( qfg(:,k) + qlg(:,k), 0. )
    
    qfg(:,k) = max( qfg(:,k), 0. )
    qlg(:,k) = max( qlg(:,k), 0. )
    
    ! account for saturation
    theta(:,k) = thetal(:,k) + sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
    temp(:) = theta(:,k)/sigkap(k)
    templ(:) = thetal(:,k)/sigkap(k)
    pres(:) = ps(:)*sig(k)
    where ( qfg(:,k)>1.e-12 )
      fice = min( qfg(:,k)/(qfg(:,k)+qlg(:,k)), 1. )
    elsewhere
      fice = 0.
    end where
    call getqsat(qsat(:,k),templ(:),pres(:),fice)
    lx(:) = lv + lf*fice(:)
    dqsdt(:) = qsat(:,k)*lx(:)/(rv*templ(:)**2)
    al(:) = cp/(cp+lx(:)*dqsdt(:))
    qc(:) = max( al(:)*(qt(:) - qsat(:,k)), qc(:) )
    where ( temp(:)>=tice )
      qfg(:,k) = max( fice(:)*qc(:), 0. )  
      qlg(:,k) = max( qc(:) - qfg(:,k), 0. )
    end where
    
    qvg(:,k) = max( qt(:) - qfg(:,k) - qlg(:,k), 0. )
    theta(:,k) = thetal(:,k) + sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
    where ( qlg(:,k)+qfg(:,k)>1.E-12 )
      stratcloud(:,k) = max( stratcloud(:,k), 1.E-8 )
    end where
  end do  
  
#ifdef scm
  uwflux(:,1) = 0.
  uwflux(:,2:kl) = -kmo(:,1:kl-1)*(uo(:,2:kl)-uo(:,1:kl-1))/dz_hl(:,1:kl-1)
  vwflux(:,1) = 0.
  vwflux(:,2:kl) = -kmo(:,1:kl-1)*(vo(:,2:kl)-vo(:,1:kl-1))/dz_hl(:,1:kl-1)
  do concurrent (k = 1:kl)
    wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                                 *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
  end do
#endif
  
end do

#ifdef scm
buoyproduction = ppb
shearproduction = pps
totaltransport = ppt
#endif

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plume rise model
    
subroutine plumerise(imax_p,lmask,cm12,                                &
                     zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,           &
                     zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,         &
                     stratcloud,km,ustar,wt0,wq0,wtv0,ps,              &
                     sig,sigkap,tke,eps,imax)

integer, intent(in) :: imax
integer, intent(in) :: imax_p
integer k, ktopmax
real, dimension(imax,kl), intent(inout) :: mflx, tlup, qvup, qlup, qfup, cfup
real, dimension(imax), intent(inout) :: zi, wstar
real, dimension(:,:), intent(in) :: theta, qvg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: zz, thetal, thetav, km 
real, dimension(imax,kl-1), intent(in) :: dz_hl
real, dimension(imax), intent(in) :: ustar, wt0, wq0, wtv0, ps
real, dimension(kl), intent(in) :: sig, sigkap
real, intent(in) :: cm12
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(in) :: eps
real, dimension(imax_p,kl) :: mflx_p, tlup_p, qvup_p, qlup_p, qfup_p, cfup_p
real, dimension(imax_p,kl) :: zz_p
real, dimension(imax_p,kl-1) :: dz_hl_p
real, dimension(imax_p,kl) ::  w2up, nn, cxup
real, dimension(imax_p) :: zi_p, tke_p, eps_p, km_p, thetal_p, theta_p, thetav_p
real, dimension(imax_p) :: qvg_p, qlg_p, qfg_p, stratcloud_p 
real, dimension(imax_p) :: ustar_p, wstar_p, wt0_p, wq0_p, wtv0_p, ps_p
real, dimension(imax_p) :: tke1, dzht, ent, templ, pres, upf, qxup, qupsat
real, dimension(imax_p) :: tempd, fice, lx, qcup, dqsdt, al, xp, as, bs, cs
real, dimension(imax_p) :: thup, tvup, qtup
logical, dimension(imax), intent(in) :: lmask
#ifdef offline
real, dimension(imax_p) :: dtr
#endif

if ( imax_p==0 ) then
  zi = zz(:,1)
  mflx(:,:) = 0.
  tlup(:,:) = thetal(:,:)
  qvup(:,:) = qvg(:,:)
  qlup(:,:) = qlg(:,:)
  qfup(:,:) = qfg(:,:)
  cfup(:,:) = stratcloud(:,:)
  return
end if

! packing
zi_p    = pack( zi, lmask )
ustar_p = pack( ustar, lmask )
wstar_p = pack( wstar, lmask )
wt0_p   = pack( wt0, lmask )
wq0_p   = pack( wq0, lmask )
wtv0_p  = pack( wtv0, lmask )
ps_p    = pack( ps, lmask )
do concurrent (k = 1:kl)
  mflx_p(:,k) = 0.
  tlup_p(:,k) = pack( thetal(:,k), lmask )
  qvup_p(:,k) = pack( qvg(:,k), lmask )
  qlup_p(:,k) = pack( qlg(:,k), lmask )
  qfup_p(:,k) = pack( qfg(:,k), lmask )
  cfup_p(:,k) = pack( stratcloud(:,k), lmask )
  zz_p(:,k) = pack( zz(:,k), lmask )
end do
do concurrent (k = 1:kl-1)
  dz_hl_p(:,k) = pack( dz_hl(:,k), lmask )  
end do

! Initialise updraft
tke1 = cm12*ustar_p*ustar_p + ce3*wstar_p*wstar_p
tke1 = max(tke1, mintke)
w2up = 0.
nn = 0.
dzht = zz_p(:,1)
ktopmax = 0

! Entrainment and detrainment rates
ent = entfn(zz_p(:,1),zi_p)

theta_p  = pack( theta(:,1), lmask )
thetal_p = pack( thetal(:,1), lmask )
thetav_p = pack( thetav(:,1), lmask )
qvg_p    = pack( qvg(:,1), lmask )
qlg_p    = pack( qlg(:,1), lmask )
qfg_p    = pack( qfg(:,1), lmask )
stratcloud_p = pack ( stratcloud(:,1), lmask )

! first level -----------------
! initial thermodynamic state
! split qtot into components (conservation of thetal and qtot is maintained)
tlup_p(:,1) = thetal_p + be*wt0_p/sqrt(tke1)       ! Hurley 2007
qvup_p(:,1) = qvg_p    + be*wq0_p/sqrt(tke1)       ! Hurley 2007
qlup_p(:,1) = qlg_p
qfup_p(:,1) = qfg_p
cfup_p(:,1) = stratcloud_p
! calculate thermodynamic variables assuming no condensation
qtup = qvup_p(:,1) + qlup_p(:,1) + qfup_p(:,1)                    ! qtot,up
thup = tlup_p(:,1) + sigkap(1)*(lv*qlup_p(:,1)+ls*qfup_p(:,1))/cp ! theta,up
tvup = thup + theta_p*0.61*qtup                                   ! thetav,up
templ = tlup_p(:,1)/sigkap(1)                                     ! templ,up
! update updraft velocity and mass flux
nn(:,1) = grav*be*wtv0_p/(thetav_p*sqrt(tke1))          ! Hurley 2007
w2up(:,1) = 2.*dzht*b2*nn(:,1)/(1.+2.*dzht*b1*ent)      ! Hurley 2007
cxup(:,1) = 0.


! updraft with condensation
do k = 2,kl
  dzht = dz_hl_p(:,k-1)
  ! Entrainment and detrainment rates
  ent = entfn(zz_p(:,k), zi_p(:))
  theta_p  = pack( theta(:,k), lmask )
  thetal_p = pack( thetal(:,k), lmask )
  thetav_p = pack( thetav(:,k), lmask )
  qvg_p = pack( qvg(:,k), lmask )
  qlg_p = pack( qlg(:,k), lmask )
  qfg_p = pack( qfg(:,k), lmask )
  stratcloud_p = pack( stratcloud(:,k), lmask )
  tke_p = pack( tke(:,k), lmask )
  eps_p = pack( eps(:,k), lmask )
  km_p  = pack( km(:,k), lmask )
  where ( w2up(:,k-1)>0. )
    ! entrain air into plume
    ! split qtot into components (conservation of qtot is maintained)
    tlup_p(:,k) = (tlup_p(:,k-1)+dzht*ent*thetal_p)/(1.+dzht*ent)
    qvup_p(:,k) = (qvup_p(:,k-1)+dzht*ent*qvg_p   )/(1.+dzht*ent)
    qlup_p(:,k) = (qlup_p(:,k-1)+dzht*ent*qlg_p   )/(1.+dzht*ent)
    qfup_p(:,k) = (qfup_p(:,k-1)+dzht*ent*qfg_p   )/(1.+dzht*ent)
    cfup_p(:,k) = (cfup_p(:,k-1)+dzht*ent*stratcloud_p )/(1.+dzht*ent)
  elsewhere
    tlup_p(:,k) = thetal_p
    qvup_p(:,k) = qvg_p
    qlup_p(:,k) = qlg_p
    qfup_p(:,k) = qfg_p
    cfup_p(:,k) = stratcloud_p
  end where
  ! calculate conserved variables
  qtup = qvup_p(:,k) + qlup_p(:,k) + qfup_p(:,k)    ! qtot,up
  templ = tlup_p(:,k)/sigkap(k)                     ! templ,up
  pres = ps_p(:)*sig(k)
  where ( qfup_p(:,k)>1.e-12 )
    fice = min( qfup_p(:,k)/(qfup_p(:,k)+qlup_p(:,k)), 1. )
  elsewhere
    fice = 0.
  end where
  call getqsat(qupsat,templ,pres,fice)
  where ( qtup > qupsat )
    qxup = qupsat
    cxup(:,k) = 1.
  elsewhere
    qxup = qtup
    cxup(:,k) = 0.
  end where
  where ( w2up(:,k-1)<=0. )
    qxup = qtup
    cxup(:,k) = 0.
  end where
  lx = lv + lf*fice
  dqsdt = qupsat*lx/(rv*templ**2)
  al = cp/(cp+lx*dqsdt)
  qcup = max(al*(qtup-qxup), 0.)                             ! qcondensate,up after redistribution
  qcup = min(qcup, qcmf)                                     ! limit condensation with simple autoconversion
  thup = tlup_p(:,k) + sigkap(k)*qcup*lx/cp                  ! theta,up after redistribution
  tvup = thup + theta_p*(0.61*qxup-qcup)                     ! thetav,up after redistribution
  where ( w2up(:,k-1)>0. )
    nn(:,k) = grav*(tvup-thetav_p)/thetav_p                          ! calculate buoyancy
    w2up(:,k) = (w2up(:,k-1)+2.*dzht*b2*nn(:,k))/(1.+2.*dzht*b1*ent) ! update updraft velocity
  elsewhere
    nn(:,k) = 0.  
    w2up(:,k) = 0.
  end where
  ! test if maximum plume height is reached
  where ( w2up(:,k)<=0. .and. w2up(:,k-1)>0. )
    as = 2.*b2*(nn(:,k)-nn(:,k-1))/dzht
    bs = 2.*b2*nn(:,k-1)
    cs = w2up(:,k-1)
    xp = -2.*cs/(bs-sqrt(max(bs*bs-4.*as*cs,0.)))
    xp = min(max(xp,0.),dzht)
    zi_p(:) = xp + zz_p(:,k-1)
  end where
  ktopmax = k - 1
  if ( all(w2up(:,k)<=0.) ) exit
end do
          
thetav_p = pack( thetav(:,1), lmask )
wstar_p = (grav*min(zi_p,zimax)*max(wtv0_p,0.)/thetav_p)**(1./3.)
          
! update mass flux
mflx_p(:,1) = m0*sqrt(max(w2up(:,1), 0.))
do concurrent (k = 2:ktopmax)
  dzht = dz_hl_p(:,k-1)
  upf = mflx_p(:,k-1)/sqrt(max(w2up(:,k-1), 1.e-8))
  where ( w2up(:,k)>0. )
    mflx_p(:,k) = (1.-cxup(:,k))*m0*sqrt(max(w2up(:,k), 0.))         &
                + cxup(:,k)*mflx_p(:,k-1)/(1.+dzht*(dtrc0-entc0))
    mflx_p(:,k) = min( mflx_p(:,k), upf*sqrt(max(w2up(:,k), 0.)) )
  elsewhere
    mflx_p(:,k) = 0.
  end where
end do

! unpacking
zi = unpack( zi_p, lmask, zz(:,1) )
tke(:,1) = unpack( tke1, lmask, tke(:,1) )
wstar = unpack( wstar_p, lmask, wstar )
do concurrent (k = 1:ktopmax)
  mflx(:,k) = unpack( mflx_p(:,k), lmask, 0. )
  tlup(:,k) = unpack( tlup_p(:,k), lmask, thetal(:,k) )
  qvup(:,k) = unpack( qvup_p(:,k), lmask, qvg(:,k) )
  qlup(:,k) = unpack( qlup_p(:,k), lmask, qlg(:,k) )
  qfup(:,k) = unpack( qfup_p(:,k), lmask, qfg(:,k) )
  cfup(:,k) = unpack( cfup_p(:,k), lmask, stratcloud(:,k) )
end do

#ifdef offline
do concurrent (k = 1:ktopmax)
  mf(:,k) = unpack( mflx_p(:,k), lmask, 0. )
  w_up(:,k) = unpack( sqrt(w2up(:,k)), lmask, 0. )
  tl_up(:,k) = unpack( tlup_p(:,k), lmask, thetal(:,k) )
  qv_up(:,k) = unpack( qvup_p(:,k), lmask, qvg(:,k) )
  ql_up(:,k) = unpack( qlup_p(:,k), lmask, qlg(:,k) )
  qf_up(:,k) = unpack( qfup_p(:,k), lmask, qfg(:,k) )
  cf_up(:,k) = unpack( cxup(:,k)*min(mflx_p(:,k)/sqrt(max(w2up(:,k),1.e-8)),1.), lmask, 0. )
end do

! Entrainment and detrainment rates
dzht = zz_p(:,1)
ent = entfn(zz_p(:,1),zi_p(:))
dtr = -1./dzht + ent
dtr = max( dtr, 0. )
ents(:,1) = unpack( ent, lmask, 0. )
dtrs(:,1) = unpack( dtr, lmask, 0. )
do concurrent (k = 2:ktopmax)
  dzht = dz_hl_p(:,k-1)
  ! Entrainment and detrainment rates
  where ( w2up(:,k)>0. )
    ent = entfn(zz_p(:,k),zi_p(:))
    dtr = mflx_p(:,k-1)/(mflx_p(:,k)*dzht) - 1./dzht + ent
    dtr = max( dtr, 0. )
  elsewhere
    ent = 0.
    dtr = 0.
  end where
  ents(:,k) = unpack( ent, lmask, 0. )
  dtrs(:,k) = unpack( dtr, lmask, 0. )
end do
#endif

return
end subroutine plumerise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stablepbl(lmask,zi,zzh,dz_hl,thetav,kmo,imax)

implicit none

integer, intent(in) :: imax
integer i, k
real, dimension(imax), intent(inout) :: zi
real, dimension(imax,kl-1), intent(in) :: dz_hl
real, dimension(imax,kl), intent(in) :: zzh, thetav, kmo 
real, dimension(kl) :: wpv_flux
real xp
logical, dimension(imax), intent(in) :: lmask

do concurrent (i = 1:imax)
  if ( lmask(i) ) then
    !wpv_flux is calculated at half levels  
    wpv_flux(1) = -kmo(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1) !+gamt_hl(i,1)
    do k = 2,kl-1
      wpv_flux(k) = -kmo(i,k)*(thetav(i,k+1)-thetav(i,k))/dz_hl(i,k) !+gamt_hl(i,k)
      if ( wpv_flux(k)*wpv_flux(1)<0. ) then ! change in sign
        xp = (0.05*wpv_flux(1)-wpv_flux(k-1))/(wpv_flux(k)-wpv_flux(k-1))
        xp = min( max( xp, 0. ), 1. )
        zi(i) = zzh(i,k-1) + xp*(zzh(i,k)-zzh(i,k-1))
        exit
      else if ( abs(wpv_flux(k))<0.05*abs(wpv_flux(1)) ) then
        xp = (0.05*abs(wpv_flux(1))-abs(wpv_flux(k-1)))/(abs(wpv_flux(k))-abs(wpv_flux(k-1)))
        xp = min( max( xp, 0. ), 1. )
        zi(i) = zzh(i,k-1) + xp*(zzh(i,k)-zzh(i,k-1))
        exit
      end if    
    end do    
  end if
end do

return
end subroutine stablepbl
                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

pure subroutine thomas(outdat,aai,bbi,cci,ddi,imax)

implicit none

integer, intent(in) :: imax
real, dimension(:,2:), intent(in) :: aai
real, dimension(:,:), intent(in) :: bbi,ddi
real, dimension(:,:), intent(in) :: cci
real, dimension(:,:), intent(out) :: outdat
real, dimension(imax,size(outdat,2)) :: cc,dd
real, dimension(imax) :: n
integer k,klin

klin=size(outdat,2)
cc(:,1)=cci(:,1)/bbi(:,1)
dd(:,1)=ddi(:,1)/bbi(:,1)

do k=2,klin-1
  n(:)=bbi(:,k)-cc(:,k-1)*aai(:,k)
  cc(:,k)=cci(:,k)/n(:)
  dd(:,k)=(ddi(:,k)-dd(:,k-1)*aai(:,k))/n(:)
end do
n(:)=bbi(:,klin)-cc(:,klin-1)*aai(:,klin)
dd(:,klin)=(ddi(:,klin)-dd(:,klin-1)*aai(:,klin))/n(:)
outdat(:,klin)=dd(:,klin)
do k=klin-1,1,-1
  outdat(:,k)=dd(:,k)-cc(:,k)*outdat(:,k+1)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

pure subroutine getqsat(qsat,templ,ps,fice)

implicit none

real, dimension(:), intent(in) :: templ
real, dimension(size(templ)), intent(in) :: ps, fice
real, dimension(size(templ)), intent(out) :: qsat
real, dimension(size(templ)) :: estafi, tdiff, tdiffx, rx, rxx
real, dimension(size(templ)) :: qsatl, qsati, deles
integer, dimension(size(templ)) :: ix, ixx

tdiff = min(max( templ(:)-123.16, 0.), 219.)
rx = tdiff - aint(tdiff)
ix = int(tdiff)
estafi = (1.-rx)*tablei(ix) + rx*tablei(ix+1)
qsati = 0.622*estafi/max(ps(:)-estafi,0.1)

tdiffx = min(max( templ(:)-273.1, -40.), 1.)
rxx = tdiffx - aint(tdiffx)
ixx = int(tdiffx)
deles = (1.-rxx)*esdiff(ixx) + rxx*esdiff(ixx+1)
qsatl = qsati + 0.622*deles/ps

qsat = fice*qsati + (1.-fice)*qsatl

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Estimate ice fraction
!
!subroutine getfice(fice,templ,qlg,qfg)
!
!implicit none
!
!real, dimension(:), intent(in) :: templ
!real, dimension(size(templ)), intent(in) :: qlg, qfg
!real, dimension(size(templ)), intent(out) :: fice
!real, dimension(size(templ)) :: temp
!
!temp = templ + (lv*qlg+ls*qfg)/cp
!
!where ( temp>=tfrz )
!  fice = 0.
!elsewhere ( temp>=tice .and. qfg>1.e-12 )
!  fice = min( qfg/(qfg+qlg), 1. )
!elsewhere ( temp>=tice )
!  fice = 0.
!elsewhere
!  fice = 1.
!end where
!
!return
!end subroutine getfice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

pure subroutine updatekmo(kmo,km,fzhl,imax)

implicit none

integer, intent(in) :: imax
real, dimension(:,:), intent(in) :: km
real, dimension(imax,kl-1), intent(in) :: fzhl
real, dimension(imax,kl), intent(out) :: kmo

kmo(:,1:kl-1)=km(:,1:kl-1)+fzhl(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))
! These terms are never used
kmo(:,kl)=0.

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

pure function entfn(zht,zi) result(ans)

implicit none

real, dimension(:), intent(in) :: zht, zi
real, dimension(size(zht)) :: ans

!ans=0.002                                            ! Angevine (2005)
!ans=2./max(100.,zi)                                  ! Angevine et al (2010)
!ans=1./zht                                           ! Siebesma et al (2003)
!ans=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))    ! Soares et al (2004)
ans = max( ent0/max( zht, 1. ) + ent1/max( min(zi,zimax)-zht, ezmin ), ent_min )

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral detrainment

!real function dtrfn(zht,zi,rat)
!
!implicit none
!
!real, intent(in) :: zht,zi,rat
!
!!dtrfn=ent+0.05/max(zi-zht,zmin)   ! Angevine et al (2010)
!dtrfn=rat/max(zi-zht,1.)+ent1/max(zi-zht,ezmin)
!
!! results in analytic solution
!!mflx(k)=A*(zht**ent0)*((zi-zht)**rat)
!
!
!return
!end function dtrfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phi_m for atmospheric stability

pure subroutine calc_phi(phim,z_on_l)

implicit none

real, dimension(:), intent(in) :: z_on_l
real, dimension(:), intent(out) :: phim

select case(stabmeth)
  case(0)
    where ( z_on_l<0. )
      phim = (1.-16.*z_on_l)**(-0.25)
    elsewhere
      phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
    end where
  case(1)
    where ( z_on_l<0. )
      phim = (1.-16.*z_on_l)**(-0.25)
    elsewhere (z_on_l<=0.4)
      phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
    elsewhere
      phim = aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar
    end where
end select

return
end subroutine calc_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Diagnose cloud fraction
!
!subroutine calc_smithcloud(cfup,qxup,sigqtup,qtup,qupsat)
!
!implicit none
!
!real, dimension(:), intent(in) :: sigqtup,qtup,qupsat
!real, dimension(:), intent(out) :: cfup, qxup
!real, dimension(size(cfup)) :: rng, dqdash
!
!rng = sqrt(6.)*sigqtup           ! variance of triangle distribution
!dqdash = (qtup-qupsat)/rng  ! scaled variance
!where ( dqdash<-1. )
!  ! gridbox all unsaturated
!  qxup = qtup
!  cfup = 0.
!elsewhere ( dqdash<0. )
!  ! gridbox minority saturated
!  qxup = qtup + 0.5*rng*(-1./3.-dqdash-dqdash**2-1./3.*dqdash**3)
!  cfup = 0.5*(dqdash+1.)**2
!elsewhere ( dqdash<1. )
!  ! gridbox majority saturated
!  qxup = qtup + 0.5*rng*(-1./3.-dqdash-dqdash**2+1./3.*dqdash**3)
!  cfup = 1. - 0.5*(dqdash-1.)**2
!elsewhere
!  ! gridbox all saturated
!  qxup = qupsat
!  cfup = 1.
!end where
!
!return
!end subroutine calc_smithcloud

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

