! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public cm0,ce0,ce1,ce2,ce3,be,ent0,ent1,entc0,ezmin,dtrc0
public m0,b1,b2,qcmf,ent_min,mfbeta
public buoymeth,maxdts,mintke,mineps,minl,maxl,stabmeth
public tkemeth
public tke_timeave_length, update_ema
public u_ema, v_ema, w_ema
public thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
public tke_ema

real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:,:), allocatable, save :: u_ema, v_ema, w_ema
real, dimension(:,:), allocatable, save :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(:,:), allocatable, save :: tke_ema

! model ED constants
real, save :: cm0      = 0.09     ! Hurley (2007) 0.09, Duynkerke (1988) 0.03, Duynkerke (1987) 0.09
real, save :: ce0      = 0.69     ! Hurley (2007) 0.69, Duynkerke (1988) 0.42, Duynkerke (1987) 0.77
real, save :: ce1      = 1.46
real, save :: ce2      = 1.83
real, save :: ce3      = 0.45     ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, save :: mintke   = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps   = 1.E-11   ! min value for eps (1.0e-6 in TAPM)
real, save :: minl     = 1.       ! min value for L   (5. in TAPM)
real, save :: maxl     = 1000.    ! max value for L   (500. in TAPM)
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
! wind gusts
real, parameter :: cs1 = 2.2
real, parameter :: cs2 = 1.63
real, parameter :: cs3 = 0.73
real, parameter :: cw1 = 1.
real, parameter :: cw2 = 0.24
real, parameter :: cw3 = 0.

#ifdef GPUPHYSICS
!!$acc declare create(cm0,ce0,ce1,ce2,maxl,minl,mintke,mineps,buoymeth,tkemeth)
!!$acc declare create(ce3,mfbeta,stabmeth,m0,qcmf,be,b1,b2)
!!$acc declare create(maxdts,dtrc0,entc0,ent_min,ent0,ent1,ezmin)
#endif

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
real, parameter :: a_1 = 1.
real, parameter :: b_1 = 2./3.
real, parameter :: c_1 = 5.
real, parameter :: d_1 = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

! Time averaging
real, save :: tke_timeave_length = 0. ! Time period for averaging source terms (Ps, Pb, Pt) in seconds
                                      ! 0 indicates alpha=2/3

#ifdef GPUPHYSICS
!!$acc declare create(tke_timeave_length)
#endif

real, dimension(0:220), parameter :: tablei = &
(/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                & !-146C
   6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                             & !-141C
   36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                          & !-136C
   0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,  & !-131C
   0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,  & !-126C
   0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,   & !-121C
   0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,       & !-116C
   0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,         & !-111C
   0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,            & !-106C
   0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,             & !-101C
   .001403, .001719, .002101, .002561, .003117, .003784,             & !-95C
   .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658, & !-87C
   .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577, & !-78C
   .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032,   & !-69C
   .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080,    & !-60C
   1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,    & !-51C
   3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,            & !-43C
   10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,    & !-34C
   27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85, & !-24C 
   77.09, 85.02, 93.70, 103.06, 113.40, 124.68, 136.98, 150.39,      & !-16C
   164.99, 180.88, 198.16, 216.94, 237.34, 259.47, 283.49, 309.51,   & !-8C
   337.71, 368.23, 401.25, 436.96, 475.54, 517.21, 562.19, 610.70,   & !0C
   656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,   & !8C
   1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,   & !16C
   1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,   & !24C
   3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,   & !32C
   5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,   & !40C
   7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,         & !47C
   11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,    & !54C
   15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,    & !61C
   21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,    & !68C
   29845.0, 31169.0 /)                                                 !70C 
real, dimension(-40:2), parameter :: esdiff = &
(/ 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,  &
   13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61, &
   22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27, &
   26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65, &
   0.08, 0., 0. /)


interface solve_sherman_morrison
  module procedure solve_sherman_morrison_2, solve_sherman_morrison_3
end interface

interface thomas
  module procedure thomas1, thomas2
end interface

interface getqsat
  module procedure getqsat2, getqsat3
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(shear(ifull,kl))
allocate(u_ema(ifull,kl), v_ema(ifull,kl), w_ema(ifull+iextra,kl))
allocate(thetal_ema(ifull,kl), qv_ema(ifull,kl), ql_ema(ifull,kl), qf_ema(ifull,kl))
allocate(cf_ema(ifull,kl))
allocate(tke_ema(ifull,kl))

tke(1:ifull+iextra,1:kl)=mintke
eps(1:ifull+iextra,1:kl)=mineps
shear(1:ifull,1:kl)=0.
u_ema(1:ifull,1:kl)=0.
v_ema(1:ifull,1:kl)=0.
w_ema(1:ifull,1:kl)=0.
thetal_ema(1:ifull,1:kl)=0.
qv_ema(1:ifull,1:kl)=0.
ql_ema(1:ifull,1:kl)=0.
qf_ema(1:ifull,1:kl)=0.
cf_ema(1:ifull,1:kl)=0.
tke_ema(1:ifull,1:kl)=0.

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux
! mode=2 combined atmosphere and ocean with mass flux
! mode=3 combined atm-ocn and no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,stratcloud,ua,va,zi,fg,eg,cduv,ps,zz,zzh,sig,rhos, &
                  ustar_ave,dt,qgmin,mode,tke,eps,shear,dx,thetal_ema,qv_ema,ql_ema,       &
                  qf_ema,cf_ema,tke_ema,                                                   &
#ifdef scm
                  wthflux,wqvflux,uwflux,vwflux,mfout,buoyproduction,                      &
                  shearproduction,totaltransport,                                          &
#endif
                  deptho_fl,deptho_hl,cd_water,cdh_water,cdbot_water,wt0rad_o,             &
                  wt0melt_o,wt0eg_o,wt0fb_o,ws0_o,ws0subsurf_o,w_t,w_s,w_u,w_v,rad_o,ibot, &
                  i_u,i_v,fracice,icefg_a,imass,cd_ice,cdbot_ice,w_tke,w_eps,              &
                  w_t_ema,w_s_ema,zo_o,shear_o,                                            &
                  land,ugs_var,sigkap,imax,kl,wlev)
#ifdef GPUPHYSICS
!!$acc routine vector                  
#endif

use mlo, only : mlo_updatekm, wrtrho, wrtemp, cp0, minsfc
     
implicit none

integer, intent(in) :: mode, imax, kl, wlev
integer k, iq
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(:,:), intent(inout) :: theta,stratcloud,ua,va
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg
real, dimension(:,:), intent(out) :: kmo
real, dimension(:,:), intent(in) :: zz, zzh
real, dimension(:,:), intent(in) :: shear
real, dimension(:,:), intent(inout) :: tke
real, dimension(:,:), intent(inout) :: eps
real, dimension(:,:), intent(inout) :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(:,:), intent(inout) :: tke_ema
real, dimension(:), intent(inout) :: zi, fg, eg
real, dimension(:), intent(in) :: cduv, ps, rhos, dx
real, dimension(:), intent(in) :: sig, sigkap
real, dimension(:), intent(out) :: ustar_ave, ugs_var
real, dimension(imax,kl) :: km, thetal
real, dimension(imax,kl) :: rhoa
real, dimension(imax,kl) :: tlup,qvup,qlup,qfup
real, dimension(imax,kl) :: cfup, mflx
real, dimension(imax,2:kl) :: idzm
real, dimension(imax,kl-1) :: idzp
real, dimension(imax,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(imax,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(imax,kl-1) :: fzzh
real, dimension(imax,kl) :: templ, fice, qsat, pres
real, dimension(imax,kl) :: ppb, pps, ppt
real, dimension(imax) :: wt0, wq0, wtv0, thetav1
real, dimension(imax) :: wstar,z_on_l,phim
real, dimension(imax) :: ustar, fg_ave
real, dimension(imax) :: zi_save, cgmap
real tff, cm12, cm34, ddts
real temp, lx, dqsdt, al, wdash_sq, qc, qt
real nstep, alpha, zturb, rhoahl
logical, dimension(imax) :: mask

integer, dimension(:), intent(in) :: ibot
real, dimension(:,:), intent(in) :: rad_o
real, dimension(:,:), intent(in) :: deptho_fl
real, dimension(:,:), intent(in) :: deptho_hl
real, dimension(:), intent(in) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o, wt0eg_o
real, dimension(:,:), intent(inout) :: w_t, w_s, w_u, w_v
real, dimension(:,:), intent(inout) :: w_t_ema, w_s_ema
real, dimension(:,2:), intent(in) :: shear_o
real, dimension(:), intent(in) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o
real, dimension(:), intent(in) :: fracice, imass, cd_ice, cdbot_ice, zo_o
real, dimension(:,:), intent(inout) :: w_tke, w_eps
real, dimension(:), intent(inout) :: i_u, i_v
logical, dimension(:), intent(in) :: land
real, dimension(imax,wlev) :: km_o, ks_o, rhs_o, gammas_o
real, dimension(imax) :: ustar_o, wu0_o, wv0_o
real deptho_dz

#ifdef scm
real, dimension(:,:), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(:,:), intent(out) :: buoyproduction, shearproduction
real, dimension(:,:), intent(out) :: totaltransport
real, dimension(:,:), intent(out) :: mfout
#endif

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

do k = 1,kl
  do iq = 1,imax
    ! Impose limits on tke and eps after being advected by the host model
    tke(iq,k) = max(tke(iq,k), mintke)
    tff       = cm34*tke(iq,k)*sqrt(tke(iq,k))
    eps(iq,k) = min(eps(iq,k), tff/minl)
    eps(iq,k) = max(eps(iq,k), tff/maxl, mineps)
  
    ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
    pres(iq,k) = ps(iq)*sig(k) ! pressure
    ! Density must be updated when dz is updated so that rho*dz is conserved
    thetav1(iq) = theta(iq,k)*(1.+0.61*qvg(iq,k)-qlg(iq,k)-qfg(iq,k))
    rhoa(iq,k) = sigkap(k)*pres(iq,k)/(rd*thetav1(iq))

    ! Calculate to thetal as it is the conserved variable
    thetal(iq,k) = theta(iq,k) - sigkap(k)*(lv*qlg(iq,k)+ls*qfg(iq,k))/cp
  
    ! Calculate first approximation to diffusion coeffs
    km(iq,k) = cm0*tke(iq,k)**2/eps(iq,k)
  end do
end do
  
! Fraction for interpolation from full levels to half levels
fzzh(:,1:kl-1) = (zzh(:,1:kl-1)-zz(:,1:kl-1))/(zz(:,2:kl)-zz(:,1:kl-1))

! Calculate dz at half levels
dz_hl(:,1:kl-1) = max( zz(:,2:kl) - zz(:,1:kl-1), 1. )

! Calculate dz at full levels
dz_fl(:,1)    = zzh(:,1)
dz_fl(:,2:kl) = zzh(:,2:kl) - zzh(:,1:kl-1)
dz_fl(:,1:kl) = max( dz_fl(:,1:kl), 1. )

! interpolate diffusion coeff to half levels
kmo(:,1:kl-1)=km(:,1:kl-1)+fzzh(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))

do k = 1,kl-1
  do iq = 1,imax  
    ! eddy diffusion terms to account for air density with level thickness
    rhoahl = rhoa(iq,k)+fzzh(iq,k)*(rhoa(iq,k+1)-rhoa(iq,k))
    idzm(iq,k+1)   = rhoahl/(rhoa(iq,k+1)*dz_fl(iq,k+1))
    idzp(iq,k) = rhoahl/(rhoa(iq,k)*dz_fl(iq,k))
  end do
end do

ustar_ave(:) = 0.
fg_ave(:) = 0.
ugs_var(:) = 0.


! Main loop to prevent time splitting errors
nstep = max( tke_timeave_length/dt, 1. ) ! this is a real value
alpha = 2./(nstep + 1.) 
mcount = int(dt/(min(maxdts,12./max(m0,0.1))+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  do k = 1,kl
    do iq = 1,imax
      ! Set-up thermodynamic variables temp, theta_v, qtot and saturated mixing ratio
      theta(iq,k) = thetal(iq,k) + sigkap(k)*(lv*qlg(iq,k)+ls*qfg(iq,k))/cp
    end do
  end do

  thetav1(:) = theta(:,1)*(1.+0.61*qvg(:,1)-qlg(:,1)-qfg(:,1))
  
  ! Calculate surface fluxes
  wt0(:) = fg(:)/(rhos(:)*cp)  ! theta flux
  wq0(:) = eg(:)/(rhos(:)*lv)  ! qtot flux  
  wtv0(:) = wt0(:) + theta(:,1)*0.61*wq0(:) ! thetav flux
  
  ! Update momentum flux
  ustar(:) = sqrt(cduv(:)*sqrt(ua(:,1)**2+va(:,1)**2))  
  wstar(:) = max(grav*zi(:)*max(wtv0(:),0.)/thetav1(:),1.e-10)**(1./3.)   

  ! update surface boundary condition
  tke(:,1) = cm12*ustar(:)**2 + ce3*wstar(:)**2
  tke(:,1) = max( tke(:,1), mintke )
  

  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  if ( mode/=1 .and. mode/=3 ) then ! mass flux

      zi_save(:) = zi(:)  

      ! plume rise model
      mask = wtv0(:)>0.

      call plumerise(mask,zi,wstar,mflx,tlup,qvup,           &
                     qlup,qfup,cfup,zz,dz_hl,thetal,qvg,     &
                     qlg,qfg,stratcloud,wt0,wq0,             &
                     ps,ustar,sig,tke,ua,va,sigkap,imax,kl)

#ifndef scm
      ! Turn off MF term if small grid spacing (mfbeta=0 implies MF is always non-zero)
      ! Based on Boutle et al 2014
      do iq = 1,imax
        zturb = 0.5*(zi_save(iq) + zi(iq))
        cgmap(iq) = 1. - tanh(mfbeta*zturb/dx(iq))*max(0.,1.-0.25*dx(iq)/zturb)
      end do
      do k = 1,kl
        do iq = 1,imax
          mflx(iq,k) = mflx(iq,k)*cgmap(iq)
        end do
      end do
#endif

  else

    mflx(:,:) = 0.
    tlup(:,:) = thetal(:,:)
    qvup(:,:) = qvg(:,:)
    qlup(:,:) = qlg(:,:)
    qfup(:,:) = qfg(:,:)
    cfup(:,:) = stratcloud(:,:)

  end if


#ifdef scm  
  do k = 1,kl-1
    mfout(:,k) = mflx(:,k)*(1.-fzzh(:,k)) &
               + mflx(:,k+1)*fzzh(:,k)
  end do  
  mfout(:,kl)=0.  
#endif


  ! calculate tke and eps boundary condition at 1st vertical level
  z_on_l(:) = -vkar*zz(:,1)*grav*wtv0(:)/(thetav1(:)*max(ustar(:)**3,1.E-10))
  z_on_l(:) = min(z_on_l(:),10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  call calc_phi(phim(:),z_on_l(:),imax)
  do iq = 1,imax
    tke(iq,1) = cm12*ustar(iq)**2+ce3*wstar(iq)**2
    eps(iq,1) = ustar(iq)**3*phim(iq)/(vkar*zz(iq,1))+grav*wtv0(iq)/thetav1(iq)
    tke(iq,1) = max( tke(iq,1), mintke )
    tff = cm34*tke(iq,1)*sqrt(tke(iq,1))
    eps(iq,1) = min( eps(iq,1), tff/minl )
    eps(iq,1) = max( eps(iq,1), tff/maxl, mineps )
  end do


  ! time averaging
  do k = 1,kl
    do iq = 1,imax
      thetal_ema(iq,k) = alpha*thetal(iq,k) + (1.-alpha)*thetal_ema(iq,k)
      qv_ema(iq,k) = alpha*qvg(iq,k) + (1.-alpha)*qv_ema(iq,k)
      ql_ema(iq,k) = alpha*qlg(iq,k) + (1.-alpha)*ql_ema(iq,k)
      qf_ema(iq,k) = alpha*qfg(iq,k) + (1.-alpha)*qf_ema(iq,k)
      cf_ema(iq,k) = alpha*stratcloud(iq,k) + (1.-alpha)*cf_ema(iq,k)
      tke_ema(iq,k) = alpha*tke(iq,k) + (1.-alpha)*tke_ema(iq,k)
    end do
  end do    


  ! Calculate shear term on full levels
  pps(:,1:kl-1) = km(:,1:kl-1)*shear(:,1:kl-1)

  ! Update TKE and eps terms
  call update_tkeeps(tke,eps,ppb,pps,ppt,qvg,qlg,qfg,thetal,         &
                     zzh,zz,idzp,idzm,wstar,zi,                      &
                     thetal_ema,qv_ema,ql_ema,qf_ema,cf_ema,tke_ema, &
                     ddts,qgmin,cm34,sigkap,imax,kl)

  
#ifdef scm
    buoyproduction(:,:) = ppb
    shearproduction(:,:) = pps
    totaltransport(:,:) = ppt
#endif
  

  if ( mode==2 .or. mode==3 ) then

    ! calculate eddy diffusivity for ocean

    ! update ocean eddy diffusivity possibly including prognostic tke and eps
    wu0_o(:) = -(1.-fracice(:))*rhoa(:,1)*cd_water(:)*(ua(:,1)-w_u(:,1))/wrtrho
    wv0_o(:) = -(1.-fracice(:))*rhoa(:,1)*cd_water(:)*(va(:,1)-w_v(:,1))/wrtrho
    wu0_o(:) = wu0_o(:) + fracice(:)*cdbot_ice(:)*(w_u(:,1)-i_u(:))
    wv0_o(:) = wv0_o(:) + fracice(:)*cdbot_ice(:)*(w_v(:,1)-i_v(:))
    ustar_o(:) = max(sqrt(sqrt(wu0_o(:)**2+wv0_o(:)**2)),1.e-6)
    call mlo_updatekm(km_o,ks_o,gammas_o,ddts,0,                         &
                      ustar_o,zo_o,deptho_fl,deptho_hl,w_t,w_s,w_u,w_v,  &
                      w_t_ema,w_s_ema,shear_o,ibot,w_tke,w_eps,imax,wlev)

    do iq = 1,imax
      deptho_dz = deptho_hl(iq,2) - deptho_hl(iq,1)
      rhs_o(iq,1) = 0.
      if ( deptho_dz>1.e-4 ) then
        rhs_o(iq,1) = ks_o(iq,2)*gammas_o(iq,2)/deptho_dz
      end if
    end do  
    do k = 2,wlev-1
      do iq = 1,imax
        rhs_o(iq,k) = 0.  
        deptho_dz = deptho_hl(iq,k+1) - deptho_hl(iq,k) 
        if ( deptho_dz>1.e-4 ) then
          rhs_o(iq,k) = (ks_o(iq,k+1)*gammas_o(iq,k+1)-ks_o(iq,k)*gammas_o(iq,k))/deptho_dz
        end if
      end do
    end do
    do iq = 1,imax    
      rhs_o(iq,wlev) = 0.
      deptho_dz = deptho_hl(iq,wlev+1) - deptho_hl(iq,wlev)
      if ( deptho_dz>1.e-4 ) then
        rhs_o(iq,wlev) = -ks_o(iq,wlev)*gammas_o(iq,wlev)/deptho_dz
      end if
    end do

    call update_coupled(thetal,qvg,qlg,qfg,stratcloud,ua,va,           &
                        tlup,qvup,qlup,qfup,cfup,fg,eg,rhos,ustar,     &
                        cduv,tke,eps,mflx,fzzh,idzp,idzm,              &
                        dz_hl,rhoa(:,1),dz_fl(:,1),                    &
                        rad_o,w_t,w_s,w_u,w_v,                         &
                        deptho_fl,deptho_hl,cd_water,                  &
                        cdh_water,cdbot_water,wt0rad_o,wt0melt_o,      &
                        wt0eg_o,icefg_a,wt0fb_o,ws0_o,ws0subsurf_o,    &
                        i_u,i_v,imass,fracice,cd_ice,cdbot_ice,        &
                        ibot,land,sigkap,km_o,ks_o,rhs_o,              &
#ifdef scm
                        wthflux,wqvflux,uwflux,vwflux,                 &
#endif
                        ddts,wrtemp,wrtrho,cp0,minsfc,                 &
                        imax,kl,wlev)

  else

    call update_atmosphere(thetal,qvg,qlg,qfg,stratcloud,ua,va,        &
                           tlup,qvup,qlup,qfup,cfup,fg,                &
                           eg,rhos,ustar,cduv,tke,eps,mflx,fzzh,idzp,  &
                           idzm,dz_hl,rhoa(:,1),dz_fl(:,1),sigkap,     &
#ifdef scm
                           wthflux,wqvflux,uwflux,vwflux,              &
#endif
                           ddts,imax,kl)

  end if

  
  ustar_ave(:) = ustar_ave(:) + ustar(:)/real(mcount)
  fg_ave(:) = fg_ave(:) + fg(:)/real(mcount)
  
  do k = 1,kl
    do iq = 1,imax
      templ(iq,k) = thetal(iq,k)/sigkap(k)
      pres(iq,k) = ps(iq)*sig(k)
      fice(iq,k) = min( qfg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-8), 1. )
    end do
  end do
  call getqsat(qsat,templ,pres,fice)
  
  ! Account for phase transistions
  do k = 1,kl
    do iq = 1,imax      
      ! Check for -ve values  
      qt = max( qfg(iq,k) + qlg(iq,k) + qvg(iq,k), 0. )
      qc = max( qfg(iq,k) + qlg(iq,k), 0. )
      
      qfg(iq,k) = max( qfg(iq,k), 0. )
      qlg(iq,k) = max( qlg(iq,k), 0. )
    
      ! account for saturation  
      lx = lv + lf*fice(iq,k)
      dqsdt = qsat(iq,k)*lx/(rv*templ(iq,k)**2)
      al = cp/(cp+lx*dqsdt)
      qc = max( al*(qt - qsat(iq,k)), qc )
      theta(iq,k) = thetal(iq,k) + sigkap(k)*(lv*qlg(iq,k)+ls*qfg(iq,k))/cp
      temp = theta(iq,k)/sigkap(k)
      !if ( temp>=tice ) then
        qfg(iq,k) = max( fice(iq,k)*qc, 0. )  
        qlg(iq,k) = max( qc - qfg(iq,k), 0. )
      !end if
      qvg(iq,k) = max( qt - qfg(iq,k) - qlg(iq,k), 0. )
      if ( qlg(iq,k)+qfg(iq,k)>1.E-12 ) then
        stratcloud(iq,k) = max( stratcloud(iq,k), 1.E-8 )
      end if
    end do
  end do

 
end do ! kcount loop

 
do k = 1,kl
  do iq = 1,imax
    theta(iq,k) = thetal(iq,k) + sigkap(k)*(lv*qlg(iq,k)+ls*qfg(iq,k))/cp
  end do
end do
  
! variance for wind gusts
! (from TAPM)
do iq = 1,imax
  wdash_sq = (1.2*ustar(iq))**2   ! Hurley
  !l_on_kz = cm34*tke(iq,1)**(3./2.)/eps(iq,1)/(vkar*zz(iq,1))
  !wdash_sq = ((2./3.)*tke(iq,1) + tke(iq,1)/(cs1*eps(iq,1))*( &
  !            (2.-cs2-cw2*l_on_kz)*pps(iq,1)                  &
  !          + (2.-cs3-cw3*l_on_kz)*ppb(iq,1)                  &
  !          - (2./3.)*eps(iq,1) ) )                           &
  !          / (1.+(cw1/cs1)*l_on_kz)
  ugs_var(iq) = tke(iq,1) - 0.5*wdash_sq      ! = (cm12-0.5*1.2)*ustar(iq)**2 + ce3*wstar(iq)**2
  !usg_var(iq) = (2.185*ustar(iq))**2         ! Schreur et al (2008)
end do

if ( mode==2 .or. mode==3 ) then
  fg(:) = fg_ave(:)
end if

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plume rise model
    
subroutine plumerise(mask,                                        &
                     zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,      &
                     zz,dz_hl,thetal,qvg,qlg,qfg,                 &
                     stratcloud,wt0,wq0,ps,ustar,                 &
                     sig,tke,ua,va,sigkap,imax,kl)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif

integer, intent(in) :: imax, kl
integer k, iq
real, dimension(imax,kl), intent(out) :: mflx, tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: qvg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: zz, thetal, ua, va 
real, dimension(imax,kl-1), intent(in) :: dz_hl
real, dimension(imax,kl), intent(in) :: tke
real, dimension(imax), intent(in) :: wt0, wq0, ps, ustar
real, dimension(kl), intent(in) :: sig, sigkap
real, dimension(imax), intent(inout) :: zi, wstar
real, dimension(imax,kl) ::  w2up, nn, cxup, rino
real, dimension(imax,kl) :: theta, thetav
real, dimension(imax) :: ent, qupsat
real, dimension(imax) :: wtv0, fice, pres, templ
real qtup, qxup, lx, dqsdt, al, qcup, thup, tvup
real vvk, as, bs, cs, xp, upf, dzht
real, parameter :: fac = 10. ! originally fac=100.
real, parameter :: ricr = 0.3
logical, dimension(imax), intent(in) :: mask

! Update theta and thetav
do k = 1,kl
  do iq = 1,imax
    theta(iq,k) = thetal(iq,k) + sigkap(k)*(lv*qlg(iq,k)+ls*qfg(iq,k))/cp
    thetav(iq,k) = theta(iq,k)*(1.+0.61*qvg(iq,k)-qlg(iq,k)-qfg(iq,k))
  end do
end do 

wtv0(:) = wt0(:) + theta(:,1)*0.61*wq0(:) ! thetav flux

! Initialise updraft
do k = 1,kl
  do iq = 1,imax
    mflx(iq,k) = 0.  
    w2up(iq,k) = 0.
    nn(iq,k) = 0.
    tlup(iq,k) = thetal(iq,k)
    qvup(iq,k) = qvg(iq,k)
    qlup(iq,k) = qlg(iq,k)
    qfup(iq,k) = qfg(iq,k)
    cfup(iq,k) = stratcloud(iq,k)
  end do
end do

! first level -----------------

! Entrainment rates
ent = entfn(zz(:,1), zi(:), imax)

! initial thermodynamic state
! split qtot into components (conservation of thetal and qtot is maintained)
do iq = 1,imax
  if ( mask(iq) ) then
    tlup(iq,1) = thetal(iq,1) + be*wt0(iq)/sqrt(max(tke(iq,1),1.5e-4))   ! Hurley 2007
    qvup(iq,1) = qvg(iq,1)    + be*wq0(iq)/sqrt(max(tke(iq,1),1.5e-4))   ! Hurley 2007
    qlup(iq,1) = qlg(iq,1)
    qfup(iq,1) = qfg(iq,1)
    cfup(iq,1) = stratcloud(iq,1)
  end if
  ! update updraft velocity and mass flux
  nn(iq,1) = grav*be*wtv0(iq)/(thetav(iq,1)*sqrt(max(tke(iq,1),1.5e-4))) ! Hurley 2007
  dzht = zz(iq,1)
  w2up(iq,1) = 2.*dzht*b2*nn(iq,1)/(1.+2.*dzht*b1*ent(iq))        ! Hurley 2007
  cxup(iq,1) = 0.
  rino(iq,1) = 0.
end do

! updraft with condensation
do k = 2,kl
  ! Entrainment rates
  ent = entfn(zz(:,k), zi(:), imax)
  templ = tlup(:,k)/sigkap(k)     ! templ,up
  pres = ps(:)*sig(k)
  fice = min( qfup(:,k)/max(qfup(:,k)+qlup(:,k),1.e-8), 1. )
  call getqsat(qupsat,templ,pres,fice)
  do iq = 1,imax
    dzht = dz_hl(iq,k-1)
    if ( w2up(iq,k-1)>0. .and. mask(iq) ) then
      ! entrain air into plume
      ! split qtot into components (conservation of qtot is maintained)
      tlup(iq,k) = (tlup(iq,k-1)+dzht*ent(iq)*thetal(iq,k))/(1.+dzht*ent(iq))
      qvup(iq,k) = (qvup(iq,k-1)+dzht*ent(iq)*qvg(iq,k)   )/(1.+dzht*ent(iq))
      qlup(iq,k) = (qlup(iq,k-1)+dzht*ent(iq)*qlg(iq,k)   )/(1.+dzht*ent(iq))
      qfup(iq,k) = (qfup(iq,k-1)+dzht*ent(iq)*qfg(iq,k)   )/(1.+dzht*ent(iq))
      cfup(iq,k) = (cfup(iq,k-1)+dzht*ent(iq)*stratcloud(iq,k))/(1.+dzht*ent(iq))
    end if
    ! calculate conserved variables
    qtup = qvup(iq,k) + qlup(iq,k) + qfup(iq,k)    ! qtot,up
    if ( qtup>qupsat(iq) .and. w2up(iq,k-1)>0. ) then
      qxup = qupsat(iq)
      cxup(iq,k) = 1.
    else
      qxup = qtup
      cxup(iq,k) = 0.
    end if
    lx = lv + lf*fice(iq)
    dqsdt = qupsat(iq)*lx/(rv*templ(iq)**2)
    al = cp/(cp+lx*dqsdt)
    qcup = max(al*(qtup-qxup), 0.)                        ! qcondensate,up after redistribution
    qcup = min(qcup, qcmf)                                ! limit condensation with simple autoconversion
    thup = tlup(iq,k) + sigkap(k)*qcup*lx/cp              ! theta,up after redistribution
    tvup = thup + theta(iq,k)*(0.61*qxup-qcup)            ! thetav,up after redistribution
    if ( w2up(iq,k-1)>0. ) then
      nn(iq,k) = grav*(tvup-thetav(iq,k))/thetav(iq,k)                        ! calculate buayancy
      w2up(iq,k) = (w2up(iq,k-1)+2.*dzht*b2*nn(iq,k))/(1.+2.*dzht*b1*ent(iq)) ! update updraft velocity
    else
      nn(iq,k) = 0.  
      w2up(iq,k) = 0.
    end if
    vvk = (ua(iq,k)-ua(iq,1))**2 + (va(iq,k)-va(iq,1))**2 + fac*ustar(iq)**2  
    rino(iq,k) = grav*(thetav(iq,k)-thetav(iq,1))*(zz(iq,k)-zz(iq,1))/max(thetav(iq,1)*vvk,1.e-20)
    ! test if maximum plume height is reached
    if ( w2up(iq,k)<=0. .and. w2up(iq,k-1)>0. .and. mask(iq) ) then ! unstable
      as = 2.*b2*(nn(iq,k)-nn(iq,k-1))/dzht
      bs = 2.*b2*nn(iq,k-1)
      cs = w2up(iq,k-1)
      xp = -2.*cs/(bs-sqrt(max(bs**2-4.*as*cs,0.)))
      xp = min(max(xp,0.),dzht)
      zi(iq) = xp + zz(iq,k-1)
    else if ( rino(iq,k)>ricr .and. rino(iq,k-1)<=ricr .and. .not.mask(iq) ) then ! stable
      xp = (ricr-rino(iq,k-1))/(rino(iq,k)-rino(iq,k-1))
      xp = min( max(xp, 0.), 1.)
      zi(iq) = zz(iq,k-1) + xp*(zz(iq,k)-zz(iq,k-1))
    end if
  end do
end do

! update wstar with new zi  
wstar = (grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)
          
! update mass flux
mflx(:,1) = m0*sqrt(max(w2up(:,1), 0.))
do k = 2,kl
  do iq = 1,imax
    dzht = dz_hl(iq,k-1)
    mflx(iq,k) = (1.-cxup(iq,k))*m0*sqrt(max(w2up(iq,k), 0.))         &
              + cxup(iq,k)*mflx(iq,k-1)/(1.+dzht*(dtrc0-entc0))
    upf = mflx(iq,k-1)/sqrt(max(w2up(iq,k-1), 1.e-8))
    mflx(iq,k) = min( mflx(iq,k), upf*sqrt(max(w2up(iq,k), 0.)) )
  end do
end do

return
end subroutine plumerise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update TKE and EPS

subroutine update_tkeeps(tke,eps,ppb,pps,ppt,qvg,qlg,qfg,thetal,                 &
                         zzh,zz,idzp,idzm,wstar,zi,thetal_ema,qv_ema,            &
                         ql_ema,qf_ema,cf_ema,tke_ema,ddts,qgmin,cm34,sigkap,    &
                         imax,kl)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif

implicit none

integer, intent(in) :: imax, kl
integer k, iq
real, dimension(imax,kl), intent(inout) :: tke, eps
real, dimension(imax,kl), intent(inout) :: ppb, pps, ppt
real, dimension(imax,kl), intent(in) :: qvg, qlg, qfg
real, dimension(imax,kl), intent(in) :: thetal
real, dimension(imax,2:kl), intent(in) :: idzm
real, dimension(imax,kl), intent(in) :: zz, zzh
real, dimension(imax,kl), intent(inout) :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(imax,kl), intent(inout) :: tke_ema
real, dimension(imax,kl-1), intent(in) :: idzp
real, dimension(imax,kl) :: qsatc, qgnc, thetalhl, quhl, qshl, qlhl, qfhl
real, dimension(imax,kl) :: qthl, thetavhl, kmo, qtot
real, dimension(imax,kl) :: dz_fl, theta, thetav, km
real, dimension(imax,kl) :: rr, ff, gg
real, dimension(imax,kl) :: theta_ema, thetav_ema, qtot_ema
real, dimension(imax,kl-1) :: fzzh, dz_hl
real, dimension(imax,2:kl) :: qq
real, dimension(imax,3:kl-1) :: aa
real, dimension(imax,2:kl-2) :: cc
real, dimension(imax,2:kl-1) :: bb, dd
real, dimension(imax), intent(in) :: wstar, zi
real, dimension(kl), intent(in) :: sigkap
real, intent(in) :: qgmin, ddts, cm34
real tbb, tcc, tff, tqq, thetac, tempc, tempt, tempv, rvar, bvf, mc, dc, fc
logical, dimension(imax,kl) :: lta

! set top boundary condition for TKE-eps source terms
pps(:,kl) = 0.
ppb(:,kl) = 0.
ppt(:,1)  = 0.
ppt(:,kl) = 0.

! Fraction for interpolation from full levels to half levels
fzzh(:,1:kl-1) = (zzh(:,1:kl-1)-zz(:,1:kl-1))/(zz(:,2:kl)-zz(:,1:kl-1))

! Calculate dz at half levels
dz_hl(:,1:kl-1) = max( zz(:,2:kl) - zz(:,1:kl-1), 1. )

! Calculate dz at full levels
dz_fl(:,1)    = zzh(:,1)
dz_fl(:,2:kl) = zzh(:,2:kl) - zzh(:,1:kl-1)
dz_fl = max( dz_fl, 1. )

! Update qtot, theta and thetav
qtot(:,:) = qvg(:,:) + qlg(:,:) + qfg(:,:)
do k = 1,kl
  theta(:,k) = thetal(:,k) + sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
end do
thetav(:,:) = theta(:,:)*(1.+0.61*qvg(:,:)-qlg(:,:)-qfg(:,:))
qtot_ema(:,:) = qv_ema(:,:) + ql_ema(:,:) + qf_ema(:,:)  
do k = 1,kl
  theta_ema(:,k) = thetal_ema(:,k) + sigkap(k)*(lv*ql_ema(:,k)+ls*qf_ema(:,k))/cp
end do
thetav_ema(:,:) = theta_ema(:,:)*(1.+0.61*qv_ema(:,:)-ql_ema(:,:)-qf_ema(:,:))

! interpolate diffusion coeffs to half levels
km = cm0*tke(:,:)**2/eps(:,:)
kmo(:,1:kl-1)=km(:,1:kl-1)+fzzh(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))

! top boundary condition to avoid unphysical behaviour at the top of the model
tke(:,kl) = mintke
eps(:,kl) = mineps

! Calculate buoyancy term
select case(buoymeth)
  case(0) ! Blend staturated and unsaturated terms - saturated method from Durran and Klemp JAS 1982 (see also WRF)
    qsatc = qv_ema(:,:)                       ! assume qvg is saturated inside cloud
    ff = qf_ema(:,:)/max(cf_ema,1.E-8)        ! inside cloud value assuming max overlap
    gg = ql_ema(:,:)/max(cf_ema,1.E-8)        ! inside cloud value assuming max overlap
    do k = 1,kl
      do iq = 1,imax
          tbb = max(1.-cf_ema(iq,k),1.E-8)
          qgnc(iq,k) = (qv_ema(iq,k)-(1.-tbb)*qsatc(iq,k))/tbb    ! outside cloud value
          qgnc(iq,k) = min(max(qgnc(iq,k),qgmin),qsatc(iq,k))
        end do
      end do
      thetalhl(:,1:kl-1)=thetal_ema(:,1:kl-1)+fzzh(:,1:kl-1)*(thetal_ema(:,2:kl)-thetal_ema(:,1:kl-1)) ! outside cloud value
      quhl(:,1:kl-1)=qgnc(:,1:kl-1)+fzzh(:,1:kl-1)*(qgnc(:,2:kl)-qgnc(:,1:kl-1))                       ! outside cloud value
      qshl(:,1:kl-1)=qsatc(:,1:kl-1)+fzzh(:,1:kl-1)*(qsatc(:,2:kl)-qsatc(:,1:kl-1))                    ! inside cloud value
      qlhl(:,1:kl-1)=gg(:,1:kl-1)+fzzh(:,1:kl-1)*(gg(:,2:kl)-gg(:,1:kl-1))                             ! inside cloud value
      qfhl(:,1:kl-1)=ff(:,1:kl-1)+fzzh(:,1:kl-1)*(ff(:,2:kl)-ff(:,1:kl-1))                             ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl) = cf_ema(:,2:kl)<=1.E-6
      do k = 2,kl-1
        where( lta(:,k) .and. .not.lta(:,k+1) )
          qlhl(:,k) = gg(:,k+1)
          qfhl(:,k) = ff(:,k+1)
        elsewhere ( .not.lta(:,k) .and. lta(:,k+1) )
          qlhl(:,k) = gg(:,k)
          qfhl(:,k) = ff(:,k)
        end where
      end do
      do k = 2,kl-1
        do iq = 1,imax
          ! saturated
          thetac = thetal_ema(iq,k)+sigkap(k)*(lv*gg(iq,k)+ls*ff(iq,k))/cp   ! inside cloud value
          tempc = thetac/sigkap(k)                                           ! inside cloud value          
          tqq = (1.+lv*qsatc(iq,k)/(rd*tempc))/(1.+lv*lv*qsatc(iq,k)/(cp*rv*tempc**2))
          tbb = -grav*km(iq,k)*(tqq*((thetalhl(iq,k)-thetalhl(iq,k-1)+sigkap(k)/cp*(lv*(qlhl(iq,k)-qlhl(iq,k-1))  &
                +ls*(qfhl(iq,k)-qfhl(iq,k-1))))/thetac+lv/cp*(qshl(iq,k)-qshl(iq,k-1))/tempc)                     &
                -qshl(iq,k)-qlhl(iq,k)-qfhl(iq,k)+qshl(iq,k-1)+qlhl(iq,k-1)+qfhl(iq,k-1))/dz_fl(iq,k)
          ! unsaturated
          tcc = -grav*km(iq,k)*(thetalhl(iq,k)-thetalhl(iq,k-1)+thetal_ema(iq,k)*0.61*(quhl(iq,k)-quhl(iq,k-1)))  &
                           /(thetal_ema(iq,k)*dz_fl(iq,k))
          ppb(iq,k) = (1.-cf_ema(iq,k))*tcc+cf_ema(iq,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)
        end do
      end do
      ! saturated
      do iq = 1,imax
        thetac = thetal_ema(iq,1)+sigkap(1)*(lv*gg(iq,1)+ls*ff(iq,1))/cp        ! inside cloud value
        tempc = thetac/sigkap(1)                                                ! inside cloud value          
        tqq = (1.+lv*qsatc(iq,1)/(rd*tempc))/(1.+lv*lv*qsatc(iq,1)/(cp*rv*tempc*tempc))
        tbb = -grav*km(iq,1)*(tqq*((thetalhl(iq,1)-thetal_ema(iq,1)+sigkap(1)/cp*(lv*(qlhl(iq,1)-ql_ema(iq,1))  &
              +ls*(qfhl(iq,1)-qf_ema(iq,1))))/thetac+lv/cp*(qshl(iq,1)-qsatc(iq,1))/tempc)                      &
              -qshl(iq,1)-qlhl(iq,1)-qfhl(iq,1)+qsatc(iq,1)+ql_ema(iq,1)+qf_ema(iq,1))/(zzh(iq,1)-zz(iq,1))
        ! unsaturated
        tcc = -grav*km(iq,1)*(thetalhl(iq,1)-thetal_ema(iq,1)+thetal_ema(iq,1)*0.61*(quhl(iq,1)-qgnc(iq,1)))    &
                       /(thetal_ema(iq,1)*(zzh(iq,1)-zz(iq,1)))
        ppb(iq,1) = (1.-cf_ema(iq,1))*tcc+cf_ema(iq,1)*tbb ! cloud fraction weighted (e.g., Smith 1990)
      end do

      
    case(1) ! Marquet and Geleyn QJRMS (2012) for partially saturated
      thetalhl(:,1:kl-1)=km(:,1:kl-1)+fzzh(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1)) 
      qthl(:,1:kl-1)=qtot(:,1:kl-1)+fzzh(:,1:kl-1)*(qtot(:,2:kl)-qtot(:,1:kl-1))
      do k = 2,kl-1
        do iq = 1,imax
          tempt = theta_ema(iq,k)/sigkap(k)
          tempv = thetav_ema(iq,k)/sigkap(k)
          rvar = rd*tempv/tempt ! rvar = qd*rd+qv*rv
          fc = (1.-cf_ema(iq,k))+cf_ema(iq,k)*(lv*rvar/(cp*rv*tempt))
          dc = (1.+0.61*qv_ema(iq,k))*lv*qv_ema(iq,k)/(rd*tempv)
          mc = (1.+dc)/(1.+lv*ql_ema(iq,k)/(cp*tempt)+dc*fc)
          bvf = grav*mc*(thetalhl(iq,k)-thetalhl(iq,k-1))/(thetal_ema(iq,k)*dz_fl(iq,k))           &
               +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(iq,k)-qthl(iq,k-1))/dz_fl(iq,k)
          ppb(iq,k) = -km(iq,k)*bvf
        end do
      end do
      do iq = 1,imax
        tempt = theta_ema(iq,1)/sigkap(1)
        tempv = thetav_ema(iq,1)/sigkap(1)
        rvar = rd*tempv/tempt ! rvar = qd*rd+qv*rv
        fc = (1.-cf_ema(iq,1))+cf_ema(iq,1)*(lv*rvar/(cp*rv*tempt))
        dc = (1.+0.61*qv_ema(iq,1))*lv*qv_ema(iq,1)/(rd*tempv)
        mc = (1.+dc)/(1.+lv*ql_ema(iq,1)/(cp*tempt)+dc*fc)
        bvf = grav*mc*(thetalhl(iq,1)-thetal_ema(iq,1))/(thetal_ema(iq,1)*(zzh(iq,1)-zz(iq,1)))         &
             +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(iq,1)-qtot_ema(iq,1))/(zzh(iq,1)-zz(iq,1))
        ppb(iq,1) = -km(iq,1)*bvf
      end do
      
    case(2) ! dry convection from Hurley 2007
      thetavhl(:,1:kl-1)=thetav_ema(:,1:kl-1)+fzzh(:,1:kl-1)*(thetav_ema(:,2:kl)-thetav_ema(:,1:kl-1))
      do k = 2,kl-1
        do iq = 1,imax
          tcc = -grav*km(iq,k)*(thetavhl(iq,k)-thetavhl(iq,k-1))/(thetav_ema(iq,k)*dz_fl(iq,k))
          ppb(iq,k) = tcc
        end do
      end do
      do iq = 1,imax
        tcc = -grav*km(iq,1)*(thetavhl(iq,1)-thetav(iq,1))/(thetav_ema(iq,1)*(zzh(iq,1)-zz(iq,1)))
        ppb(iq,1) = tcc 
      end do

  end select
    
  ! Calculate transport source term on full levels
  ppt(:,2:kl-1)= kmo(:,2:kl-1)*idzp(:,2:kl-1)*(tke_ema(:,3:kl)-tke_ema(:,2:kl-1))/dz_hl(:,2:kl-1)  &
                -kmo(:,1:kl-2)*idzm(:,2:kl-1)*(tke_ema(:,2:kl-1)-tke_ema(:,1:kl-2))/dz_hl(:,1:kl-2)
  
  ! Pre-calculate eddy diffusivity mixing terms
  ! -ve because gradient is calculated at t+1
  qq(:,2:kl-1)=-ddts*idzm(:,2:kl-1)/dz_hl(:,1:kl-2)
  rr(:,2:kl-1)=-ddts*idzp(:,2:kl-1)/dz_hl(:,2:kl-1)
  
  ! eps vertical mixing
  aa(:,3:kl-1)=ce0*kmo(:,2:kl-2)*qq(:,3:kl-1)
  cc(:,2:kl-2)=ce0*kmo(:,2:kl-2)*rr(:,2:kl-2)
  ! follow Hurley 2007 to make scheme more numerically stable
  bb(:,2:kl-1)=1.-ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)-ce0*kmo(:,2:kl-1)*rr(:,2:kl-1) &
               +ddts*ce2*eps(:,2:kl-1)/tke(:,2:kl-1)
  dd(:,2:kl-1)=eps(:,2:kl-1)+ddts*eps(:,2:kl-1)/tke(:,2:kl-1)                    &
              *ce1*(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
  dd(:,2)     =dd(:,2)   -ce0*kmo(:,1)*qq(:,2)*eps(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-ce0*kmo(:,kl-1)*rr(:,kl-1)*mineps
  call thomas(eps(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax,kl-2)
  
  ! TKE vertical mixing
  aa(:,3:kl-1)=kmo(:,2:kl-2)*qq(:,3:kl-1)
  cc(:,2:kl-2)=kmo(:,2:kl-2)*rr(:,2:kl-2)
  bb(:,2:kl-1)=1.-kmo(:,1:kl-2)*qq(:,2:kl-1)-kmo(:,2:kl-1)*rr(:,2:kl-1)
  dd(:,2:kl-1)=tke(:,2:kl-1)+ddts*(pps(:,2:kl-1)+ppb(:,2:kl-1)-eps(:,2:kl-1))
  dd(:,2)     =dd(:,2)   -kmo(:,1)*qq(:,2)*tke(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-kmo(:,kl-1)*rr(:,kl-1)*mintke
  call thomas(tke(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax,kl-2)

  ! limit decay of TKE and EPS with coupling to mass flux term
  if ( tkemeth==1 ) then
    do k = 2,kl-1
      do iq = 1,imax
        if ( wstar(iq)>0.5 .and. zz(iq,k)>0.5*zi(iq) .and. zz(iq,k)<0.95*zi(iq) ) then
          tbb = max(1.-0.05*dz_hl(iq,k-1)/250.,0.)          
          tke(iq,k) = max( tke(iq,k), tbb*tke(iq,k-1) )
          eps(iq,k) = max( eps(iq,k), tbb*eps(iq,k-1) )
        end if
      end do
    end do
  end if
  
  ! apply limits and corrections to tke and eps terms
  do k = 2,kl-1
    do iq = 1,imax
      tke(iq,k) = max(tke(iq,k),mintke)
      tff = cm34*tke(iq,k)*sqrt(tke(iq,k))
      eps(iq,k) = min(eps(iq,k),tff/minl)
      eps(iq,k) = max(eps(iq,k),tff/maxl,mineps)
    end do
  end do

return
end subroutine update_tkeeps
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update atmosphere only

subroutine update_atmosphere(thetal,qvg,qlg,qfg,stratcloud,ua,va, &
                             tlup,qvup,qlup,qfup,cfup,fg,eg,rhos, &
                             ustar,cduv,                          &
                             tke,eps,mflx,fzzh,idzp,idzm,dz_hl,   &
                             rhoa1,dz_fl1,sigkap,                 &
#ifdef scm
                             wthflux,wqvflux,uwflux,vwflux,       &
#endif
                             ddts,imax,kl)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif            

implicit none

integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(inout) :: thetal, qvg, qlg, qfg, stratcloud, ua, va
real, dimension(imax,kl), intent(in) :: tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: mflx, tke, eps
real, dimension(imax,2:kl), intent(in) :: idzm
real, dimension(imax,kl-1), intent(in) :: fzzh, idzp, dz_hl
real, dimension(imax,kl) :: km, kmo
real, dimension(imax,kl) :: rr, bb, cc
real, dimension(imax,kl,5) :: dd, ans
real, dimension(imax,2:kl) :: qq, aa
real, dimension(imax), intent(inout) :: fg, eg, ustar
real, dimension(imax), intent(in) :: rhos, rhoa1, dz_fl1, cduv
real, dimension(imax) :: wt0, wq0
real, dimension(kl), intent(in) :: sigkap
real, intent(in) :: ddts

#ifdef scm
integer k
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl) :: wthlflux, wqlflux, wqfflux
#endif

! Calculate surface fluxes
wt0 = fg/(rhos*cp)  ! theta flux
wq0 = eg/(rhos*lv)  ! qtot flux

! estimate eddy diffusivity mixing coeff
km(:,:) = cm0*tke(:,:)**2/eps(:,:)
kmo(:,1:kl-1)=km(:,1:kl-1)+fzzh(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))
  
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
bb(:,2:kl-1)=1.-qq(:,2:kl-1)-rr(:,2:kl-1)+ddts*(mflx(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)         &
                                               -mflx(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1))
aa(:,kl)=qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
bb(:,kl)=1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)


! thetal vertical mixing
dd(:,1,1)=thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                               &
                           +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                   &
                     +ddts*rhos*wt0/(rhoa1(:)*dz_fl1(:))
dd(:,2:kl-1,1)=thetal(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*tlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1) &
                                     +mflx(:,2:kl-1)*tlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)      &
                                     -mflx(:,2:kl-1)*tlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1) &
                                     -mflx(:,3:kl)*tlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl,1)=thetal(:,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                   &
                             +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
! qv (part of qtot) vertical mixing
dd(:,1,2)=qvg(:,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                  &
                        +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                                      &
                        +ddts*rhos*wq0/(rhoa1(:)*dz_fl1(:))
dd(:,2:kl-1,2)=qvg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qvup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)    &
                                  +mflx(:,2:kl-1)*qvup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)         &
                                  -mflx(:,2:kl-1)*qvup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)    &
                                  -mflx(:,3:kl)*qvup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl,2)=qvg(:,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                      &
                          +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
! ql (part of qtot) vertical mixing
dd(:,1,3)=qlg(:,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                  &
                        +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
dd(:,2:kl-1,3)=qlg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)    &
                                  +mflx(:,2:kl-1)*qlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)         &
                                  -mflx(:,2:kl-1)*qlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)    &
                                  -mflx(:,3:kl)*qlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl,3)=qlg(:,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                      &
                          +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
! qf (part of qtot) vertical mixing
dd(:,1,4)=qfg(:,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                 &
                        +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
dd(:,2:kl-1,4)=qfg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)   &
                                  +mflx(:,2:kl-1)*qfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                  -mflx(:,2:kl-1)*qfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                  -mflx(:,3:kl)*qfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl,4)=qfg(:,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                     &
                          +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
! cloud fraction vertical mixing
dd(:,1,5) = stratcloud(:,1) - ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                 &
                                   +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))
dd(:,2:kl-1,5) = stratcloud(:,2:kl-1) + ddts*(mflx(:,1:kl-2)*cfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)   &
                                             +mflx(:,2:kl-1)*cfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                             -mflx(:,2:kl-1)*cfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                             -mflx(:,3:kl)*cfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl,5) = stratcloud(:,kl) + ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                     &
                                     +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))

call thomas(ans(:,:,1:5),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,:,1:5),imax,kl,5)
thetal(:,:)     = ans(:,:,1)
qvg(:,:)        = ans(:,:,2)
qlg(:,:)        = ans(:,:,3)
qfg(:,:)        = ans(:,:,4)
stratcloud(:,:) = min( max( ans(:,:,5), 0. ), 1. )
    
#ifdef scm  
wthlflux(:,1)=wt0(:)
wthlflux(:,2:kl)=-kmo(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)                    &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))               &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
wqvflux(:,1)=wq0(:)
wqvflux(:,2:kl)=-kmo(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)                           &
                +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                   &
                +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
wqlflux(:,1)=0.
wqlflux(:,2:kl)=-kmo(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)                                &
                +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                        &
                +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
wqfflux(:,1)=0.
wqfflux(:,2:kl)=-kmo(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)                                &
                +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                        &
                +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif

  
! momentum vertical mixing
aa(:,2:kl)   = qq(:,2:kl)
cc(:,1:kl-1) = rr(:,1:kl-1)
bb(:,1) = 1. - cc(:,1) + ddts*rhos*cduv/(rhoa1(:)*dz_fl1(:)) ! implicit  
bb(:,2:kl-1) = 1. - aa(:,2:kl-1) - cc(:,2:kl-1)
bb(:,kl) = 1. - aa(:,kl)

dd(:,1:kl,1) = ua(:,1:kl)
dd(:,1:kl,2) = va(:,1:kl)
call thomas(ans(:,:,1:2),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,:,1:2),imax,kl,2)
ua(:,:) = ans(:,:,1)
va(:,:) = ans(:,:,2)
  
  
! update surface momentum flux
ustar = sqrt(cduv*sqrt(ua(:,1)**2+va(:,1)**2))


#ifdef scm
uwflux(:,1) = 0.
uwflux(:,2:kl) = -kmo(:,1:kl-1)*(ua(:,2:kl)-ua(:,1:kl-1))/dz_hl(:,1:kl-1)
vwflux(:,1) = 0.
vwflux(:,2:kl) = -kmo(:,1:kl-1)*(va(:,2:kl)-va(:,1:kl-1))/dz_hl(:,1:kl-1)
wthflux(:,1) = wthlflux(:,1) - ((1.-fzzh(:,1)+sigkap(1)*fzzh(:,1)) &
                             *(lv*wqlflux(:,1)+ls*wqfflux(:,1)))
do k = 2,kl
  wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                               *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
end do
#endif

return
end subroutine update_atmosphere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update coupled

subroutine update_coupled(thetal,qvg,qlg,qfg,stratcloud,ua,va,    &
                          tlup,qvup,qlup,qfup,cfup,fg,eg,rhos,    &
                          ustar,cduv,                             &    
                          tke,eps,mflx,fzzh,idzp,idzm,dz_hl,      &
                          rhoa1,dz_fl1,                           &
                          rad_o,w_t,w_s,w_u,w_v,                  &
                          deptho_fl,deptho_hl,                    &
                          cd_water,cdh_water,cdbot_water,         &
                          wt0rad_o,wt0melt_o,                     &
                          wt0eg_o,icefg_a,wt0fb_o,ws0_o,          &
                          ws0subsurf_o,i_u,i_v,imass,fracice,     &
                          cd_ice,cdbot_ice,ibot,land,sigkap,      &
                          km_o,ks_o,rhs_o,                        &
#ifdef scm
                          wthflux,wqvflux,uwflux,vwflux,          &
#endif
                          ddts,wrtemp,wrtrho,cp0,minsfc,          &
                          imax,kl,wlev)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif
         
implicit none

integer, intent(in) :: imax, kl, wlev
integer, dimension(:), intent(in) :: ibot
integer k, kn, iq
real, dimension(:,:), intent(inout) :: thetal, qvg, qlg, qfg, stratcloud, ua, va
real, dimension(:,:), intent(in) :: tlup, qvup, qlup, qfup, cfup
real, dimension(:,:), intent(in) :: mflx, tke, eps
real, dimension(:,2:), intent(in) :: idzm
real, dimension(:,:), intent(in) :: fzzh, idzp, dz_hl
real, dimension(imax,kl) :: km_a, kmo_a
real, dimension(imax,kl-1) :: rr
real, dimension(imax,2:kl) :: qq
real, dimension(imax,kl) :: bb_a
real, dimension(imax,kl-1) :: cc_a
real, dimension(imax,kl,4) :: dd_a, tt_a
real, dimension(imax,2:kl) :: aa_a
real, dimension(:), intent(in) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o
real, dimension(:,:), intent(in) :: rad_o
real, dimension(:), intent(in) :: wt0eg_o
real, dimension(:,:), intent(inout) :: w_t, w_s, w_u, w_v
real, dimension(imax,wlev), intent(in) :: km_o, ks_o, rhs_o
real, dimension(imax,wlev) :: bb_o
real, dimension(imax,wlev-1) :: cc_o
real, dimension(imax,2:wlev) :: aa_o
real, dimension(imax,wlev,2) :: dd_o, w_tt
real, dimension(imax,wlev), intent(in) :: deptho_fl
real, dimension(imax,wlev+1), intent(in) :: deptho_hl
real, dimension(imax,wlev) :: deptho_dz
real, dimension(imax,2:wlev) :: deptho_dz_hl
real, dimension(imax) :: bb_i
real, dimension(imax,2) :: dd_i
real, dimension(:), intent(inout) :: fg, eg, ustar
real, dimension(:), intent(in) :: rhos, rhoa1, dz_fl1, cduv
real, dimension(:), intent(in) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o
real, dimension(:), intent(in) :: fracice, imass, cd_ice, cdbot_ice
real, dimension(imax) :: wq0_a, wt0_o
real, dimension(:), intent(inout) :: i_u, i_v
real, dimension(imax) :: wt0_a, t1, t3, t4
real, dimension(imax) :: fulldeptho, targetdeptho
real, dimension(imax) :: f_ao, f_oa, f_ai, f_ia, f_oi, f_io
real, dimension(:), intent(in) :: sigkap
real, intent(in) :: ddts, wrtemp, wrtrho, cp0, minsfc
real deltazo
logical, dimension(:), intent(in) :: land

#ifdef scm
real, dimension(:,:), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl) :: wthlflux, wqlflux, wqfflux
#endif


deptho_dz_hl(:,2:wlev) = deptho_fl(:,2:wlev) - deptho_fl(:,1:wlev-1)
deptho_dz(:,1:wlev) = deptho_hl(:,2:wlev+1) - deptho_hl(:,1:wlev)

! Calculate surface fluxes
wq0_a(:) = eg(:)/(rhos(:)*lv)  ! qtot flux
wt0_a(:) = fg(:)/(rhos(:)*cp)  ! theta flux

t1(:) = rhos(:)*cdh_water(:)*cp
where ( .not.land(:) )
  wt0_o(:) = wt0rad_o(:) + wt0melt_o(:) + wt0eg_o(:) + wt0fb_o(:) &
         + (1.-fracice(:))*t1*(w_t(:,1)+wrtemp-thetal(:,1)-sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp)/(wrtrho*cp0)
elsewhere
  wt0_o(:) = 0.    
end where 


! estimate eddy diffusivity mixing coeff for atmosphere
km_a(:,:) = cm0*tke(:,:)**2/eps(:,:)
! interpolate diffusion coeffs to half levels
kmo_a(:,1:kl-1)=km_a(:,1:kl-1)+fzzh(:,1:kl-1)*(km_a(:,2:kl)-km_a(:,1:kl-1))

! Pre-calculate eddy diffiusivity mixing terms using updated kmo values
! -ve because gradient is calculated at t+1
qq(:,2:kl)  =-ddts*kmo_a(:,1:kl-1)*idzm(:,2:kl)/dz_hl(:,1:kl-1)
rr(:,1:kl-1)=-ddts*kmo_a(:,1:kl-1)*idzp(:,1:kl-1)/dz_hl(:,1:kl-1)


!------------------------------------------------------------------------------

! thetal, thetao - coupled ----

t1(:) = rhos(:)*cdh_water(:)*cp

! k=1 is top of atmosphere and k=kl is bottom of atmosphere
cc_a(:,1) = qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
bb_a(:,1) = 1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)
aa_a(:,2:kl-1) = rr(:,kl-1:2:-1)-ddts*mflx(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1)
cc_a(:,2:kl-1) = qq(:,kl-1:2:-1)+ddts*mflx(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1)
bb_a(:,2:kl-1) = 1.-qq(:,kl-1:2:-1)-rr(:,kl-1:2:-1)+ddts*(mflx(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)               &
                                                         -mflx(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1))
aa_a(:,kl) = rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)
bb_a(:,kl) = 1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)

where ( deptho_dz(:,2)*deptho_dz(:,1)>1.e-4 )
  cc_o(:,1) = -ddts*ks_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
elsewhere
  cc_o(:,1) = 0.
end where
bb_o(:,1) = 1. - cc_o(:,1)
where ( deptho_dz(:,2:wlev-1)*deptho_dz(:,2:wlev-1)>1.e-4 )
  aa_o(:,2:wlev-1) = -ddts*ks_o(:,2:wlev-1)/(deptho_dz_hl(:,2:wlev-1)*deptho_dz(:,2:wlev-1))
elsewhere
  aa_o(:,2:wlev-1) = 0.
end where
where ( deptho_dz(:,3:wlev)*deptho_dz(:,2:wlev-1)>1.e-4 )
  cc_o(:,2:wlev-1) = -ddts*ks_o(:,3:wlev)/(deptho_dz_hl(:,3:wlev)*deptho_dz(:,2:wlev-1))
elsewhere
  cc_o(:,2:wlev-1) = 0.
end where
bb_o(:,2:wlev-1) = 1. - aa_o(:,2:wlev-1) - cc_o(:,2:wlev-1)
where ( deptho_dz(:,wlev)*deptho_dz(:,wlev)>1.e-4 )
  aa_o(:,wlev) = -ddts*ks_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
elsewhere
  aa_o(:,wlev) = 0.
end where
bb_o(:,wlev) = 1. - aa_o(:,wlev)


! k=1 is top of atmosphere and k=kl is bottom of atmosphere
dd_a(:,1,1) = thetal(:,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)   &
                                +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
dd_a(:,2:kl-1,1) = thetal(:,kl-1:2:-1)+ddts*(mflx(:,kl-2:1:-1)*tlup(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1) &
                                            +mflx(:,kl-1:2:-1)*tlup(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)      &
                                            -mflx(:,kl-1:2:-1)*tlup(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1) &
                                            -mflx(:,kl:3:-1)*tlup(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1))
dd_a(:,kl,1) = thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)             &
                                +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))

where ( land(:) )
  dd_a(:,kl,1) = dd_a(:,kl,1)+ddts*rhos(:)*wt0_a(:)/(rhoa1(:)*dz_fl1(:))
elsewhere
  bb_a(:,kl) = bb_a(:,kl)+ddts*t1(:)*(1.-fracice(:))/cp/(rhoa1(:)*dz_fl1(:)) 
  dd_a(:,kl,1) = dd_a(:,kl,1)-ddts*t1(:)*(1.-fracice(:))/cp*sigkap(1)                                     & 
                *(lv*qlg(:,1)+ls*qfg(:,1))/cp/(rhoa1(:)*dz_fl1(:))
  dd_a(:,kl,1) = dd_a(:,kl,1)+ddts*t1(:)*(1.-fracice(:))*wrtemp/cp/(rhoa1(:)*dz_fl1(:))
  dd_a(:,kl,1) = dd_a(:,kl,1)+ddts*(fracice(:)*icefg_a(:)/cp)/(rhoa1(:)*dz_fl1(:))
end where  

dd_o(:,1,1) = w_t(:,1) + ddts*rhs_o(:,1)*wt0_o(:)
where ( deptho_dz(:,1)>1.e-4 )
  dd_o(:,1,1) = dd_o(:,1,1) - ddts*rad_o(:,1)/deptho_dz(:,1)
end where
do k = 2,wlev-1
  dd_o(:,k,1) = w_t(:,k) + ddts*rhs_o(:,k)*wt0_o(:)
end do
where ( deptho_dz(:,2:wlev-1)>1.e-4 )
  dd_o(:,2:wlev-1,1) = dd_o(:,2:wlev-1,1) - ddts*rad_o(:,2:wlev-1)/deptho_dz(:,2:wlev-1)
end where
dd_o(:,wlev,1) = w_t(:,wlev) + ddts*rhs_o(:,wlev)*wt0_o(:)
where ( deptho_dz(:,wlev)>1.e-4 )
  dd_o(:,wlev,1) = dd_o(:,wlev,1) - ddts*rad_o(:,wlev)/deptho_dz(:,wlev)
end where
where ( .not.land(:) )
  bb_o(:,1) = bb_o(:,1) + ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0*deptho_dz(:,1))
  dd_o(:,1,1) = dd_o(:,1,1) - ddts*t1(:)*(1.-fracice(:))*wrtemp/(wrtrho*cp0*deptho_dz(:,1))
  dd_o(:,1,1) = dd_o(:,1,1) - ddts*(wt0rad_o(:)+wt0melt_o(:)+wt0eg_o(:))/deptho_dz(:,1)
  dd_o(:,1,1) = dd_o(:,1,1) + ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0)*sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1)) &
                              /cp/deptho_dz(:,1)
  dd_o(:,1,1) = dd_o(:,1,1) - ddts*wt0fb_o(:)/deptho_dz(:,1)
end where
where ( .not.land(:) )
  f_ao(:) = -ddts*t1(:)*(1.-fracice(:))/cp/(rhoa1(:)*dz_fl1(:))
  f_oa(:) = -ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0*deptho_dz(:,1))
elsewhere
  f_ao = 0.
  f_oa = 0.  
end where  

call solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a(:,:,1),tt_a(:,:,1),  &
                            aa_o,bb_o,cc_o,dd_o(:,:,1),w_t,          &
                            f_ao,f_oa,imax,kl,wlev) 

thetal(:,1:kl) = tt_a(:,kl:1:-1,1)

! Derive sensible heat flux for water gridpoints
!wt0_o(:) = 0.
where ( .not.land(:) )
  fg(:) = (1.-fracice(:))*t1(:)*(w_t(:,1)+wrtemp-thetal(:,1) &
              -sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp) &
              + fracice(:)*icefg_a(:)
  !wt0_o(:) = wt0rad_o(:) + wt0melt_o(:) + wt0eg_o(:) + wt0fb_o(:) &
  !       + (1.-fracice(:))*t1*(w_t(:,1)+wrtemp-thetal(:,1)-sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp)/(wrtrho*cp0)
end where 

#ifdef scm
wthlflux(:,1)=fg(:)/(rhos(:)*cp)  ! theta flux
wthlflux(:,2:kl)=-kmo_a(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)    &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1)) &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif


!------------------------------------------------------------------------------

! set-up decoupled atmosphere matrices

cc_a(:,1) = qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
bb_a(:,1) = 1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)
aa_a(:,2:kl-1) = rr(:,kl-1:2:-1)-ddts*mflx(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1) 
cc_a(:,2:kl-1) = qq(:,kl-1:2:-1)+ddts*mflx(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1)
bb_a(:,2:kl-1) = 1.-qq(:,kl-1:2:-1)-rr(:,kl-1:2:-1)+ddts*(mflx(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)       &
                                                         -mflx(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1))
aa_a(:,kl) = rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)
bb_a(:,kl) = 1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)


! Note that vertical interpolation is linear so that qtot can be
! decomposed into qv, ql and qf.

! qv (part of qtot) - atmosphere
dd_a(:,1,1)=qvg(:,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                           +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
dd_a(:,2:kl-1,1)=qvg(:,kl-1:2:-1)+ddts*(mflx(:,kl-2:1:-1)*qvup(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1)   &
                                       +mflx(:,kl-1:2:-1)*qvup(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)        &
                                       -mflx(:,kl-1:2:-1)*qvup(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1)   &
                                       -mflx(:,kl:3:-1)*qvup(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1))
dd_a(:,kl,1)=qvg(:,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                           +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                   &
                     +ddts*rhos(:)*wq0_a(:)/(rhoa1(:)*dz_fl1(:))

! ql (part of qtot) - atmosphere
dd_a(:,1,2)=qlg(:,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                           +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
dd_a(:,2:kl-1,2)=qlg(:,kl-1:2:-1)+ddts*(mflx(:,kl-2:1:-1)*qlup(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1) &
                                       +mflx(:,kl-1:2:-1)*qlup(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)      &
                                       -mflx(:,kl-1:2:-1)*qlup(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1) &
                                       -mflx(:,kl:3:-1)*qlup(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1))
dd_a(:,kl,2)=qlg(:,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                           +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))

! qf (part of qtot) - atmosphere
dd_a(:,1,3)=qfg(:,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                           +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
dd_a(:,2:kl-1,3)=qfg(:,kl-1:2:-1)+ddts*(mflx(:,kl-2:1:-1)*qfup(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1) &
                                       +mflx(:,kl-1:2:-1)*qfup(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)      &
                                       -mflx(:,kl-1:2:-1)*qfup(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1) &
                                       -mflx(:,kl:3:-1)*qfup(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1))
dd_a(:,kl,3)=qfg(:,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                           +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))

! stratcloud - atmosphere
dd_a(:,1,4)=stratcloud(:,kl)+ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                                  +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
dd_a(:,2:kl-1,4)=stratcloud(:,kl-1:2:-1)+ddts*(mflx(:,kl-2:1:-1)*cfup(:,kl-2:1:-1)*(1.-fzzh(:,kl-2:1:-1))*idzm(:,kl-1:2:-1) &
                                              +mflx(:,kl-1:2:-1)*cfup(:,kl-1:2:-1)*fzzh(:,kl-2:1:-1)*idzm(:,kl-1:2:-1)      &
                                              -mflx(:,kl-1:2:-1)*cfup(:,kl-1:2:-1)*(1.-fzzh(:,kl-1:2:-1))*idzp(:,kl-1:2:-1) &
                                              -mflx(:,kl:3:-1)*cfup(:,kl:3:-1)*fzzh(:,kl-1:2:-1)*idzp(:,kl-1:2:-1))
dd_a(:,kl,4)=stratcloud(:,1)-ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                                  +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))

call thomas(tt_a(:,:,1:4),aa_a,bb_a,cc_a,dd_a(:,:,1:4),imax,kl,4)

qvg(:,:) = tt_a(:,kl:1:-1,1)
qlg(:,:) = tt_a(:,kl:1:-1,2)
qfg(:,:) = tt_a(:,kl:1:-1,3)
stratcloud(:,:) = min( max( tt_a(:,kl:1:-1,4), 0. ), 1. )


#ifdef scm
wqvflux(:,1)=wq0_a(:)
wqvflux(:,2:kl)=-kmo_a(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)         &
                +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))   &
                +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
wqlflux(:,1)=0.
wqlflux(:,2:kl)=-kmo_a(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)         &
                +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))   &
                +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
wqfflux(:,1)=0.
wqfflux(;,2:kl)=-kmo_a(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)         &
                +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))   &
                +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif


!------------------------------------------------------------------------------

! set-up decoupled ocean matrices

where ( deptho_dz(:,2)*deptho_dz(:,1)>1.e-4 )
  cc_o(:,1) = -ddts*ks_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
elsewhere
  cc_o(:,1) = 0.
end where
bb_o(:,1) = 1. - cc_o(:,1)
where ( deptho_dz(:,1:wlev-2)*deptho_dz(:,2:wlev-1)>1.e-4 ) 
  aa_o(:,2:wlev-1) = -ddts*ks_o(:,2:wlev-1)/(deptho_dz_hl(:,2:wlev-1)*deptho_dz(:,2:wlev-1))
elsewhere
  aa_o(:,2:wlev-1) = 0. 
end where
where ( deptho_dz(:,3:wlev)*deptho_dz(:,2:wlev-1)>1.e-4 )
  cc_o(:,2:wlev-1) = -ddts*ks_o(:,3:wlev)/(deptho_dz_hl(:,3:wlev)*deptho_dz(:,2:wlev-1))
elsewhere
  cc_o(:,2:wlev-1) = 0.
end where
bb_o(:,2:wlev-1) = 1. - aa_o(:,2:wlev-1) - cc_o(:,2:wlev-1)
where ( deptho_dz(:,wlev-1)*deptho_dz(:,wlev)>1.e-4 )
  aa_o(:,wlev) = -ddts*ks_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
elsewhere
  aa_o(:,wlev) = 0.
end where
bb_o(:,wlev) = 1. - aa_o(:,wlev)


! sal - ocean
fulldeptho(:) = 0.
targetdeptho(:) = 0.
do k = 1,wlev
  targetdeptho(:) = min( targetdeptho(:) + deptho_dz(:,k),  minsfc )
end do
do k = 1,wlev
  do iq = 1,imax
    dd_o(iq,k,1) = w_s(iq,k) + ddts*rhs_o(iq,k)*ws0_o(iq)
    if ( deptho_dz(iq,k)>1.e-4 ) then
      deltazo = max( min( deptho_dz(iq,k), targetdeptho(iq)-fulldeptho(iq) ), 0. )
      fulldeptho(iq) = fulldeptho(iq) + deltazo
      dd_o(iq,k,1) = dd_o(iq,k,1) - ddts*ws0subsurf_o(iq)*deltazo/max(deptho_dz(iq,k)*targetdeptho(iq),1.e-4)
    end if
  end do
end do
where ( .not.land(:) )
  dd_o(:,1,1) = dd_o(:,1,1) - ddts*ws0_o(:)/deptho_dz(:,1)
end where

call thomas(w_s,aa_o,bb_o,cc_o,dd_o(:,:,1),imax,wlev)


!------------------------------------------------------------------------------

! momentum for atmosphere, ocean and sea-ice ----

t1(:) = rhos(:)*cd_water(:)
t3(:) = rhos(:)*cd_ice(:)
t4(:) = wrtrho*cdbot_ice(:)

cc_a(:,1) = qq(:,kl)
bb_a(:,1) = 1.-qq(:,kl)
aa_a(:,2:kl-1) = rr(:,kl-1:2:-1)
cc_a(:,2:kl-1) = qq(:,kl-1:2:-1)
bb_a(:,2:kl-1) = 1.-qq(:,kl-1:2:-1)-rr(:,kl-1:2:-1)
aa_a(:,kl) = rr(:,1)
bb_a(:,kl) = 1.-rr(:,1)
where ( land(:) )
  bb_a(:,kl) = bb_a(:,kl) + ddts*rhos(:)*cduv(:)/(rhoa1(:)*dz_fl1(:)) ! implicit  
elsewhere
  bb_a(:,kl) = bb_a(:,kl) + ddts*(t1(:)*(1.-fracice(:))+t3(:)*fracice(:))/(rhoa1(:)*dz_fl1(:))
end where  

where ( deptho_dz(:,2)*deptho_dz(:,1) > 1.e-4 )
  cc_o(:,1) = -ddts*km_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
elsewhere
  cc_o(:,1) = 0.
end where
bb_o(:,1) = 1. - cc_o(:,1)
where ( deptho_dz(:,1:wlev-2)*deptho_dz(:,2:wlev-1) > 1.e-4 ) 
  aa_o(:,2:wlev-1) = -ddts*km_o(:,2:wlev-1)/(deptho_dz_hl(:,2:wlev-1)*deptho_dz(:,2:wlev-1))
elsewhere
  aa_o(:,2:wlev-1) = 0. 
end where
where ( deptho_dz(:,3:wlev)*deptho_dz(:,2:wlev-1) > 1.e-4 )
  cc_o(:,2:wlev-1) = -ddts*km_o(:,3:wlev)/(deptho_dz_hl(:,3:wlev)*deptho_dz(:,2:wlev-1))
elsewhere
  cc_o(:,2:wlev-1) = 0.
end where
bb_o(:,2:wlev-1) = 1. - aa_o(:,2:wlev-1) - cc_o(:,2:wlev-1)
where ( deptho_dz(:,wlev-1)*deptho_dz(:,wlev) > 1.e-4 )
  aa_o(:,wlev) = -ddts*km_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
elsewhere
  aa_o(:,wlev) = 0.
end where
bb_o(:,wlev) = 1. - aa_o(:,wlev)

where ( .not.land(:) )
  bb_o(:,1) = bb_o(:,1) + ddts*(t1(:)*(1.-fracice(:))+t4(:)*fracice(:))/(wrtrho*deptho_dz(:,1))
end where
! bottom drag
do iq = 1,imax
  if ( .not.land(iq) ) then  
    k = ibot(iq)
    bb_o(iq,k) = bb_o(iq,k) + ddts*cdbot_water(iq)/deptho_dz(iq,k) 
  end if
end do

where ( .not.land(:) )
  bb_i(:) = 1. + ddts*(t3(:)+t4(:))/imass(:)
elsewhere
  bb_i(:) = 1.
end where

where ( .not.land(:) )
  f_ao(:) = -ddts*t1(:)*(1.-fracice(:))/(rhoa1(:)*dz_fl1(:))
  f_oa(:) = -ddts*t1(:)*(1.-fracice(:))/(wrtrho*deptho_dz(:,1))
  f_ai(:) = -ddts*t3(:)*fracice(:)/(rhoa1(:)*dz_fl1(:))
  f_ia(:) = -ddts*t3(:)/imass(:)
  f_oi(:) = -ddts*t4(:)*fracice(:)/(wrtrho*deptho_dz(:,1))
  f_io(:) = -ddts*t4(:)/imass(:)
elsewhere
  f_ao = 0.
  f_oa = 0.
  f_ai = 0.
  f_ia = 0.
  f_oi = 0.
  f_io = 0.
end where  

! ua, uo, ui - coupled ----
! va, vo, vi - coupled ----
! k=1 is top of atmosphere and k=kl is bottom of atmosphere
dd_a(:,:,1) = ua(:,kl:1:-1)
dd_a(:,:,2) = va(:,kl:1:-1)
dd_o(:,:,1) = w_u(:,:)
dd_o(:,:,2) = w_v(:,:)
dd_i(:,1) = i_u(:)
dd_i(:,2) = i_v(:)

call solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a(:,:,1),dd_a(:,:,2),  &
                            tt_a(:,:,1),tt_a(:,:,2),                 &
                            aa_o,bb_o,cc_o,dd_o(:,:,1),dd_o(:,:,2),  &
                            w_u,w_v,                                 &
                            bb_i,dd_i(:,1),dd_i(:,2),                &
                            i_u(:),i_v(:),                           &
                            f_ao,f_oa,f_ai,f_ia,f_oi,f_io,           &
                            imax,kl,wlev) 

ua(:,:) = tt_a(:,kl:1:-1,1)
va(:,:) = tt_a(:,kl:1:-1,2)


! update surface momentum flux

where ( land(:) )
  ustar(:) = sqrt(cduv(:)*sqrt(ua(:,1)**2+va(:,1)**2))
elsewhere
  ustar(:) = sqrt((1.-fracice(:))*cd_water(:)*sqrt((ua(:,1)-w_u(:,1))**2+(va(:,1)-w_v(:,1))**2)  &
           +fracice(:)*cd_ice(:)*sqrt((ua(:,1)-i_u(:))**2+(va(:,1)-i_v(:))**2))
end where


#ifdef scm
uwflux(:,1) = 0.
uwflux(:,2:kl) = -kmo_a(:,1:kl-1)*(ua(:,2:kl)-ua(:,1:kl-1))/dz_hl(:,1:kl-1)
vwflux(:,1) = 0.
vwflux(:,2:kl) = -kmo_a(:,1:kl-1)*(va(:,2:kl)-va(:,1:kl-1))/dz_hl(:,1:kl-1)
wthflux(:,1) = wthlflux(:,1) - ((1.-fzzh(:,1)+sigkap(1)*fzzh(:,1)) &
                               *(lv*wqlflux(:,1)+ls*wqfflux(:,1)))
do k = 2,kl
  wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                                *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
end do
#endif


return
end subroutine update_coupled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve Sherman-Morrison formula

subroutine solve_sherman_morrison_2(aa_a,bb_a,cc_a,dd_a,tt_a, &
                                    aa_o,bb_o,cc_o,dd_o,tt_o, &
                                    f_ao,f_oa,                &
                                    imax,kl,wlev)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif

implicit none

integer, intent(in) :: imax, kl, wlev
integer k, iq
real, dimension(imax,2:kl), intent(in) :: aa_a
real, dimension(imax,kl), intent(in) :: bb_a, dd_a
real, dimension(imax,kl-1), intent(in) :: cc_a
real, dimension(imax,kl), intent(out) :: tt_a
real, dimension(imax,2:wlev), intent(in) :: aa_o
real, dimension(imax,wlev), intent(in) :: bb_o, dd_o
real, dimension(imax,wlev-1), intent(in) :: cc_o
real, dimension(imax,wlev), intent(out) :: tt_o
real, dimension(imax,2:kl+wlev) :: aad
real, dimension(imax,kl+wlev) :: bbd, ddd
real, dimension(imax,kl+wlev-1) :: ccd
real, dimension(imax,kl+wlev) :: xx
real, dimension(imax), intent(in) :: f_ao, f_oa

! Atmosphere and ocean (no sea-ice) large tridiagonal matrix

! Coupling matrix (inverted levels)
! [ bb_a cc_a                               ] [ t_a ]   [ dd_a ]
! [ aa_a bb_a cc_a                          ] [ t_a ]   [ dd_a ]
! [      .... .... ....                     ] [ ... ]   [ .... ]
! [           aa_a bb_a f_ao                ] [ t_a ]   [ dd_a ]
! [                f_oa bb_o cc_o           ] [ t_o ]   [ dd_i ]
! [                     aa_o bb_o cc_o      ] [ t_o ]   [ dd_i ]
! [                          .... .... .... ] [ ... ] = [ .... ]
! [                               aa_o bb_o ] [ t_o ]   [ dd_o ]

! Construct A
aad(:,2:kl) = aa_a(:,2:kl)
bbd(:,1:kl) = bb_a(:,1:kl)
ccd(:,1:kl-1) = cc_a(:,1:kl-1)
ddd(:,1:kl) = dd_a(:,1:kl)
ccd(:,kl) = f_ao
aad(:,kl+1) = f_oa
aad(:,kl+2:kl+wlev) = aa_o(:,2:wlev)
bbd(:,kl+1:kl+wlev) = bb_o(:,1:wlev)
ccd(:,kl+1:kl+wlev-1) = cc_o(:,1:wlev-1)
ddd(:,kl+1:kl+wlev) = dd_o(:,1:wlev)

! pure tridiagonal matrix for atmosphere and ocean (no sea-ice)
call thomas(xx(:,1:kl+wlev),aad(:,2:kl+wlev),bbd(:,1:kl+wlev),ccd(:,1:kl+wlev-1),ddd(:,1:kl+wlev),imax,kl+wlev)

! Unpack solution
tt_a(:,:) = xx(:,1:kl)
tt_o(:,:) = xx(:,kl+1:kl+wlev)

return
end subroutine solve_sherman_morrison_2

subroutine solve_sherman_morrison_3(aa_a,bb_a,cc_a,dd_au,dd_av, &
                                    tt_au,tt_av,                &
                                    aa_o,bb_o,cc_o,dd_ou,dd_ov, &
                                    tt_ou,tt_ov,                &
                                    bb_i,dd_iu,dd_iv,           &
                                    tt_iu,tt_iv,                &
                                    f_ao,f_oa,f_ai,f_ia,        &
                                    f_oi,f_io,                  &
                                    imax,kl,wlev)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif

implicit none

integer, intent(in) :: imax, kl, wlev
integer k, n, iq
real, dimension(imax,2:kl), intent(in) :: aa_a
real, dimension(imax,kl), intent(in) :: bb_a, dd_au, dd_av
real, dimension(imax,kl-1), intent(in) :: cc_a
real, dimension(imax,kl), intent(out) :: tt_au, tt_av
real, dimension(imax,2:wlev), intent(in) :: aa_o
real, dimension(imax,wlev), intent(in) :: bb_o, dd_ou, dd_ov
real, dimension(imax,wlev-1), intent(in) :: cc_o
real, dimension(imax,wlev), intent(out) :: tt_ou, tt_ov
real, dimension(imax), intent(in) :: bb_i, dd_iu, dd_iv
real, dimension(imax), intent(out) :: tt_iu, tt_iv
real, dimension(imax,2:kl+wlev+1) :: aad
real, dimension(imax,kl+wlev+1) :: bbd
real, dimension(imax,kl+wlev) :: ccd
real, dimension(imax,kl+wlev+1,3) :: ddd
real, dimension(imax,kl+wlev+1) :: uu, vv
real, dimension(imax,kl+wlev+1,3) :: yy
real, dimension(imax), intent(in) :: f_ao, f_oa, f_ai, f_ia, f_oi, f_io
real, dimension(imax,3) :: tmp

! Coupling matrix (inverted levels)
! [ bb_a cc_a                                    ] [ t_a ]   [ dd_a ]
! [ aa_a bb_a cc_a                               ] [ t_a ]   [ dd_a ]
! [      .... .... ....                          ] [ ... ]   [ .... ]
! [           aa_a bb_a f_ao                f_ai ] [ t_a ]   [ dd_a ]
! [                f_oa bb_o cc_o           f_oi ] [ t_o ]   [ dd_o ]
! [                     aa_o bb_o cc_o           ] [ t_o ]   [ dd_o ]
! [                          .... .... ....      ] [ ... ] = [ .... ]
! [                               aa_o bb_o 0.   ] [ t_o ]   [ dd_o ]
! [                f_ia f_io           0.   bb_i ] [ t_i ]   [ dd_i ]

! u^T = [ ..  f_ai  f_oi .. -1 ]
! v^T = [ .. -f_ia -f_io ..  1 ]


! Use Sherman-Morrison formula
! A t = d
! (A' + u v^T) t = d
! A' y = d
! A' q = u
! t = y - {(v^t y)/(1 + (v^t q))} q


! Construct A'
aad(:,2:kl) = aa_a(:,2:kl)
bbd(:,1:kl-1) = bb_a(:,1:kl-1)
ccd(:,1:kl-1) = cc_a(:,1:kl-1)
ddd(:,1:kl,1) = dd_au(:,1:kl)
ddd(:,1:kl,2) = dd_av(:,1:kl)
bbd(:,kl) = bb_a(:,kl)  + f_ai*f_ia ! due to u v^T matrix
ccd(:,kl) = f_ao        + f_ai*f_io ! due to u v^T matrix
aad(:,kl+1) = f_oa      + f_oi*f_ia ! due to u v^T matrix
bbd(:,kl+1) = bb_o(:,1) + f_oi*f_io ! due to u v^T matrix
aad(:,kl+2:kl+wlev) = aa_o(:,2:wlev)
bbd(:,kl+2:kl+wlev) = bb_o(:,2:wlev)
ccd(:,kl+1:kl+wlev-1) = cc_o(:,1:wlev-1)
ddd(:,kl+1:kl+wlev,1) = dd_ou(:,1:wlev)
ddd(:,kl+1:kl+wlev,2) = dd_ov(:,1:wlev)
ccd(:,kl+wlev) = 0.
aad(:,kl+wlev+1) = 0.
bbd(:,kl+wlev+1) = bb_i(:) + 1. ! note +1. due to u v^T matrix
ddd(:,kl+wlev+1,1) = dd_iu(:)
ddd(:,kl+wlev+1,2) = dd_iv(:)

! Construct u and v
uu(:,:) = 0.
uu(:,kl) = f_ai
uu(:,kl+1) = f_oi
uu(:,kl+wlev+1) = -1.
vv(:,:) = 0.
vv(:,kl) = -f_ia
vv(:,kl+1) = -f_io
vv(:,kl+wlev+1) = 1.
ddd(:,:,3) = uu(:,:)

! Solve A' y = d  and  A' q = u for combined atmosphere, ocean and sea-ice
call thomas(yy(:,1:kl+wlev+1,1:3),aad(:,2:kl+wlev+1),bbd(:,1:kl+wlev+1),ccd(:,1:kl+wlev),ddd(:,1:kl+wlev+1,1:3),imax,kl+wlev+1,3)

! Solve for x = y - {(v^t y)/(1 + (v^t q))} q
do n = 1,3
  do iq = 1,imax
    tmp(iq,n) = vv(iq,kl)*yy(iq,kl,n) + vv(iq,kl+1)*yy(iq,kl+1,n) &
              + vv(iq,kl+wlev+1)*yy(iq,kl+wlev+1,n)
  end do
end do

! Unpack solution
do k = 1,kl
  do iq = 1,imax
    tt_au(iq,k) = yy(iq,k,1) - yy(iq,k,3)*tmp(iq,1)/(1.+tmp(iq,3))
    tt_av(iq,k) = yy(iq,k,2) - yy(iq,k,3)*tmp(iq,2)/(1.+tmp(iq,3))
  end do
end do
do k = 1,wlev
  do iq = 1,imax
    tt_ou(iq,k) = yy(iq,kl+k,1) - yy(iq,kl+k,3)*tmp(iq,1)/(1.+tmp(iq,3))
    tt_ov(iq,k) = yy(iq,kl+k,2) - yy(iq,kl+k,3)*tmp(iq,2)/(1.+tmp(iq,3))
  end do
end do
tt_iu(:) = yy(:,kl+wlev+1,1) - yy(:,kl+wlev+1,3)*tmp(:,1)/(1.+tmp(:,3))
tt_iv(:) = yy(:,kl+wlev+1,2) - yy(:,kl+wlev+1,3)*tmp(:,2)/(1.+tmp(:,3))

return
end subroutine solve_sherman_morrison_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

pure subroutine thomas1(outdat,aai,bbi,cci,ddi,imax,klin)
#ifdef GPUPHYSICS
!!$acc routine vector
#endif

implicit none

integer, intent(in) :: imax, klin
real, dimension(imax,2:klin), intent(in) :: aai
real, dimension(imax,klin), intent(in) :: bbi,ddi
real, dimension(imax,klin-1), intent(in) :: cci
real, dimension(imax,klin), intent(out) :: outdat
real, dimension(imax,klin) :: cc,dd
real n_s
integer k, iq

do iq = 1,imax
  cc(iq,1) = cci(iq,1)/bbi(iq,1)
  dd(iq,1) = ddi(iq,1)/bbi(iq,1)
end do

do k = 2,klin-1
  do iq = 1,imax
    n_s = 1./(bbi(iq,k)-cc(iq,k-1)*aai(iq,k))
    cc(iq,k) = cci(iq,k)*n_s
    dd(iq,k) = (ddi(iq,k)-dd(iq,k-1)*aai(iq,k))*n_s
  end do
end do
do iq = 1,imax
  n_s = 1./(bbi(iq,klin)-cc(iq,klin-1)*aai(iq,klin))
  outdat(iq,klin) = (ddi(iq,klin)-dd(iq,klin-1)*aai(iq,klin))*n_s
end do
do k = klin-1,1,-1
  outdat(:,k) = dd(:,k)-cc(:,k)*outdat(:,k+1)
end do

return
end subroutine thomas1

pure subroutine thomas2(outdat,aai,bbi,cci,ddi,imax,klin,ndim)
!!$acc routine vector

implicit none

integer, intent(in) :: imax, klin, ndim
real, dimension(imax,2:klin), intent(in) :: aai
real, dimension(imax,klin), intent(in) :: bbi
real, dimension(imax,klin-1), intent(in) :: cci
real, dimension(imax,klin,ndim), intent(in) :: ddi
real, dimension(imax,klin,ndim), intent(out) :: outdat
integer n, iq, k
real, dimension(imax,ndim,klin) :: cc,dd
real n_s

! Thomas
do n = 1,ndim
  do iq = 1,imax
    cc(iq,n,1) = cci(iq,1)/bbi(iq,1)
    dd(iq,n,1) = ddi(iq,1,n)/bbi(iq,1)
  end do
end do

do k = 2,klin-1
  do n = 1,ndim
    do iq = 1,imax
      n_s = 1./(bbi(iq,k)-cc(iq,n,k-1)*aai(iq,k))
      cc(iq,n,k) = cci(iq,k)*n_s
      dd(iq,n,k) = (ddi(iq,k,n)-dd(iq,n,k-1)*aai(iq,k))*n_s
    end do
  end do  
end do
do n = 1,ndim
  do iq = 1,imax
    n_s = 1./(bbi(iq,klin)-cc(iq,n,klin-1)*aai(iq,klin))
    outdat(iq,klin,n) = (ddi(iq,klin,n)-dd(iq,n,klin-1)*aai(iq,klin))*n_s
  end do  
end do
do k = klin-1,1,-1
  do n = 1,ndim
    do iq = 1,imax
      outdat(iq,k,n) = dd(iq,n,k)-cc(iq,n,k)*outdat(iq,k+1,n)
    end do
  end do
end do

return
end subroutine thomas2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

pure subroutine getqsat2(qsat,templ,ps,fice)
!!$acc routine vector

implicit none

integer iq, ix, ixx
real, dimension(:), intent(in) :: templ
real, dimension(:), intent(in) :: ps, fice
real, dimension(:), intent(out) :: qsat
real estafi, tdiff, tdiffx, rx, rxx
real qsatl, qsati, deles

do iq = 1,size(qsat)
  tdiff = min(max( templ(iq)-123.16, 0.), 219.)
  rx = tdiff - aint(tdiff)
  ix = int(tdiff)
  estafi = (1.-rx)*tablei(ix) + rx*tablei(ix+1)
  qsati = 0.622*estafi/max(ps(iq)-estafi,0.1)

  tdiffx = min(max( templ(iq)-273.1, -40.), 1.)
  rxx = tdiffx - aint(tdiffx)
  ixx = int(tdiffx)
  deles = (1.-rxx)*esdiff(ixx) + rxx*esdiff(ixx+1)
  qsatl = qsati + 0.622*deles/ps(iq)

  qsat(iq) = fice(iq)*qsati + (1.-fice(iq))*qsatl
end do

return
end subroutine getqsat2

pure subroutine getqsat3(qsat,templ,ps,fice)
!!$acc routine vector

implicit none

integer iq, ix, ixx, k
real, dimension(:,:), intent(in) :: templ
real, dimension(:,:), intent(in) :: ps, fice
real, dimension(:,:), intent(out) :: qsat
real estafi, tdiff, tdiffx, rx, rxx
real qsatl, qsati, deles

do k = 1,size(qsat,2)
  do iq = 1,size(qsat,1)
    tdiff = min(max( templ(iq,k)-123.16, 0.), 219.)
    rx = tdiff - aint(tdiff)
    ix = int(tdiff)
    estafi = (1.-rx)*tablei(ix) + rx*tablei(ix+1)
    qsati = 0.622*estafi/max(ps(iq,k)-estafi,0.1)

    tdiffx = min(max( templ(iq,k)-273.1, -40.), 1.)
    rxx = tdiffx - aint(tdiffx)
    ixx = int(tdiffx)
    deles = (1.-rxx)*esdiff(ixx) + rxx*esdiff(ixx+1)
    qsatl = qsati + 0.622*deles/ps(iq,k)

    qsat(iq,k) = fice(iq,k)*qsati + (1.-fice(iq,k))*qsatl
  end do  
end do

return
end subroutine getqsat3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

function entfn(zht,zi,imax) result(ans)
!!$acc routine vector

implicit none

integer, intent(in) :: imax
real, dimension(imax), intent(in) :: zht, zi
real, dimension(imax) :: ans

!ans=0.002                                            ! Angevine (2005)
!ans=2./max(100.,zi)                                  ! Angevine et al (2010)
!ans=1./zht                                           ! Siebesma et al (2003)
!ans=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))    ! Soares et al (2004)
ans = max( ent0/max( zht, 1. ) + ent1/max( zi-zht, ezmin ), ent_min )

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phi_m for atmospheric stability

pure subroutine calc_phi(phim,z_on_l,imax)
!!$acc routine vector

implicit none

integer, intent(in) :: imax
real, dimension(imax), intent(in) :: z_on_l
real, dimension(imax), intent(out) :: phim

select case(stabmeth)
  case default
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update exponential weighted moving average

subroutine update_ema(field,ema,dt)

implicit none

integer imax, kx, k, iq
real, intent(in) :: dt
real alpha, nstep
real, dimension(:,:), intent(in) :: field
real, dimension(:,:), intent(inout) :: ema

imax = min( size(field,1), size(ema,1) )
kx = size(field,2)

nstep = max( tke_timeave_length/dt, 1. ) ! this is a real value
alpha = 2./(nstep + 1.) 

do k = 1,kx
  do iq = 1,imax
    ema(iq,k) = alpha*field(iq,k) + (1.-alpha)*ema(iq,k)
  end do
end do
  
return
end subroutine update_ema

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE-eps

subroutine tkeend

implicit none

deallocate(tke,eps)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps

