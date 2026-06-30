! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! clouds and saturated air following Marquet and Geleyn 2012.  Scale aware mass flux is based on
! Boutle et al 2014.

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
public tkemeth, tcalmeth

real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps

! model ED constants
real, save :: cm0      = 0.09      ! Hurley (2007) 0.09, Duynkerke (1988) 0.03, Duynkerke (1987) 0.09
real, save :: ce0      = 0.69      ! Hurley (2007) 0.69, Duynkerke (1988) 0.42, Duynkerke (1987) 0.77
real, save :: ce1      = 1.46
real, save :: ce2      = 1.83
real, save :: ce3      = 0.45      ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, save :: mintke   = 1.5e-4    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps   = 1.e-6     ! min value for eps (1.0e-6 in TAPM)
real, save :: minl     = 5.        ! min value for L   (5. in TAPM)
real, save :: maxl     = 500.      ! max value for L   (500. in TAPM)
! model MF constants
real, save :: be       = 1.        ! Surface boundary condition (Hurley (2007) 1., Soares et al (2004) 0.3)
real, save :: ent0     = 0.5       ! Entrainment constant (Controls height of boundary layer) (Hurley (2007) 0.5)
real, save :: ent1     = 0.  
real, save :: ent_min  = 0.001     ! Minimum entrainment
real, save :: ezmin    = 10.       ! Limits entrainment at plume top
real, save :: entc0    = 2.e-3     ! Saturated entrainment constant for mass flux
real, save :: dtrc0    = 3.e-3     ! Saturated detrainment constant for mass flux
real, save :: m0       = 0.1       ! Mass flux area constant (Hurley (2007) 0.1)
real, save :: b1       = 2.        ! Updraft entrainment coeff (Soares et al (2004) 1., Siebesma et al (2003) 2.)
real, save :: b2       = 0.3333    ! Updraft buoyancy coeff (Soares et al (2004) 2., Siebesma et al (2003) 1./3.)
real, save :: qcmf     = 1.e-4     ! Critical mixing ratio of liquid water before autoconversion
real, save :: mfbeta   = 0.15      ! Horizontal scale factor
! generic constants
integer, save :: buoymeth = 1      ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0      ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: tkemeth  = 1      ! Method for TKE calculation (0=D&K84, 1=Hurley)
integer, save :: tcalmeth = 1      ! Method for correcting saturated air (0=Remove, 1=Retain, 2=Remove below pbl)
real, save :: maxdts      = 120.   ! max timestep for split
! wind gusts
real, parameter :: cs1 = 2.2
real, parameter :: cs2 = 1.63
real, parameter :: cs3 = 0.73
real, parameter :: cw1 = 1.
real, parameter :: cw2 = 0.24
real, parameter :: cw3 = 0.

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

real, dimension(0:220), save :: tablei 
real, dimension(-40:2), save :: esdiff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(shear(ifull,kl))

tke(1:ifull+iextra,1:kl)=mintke
eps(1:ifull+iextra,1:kl)=mineps
shear(1:ifull,1:kl)=0.

tablei = &
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
esdiff = &
(/ 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,  &
   13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61, &
   22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27, &
   26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65, &
   0.08, 0., 0. /)

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux
! mode=2 combined atmosphere and ocean with mass flux
! mode=3 combined atm-ocn and no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,ni,cf,ua,va,zi,fg,eg,cduv,                &
                  ps,zz,zzh,sig,rhos,ustar_ave,dt,mode,tke,eps,shear,dx,land,     &
                  w_t,w_s,w_u,w_v,w_tke,w_eps,deptho_depth,deptho_depth_hl,       &
                  deptho_dz,deptho_dz_hl,dwdx_o,dwdy_o,ibot,ubot_o,vbot_o,        &
                  utop_o,vtop_o,cd_water,cdh_water,cdbot_water,wt0_o,wt0rad_o,    &
                  wt0melt_o,wt0eg_o,wt0fb_o,ws0_o,ws0subsurf_o,wu0_o,wv0_o,zo_o,  &
                  rad_o,i_u,i_v,fracice,icefg_a,imass,cd_ice,cdbot_ice,imax,kl,ol)

use mlo_ctrl, only : mlocheck, kemaxdt, maxsal, wrtrho, fluxwgt, omink,  &
                     omineps, limitL, omaxl, ominl, fixedstabfunc, nops, &
                     nopb, fixedce3, eps_mode, k_mode, minsfc, cp0,      &
                     wrtemp

implicit none

integer, intent(in) :: imax, kl, ol, mode
real, intent(in) :: dt
real, dimension(:,:), intent(inout) :: theta,cf,ua,va
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg,ni
real, dimension(:,:), intent(out) :: kmo
real, dimension(:,:), intent(in) :: zz,zzh
real, dimension(:,:), intent(in) :: shear
real, dimension(:,:), intent(inout) :: tke, eps
real, dimension(:), intent(inout) :: zi,fg,eg
real, dimension(:), intent(in) :: cduv,ps,rhos,dx
real, dimension(:), intent(out) :: ustar_ave
real, dimension(:), intent(in) :: sig

integer, dimension(:), intent(in) :: ibot
real, dimension(:,:), intent(in) :: deptho_depth
real, dimension(:,:), intent(in) :: deptho_depth_hl
real, dimension(:,:), intent(in) :: deptho_dz
real, dimension(:,2:), intent(in) :: deptho_dz_hl
real, dimension(:,:), intent(in) :: rad_o
real, dimension(:,:), intent(inout) :: w_t, w_s, w_u, w_v
real, dimension(:), intent(inout) :: ubot_o, vbot_o, utop_o, vtop_o
real, dimension(:), intent(in) :: cd_water, cdh_water, cdbot_water
real, dimension(:), intent(in) :: wt0rad_o, wt0melt_o, wt0eg_o, wt0fb_o, ws0subsurf_o
real, dimension(:), intent(inout) :: wu0_o, wv0_o, wt0_o
real, dimension(:), intent(in) :: ws0_o, zo_o
real, dimension(:,:), intent(in) :: dwdx_o, dwdy_o
real, dimension(:), intent(inout) :: i_u, i_v
real, dimension(:), intent(in) :: fracice
real, dimension(:), intent(in) :: icefg_a, imass, cd_ice, cdbot_ice
real, dimension(:,:), intent(inout) :: w_tke, w_eps
logical, dimension(:), intent(in) :: land

integer k, i, kcount, mcount
real, dimension(imax,kl) :: km, rhoa, pps, ppb
real, dimension(imax,kl) :: thetav
real, dimension(imax,kl) :: mflx, tlup, qvup, qlup, qfup, niup, cfup
real, dimension(imax,kl) :: thetal
real, dimension(imax,kl) :: idzm
real, dimension(imax,kl-1) :: dz_hl, idzp, fzzh
real, dimension(imax) :: z_on_l, phim, pres, wt0, wq0, wtv0
real, dimension(imax) :: cgmap, templ, fice, qc, qt, temp, zi_save
real, dimension(imax) :: wstar, qsat, ustar, fg_ave
real, dimension(kl) :: sigkap
real dz_fl     ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real tff, cm12, cm34, ddts, zturb
real rhoahl, rhoahlm1
real wdash_sq, l_on_kz
real lx, dqsdt, al
logical, dimension(imax) :: mask

real, dimension(imax,ol) :: pps_o, ppb_o, n2_o
real, dimension(imax,ol) :: km_o, ks_o
real, dimension(imax,ol) :: kmo_hl, kso_hl
real, dimension(imax,ol) :: fdeptho_hl
real ustar_o, uoave, voave, umag, zrough
real, parameter :: cu0 = 0.5562
real, parameter :: utide   = 0.05    ! Tide velocity for bottom drag (m/s)
real, parameter :: cdbot   = 2.4e-3  ! bottom drag coefficent
logical use_ocean

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))
use_ocean = mode==2 .or. mode==3
use_ocean = use_ocean .and. any(.not.land(1:imax))

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

call mlocheck("Before coupled-mixing",water_temp=w_t,water_sal=w_s,water_u=w_u, &
              water_v=w_v,ice_u=i_u,ice_v=i_v)

do k = 1,kl
  sigkap(k) = sig(k)**(-rd/cp)
end do
  
do k = 1,kl
  do i = 1,imax
     
    ! Transform to thetal as it is the conserved variable
    thetal(i,k) = theta(i,k) - sigkap(k)*(lv*qlg(i,k)+ls*qfg(i,k))/cp
            
    ! Impose limits on tke and eps after being advected by the host model
    tke(i,k) = max(tke(i,k), mintke)
    tff = cm34*tke(i,k)*sqrt(tke(i,k))
    eps(i,k) = min(eps(i,k), tff/minl)
    eps(i,k) = max(eps(i,k), tff/maxl, mineps)
  
    ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
    pres(i) = ps(i)*sig(k) ! pressure
    ! Density must be updated when dz is updated so that rho*dz is conserved
    thetav(i,k) = theta(i,k)*(1.+0.61*qvg(i,k)-qlg(i,k)-qfg(i,k))
    rhoa(i,k) = sigkap(k)*pres(i)/(rd*thetav(i,k))

    ! Calculate first approximation to diffusion coeffs
    km(i,k) = cm0*tke(i,k)**2/eps(i,k)
  end do
end do
  
do k = 1,kl-1
  do i = 1,imax
    ! Fraction for interpolation from full levels to half levels
    fzzh(i,k) = (zzh(i,k)-zz(i,k))/(zz(i,k+1)-zz(i,k))

    ! Calculate dz at half levels
    dz_hl(i,k) = max( zz(i,k+1) - zz(i,k), 1. )
  end do  
end do

do k = 1,kl-1
  do i = 1,imax
    ! interpolate diffusion coeff to half levels
    kmo(i,k) = km(i,k) + fzzh(i,k)*(km(i,k+1)-km(i,k))
  end do
end do
do i = 1,imax
  kmo(i,kl) = 0.
end do  

! eddy diffusion terms to account for air density with level thickness
do k = 2,kl
  do i = 1,imax 
    dz_fl = zzh(i,k) - zzh(i,k-1)
    !rhoahlm1 = rhoahl(:,k-1)
    rhoahlm1 = rhoa(i,k-1) + fzzh(i,k-1)*(rhoa(i,k)-rhoa(i,k-1))
    idzm(i,k) = rhoahlm1/(rhoa(i,k)*dz_fl)
  end do  
end do
do i = 1,imax
  rhoahl = rhoa(i,1) + fzzh(i,1)*(rhoa(i,2)-rhoa(i,1))
  idzp(i,1) = rhoahl/(rhoa(i,1)*zzh(i,1))
end do  
do k = 2,kl-1
  do i = 1,imax
    dz_fl = zzh(i,k) - zzh(i,k-1)
    rhoahl = rhoa(i,k) + fzzh(i,k)*(rhoa(i,k+1)-rhoa(i,k))
    idzp(i,k) = rhoahl/(rhoa(i,k)*dz_fl)
  end do  
end do

! additional ocean terms
do k = 2,ol
  do i = 1,imax
    fdeptho_hl(i,k) = (deptho_depth_hl(i,k)-deptho_depth(i,k-1)) &
                     /max(deptho_depth(i,k)-deptho_depth(i,k-1),1.e-8)
  end do
end do

! reset average ustar and fg
do i = 1,imax
  ustar_ave(i) = 0.
  fg_ave(i) = 0.
end do  


! Main loop to prevent time splitting errors
mcount = int(dt/(min(maxdts,12./max(m0,0.1))+0.01)) + 1
mcount = max( mcount, int(dt/(kemaxdt+0.01)) + 1 )
ddts   = dt/real(mcount)
do kcount = 1,mcount
    
  ! calculate theta and theta_v at k=1
  do k = 1,kl
    do i = 1,imax
      theta(i,k) = thetal(i,k) + sigkap(k)*(lv*qlg(i,k)+ls*qfg(i,k))/cp
    end do
  end do
  thetav(:,1) = theta(:,1)*(1.+0.61*qvg(:,1)-qlg(:,1)-qfg(:,1))
  
  ! Calculate surface fluxes
  wt0(:) = fg(:)/(rhos(:)*cp)  ! theta flux
  wq0(:) = eg(:)/(rhos(:)*lv)  ! qtot flux  
  wtv0(:) = wt0(:) + theta(:,1)*0.61*wq0(:) ! thetav flux
  
  ! Update momentum flux
  ustar(:) = sqrt(cduv(:)*sqrt(ua(:,1)**2+va(:,1)**2))  
  wstar(:) = max(grav*zi(:)*max(wtv0(:),0.)/thetav(:,1),1.e-10)**(1./3.)

    
  !=============================================================================
  ! Start: Plumerise
  !=============================================================================
        
  ! Calculate non-local mass-flux terms for theta and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  if ( mode/=1 .and. mode/=3 ) then ! mass flux

    zi_save(:) = zi(:)
    mask(:) = wtv0(:)>0.

      
    call plumerise(thetal,qvg,qlg,qfg,ni,cf,ua,va,                 &
                   zz,dz_hl,tke,theta,thetav,ps,mask,              &
                   zi,ustar,wstar,wt0,wq0,wtv0,sig,                &
                   sigkap,mflx,tlup,qvup,qlup,qfup,niup,cfup,cm12, &
                   imax,kl)
      

    ! Turn off MF term if small grid spacing (mfbeta=0 implies MF is always non-zero)
    ! Based on Boutle et al 2014
    do i = 1,imax
      zturb = 0.5*(zi_save(i) + zi(i))
      cgmap(i) = 1. - tanh(mfbeta*zturb/dx(i))*max(0.,1.-0.25*dx(i)/zturb)
    end do  
    do k = 1,kl
      do i = 1,imax  
        mflx(i,k) = mflx(i,k)*cgmap(i)
      end do  
    end do

  else
      
    do k = 1,kl
      do i = 1,imax
        mflx(i,k) = 0.
        tlup(i,k) = thetal(i,k)
        qvup(i,k) = qvg(i,k)
        qlup(i,k) = qlg(i,k)
        qfup(i,k) = qfg(i,k)
        niup(i,k) = ni(i,k)
        cfup(i,k) = cf(i,k)
      end do
    end do
      
  end if

  !==============================================================================
  ! End: Plumerise
  !==============================================================================

      
  ! calculate tke and eps boundary condition at 1st vertical level
  z_on_l(:) = -vkar*zz(:,1)*grav*wtv0(:)/(thetav(:,1)*max(ustar(:)**3,1.E-10))
  z_on_l(:) = min(z_on_l(:),10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  select case(stabmeth)
    case default
      where ( z_on_l(:)<0. )
        phim = (1.-16.*z_on_l)**(-0.25)
      elsewhere
        phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
      end where
    case(1)
      where ( z_on_l(:)<0. )
        phim = (1.-16.*z_on_l)**(-0.25)
      elsewhere (z_on_l<=0.4)
        phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
      elsewhere
        phim = aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar
      end where
  end select

  do i = 1,imax
    tke(i,1) = cm12*ustar(i)**2+ce3*wstar(i)**2
    eps(i,1) = ustar(i)**3*phim(i)/(vkar*zz(i,1))+grav*wtv0(i)/thetav(i,1)
    tke(i,1) = max( tke(i,1), mintke )
    tff = cm34*tke(i,1)*sqrt(tke(i,1))
    eps(i,1) = min( eps(i,1), tff/minl )
    eps(i,1) = max( eps(i,1), tff/maxl, mineps )
  end do
    
  do k = 1,kl
    do i = 1,imax
      ! update eddy diffusivity
      km(i,k) = cm0*tke(i,k)**2/eps(i,k)
    end do  
  end do

    
  !=============================================================================
  ! Start: calculate_buoyancy
  !=============================================================================

  call calculate_buoyancy(ppb,thetal,qvg,qlg,qfg,cf,   &
                          theta,thetav,km,zz,zzh,fzzh, &
                          sigkap,imax,kl)
    
  !=============================================================================
  ! End: calculate_buoyancy
  !=============================================================================

    
  do k = 1,kl
    do i = 1,imax
      pps(i,k) = km(i,k)*shear(i,k)
    end do
  end do
  do i = 1,imax
    pps(i,kl) = 0.
  end do


  !=============================================================================
  ! Start: update_tkeeps
  !=============================================================================

  call update_tkeeps(tke,eps,pps,ppb,wstar,zi(:), &
                     km,zz,dz_hl,idzm,idzp,fzzh,  &
                     ddts,cm34,imax,kl)
    
  !=============================================================================
  ! End: update_tkeeps
  !=============================================================================

    
  if ( use_ocean ) then
      

    !=============================================================================
    ! Start: mlo_update_keps
    !=============================================================================

    call mlo_update_n2(n2_o,w_t,w_s,deptho_depth,deptho_depth_hl, &
                       deptho_dz,fdeptho_hl,wrtrho,imax,ol)

    !update arrays based on current state
    do k = 1,ol
      do i = 1,imax
        if ( deptho_dz(i,k)<=1.e-4 ) then
          w_tke(i,k) = omink
          w_eps(i,k) = omineps
        end if
      end do  
    end do

    !boundary conditions
    do i = 1,imax
      ustar_o = max(sqrt(sqrt(wu0_o(i)**2+wv0_o(i)**2)),1.e-6) 
      w_tke(i,1) = (ustar_o/cu0)**2
      w_eps(i,1) = min((cu0)**3*w_tke(i,1)**1.5,1.e9)/(vkar*(0.5*max(deptho_dz(i,1),1.e-4)+zo_o(i)))
      k = ibot(i)
      uoave = fluxwgt*w_u(i,k) + (1.-fluxwgt)*ubot_o(i)
      voave = fluxwgt*w_v(i,k) + (1.-fluxwgt)*vbot_o(i)
      umag = sqrt(max(uoave**2+voave**2,1.e-4))
      umag = max( umag, utide )
      zrough = 0.5*max(deptho_dz(i,k),1.e-4)/exp(vkar/sqrt(cdbot))
      if ( k==1 ) then
        w_tke(i,1) = 0.5*( w_tke(i,1) + cdbot*(umag/cu0)**2 )
        w_eps(i,1) = 0.5*( w_eps(i,1) + min((cu0)**3*max(w_tke(i,1),1.e-8)**1.5,1.e9) &
                           /(vkar*(0.5*max(deptho_dz(i,1),1.e-4)+zrough)) )
      else
        w_tke(i,k) = cdbot*(umag/cu0)**2
        w_eps(i,k) = min((cu0)**3*max(w_tke(i,k),1.e-8)**1.5,1.e9)/(vkar*(0.5*max(deptho_dz(i,k),1.e-4)+zrough))
      end if
    end do
        
    call mlo_limit_tke(w_tke,w_eps,km_o,ks_o,n2_o,omink,omineps,omaxL,  &
                       ominL,limitL,fixedstabfunc,cu0,imax,ol)

    call mlo_update_buoyancy(pps_o,ppb_o,km_o,ks_o,w_u,w_v,             &
                             n2_o,dwdx_o,dwdy_o,deptho_dz,fdeptho_hl,   &
                             nops,nopb,imax,ol)
        
    !km & ks at half levels
    do i = 1,imax
      kmo_hl(i,1) = 1.e-6
      kso_hl(i,1) = 1.e-6
    end do
    do k = 2,ol
      do i = 1,imax  
        kmo_hl(i,k) = km_o(i,k-1) + fdeptho_hl(i,k)*(km_o(i,k)-km_o(i,k-1))
        kso_hl(i,k) = ks_o(i,k-1) + fdeptho_hl(i,k)*(ks_o(i,k)-ks_o(i,k-1))
      end do
    end do      
      
    call mlo_update_keps(kmo_hl,kso_hl,w_tke,w_eps,pps_o,ppb_o,   &
                         deptho_depth,deptho_depth_hl,deptho_dz,  &
                         deptho_dz_hl,ibot,                       &
                         ddts,fixedce3,eps_mode,k_mode,imax,ol)      
      
    call mlo_limit_tke(w_tke,w_eps,km_o,ks_o,n2_o,omink,omineps,omaxL,  &
                       ominL,limitL,fixedstabfunc,cu0,imax,ol)

    !km & ks at half levels
    do i = 1,imax
      kmo_hl(i,1) = 1.e-6
      kso_hl(i,1) = 1.e-6
    end do
    do k = 2,ol
      do i = 1,imax  
        kmo_hl(i,k) = km_o(i,k-1) + fdeptho_hl(i,k)*(km_o(i,k)-km_o(i,k-1))
        kso_hl(i,k) = ks_o(i,k-1) + fdeptho_hl(i,k)*(ks_o(i,k)-ks_o(i,k-1))
      end do
    end do
      
    ! store currents for time averaging at next time-step  
    ! Assume ocean mixing occurs after this routine is called
    do i = 1,imax
      utop_o(i) = w_u(i,1)  
      vtop_o(i) = w_v(i,1)
      k = ibot(i)
      ubot_o(i) = w_u(i,k)
      vbot_o(i) = w_v(i,k)
    end do  

    !=============================================================================
    ! End: mlo_update_keps
    !=============================================================================
 

  end if  

    
  !=============================================================================
  ! Start: update_coupled
  !=============================================================================

  call update_coupled(thetal,qvg,qlg,qfg,ni,cf,ua,va,            &
                      tlup,qvup,qlup,qfup,niup,cfup,fg,eg,       &
                      rhos(:),ustar,cduv(:),tke,eps,mflx,fzzh,   &
                      idzp,idzm,dz_hl,rhoa(:,1),zzh(:,1),        &
                      deptho_dz,deptho_dz_hl,rad_o,              &
                      w_t,w_s,w_u,w_v,                           &
                      cd_water(:),cdh_water(:),cdbot_water(:),   &
                      wt0rad_o(:),wt0melt_o(:),wt0eg_o(:),       &
                      icefg_a(:),wt0fb_o(:),ws0_o(:),            &
                      ws0subsurf_o(:),i_u(:),i_v(:),             &
                      imass(:),fracice(:),cd_ice(:),             &
                      cdbot_ice(:),ibot(:),land(:),wu0_o(:),     &
                      wv0_o(:),wt0_o(:),kmo_hl,kso_hl,sigkap,    &
                      ddts,minsfc,cp0,wrtemp,wrtrho,use_ocean,   &
                      imax,kl,ol)

  !=============================================================================
  ! End: update_coupled
  !=============================================================================
  
    
  do i = 1,imax
    ustar_ave(i) = ustar_ave(i) + ustar(i)/real(mcount)
    fg_ave(i) = fg_ave(i) + fg(i)/real(mcount)
  end do  

  
  ! Account for phase transistions
  select case(tcalmeth)
      
    case(0) ! correct saturated air
      do k = 1,kl
        ! Check for -ve values  
        do i = 1,imax  
          qt(i) = max( qfg(i,k) + qlg(i,k) + qvg(i,k), 0. )
          qc(i) = max( qfg(i,k) + qlg(i,k), 0. ) 
          qfg(i,k) = max( qfg(i,k), 0. )
          qlg(i,k) = max( qlg(i,k), 0. )
   
          ! account for saturation  
          templ(i) = thetal(i,k)/sigkap(k)
          temp(i) = templ(i) + (lv*qlg(i,k)+ls*qfg(i,k))/cp
          pres(i) = ps(i)*sig(k)
          if ( qfg(i,k)+qlg(i,k)>0. ) then
            fice(i) = qfg(i,k)/(qfg(i,k)+qlg(i,k))
          else
            fice(i) = 0.
          end if
        end do  
         
        qsat(:) = getqsat(templ(:),pres(:),fice)
          
        do i = 1,imax
          lx = lv + lf*fice(i)
          dqsdt = qsat(i)*lx/(rv*templ(i)**2)
          al = cp/(cp+lx*dqsdt)
          qc(i) = max( al*(qt(i) - qsat(i)), qc(i) )
          if ( temp(i)>=tice ) then
            qfg(i,k) = max( fice(i)*qc(i), 0. )  
            qlg(i,k) = max( qc(i) - qfg(i,k), 0. )
          end if
   
          qvg(i,k) = max( qt(i) - qfg(i,k) - qlg(i,k), 0. )
          if ( qlg(i,k)+qfg(i,k)>1.e-8 ) then
            cf(i,k) = max( cf(i,k), 1.e-10 )
          end if
        end do  
      end do  

    case(1) ! fix rounding errors only
      do k = 1,kl
        ! Check for -ve values  
        do i = 1,imax
          qt(i) = max( qfg(i,k) + qlg(i,k) + qvg(i,k), 0. )
          qfg(i,k) = max( qfg(i,k), 0. )
          qlg(i,k) = max( qlg(i,k), 0. )
          qvg(i,k) = max( qt(i) - qfg(i,k) - qlg(i,k), 0. )
          if ( qlg(i,k)+qfg(i,k)>1.e-8 ) then
            cf(i,k) = max( cf(i,k), 1.e-10 )
          end if
        end do
      end do  
      
    case(2) ! correct saturated air (below pbl)
      do k = 1,kl
        do i = 1,imax  
          templ(i) = thetal(i,k)/sigkap(k)
          pres(i) = ps(i)*sig(k)
          if ( qfg(i,k)+qlg(i,k)>0. ) then
            fice(i) = qfg(i,k)/(qfg(i,k)+qlg(i,k))
          else
            fice(i) = 0.
          end if
        end do  

        qsat(:) = getqsat(templ(:),pres(:),fice)
          
        do i = 1,imax
          if ( zz(i,k)<zi(i) ) then
        
            ! Check for -ve values  
            qt(i) = max( qfg(i,k) + qlg(i,k) + qvg(i,k), 0. )
            qc(i) = max( qfg(i,k) + qlg(i,k), 0. )
        
            qfg(i,k) = max( qfg(i,k), 0. )
            qlg(i,k) = max( qlg(i,k), 0. )
    
            ! account for saturation  
            temp(i) = thetal(i,k)/sigkap(k) + (lv*qlg(i,k)+ls*qfg(i,k))/cp
            lx = lv + lf*fice(i)
            dqsdt = qsat(i)*lx/(rv*templ(i)**2)
            al = cp/(cp+lx*dqsdt)
            qc(i) = max( al*(qt(i) - qsat(i)), qc(i) )
            if ( temp(i)>=tice ) then
              qfg(i,k) = max( fice(i)*qc(i), 0. )  
              qlg(i,k) = max( qc(i) - qfg(i,k), 0. )
            end if
   
            qvg(i,k) = max( qt(i) - qfg(i,k) - qlg(i,k), 0. )
            if ( qlg(i,k)+qfg(i,k)>1.e-8 ) then
              cf(i,k) = max( cf(i,k), 1.e-10 )
            end if
            
          end if ! zz<zi
        end do   ! i
      end do     ! k
      
  end select

end do ! kcount loop
    

! update theta (output) from thetal (conserved)
do k = 1,kl
  do i = 1,imax  
    theta(i,k) = thetal(i,k) + sigkap(k)*(lv*qlg(i,k)+ls*qfg(i,k))/cp
  end do  
end do  

if ( use_ocean ) then
  fg(:) = fg_ave(:)
end if ! use_ocean
  

call mlocheck("Coupled-mixing",water_temp=w_t,water_sal=w_s,water_u=w_u, &
              water_v=w_v,ice_u=i_u,ice_v=i_v)

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plumerise
                  
subroutine plumerise(thetal,qvg,qlg,qfg,ni,cf,ua,va,     &
                     zz,dz_hl,tke,theta,thetav,ps,       &
                     mask,zi,ustar,wstar,wt0,wq0,wtv0,   &
                     sig,sigkap,                         &
                     mflx,tlup,qvup,qlup,qfup,niup,cfup, &
                     cm12,imax,kl)

implicit none

integer, intent(in) :: imax,kl
integer i, k
real, dimension(imax,kl), intent(in) :: thetal, qvg, qlg, qfg, ni, cf
real, dimension(imax,kl), intent(in) :: ua, va
real, dimension(imax,kl), intent(in) :: zz, dz_hl, theta
real, dimension(kl), intent(in) :: sig, sigkap
real, dimension(imax,kl), intent(inout) :: tke, thetav
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: ustar, wstar, wt0, wq0, wtv0
real, dimension(imax), intent(inout) :: zi
real, dimension(imax,kl), intent(out) :: mflx, tlup, qvup, qlup, qfup, niup, cfup
real, intent(in) :: cm12
logical, dimension(imax), intent(in) :: mask

real, dimension(imax,kl) :: nn, w2up, cxup, rino
real, dimension(imax) :: qupsat, templ, pres, fice
real ent, dzht, qtup, qxup, qcup, thup, tvup, vvk
real as, bs, cs, xp, upf
real lx, dqsdt, al
real, parameter :: fac = 10. ! originally fac=100.
real, parameter :: ricr = 0.3

! Initialise updraft
do k = 1,kl
  do i = 1,imax  
    mflx(i,k) = 0.  
    tlup(i,k) = thetal(i,k)
    qvup(i,k) = qvg(i,k)
    qlup(i,k) = qlg(i,k)
    qfup(i,k) = qfg(i,k)
    niup(i,k) = ni(i,k)
    cfup(i,k) = cf(i,k)
    w2up(i,k) = 0.
    nn(i,k) = 0.
  end do  
end do

do i = 1,imax
  tke(i,1) = cm12*ustar(i)**2 + ce3*wstar(i)**2
  tke(i,1) = max(tke(i,1), mintke)
end do


! first level -----------------

do i = 1,imax
    
  ! Entrainment rates  
  !ent(:) = entfn(zz(:,1), zi(:), imax)
  ent = max( ent0/max( zz(i,1), 1. ) + ent1/max( zi(i)-zz(i,1), ezmin ), ent_min )

  ! initial thermodynamic state
  ! split qtot into components (conservation of thetal and qtot is maintained)
  if ( mask(i) ) then
    tlup(i,1) = thetal(i,1) + be*wt0(i)/sqrt(tke(i,1)) ! Hurley 2007
    qvup(i,1) = qvg(i,1)    + be*wq0(i)/sqrt(tke(i,1)) ! Hurley 2007
    qlup(i,1) = qlg(i,1)
    qfup(i,1) = qfg(i,1)
    niup(i,1) = ni(i,1)
    cfup(i,1) = cf(i,1)
  end if
  qvup(i,1) = max( qvup(i,1), 0. ) ! in case of bad land surface flux
  ! update updraft velocity and mass flux
  thetav(i,1) = theta(i,1)*(1.+0.61*qvg(i,1)-qlg(i,1)-qfg(i,1))
  nn(i,1) = grav*be*wtv0(i)/(thetav(i,1)*sqrt(tke(i,1))) ! Hurley 2007
  dzht = zz(i,1)
  if ( mask(i) ) then
    w2up(i,1) = 2.*dzht*b2*nn(i,1)/(1.+2.*dzht*b1*ent)    ! Hurley 2007
  end if
  cxup(i,1) = 0.
  rino(i,1) = 0.
  
end do ! i loop

! updraft with condensation
do k = 2,kl

  do i = 1,imax  
    templ(i) = tlup(i,k)/sigkap(k)     ! templ,up
    pres(i) = ps(i)*sig(k)
    if ( qfg(i,k)+qlg(i,k)>0. ) then
      fice(i) = qfg(i,k)/(qfg(i,k)+qlg(i,k))
    else
      fice(i) = 0.
    end if 
  end do ! i loop
  
  qupsat(:) = getqsat(templ(:),pres(:),fice(:))
  
  do i = 1,imax
    ! Entrainment rates
    !ent(:) = entfn(zz(:,k), zi(:), imax)
    ent = max( ent0/max( zz(i,k), 1. ) + ent1/max( zi(i)-zz(i,k), ezmin ), ent_min )  
    dzht = dz_hl(i,k-1)
    if ( w2up(i,k-1)>0. .and. mask(i) ) then
      ! entrain air into plume
      ! split qtot into components (conservation of qtot is maintained)
      tlup(i,k) = (tlup(i,k-1)+dzht*ent*thetal(i,k))/(1.+dzht*ent)
      qvup(i,k) = (qvup(i,k-1)+dzht*ent*qvg(i,k)   )/(1.+dzht*ent)
      qlup(i,k) = (qlup(i,k-1)+dzht*ent*qlg(i,k)   )/(1.+dzht*ent)
      qfup(i,k) = (qfup(i,k-1)+dzht*ent*qfg(i,k)   )/(1.+dzht*ent)
      niup(i,k) = (niup(i,k-1)+dzht*ent*ni(i,k)    )/(1.+dzht*ent)
      cfup(i,k) = (cfup(i,k-1)+dzht*ent*cf(i,k)    )/(1.+dzht*ent)
    end if
    ! calculate conserved variables
    thetav(i,k) = theta(i,k)*(1.+0.61*qvg(i,k)-qlg(i,k)-qfg(i,k))
    qtup = qvup(i,k) + qlup(i,k) + qfup(i,k)    ! qtot,up
    if ( qtup>qupsat(i) .and. w2up(i,k-1)>0. ) then
      qxup = qupsat(i)
      cxup(i,k) = 1.
    else
      qxup = qtup
      cxup(i,k) = 0.
    end if
    lx = lv + lf*fice(i)
    dqsdt = qupsat(i)*lx/(rv*templ(i)**2)
    al = cp/(cp+lx*dqsdt)
    qcup = max(al*(qtup-qxup), 0.)             ! qcondensate,up after redistribution
    qcup = min(qcup, qcmf)                     ! limit condensation with simple autoconversion
    thup = tlup(i,k) + sigkap(k)*qcup*lx/cp    ! theta,up after redistribution
    tvup = thup + theta(i,k)*(0.61*qxup-qcup)  ! thetav,up after redistribution
    if ( w2up(i,k-1)>0. .and. mask(i) ) then
      nn(i,k) = grav*(tvup-thetav(i,k))/thetav(i,k) ! calculate buayancy
      w2up(i,k) = (w2up(i,k-1)+2.*dzht*b2*nn(i,k))/(1.+2.*dzht*b1*ent) ! update updraft velocity
    else
      nn(i,k) = 0.  
      w2up(i,k) = 0.
    end if
    vvk = (ua(i,k)-ua(i,1))**2 + (va(i,k)-va(i,1))**2 + fac*ustar(i)**2
    rino(i,k) = grav*(thetav(i,k)-thetav(i,1))*(zz(i,k)-zz(i,1))/max(thetav(i,k)*vvk,1.e-20)
    ! test if maximum plume height is reached
    if ( w2up(i,k)<=0. .and. w2up(i,k-1)>0. .and. mask(i) ) then ! unstable
      as = 2.*b2*(nn(i,k)-nn(i,k-1))/dzht
      bs = 2.*b2*nn(i,k-1)
      cs = w2up(i,k-1)
      xp = -2.*cs/(bs-sqrt(max(bs**2-4.*as*cs,0.)))
      xp = min(max(xp,0.),dzht)
      zi(i) = xp + zz(i,k-1)
    else if ( rino(i,k)>ricr .and. rino(i,k-1)<=ricr .and. .not.mask(i) ) then ! stable
      xp = (ricr-rino(i,k-1))/(rino(i,k)-rino(i,k-1))
      xp = min( max(xp, 0.), 1.)
      zi(i) = zz(i,k-1) + xp*(zz(i,k)-zz(i,k-1))
    end if
  end do ! i loop
end do   ! k loop
  
! update mass flux
do i = 1,imax
  mflx(i,1) = m0*sqrt(max(w2up(i,1), 0.))
end do  
  
do k = 2,kl
  do i = 1,imax
    dzht = dz_hl(i,k-1)
    mflx(i,k) = (1.-cxup(i,k))*m0*sqrt(max(w2up(i,k), 0.))         &
              + cxup(i,k)*mflx(i,k-1)/(1.+dzht*(dtrc0-entc0))
    upf = mflx(i,k-1)/sqrt(max(w2up(i,k-1), 1.e-8))
    mflx(i,k) = min( mflx(i,k), upf*sqrt(max(w2up(i,k), 0.)) )
  end do
end do    

return
end subroutine plumerise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate buoyancy

subroutine calculate_buoyancy(ppb,thetal,qvg,qlg,qfg,cf,    &
                              theta,thetav,km,zz,zzh,fzzh,  &
                              sigkap,imax,kl)

implicit none

integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(in) :: thetal, qvg, qlg, qfg, cf
real, dimension(imax,kl), intent(in) :: theta, km, zz, zzh
real, dimension(imax,kl-1), intent(in) :: fzzh
real, dimension(imax,kl), intent(inout) :: thetav
real, dimension(imax,kl), intent(out) :: ppb
real, dimension(kl), intent(in) :: sigkap

integer i, k
real, dimension(imax,kl) :: thetalhl, thetavhl, qthl
real tempt, tempv, rvar, fc, dc, mc, bvf, qtot, tcc
real dz_fl


! set top boundary condition for TKE-eps source terms
do i = 1,imax
  ppb(i,kl) = 0.
end do

! Calculate buoyancy term
select case(buoymeth)
  !case(0) ! Blend staturated and unsaturated terms - saturated method from Durran and Klemp JAS 1982 (see also WRF)
  !  qsatc(:,1:kl) = qvg(:,1:kl)                      ! assume qvg is saturated inside cloud
  !  ff(:,1:kl) = qfg(:,1:kl)/max(cf(:,1:kl),1.e-10)  ! inside cloud value assuming max overlap
  !  gg(:,1:kl) = qlg(:,1:kl)/max(cf(:,1:kl),1.e-10)  ! inside cloud value assuming max overlap
  !  do k = 1,kl
  !    do i = 1,imax
  !      tbb = max(1.-cf(i,k),1.e-10)
  !      qgnc(i,k) = (qvg(i,k)-(1.-tbb)*qsatc(i,k))/tbb  ! outside cloud value
  !      qgnc(i,k) = min(max(qgnc(i,k),0.),qsatc(i,k))
  !    end do
  !  end do
  !  do k = 1,kl-1
  !    thetalhl(:,k) = thetal(:,k) + fzzh(:,k)*(thetal(:,k+1)-thetal(:,k))
  !    quhl(:,k) = qgnc(:,k) + fzzh(:,k)*(qgnc(:,k+1)-qgnc(:,k))
  !    qshl(:,k) = qsatc(:,k) + fzzh(:,k)*(qsatc(:,k+1)-qsatc(:,k))
  !    qlhl(:,k) = gg(:,k) + fzzh(:,k)*(gg(:,k+1)-gg(:,k))
  !    qfhl(:,k) = ff(:,k) + fzzh(:,k)*(ff(:,k+1)-ff(:,k))
  !  end do
  !  do k = 2,kl-1
  !    ! saturated
  !    do i = 1,imax
  !      dz_fl = zzh_t(i,k,tile) - zzh_t(i,k-1,tile)
  !      thetac = thetal(i,k)+sigkap(k)*(lv*gg(i,k)+ls*ff(i,k))/cp   ! inside cloud value
  !      tempc = thetac/sigkap(k)                                    ! inside cloud value          
  !      tqq = (1.+lv*qsatc(i,k)/(rd*tempc))/(1.+lv*lv*qsatc(i,k)/(cp*rv*tempc**2))
  !      tbb = -grav*km(i,k)*(tqq*((thetalhl(i,k)-thetalhl(i,k-1)+sigkap(k)/cp*(lv*(qlhl(i,k)-qlhl(i,k-1))  &
  !            +ls*(qfhl(i,k)-qfhl(i,k-1))))/thetac+lv/cp*(qshl(i,k)-qshl(i,k-1))/tempc)                    &
  !            -qshl(i,k)-qlhl(i,k)-qfhl(i,k)+qshl(i,k-1)+qlhl(i,k-1)+qfhl(i,k-1))/dz_fl(i,k)
  !      ! unsaturated
  !      tcc = -grav*km(i,k)*(thetalhl(i,k)-thetalhl(i,k-1)+thetal(i,k)*0.61*(quhl(i,k)-quhl(i,k-1)))  &
  !                       /(thetal(i,k)*dz_fl)
  !      ppb(i,k) = (1.-cf(i,k))*tcc+cf(i,k)*tbb ! cloud fraction weighted
  !    end do
  !  end do
  !  ! saturated
  !  do i = 1,imax
  !    thetac = thetal(i,1)+sigkap(1)*(lv*gg(i,1)+ls*ff(i,1))/cp        ! inside cloud value
  !    tempc = thetac/sigkap(1)                                         ! inside cloud value          
  !    tqq = (1.+lv*qsatc(i,1)/(rd*tempc))/(1.+lv*lv*qsatc(i,1)/(cp*rv*tempc*tempc))
  !    tbb = -grav*km(i,1)*(tqq*((thetalhl(i,1)-thetal(i,1)+sigkap(1)/cp*(lv*(qlhl(i,1)-qlg(i,1))  &
  !          +ls*(qfhl(i,1)-qfg(i,1))))/thetac+lv/cp*(qshl(i,1)-qsatc(i,1))/tempc)                 &
  !          -qshl(i,1)-qlhl(i,1)-qfhl(i,1)+qsatc(i,1)+qlg(i,1)+qfg(i,1))/(zzh(i,1)-zz(i,1))     ! unsaturated
  !    tcc = -grav*km(i,1)*(thetalhl(i,1)-thetal(i,1)+thetal(i,1)*0.61*(quhl(i,1)-qgnc(i,1)))      &
  !                   /(thetal(i,1)*(zzh(i,1)-zz(i,1)))
  !    ppb(i,1) = (1.-cf(i,1))*tcc+cf(i,1)*tbb ! cloud fraction weighted
  !  end do

  case(1) ! Marquet and Geleyn QJRMS (2012) for partially saturated
    do k = 1,kl-1
      do i = 1,imax
        thetalhl(i,k) = thetal(i,k) + fzzh(i,k)*(thetal(i,k+1)-thetal(i,k))
        ! calculate qtot at half levels
        qthl(i,k) = qvg(i,k) + qlg(i,k) + qfg(i,k)              &
                  + fzzh(i,k)*(qvg(i,k+1)+qlg(i,k+1)+qfg(i,k+1) &
                              -qvg(i,k)-qlg(i,k)-qfg(i,k))
        thetav(i,k) = theta(i,k)*(1.+0.61*qvg(i,k)-qlg(i,k)-qfg(i,k)) 
      end do  
    end do
    do k = 2,kl-1
      do i = 1,imax
        dz_fl = zzh(i,k) - zzh(i,k-1)
        tempt = thetal(i,k)/sigkap(k) + (lv*qlg(i,k)+ls*qfg(i,k))/cp
        tempv = thetav(i,k)/sigkap(k)
        rvar = rd*tempv/tempt ! rvar = qd*rd+qv*rv
        fc = (1.-cf(i,k))+cf(i,k)*(lv*rvar/(cp*rv*tempt))
        dc = (1.+0.61*qvg(i,k))*lv*qvg(i,k)/(rd*tempv)
        mc = (1.+dc)/(1.+lv*qlg(i,k)/(cp*tempt)+dc*fc)
        bvf = grav*mc*(thetalhl(i,k)-thetalhl(i,k-1))/(thetal(i,k)*dz_fl)     &
             +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(i,k)-qthl(i,k-1))/dz_fl
        ppb(i,k) = -km(i,k)*bvf
      end do
    end do
    do i = 1,imax
      tempt = thetal(i,1)/sigkap(1) + (lv*qlg(i,1)+ls*qfg(i,1))/cp
      tempv = thetav(i,1)/sigkap(1)
      rvar = rd*tempv/tempt ! rvar = qd*rd+qv*rv
      fc = (1.-cf(i,1))+cf(i,1)*(lv*rvar/(cp*rv*tempt))
      dc = (1.+0.61*qvg(i,1))*lv*qvg(i,1)/(rd*tempv)
      mc = (1.+dc)/(1.+lv*qlg(i,1)/(cp*tempt)+dc*fc)
      qtot = qvg(i,1) + qlg(i,1) + qfg(i,1)
      bvf = grav*mc*(thetalhl(i,1)-thetal(i,1))/(thetal(i,1)*(zzh(i,1)-zz(i,1))) &
           +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(i,1)-qtot)/(zzh(i,1)-zz(i,1))
      ppb(i,1) = -km(i,1)*bvf
    end do
      
  case(2) ! dry convection from Hurley 2007
    do k = 1,kl
      do i = 1,imax
        thetav(i,k) = theta(i,k)*(1.+0.61*qvg(i,k)-qlg(i,k)-qfg(i,k))
      end do  
    end do
    do k = 1,kl-1
      do i = 1,imax  
        thetavhl(i,k) = thetav(i,k) + fzzh(i,k)*(thetav(i,k+1)-thetav(i,k))
      end do
    end do
    do k = 2,kl-1
      do i = 1,imax
        dz_fl = zzh(i,k) - zzh(i,k-1)
        tcc = -grav*km(i,k)*(thetavhl(i,k)-thetavhl(i,k-1))/(thetav(i,k)*dz_fl)
        ppb(i,k) = tcc
      end do
    end do
    do i = 1,imax
      tcc = -grav*km(i,1)*(thetavhl(i,1)-thetav(i,1))/(thetav(i,1)*(zzh(i,1)-zz(i,1)))
      ppb(i,1) = tcc 
    end do

end select

return
end subroutine calculate_buoyancy
                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! update_tkeeps
         
subroutine update_tkeeps(tke,eps,pps,ppb,wstar,zi,km,zz,dz_hl,idzm,idzp,fzzh, &
                         ddts,cm34,imax,kl)

implicit none

integer, intent(in) :: imax, kl
real, intent(in) :: ddts, cm34
real, dimension(imax,kl), intent(inout) :: tke, eps
real, dimension(imax,kl), intent(in) :: km, idzm, zz
real, dimension(imax,kl-1), intent(in) :: idzp, fzzh, dz_hl
real, dimension(imax,kl), intent(in) :: pps, ppb
real, dimension(imax), intent(in) :: wstar, zi

integer i, k
real, dimension(imax,kl) :: kmo, rr
real, dimension(imax,2:kl) :: qq
real, dimension(imax,kl-1) :: ccd    ! working
real, dimension(imax,kl) :: bbd, ddd ! working
real, dimension(imax,2:kl) :: aad    ! working
real ppt, tbb, tff

do k = 1,kl-1
  do i = 1,imax
    kmo(i,k) = km(i,k) + fzzh(i,k)*(km(i,k+1)-km(i,k))
  end do
end do  

! top boundary condition to avoid unphysical behaviour at the top of the model
do i = 1,imax
  tke(i,kl) = mintke
  eps(i,kl) = mineps
end do  
  
! Pre-calculate eddy diffusivity mixing terms
! -ve because gradient is calculated at t+1
do k = 2,kl-1
  do i = 1,imax
    qq(i,k) = -ddts*idzm(i,k)/dz_hl(i,k-1)
    rr(i,k) = -ddts*idzp(i,k)/dz_hl(i,k)
  end do
end do
  
! eps vertical mixing
do k = 2,kl-1
  do i = 1,imax
    ! Calculate transport source term on full levels
    ppt = kmo(i,k)*idzp(i,k)*(tke(i,k+1)-tke(i,k))/dz_hl(i,k)  &
         -kmo(i,k-1)*idzm(i,k)*(tke(i,k)-tke(i,k-1))/dz_hl(i,k-1)
    
    aad(i,k) = ce0*kmo(i,k-1)*qq(i,k)
    ccd(i,k) = ce0*kmo(i,k)*rr(i,k)
    ! follow Hurley 2007 to make scheme more numerically stable
    bbd(i,k) = 1.-aad(i,k)-ccd(i,k)+ddts*ce2*eps(i,k)/tke(i,k)
    ddd(i,k) = eps(i,k)+ddts*eps(i,k)/tke(i,k)  &
              *ce1*(pps(i,k)+max(ppb(i,k),0.)+max(ppt,0.))
  end do  
end do
do i = 1,imax
  ddd(i,2)    = ddd(i,2)    - aad(i,2)*eps(i,1)
  ddd(i,kl-1) = ddd(i,kl-1) - ccd(i,kl-1)*mineps
end do  
call thomas(eps(:,2:kl-1),aad(:,3:kl-1),bbd(:,2:kl-1),ccd(:,2:kl-2),ddd(:,2:kl-1))
  
! TKE vertical mixing
do k = 2,kl-1
  do i = 1,imax
    aad(i,k) = kmo(i,k-1)*qq(i,k)
    ccd(i,k) = kmo(i,k)*rr(i,k)
    bbd(i,k) = 1. - aad(i,k) - ccd(i,k)
    ddd(i,k) = tke(i,k) + ddts*(pps(i,k)+ppb(i,k)-eps(i,k))
  end do
end do
do i = 1,imax
  ddd(i,2)    = ddd(i,2)    - aad(i,2)*tke(i,1)
  ddd(i,kl-1) = ddd(i,kl-1) - ccd(i,kl-1)*mintke
end do
call thomas(tke(:,2:kl-1),aad(:,3:kl-1),bbd(:,2:kl-1),ccd(:,2:kl-2),ddd(:,2:kl-1))
  
! limit decay of TKE and EPS with coupling to mass flux term
if ( tkemeth==1 ) then
  do k = 2,kl-1
    do i = 1,imax
      if ( wstar(i)>0.5 .and. zz(i,k)>0.5*zi(i) .and. zz(i,k)<0.95*zi(i) ) then
        tbb = max(1.-0.05*dz_hl(i,k-1)/250.,0.)        
        tke(i,k) = max( tke(i,k), tbb*tke(i,k-1) )
        eps(i,k) = max( eps(i,k), tbb*eps(i,k-1) )
      end if
    end do
  end do
end if

! apply limits and corrections to tke and eps terms
do k = 2,kl-1
  do i = 1,imax
    tke(i,k) = max(tke(i,k),mintke)
    tff = cm34*tke(i,k)*sqrt(tke(i,k))
    eps(i,k) = min(eps(i,k),tff/minl)
    eps(i,k) = max(eps(i,k),tff/maxl,mineps)
  end do
end do

return
end subroutine update_tkeeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update ocean buoyancy
                         
subroutine mlo_update_n2(n2_o,w_t,w_s,deptho_depth,deptho_depth_hl,deptho_dz, &
                         fdeptho_hl,wrtrho,imax,ol)

use mlo_ctrl, only : calcdensity

implicit none

integer, intent(in) :: imax, ol
real, intent(in) :: wrtrho
real, dimension(imax,ol), intent(out) :: n2_o
real, dimension(imax,ol), intent(in) :: w_t, w_s
real, dimension(imax,ol), intent(in) :: deptho_depth
real, dimension(imax,ol+1), intent(in) :: deptho_depth_hl
real, dimension(imax,ol), intent(in) :: deptho_dz
real, dimension(imax,ol), intent(in) :: fdeptho_hl

integer i, k
real, dimension(imax,ol) :: rho_o, alpha_o, beta_o
real, dimension(imax,ol) :: rho_hl_o
real, dimension(imax) :: pxtr_o

! calculate density
do i = 1,imax
  pxtr_o(i) = 0. ! neglect surface pressure as rho is only used for n2 calculation
end do  
call calcdensity(rho_o,alpha_o,beta_o,w_t,w_s,deptho_dz,pxtr_o)

do i = 1,imax
  rho_hl_o(i,1) = 0.  
end do
do k = 2,ol
  do i = 1,imax  
    rho_hl_o(i,k) = rho_o(i,k-1) + fdeptho_hl(i,k)*(rho_o(i,k)-rho_o(i,k-1))
  end do
end do
      
!n2 (full levels)
do k = 1,ol
  do i = 1,imax  
    n2_o(i,k) = 0.
  end do
end do

do k = 2,ol-1
  do i = 1,imax
    if ( deptho_dz(i,k)>1.e-4 ) then
      n2_o(i,k) = -grav/wrtrho*(rho_hl_o(i,k)-rho_hl_o(i,k+1))/deptho_dz(i,k)
    end if
  end do  
end do

do i = 1,imax
  if ( deptho_dz(i,1)>1.e-4 ) then
    n2_o(i,1) = -grav/wrtrho*(rho_o(i,1)-rho_hl_o(i,2))/(deptho_depth_hl(i,2)-deptho_depth(i,1))
  end if
  if ( deptho_dz(i,ol)>1.e-4 ) then
    n2_o(i,ol) = -grav/wrtrho*(rho_hl_o(i,ol-1)-rho_o(i,ol))/(deptho_depth(i,ol)-deptho_depth_hl(i,ol-1))
  end if  
end do

return
end subroutine mlo_update_n2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate buoyancy and shear for mlo

subroutine mlo_update_buoyancy(pps_o,ppb_o,km_o,ks_o,w_u,w_v,n2_o,dwdx_o,dwdy_o, &
                               deptho_dz,fdeptho_hl,nops,nopb,imax,ol)

implicit none

integer, intent(in) :: imax, ol
integer, intent(in) :: nops, nopb
real, dimension(imax,ol), intent(out) :: pps_o, ppb_o
real, dimension(imax,ol), intent(in) :: km_o, ks_o
real, dimension(imax,ol), intent(in) :: w_u, w_v
real, dimension(imax,ol), intent(in) :: n2_o
real, dimension(imax,ol), intent(in) :: dwdx_o, dwdy_o
real, dimension(imax,ol), intent(in) :: deptho_dz
real, dimension(imax,ol), intent(in) :: fdeptho_hl

integer i, k
real, dimension(imax,ol) :: u_hl_o, v_hl_o
real dudz, dvdz, shear_o

!shear production
do k = 2,ol-1
  do i = 1,imax
    pps_o(i,k) = 0.
  end do  
end do

if ( nops==0 ) then    
  do k = 2,ol
    do i = 1,imax  
      u_hl_o(i,k) = w_u(i,k-1) + fdeptho_hl(i,k)*(w_u(i,k)-w_u(i,k-1))
      v_hl_o(i,k) = w_v(i,k-1) + fdeptho_hl(i,k)*(w_v(i,k)-w_v(i,k-1))
    end do
  end do
  do k = 2,ol-1
    do i = 1,imax
      dudz = 0.
      dvdz = 0.
      if ( deptho_dz(i,k)>1.e-4 ) then  
        dudz = (u_hl_o(i,k+1)-u_hl_o(i,k))/deptho_dz(i,k)
        dvdz = (v_hl_o(i,k+1)-v_hl_o(i,k))/deptho_dz(i,k)
      end if  
      shear_o = (dudz+dwdx_o(i,k))**2 &
              + (dvdz+dwdy_o(i,k))**2
      pps_o(i,k) = km_o(i,k)*shear_o
    end do  
  end do
end if

!buoyancy production
do k = 2,ol-1
  do i = 1,imax
    ppb_o(i,k) = 0.
  end do  
end do

if ( nopb==0 ) then
  do k = 2,ol-1
    do i = 1,imax
      ppb_o(i,k) = -ks_o(i,k)*n2_o(i,k)
    end do  
  end do
end if

return
end subroutine mlo_update_buoyancy
                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update ocean tke and eps 
                         
subroutine mlo_update_keps(kmo_hl,kso_hl,w_tke,w_eps,pps_o,ppb_o,deptho_depth, &
                           deptho_depth_hl,deptho_dz,deptho_dz_hl,ibot,        &
                           ddts,fixedce3,eps_mode,k_mode,imax,ol)

implicit none

integer, intent(in) :: imax, ol
integer, intent(in) :: fixedce3, eps_mode, k_mode
integer, dimension(imax), intent(in) :: ibot
real, intent(in) :: ddts
real, dimension(imax,ol), intent(in) :: kmo_hl, kso_hl
real, dimension(imax,ol), intent(inout) :: w_tke, w_eps
real, dimension(imax,ol), intent(in) :: deptho_depth
real, dimension(imax,ol+1), intent(in) :: deptho_depth_hl
real, dimension(imax,ol), intent(in) :: deptho_dz
real, dimension(imax,2:ol), intent(in) :: deptho_dz_hl
real, dimension(imax,ol), intent(in) :: pps_o, ppb_o

integer i, k
real, dimension(imax,ol-1) :: ccd    ! working
real, dimension(imax,ol) :: bbd, ddd ! working
real, dimension(imax,2:ol) :: aad    ! working
real, dimension(imax,ol) :: ce3_o
real, parameter :: ce1 = 1.44        ! eps production coefficient
real, parameter :: ce2 = 1.92        ! eps sink coefficient
real, parameter :: ce3stable = -0.4  ! eps buoyancy coefficient for stable stratification
real, parameter :: ce3unstable = 1.0 ! eps buoyancy coefficient for unstable stratification
real, parameter :: sigmaeps = 1.08   ! eps Schmidt number

!calculate ce3
if ( fixedce3==1 ) then
  do k = 2,ol-1
    do i = 1,imax
      ce3_o(i,k) = ce3stable
    end do
  end do
else if ( fixedce3==0 ) then
  do k = 2,ol-1
    do i = 1,imax
      if ( ppb_o(i,k) < 0. ) then
        ce3_o(i,k) = ce3unstable
      else
        ce3_o(i,k) = ce3stable
      end if
    end do  
  end do
end if

!solve eps
!setup diagonals
do k = 2,ol
  do i = 1,imax
    aad(i,k) = 0.
  end do
end do
do k = 1,ol-1
  do i = 1,imax
    ccd(i,k) = 0.
  end do
end do
do k = 1,ol
  do i = 1,imax
    bbd(i,k) = 1.
    ddd(i,k) = w_eps(i,k)
  end do
end do
do k = 2,ol-1
  do i = 1,imax
    if ( deptho_dz(i,k)*deptho_dz(i,k-1)>1.e-4 .and. ibot(i)>k ) then
      aad(i,k) = -ddts*kmo_hl(i,k)/(deptho_dz(i,k)*deptho_dz_hl(i,k)*sigmaeps)
    end if
    if ( deptho_dz(i,k)*deptho_dz(i,k+1)>1.e-4 .and. ibot(i)>k ) then
      ccd(i,k) = -ddts*kmo_hl(i,k+1)/(deptho_dz(i,k)*deptho_dz_hl(i,k+1)*sigmaeps)
    end if
    bbd(i,k) = 1. - aad(i,k) - ccd(i,k)
  end do  
end do
if ( eps_mode==0 ) then !explicit eps
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        ddd(i,k) = ddd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)*(ce1*pps_o(i,k) &
                            + ce3_o(i,k)*ppb_o(i,k) - ce2*w_eps(i,k))
      end if
    end do  
  end do
else if ( eps_mode==1 ) then !quasi implicit for eps, Patanker (1980)
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts*ce2*w_eps(i,k)/w_tke(i,k)
        ddd(i,k) = ddd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)*(ce1*pps_o(i,k) &
                            + ce3_o(i,k)*ppb_o(i,k))
      end if
    end do  
  end do
else if ( eps_mode==2 ) then !quasi implicit for eps & pb, Patanker (1980) & Burchard et al sect 4 (1998)
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. (ce1*pps_o(i,k)+ce3_o(i,k)*ppb_o(i,k))>0. .and. &
           ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts*ce2*w_eps(i,k)/w_tke(i,k)
        ddd(i,k) = ddd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)*(ce1*pps_o(i,k) &
                            + ce3_o(i,k)*ppb_o(i,k))
      else if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)*(ce2 - ce3_o(i,k)*ppb_o(i,k)/w_eps(i,k))
        ddd(i,k) = ddd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)*(ce1*pps_o(i,k))
      end if
    end do  
  end do  
end if
do i = 1,imax
  ddd(i,2) = ddd(i,2) - aad(i,2)*w_eps(i,1)
  ddd(i,ol-1) = ddd(i,ol-1) - ccd(i,ol-1)*w_eps(i,ol)
end do
  
!solve using thomas algorithm
call thomas(w_eps(:,2:ol-1),aad(:,3:ol-1),bbd(:,2:ol-1),ccd(:,2:ol-2),ddd(:,2:ol-1))
  
!solve k
!setup diagonals
do k = 2,ol
  do i = 1,imax
    aad(i,k) = 0.
  end do
end do
do k = 1,ol-1
  do i = 1,imax
    ccd(i,k) = 0.
  end do
end do
do k = 1,ol
  do i = 1,imax
    bbd(i,k) = 1.
    ddd(i,k) = w_tke(i,k)
  end do
end do

do k = 2,ol-1
  do i = 1,imax
    if ( deptho_dz(i,k)*deptho_dz(i,k-1)>1.e-4 .and. ibot(i)>k ) then
      aad(i,k) = -ddts*kmo_hl(i,k)/(deptho_dz(i,k)*deptho_dz_hl(i,k))
    end if
    if ( deptho_dz(i,k)*deptho_dz(i,k+1)>1.e-4 .and. ibot(i)>k ) then
      ccd(i,k) = -ddts*kmo_hl(i,k+1)/(deptho_dz(i,k)*deptho_dz_hl(i,k+1))
    end if
    bbd(i,k) = 1. - aad(i,k) - ccd(i,k)
  end do  
end do
if ( k_mode==0 ) then !explicit eps
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        ddd(i,k) = ddd(i,k) + ddts*(pps_o(i,k) + ppb_o(i,k) - w_eps(i,k))
      end if
    end do  
  end do    
else if ( k_mode==1 ) then !quasi impliciit for eps, Patanker (1980)
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)
        ddd(i,k) = ddd(i,k) + ddts*(pps_o(i,k) + ppb_o(i,k))
      end if
    end do  
  end do    
else if ( k_mode==2 ) then !quasi implicit for eps & pb, Patanker (1980) & Burchard et al sect 4 (1998)
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)>1.e-4 .and. (pps_o(i,k)+ppb_o(i,k))>0. .and. &
           ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts*w_eps(i,k)/w_tke(i,k)
        ddd(i,k) = ddd(i,k) + ddts*(pps_o(i,k) + ppb_o(i,k))
      else if ( deptho_dz(i,k)>1.e-4 .and. ibot(i)>k ) then
        bbd(i,k) = bbd(i,k) + ddts/w_tke(i,k)*(w_eps(i,k) - ppb_o(i,k))
        ddd(i,k) = ddd(i,k) + ddts*pps_o(i,k)
      end if
    end do  
  end do
end if
do i = 1,imax
  ddd(i,2) = ddd(i,2) - aad(i,2)*w_tke(i,1)
  ddd(i,ol-1) = ddd(i,ol-1) - ccd(i,ol-1)*w_tke(i,ol)
end do  
  
!solve using thomas algorithm
call thomas(w_tke(:,2:ol-1),aad(:,3:ol-1),bbd(:,2:ol-1),ccd(:,2:ol-2),ddd(:,2:ol-1))

return
end subroutine mlo_update_keps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Limit tke and eps

subroutine mlo_limit_tke(w_tke,w_eps,km_o,ks_o,n2_o,omink,omineps,omaxL,ominL, &
                         limitL,fixedstabfunc,cu0,imax,ol)

implicit none

integer, intent(in) :: imax, ol
integer, intent(in) :: limitL, fixedstabfunc
real, intent(in) :: omink, omineps, omaxL, ominL, cu0
real, dimension(imax,ol), intent(inout) :: w_tke, w_eps
real, dimension(imax,ol), intent(out) :: km_o, ks_o
real, dimension(imax,ol), intent(in) :: n2_o

integer i, k
real minL, alpha, cu, cud
real, dimension(imax,ol) :: L_o

!limit k & eps
do k = 1,ol
  do i = 1,imax
    w_tke(i,k) = max( w_tke(i,k), omink )
    w_eps(i,k) = max( w_eps(i,k), omineps )
  end do
end do

!limit length scale
do k = 1,ol
  do i = 1,imax
    L_o(i,k) = cu0**3*w_tke(i,k)**1.5/w_eps(i,k)
  end do
end do
if ( limitL==1 ) then
  minL = cu0**3*omink**1.5/omineps
  do k = 2,ol-1
    do i = 1,imax
      L_o(i,k) = max( L_o(i,k), minL )
      if ( n2_o(i,k) > 0. ) then
        L_o(i,k) = min( L_o(i,k), sqrt(0.56*w_tke(i,k)/n2_o(i,k)) )
      end if
    end do  
  end do
end if
do k = 1,ol
  do i = 1,imax
    L_o(i,k) = max( min( L_o(i,k), omaxl), ominl )
  end do
end do

!stability functions
if ( fixedstabfunc==1 ) then
  alpha = 0.
  cu = cu0
  cud = 0.6985
  do k = 1,ol
    do i = 1,imax
      km_o(i,k) = max( cu*sqrt(w_tke(i,k))*L_o(i,k), 1.e-6 )
      ks_o(i,k) = max( cud*sqrt(w_tke(i,k))*L_o(i,k), 1.e-6 )
    end do
  end do
else
  do k = 1,ol
    do i = 1,imax
      alpha = L_o(i,k)**2*n2_o(i,k)/w_tke(i,k)
      alpha = max( min( alpha, 0.56), -0.0466 ) ! SHOC (before eq 6.7.5)
      cu = (cu0 + 2.182*alpha)/(1. + 20.4*alpha + 53.12*alpha**2)
      cud = 0.6985/(1. + 17.34*alpha)
      km_o(i,k) = max( cu*sqrt(w_tke(i,k))*L_o(i,k), 1.e-6 )
      ks_o(i,k) = max( cud*sqrt(w_tke(i,k))*L_o(i,k), 1.e-6 )
    end do
  end do
end if

return
end subroutine mlo_limit_tke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update coupled

subroutine update_coupled(thetal,qvg,qlg,qfg,ni,cf,                &
                          ua,va,tlup,qvup,qlup,qfup,niup,          &
                          cfup,fg,eg,rhos,ustar,cduv,tke,eps,      &
                          mflx,fzzh,idzp,idzm,dz_hl,rhoa1,         &
                          dz1,deptho_dz,deptho_dz_hl,              &
                          rad_o,w_t,w_s,w_u,w_v,cd_water,          &
                          cdh_water,cdbot_water,wt0rad_o,          &
                          wt0melt_o,wt0eg_o,icefg_a,wt0fb_o,       &
                          ws0_o,ws0subsurf_o,i_u,i_v,imass,        &
                          fracice,cd_ice,cdbot_ice,ibot,land,      &
                          wu0_o,wv0_o,wt0_o,km_o,ks_o,sigkap,      &
                          ddts,minsfc,cp0,wrtemp,wrtrho,use_ocean, &
                          imax,kl,ol)

implicit none

integer, intent(in) :: imax, kl, ol
integer, dimension(imax), intent(in) :: ibot
integer k, kn, i
real, intent(in) :: minsfc, cp0, wrtemp, wrtrho
real, dimension(imax,kl), intent(inout) :: thetal, qvg, qlg, qfg, ni
real, dimension(imax,kl), intent(inout) :: cf, ua, va
real, dimension(imax,kl), intent(in) :: tke, eps
real, dimension(imax,kl), intent(in) :: tlup, qvup, qlup, qfup
real, dimension(imax,kl), intent(in) :: niup, cfup, mflx
real, dimension(imax,kl), intent(in) :: idzm
real, dimension(imax,kl-1), intent(in) :: fzzh, dz_hl, idzp
real, dimension(imax), intent(inout) :: fg, eg, ustar
real, dimension(kl), intent(in) :: sigkap
real, dimension(imax,ol), intent(inout) :: w_t, w_s, w_u, w_v
real, dimension(imax,ol), intent(in) :: deptho_dz
real, dimension(imax,2:ol), intent(in) :: deptho_dz_hl
real, dimension(imax,ol), intent(in) :: km_o, ks_o
real, dimension(imax,ol), intent(in) :: rad_o
real, dimension(imax), intent(inout) :: wu0_o, wv0_o, wt0_o
real, dimension(imax), intent(in) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o
real, dimension(imax), intent(in) :: wt0eg_o
real, dimension(imax), intent(in) :: rhos, rhoa1, dz1, cduv
real, dimension(imax), intent(in) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o
real, dimension(imax), intent(in) :: fracice, imass, cd_ice, cdbot_ice
real, dimension(imax), intent(inout) :: i_u, i_v
real, intent(in) :: ddts
real, dimension(imax,kl) :: km_a, rr
real, dimension(imax,2:kl) :: qq
real, dimension(imax,2:kl+ol+1) :: aad       ! working
real, dimension(imax,kl+ol+1) :: bbd, ddd    ! working
real, dimension(imax,kl+ol) :: ccd           ! working
real, dimension(imax,kl+ol+1) :: uu, yy, yyu ! working
real, dimension(imax) :: tmp, tmpu
real, dimension(imax) :: fulldeptho, targetdeptho
real, dimension(imax) :: f_ao, f_oa, f_ai, f_ia, f_oi, f_io
real, dimension(imax) :: t1, t3, t4
real deltazo, wq0_a, wt0_a, kmo_a
logical, intent(in) :: use_ocean
logical, dimension(imax), intent(in) :: land


! estimate eddy diffusivity mixing coeff for atmosphere
do k = 1,kl
  do i = 1,imax
    km_a(i,k) = cm0*tke(i,k)**2/eps(i,k)
  end do
end do

! Pre-calculate eddy diffiusivity mixing terms using updated kmo values
! -ve because gradient is calculated at t+1
do k = 1,kl-1
  do i = 1,imax
    kmo_a = km_a(i,k) + fzzh(i,k)*(km_a(i,k+1)-km_a(i,k))
    qq(i,k+1) = -ddts*kmo_a*idzm(i,k+1)/dz_hl(i,k)
    rr(i,k)   = -ddts*kmo_a*idzp(i,k)/dz_hl(i,k)
  end do
end do

! default values for working arrays
do k = 2,kl+ol+1
  do i = 1,imax
    aad(i,k) = 0.
  end do
end do
do k = 1,kl+ol
  do i = 1,imax
    ccd(i,k) = 0.
  end do
end do
do k = 1,kl+ol+1
  do i = 1,imax
    bbd(i,k) = 1.
    ddd(i,k) = 0.
  end do
end do
  
! zero coupling terms for land
do i = 1,imax
  f_ao(i) = 0.
  f_oa(i) = 0.
  f_ai(i) = 0.
  f_ia(i) = 0.
  f_io(i) = 0.
  f_oi(i) = 0.
end do
  



!------------------------------------------------------------------------------
! thetal, thetao - coupled ----

! Atmosphere and ocean (no sea-ice) large tridiagonal matrix

! Original coupling matrix (inverted levels)
! [ bb_a cc_a                               ] [ t_a ]   [ dd_a ]
! [ aa_a bb_a cc_a                          ] [ t_a ]   [ dd_a ]
! [      .... .... ....                     ] [ ... ]   [ .... ]
! [           aa_a bb_a f_ao                ] [ t_a ]   [ dd_a ]
! [                f_oa bb_o cc_o           ] [ t_o ]   [ dd_i ]
! [                     aa_o bb_o cc_o      ] [ t_o ]   [ dd_i ]
! [                          .... .... .... ] [ ... ] = [ .... ]
! [                               aa_o bb_o ] [ t_o ]   [ dd_o ]
    
! k=1 is top of atmosphere and k=kl is bottom of atmosphere
do i = 1,imax
  ccd(i,1) = qq(i,kl) + ddts*mflx(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    aad(i,k) = rr(i,kn)-ddts*mflx(i,kn+1)*fzzh(i,kn)*idzp(i,kn)
    ccd(i,k) = qq(i,kn)+ddts*mflx(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)
  end do  
end do
do i = 1,imax
  aad(i,kl) = rr(i,1) - ddts*mflx(i,2)*fzzh(i,1)*idzp(i,1)
end do

if ( use_ocean ) then
  ! ocean data at kl+k
  do i = 1,imax
    if ( deptho_dz(i,2)*deptho_dz(i,1)>1.e-4 ) then
      ccd(i,kl+1) = -ddts*ks_o(i,2)/(deptho_dz_hl(i,2)*deptho_dz(i,1))
    end if
  end do
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k)*deptho_dz(i,k)>1.e-4 ) then
        aad(i,kl+k) = -ddts*ks_o(i,k)/(deptho_dz_hl(i,k)*deptho_dz(i,k))
      end if
      if ( deptho_dz(i,k+1)*deptho_dz(i,k)>1.e-4 ) then
        ccd(i,kl+k) = -ddts*ks_o(i,k+1)/(deptho_dz_hl(i,k+1)*deptho_dz(i,k))
      end if
    end do
  end do
  do i = 1,imax
    if ( deptho_dz(i,ol)*deptho_dz(i,ol)>1.e-4 ) then
      aad(i,kl+ol) = -ddts*ks_o(i,ol)/(deptho_dz_hl(i,ol)*deptho_dz(i,ol))
    end if
  end do
  
end if
  

do i = 1,imax
  ! k=1 is top of atmosphere and k=kl is bottom of atmosphere
  bbd(i,1) = 1.-qq(i,kl)+ddts*mflx(i,kl)*fzzh(i,kl-1)*idzm(i,kl)
  ddd(i,1) = thetal(i,kl)+ddts*(mflx(i,kl-1)*tlup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)   &
                               +mflx(i,kl)*tlup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    bbd(i,k) = 1.-qq(i,kn)-rr(i,kn)+ddts*(mflx(i,kn)*fzzh(i,kn-1)*idzm(i,kn)     &
                                         -mflx(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn))
    ddd(i,k) = thetal(i,kn)+ddts*(mflx(i,kn-1)*tlup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn) &
                                 +mflx(i,kn)*tlup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)          &
                                 -mflx(i,kn)*tlup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)       &
                                 -mflx(i,kn+1)*tlup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do
end do
do i = 1,imax
  bbd(i,kl) = 1.-rr(i,1)-ddts*mflx(i,1)*(1.-fzzh(i,1))*idzp(i,1)
  ddd(i,kl) = thetal(i,1)-ddts*(mflx(i,1)*tlup(i,1)*(1.-fzzh(i,1))*idzp(i,1) &
                               +mflx(i,2)*tlup(i,2)*fzzh(i,1)*idzp(i,1))
end do  


if ( use_ocean ) then

  do i = 1,imax
    t1(i) = rhos(i)*cdh_water(i)*cp
    if ( land(i) ) then
      wt0_a = fg(i)/(rhos(i)*cp)  ! theta flux
      ddd(i,kl) = ddd(i,kl) + ddts*rhos(i)*wt0_a/(rhoa1(i)*dz1(i))
    else
      bbd(i,kl) = bbd(i,kl)+ddts*t1(i)*(1.-fracice(i))/cp/(rhoa1(i)*dz1(i))
      ddd(i,kl) = ddd(i,kl)-ddts*t1(i)*(1.-fracice(i))/cp*sigkap(1)           &
                           *(lv*qlg(i,1)+ls*qfg(i,1))/cp/(rhoa1(i)*dz1(i))
      ddd(i,kl) = ddd(i,kl)+ddts*t1(i)*(1.-fracice(i))*wrtemp/cp/(rhoa1(i)*dz1(i))
      ddd(i,kl) = ddd(i,kl)+ddts*(fracice(i)*icefg_a(i)/cp)/(rhoa1(i)*dz1(i))
    end if
  end do
    
  ! ocean data at kl+k
  do i = 1,imax
    bbd(i,kl+1) = 1. - ccd(i,kl+1)
    ddd(i,kl+1) = w_t(i,1)
    if ( deptho_dz(i,1)>1.e-4 ) then
      ddd(i,kl+1) = ddd(i,kl+1) - ddts*rad_o(i,1)/deptho_dz(i,1)
    end if
  end do  
  do k = 2,ol-1
    do i = 1,imax
      bbd(i,kl+k) = 1. - aad(i,kl+k) - ccd(i,kl+k)
      ddd(i,kl+k) = w_t(i,k)
      if ( deptho_dz(i,k)>1.e-4 ) then
        ddd(i,kl+k) = ddd(i,kl+k) - ddts*rad_o(i,k)/deptho_dz(i,k)
      end if
    end do  
  end do
  do i = 1,imax
    bbd(i,kl+ol) = 1. - aad(i,kl+ol)
    ddd(i,kl+ol) = w_t(i,ol)
    if ( deptho_dz(i,ol)>1.e-4 ) then
      ddd(i,kl+ol) = ddd(i,kl+ol) - ddts*rad_o(i,ol)/deptho_dz(i,ol)
    end if
    if ( .not.land(i) ) then
      bbd(i,kl+1) = bbd(i,kl+1) + ddts*t1(i)*(1.-fracice(i))/(wrtrho*cp0*deptho_dz(i,1))
      ddd(i,kl+1) = ddd(i,kl+1) - ddts*t1(i)*(1.-fracice(i))*wrtemp/(wrtrho*cp0*deptho_dz(i,1))
      ddd(i,kl+1) = ddd(i,kl+1) - ddts*(wt0rad_o(i)+wt0melt_o(i)+wt0eg_o(i))/deptho_dz(i,1)
      ddd(i,kl+1) = ddd(i,kl+1) + ddts*t1(i)*(1.-fracice(i))/(wrtrho*cp0)*sigkap(1)*(lv*qlg(i,1)+ls*qfg(i,1)) &
                                 /cp/deptho_dz(i,1)
      ddd(i,kl+1) = ddd(i,kl+1) - ddts*wt0fb_o(i)/deptho_dz(i,1)
    end if
  end do  
  
  do i = 1,imax
    if ( .not.land(i) ) then
      f_ao(i) = -ddts*t1(i)*(1.-fracice(i))/cp/(rhoa1(i)*dz1(i))
      f_oa(i) = -ddts*t1(i)*(1.-fracice(i))/(wrtrho*cp0*deptho_dz(i,1))
      ccd(i,kl) = f_ao(i)
      aad(i,kl+1) = f_oa(i)
    end if
  end do

  ! pure tridiagonal matrix for atmosphere and ocean (no sea-ice)
  call thomas(yy(:,1:kl+ol),aad(:,2:kl+ol),bbd(:,1:kl+ol),ccd(:,1:kl+ol-1),ddd(:,1:kl+ol))

  do k = 1,ol
    do i = 1,imax
      w_t(i,k) = yy(i,kl+k)
    end do
  end do

  ! Derive sensible heat flux for water gridpoints
  do i = 1,imax
    wt0_o(i) = 0.
    if ( .not.land(i) ) then
      fg(i) = (1.-fracice(i))*t1(i)*(yy(i,kl+1)+wrtemp-yy(i,kl) &
              -sigkap(1)*(lv*qlg(i,1)+ls*qfg(i,1))/cp)          &
            + fracice(i)*icefg_a(i)
      ! wt0fb_o has already been multiplied by fracice
      wt0_o(i) = wt0rad_o(i) + wt0melt_o(i) + wt0eg_o(i) + wt0fb_o(i)    &
                + (1.-fracice(i))*t1(i)*(yy(i,kl+1)+wrtemp-yy(i,kl)      &
                  -sigkap(1)*(lv*qlg(i,1)+ls*qfg(i,1))/cp)               &
                  /(wrtrho*cp0)
    end if
  end do
  
else ! use_ocean ..else.

  do i = 1,imax
    wt0_a = fg(i)/(rhos(i)*cp)  ! theta flux
    ddd(i,kl) = ddd(i,kl) + ddts*rhos(i)*wt0_a/(rhoa1(i)*dz1(i))
  end do
    
  ! pure tridiagonal matrix for atmosphere and ocean (no sea-ice)
  call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
    
end if ! use_ocean ..else..


do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    thetal(i,k) = yy(i,kn)
  end do
end do

  
  


!------------------------------------------------------------------------------
! set-up decoupled atmosphere matrices (qv, ql, qf, ni, cf)

do i = 1,imax
  bbd(i,1) = 1.-qq(i,kl)+ddts*mflx(i,kl)*fzzh(i,kl-1)*idzm(i,kl)
  ccd(i,1) = qq(i,kl)+ddts*mflx(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)
end do  
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    aad(i,k) = rr(i,kn)-ddts*mflx(i,kn+1)*fzzh(i,kn)*idzp(i,kn)
    bbd(i,k) = 1.-qq(i,kn)-rr(i,kn)+ddts*(mflx(i,kn)*fzzh(i,kn-1)*idzm(i,kn)      &
                                         -mflx(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn))
    ccd(i,k) = qq(i,kn)+ddts*mflx(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)
  end do  
end do
do i = 1,imax
  aad(i,kl) = rr(i,1)-ddts*mflx(i,2)*fzzh(i,1)*idzp(i,1)
  bbd(i,kl) = 1.-rr(i,1)-ddts*mflx(i,1)*(1.-fzzh(i,1))*idzp(i,1)
end do


! Note that vertical interpolation is linear so that qtot can be
! decomposed into qv, ql and qf.

! qv (part of qtot) - atmosphere
do i = 1,imax
  ddd(i,1)=qvg(i,kl)+ddts*(mflx(i,kl-1)*qvup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)     &
                          +mflx(i,kl)*qvup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k)=qvg(i,kn)+ddts*(mflx(i,kn-1)*qvup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)   &
                            +mflx(i,kn)*qvup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)            &
                            -mflx(i,kn)*qvup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)         &
                            -mflx(i,kn+1)*qvup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do  
end do
do i = 1,imax
  wq0_a = eg(i)/(rhos(i)*lv)
  ddd(i,kl)=qvg(i,1)-ddts*(mflx(i,1)*qvup(i,1)*(1.-fzzh(i,1))*idzp(i,1) &
                          +mflx(i,2)*qvup(i,2)*fzzh(i,1)*idzp(i,1))     &
                    +ddts*rhos(i)*wq0_a/(rhoa1(i)*dz1(i))
end do
call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    qvg(i,k) = yy(i,kn)
  end do
end do

! ql (part of qtot) - atmosphere
do i = 1,imax
  ddd(i,1)=qlg(i,kl)+ddts*(mflx(i,kl-1)*qlup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)     &
                          +mflx(i,kl)*qlup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k)=qlg(i,kn)+ddts*(mflx(i,kn-1)*qlup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)  &
                            +mflx(i,kn)*qlup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)           &
                            -mflx(i,kn)*qlup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)        &
                            -mflx(i,kn+1)*qlup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do  
end do
do i = 1,imax
  ddd(i,kl)=qlg(i,1)-ddts*(mflx(i,1)*qlup(i,1)*(1.-fzzh(i,1))*idzp(i,1)         &
                          +mflx(i,2)*qlup(i,2)*fzzh(i,1)*idzp(i,1))
end do  
call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    qlg(i,k) = yy(i,kn)
  end do  
end do

! qf (part of qtot) - atmosphere
do i = 1,imax
  ddd(i,1)=qfg(i,kl)+ddts*(mflx(i,kl-1)*qfup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)     &
                          +mflx(i,kl)*qfup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do  
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k)=qfg(i,kn)+ddts*(mflx(i,kn-1)*qfup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)   &
                            +mflx(i,kn)*qfup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)            &
                            -mflx(i,kn)*qfup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)         &
                            -mflx(i,kn+1)*qfup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do  
end do
do i = 1,imax
  ddd(i,kl)=qfg(i,1)-ddts*(mflx(i,1)*qfup(i,1)*(1.-fzzh(i,1))*idzp(i,1)         &
                          +mflx(i,2)*qfup(i,2)*fzzh(i,1)*idzp(i,1))
end do
call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    qfg(i,k) = yy(i,kn)
  end do  
end do

! ni - atmosphere
do i = 1,imax
  ddd(i,1)=ni(i,kl)+ddts*(mflx(i,kl-1)*niup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)     &
                         +mflx(i,kl)*niup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k)=ni(i,kn)+ddts*(mflx(i,kn-1)*niup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)   &
                           +mflx(i,kn)*niup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)            &
                           -mflx(i,kn)*niup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)         &
                           -mflx(i,kn+1)*niup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do
end do
do i = 1,imax
  ddd(i,kl)=ni(i,1)-ddts*(mflx(i,1)*niup(i,1)*(1.-fzzh(i,1))*idzp(i,1)         &
                         +mflx(i,2)*niup(i,2)*fzzh(i,1)*idzp(i,1))
end do
call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    ni(i,k) = yy(i,kn)
  end do  
end do

! cf - atmosphere
do i = 1,imax
  ddd(i,1)=cf(i,kl)+ddts*(mflx(i,kl-1)*cfup(i,kl-1)*(1.-fzzh(i,kl-1))*idzm(i,kl)     &
                         +mflx(i,kl)*cfup(i,kl)*fzzh(i,kl-1)*idzm(i,kl))
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k)=cf(i,kn)+ddts*(mflx(i,kn-1)*cfup(i,kn-1)*(1.-fzzh(i,kn-1))*idzm(i,kn)   &
                           +mflx(i,kn)*cfup(i,kn)*fzzh(i,kn-1)*idzm(i,kn)            &
                           -mflx(i,kn)*cfup(i,kn)*(1.-fzzh(i,kn))*idzp(i,kn)         &
                           -mflx(i,kn+1)*cfup(i,kn+1)*fzzh(i,kn)*idzp(i,kn))
  end do
end do
do i = 1,imax
  ddd(i,kl)=cf(i,1)-ddts*(mflx(i,1)*cfup(i,1)*(1.-fzzh(i,1))*idzp(i,1)        &
                         +mflx(i,2)*cfup(i,2)*fzzh(i,1)*idzp(i,1))
end do
call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    cf(i,k) = min( max( yy(i,kn), 0. ), 1. )
  end do  
end do




if ( use_ocean ) then

  !------------------------------------------------------------------------------
  ! set-up decoupled ocean matrices (salinity)

  ! ocean data at kl+k
  do i = 1,imax
    if ( deptho_dz(i,2)*deptho_dz(i,1)>1.e-4 ) then
      ccd(i,kl+1) = -ddts*ks_o(i,2)/(deptho_dz_hl(i,2)*deptho_dz(i,1))
    end if
    bbd(i,kl+1) = 1. - ccd(i,kl+1)
  end do
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k-1)*deptho_dz(i,k)>1.e-4 ) then
        aad(i,kl+k) = -ddts*ks_o(i,k)/(deptho_dz_hl(i,k)*deptho_dz(i,k))
      end if
      if ( deptho_dz(i,k+1)*deptho_dz(i,k)>1.e-4 ) then
        ccd(i,kl+k) = -ddts*ks_o(i,k+1)/(deptho_dz_hl(i,k+1)*deptho_dz(i,k))
      end if
      bbd(i,kl+k) = 1. - aad(i,kl+k) - ccd(i,kl+k)
    end do
  end do
  do i = 1,imax
    if ( deptho_dz(i,ol-1)*deptho_dz(i,ol)>1.e-4 ) then
      aad(i,kl+ol) = -ddts*ks_o(i,ol)/(deptho_dz_hl(i,ol)*deptho_dz(i,ol))
    end if
    bbd(i,kl+ol) = 1. - aad(i,kl+ol)
  end do

  ! sal - ocean
  do i = 1,imax
    fulldeptho(i) = 0.
    targetdeptho(i) = 0.
  end do
  do k = 1,ol
    do i = 1,imax
      targetdeptho(i) = targetdeptho(i) + deptho_dz(i,k)
      targetdeptho(i) = min( minsfc, targetdeptho(i) )
    end do
  end do
  do k = 1,ol
    do i = 1,imax
      !ddd(i,kl+k) = w_s(i,k)
      ddd(i,kl+k) = w_s(i,k)
      if ( deptho_dz(i,k)>1.e-4 ) then
        deltazo = max( min( deptho_dz(i,k), targetdeptho(i)-fulldeptho(i) ), 0. )
        fulldeptho(i) = fulldeptho(i) + deltazo
        !ddd(i,kl+k) = ddd(i,kl+k) - ddts*ws0subsurf_o(iq)*deltazo/max(deptho_dz(i,k)*targetdeptho(i),1.e-4)
        bbd(i,kl+k) = bbd(i,kl+k) + ddts*ws0subsurf_o(i)*deltazo/max(deptho_dz(i,k)*targetdeptho(i),1.e-4)
      end if
    end do  
  end do
  do i = 1,imax
    if ( .not.land(i) ) then
      !ddd(i,kl+1) = ddd(i,kl+1) - ddts*ws0_o(i)/deptho_dz_t(i,1)
      bbd(i,kl+1) = bbd(i,kl+1) + ddts*ws0_o(i)/deptho_dz(i,1)
    end if
  end do
  call thomas(w_s(:,1:ol),aad(:,kl+2:kl+ol),bbd(:,kl+1:kl+ol),ccd(:,kl+1:kl+ol-1),ddd(:,kl+1:kl+ol))

end if ! use_ocean




!------------------------------------------------------------------------------
! momentum for atmosphere, ocean and sea-ice ----

! Original coupling matrix (inverted levels)
! [ bb_a cc_a                                                   ] [ t_a ]   [ dd_a ]
! [ aa_a bb_a cc_a                                              ] [ t_a ]   [ dd_a ]
! [      .... .... ....                                         ] [ ... ]   [ .... ]
! [           aa_a bb_a f_ao                f_ai                ] [ t_a ]   [ dd_a ]
! [                f_oa bb_o cc_o                          f_oi ] [ t_o ]   [ dd_o ]
! [                     aa_o bb_o cc_o                          ] [ t_o ]   [ dd_o ]
! [                          .... .... ....                     ] [ ... ] = [ .... ]
! [                               aa_o bb_o 0.                  ] [ t_o ]   [ dd_o ]
! [                f_ia                0.   bb_i cc_i           ] [ t_i ]   [ dd_i ]
! [                                         aa_i bb_i cc_i      ] [ t_i ]   [ dd_i ]
! [                                              .... .... .... ] [ ... ]   [ .... ]
! [                     f_io                          aa_i bb_i ] [ t_i ]   [ dd_i ]

! u^T = [ .. -f_ai ..  1 .. -f_io ]
! v^T = [ ..  f_ai .. -1 ..  f_io ]

! Improved version (momentum)
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
! t = y - {(v^T y)/(1 + (v^T q))} q
  

do i = 1,imax
  bbd(i,1) = 1.-qq(i,kl)
  ccd(i,1) = qq(i,kl)
end do
do k = 2,kl-1
  kn = kl - k + 1
  do i = 1,imax
    aad(i,k) = rr(i,kn)
    bbd(i,k) = 1.-qq(i,kn)-rr(i,kn)
    ccd(i,k) = qq(i,kn)
  end do
end do
do i = 1,imax
  aad(i,kl) = rr(i,1)
  bbd(i,kl) = 1.-rr(i,1)
end do


! coupling terms
if ( use_ocean ) then
  
  do i = 1,imax
    t1(i) = rhos(i)*cd_water(i)
    t3(i) = rhos(i)*cd_ice(i)
    t4(i) = wrtrho*cdbot_ice(i)
  end do
  
  do i = 1,imax
    if ( land(i) ) then
      bbd(i,kl) = bbd(i,kl) + ddts*rhos(i)*cduv(i)/(rhoa1(i)*dz1(i)) ! implicit
    else
      bbd(i,kl) = bbd(i,kl) + ddts*(t1(i)*(1.-fracice(i))+t3(i)*fracice(i))/(rhoa1(i)*dz1(i))
    end if
  end do

  ! ocean is kl+k
  do i = 1,imax
    if ( deptho_dz(i,2)*deptho_dz(i,1) > 1.e-4 ) then
      ccd(i,kl+1) = -ddts*km_o(i,2)/(deptho_dz_hl(i,2)*deptho_dz(i,1))
    end if
    bbd(i,kl+1) = 1. - ccd(i,kl+1)
    if ( .not.land(i) ) then
      bbd(i,kl+1) = bbd(i,kl+1) + ddts*(t1(i)*(1.-fracice(i))+t4(i)*fracice(i))/(wrtrho*deptho_dz(i,1))
    end if
  end do
  
  do k = 2,ol-1
    do i = 1,imax
      if ( deptho_dz(i,k-1)*deptho_dz(i,k) > 1.e-4 ) then
        aad(i,kl+k) = -ddts*km_o(i,k)/(deptho_dz_hl(i,k)*deptho_dz(i,k))
      end if
      if ( deptho_dz(i,k+1)*deptho_dz(i,k) > 1.e-4 ) then
        ccd(i,kl+k) = -ddts*km_o(i,k+1)/(deptho_dz_hl(i,k+1)*deptho_dz(i,k))
      end if
      bbd(i,kl+k) = 1. - aad(i,kl+k) - ccd(i,kl+k)
    end do
  end do
  do i = 1,imax
    if ( deptho_dz(i,ol-1)*deptho_dz(i,ol) > 1.e-4 ) then
      aad(i,kl+ol) = -ddts*km_o(i,ol)/(deptho_dz_hl(i,ol)*deptho_dz(i,ol))
    end if
    bbd(i,kl+ol) = 1. - aad(i,kl+ol)
  end do  
  ! bottom drag
  do i = 1,imax
    if ( .not.land(i) ) then
      k = ibot(i)
      if ( deptho_dz(i,k)>=1.e-4 ) then
        ! cdbot_water is cd*|u|
        bbd(i,kl+k) = bbd(i,kl+k) + ddts*cdbot_water(i)/deptho_dz(i,k)
      end if
    end if
  end do

  ! sea-ice is kl+ol+k
  do i = 1,imax
    bbd(i,kl+ol+1) = 1.
    if ( .not.land(i) ) then
      bbd(i,kl+ol+1) = 1. + ddts*(t3(i)+t4(i))/imass(i)
    end if
  end do
  
  do k = 1,kl+ol+1
    do i = 1,imax  
      uu(i,k) = 0.
      !vv(i,k) = 0.
    end do
  end do
  do i = 1,imax
    if ( .not.land(i) ) then
      f_ao(i) = -ddts*t1(i)*(1.-fracice(i))/(rhoa1(i)*dz1(i))
      f_oa(i) = -ddts*t1(i)*(1.-fracice(i))/(wrtrho*deptho_dz(i,1))
      f_ai(i) = -ddts*t3(i)*fracice(i)/(rhoa1(i)*dz1(i))
      f_ia(i) = -ddts*t3(i)/imass(i)
      f_oi(i) = -ddts*t4(i)*fracice(i)/(wrtrho*deptho_dz(i,1))
      f_io(i) = -ddts*t4(i)/imass(i)
      bbd(i,kl) = bbd(i,kl) + f_ai(i)*f_ia(i)     ! due to u v^T matrix
      ccd(i,kl) = f_ao(i) + f_ai(i)*f_io(i)       ! due to u v^T matrix
      aad(i,kl+1) = f_oa(i) + f_oi(i)*f_ia(i)     ! due to u v^T matrix
      bbd(i,kl+1) = bbd(i,kl+1) + f_oi(i)*f_io(i) ! due to u v^T matrix
      bbd(i,kl+ol+1) = bbd(i,kl+ol+1) + 1.        ! due to u v^T matrix
      ! Construct u and v
      uu(i,kl) = f_ai(i)
      uu(i,kl+1) = f_oi(i)
      uu(i,kl+ol+1) = -1.
      !vv(i,kl) = -f_ia(i)
      !vv(i,kl+1) = -f_io(i)
      !vv(i,kl+ol+1) = 1.
    end if
  end do
 
  ! Solve A' q = u for combined atmosphere, ocean and sea-ice
  call thomas(yyu(:,1:kl+ol+1),aad(:,2:kl+ol+1),bbd(:,1:kl+ol+1),ccd(:,1:kl+ol),uu(:,1:kl+ol+1))
  
  ! Solve  v^t q
  do i = 1,imax
    !tmp(i) = vv(i,kl)*yyu(i,kl) + vv(i,kl+1)*yyu(i,kl+1) + vv(i,kl+ol+1)*yyu(i,kl+ol+1)
    tmpu(i) = (-f_ia(i))*yyu(i,kl) + (-f_io(i))*yyu(i,kl+1) + (1.)*yyu(i,kl+ol+1)
  end do
  
end if ! use_ocean  
  



! ua, uo, ui - coupled ----

! k=1 is top of atmosphere and k=kl is bottom of atmosphere
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k) = ua(i,kn)
  end do
end do


if ( use_ocean ) then

  do k = 1,ol ! ocean at kl+k
    do i = 1,imax
      ddd(i,kl+k) = w_u(i,k)
    end do
  end do
  do i = 1,imax ! ice at kl+ol+k
    ddd(i,kl+ol+1) = 0.
    if ( .not.land(i) ) then
      ddd(i,kl+ol+1) = i_u(i)
    end if
  end do  
  
  ! Solve A' y = d (for u and v) for combined atmosphere, ocean and sea-ice
  call thomas(yy(:,1:kl+ol+1),aad(:,2:kl+ol+1),bbd(:,1:kl+ol+1),ccd(:,1:kl+ol),ddd(:,1:kl+ol+1))

  ! Solve  v^t y
  do i = 1,imax
    !tmp(i) = vv(i,kl)*yy(i,kl) + vv(i,kl+1)*yy(i,kl+1) + vv(i,kl+ol+1)*yy(i,kl+ol+1)
    tmp(i) = (-f_ia(i))*yy(i,kl) + (-f_io(i))*yy(i,kl+1) + (1.)*yy(i,kl+ol+1)
  end do

  ! Solve for x = y - {(v^t y)/(1 + (v^t q))} q
  do k = 1,kl
    kn = kl - k + 1
    do i = 1,imax
      ua(i,k) = yy(i,kn) - yyu(i,kn)*tmp(i)/(1.+tmpu(i))
    end do  
  end do
  do k = 1,ol ! ocean at kl+k
    do i = 1,imax
      w_u(i,k) = yy(i,kl+k) - yyu(i,kl+k)*tmp(i)/(1.+tmpu(i))
    end do  
  end do
  do i = 1,imax ! ice at kl+ol+k
    i_u(i) = yy(i,kl+ol+1) - yyu(i,kl+ol+1)*tmp(i)/(1.+tmpu(i))
  end do

else ! use_ocean ..else..
    
  do i = 1,imax
    bbd(i,kl) = bbd(i,kl) + ddts*rhos(i)*cduv(i)/(rhoa1(i)*dz1(i)) ! implicit
  end do
  call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
  do k = 1,kl
    kn = kl - k + 1
    do i = 1,imax
      ua(i,k) = yy(i,kn)
    end do  
  end do
  
end if ! use_ocean ..else..



  
! va, vo, vi - coupled ----

! k=1 is top of atmosphere and k=kl is bottom of atmosphere
do k = 1,kl
  kn = kl - k + 1
  do i = 1,imax
    ddd(i,k) = va(i,kn)
  end do
end do
  

if ( use_ocean ) then
  
  do k = 1,ol ! ocean at kl+k
    do i = 1,imax
      ddd(i,kl+k) = w_v(i,k)
    end do
  end do
  do i = 1,imax ! ice at kl+ol+k
    ddd(i,kl+ol+1) = 0.
    if ( .not.land(i) ) then
      ddd(i,kl+ol+1) = i_v(i)
    end if
  end do  
  
  ! Solve A' y = d (for u and v) for combined atmosphere, ocean and sea-ice
  call thomas(yy(:,1:kl+ol+1),aad(:,2:kl+ol+1),bbd(:,1:kl+ol+1),ccd(:,1:kl+ol),ddd(:,1:kl+ol+1))

  ! Solve  v^t y
  do i = 1,imax
    !tmp(i) = vv(i,kl)*yy(i,kl) + vv(i,kl+1)*yy(i,kl+1) + vv(i,kl+ol+1)*yy(i,kl+ol+1)
    tmp(i) = (-f_ia(i))*yy(i,kl) + (-f_io(i))*yy(i,kl+1) + (1.)*yy(i,kl+ol+1)
  end do

  ! Solve for x = y - {(v^t y)/(1 + (v^t q))} q
  do k = 1,kl
    kn = kl - k + 1
    do i = 1,imax
      va(i,k) = yy(i,kn) - yyu(i,kn)*tmp(i)/(1.+tmpu(i))
    end do  
  end do
  do k = 1,ol ! ocean at kl+k
    do i = 1,imax
      w_v(i,k) = yy(i,kl+k) - yyu(i,kl+k)*tmp(i)/(1.+tmpu(i))
    end do  
  end do
  do i = 1,imax ! ice at kl+ol+k
    i_v(i) = yy(i,kl+ol+1) - yyu(i,kl+ol+1)*tmp(i)/(1.+tmpu(i))
  end do
  
else ! use_ocean ..else..

  ! Solve A' y = d (for u and v) for combined atmosphere, ocean and sea-ice
  call thomas(yy(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl))
  do k = 1,kl
    kn = kl - k + 1
    do i = 1,imax
      va(i,k) = yy(i,kn)
    end do  
  end do

end if ! use_ocean ..else..


if ( use_ocean ) then

  do i = 1,imax
    ! update surface momentum flux
    if ( land(i) ) then
      ustar(i) = sqrt(cduv(i)*sqrt(ua(i,1)**2+va(i,1)**2))
    else
      ustar(i) = sqrt((1.-fracice(i))*cd_water(i)*sqrt((ua(i,1)-w_u(i,1))**2+(va(i,1)-w_v(i,1))**2)  &
                +fracice(i)*cd_ice(i)*sqrt((ua(i,1)-i_u(i))**2+(va(i,1)-i_v(i))**2))
    end if
    
    ! update wu0 and wv0
    wu0_o(i) = -(1.-fracice(i))*rhoa1(i)*cd_water(i)*(ua(i,1)-w_u(i,1))/wrtrho
    wv0_o(i) = -(1.-fracice(i))*rhoa1(i)*cd_water(i)*(va(i,1)-w_v(i,1))/wrtrho
    wu0_o(i) = wu0_o(i) + fracice(i)*cdbot_ice(i)*(w_u(i,1)-i_u(i))
    wv0_o(i) = wv0_o(i) + fracice(i)*cdbot_ice(i)*(w_v(i,1)-i_v(i))
  end do
  
else ! use_ocean ..else..

  do i = 1,imax
    ! update surface momentum flux
    ustar(i) = sqrt(cduv(i)*sqrt(ua(i,1)**2+va(i,1)**2))    
  end do
    
end if ! use_ocean ..else..

return
end subroutine update_coupled
                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

pure subroutine thomas(outdat,aai,bbi,cci,ddi)

implicit none

real, dimension(:,2:), intent(in) :: aai
real, dimension(:,:), intent(in) :: bbi,ddi
real, dimension(:,:), intent(in) :: cci
real, dimension(:,:), intent(out) :: outdat
real, dimension(size(outdat,1),size(outdat,2)) :: cc,dd
real n_s
integer k, iq, nx, kx

nx = size(bbi,1)
kx = size(bbi,2)

do iq = 1,nx
  cc(iq,1) = cci(iq,1)/bbi(iq,1)
  dd(iq,1) = ddi(iq,1)/bbi(iq,1)
end do

do k = 2,kx-1
  do iq = 1,nx
    n_s = 1./(bbi(iq,k)-cc(iq,k-1)*aai(iq,k))
    cc(iq,k) = cci(iq,k)*n_s
    dd(iq,k) = (ddi(iq,k)-dd(iq,k-1)*aai(iq,k))*n_s
  end do
end do

do iq = 1,nx
  n_s = 1./(bbi(iq,kx)-cc(iq,kx-1)*aai(iq,kx))
  outdat(iq,kx) = (ddi(iq,kx)-dd(iq,kx-1)*aai(iq,kx))*n_s
end do

do k = kx-1,1,-1
  do iq = 1,nx
    outdat(iq,k) = dd(iq,k)-cc(iq,k)*outdat(iq,k+1)
  end do  
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

pure function getqsat(templ,ps,fice) result(qsat)

implicit none

integer iq, ix, ixx
real, dimension(:), intent(in) :: templ
real, dimension(:), intent(in) :: ps, fice
real, dimension(size(templ)) :: qsat
real estafi, tdiff, tdiffx, rx, rxx
real qsatl, qsati, deles

do iq = 1,size(templ)
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
end function getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate lateral entrainment
!
!pure function entfn(zht,zi,imax) result(ans)
!
!implicit none
!
!integer, intent(in) :: imax
!real, dimension(imax), intent(in) :: zht, zi
!real, dimension(imax) :: ans
!
!!ans=0.002                                            ! Angevine (2005)
!!ans=2./max(100.,zi)                                  ! Angevine et al (2010)
!!ans=1./zht                                           ! Siebesma et al (2003)
!!ans=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))    ! Soares et al (2004)
!ans = max( ent0/max( zht, 1. ) + ent1/max( zi-zht, ezmin ), ent_min )
!
!return
!end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE-eps

subroutine tkeend

implicit none

deallocate(tke,eps)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps

