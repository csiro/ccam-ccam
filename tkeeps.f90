! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
#ifdef CCAM
                  land,tile,                                                               &
#endif
                  imax,kl)

#ifdef CCAM
use mlo, only : wlev
#endif
                  
implicit none

integer, intent(in) :: imax, kl, mode
integer k, iq
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(imax,kl), intent(inout) :: theta,stratcloud,ua,va
real, dimension(imax,kl), intent(inout) :: qvg,qlg,qfg
real, dimension(imax,kl), intent(out) :: kmo
real, dimension(imax,kl), intent(in) :: zz,zzh
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(inout) :: eps
real, dimension(imax,kl), intent(inout) :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(imax,kl), intent(inout) :: tke_ema
real, dimension(imax), intent(inout) :: zi,fg,eg
real, dimension(imax), intent(in) :: cduv,ps,rhos,dx
real, dimension(imax), intent(out) :: ustar_ave
real, dimension(imax), intent(in) :: sig
real, dimension(imax,kl) :: km, thetav, thetal, qsat
real, dimension(imax,kl) :: rr
real, dimension(imax,kl) :: rhoa,rhoahl
real, dimension(imax,kl) :: tlup,qvup,qlup,qfup
real, dimension(imax,kl) :: cfup,mflx
real, dimension(imax,kl) :: pps,ppt,ppb
real, dimension(imax,kl) :: idzm
real, dimension(imax,kl-1) :: idzp
real, dimension(imax,2:kl) :: qq
real, dimension(imax,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(imax,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(imax,kl-1) :: fzzh
real, dimension(imax) :: wt0,wq0,wtv0
real, dimension(imax) :: wstar,z_on_l,phim
real, dimension(imax) :: pres, temp
real, dimension(imax) :: ustar, fg_ave
real, dimension(imax) :: zi_save, zturb, cgmap
real, dimension(imax) :: templ, fice, al, dqsdt, qc, qt, lx
real, dimension(kl) :: sigkap
real tff, cm12, cm34, ddts
logical, dimension(imax) :: mask

#ifdef CCAM
integer, intent(in) :: tile
integer, dimension(imax) :: ibot
real, dimension(imax,wlev) :: deptho_dz
real, dimension(imax,2:wlev) :: deptho_dz_hl
real, dimension(imax,wlev) :: rad_o
real, dimension(imax) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o, wt0eg_o
real, dimension(imax,wlev) :: w_t, w_s, w_u, w_v
real, dimension(imax) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o
real, dimension(imax) :: fracice, i_u, i_v, imass, cd_ice, cdbot_ice
logical, dimension(imax), intent(in) :: land
#endif

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl), intent(out) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(out) :: totaltransport
real, dimension(imax,kl), intent(out) :: mfout
#endif

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

sigkap(1:kl) = sig(1:kl)**(-rd/cp)

do k = 1,kl
  do iq = 1,imax
    ! Impose limits on tke and eps after being advected by the host model
    tke(iq,k) = max(tke(iq,k), mintke)
    tff       = cm34*tke(iq,k)*sqrt(tke(iq,k))
    eps(iq,k) = min(eps(iq,k), tff/minl)
    eps(iq,k) = max(eps(iq,k), tff/maxl, mineps)
  
    ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
    pres(iq) = ps(iq)*sig(k) ! pressure
    ! Density must be updated when dz is updated so that rho*dz is conserved
    thetav(iq,k) = theta(iq,k)*(1.+0.61*qvg(iq,k)-qlg(iq,k)-qfg(iq,k))
    rhoa(iq,k) = sigkap(k)*pres(iq)/(rd*thetav(iq,k))

    ! Transform to thetal as it is the conserved variable
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
dz_fl = max( dz_fl, 1. )

! Calculate shear term on full levels
pps(:,1:kl-1) = km(:,1:kl-1)*shear(:,1:kl-1)

! interpolate diffusion coeff and air density to half levels
call updatekmo(kmo,   km,  fzzh,imax,kl)
call updatekmo(rhoahl,rhoa,fzzh,imax,kl)
! eddy diffusion terms to account for air density with level thickness
idzm(:,2:kl)   = rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1) = rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

ustar_ave(:) = 0.
fg_ave(:) = 0.

#ifdef CCAM
if ( mode==2 .or. mode==3 ) then
  call unpack_coupled(deptho_dz,deptho_dz_hl,rad_o,                   &
                      w_t,w_s,w_u,w_v,cd_water,cdh_water,cdbot_water, &
                      wt0rad_o,wt0melt_o,wt0eg_o,                     &
                      icefg_a,wt0fb_o,ws0_o,ws0subsurf_o,i_u,i_v,     &
                      imass,fracice,cd_ice,cdbot_ice,ibot,imax,tile)
end if
#endif

! Main loop to prevent time splitting errors
mcount = int(dt/(min(maxdts,12./max(m0,0.1))+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  ! Set-up thermodynamic variables temp, theta_v, qtot and saturated mixing ratio
  do k = 1,kl
    thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
  end do
  
  ! time averaging
  call update_ema(thetal,thetal_ema,ddts)
  call update_ema(qvg,qv_ema,ddts)
  call update_ema(qlg,ql_ema,ddts)
  call update_ema(qfg,qf_ema,ddts)
  call update_ema(stratcloud,cf_ema,ddts)
  call update_ema(tke,tke_ema,ddts)
  
  ! Calculate surface fluxes
  wt0 = fg/(rhos*cp)  ! theta flux
  wq0 = eg/(rhos*lv)  ! qtot flux  
  wtv0 = wt0 + theta(:,1)*0.61*wq0 ! thetav flux
  
  ! Update momentum flux
  ustar = sqrt(cduv*sqrt(ua(:,1)**2+va(:,1)**2))  
  wstar = (grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)   
  
  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  if ( mode/=1 .and. mode/=3 ) then ! mass flux
      
    zi_save = zi  

    ! plume rise model
    mask = wtv0>0.
    tke(:,1) = cm12*ustar**2 + ce3*wstar**2
    tke(:,1) = max(tke(:,1), mintke)
    call plumerise(mask,                                         &
                   zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,       &
                   zz,dz_hl,thetal,qvg,qlg,qfg,                  &
                   stratcloud,wt0,wq0,ps,ustar,                  &
                   sig,sigkap,tke,eps,ua,va,imax,kl)

#ifndef scm
    ! Turn off MF term if small grid spacing (mfbeta=0 implies MF is always non-zero)
    ! Based on Boutle et al 2014
    zturb = 0.5*(zi_save + zi)
    cgmap(:) = 1. - tanh(mfbeta*zturb/dx)*max(0.,1.-0.25*dx/zturb)
    do k = 1,kl
      mflx(:,k) = mflx(:,k)*cgmap(:)
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
  z_on_l = -vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar**3,1.E-10))
  z_on_l = min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  call calc_phi(phim,z_on_l,imax)
  do iq = 1,imax
    tke(iq,1) = cm12*ustar(iq)**2+ce3*wstar(iq)**2
    eps(iq,1) = ustar(iq)**3*phim(iq)/(vkar*zz(iq,1))+grav*wtv0(iq)/thetav(iq,1)
    tke(iq,1) = max( tke(iq,1), mintke )
    tff = cm34*tke(iq,1)*sqrt(tke(iq,1))
    eps(iq,1) = min( eps(iq,1), tff/minl )
    eps(iq,1) = max( eps(iq,1), tff/maxl, mineps )
  end do


  ! Update TKE and eps terms
  call update_tkeeps(tke,eps,ppb,pps,ppt,qvg,qlg,qfg,stratcloud,thetal,         &
                     ua,va,zzh,zz,sigkap,idzp,idzm,wstar,zi,thetal_ema,qv_ema,  &
                     ql_ema,qf_ema,cf_ema,tke_ema,ddts,qgmin,cm34,imax,kl)
  
#ifdef CCAM
  if ( mode==2 .or. mode==3 ) then
    call update_coupled(thetal,qvg,qlg,qfg,stratcloud,ua,va,    &
                        tlup,qvup,qlup,qfup,cfup,fg,eg,rhos,    &
                        ustar,cduv,ps,                          &
                        tke,eps,mflx,fzzh,idzp,idzm,dz_hl,      &
                        rhoa(:,1),dz_fl(:,1),                   &
                        deptho_dz,deptho_dz_hl,                 &
                        rad_o,w_t,w_s,w_u,w_v,                  &
                        cd_water,cdh_water,cdbot_water,         &
                        wt0rad_o,wt0melt_o,                     &
                        wt0eg_o,icefg_a,wt0fb_o,ws0_o,          &
                        ws0subsurf_o,i_u,i_v,imass,fracice,     &
                        cd_ice,cdbot_ice,ibot,land,sigkap,      &
#ifdef scm
                        wthflux,wqvflux,uwflux,vwflux,          &
#endif
                        ddts,imax,kl,tile)  
  else
    call update_atmosphere(thetal,qvg,qlg,qfg,stratcloud,ua,va, &
                           tlup,qvup,qlup,qfup,cfup,fg,eg,rhos, &
                           ustar,cduv,                          &
                           tke,eps,mflx,fzzh,idzp,idzm,dz_hl,   &
                           rhoa(:,1),dz_fl(:,1),sigkap,         &
#ifdef scm
                           wthflux,wqvflux,uwflux,vwflux,       &
#endif
                           ddts,imax,kl)
  end if
#else
  call update_atmosphere(thetal,qvg,qlg,qfg,stratcloud,ua,va, &
                         tlup,qvup,qlup,qfup,cfup,fg,eg,rhos, &
                         ustar,cduv,                          &
                         tke,eps,mflx,fzzh,idzp,idzm,dz_hl,   &
                         rhoa(:,1),dz_fl(:,1),sigkap,         &
#ifdef scm
                         wthflux,wqvflux,uwflux,vwflux,       &
#endif
                         ddts,imax,kl)
#endif
  
  ustar_ave = ustar_ave + ustar/real(mcount)
  fg_ave = fg_ave + fg/real(mcount)

  ! Account for phase transistions
  do k = 1,kl
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
    fice = min( qfg(:,k)/max(qfg(:,k)+qlg(:,k),1.e-8), 1. )
    call getqsat(qsat(:,k),templ(:),pres(:),fice,imax)
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

end do

#ifdef CCAM
if ( mode==2 .or. mode==3 ) then
  call pack_coupled(w_t,w_s,w_u,w_v,i_u,i_v,imax,tile)
  fg = fg_ave
end if
#endif

#ifdef scm
do k = 1,kl
  buoyproduction(:,k) = ppb(:,k)
  shearproduction(:,k) = pps(:,k)
  totaltransport(:,k) = ppt(:,k)
end do
#endif

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plume rise model
    
pure subroutine plumerise(mask,                                        &
                     zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,           &
                     zz,dz_hl,thetal,qvg,qlg,qfg,                      &
                     stratcloud,wt0,wq0,ps,ustar,                      &
                     sig,sigkap,tke,eps,ua,va,imax,kl)

integer, intent(in) :: imax, kl
integer k
real, dimension(imax,kl), intent(inout) :: mflx, tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: qvg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: zz, thetal, ua, va 
real, dimension(imax,kl), intent(in) :: dz_hl
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(in) :: eps
real, dimension(imax), intent(in) :: wt0, wq0, ps, ustar
real, dimension(kl), intent(in) :: sig, sigkap
real, dimension(imax), intent(inout) :: zi, wstar
real, dimension(imax,kl) ::  w2up, nn, cxup, rino
real, dimension(imax,kl) :: theta, thetav, km
real, dimension(imax) :: dzht, ent, templ, pres, upf, qxup, qupsat
real, dimension(imax) :: fice, lx, qcup, dqsdt, al, xp, as, bs, cs
real, dimension(imax) :: thup, tvup, qtup, vvk, wtv0
real, parameter :: fac = 10. ! originally fac=100.
real, parameter :: ricr = 0.3
logical, dimension(imax), intent(in) :: mask

! Update theta and thetav
do k = 1,kl
  theta(:,k) = thetal(:,k) + sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
  thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
  km(:,k) = cm0*tke(:,k)**2/eps(:,k)
end do 

wtv0(:) = wt0(:) + theta(:,1)*0.61*wq0(:) ! thetav flux

! Initialise updraft
do k = 1,kl
  mflx(:,k) = 0.  
  w2up(:,k) = 0.
  nn(:,k) = 0.
  tlup(:,k) = thetal(:,k)
  qvup(:,k) = qvg(:,k)
  qlup(:,k) = qlg(:,k)
  qfup(:,k) = qfg(:,k)
  cfup(:,k) = stratcloud(:,k)
end do

! first level -----------------

dzht = zz(:,1)
! Entrainment rates
ent = entfn(zz(:,1), zi(:), imax)

! initial thermodynamic state
! split qtot into components (conservation of thetal and qtot is maintained)
where ( mask )
  tlup(:,1) = thetal(:,1) + be*wt0/sqrt(max(tke(:,1),1.5e-4))       ! Hurley 2007
  qvup(:,1) = qvg(:,1)    + be*wq0/sqrt(max(tke(:,1),1.5e-4))       ! Hurley 2007
  qlup(:,1) = qlg(:,1)
  qfup(:,1) = qfg(:,1)
  cfup(:,1) = stratcloud(:,1)
end where
! update updraft velocity and mass flux
nn(:,1) = grav*be*wtv0/(thetav(:,1)*sqrt(max(tke(:,1),1.5e-4))) ! Hurley 2007
w2up(:,1) = 2.*dzht*b2*nn(:,1)/(1.+2.*dzht*b1*ent)              ! Hurley 2007
cxup(:,1) = 0.
rino(:,1) = 0.


! updraft with condensation
do k = 2,kl
  dzht = dz_hl(:,k-1)
  ! Entrainment rates
  ent = entfn(zz(:,k), zi(:), imax)
  where ( w2up(:,k-1)>0. .and. mask )
    ! entrain air into plume
    ! split qtot into components (conservation of qtot is maintained)
    tlup(:,k) = (tlup(:,k-1)+dzht*ent*thetal(:,k))/(1.+dzht*ent)
    qvup(:,k) = (qvup(:,k-1)+dzht*ent*qvg(:,k)   )/(1.+dzht*ent)
    qlup(:,k) = (qlup(:,k-1)+dzht*ent*qlg(:,k)   )/(1.+dzht*ent)
    qfup(:,k) = (qfup(:,k-1)+dzht*ent*qfg(:,k)   )/(1.+dzht*ent)
    cfup(:,k) = (cfup(:,k-1)+dzht*ent*stratcloud(:,k))/(1.+dzht*ent)
  end where
  ! calculate conserved variables
  qtup = qvup(:,k) + qlup(:,k) + qfup(:,k)    ! qtot,up
  templ = tlup(:,k)/sigkap(k)                 ! templ,up
  pres = ps(:)*sig(k)
  fice = min( qfup(:,k)/max(qfup(:,k)+qlup(:,k),1.e-8), 1. )
  call getqsat(qupsat,templ,pres,fice,imax)
  where ( qtup>qupsat .and. w2up(:,k-1)>0. )
    qxup = qupsat
    cxup(:,k) = 1.
  elsewhere
    qxup = qtup
    cxup(:,k) = 0.
  end where
  lx = lv + lf*fice
  dqsdt = qupsat*lx/(rv*templ**2)
  al = cp/(cp+lx*dqsdt)
  qcup = max(al*(qtup-qxup), 0.)                           ! qcondensate,up after redistribution
  qcup = min(qcup, qcmf)                                   ! limit condensation with simple autoconversion
  thup = tlup(:,k) + sigkap(k)*qcup*lx/cp                  ! theta,up after redistribution
  tvup = thup + theta(:,k)*(0.61*qxup-qcup)                ! thetav,up after redistribution
  where ( w2up(:,k-1)>0. )
    nn(:,k) = grav*(tvup-thetav(:,k))/thetav(:,k)                    ! calculate buayancy
    w2up(:,k) = (w2up(:,k-1)+2.*dzht*b2*nn(:,k))/(1.+2.*dzht*b1*ent) ! update updraft velocity
  elsewhere
    nn(:,k) = 0.  
    w2up(:,k) = 0.
  end where
  vvk(:) = (ua(:,k)-ua(:,1))**2 + (va(:,k)-va(:,1))**2 + fac*ustar(:)**2  
  rino(:,k) = grav*(thetav(:,k)-thetav(:,1))*(zz(:,k)-zz(:,1))/max(thetav(:,1)*vvk,1.e-30)
  ! test if maximum plume height is reached
  where ( w2up(:,k)<=0. .and. w2up(:,k-1)>0. .and. mask ) ! unstable
    as = 2.*b2*(nn(:,k)-nn(:,k-1))/dzht
    bs = 2.*b2*nn(:,k-1)
    cs = w2up(:,k-1)
    xp = -2.*cs/(bs-sqrt(max(bs**2-4.*as*cs,0.)))
    xp = min(max(xp,0.),dzht)
    zi(:) = xp + zz(:,k-1)
  elsewhere ( rino(:,k)>ricr .and. rino(:,k-1)<=ricr .and. .not.mask ) ! stable
    xp(:) = (ricr-rino(:,k-1))/(rino(:,k)-rino(:,k-1))
    xp(:) = min( max(xp(:), 0.), 1.)
    zi(:) = zz(:,k-1) + xp(:)*(zz(:,k)-zz(:,k-1))
  end where
end do

! update wstar with new zi          
wstar = (grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)
          
! update mass flux
mflx(:,1) = m0*sqrt(max(w2up(:,1), 0.))
do k = 2,kl
  dzht = dz_hl(:,k-1)
  upf = mflx(:,k-1)/sqrt(max(w2up(:,k-1), 1.e-8))
  mflx(:,k) = (1.-cxup(:,k))*m0*sqrt(max(w2up(:,k), 0.))         &
            + cxup(:,k)*mflx(:,k-1)/(1.+dzht*(dtrc0-entc0))
  mflx(:,k) = min( mflx(:,k), upf*sqrt(max(w2up(:,k), 0.)) )
end do

return
end subroutine plumerise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update TKE and EPS

subroutine update_tkeeps(tke,eps,ppb,pps,ppt,qvg,qlg,qfg,stratcloud,thetal,         &
                         ua,va,zzh,zz,sigkap,idzp,idzm,wstar,zi,thetal_ema,qv_ema,  &
                         ql_ema,qf_ema,cf_ema,tke_ema,ddts,qgmin,cm34,imax,kl)

implicit none

integer, intent(in) :: imax, kl
integer k, iq
real, dimension(imax,kl), intent(inout) :: tke, eps
real, dimension(imax,kl), intent(inout) :: ppb, pps, ppt
real, dimension(imax,kl), intent(in) :: qvg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: thetal
real, dimension(imax,kl), intent(in) :: ua, va
real, dimension(imax,kl), intent(in) :: idzm
real, dimension(imax,kl), intent(in) :: zzh, zz
real, dimension(imax,kl), intent(inout) :: thetal_ema, qv_ema, ql_ema, qf_ema, cf_ema
real, dimension(imax,kl), intent(inout) :: tke_ema
real, dimension(imax,kl-1), intent(in) :: idzp
real, dimension(imax,kl) :: qsatc, qgnc, thetalhl, quhl, qshl, qlhl, qfhl
real, dimension(imax,kl) :: qthl, thetavhl, kmo, ua_hl, va_hl, qtot
real, dimension(imax,kl) :: dz_fl, theta, thetav, km
real, dimension(imax,kl) :: rr, bb, cc, dd, ff
real, dimension(imax,kl) :: theta_ema, thetav_ema, qtot_ema
real, dimension(imax,kl-1) :: fzzh, dz_hl
real, dimension(imax,2:kl) :: qq, aa
real, dimension(imax), intent(in) :: wstar, zi
real, dimension(kl), intent(in) :: sigkap
real, intent(in) :: qgmin, ddts, cm34
real :: tbb, tcc, tff, tqq, thetac, tempc, tempt, tempv, rvar, bvf, mc, dc, fc
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
do k = 1,kl
  qtot(:,k) = qvg(:,k) + qlg(:,k) + qfg(:,k)
  theta(:,k) = thetal(:,k) + sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
  thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
  qtot_ema(:,k) = qv_ema(:,k) + ql_ema(:,k) + qf_ema(:,k)  
  theta_ema(:,k) = thetal_ema(:,k) + sigkap(k)*(lv*ql_ema(:,k)+ls*qf_ema(:,k))/cp
  thetav_ema(:,k) = theta_ema(:,k)*(1.+0.61*qv_ema(:,k)-ql_ema(:,k)-qf_ema(:,k))    
end do

! interpolate diffusion coeffs to half levels
km = cm0*tke(:,:)**2/eps(:,:)
call updatekmo(kmo,km,fzzh,imax,kl) 

! top boundary condition to avoid unphysical behaviour at the top of the model
tke(:,kl) = mintke
eps(:,kl) = mineps

! Calculate buoyancy term
select case(buoymeth)
  case(0) ! Blend staturated and unsaturated terms - saturated method from Durran and Klemp JAS 1982 (see also WRF)
    qsatc = qv_ema(:,:)                       ! assume qvg is saturated inside cloud
    ff = qf_ema(:,:)/max(cf_ema,1.E-8)        ! inside cloud value assuming max overlap
    dd = ql_ema(:,:)/max(cf_ema,1.E-8)        ! inside cloud value assuming max overlap
    do k = 1,kl
      do iq = 1,imax
          tbb = max(1.-cf_ema(iq,k),1.E-8)
          qgnc(iq,k) = (qv_ema(iq,k)-(1.-tbb)*qsatc(iq,k))/tbb    ! outside cloud value
          qgnc(iq,k) = min(max(qgnc(iq,k),qgmin),qsatc(iq,k))
        end do
      end do
      call updatekmo(thetalhl,thetal_ema,fzzh,imax,kl)          ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh,imax,kl)                    ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh,imax,kl)                   ! inside cloud value
      call updatekmo(qlhl,dd,fzzh,imax,kl)                      ! inside cloud value
      call updatekmo(qfhl,ff,fzzh,imax,kl)                      ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl) = cf_ema(:,2:kl)<=1.E-6
      do k = 2,kl-1
        where( lta(:,k) .and. .not.lta(:,k+1) )
          qlhl(:,k) = dd(:,k+1)
          qfhl(:,k) = ff(:,k+1)
        elsewhere ( .not.lta(:,k) .and. lta(:,k+1) )
          qlhl(:,k) = dd(:,k)
          qfhl(:,k) = ff(:,k)
        end where
      end do
      do k = 2,kl-1
        ! saturated
        do iq = 1,imax
          thetac = thetal_ema(iq,k)+sigkap(k)*(lv*dd(iq,k)+ls*ff(iq,k))/cp   ! inside cloud value
          tempc = thetac/sigkap(k)                                          ! inside cloud value          
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
        thetac = thetal_ema(iq,1)+sigkap(1)*(lv*dd(iq,1)+ls*ff(iq,1))/cp        ! inside cloud value
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
      call updatekmo(thetalhl,thetal_ema,fzzh,imax,kl)
      call updatekmo(qthl,qtot_ema,fzzh,imax,kl)
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
      call updatekmo(thetavhl,thetav_ema,fzzh,imax,kl)
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
  aa(:,2:kl-1)=ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=ce0*kmo(:,2:kl-1)*rr(:,2:kl-1)
  ! follow Hurley 2007 to make scheme more numerically stable
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)+ddts*ce2*eps(:,2:kl-1)/tke(:,2:kl-1)
  dd(:,2:kl-1)=eps(:,2:kl-1)+ddts*eps(:,2:kl-1)/tke(:,2:kl-1)                    &
              *ce1*(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
  dd(:,2)     =dd(:,2)   -aa(:,2)*eps(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mineps
  call thomas(eps(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax,kl-2)
  
  ! TKE vertical mixing
  aa(:,2:kl-1)=kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=kmo(:,2:kl-1)*rr(:,2:kl-1)
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  dd(:,2:kl-1)=tke(:,2:kl-1)+ddts*(pps(:,2:kl-1)+ppb(:,2:kl-1)-eps(:,2:kl-1))
  dd(:,2)     =dd(:,2)   -aa(:,2)*tke(:,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mintke
  call thomas(tke(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax,kl-2)

  ! limit decay of TKE and EPS with coupling to mass flux term
  if ( tkemeth==1 ) then
    do k = 2,kl-1
      do iq = 1,imax
        tbb = max(1.-0.05*dz_hl(iq,k-1)/250.,0.)
        if ( wstar(iq)>0.5 .and. zz(iq,k)>0.5*zi(iq) .and. zz(iq,k)<0.95*zi(iq) ) then
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

implicit none

integer, intent(in) :: imax, kl
integer k
real, dimension(imax,kl), intent(inout) :: thetal, qvg, qlg, qfg, stratcloud, ua, va
real, dimension(imax,kl), intent(in) :: tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: mflx, tke, eps
real, dimension(imax,kl), intent(in) :: idzm
real, dimension(imax,kl-1), intent(in) :: fzzh, idzp, dz_hl
real, dimension(imax,kl) :: km, kmo
real, dimension(imax,kl) :: rr, bb, cc, dd
real, dimension(imax,2:kl) :: qq, aa
real, dimension(imax), intent(inout) :: fg, eg, ustar
real, dimension(imax), intent(in) :: rhos, rhoa1, dz_fl1, cduv
real, dimension(imax) :: wt0, wq0
real, dimension(kl), intent(in) :: sigkap
real, intent(in) :: ddts

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl) :: wthlflux, wqlflux, wqfflux
#endif

! Calculate surface fluxes
wt0 = fg/(rhos*cp)  ! theta flux
wq0 = eg/(rhos*lv)  ! qtot flux

! estimate eddy diffusivity mixing coeff
km(:,:) = cm0*tke(:,:)**2/eps(:,:)
call updatekmo(kmo,km,fzzh,imax,kl) ! interpolate diffusion coeffs to half levels
  
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
dd(:,1)=thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                               &
                         +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                   &
                   +ddts*rhos*wt0/(rhoa1(:)*dz_fl1(:))
dd(:,2:kl-1)=thetal(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*tlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1) &
                                   +mflx(:,2:kl-1)*tlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)      &
                                   -mflx(:,2:kl-1)*tlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1) &
                                   -mflx(:,3:kl)*tlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl)=thetal(:,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                   &
                           +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)

    
#ifdef scm  
wthlflux(:,1)=wt0(:)
wthlflux(:,2:kl)=-kmo(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)                    &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))               &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif


! qv (part of qtot) vertical mixing
dd(:,1)=qvg(:,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                  &
                         +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                                   &
                         +ddts*rhos*wq0/(rhoa1(:)*dz_fl1(:))
dd(:,2:kl-1)=qvg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qvup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)    &
                                 +mflx(:,2:kl-1)*qvup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                 -mflx(:,2:kl-1)*qvup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                 -mflx(:,3:kl)*qvup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl)=qvg(:,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                      &
                        +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
#ifdef scm
wqvflux(:,1)=wq0(:)
wqvflux(:,2:kl)=-kmo(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)                           &
                +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                   &
                +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif


! ql (part of qtot) vertical mixing
dd(:,1)=qlg(:,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                            +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
dd(:,2:kl-1)=qlg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                      +mflx(:,2:kl-1)*qlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                      -mflx(:,2:kl-1)*qlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                      -mflx(:,3:kl)*qlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl)=qlg(:,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                              +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
#ifdef scm
wqlflux(:,1)=0.
wqlflux(:,2:kl)=-kmo(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)                                &
                +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                        &
                +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif


! qf (part of qtot) vertical mixing
dd(:,1)=qfg(:,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                            +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
dd(:,2:kl-1)=qfg(:,2:kl-1)+ddts*(mflx(:,1:kl-2)*qfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                      +mflx(:,2:kl-1)*qfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)        &
                                      -mflx(:,2:kl-1)*qfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)   &
                                      -mflx(:,3:kl)*qfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
dd(:,kl)=qfg(:,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                              +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
#ifdef scm
wqfflux(:,1)=0.
wqfflux(:,2:kl)=-kmo(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)                                &
                +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                        &
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
call thomas(stratcloud,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
stratcloud(:,:) = min( max( stratcloud(:,:), 0. ), 1. )

  
! momentum vertical mixing
aa(:,2:kl)   = qq(:,2:kl)
cc(:,1:kl-1) = rr(:,1:kl-1)
bb(:,1) = 1. - cc(:,1) + ddts*rhos*cduv/(rhoa1(:)*dz_fl1(:)) ! implicit  
bb(:,2:kl-1) = 1. - aa(:,2:kl-1) - cc(:,2:kl-1)
bb(:,kl) = 1. - aa(:,kl)
dd(:,1:kl) = ua(:,1:kl)
call thomas(ua,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
dd(:,1:kl) = va(:,1:kl)
call thomas(va,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
  
  
! update surface momentum flux
ustar = sqrt(cduv*sqrt(ua(:,1)**2+va(:,1)**2))


#ifdef scm
uwflux(:,1) = 0.
uwflux(:,2:kl) = -kmo(:,1:kl-1)*(ua(:,2:kl)-ua(:,1:kl-1))/dz_hl(:,1:kl-1)
vwflux(:,1) = 0.
vwflux(:,2:kl) = -kmo(:,1:kl-1)*(va(:,2:kl)-va(:,1:kl-1))/dz_hl(:,1:kl-1)
do k = 1,kl
  wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                               *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
end do
#endif

return
end subroutine update_atmosphere

#ifdef CCAM
subroutine unpack_coupled(deptho_dz,deptho_dz_hl,rad_o,                   &
                          w_t,w_s,w_u,w_v,cd_water,cdh_water,cdbot_water, &
                          wt0rad_o,wt0melt_o,wt0eg_o,                     &
                          icefg_a,wt0fb_o,ws0_o,ws0subsurf_o,i_u,i_v,     &
                          imass,fracice,cd_ice,cdbot_ice,ibot,imax,tile)

use mlo, only : mloexport, mloexpdep, mlodiag, mloexpice, mlodiagice, wlev, &
                water_g, dgwater_g, ice_g, dgice_g, wpack_g, wfull_g,       &
                depth_g, minwater

implicit none

integer, intent(in) :: imax, tile
integer, dimension(imax), intent(out) :: ibot
integer k
real, dimension(imax,wlev), intent(out) :: deptho_dz
real, dimension(imax,2:wlev), intent(out) :: deptho_dz_hl
real, dimension(imax,wlev), intent(out) :: rad_o
real, dimension(imax), intent(out) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o, wt0eg_o
real, dimension(imax,wlev), intent(out) :: w_t, w_s, w_u, w_v
real, dimension(imax), intent(out) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o, fracice
real, dimension(imax), intent(out) :: i_u, i_v, imass, cd_ice, cdbot_ice
real, dimension(imax) :: d_zcr, rbot, neta, deptho_max

neta = 0.
call mloexport("eta",neta,0,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
deptho_max = 0.
call mloexpdep("depth_hl",deptho_max,wlev+1,0,depth_g(tile),wpack_g(:,tile),wfull_g(tile))

d_zcr = max( 1. + neta/max(deptho_max,1.e-4), minwater/max(deptho_max,1.e-4) )

do k = 1,wlev
  call mloexpdep("depth_dz_fl",deptho_dz(:,k),k,0,depth_g(tile),wpack_g(:,tile),wfull_g(tile))
  deptho_dz(:,k) = deptho_dz(:,k)*d_zcr
end do
do k = 2,wlev
  call mloexpdep("depth_dz_hl",deptho_dz_hl(:,k),k,0,depth_g(tile),wpack_g(:,tile),wfull_g(tile))
  deptho_dz_hl(:,k) = deptho_dz_hl(:,k)*d_zcr
end do

cd_water = 0.
cdh_water = 0.
cdbot_water = 0.
wt0rad_o = 0.
wt0melt_o = 0.
wt0eg_o = 0.
wt0fb_o = 0.
ws0_o = 0.
ws0subsurf_o = 0.
call mlodiag("cd",cd_water,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("cdh",cdh_water,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("cd_bot",cdbot_water,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("wt0_rad",wt0rad_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("wt0_melt",wt0melt_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("wt0_eg",wt0eg_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("wt0_fb",wt0fb_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("ws0",ws0_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("ws0_subsurf",ws0subsurf_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))

w_t = 0.
w_s = 0.
w_u = 0.
w_v = 0.
do k = 1,wlev
  call mloexport("temp",w_t(:,k),k,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloexport("sal",w_s(:,k),k,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloexport("u",w_u(:,k),k,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloexport("v",w_v(:,k),k,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mlodiag("rad",rad_o(:,k),k,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
end do

rbot = 0.
call mloexport("ibot",rbot,0,0,water_g(tile),wpack_g(:,tile),wfull_g(tile))
ibot = nint(rbot)

i_u = 0.
i_v = 0.
fracice = 0.
call mloexpice("u",i_u,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mloexpice("v",i_v,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mloexpice("fracice",fracice,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))

icefg_a = 0.
imass = 100.
cd_ice = 0.
cdbot_ice = 0.
call mlodiagice("fg",icefg_a,0,dgice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiagice("mass",imass,0,dgice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiagice("cd",cd_ice,0,dgice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiagice("cd_bot",cdbot_ice,0,dgice_g(tile),wpack_g(:,tile),wfull_g(tile))

return
end subroutine unpack_coupled

subroutine pack_coupled(w_t,w_s,w_u,w_v,i_u,i_v,imax,tile)

use mlo, only : mloimport, mloimpice, wlev, water_g, ice_g, depth_g, wpack_g, wfull_g

implicit none

integer, intent(in) :: tile, imax
integer k
real, dimension(imax,wlev), intent(in) :: w_t, w_s, w_u, w_v
real, dimension(imax), intent(in) :: i_u, i_v

do k = 1,wlev
  call mloimport("temp",w_t(:,k),k,0,water_g(tile),depth_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloimport("sal",w_s(:,k),k,0,water_g(tile),depth_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloimport("u",w_u(:,k),k,0,water_g(tile),depth_g(tile),wpack_g(:,tile),wfull_g(tile))
  call mloimport("v",w_v(:,k),k,0,water_g(tile),depth_g(tile),wpack_g(:,tile),wfull_g(tile))
end do

call mloimpice("u",i_u,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))
call mloimpice("v",i_v,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile))

return
end subroutine pack_coupled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update coupled

subroutine update_coupled(thetal,qvg,qlg,qfg,stratcloud,ua,va,    &
                          tlup,qvup,qlup,qfup,cfup,fg,eg,rhos,    &
                          ustar,cduv,ps,                          &    
                          tke,eps,mflx,fzzh,idzp,idzm,dz_hl,      &
                          rhoa1,dz_fl1,                           &
                          deptho_dz,deptho_dz_hl,                 &
                          rad_o,w_t,w_s,w_u,w_v,                  &
                          cd_water,cdh_water,cdbot_water,         &
                          wt0rad_o,wt0melt_o,                     &
                          wt0eg_o,icefg_a,wt0fb_o,ws0_o,          &
                          ws0subsurf_o,i_u,i_v,imass,fracice,     &
                          cd_ice,cdbot_ice,ibot,land,sigkap,      &
#ifdef scm
                          wthflux,wqvflux,uwflux,vwflux,          &
#endif
                          ddts,imax,kl,tile)

use mlo, only : wlev, mlo_updatekm, mlodiag, dgwater_g, wpack_g, wfull_g, &
                mlo_updatediag, water_g, turb_g, wrtemp, cp0, depth_g,    &
                cdbot, wrtrho, ice_g, mlocheck, minsfc, mlosalfix
                          
implicit none

integer, intent(in) :: imax, kl, tile
integer, dimension(imax), intent(in) :: ibot
integer k, kn, iq
real, dimension(imax,kl), intent(inout) :: thetal, qvg, qlg, qfg, stratcloud, ua, va
real, dimension(imax,kl), intent(in) :: tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: mflx, tke, eps
real, dimension(imax,kl), intent(in) :: idzm
real, dimension(imax,kl-1), intent(in) :: fzzh, idzp, dz_hl
real, dimension(imax,kl) :: km_a, kmo_a
real, dimension(imax,kl) :: bb_a, cc_a, dd_a, rr, tt_a
real, dimension(imax,2:kl) :: qq, aa_a
real, dimension(imax), intent(in) :: cd_water, cdh_water, cdbot_water, wt0rad_o, wt0melt_o
real, dimension(imax,wlev), intent(in) :: deptho_dz
real, dimension(imax,2:wlev), intent(in) :: deptho_dz_hl
real, dimension(imax,wlev), intent(in) :: rad_o
real, dimension(imax) :: t0_o, wt0eg_o
real, dimension(imax,wlev), intent(inout) :: w_t, w_s, w_u, w_v
real, dimension(imax,wlev) :: km_o, ks_o, rhs_o
real, dimension(imax,wlev) :: bb_o, cc_o, dd_o
real, dimension(imax,2:wlev) :: aa_o
real, dimension(imax,wlev) :: gammas_o
real, dimension(imax,1) :: bb_i, dd_i, tt_i
real, dimension(imax), intent(inout) :: fg, eg, ustar
real, dimension(imax), intent(in) :: rhos, rhoa1, dz_fl1, cduv, ps
real, dimension(imax), intent(in) :: icefg_a, wt0fb_o, ws0_o, ws0subsurf_o
real, dimension(imax), intent(in) :: fracice, imass, cd_ice, cdbot_ice
real, dimension(imax) :: wt0_a, wq0_a, f_ao, f_oa, f_ai, f_ia, f_oi, f_io
real, dimension(imax) :: wt0_o, wu0_o, wv0_o
real, dimension(imax) :: fulldeptho, targetdeptho, deltazo
real, dimension(imax), intent(inout) :: i_u, i_v
real, dimension(imax) :: t1, t3, t4
real, dimension(kl), intent(in) :: sigkap
real, intent(in) :: ddts
logical, dimension(imax), intent(in) :: land

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl) :: wthlflux, wqlflux, wqfflux
#endif

call mlocheck("Before coupled-mixing",water_temp=w_t,water_sal=w_s,water_u=w_u, &
              water_v=w_v,ice_u=i_u,ice_v=i_v)

! Calculate surface fluxes
wq0_a = eg/(rhos*lv)  ! qtot flux

! estimate eddy diffusivity mixing coeff for atmosphere
km_a(:,:) = cm0*tke(:,:)**2/eps(:,:)
call updatekmo(kmo_a,km_a,fzzh,imax,kl) ! interpolate diffusion coeffs to half levels

! Pre-calculate eddy diffiusivity mixing terms using updated kmo values
! -ve because gradient is calculated at t+1
qq(:,2:kl)  =-ddts*kmo_a(:,1:kl-1)*idzm(:,2:kl)/dz_hl(:,1:kl-1)
rr(:,1:kl-1)=-ddts*kmo_a(:,1:kl-1)*idzp(:,1:kl-1)/dz_hl(:,1:kl-1)

! k=1 is top of atmosphere and k=kl is bottom of atmosphere
cc_a(:,1) = qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
do k = 2,kl-1
  kn = kl - k + 1  
  aa_a(:,k) = rr(:,kn)-ddts*mflx(:,kn+1)*fzzh(:,kn)*idzp(:,kn)
  cc_a(:,k) = qq(:,kn)+ddts*mflx(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)
end do
aa_a(:,kl) = rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)


! calculate eddy diffusivity for ocean
km_o = 0.
ks_o = 0.
gammas_o = 0.
aa_o = 0.
bb_o = 1.
cc_o = 0.
dd_o = 0.
rhs_o = 0.

! update ocean eddy diffusivity possibly including prognostic tke and eps
call mlo_updatekm(km_o,ks_o,gammas_o,ddts,ps,0,depth_g(tile),ice_g(tile),dgwater_g(tile),water_g(tile), &
                  turb_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlodiag("t0",t0_o,0,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
where ( deptho_dz(:,1)>1.e-4 )
  rhs_o(:,1) = ks_o(:,2)*gammas_o(:,2)/deptho_dz(:,1)
end where
do k = 2,wlev-1
  where ( deptho_dz(:,k)>1.e-4 )  
    rhs_o(:,k) = (ks_o(:,k+1)*gammas_o(:,k+1)-ks_o(:,k)*gammas_o(:,k))/deptho_dz(:,k)
  end where  
end do
where ( deptho_dz(:,wlev)>1.e-4 )
  rhs_o(:,wlev) = -ks_o(:,wlev)*gammas_o(:,wlev)/deptho_dz(:,wlev)
end where

where ( deptho_dz(:,2)*deptho_dz(:,1)>1.e-4 )
  cc_o(:,1) = -ddts*ks_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
end where
do k = 2,wlev-1
  where ( deptho_dz(:,k-1)*deptho_dz(:,k)>1.e-4 )
    aa_o(:,k) = -ddts*ks_o(:,k)/(deptho_dz_hl(:,k)*deptho_dz(:,k))
  end where
  where ( deptho_dz(:,k+1)*deptho_dz(:,k)>1.e-4 )
    cc_o(:,k) = -ddts*ks_o(:,k+1)/(deptho_dz_hl(:,k+1)*deptho_dz(:,k))
  end where  
end do
where ( deptho_dz(:,wlev-1)*deptho_dz(:,wlev)>1.e-4 )
  aa_o(:,wlev) = -ddts*ks_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
end where


bb_i = 1.
dd_i = 0.


! zero coupling terms for land
f_ao(:) = 0.
f_oa(:) = 0.
f_ai(:) = 0.
f_ia(:) = 0.
f_io(:) = 0.
f_oi(:) = 0.


! thetal, thetao - coupled ----
t1 = rhos(:)*cdh_water*cp
! k=1 is top of atmosphere and k=kl is bottom of atmosphere
bb_a(:,1) = 1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)
dd_a(:,1) = thetal(:,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)   &
                              +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
do k = 2,kl-1
  kn = kl - k + 1  
  bb_a(:,k) = 1.-qq(:,kn)-rr(:,kn)+ddts*(mflx(:,kn)*fzzh(:,kn-1)*idzm(:,kn)             &
                                        -mflx(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn))
  dd_a(:,k) = thetal(:,kn)+ddts*(mflx(:,kn-1)*tlup(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn) &
                                +mflx(:,kn)*tlup(:,kn)*fzzh(:,kn-1)*idzm(:,kn)          &
                                -mflx(:,kn)*tlup(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn)       &
                                -mflx(:,kn+1)*tlup(:,kn+1)*fzzh(:,kn)*idzp(:,kn))
end do
bb_a(:,kl) = 1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)
dd_a(:,kl) = thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)             &
                              +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))
where ( land(1:imax) )
  wt0_a = fg/(rhos*cp)  ! theta flux    
  dd_a(:,kl) = dd_a(:,kl)+ddts*rhos*wt0_a/(rhoa1(:)*dz_fl1(:))
elsewhere
  bb_a(:,kl) = bb_a(:,kl)+ddts*t1(:)*(1.-fracice(:))/cp/(rhoa1(:)*dz_fl1(:)) 
  dd_a(:,kl) = dd_a(:,kl)-ddts*t1(:)*(1.-fracice(:))/cp*sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp/(rhoa1(:)*dz_fl1(:))
  dd_a(:,kl) = dd_a(:,kl)+ddts*t1(:)*(1.-fracice(:))*wrtemp/cp/(rhoa1(:)*dz_fl1(:))
  dd_a(:,kl) = dd_a(:,kl)+ddts*(fracice*icefg_a/cp)/(rhoa1(:)*dz_fl1(:))
end where  

bb_o(:,1) = 1. - cc_o(:,1)
dd_o(:,1) = w_t(:,1) + ddts*rhs_o(:,1)*t0_o
where ( deptho_dz(:,1)>1.e-4 )
  dd_o(:,1) = dd_o(:,1) - ddts*rad_o(:,1)/deptho_dz(:,1)
end where
do k = 2,wlev-1
  bb_o(:,k) = 1. - aa_o(:,k) - cc_o(:,k)
  dd_o(:,k) = w_t(:,k) + ddts*rhs_o(:,k)*t0_o
  where ( deptho_dz(:,k)>1.e-4 )
    dd_o(:,k) = dd_o(:,k) - ddts*rad_o(:,k)/deptho_dz(:,k)
  end where
end do
bb_o(:,wlev) = 1. - aa_o(:,wlev)
dd_o(:,wlev) = w_t(:,wlev) + ddts*rhs_o(:,wlev)*t0_o
where ( deptho_dz(:,wlev)>1.e-4 )
  dd_o(:,wlev) = dd_o(:,wlev) - ddts*rad_o(:,wlev)/deptho_dz(:,wlev)
end where
where ( .not.land(1:imax) )
  bb_o(:,1) = bb_o(:,1) + ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0*deptho_dz(:,1))
  dd_o(:,1) = dd_o(:,1) - ddts*t1(:)*(1.-fracice(:))*wrtemp/(wrtrho*cp0*deptho_dz(:,1))
  dd_o(:,1) = dd_o(:,1) - ddts*(wt0rad_o+wt0melt_o+wt0eg_o)/deptho_dz(:,1)
  dd_o(:,1) = dd_o(:,1) + ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0)*sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp/deptho_dz(:,1)
  dd_o(:,1) = dd_o(:,1) - ddts*wt0fb_o/deptho_dz(:,1)
end where  
where ( .not.land(1:imax) )
  f_ao(:) = -ddts*t1(:)*(1.-fracice(:))/cp/(rhoa1(:)*dz_fl1(:))
  f_oa(:) = -ddts*t1(:)*(1.-fracice(:))/(wrtrho*cp0*deptho_dz(:,1))
end where  
call solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a,tt_a,  &
                            aa_o,bb_o,cc_o,dd_o,w_t,   &
                            bb_i,dd_i,tt_i,            &
                            f_ao,f_oa,f_ai,f_ia,       &
                            f_oi,f_io,                 &
                            1,imax,kl,wlev) 
do k = 1,kl
  kn = kl - k + 1
  thetal(:,k) = tt_a(:,kn)
end do

! Derive sensible heat flux for water gridpoints
wt0_o = 0.
where ( .not.land(1:imax) )
  fg = (1.-fracice)*t1*(w_t(:,1)+wrtemp-thetal(:,1)-sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp) &
      + fracice*icefg_a
  wt0_o = wt0rad_o + wt0melt_o + wt0eg_o + wt0fb_o &
         + (1.-fracice)*t1*(w_t(:,1)-thetal(:,1)-sigkap(1)*(lv*qlg(:,1)+ls*qfg(:,1))/cp)/(wrtrho*cp0)
end where 
call mlo_updatediag("wt0",wt0_o,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))

#ifdef scm  
wthlflux(:,1)=fg/(rhos*cp)  ! theta flux
wthlflux(:,2:kl)=-kmo_a(:,1:kl-1)*(thetal(:,2:kl)-thetal(:,1:kl-1))/dz_hl(:,1:kl-1)    &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1)) &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif

bb_a(:,1) = 1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)
cc_a(:,1) = qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
do k = 2,kl-1
  kn = kl - k + 1  
  aa_a(:,k) = rr(:,kn)-ddts*mflx(:,kn+1)*fzzh(:,kn)*idzp(:,kn)
  
  bb_a(:,k) = 1.-qq(:,kn)-rr(:,kn)+ddts*(mflx(:,kn)*fzzh(:,kn-1)*idzm(:,kn)         &
                                        -mflx(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn))
  cc_a(:,k) = qq(:,kn)+ddts*mflx(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)
end do
aa_a(:,kl) = rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)
bb_a(:,kl) = 1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)


where ( deptho_dz(:,2)*deptho_dz(:,1)>1.e-4 )
  cc_o(:,1) = -ddts*ks_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
end where
bb_o(:,1) = 1. - cc_o(:,1)
do k = 2,wlev-1
  where ( deptho_dz(:,k-1)*deptho_dz(:,k)>1.e-4 )  
    aa_o(:,k) = -ddts*ks_o(:,k)/(deptho_dz_hl(:,k)*deptho_dz(:,k)) 
  end where
  where ( deptho_dz(:,k+1)*deptho_dz(:,k)>1.e-4 )
    cc_o(:,k) = -ddts*ks_o(:,k+1)/(deptho_dz_hl(:,k+1)*deptho_dz(:,k))
  end where
  bb_o(:,k) = 1. - aa_o(:,k) - cc_o(:,k)
end do
where ( deptho_dz(:,wlev-1)*deptho_dz(:,wlev)>1.e-4 )
  aa_o(:,wlev) = -ddts*ks_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
end where
bb_o(:,wlev) = 1. - aa_o(:,wlev)


! Note that vertical interpolation is linear so that qtot can be
! decomposed into qv, ql and qf.

! qv (part of qtot) - atmosphere
dd_a(:,1)=qvg(:,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                         +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
do k = 2,kl-1
  kn = kl - k + 1
  dd_a(:,k)=qvg(:,kn)+ddts*(mflx(:,kn-1)*qvup(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)   &
                           +mflx(:,kn)*qvup(:,kn)*fzzh(:,kn-1)*idzm(:,kn)            &
                           -mflx(:,kn)*qvup(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn)         &
                           -mflx(:,kn+1)*qvup(:,kn+1)*fzzh(:,kn)*idzp(:,kn))
end do
dd_a(:,kl)=qvg(:,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                         +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                   &
                   +ddts*rhos*wq0_a/(rhoa1(:)*dz_fl1(:))
call thomas(tt_a,aa_a(:,2:kl),bb_a(:,1:kl),cc_a(:,1:kl-1),dd_a(:,1:kl),imax,kl)
do k = 1,kl
  kn = kl - k + 1  
  qvg(:,k) = tt_a(:,kn)
end do

#ifdef scm
wqvflux(:,1)=wq0_a(:)
wqvflux(:,2:kl)=-kmo_a(:,1:kl-1)*(qvg(:,2:kl)-qvg(:,1:kl-1))/dz_hl(:,1:kl-1)         &
                +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))   &
                +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif

! ql (part of qtot) - atmosphere
dd_a(:,1)=qlg(:,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                         +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
do k = 2,kl-1
  kn = kl - k + 1
  dd_a(:,k)=qlg(:,kn)+ddts*(mflx(:,kn-1)*qlup(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)   &
                           +mflx(:,kn)*qlup(:,kn)*fzzh(:,kn-1)*idzm(:,kn)            &
                           -mflx(:,kn)*qlup(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn)         &
                           -mflx(:,kn+1)*qlup(:,kn+1)*fzzh(:,kn)*idzp(:,kn))
end do
dd_a(:,kl)=qlg(:,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                         +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
call thomas(tt_a,aa_a(:,2:kl),bb_a(:,1:kl),cc_a(:,1:kl-1),dd_a(:,1:kl),imax,kl)
do k = 1,kl
  kn = kl - k + 1  
  qlg(:,k) = tt_a(:,kn)
end do

#ifdef scm
wqlflux(:,1)=0.
wqlflux(:,2:kl)=-kmo_a(:,1:kl-1)*(qlg(:,2:kl)-qlg(:,1:kl-1))/dz_hl(:,1:kl-1)         &
                +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))   &
                +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif

! qf (part of qtot) - atmosphere
dd_a(:,1)=qfg(:,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                         +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
do k = 2,kl-1
  kn = kl - k + 1
  dd_a(:,k)=qfg(:,kn)+ddts*(mflx(:,kn-1)*qfup(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)   &
                           +mflx(:,kn)*qfup(:,kn)*fzzh(:,kn-1)*idzm(:,kn)            &
                           -mflx(:,kn)*qfup(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn)         &
                           -mflx(:,kn+1)*qfup(:,kn+1)*fzzh(:,kn)*idzp(:,kn))
end do
dd_a(:,kl)=qfg(:,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                         +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
call thomas(tt_a,aa_a(:,2:kl),bb_a(:,1:kl),cc_a(:,1:kl-1),dd_a(:,1:kl),imax,kl)
do k = 1,kl
  kn = kl - k + 1  
  qfg(:,k) = tt_a(:,kn)
end do

#ifdef scm
wqfflux(:,1)=0.
wqfflux(:,2:kl)=-kmo_a(:,1:kl-1)*(qfg(:,2:kl)-qfg(:,1:kl-1))/dz_hl(:,1:kl-1)                &
                +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))          &
                +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif

! stratcloud - atmosphere
dd_a(:,1)=stratcloud(:,kl)+ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)     &
                                +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
do k = 2,kl-1
  kn = kl - k + 1
  dd_a(:,k)=stratcloud(:,kn)+ddts*(mflx(:,kn-1)*cfup(:,kn-1)*(1.-fzzh(:,kn-1))*idzm(:,kn)   &
                                  +mflx(:,kn)*cfup(:,kn)*fzzh(:,kn-1)*idzm(:,kn)            &
                                  -mflx(:,kn)*cfup(:,kn)*(1.-fzzh(:,kn))*idzp(:,kn)         &
                                  -mflx(:,kn+1)*cfup(:,kn+1)*fzzh(:,kn)*idzp(:,kn))
end do
dd_a(:,kl)=stratcloud(:,1)-ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)               &
                                +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))
call thomas(tt_a,aa_a(:,2:kl),bb_a(:,1:kl),cc_a(:,1:kl-1),dd_a(:,1:kl),imax,kl)
do k = 1,kl
  kn = kl - k + 1  
  stratcloud(:,k) = min( max( tt_a(:,kn), 0. ), 1. )
end do

! sal - ocean
fulldeptho = 0.
targetdeptho = max( min( minsfc, sum( deptho_dz, dim=2 ) ), 1.e-4 )
do k = 1,wlev
  dd_o(:,k) = w_s(:,k) + ddts*rhs_o(:,k)*ws0_o
  where ( deptho_dz(:,k)>1.e-4 )
    deltazo = max( min( deptho_dz(:,k), targetdeptho-fulldeptho ), 0. )
    fulldeptho = fulldeptho + deltazo
    dd_o(:,k) = dd_o(:,k) - ddts*ws0subsurf_o*deltazo/(deptho_dz(:,k)*targetdeptho)
  end where
end do
where ( .not.land(1:imax) )
  dd_o(:,1) = dd_o(:,1) - ddts*ws0_o/deptho_dz(:,1)
end where
call thomas(w_s,aa_o(:,2:wlev),bb_o(:,1:wlev),cc_o(:,1:wlev-1),dd_o(:,1:wlev),imax,wlev)

! momentum ----

t1 = rhos(:)*cd_water
t3 = rhos(:)*cd_ice
t4 = wrtrho*cdbot_ice

bb_a(:,1) = 1.-qq(:,kl)
cc_a(:,1) = qq(:,kl)
do k = 2,kl-1
  kn = kl - k + 1  
  aa_a(:,k) = rr(:,kn)
  bb_a(:,k) = 1.-qq(:,kn)-rr(:,kn)
  cc_a(:,k) = qq(:,kn)
end do
aa_a(:,kl) = rr(:,1)
bb_a(:,kl) = 1.-rr(:,1)
where ( land(1:imax) )
  bb_a(:,kl) = bb_a(:,kl) + ddts*rhos*cduv/(rhoa1(:)*dz_fl1(:)) ! implicit  
elsewhere
  bb_a(:,kl) = bb_a(:,kl) + ddts*(t1(:)*(1.-fracice(:))+t3(:)*fracice(:))/(rhoa1(:)*dz_fl1(:))
end where  

where ( deptho_dz(:,2)*deptho_dz(:,1) > 1.e-4 )
  cc_o(:,1) = -ddts*km_o(:,2)/(deptho_dz_hl(:,2)*deptho_dz(:,1))
end where
bb_o(:,1) = 1. - cc_o(:,1)
do k = 2,wlev-1
  where ( deptho_dz(:,k-1)*deptho_dz(:,k) > 1.e-4 )  
    aa_o(:,k) = -ddts*km_o(:,k)/(deptho_dz_hl(:,k)*deptho_dz(:,k)) 
  end where
  where ( deptho_dz(:,k+1)*deptho_dz(:,k) > 1.e-4 )
    cc_o(:,k) = -ddts*km_o(:,k+1)/(deptho_dz_hl(:,k+1)*deptho_dz(:,k))
  end where
  bb_o(:,k) = 1. - aa_o(:,k) - cc_o(:,k)
end do
where ( deptho_dz(:,wlev-1)*deptho_dz(:,wlev) > 1.e-4 )
  aa_o(:,wlev) = -ddts*km_o(:,wlev)/(deptho_dz_hl(:,wlev)*deptho_dz(:,wlev))
end where
bb_o(:,wlev) = 1. - aa_o(:,wlev)
where ( .not.land(1:imax) )
  bb_o(:,1) = bb_o(:,1) + ddts*(t1(:)*(1.-fracice(:))+t4(:)*fracice(:))/(wrtrho*deptho_dz(:,1))
end where
! bottom drag
do iq = 1,imax
  if ( .not.land(iq) ) then  
    k = ibot(iq)
    bb_o(iq,k) = bb_o(iq,k) + ddts*cdbot_water(iq)/deptho_dz(iq,k) 
  end if
end do

bb_i(:,1) = 1.
where ( .not.land(1:imax) )
  bb_i(:,1) = 1. + ddts*(t3(:)+t4(:))/imass
end where

where ( .not.land(1:imax) )
  f_ao(:) = -ddts*t1(:)*(1.-fracice(:))/(rhoa1(:)*dz_fl1(:))
  f_oa(:) = -ddts*t1(:)*(1.-fracice(:))/(wrtrho*deptho_dz(:,1))
  f_ai(:) = -ddts*t3(:)*fracice(:)/(rhoa1(:)*dz_fl1(:))
  f_ia(:) = -ddts*t3(:)/imass
  f_oi(:) = -ddts*t4(:)*fracice(:)/(wrtrho*deptho_dz(:,1))
  f_io(:) = -ddts*t4(:)/imass
end where  

! ua, uo - coupled ----
! k=1 is top of atmosphere and k=kl is bottom of atmosphere
do k = 1,kl
  kn = kl - k + 1  
  dd_a(:,k) = ua(:,kn)
end do
do k = 1,wlev
  dd_o(:,k) = w_u(:,k)
end do
dd_i(:,1) = 0.
where ( .not.land(1:imax) )
  dd_i(:,1) = i_u(:)
end where
call solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a,tt_a,  &
                            aa_o,bb_o,cc_o,dd_o,w_u,   &
                            bb_i,dd_i,tt_i,            &
                            f_ao,f_oa,f_ai,f_ia,       &
                            f_oi,f_io,                 &
                            2,imax,kl,wlev) 
do k = 1,kl
  kn = kl - k + 1
  ua(:,k) = tt_a(:,kn)
end do
i_u(:) = tt_i(:,1)

! va, vo - coupled ----
! k=1 is top of atmosphere and k=kl is bottom of atmosphere
do k = 1,kl
  kn = kl - k + 1  
  dd_a(:,k) = va(:,kn)
end do
do k = 1,wlev
  dd_o(:,k) = w_v(:,k)
end do
dd_i(:,1) = 0.
where ( .not.land(1:imax) )
  dd_i(:,1) = i_v(:)
end where
call solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a,tt_a,  &
                            aa_o,bb_o,cc_o,dd_o,w_v,   &
                            bb_i,dd_i,tt_i,            &
                            f_ao,f_oa,f_ai,f_ia,       &
                            f_oi,f_io,                 &
                            2,imax,kl,wlev) 
do k = 1,kl
  kn = kl - k + 1
  va(:,k) = tt_a(:,kn)
end do
i_v(:) = tt_i(:,1)

! update surface momentum flux
where ( land(1:imax) )
  ustar = sqrt(cduv*sqrt(ua(:,1)**2+va(:,1)**2))
elsewhere
  ustar = sqrt((1.-fracice)*cd_water*sqrt((ua(:,1)-w_u(:,1))**2+(va(:,1)-w_v(:,1))**2)  &
              +fracice*cd_ice*sqrt((ua(:,1)-i_u(:))**2+(va(:,1)-i_v(:))**2))
end where

! update wu0 and wv0
wu0_o = -(1.-fracice)*rhoa1*cd_water*(ua(:,1)-w_u(:,1))/wrtrho
wv0_o = -(1.-fracice)*rhoa1*cd_water*(va(:,1)-w_v(:,1))/wrtrho
wu0_o = wu0_o + fracice*cdbot_ice*(w_u(:,1)-i_u(:))
wv0_o = wv0_o + fracice*cdbot_ice*(w_v(:,1)-i_v(:))
call mlo_updatediag("wu0",wu0_o,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))
call mlo_updatediag("wv0",wv0_o,0,dgwater_g(tile),wpack_g(:,tile),wfull_g(tile))

call mlosalfix(w_s)
call mlocheck("Coupled-mixing",water_temp=w_t,water_sal=w_s,water_u=w_u, &
              water_v=w_v,ice_u=i_u,ice_v=i_v)

#ifdef scm
uwflux(:,1) = 0.
uwflux(:,2:kl) = -kmo_a(:,1:kl-1)*(ua(:,2:kl)-ua(:,1:kl-1))/dz_hl(:,1:kl-1)
vwflux(:,1) = 0.
vwflux(:,2:kl) = -kmo_a(:,1:kl-1)*(va(:,2:kl)-va(:,1:kl-1))/dz_hl(:,1:kl-1)
do k = 1,kl
  wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                               *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
end do
#endif

return
end subroutine update_coupled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve Sherman-Morrison formula

subroutine solve_sherman_morrison(aa_a,bb_a,cc_a,dd_a,tt_a, &
                                  aa_o,bb_o,cc_o,dd_o,tt_o, &
                                  bb_i,dd_i,tt_i,           &
                                  f_ao,f_oa,f_ai,f_ia,      &
                                  f_oi,f_io,                &
                                  aoi_mode,imax,kl,wlev)

implicit none

integer, intent(in) :: imax, kl, wlev, aoi_mode
integer k
real, dimension(imax,2:kl), intent(in) :: aa_a
real, dimension(imax,kl), intent(in) :: bb_a, dd_a
real, dimension(imax,kl-1), intent(in) :: cc_a
real, dimension(imax,kl), intent(out) :: tt_a
real, dimension(imax,2:wlev), intent(in) :: aa_o
real, dimension(imax,wlev), intent(in) :: bb_o, dd_o
real, dimension(imax,wlev-1), intent(in) :: cc_o
real, dimension(imax,wlev), intent(out) :: tt_o
real, dimension(imax,1), intent(in) :: bb_i, dd_i
real, dimension(imax,1), intent(out) :: tt_i
real, dimension(imax,2:kl+wlev+1) :: aad
real, dimension(imax,kl+wlev+1) :: bbd, ddd
real, dimension(imax,kl+wlev) :: ccd
real, dimension(imax,kl+wlev+1) :: uu, vv, yy, qq
real, dimension(imax,kl+wlev+1) :: xx
real, dimension(imax), intent(in) :: f_ao, f_oa, f_ai, f_ia, f_oi, f_io
real, dimension(imax) :: t1, t2

! Original coupling matrix (inverted levels)
! [ bb_a cc_a                                                   ] [ t_a ]   [ dd_a ]
! [ aa_a bb_a cc_a                                              ] [ t_a ]   [ dd_a ]
! [      .... .... ....                                         ] [ ... ]   [ .... ]
! [           aa_a bb_a f_ai                f_ao                ] [ t_a ]   [ dd_a ]
! [                f_ia bb_i cc_i                               ] [ t_i ]   [ dd_i ]
! [                     aa_i bb_i cc_i                          ] [ t_i ]   [ dd_i ]
! [                          .... .... ....                     ] [ ... ] = [ .... ]
! [                               aa_i bb_i f_io                ] [ t_i ]   [ dd_i ]
! [                f_oa                f_oi bb_o cc_o           ] [ t_o ]   [ dd_o ]
! [                                         aa_o bb_o cc_o      ] [ t_o ]   [ dd_o ]
! [                                              .... .... .... ] [ ... ]   [ .... ]
! [                                                   aa_o bb_o ] [ t_o ]   [ dd_o ]

! u^T = [ .. -bb_a ..  f_oa      .. ]   
! v^T = [ .. 1     .. -f_ao/bb_a .. ]


! Improved version (scalar)
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
! x = y - {(v^t y)/(1 + (v^t q))} q


! Construct A'
bbd(:,1) = bb_a(:,1)
ccd(:,1) = cc_a(:,1)
ddd(:,1) = dd_a(:,1)
do k = 2,kl-1
  aad(:,k) = aa_a(:,k)
  bbd(:,k) = bb_a(:,k)
  ccd(:,k) = cc_a(:,k)
  ddd(:,k) = dd_a(:,k)
end do
aad(:,kl) = aa_a(:,kl)
bbd(:,kl) = bb_a(:,kl)  + f_ai*f_ia ! due to u v^T matrix
ccd(:,kl) = f_ao        + f_ai*f_io ! due to u v^T matrix
ddd(:,kl) = dd_a(:,kl)
aad(:,kl+1) = f_oa      + f_oi*f_ia ! due to u v^T matrix
bbd(:,kl+1) = bb_o(:,1) + f_oi*f_io ! due to u v^T matrix
ccd(:,kl+1) = cc_o(:,1)
ddd(:,kl+1) = dd_o(:,1)
do k = 2,wlev-1
  aad(:,kl+k) = aa_o(:,k)
  bbd(:,kl+k) = bb_o(:,k)
  ccd(:,kl+k) = cc_o(:,k)
  ddd(:,kl+k) = dd_o(:,k)
end do
aad(:,kl+wlev) = aa_o(:,wlev)
bbd(:,kl+wlev) = bb_o(:,wlev)
ccd(:,kl+wlev) = 0.
ddd(:,kl+wlev) = dd_o(:,wlev)
aad(:,kl+wlev+1) = 0.
bbd(:,kl+wlev+1) = bb_i(:,1) + 1. ! note +1. due to u v^T matrix
ddd(:,kl+wlev+1) = dd_i(:,1)

! Construct u and v
uu(:,:) = 0.
uu(:,kl) = f_ai
uu(:,kl+1) = f_oi
uu(:,kl+wlev+1) = -1.
vv(:,:) = 0.
vv(:,kl) = -f_ia
vv(:,kl+1) = -f_io
vv(:,kl+wlev+1) = 1.

select case(aoi_mode)
  case(0)
    ! pure tridiagonal matrix for atmosphere only
    xx(:,:) = 0.  
    call thomas(xx(:,1:kl),aad(:,2:kl),bbd(:,1:kl),ccd(:,1:kl-1),ddd(:,1:kl),imax,kl)
  case(1)
    ! pure tridiagonal matrix for atmosphere and ocean (no sea-ice)
    xx(:,:) = 0.
    call thomas(xx(:,1:kl+wlev),aad(:,2:kl+wlev),bbd(:,1:kl+wlev),ccd(:,1:kl+wlev-1),ddd(:,1:kl+wlev),imax,kl+wlev)
  case default
    ! Solve A' y = d  and  A' q = u for combined atmosphere, ocean and sea-ice
    call thomas(yy(:,1:kl+wlev+1),aad(:,2:kl+wlev+1),bbd(:,1:kl+wlev+1),ccd(:,1:kl+wlev),ddd(:,1:kl+wlev+1),imax,kl+wlev+1)
    call thomas(qq(:,1:kl+wlev+1),aad(:,2:kl+wlev+1),bbd(:,1:kl+wlev+1),ccd(:,1:kl+wlev),uu(:,1:kl+wlev+1),imax,kl+wlev+1)
    ! Solve for x = y - {(v^t y)/(1 + (v^t q))} q
    t1(:) = vv(:,kl)*yy(:,kl) + vv(:,kl+1)*yy(:,kl+1) + vv(:,kl+wlev+1)*yy(:,kl+wlev+1)
    t2(:) = vv(:,kl)*qq(:,kl) + vv(:,kl+1)*qq(:,kl+1) + vv(:,kl+wlev+1)*qq(:,kl+wlev+1)
    do k = 1,kl+wlev+1
      xx(:,k) = yy(:,k) - qq(:,k)*t1(:)/(1.+t2(:))
    end do  
end select

! Unpack solution
do k = 1,kl
  tt_a(:,k) = xx(:,k)
end do
do k = 1,wlev
  tt_o(:,k) = xx(:,kl+k)
end do
tt_i(:,1) = xx(:,kl+wlev+1)

return
end subroutine solve_sherman_morrison
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

pure subroutine thomas(outdat,aai,bbi,cci,ddi,imax,klin)

implicit none

integer, intent(in) :: imax, klin
real, dimension(imax,2:klin), intent(in) :: aai
real, dimension(imax,klin), intent(in) :: bbi,ddi
real, dimension(imax,klin-1), intent(in) :: cci
real, dimension(imax,klin), intent(out) :: outdat
real, dimension(imax,klin) :: cc,dd
real :: n_s
integer k, iq

cc(:,1)=cci(:,1)/bbi(:,1)
dd(:,1)=ddi(:,1)/bbi(:,1)

do k=2,klin-1
  do iq = 1,imax
    n_s=bbi(iq,k)-cc(iq,k-1)*aai(iq,k)
    cc(iq,k)=cci(iq,k)/n_s
    dd(iq,k)=(ddi(iq,k)-dd(iq,k-1)*aai(iq,k))/n_s
  end do
end do
do iq = 1,imax
  n_s=bbi(iq,klin)-cc(iq,klin-1)*aai(iq,klin)
  dd(iq,klin)=(ddi(iq,klin)-dd(iq,klin-1)*aai(iq,klin))/n_s
  outdat(iq,klin)=dd(iq,klin)
end do
do k=klin-1,1,-1
  outdat(:,k)=dd(:,k)-cc(:,k)*outdat(:,k+1)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

pure subroutine getqsat(qsat,templ,ps,fice,imax)

implicit none

integer, intent(in) :: imax
real, dimension(imax), intent(in) :: templ
real, dimension(imax), intent(in) :: ps, fice
real, dimension(imax), intent(out) :: qsat
real, dimension(imax) :: estafi, tdiff, tdiffx, rx, rxx
real, dimension(imax) :: qsatl, qsati, deles
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

pure subroutine updatekmo(kmo,km,fzhl,imax,kl)

implicit none

integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(out) :: kmo
real, dimension(imax,kl), intent(in) :: km
real, dimension(imax,kl-1), intent(in) :: fzhl

kmo(:,1:kl-1)=km(:,1:kl-1)+fzhl(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))
! These terms are never used
kmo(:,kl)=0.

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

pure function entfn(zht,zi,imax) result(ans)

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

subroutine calc_phi(phim,z_on_l,imax)

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

