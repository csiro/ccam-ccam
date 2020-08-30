! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public tkemeth,upshear

real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps

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
!$acc declare create(cm0,minl,maxl,mintke,mineps,ce0,ce1,ce2,ce3)
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
!$acc declare create(be,ent0,ent1,ent_min,ezmin,entc0,dtrc0,m0,b1,b2,qcmf,mfbeta)
! generic constants
integer, save :: buoymeth = 1     ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0     ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: tkemeth  = 1     ! Method for TKE calculation (0=D&K84, 1=Hurley)
integer, save :: upshear  = 0     ! Method for updating shear (0=Host, 1=Sub-timestep)
real, save :: maxdts      = 120.  ! max timestep for split
!$acc declare create(buoymeth,tkemeth,stabmeth,maxdts,upshear)

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

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(shear(ifull,kl))

tke=mintke
eps=mineps
shear=0.

!$acc update device(cm0,minl,maxl,mintke,mineps,ce0,ce1,ce2,ce3)
!$acc update device(be,ent0,ent1,ent_min,ezmin,entc0,dtrc0,m0,b1,b2,qcmf,mfbeta)
!$acc update device(buoymeth,tkemeth,stabmeth,maxdts,upshear)

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,stratcloud,uo,vo,zi,fg,eg,cduv,ps,zz,zzh,sig,rhos, &
                  ustar_ave,dt,qgmin,mode,tke,eps,shear,dx                                 &
#ifdef scm
                  ,wthflux,wqvflux,uwflux,vwflux,mfout,buoyproduction                      &
                  ,shearproduction,totaltransport                                          &
#endif
                  ,imax,kl)
!$acc routine vector
                  
implicit none

integer, intent(in) :: imax, kl, mode
integer k, iq
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(imax,kl), intent(inout) :: theta,stratcloud,uo,vo
real, dimension(imax,kl), intent(inout) :: qvg,qlg,qfg
real, dimension(imax,kl), intent(out) :: kmo
real, dimension(imax,kl), intent(in) :: zz,zzh
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(inout) :: eps
real, dimension(imax), intent(inout) :: zi
real, dimension(imax), intent(in) :: fg,eg,cduv,ps,rhos,dx
real, dimension(imax), intent(out) :: ustar_ave
real, dimension(imax), intent(in) :: sig
real, dimension(imax,kl) :: km,thetav,thetal,qsat
real, dimension(imax,kl) :: qsatc,qgnc,ff
real, dimension(imax,kl) :: thetalhl,thetavhl,uo_hl,vo_hl
real, dimension(imax,kl) :: quhl,qshl,qlhl,qfhl
real, dimension(imax,kl) :: bb,cc,dd,rr
real, dimension(imax,kl) :: rhoa,rhoahl
real, dimension(imax,kl) :: qtot,qthl
real, dimension(imax,kl) :: tlup,qvup,qlup,qfup
real, dimension(imax,kl) :: cfup,mflx
real, dimension(imax,kl) :: pps,ppt,ppb
real, dimension(imax,kl) :: idzm
real, dimension(imax,kl-1) :: idzp
real, dimension(imax,2:kl) :: aa,qq
real, dimension(imax,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(imax,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(imax,kl-1) :: fzzh
real, dimension(imax) :: wt0,wq0,wtv0
real, dimension(imax) :: wstar,z_on_l,phim
real, dimension(imax) :: pres,temp
real :: tff,tempc,thetac,tempt
real, dimension(imax) :: ustar
real :: tbb,tqq,tcc,tempv,fc,dc,mc,bvf,rvar
real, dimension(imax) :: zi_save, zturb, cgmap
real, dimension(imax) :: templ, dqsdt, al, fice, qc, qt, lx
real, dimension(kl) :: sigkap
real cm12, cm34, ddts
logical, dimension(imax,kl) :: lta
logical, dimension(imax) :: mask

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl), intent(out) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(out) :: totaltransport
real, dimension(imax,kl), intent(out) :: mfout
real, dimension(imax,kl) :: wthlflux, wqlflux
real, dimension(imax,kl) :: wqfflux
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
    tff      = cm34*tke(iq,k)*sqrt(tke(iq,k))
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
call updatekmo(kmo,   km,  fzzh,imax,kl)
call updatekmo(rhoahl,rhoa,fzzh,imax,kl)
! eddy diffusion terms to account for air density with level thickness
idzm(:,2:kl)   = rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1) = rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

ustar_ave(:) = 0.

! Main loop to prevent time splitting errors
mcount = int(dt/(min(maxdts,12./max(m0,0.1))+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  ! Set-up thermodynamic variables temp, theta_v, qtot and saturated mixing ratio
  do k = 1,kl
    templ(:) = thetal(:,k)/sigkap(k)
    pres(:) = ps(:)*sig(k)
    ! calculate saturated air mixing ratio
    where ( qfg(:,k)>1.e-12 )
      fice = min( qfg(:,k)/(qfg(:,k)+qlg(:,k)), 1. )
    elsewhere
      fice = 0.
    end where
    call getqsat(qsat(:,k),templ(:),pres(:),fice,imax)
    thetav(:,k) = theta(:,k)*(1.+0.61*qvg(:,k)-qlg(:,k)-qfg(:,k))
    qtot(:,k) = qvg(:,k) + qlg(:,k) + qfg(:,k)
  end do
  
  ! Update thetav flux
  wtv0 = wt0 + theta(:,1)*0.61*wq0 ! thetav flux
  
  ! Update momentum flux
  ustar = sqrt(cduv*sqrt(uo(:,1)**2+vo(:,1)**2))  
  wstar = (grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)   
  
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

  if ( mode/=1 ) then ! mass flux
      
    zi_save = zi  

    ! plume rise model
    mask = wtv0>0.
    tke(:,1) = cm12*ustar**2 + ce3*wstar**2
    tke(:,1) = max(tke(:,1), mintke)
    call plumerise(mask,                                         &
                   zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,       &
                   zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,     &
                   stratcloud,km,wt0,wq0,wtv0,ps,ustar,          &
                   sig,sigkap,tke,eps,uo,vo,imax,kl)

#ifndef scm
    ! Turn off MF term if small grid spacing (mfbeta=0 implies MF is always non-zero)
    ! Based on Boutle et al 2014
    zturb = 0.5*(zi_save + zi)
    cgmap(:) = 1. - tanh(mfbeta*zturb/dx)*max(0.,1.-0.25*dx/zturb)
    do k = 1,kl
      mflx(:,k) = mflx(:,k)*cgmap(:)
    end do
#endif
    
  end if

#ifdef scm  
  do k = 1,kl-1
    mfout(:,k) = mflx(:,k)*(1.-fzzh(:,k)) &
               + mflx(:,k+1)*fzzh(:,k)
  end do  
#endif
  
  
  ! calculate tke and eps boundary condition at 1st vertical level
  z_on_l = -vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-10))
  z_on_l = min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  call calc_phi(phim,z_on_l,imax)
  do iq = 1,imax
    tke(iq,1) = cm12*ustar(iq)*ustar(iq)+ce3*wstar(iq)*wstar(iq)
    eps(iq,1) = ustar(iq)*ustar(iq)*ustar(iq)*phim(iq)/(vkar*zz(iq,1))+grav*wtv0(iq)/thetav(iq,1)
    tke(iq,1) = max( tke(iq,1), mintke )
    tff = cm34*tke(iq,1)*sqrt(tke(iq,1))
    eps(iq,1) = min( eps(iq,1), tff/minl )
    eps(iq,1) = max( eps(iq,1), tff/maxl, mineps )
  end do


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
      do k=1,kl
        do iq =1,imax
          tbb=max(1.-stratcloud(iq,k),1.E-8)
          qgnc(iq,k)=(qvg(iq,k)-(1.-tbb)*qsatc(iq,k))/tbb    ! outside cloud value
          qgnc(iq,k)=min(max(qgnc(iq,k),qgmin),qsatc(iq,k))
        end do
      end do
      call updatekmo(thetalhl,thetal,fzzh,imax,kl)              ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh,imax,kl)                    ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh,imax,kl)                   ! inside cloud value
      call updatekmo(qlhl,dd,fzzh,imax,kl)                      ! inside cloud value
      call updatekmo(qfhl,ff,fzzh,imax,kl)                      ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl)=stratcloud(:,2:kl)<=1.E-6
      do k=2,kl-1
        where(lta(:,k).and..not.lta(:,k+1))
          qlhl(:,k)=dd(:,k+1)
          qfhl(:,k)=ff(:,k+1)
        elsewhere (.not.lta(:,k).and.lta(:,k+1))
          qlhl(:,k)=dd(:,k)
          qfhl(:,k)=ff(:,k)
        end where
      end do
      do k=2,kl-1
        ! saturated
        do iq = 1,imax
          thetac=thetal(iq,k)+sigkap(k)*(lv*dd(iq,k)+ls*ff(iq,k))/cp              ! inside cloud value
          tempc=thetac/sigkap(k)                                            ! inside cloud value          
          tqq=(1.+lv*qsatc(iq,k)/(rd*tempc))/(1.+lv*lv*qsatc(iq,k)/(cp*rv*tempc*tempc))
          tbb=-grav*km(iq,k)*(tqq*((thetalhl(iq,k)-thetalhl(iq,k-1)+sigkap(k)/cp*(lv*(qlhl(iq,k)-qlhl(iq,k-1))  &
              +ls*(qfhl(iq,k)-qfhl(iq,k-1))))/thetac+lv/cp*(qshl(iq,k)-qshl(iq,k-1))/tempc)              &
              -qshl(iq,k)-qlhl(iq,k)-qfhl(iq,k)+qshl(iq,k-1)+qlhl(iq,k-1)+qfhl(iq,k-1))/dz_fl(iq,k)
          ! unsaturated
          tcc=-grav*km(iq,k)*(thetalhl(iq,k)-thetalhl(iq,k-1)+thetal(iq,k)*0.61*(quhl(iq,k)-quhl(iq,k-1)))  &
                           /(thetal(iq,k)*dz_fl(iq,k))
          ppb(iq,k)=(1.-stratcloud(iq,k))*tcc+stratcloud(iq,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)
        end do
      end do
      ! saturated
      do iq = 1,imax
        thetac=thetal(iq,1)+sigkap(1)*(lv*dd(iq,1)+ls*ff(iq,1))/cp              ! inside cloud value
        tempc=thetac/sigkap(1)                                            ! inside cloud value          
        tqq=(1.+lv*qsatc(iq,1)/(rd*tempc))/(1.+lv*lv*qsatc(iq,1)/(cp*rv*tempc*tempc))
        tbb=-grav*km(iq,1)*(tqq*((thetalhl(iq,1)-thetal(iq,1)+sigkap(1)/cp*(lv*(qlhl(iq,1)-qlg(iq,1))         &
            +ls*(qfhl(iq,1)-qfg(iq,1))))/thetac+lv/cp*(qshl(iq,1)-qsatc(iq,1))/tempc)                  &
            -qshl(iq,1)-qlhl(iq,1)-qfhl(iq,1)+qsatc(iq,1)+qlg(iq,1)+qfg(iq,1))/(zzh(iq,1)-zz(iq,1))
        ! unsaturated
        tcc=-grav*km(iq,1)*(thetalhl(iq,1)-thetal(iq,1)+thetal(iq,1)*0.61*(quhl(iq,1)-qgnc(iq,1)))        &
                       /(thetal(iq,1)*(zzh(iq,1)-zz(iq,1)))
        ppb(iq,1)=(1.-stratcloud(iq,1))*tcc+stratcloud(iq,1)*tbb ! cloud fraction weighted (e.g., Smith 1990)
      end do

      
    case(1) ! Marquet and Geleyn QJRMS (2012) for partially saturated
      call updatekmo(thetalhl,thetal,fzzh,imax,kl)
      call updatekmo(qthl,qtot,fzzh,imax,kl)
      do k=2,kl-1
        do iq = 1,imax
          tempt  = theta(iq,k)/sigkap(k)
          tempv = thetav(iq,k)/sigkap(k)
          rvar=rd*tempv/tempt ! rvar = qd*rd+qv*rv
          fc=(1.-stratcloud(iq,k))+stratcloud(iq,k)*(lv*rvar/(cp*rv*tempt))
          dc=(1.+0.61*qvg(iq,k))*lv*qvg(iq,k)/(rd*tempv)
          mc=(1.+dc)/(1.+lv*qlg(iq,k)/(cp*tempt)+dc*fc)
          bvf=grav*mc*(thetalhl(iq,k)-thetalhl(iq,k-1))/(thetal(iq,k)*dz_fl(iq,k))           &
             +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(iq,k)-qthl(iq,k-1))/dz_fl(iq,k)
          ppb(iq,k)=-km(iq,k)*bvf
        end do
      end do
      do iq = 1,imax
        tempt  = theta(iq,1)/sigkap(1)
        tempv = thetav(iq,1)/sigkap(1)
        rvar=rd*tempv/tempt ! rvar = qd*rd+qv*rv
        fc=(1.-stratcloud(iq,1))+stratcloud(iq,1)*(lv*rvar/(cp*rv*tempt))
        dc=(1.+0.61*qvg(iq,1))*lv*qvg(iq,1)/(rd*tempv)
        mc=(1.+dc)/(1.+lv*qlg(iq,1)/(cp*tempt)+dc*fc)
        bvf=grav*mc*(thetalhl(iq,1)-thetal(iq,1))/(thetal(iq,1)*(zzh(iq,1)-zz(iq,1)))         &
           +grav*(mc*fc*1.61-1.)*(tempt/tempv)*(qthl(iq,1)-qtot(iq,1))/(zzh(iq,1)-zz(iq,1))
        ppb(iq,1) = -km(iq,1)*bvf
      end do
      
    case(2) ! dry convection from Hurley 2007
      call updatekmo(thetavhl,thetav,fzzh,imax,kl)
      do k=2,kl-1
        do iq = 1,imax
          tcc=-grav*km(iq,k)*(thetavhl(iq,k)-thetavhl(iq,k-1))/(thetav(iq,k)*dz_fl(iq,k))
          ppb(iq,k)=tcc
        end do
      end do
      do iq = 1,imax
        tcc=-grav*km(iq,1)*(thetavhl(iq,1)-thetav(iq,1))/(thetav(iq,1)*(zzh(iq,1)-zz(iq,1)))
        ppb(iq,1)=tcc 
      end do

  end select

  ! Calculate transport source term on full levels
  ppt(:,2:kl-1)= kmo(:,2:kl-1)*idzp(:,2:kl-1)*(tke(:,3:kl)-tke(:,2:kl-1))/dz_hl(:,2:kl-1)  &
               -kmo(:,1:kl-2)*idzm(:,2:kl-1)*(tke(:,2:kl-1)-tke(:,1:kl-2))/dz_hl(:,1:kl-2)
  
  if ( upshear==1 ) then
    call updatekmo(uo_hl,uo,fzzh,imax,kl)
    call updatekmo(vo_hl,vo,fzzh,imax,kl)
    ! pps(:,1) is not used
    pps(:,1) = kmo(:,1)*(((uo_hl(:,1)-uo(:,1))/(zzh(:,1)-zz(:,1)))**2 + &
                         ((vo_hl(:,1)-vo(:,1))/(zzh(:,1)-zz(:,1)))**2)
    do k = 2,kl-1  
      pps(:,k) = kmo(:,k)*(((uo_hl(:,k)-uo_hl(:,k-1))/dz_fl(:,k))**2 + &
                           ((vo_hl(:,k)-vo_hl(:,k-1))/dz_fl(:,k))**2)
    end do  
  end if
  
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
    
  ! estimate eddy diffusivity mixing coeff
  km = cm0*tke(:,:)**2/eps(:,:)
  call updatekmo(kmo,km,fzzh,imax,kl) ! interpolate diffusion coeffs to half levels
  
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
  bb(:,2:kl-1)=1.-qq(:,2:kl-1)-rr(:,2:kl-1)+ddts*(mflx(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)         &
                                                 -mflx(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1))
  aa(:,kl)=qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
  bb(:,kl)=1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)


  ! theta vertical mixing
  dd(:,1)=thetal(:,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                               &
                           +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                   &
                     +ddts*rhos*wt0/(rhoa(:,1)*dz_fl(:,1))
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
                           +ddts*rhos*wq0/(rhoa(:,1)*dz_fl(:,1))
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
  bb(:,1) = 1. - cc(:,1) + ddts*rhos*cduv/(rhoa(:,1)*dz_fl(:,1)) ! implicit  
  bb(:,2:kl-1) = 1. - aa(:,2:kl-1) - cc(:,2:kl-1)
  bb(:,kl) = 1. - aa(:,kl)
  dd(:,1:kl) = uo(:,1:kl)
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
  dd(:,1:kl) = vo(:,1:kl)
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax,kl)
  
  
  ! update surface momentum flux
  ustar = sqrt(cduv*sqrt(uo(:,1)**2+vo(:,1)**2))
  ustar_ave = ustar_ave + ustar/real(mcount)

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
    where ( qfg(:,k)>1.e-12 )
      fice = min( qfg(:,k)/(qfg(:,k)+qlg(:,k)), 1. )
    elsewhere
      fice = 0.
    end where
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
  
#ifdef scm
  uwflux(:,1) = 0.
  uwflux(:,2:kl) = -kmo(:,1:kl-1)*(uo(:,2:kl)-uo(:,1:kl-1))/dz_hl(:,1:kl-1)
  vwflux(:,1) = 0.
  vwflux(:,2:kl) = -kmo(:,1:kl-1)*(vo(:,2:kl)-vo(:,1:kl-1))/dz_hl(:,1:kl-1)
  do k = 1,kl
    wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                                 *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
  end do
#endif
  
end do

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
                     zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,         &
                     stratcloud,km,wt0,wq0,wtv0,ps,ustar,              &
                     sig,sigkap,tke,eps,uo,vo,imax,kl)
!$acc routine vector

integer, intent(in) :: imax, kl
integer k
real, dimension(imax,kl), intent(inout) :: mflx, tlup, qvup, qlup, qfup, cfup
real, dimension(imax,kl), intent(in) :: theta, qvg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: zz, thetal, thetav, km, uo, vo 
real, dimension(imax,kl), intent(in) :: dz_hl
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(in) :: eps
real, dimension(imax), intent(in) :: wt0, wq0, wtv0, ps, ustar
real, dimension(kl), intent(in) :: sig, sigkap
real, dimension(imax), intent(inout) :: zi, wstar
real, dimension(imax,kl) ::  w2up, nn, cxup, rino
real, dimension(imax) :: dzht, ent, templ, pres, upf, qxup, qupsat
real, dimension(imax) :: fice, lx, qcup, dqsdt, al, xp, as, bs, cs
real, dimension(imax) :: thup, tvup, qtup, vvk
real, parameter :: fac = 10. ! originally fac=100.
real, parameter :: ricr = 0.3
logical, dimension(imax), intent(in) :: mask

! Initialise updraft
do k = 1,kl
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
    nn(:,k) = grav*(tvup-thetav(:,K))/thetav(:,k)                    ! calculate buoyancy
    w2up(:,k) = (w2up(:,k-1)+2.*dzht*b2*nn(:,k))/(1.+2.*dzht*b1*ent) ! update updraft velocity
  elsewhere
    nn(:,k) = 0.  
    w2up(:,k) = 0.
  end where
  vvk(:) = (uo(:,k)-uo(:,1))**2 + (vo(:,k)-vo(:,1))**2 + fac*ustar(:)**2  
  rino(:,k) = grav*(thetav(:,k)-thetav(:,1))*(zz(:,k)-zz(:,1))/max(thetav(:,1)*vvk,1.e-30)
  ! test if maximum plume height is reached
  where ( w2up(:,k)<=0. .and. w2up(:,k-1)>0. .and. mask )
    as = 2.*b2*(nn(:,k)-nn(:,k-1))/dzht
    bs = 2.*b2*nn(:,k-1)
    cs = w2up(:,k-1)
    xp = -2.*cs/(bs-sqrt(max(bs**2-4.*as*cs,0.)))
    xp = min(max(xp,0.),dzht)
    zi(:) = xp + zz(:,k-1)
  elsewhere ( rino(:,k)>ricr .and. rino(:,k-1)<=ricr .and. .not.mask )
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
! Tri-diagonal solver (array version)

pure subroutine thomas(outdat,aai,bbi,cci,ddi,imax,klin)
!$acc routine vector

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
!$acc routine vector

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
!$acc routine vector

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
!$acc routine vector

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
!$acc routine vector

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
! End TKE-eps

subroutine tkeend

implicit none

deallocate(tke,eps)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps

