! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public tkeinit,tkemix,tkeend,tke,eps,shear,zidry
public cm0,ce0,ce1,ce2,ce3,cq,be,ent0,ezmin,dtrn0,dtrc0,m0,b1,b2
public buoymeth,icm1,maxdts,mintke,mineps,minl,maxl,zidrytol,stabmeth
public tke_umin
#ifdef offline
public wthl,wqv,wql,wqf
public mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
public ents,dtrs
#endif

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:), allocatable, save :: zidry
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
real, save :: ezmin   = 100.   ! Limits entrainment at plume top
real, save :: dtrn0   = 0.4    ! Unsaturated detrainment constant
real, save :: dtrc0   = 0.9    ! Saturated detrainment constant
real, save :: m0      = 0.1    ! Mass flux amplitude constant
real, save :: b1      = 2.     ! Updraft entrainment coeff (Soares et al (2004) 1., Siebesma et al (2003) 2.)
real, save :: b2      = 1./3.  ! Updraft buoyancy coeff (Soares et al (2004) 2., Siebesma et al (2003) 1./3.)
! numerical constants
integer, save :: buoymeth = 0        ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0        ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: icm1     = 5        ! max iterations for calculating pblh
real, save :: maxdts      = 120.     ! max timestep for split
real, save :: mintke      = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps      = 1.E-10   ! min value for eps (1.0e-6 in TAPM)
real, save :: minl        = 1.       ! min value for L   (5. in TAPM)
real, save :: maxl        = 1000.    ! max value for L   (500. in TAPM)
real, save :: zidrytol    = 1.       ! tolerance for calculating pblh (m)
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
allocate(shear(ifull,kl),zidry(ifull))

cm34=cm0**0.75
tke=mintke
eps=mineps
shear=0.
zidry=0.

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

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,qrg,qsg,qgrg,cfrac,cfrain,cfsnow,cfgrap,uo,vo, &
                  zi,fg,eg,ps,zom,zz,zzh,sig,rhos,dt,qgmin,mode,diag,naero,aero,cgmap)

implicit none

integer, intent(in) :: diag,mode,naero
integer k,i,j,ktopmax,klcl
integer kcount,mcount,icount
real, intent(in) :: dt,qgmin
real, dimension(:,:,:), intent(inout) :: aero
real, dimension(:,:), intent(inout) :: theta,cfrac,cfrain,cfsnow,cfgrap,uo,vo
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg,qrg,qsg,qgrg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl), intent(in) :: zz,zzh
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: fg,eg,ps,zom,rhos,cgmap
real, dimension(kl), intent(in) :: sig
real, dimension(ifull,kl,naero) :: arup
real, dimension(ifull,kl) :: gamtl,gamqv,gamql,gamqf,gamqr,gamqs,gamqgr
real, dimension(ifull,kl) :: km,thetav,thetal,temp,qsat
real, dimension(ifull,kl) :: qsatc,qgnc,ff,thetac,tempc
real, dimension(ifull,kl) :: thetalhl,thetavhl
real, dimension(ifull,kl) :: quhl,qshl,qlhl,qfhl
real, dimension(ifull,kl) :: tkenew,epsnew,bb,cc,dd,rr
real, dimension(ifull,kl) :: rhoa,rhoahl
real, dimension(ifull,kl) :: pres,qtot,qthl
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
real, dimension(ifull) :: tff,tgg
real, dimension(ifull) :: cdrag,umag,ustar
real, dimension(ifull) :: tempv,rvar,bvf,dc,mc,fc
real, dimension(ifull) :: tbb,tcc,tqq
real, dimension(ifull) :: avearray
real, dimension(kl) :: sigkap,w2up,nn,dqdash,qupsat
real, dimension(kl) :: qtup,ttup,tvup,thup
real, dimension(1) :: templ
real qrup,qsup,qgup
real xp,as,bs,cs,cm12,cm34,qcup
real dzht,ziold,ent,dtr,dtrc,dtrn,dtrx
real ddts,zlcl
real lx,tempd,fice,qxup,dqsdt,al
real sigqtup,rng
logical, dimension(ifull,kl) :: lta
logical scond

cm12=1./sqrt(cm0)
cm34=sqrt(sqrt(cm0**3))
zlcl=0.
ktopmax=0

if (diag>0) write(6,*) "Update PBL mixing with TKE-eps + MF turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

do k=1,kl
  ! Impose limits on tke and eps after being advected by the host model
  tke(1:ifull,k)=max(tke(1:ifull,k),mintke)
  tff=cm34*tke(1:ifull,k)*sqrt(tke(1:ifull,k))
  eps(1:ifull,k)=min(eps(1:ifull,k),tff/minl)
  eps(1:ifull,k)=max(eps(1:ifull,k),tff/maxl,mineps)

  ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
  sigkap(k)=sig(k)**(-rd/cp)
  pres(:,k)=ps(:)*sig(k) ! pressure
  ! Density must be updated when dz is updated so that rho*dz is conserved
  thetav(:,k)=theta(1:ifull,k)*(1.+0.61*qvg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k)-qrg(1:ifull,k) &
                                -qsg(1:ifull,k)-qgrg(1:ifull,k))
  rhoa(:,k)=sigkap(k)*pres(:,k)/(rd*thetav(:,k))

  ! Transform to thetal as it is the conserved variable
  thetal(:,k)=theta(1:ifull,k)-sigkap(k)*(lv*(qlg(1:ifull,k)+qrg(1:ifull,k)) &
              +ls*(qfg(1:ifull,k)+qsg(1:ifull,k)+qgrg(1:ifull,k)))/cp
  
  ! Calculate first approximation to diffusion coeffs
  km(:,k)=cm0*tke(1:ifull,k)*tke(1:ifull,k)/eps(1:ifull,k)
end do

! Calculate surface fluxes
wt0=fg/(rhos*cp)  ! theta flux
wq0=eg/(rhos*lv)  ! qtot flux (=qv flux)

do k=1,kl-1
  ! Fraction for interpolation
  fzzh(:,k)=(zzh(:,k)-zz(:,k))/(zz(:,k+1)-zz(:,k))

  ! Calculate dz at half levels
  dz_hl(:,k)=zz(:,k+1)-zz(:,k)
end do

! Calculate dz at full levels
dz_fl(:,1)   =zzh(:,1)
dz_fl(:,2:kl)=zzh(:,2:kl)-zzh(:,1:kl-1)

! Calculate shear term on full levels (see hordifg.f for calculation of horizontal shear)
pps(:,2:kl-1)=km(:,2:kl-1)*shear(:,2:kl-1)

! set top BC for TKE-eps source terms
pps(:,kl)=0.
ppb(:,kl)=0.
ppt(:,kl)=0.

! interpolate diffusion coeff and air density to half levels
call updatekmo(kmo,   km,  fzzh)
call updatekmo(rhoahl,rhoa,fzzh)
idzm(:,2:kl)  =rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1)=rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

! Main loop to prevent time splitting errors
mcount=int(dt/(maxdts+0.01))+1
ddts  =dt/real(mcount)
do kcount=1,mcount

  ! Update momentum flux
  wtv0 = wt0 + theta(1:ifull,1)*0.61*wq0 ! thetav flux
  umag = sqrt(max( uo(1:ifull,1)*uo(1:ifull,1)+vo(1:ifull,1)*vo(1:ifull,1), tke_umin ))
  call dyerhicks(cdrag,wtv0,zom,umag,thetav(:,1),zz(:,1))
  ustar = sqrt(cdrag)*umag
    
  ! Set-up thermodynamic variables temp, theta_v and surface fluxes
  do k=1,kl
    temp(:,k)=theta(1:ifull,k)/sigkap(k)
    ! calculate saturated air mixing ratio
    call getqsat(qsat(:,k),temp(:,k),pres(:,k))
  end do
  thetav=theta(1:ifull,:)*(1.+0.61*qvg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:)-qrg(1:ifull,:) &
         -qsg(1:ifull,:)-qgrg(1:ifull,:))
  qtot=qvg(1:ifull,:)+qlg(1:ifull,:)+qfg(1:ifull,:)+qrg(1:ifull,:)+qsg(1:ifull,:)+qgrg(1:ifull,:)

  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  mflx=0.
  tlup=thetal(1:ifull,:)
  qvup=qvg(1:ifull,:)
  qlup=qlg(1:ifull,:)
  qfup=qfg(1:ifull,:)
  cfup=cfrac(1:ifull,:)
  if (naero>0) then
    arup=aero(1:ifull,:,:)
  end if

#ifdef offline
  mf=0.
  w_up=0.
  tl_up=thetal(1:ifull,:)
  qv_up=qvg(1:ifull,:)
  ql_up=qlg(1:ifull,:)
  qf_up=qfg(1:ifull,:)
  cf_up=0.
  ents=0.
  dtrs=0.
#endif

  if (mode/=1) then ! mass flux

    ! Note that wstar is actually based on zidry, not zi
    wstar=(grav*zidry*max(wtv0,0.)/thetav(:,1))**(1./3.)
  
    do i=1,ifull
      if (wtv0(i)>0.) then ! unstable
        do icount=1,icm1
          ! Initialise updraft
          ziold=zidry(i)
          zlcl=0.
          tke(i,1)=cm12*ustar(i)*ustar(i)+ce3*wstar(i)*wstar(i)
          tke(i,1)=max(tke(i,1),mintke)
          ktopmax=0
          klcl=kl+1
          w2up=0.
          cfup(i,:)=0.
          scond=.false.
          dzht=zz(i,1)
          ! Entrainment and detrainment rates
          ent=entfn(zz(i,1),zi(i))
          
          ! first level -----------------
          ! initial thermodynamic state
          ! split qtot into components (conservation of thetal and qtot is maintained)
          tlup(i,1)=thetal(i,1)+be*wt0(i)/sqrt(tke(i,1))       ! Hurley 2007
          qvup(i,1)=qvg(i,1)   +be*wq0(i)/sqrt(tke(i,1))       ! Hurley 2007
          qlup(i,1)=qlg(i,1)
          qfup(i,1)=qfg(i,1)
          qrup=qrg(i,1)
          qsup=qsg(i,1)
          qgup=qgrg(i,1)
          ! diagnose thermodynamic variables assuming no condensation
          qtup(1)=qvup(i,1)+qlup(i,1)+qfup(i,1)+qrup+qsup+qgup     ! qtot,up
          ! state of plume after evaporation
          qvup(i,1)=qtup(1)
          qlup(i,1)=0.
          qfup(i,1)=0.
          qrup=0.
          qsup=0.
          qgup=0.
          thup(1)=tlup(i,1) !+sigkap(1)*(lv*qlup(1)+ls*qfup(1))/cp ! theta,up
          tvup(1)=thup(1)+theta(i,1)*0.61*qtup(1)                  ! thetav,up
          ttup(1)=thup(1)/sigkap(1)                                ! temp,up
          templ(1)=tlup(i,1)/sigkap(1)                             ! templ,up
          ! update updraft velocity and mass flux
          nn(1)  =grav*be*wtv0(i)/(thetav(i,1)*sqrt(tke(i,1))) ! Hurley 2007
          w2up(1)=2.*dzht*b2*nn(1)/(1.+2.*dzht*b1*ent)         ! Hurley 2007
          ! estimate variance of qtup in updraft
          call getqsat(qupsat(1:1),templ(1:1),pres(i:i,1))
          sigqtup=1.E-5
          rng=sqrt(6.)*sigqtup               ! variance of triangle distribution
          dqdash(1)=(qtup(1)-qupsat(1))/rng  ! scaled variance
          dqdash(1)=min(dqdash(1),-1.)
        
          ! updraft without condensation
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
            qrup     =(            dzht*ent*qrg(i,k)   )/(1.+dzht*ent)
            qsup     =(            dzht*ent*qsg(i,k)   )/(1.+dzht*ent)
            qgup     =(            dzht*ent*qgrg(i,k)  )/(1.+dzht*ent)
            ! calculate conserved variables
            qtup(k)=qvup(i,k)+qlup(i,k)+qfup(i,k)+qrup+qsup+qgup ! qtot,up
            ! estimate air temperature
            templ(1)=tlup(i,k)/sigkap(k)                         ! templ,up
            if (.not.scond) then
              ! MJT notes - here we determine a fraction of the updraft that has undergone condensation
              ! More correctly, we estimate a fraction of the number of updrafts at different locations
              ! in the grid-box that have undergone condesation.
              call getqsat(qupsat(k:k),templ(1:1),pres(i:i,k))
              ! estimate variance of qtup in updraft (following Hurley and TAPM)
              sigqtup=sqrt(max(1.E-10,1.6*tke(i,k)/eps(i,k)*cq*km(i,k)*((qtup(k)-qtup(k-1))/dzht)**2))
              ! MJT condensation scheme - follow Smith 1990 and assume triangle distribution for qtup.
              rng=sqrt(6.)*sigqtup               ! variance of triangle distribution
              dqdash(k)=(qtup(k)-qupsat(k))/rng  ! scaled variance
              if (dqdash(k)>=-1.) then
                scond=.true.
                klcl=k
                xp=dzht*(-1.-dqdash(k-1))/(dqdash(k)-dqdash(k-1))
                xp=min(max(xp,0.),dzht)
                zlcl=xp+zz(i,k-1)
              end if
            end if
            ! state of plume after redistribution
            qvup(i,k)=qtup(k)
            qlup(i,k)=0.
            qfup(i,k)=0.
            qrup=0.
            qsup=0.
            qgup=0.
            thup(k)=tlup(i,k) !+sigkap(k)*(lv*qlup(k)+ls*qfup(k))/cp ! theta,up
            tvup(k)=tlup(i,k)+theta(i,k)*0.61*qtup(k)                ! thetav,up after redistribution
            ttup(k)=thup(k)/sigkap(k)                                ! temp,up
            ! calculate buoyancy
            nn(k)  =grav*(tvup(k)-thetav(i,k))/thetav(i,k)
            ! update updraft velocity
            w2up(k)=(w2up(k-1)+2.*dzht*b2*nn(k))/(1.+2.*dzht*b1*ent)
            ! test if maximum plume height is reached
            if (w2up(k)<=0.) then
              as=min(2.*b2*(nn(k)-nn(k-1))/dzht,-1.E-20)
              bs=2.*b2*nn(k-1)
              cs=w2up(k-1)
              xp=0.5*(-bs-sqrt(max(bs*bs-4.*as*cs,0.)))/as
              xp=min(max(xp,0.),dzht)
              zidry(i)=xp+zz(i,k-1)
              ktopmax=max(ktopmax,k-1)
              exit
            end if
          end do
          
          if (.not.scond) then
            zi(i)=zidry(i)
            zlcl=zidry(i)
          else
            zi(i)=max(zi(i),zidry(i))
            zlcl =min(zlcl, zidry(i))

            ! updraft with condensation
            do k=klcl,kl
              dzht=dz_hl(i,k-1)
              ! Entrainment and detrainment rates
              ent=entfn(zz(i,k),zi(i))
              ! update thermodynamics of plume
              ! split qtot into components (conservation is maintained)
              tlup(i,k)=(tlup(i,k-1)+dzht*ent*thetal(i,k))/(1.+dzht*ent)
              qvup(i,k)=(qvup(i,k-1)+dzht*ent*qvg(i,k)   )/(1.+dzht*ent)
              qlup(i,k)=(qlup(i,k-1)+dzht*ent*qlg(i,k)   )/(1.+dzht*ent)
              qfup(i,k)=(qfup(i,k-1)+dzht*ent*qfg(i,k)   )/(1.+dzht*ent)
              qrup     =(            dzht*ent*qrg(i,k)   )/(1.+dzht*ent)
              qsup     =(            dzht*ent*qsg(i,k)   )/(1.+dzht*ent)
              qgup     =(            dzht*ent*qgrg(i,k)  )/(1.+dzht*ent)
              ! calculate conserved variables
              qtup(k)=qvup(i,k)+qlup(i,k)+qfup(i,k)+qrup+qsup+qgup ! qtot,up
              ! estimate air temperature
              templ(1)=tlup(i,k)/sigkap(k)                         ! templ,up
              call getqsat(qupsat(k:k),templ(1:1),pres(i:i,k))
              ! estimate variance of qtup in updraft (following Hurley and TAPM)
              sigqtup=sqrt(max(1.E-10,1.6*tke(i,k)/eps(i,k)*cq*km(i,k)*((qtup(k)-qtup(k-1))/dzht)**2))
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
              thup(k)=tlup(i,k)+sigkap(k)*(lv*(qlup(i,k)+qrup)+ls*(qfup(i,k)+qsup+qgup))/cp ! theta,up before redistribution
              tempd  =thup(k)/sigkap(k)                                                     ! temp,up before redistribution
              fice=min(max(273.16-tempd,0.),40.)/40. ! approximate ice fraction based on temperature (not templ)
              lx=lv+lf*fice
              dqsdt=qupsat(k)*lx/(rv*templ(1)*templ(1))
              al=cp/(cp+lx*dqsdt)
              qcup=al*(qtup(k)-qxup)                                 ! qcondensate,up after redistribution
              qxup=qtup(k)-qcup                                      ! qv,up after redistribution
              ttup(k)=templ(1)+lx*qcup/cp                            ! temp,up after redistribution
              thup(k)=ttup(k)*sigkap(k)                              ! theta,up after redistribution
              tvup(k)=thup(k)+theta(i,k)*(1.61*qxup-qtup(k))         ! thetav,up after redistribution
              ! state of plume after redistribution
              qvup(i,k)=qxup                                         ! qv,up after redistribution
              qlup(i,k)=qcup*(1.-fice)                               ! ql,up after redistribution
              qfup(i,k)=qcup*fice                                    ! qf,up after redistribution
              qrup=0.                                                ! qr,up after redistribution
              qsup=0.                                                ! qs,up after redistribution
              qgup=0.                                                ! qg,up after redistribution
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
          
          end if

          ! update surface boundary conditions (note zidry, rather than zi)
          wstar(i)=(grav*zidry(i)*max(wtv0(i),0.)/thetav(i,1))**(1./3.)
          
          zi(i)=max(zidry(i),zi(i))

          ! check for convergence
          if (abs(zidry(i)-ziold)<zidrytol) exit
        end do

        ! update mass flux
        mflx(i,1)=m0*sqrt(max(w2up(1),0.))*zz(i,1)**ent0*max(zi(i)-zz(i,1),1.)**dtrn0 ! MJT suggestion
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          xp  =(zz(i,k)-zlcl)/max(zidry(i)-zlcl,0.1)
          xp  =min(max(xp,0.),1.)
          ent =entfn(zz(i,k),zi(i))
          dtrn=dtrfn(zz(i,k),zidry(i),dtrn0)
          dtrc=dtrfn(zz(i,k),zi(i),   dtrn0)
          dtrx=(1.-xp)*dtrn+xp*dtrc
          dtrc=dtrfn(zz(i,k),zi(i),   dtrc0)
          dtr =(1.-cfup(i,k))*dtrx+cfup(i,k)*dtrc
          mflx(i,k)=mflx(i,k-1)/(1.+dzht*(dtr-ent))
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
          
        dzht=zz(i,1)
        ! Entrainment and detrainment rates
        xp  =(zz(i,1)-zlcl)/max(zidry(i)-zlcl,0.1)
        xp  =min(max(xp,0.),1.)
        ent =entfn(zz(i,1),zi(i))
        dtrn=dtrfn(zz(i,1),zidry(i),dtrn0)
        dtrc=dtrfn(zz(i,1),zi(i),   dtrn0)
        dtrx=(1.-xp)*dtrn+xp*dtrc
        dtrc=dtrfn(zz(i,1),zi(i),   dtrc0)
        dtr =(1.-cfup(i,1))*dtrx+cfup(i,1)*dtrc
        ents(i,1)=ent
        dtrs(i,1)=dtr
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          xp  =(zz(i,k)-zlcl)/max(zidry(i)-zlcl,0.1)
          xp  =min(max(xp,0.),1.)
          ent =entfn(zz(i,k),zi(i))
          dtrn=dtrfn(zz(i,k),zidry(i),dtrn0)
          dtrc=dtrfn(zz(i,k),zi(i),   dtrn0)
          dtrx=(1.-xp)*dtrn+xp*dtrc          
          dtrc=dtrfn(zz(i,k),zi(i),   dtrc0)
          dtr =(1.-cfup(i,k))*dtrx+cfup(i,k)*dtrc
          ents(i,k)=ent
          dtrs(i,k)=dtr
        end do
#endif

        ! update reamining scalars which are not used in the iterative loop
        do j=1,naero
          arup(i,1,j)=aero(i,1,j)
          do k=2,ktopmax
            dzht=dz_hl(i,k-1)
            ent=entfn(zz(i,k),zi(i))
            arup(i,k,j)=(arup(i,k-1,j)+dzht*ent*aero(i,k,j))/(1.+dzht*ent)
          end do
        end do
        
      else                   ! stable
        !wpv_flux is calculated at half levels
        !wpv_flux(1)=-kmo(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1) !+gamt_hl(i,k)
        !do k=2,kl-1
        !  wpv_flux(k)=-kmo(i,k)*(thetav(i,k+1)-thetav(i,k))/dz_hl(i,k) !+gamt_hl(i,k)
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
        zi(i)=zz(i,1) ! MJT suggestion
        zidry(i)=zz(i,1)
      end if
    end do
  
  else
    zidry=zi   
    ! Note that wstar is actually based on zidry, not zi
    wstar=(grav*zidry*max(wtv0,0.)/thetav(:,1))**(1./3.)   
  end if

  
  ! turn off MF term if small grid spacing
  do k=1,kl
    mflx(:,k)=mflx(:,k)*cgmap(:)
  end do
  
  
  ! explicit version of counter gradient terms
  gamtl(:,:) =mflx(:,:)*(tlup(:,:)-thetal(1:ifull,:))
  gamqv(:,:) =mflx(:,:)*(qvup(:,:)-qvg(1:ifull,:))
  gamql(:,:) =mflx(:,:)*(qlup(:,:)-qlg(1:ifull,:))
  gamqf(:,:) =mflx(:,:)*(qfup(:,:)-qfg(1:ifull,:))
  gamqr(:,:) =mflx(:,:)*(         -qrg(1:ifull,:))
  gamqs(:,:) =mflx(:,:)*(         -qsg(1:ifull,:))
  gamqgr(:,:)=mflx(:,:)*(         -qgrg(1:ifull,:))
  
  
  ! calculate tke and eps at 1st level
  z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-10))
  z_on_l=min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  select case(stabmeth)
    case(0)
      where (z_on_l<0.)
        phim=(1.-16.*z_on_l)**(-0.25)
      elsewhere
        phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
      end where
    case(1)
      where (z_on_l<0.)
        phim=(1.-16.*z_on_l)**(-0.25)
      elsewhere (z_on_l<=0.4)
        phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
      elsewhere
        phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar
      end where
    case default
      write(6,*) "ERROR: Invalid option for stabmeth in tkeeps ",stabmeth
      stop
  end select
  tke(1:ifull,1)=cm12*ustar*ustar+ce3*wstar*wstar
  eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
  tke(1:ifull,1)=max(tke(1:ifull,1),mintke)
  tff=cm34*tke(1:ifull,1)*sqrt(tke(1:ifull,1))
  eps(1:ifull,1)=min(eps(1:ifull,1),tff/minl)
  eps(1:ifull,1)=max(eps(1:ifull,1),tff/maxl,mineps)


  ! Update TKE and eps terms

  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tke(1:ifull,kl)=mintke
  eps(1:ifull,kl)=mineps
  
  ! Calculate buoyancy term
  select case(buoymeth)
    case(0) ! saturated from Durran and Klemp JAS 1982 (see also WRF)
      qsatc=max(qsat,qvg(1:ifull,:))                                             ! assume qvg is saturated inside cloud
      ff=qfg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                            &
        +qsg(1:ifull,:)/max(cfrac(1:ifull,:),cfsnow(1:ifull,:),1.E-8)          &
        +qgrg(1:ifull,:)/max(cfrac(1:ifull,:),cfgrap(1:ifull,:),1.E-8)           ! inside cloud value  assuming max overlap
      dd=qlg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                            &
        +qrg(1:ifull,:)/max(cfrac(1:ifull,:),cfrain(1:ifull,:),1.E-8)            ! inside cloud value assuming max overlap
      do k=1,kl
        tbb=max(1.-cfrac(1:ifull,k),1.E-8)
        qgnc(:,k)=(qvg(1:ifull,k)-(1.-tbb)*qsatc(:,k))/tbb                       ! outside cloud value
        qgnc(:,k)=min(max(qgnc(:,k),qgmin),qsatc(:,k))
        thetac(:,k)=thetal(:,k)+sigkap(k)*(lv*dd(:,k)+ls*ff(:,k))/cp             ! inside cloud value
        tempc(:,k)=thetac(:,k)/sigkap(k)                                         ! inside cloud value
      end do
      call updatekmo(thetalhl,thetal,fzzh)                                       ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh)                                             ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh)                                            ! inside cloud value
      call updatekmo(qlhl,dd,fzzh)                                               ! inside cloud value
      call updatekmo(qfhl,ff,fzzh)                                               ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl)=cfrac(1:ifull,2:kl)<=1.E-6
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
        tqq=(1.+lv*qsatc(:,k)/(rd*tempc(:,k)))/(1.+lv*lv*qsatc(:,k)/(cp*rv*tempc(:,k)*tempc(:,k)))
        tbb=-grav*km(:,k)*(tqq*((thetalhl(:,k)-thetalhl(:,k-1)+sigkap(k)/cp*(lv*(qlhl(:,k)-qlhl(:,k-1))  &
            +ls*(qfhl(:,k)-qfhl(:,k-1))))/thetac(:,k)+lv/cp*(qshl(:,k)-qshl(:,k-1))/tempc(:,k))          &
            -qshl(:,k)-qlhl(:,k)-qfhl(:,k)+qshl(:,k-1)+qlhl(:,k-1)+qfhl(:,k-1))/dz_fl(:,k)
        tbb=tbb+grav*(tqq*((gamtl(:,k)                                                                   &
            +sigkap(k)/cp*(lv*gamql(:,k)/max(cfrac(1:ifull,k),1.E-8)                                     &
                          +lv*gamqr(:,k)/max(cfrac(1:ifull,k),cfrain(1:ifull,k),1.E-8)                   &
                          +ls*gamqf(:,k)/max(cfrac(1:ifull,k),1.E-8)                                     &
                          +ls*gamqs(:,k)/max(cfrac(1:ifull,k),cfsnow(1:ifull,k),1.E-8)                   &
                          +ls*gamqgr(:,k)/max(cfrac(1:ifull,k),cfgrap(1:ifull,k),1.E-8)))/thetac(:,k)    &
            +lv/cp*gamqv(:,k)/tempc(:,k))                                                                &
            -gamqv(:,k)-(gamql(:,k)+gamqf(:,k))/max(cfrac(1:ifull,k),1.E-8)                              &
            -gamqr(:,k)/max(cfrac(1:ifull,k),cfrain(1:ifull,k),1.E-8)                                    &
            -gamqs(:,k)/max(cfrac(1:ifull,k),cfsnow(1:ifull,k),1.E-8)                                    &
            -gamqgr(:,k)/max(cfrac(1:ifull,k),cfrain(1:ifull,k),1.E-8))
        ! unsaturated
        tcc=-grav*km(:,k)*(thetalhl(:,k)-thetalhl(:,k-1)+thetal(1:ifull,k)*0.61*(quhl(:,k)-quhl(:,k-1))) &
                         /(thetal(1:ifull,k)*dz_fl(:,k))
        tcc=tcc+grav*(gamtl(:,k)+thetal(1:ifull,k)*0.61*gamqv(:,k))/thetal(1:ifull,k)
        ppb(:,k)=(1.-cfrac(1:ifull,k))*tcc+cfrac(1:ifull,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)
      end do
      
    case(1) ! follow Marquet and Geleyn QJRMS (2012)
      call updatekmo(thetalhl,thetal,fzzh)
      call updatekmo(qthl,qtot,fzzh)
      do k=2,kl-1
        tempv=temp(1:ifull,k)*thetav(1:ifull,k)/theta(1:ifull,k)
        rvar=rd*tempv/temp(1:ifull,k) ! rvar = qd*rd+qv*rv
        fc=(1.-cfrac(1:ifull,k))+cfrac(1:ifull,k)*(lv*rvar/(cp*rv*temp(1:ifull,k)))
        dc=(1.+0.61*qvg(1:ifull,k))*lv*qvg(1:ifull,k)/(rd*tempv)
        mc=(1.+dc)/(1.+(lv*(qlg(1:ifull,k)+qrg(1:ifull,k))+ls*(qfg(1:ifull,k)                &
           +qsg(1:ifull,k)+qgrg(1:ifull,k)))/(cp*temp(1:ifull,k))+dc*fc)
        bvf=grav*mc*(thetalhl(:,k)-thetalhl(:,k-1))/(thetal(1:ifull,k)*dz_fl(:,k))           &
           +grav*(mc*fc*1.61-1.)*(temp(1:ifull,k)/tempv)*(qthl(:,k)-qthl(:,k-1))/dz_fl(:,k)
        tcc=-grav*mc*gamtl(:,k)/thetal(1:ifull,k)                                            &
            -grav*(mc*fc*1.61-1.)*(temp(1:ifull,k)/tempv)                                    &
                 *(gamqv(:,k)+gamql(:,k)+gamqf(:,k)+gamqr(:,k)+gamqs(:,k)+gamqgr(:,k))
        ppb(:,k)=-km(:,k)*bvf-tcc
      end do
      
    case(2) ! dry convection
      call updatekmo(thetavhl,thetav,fzzh)
      do k=2,kl-1
        tcc=-grav*km(:,k)*(thetavhl(:,k)-thetavhl(:,k-1))/(thetav(:,k)*dz_fl(:,k))
        tcc=tcc+grav*(gamtl(:,k)+sigkap(k)/cp*(lv*(gamql(:,k)+gamqr(:,k))                &
            +ls*(gamqf(:,k)+gamqs(:,k)+gamqgr(:,k))))/thetav(:,k)
        tcc=tcc+grav*(theta(1:ifull,k)*(0.61*gamqv(:,k)-gamql(:,k)-gamqf(:,k)-gamqr(:,k) &
            -gamqs(:,k)-gamqgr(:,k)))/thetav(:,k)
        ppb(:,k)=tcc
      end do
      
   case default
     write(6,*) "ERROR: Unknown buoymeth option ",buoymeth
     stop
  end select

  ! Calculate transport source term on full levels
  do k=2,kl-1
    ppt(:,k)= kmo(:,k)*idzp(:,k)*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k)                &
             -kmo(:,k-1)*idzm(:,k)*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1)
  end do
  
  qq(:,2:kl-1)=-ddts*idzm(:,2:kl-1)/dz_hl(:,1:kl-2)
  rr(:,2:kl-1)=-ddts*idzp(:,2:kl-1)/dz_hl(:,2:kl-1)
  
  ! eps vertical mixing
  aa(:,2:kl-1)=ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=ce0*kmo(:,2:kl-1)*rr(:,2:kl-1)
  ! follow PH to make scheme more numerically stable
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)+ddts*ce2*eps(1:ifull,2:kl-1)/tke(1:ifull,2:kl-1)
  dd(:,2:kl-1)=eps(1:ifull,2:kl-1)+ddts*eps(1:ifull,2:kl-1)/tke(1:ifull,2:kl-1)                    &
              *ce1*(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
  dd(:,2)     =dd(:,2)   -aa(:,2)*eps(1:ifull,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mineps
  call thomas(epsnew(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1))

  ! TKE vertical mixing
  aa(:,2:kl-1)=kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=kmo(:,2:kl-1)*rr(:,2:kl-1)
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  dd(:,2:kl-1)=tke(1:ifull,2:kl-1)+ddts*(pps(:,2:kl-1)+ppb(:,2:kl-1)-epsnew(:,2:kl-1))
  dd(:,2)     =dd(:,2)   -aa(:,2)*tke(1:ifull,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mintke
  call thomas(tkenew(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1))

  do k=2,kl-1
    tke(1:ifull,k)=max(tkenew(:,k),mintke)
    tff=cm34*tke(1:ifull,k)*sqrt(tke(1:ifull,k))
    eps(1:ifull,k)=min(epsnew(:,k),tff/minl)
    eps(1:ifull,k)=max(eps(1:ifull,k),tff/maxl,mineps)
  end do
    
  km=cm0*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:)
  call updatekmo(kmo,km,fzzh) ! interpolate diffusion coeffs to half levels
  
  ! update scalars
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
  
  avearray=0.5*(maxval(thetal(1:ifull,:),dim=2)+minval(thetal(1:ifull,:),dim=2))
  do k=1,kl
    thetal(1:ifull,k)=thetal(1:ifull,k)-avearray
    tlup(:,k)=tlup(:,k)-avearray
  end do
  dd(:,1)=thetal(1:ifull,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                   &
                                 +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                       &
                           +ddts*rhos*wt0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=thetal(1:ifull,2:kl-1)+ddts*(mflx(:,1:kl-2)*tlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)     &
                                           +mflx(:,2:kl-1)*tlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)          &
                                           -mflx(:,2:kl-1)*tlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)     &
                                           -mflx(:,3:kl)*tlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=thetal(1:ifull,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                       &
                                   +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  do k=1,kl
    thetal(1:ifull,k)=thetal(1:ifull,k)+avearray
    tlup(:,k)=tlup(:,k)+avearray
  end do
#ifdef offline
  wthl(:,1:kl-1)=-kmo(:,1:kl-1)*(thetal(1:ifull,2:kl)-thetal(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)                    &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                           &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif

  avearray=0.5*(maxval(qvg(1:ifull,:),dim=2)+minval(qvg(1:ifull,:),dim=2))
  do k=1,kl
    qvg(1:ifull,k)=qvg(1:ifull,k)-avearray
    qvup(:,k)=qvup(:,k)-avearray
  end do
  dd(:,1)=qvg(1:ifull,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                      &
                              +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                                          &
                           +ddts*rhos*wq0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=qvg(1:ifull,2:kl-1)+ddts*(mflx(:,1:kl-2)*qvup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)        &
                                        +mflx(:,2:kl-1)*qvup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qvup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qvup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qvg(1:ifull,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                          &
                                +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  do k=1,kl
    qvg(1:ifull,k)=qvg(1:ifull,k)+avearray
    qvup(:,k)=qvup(:,k)+avearray
  end do  
#ifdef offline
  wqv(:,1:kl-1)=-kmo(:,1:kl-1)*(qvg(1:ifull,2:kl)-qvg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)                           &
                 +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif

  dd(:,1)=qlg(1:ifull,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                      &
                              +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qlg(1:ifull,2:kl-1)+ddts*(mflx(:,1:kl-2)*qlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)        &
                                        +mflx(:,2:kl-1)*qlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qlg(1:ifull,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                          &
                                +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wql(:,1:kl-1)=-kmo(:,1:kl-1)*(qlg(1:ifull,2:kl)-qlg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)                           &
                 +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif

  dd(:,1)=qfg(1:ifull,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                      &
                              +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qfg(1:ifull,2:kl-1)+ddts*(mflx(:,1:kl-2)*qfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)        &
                                        +mflx(:,2:kl-1)*qfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qfg(1:ifull,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                          &
                                +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wqf(:,1:kl-1)=-kmo(:,1:kl-1)*(qfg(1:ifull,2:kl)-qfg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)                           &
                 +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif

  dd(:,1:kl)=qrg(1:ifull,1:kl)
  call thomas(qrg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))

  dd(:,1:kl)=qsg(1:ifull,1:kl)
  call thomas(qsg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))

  dd(:,1:kl)=qgrg(1:ifull,1:kl)
  call thomas(qgrg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))  
  
  ! account for phase transitions
  do k=1,kl
    tgg=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k)+qrg(1:ifull,k) &
        +qsg(1:ifull,k)+qgrg(1:ifull,k),qgmin)                            ! qtot before phase transition
    qvg(1:ifull,k)=max(qvg(1:ifull,k),0.)    
    qlg(1:ifull,k)=max(qlg(1:ifull,k),0.)
    qfg(1:ifull,k)=max(qfg(1:ifull,k),0.)
    qrg(1:ifull,k)=max(qrg(1:ifull,k),0.)
    qsg(1:ifull,k)=max(qsg(1:ifull,k),0.)
    qgrg(1:ifull,k)=max(qgrg(1:ifull,k),0.)
    tff=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k)+qrg(1:ifull,k) &
        +qsg(1:ifull,k)+qgrg(1:ifull,k),qgmin)
    tgg=tgg/tff ! scale factor for conservation
    qvg(1:ifull,k)=qvg(1:ifull,k)*tgg
    qlg(1:ifull,k)=qlg(1:ifull,k)*tgg
    qfg(1:ifull,k)=qfg(1:ifull,k)*tgg
    qrg(1:ifull,k)=qrg(1:ifull,k)*tgg
    qsg(1:ifull,k)=qsg(1:ifull,k)*tgg
    qgrg(1:ifull,k)=qgrg(1:ifull,k)*tgg
    ! update theta for output or next time step
    theta(1:ifull,k)=thetal(:,k)+sigkap(k)*(lv*(qlg(1:ifull,k)+qrg(1:ifull,k)) &
                    +ls*(qfg(1:ifull,k)+qsg(1:ifull,k)+qgrg(1:ifull,k)))/cp
  end do

  ! update cloud fraction terms
  dd(:,1)=cfrac(1:ifull,1)-ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                    &
                                +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=cfrac(1:ifull,2:kl-1)+ddts*(mflx(:,1:kl-2)*cfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)      &
                                          +mflx(:,2:kl-1)*cfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)           &
                                          -mflx(:,2:kl-1)*cfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)      &
                                          -mflx(:,3:kl)*cfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=cfrac(1:ifull,kl)+ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                        &
                                  +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(cfrac,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfrac(1:ifull,:)=min(max(cfrac(1:ifull,:),0.),1.)
  where (qlg(1:ifull,:)+qfg(1:ifull,:)>1.E-12)
    cfrac(1:ifull,:)=max(cfrac(1:ifull,:),1.E-8)
  end where

  dd(:,1:kl)=cfrain(1:ifull,1:kl)
  call thomas(cfrain,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfrain(1:ifull,:)=min(max(cfrain(1:ifull,:),0.),1.)
  where (qrg(1:ifull,:)>1.E-12)
    cfrain(1:ifull,:)=max(cfrain(1:ifull,:),1.E-8)
  end where

  dd(:,1:kl)=cfsnow(1:ifull,1:kl)
  call thomas(cfsnow,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfsnow(1:ifull,:)=min(max(cfsnow(1:ifull,:),0.),1.)
  where (qsg(1:ifull,:)>1.E-12)
    cfsnow(1:ifull,:)=max(cfsnow(1:ifull,:),1.E-8)
  end where  

  dd(:,1:kl)=cfgrap(1:ifull,1:kl)
  call thomas(cfgrap,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfgrap(1:ifull,:)=min(max(cfgrap(1:ifull,:),0.),1.)
  where (qgrg(1:ifull,:)>1.E-12)
    cfgrap(1:ifull,:)=max(cfgrap(1:ifull,:),1.E-8)
  end where    
  
  ! Aerosols
  do j=1,naero
    dd(:,1)=aero(1:ifull,1,j)-ddts*(mflx(:,1)*arup(:,1,j)*(1.-fzzh(:,1))*idzp(:,1)                               &
                                   +mflx(:,2)*arup(:,2,j)*fzzh(:,1)*idzp(:,1))
    dd(:,2:kl-1)=aero(1:ifull,2:kl-1,j)+ddts*(mflx(:,1:kl-2)*arup(:,1:kl-2,j)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1) &
                                             +mflx(:,2:kl-1)*arup(:,2:kl-1,j)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)      &
                                             -mflx(:,2:kl-1)*arup(:,2:kl-1,j)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1) &
                                             -mflx(:,3:kl)*arup(:,3:kl,j)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
    dd(:,kl)=aero(1:ifull,kl,j)+ddts*(mflx(:,kl-1)*arup(:,kl-1,j)*(1.-fzzh(:,kl-1))*idzm(:,kl)                   &
                                     +mflx(:,kl)*arup(:,kl,j)*fzzh(:,kl-1)*idzm(:,kl))
    call thomas(aero(:,:,j),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  end do
  aero(:,:,:) = max( aero(:,:,:), 0. )

  ! Winds
  aa(:,2:kl)  =qq(:,2:kl)
  cc(:,1:kl-1)=rr(:,1:kl-1)
  bb(:,1)=1.-cc(:,1)+ddts*rhos*cdrag*umag/(rhoa(:,1)*dz_fl(:,1)) ! implicit
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  bb(:,kl)=1.-aa(:,kl)
  dd(:,1:kl)=uo(1:ifull,1:kl)
  ! bb(:,1)=1.-cc(:,1)                                           ! explicit
  ! dd(:,1:kl)=uo(1:ifull,1:kl)-ddts*taux/(rhoa(:,1)*dz_fl(:,1)) ! explicit
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  dd(:,1:kl)=vo(1:ifull,1:kl)
  ! dd(:,1:kl)=vo(1:ifull,1:kl)-ddts*tauy/(rhoa(:,1)*dz_fl(:,1)) ! explicit
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))

  ! umag=sqrt(max(uo(1:ifull,1)*uo(1:ifull,1)+vo(1:ifull,1)*vo(1:ifull,1),1.e-4)) ! explicit
  ! call dyerhicks(cdrag,wtv0,zom,umag,thetav(:,1),zz(:,1))                       ! explicit
  ! taux=rhos*cdrag*umag*uo(1:ifull,1)                                            ! explicit
  ! tauy=rhos*cdrag*umag*vo(1:ifull,1)                                            ! explicit
  ! ustar=sqrt(sqrt(taux*taux+tauy*tauy)/rhos)                                    ! explicit
  
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
real, dimension(0:220), save :: table
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix
logical, save :: first=.true.

if (first) then
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
  first=.false.
end if

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

kmo(1:ifull,1:kl-1)=km(1:ifull,1:kl-1)+fzhl(:,1:kl-1)*(km(1:ifull,2:kl)-km(1:ifull,1:kl-1))
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
entfn = ent0/max( zht, 1. ) + ent0/max( zi-zht, ezmin )

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral detrainment

real function dtrfn(zht,zi,rat)

implicit none

real, intent(in) :: zht,zi,rat

!dtrfn=ent+0.05/max(zi-zht,zmin)   ! Angevine et al (2010)
dtrfn=rat/max(zi-zht,1.)+ent0/max(zi-zht,ezmin)

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

ilzom      = log(zmin/zom)
ustar      = vkar*max(umag,1.e-2)/ilzom ! first guess

select case(stabmeth)
  case(0)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*zmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/zmin
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
      ustar = vkar*max(umag,1.e-2)/integralm
    end do
    
  case(1)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*zmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/zmin
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
        ustar = vkar*max(umag,1.e-2)/integralm
      elsewhere ! Luhar
        ustar = vkar*max(umag,1.e-2)/(aa1*(( z_on_l**bb1)*(1.+ cc1*z_on_l**(1.-bb1)) &
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
deallocate(shear,zidry)

return
end subroutine tkeend

end module tkeeps

