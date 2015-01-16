
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
public mintke,mineps,cm0,cq,minl,maxl
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

! model constants
real, parameter :: be      = 0.1    ! Hurley (2007) 1., Soares et al (2004) 0.3
real, parameter :: cm0     = 0.09   ! Hurley (2007) 0.09, Duynkerke 1988 0.03, Duynkerke 1987 0.09
real, parameter :: ce0     = 0.69   ! Hurley (2007) 0.69, Duynkerke 1988 0.42, Duynkerke 1987 0.77
real, parameter :: ce1     = 1.46
real, parameter :: ce2     = 1.83
real, parameter :: ce3     = 0.45   ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, parameter :: cq      = 2.5    ! Adjustment to ED in absence of MF
real, parameter :: ent0    = 0.25   ! Entrainment constant (Controls height of boundary layer)
real, parameter :: dtrn0   = 0.4    ! Unsaturated detrainment constant
real, parameter :: dtrc0   = 0.9    ! Saturated detrainment constant
real, parameter :: m0      = 0.06   ! Mass flux amplitude constant
real, parameter :: b1      = 2.     ! Soares et al (2004) 1., Siebesma et al (2003) 2.
real, parameter :: b2      = 1./3.  ! Soares et al (2004) 2., Siebesma et al (2003) 1./3.

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

integer, parameter :: icm1   = 5        ! max iterations for calculating pblh
real, parameter :: maxdts    = 120.     ! max timestep for split
real, parameter :: mintke    = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, parameter :: mineps    = 1.E-10   ! min value for eps (1.0e-6 in TAPM)
real, parameter :: minl      = 1.       ! min value for L   (5. in TAPM)
real, parameter :: maxl      = 1000.    ! max value for L   (500. in TAPM)

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

subroutine tkemix(kmo,theta,qvg,qlg,qfg,qrg,cfrac,cfrain,uo,vo,zi,fg,eg,ps, &
                  zom,zz,zzh,sig,rhos,dt,qgmin,mode,diag,naero,aero)

implicit none

integer, intent(in) :: diag,mode,naero
integer k,i,j,ktopmax,klcl
integer kcount,mcount,icount
real, intent(in) :: dt,qgmin
real, dimension(:,:,:), intent(inout) :: aero
real, dimension(:,:), intent(inout) :: theta,cfrac,cfrain,uo,vo
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg,qrg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl), intent(in) :: zz,zzh
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: fg,eg,ps,zom,rhos
real, dimension(kl), intent(in) :: sig
real, dimension(ifull,kl,naero) :: gamar
real, dimension(ifull,kl) :: km,thetav,thetal,temp,qsat
real, dimension(ifull,kl) :: gamtv,gamtl,gamqv,gamql,gamqf
real, dimension(ifull,kl) :: gamqr,gamhl,gamth,gamcf,gamcr
real, dimension(ifull,kl) :: thetavnc,qsatc,thetac,tempc
real, dimension(ifull,kl) :: tkenew,epsnew,bb,cc,dd,ff,rr
real, dimension(ifull,kl) :: rhoa,rhoahl,thetavhl,thetahl
real, dimension(ifull,kl) :: qshl,qlhl,qfhl
real, dimension(ifull,kl) :: pres
real, dimension(ifull,2:kl) :: idzm
real, dimension(ifull,1:kl-1) :: idzp
real, dimension(ifull,2:kl) :: aa,qq,pps,ppt,ppb
real, dimension(ifull,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,kl-1) :: fzzh
real, dimension(ifull,naero) :: arup
real, dimension(ifull) :: wt0,wq0,wtv0
real, dimension(ifull) :: wstar,z_on_l,phim
real, dimension(ifull) :: tff,tbb,tcc,tgg,tqq,qgnc,dum
real, dimension(ifull) :: cdrag,umag,ustar
real, dimension(kl) :: sigkap
real, dimension(kl) :: w2up,nn,dqdash
real, dimension(kl) :: qtup,qupsat,ttup,tvup,tlup,thup
real, dimension(kl) :: hup,qvup,qlup,qfup,qrup
real, dimension(kl) :: cfup,crup,mflx
real, dimension(1) :: templ
real xp,as,bs,cs,cm12,cm34,qcup
real dzht,ziold,ent,entc,entn,dtr,dtrc,dtrn,dtrx
real ddts,zlcl
real lx,tempd,fice,qxup,dqsdt,al
real sigqtup,rng
logical, dimension(ifull,kl) :: lta
logical scond

cm12=1./sqrt(cm0)
cm34=sqrt(sqrt(cm0**3))

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
  thetav(:,k)=theta(1:ifull,k)*(1.+0.61*qvg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k))
  rhoa(:,k)=sigkap(k)*pres(:,k)/(rd*thetav(:,k))

  ! Transform to thetal as it is the conserved variable
  thetal(:,k)=theta(1:ifull,k)-sigkap(k)*(lv*(qlg(1:ifull,k)+qrg(1:ifull,k))+ls*qfg(1:ifull,k))/cp
  
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
  wtv0  =wt0+theta(1:ifull,1)*0.61*wq0 ! thetav flux
  umag=sqrt(max(uo(1:ifull,1)*uo(1:ifull,1)+vo(1:ifull,1)*vo(1:ifull,1),1.e-4))
  call dyerhicks(cdrag,wtv0,zom,umag,thetav(:,1),zz(:,1))
  ustar=sqrt(cdrag)*umag
    
  ! Set-up thermodynamic variables temp, theta_v and surface fluxes
  do k=1,kl
    temp(:,k)=theta(1:ifull,k)/sigkap(k)
    ! calculate saturated air mixing ratio
    call getqsat(qsat(:,k),temp(:,k),pres(:,k))
  end do
  thetav=theta(1:ifull,:)*(1.+0.61*qvg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:)-qrg(1:ifull,:))

  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  gamtl=0.
  gamth=0.
  gamtv=0.
  gamqv=0.
  gamql=0.
  gamqf=0.
  gamqr=0.
  gamcf=0.
  gamcr=0.
  if (naero>0) then
    gamar=0.
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
          cfup(:)=0.
          scond=.false.
          dzht=zz(i,1)
          ! Entrainment and detrainment rates
          ent=entfn(zz(i,1),zi(i),zz(i,1))
          
          ! first level -----------------
          ! initial thermodynamic state
          ! split qtot into components (conservation of thetal and qtot is maintained)
          tlup(1)=thetal(i,1)+be*wt0(i)/sqrt(tke(i,1))   ! Hurley 2007
          qvup(1)=qvg(i,1)   +be*wq0(i)/sqrt(tke(i,1))   ! Hurley 2007
          qlup(1)=qlg(i,1)
          qfup(1)=qfg(i,1)
          qrup(1)=qrg(i,1)
          ! diagnose thermodynamic variables assuming no condensation
          qtup(1)=qvup(1)+qlup(1)+qfup(1)+qrup(1)                          ! qtot,up
          ! state of plume after evaporation
          qvup(1)=qtup(1)
          qlup(1)=0.
          qfup(1)=0.
          qrup(1)=0.
          thup(1)=tlup(1) !+sigkap(1)*(lv*(qlup(1)+qrup(1))+ls*qfup(1))/cp ! theta,up
          tvup(1)=thup(1)+theta(i,1)*0.61*qtup(1) ! +...                   ! thetav,up
          ttup(1)=thup(1)/sigkap(1)                                        ! temp,up
          templ(1)=tlup(1)/sigkap(1)                                       ! templ,up
          ! update updraft velocity and mass flux
          nn(1)  =grav*be*wtv0(i)/(thetav(i,1)*sqrt(tke(i,1)))             ! Hurley 2007
          w2up(1)=2.*dzht*b2*nn(1)/(1.+2.*dzht*b1*ent)                     ! Hurley 2007
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
            ent=entfn(zz(i,k),zi(i),zz(i,1))
            ! update thermodynamics of plume
            ! split qtot into components (conservation is maintained)
            tlup(k)=(tlup(k-1)+dzht*ent*thetal(i,k))/(1.+dzht*ent)
            qvup(k)=(qvup(k-1)+dzht*ent*qvg(i,k)   )/(1.+dzht*ent)
            qlup(k)=(qlup(k-1)+dzht*ent*qlg(i,k)   )/(1.+dzht*ent)
            qfup(k)=(qfup(k-1)+dzht*ent*qfg(i,k)   )/(1.+dzht*ent)
            qrup(k)=(qrup(k-1)+dzht*ent*qrg(i,k)   )/(1.+dzht*ent)
            ! calculate conserved variables
            qtup(k)=qvup(k)+qlup(k)+qfup(k)+qrup(k)                         ! qtot,up
            ! estimate air temperature
            templ(1)=tlup(k)/sigkap(k)                                      ! templ,up
            if (.not.scond) then
              ! MJT notes - here we determine a fraction of the updraft that has undergone condensation
              ! More correctly, we estimate a fraction of the number of updrafts at different locations
              ! in the grid-box that have undergone condesation.  For this reason, the variance of
              ! the moisture between updrafts is estimated from the grid-box mean.
              call getqsat(qupsat(k:k),templ(1:1),pres(i:i,k))
              ! estimate variance of qtup in updraft (following Hurley and TAPM)
              sigqtup=sqrt(max(1.E-10,1.6*tke(i,k)/eps(i,k)*cq*km(i,k)*((qtup(k)-qtup(k-1))/dzht)**2))
              ! MJT condensation scheme - follow Smith 1990 and assume triangle distribution for qtup.
              ! The average qtup is qxup after accounting for saturation
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
            qvup(k)=qtup(k)
            qlup(k)=0.
            qfup(k)=0.
            qrup(k)=0.
            thup(k)=tlup(k) !+sigkap(k)*(lv*(qlup(k)+qrup(k))+ls*qfup(k))/cp ! theta,up
            tvup(k)=tlup(k)+theta(i,k)*0.61*qtup(k) ! +...                   ! thetav,up after redistribution
            ttup(k)=thup(k)/sigkap(k)                                        ! temp,up
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
              ent=entfn(zz(i,k),zi(i),zz(i,1))
              ! update thermodynamics of plume
              ! split qtot into components (conservation is maintained)
              tlup(k)=(tlup(k-1)+dzht*ent*thetal(i,k))/(1.+dzht*ent)
              qvup(k)=(qvup(k-1)+dzht*ent*qvg(i,k)   )/(1.+dzht*ent)
              qlup(k)=(qlup(k-1)+dzht*ent*qlg(i,k)   )/(1.+dzht*ent)
              qfup(k)=(qfup(k-1)+dzht*ent*qfg(i,k)   )/(1.+dzht*ent)
              qrup(k)=(qrup(k-1)+dzht*ent*qrg(i,k)   )/(1.+dzht*ent)
              ! calculate conserved variables
              qtup(k)=qvup(k)+qlup(k)+qfup(k)+qrup(k)                         ! qtot,up
              ! estimate air temperature
              templ(1)=tlup(k)/sigkap(k)                                      ! templ,up
              call getqsat(qupsat(k:k),templ(1:1),pres(i:i,k))
              ! estimate variance of qtup in updraft (following Hurley and TAPM)
              sigqtup=sqrt(max(1.E-10,1.6*tke(i,k)/eps(i,k)*cq*km(i,k)*((qtup(k)-qtup(k-1))/dzht)**2))
              ! MJT condensation scheme -  follow Smith 1990 and assume
              ! triangle distribution for qtup.  The average qtup is qxup
              ! after accounting for saturation
              rng=sqrt(6.)*sigqtup               ! variance of triangle distribution
              dqdash(k)=(qtup(k)-qupsat(k))/rng  ! scaled variance
              if (dqdash(k)<-1.) then
                ! gridbox all unsaturated
                qxup=qtup(k)
                cfup(k)=0.
              else if (dqdash(k)<0.) then
                ! gridbox minority saturated
                qxup=qtup(k)+0.5*rng*(-1./3.-dqdash(k)-dqdash(k)**2-1./3.*dqdash(k)**3)
                cfup(k)=0.5*(dqdash(k)+1.)**2
              else if (dqdash(k)<1.) then
                ! gridbox majority saturated
                qxup=qtup(k)+0.5*rng*(-1./3.-dqdash(k)-dqdash(k)**2+1./3.*dqdash(k)**3)
                cfup(k)=1.-0.5*(dqdash(k)-1.)**2
              else
                ! gridbox all saturated
                qxup=qupsat(k)
                cfup(k)=1.
              end if
              thup(k)=tlup(k)+sigkap(k)*(lv*(qlup(k)+qrup(k))+ls*qfup(k))/cp  ! theta,up before redistribution
              tempd  =thup(k)/sigkap(k)                                       ! temp,up before redistribution
              fice=min(max(273.16-tempd,0.),40.)/40. ! approximate ice fraction based on temperature (not templ)
              lx=lv+lf*fice
              ! MJT notes - PC2 from UKMO advises that the following expression could be based on tempd instead of templ
              ! However, Smith (1990) suggests using templ.
              dqsdt=qupsat(k)*lx/(rv*templ(1)*templ(1))
              al=cp/(cp+lx*dqsdt)
              qcup=al*(qtup(k)-qxup)                                 ! qcondensate,up after redistribution
              qxup=qtup(k)-qcup                                      ! qv,up after redistribution
              ttup(k)=templ(1)+lx*qcup/cp                            ! temp,up after redistribution
              thup(k)=ttup(k)*sigkap(k)                              ! theta,up after redistribution
              tvup(k)=thup(k)+theta(i,k)*(1.61*qxup-qtup(k))         ! thetav,up after redistribution
              ! state of plume after redistribution
              qvup(k)=qxup                                           ! qv,up after redistribution
              qlup(k)=qcup*(1.-fice)                                 ! ql,up after redistribution
              qfup(k)=qcup*fice                                      ! qf,up after redistribution
              qrup(k)=0.                                             ! qr,up after redistribution
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
          if (abs(zidry(i)-ziold)<1.) exit
        end do

        ! update mass flux
        mflx(1)=m0*sqrt(max(w2up(1),0.))*zz(i,1)**ent0*max(zi(i)-zz(i,1),1.)**dtrn0 ! MJT suggestion
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          xp  =(zz(i,k)-zlcl)/max(zidry(i)-zlcl,0.1)
          xp  =min(max(xp,0.),1.)
          ent =entfn(zz(i,k),zi(i),   zz(i,1))
          dtrn=dtrfn(zz(i,k),zidry(i),zz(i,1),dtrn0)
          dtrc=dtrfn(zz(i,k),zi(i),   zz(i,1),dtrn0)
          dtrx=(1.-xp)*dtrn+xp*dtrc
          dtrc=dtrfn(zz(i,k),zi(i),   zz(i,1),dtrc0)
          dtr =(1.-cfup(k))*dtrx+cfup(k)*dtrc
          mflx(k)=mflx(k-1)/(1.+dzht*(dtr-ent))
        end do

#ifdef offline
        do k=1,ktopmax
          mf(i,k)=mflx(k)
          w_up(i,k)=sqrt(w2up(k))
          tl_up(i,k)=tlup(k)
          qv_up(i,k)=qvup(k)
          ql_up(i,k)=qlup(k)
          qf_up(i,k)=qfup(k)
          cf_up(i,k)=cfup(k)*min(mflx(k)/sqrt(w2up(k)),1.)
        end do
          
        dzht=zz(i,1)
        ! Entrainment and detrainment rates
        xp  =(zz(i,1)-zlcl)/max(zidry(i)-zlcl,0.1)
        xp  =min(max(xp,0.),1.)
        ent =entfn(zz(i,1),zi(i),   zz(i,1))
        dtrn=dtrfn(zz(i,1),zidry(i),zz(i,1),dtrn0)
        dtrc=dtrfn(zz(i,1),zi(i),   zz(i,1),dtrn0)
        dtrx=(1.-xp)*dtrn+xp*dtrc
        dtrc=dtrfn(zz(i,1),zi(i),   zz(i,1),dtrc0)
        dtr =(1.-cfup(1))*dtrx+cfup(1)*dtrc
        ents(i,1)=ent
        dtrs(i,1)=dtr
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          ! Entrainment and detrainment rates
          xp  =(zz(i,k)-zlcl)/max(zidry(i)-zlcl,0.1)
          xp  =min(max(xp,0.),1.)
          ent =entfn(zz(i,k),zi(i),   zz(i,1))
          dtrn=dtrfn(zz(i,k),zidry(i),zz(i,1),dtrn0)
          dtrc=dtrfn(zz(i,k),zi(i),   zz(i,1),dtrn0)
          dtrc=dtrfn(zz(i,k),zi(i),   zz(i,1),dtrc0)
          dtrx=(1.-xp)*dtrn+xp*dtrc
          dtr =(1.-cfup(k))*dtrx+cfup(k)*dtrc
          ents(i,k)=ent
          dtrs(i,k)=dtr
        end do
#endif

        ! update explicit counter gradient terms
        do k=1,ktopmax
          gamtl(i,k)=mflx(k)*(tlup(k)-thetal(i,k))
          gamqv(i,k)=mflx(k)*(qvup(k)-qvg(i,k) )
          gamql(i,k)=mflx(k)*(qlup(k)-qlg(i,k))
          gamqf(i,k)=mflx(k)*(qfup(k)-qfg(i,k))
          gamqr(i,k)=mflx(k)*(qrup(k)-qrg(i,k))
          gamth(i,k)=gamtl(i,k)+sigkap(k)*(lv*(gamql(i,k)+gamqr(i,k))+ls*gamqf(i,k))/cp
          gamtv(i,k)=gamth(i,k)+theta(i,k)*(0.61*gamqv(i,k)-gamql(i,k)-gamqf(i,k)-gamqr(i,k))
        end do

        ! update reamining scalars which are not used in the iterative loop
        crup(1)=cfrain(i,1)
        do j=1,naero
          arup(1,j)=aero(i,1,j)
        end do
        gamcf(i,1)=0.
        gamcr(i,1)=0.
        do j=1,naero
          gamar(i,1,j)=0.
        end do
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          ent=entfn(zz(i,k),zi(i),zz(i,1))
          crup(k)=(crup(k-1)+dzht*ent*cfrain(i,k))/(1.+dzht*ent)
          do j=1,naero
            arup(k,j)=(arup(k-1,j)+dzht*ent*aero(i,k,j))/(1.+dzht*ent)
          end do
          gamcf(i,k)=mflx(k)*(cfup(k)-cfrac(i,k))
          gamcr(i,k)=mflx(k)*(crup(k)-cfrain(i,k))
          do j=1,naero
            gamar(i,k,j)=mflx(k)*(arup(k,j)-aero(i,k,j))
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

  ! calculate tke and eps at 1st level
  z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-10))
  z_on_l=min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  where (z_on_l<0.)
    phim=(1.-16.*z_on_l)**(-0.25)
  elsewhere
    phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
  end where
  tke(1:ifull,1)=cm12*ustar*ustar+ce3*wstar*wstar
  eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
  tke(1:ifull,1)=max(tke(1:ifull,1),mintke)
  tff=cm34*tke(1:ifull,1)*sqrt(tke(1:ifull,1))
  eps(1:ifull,1)=min(eps(1:ifull,1),tff/minl)
  eps(1:ifull,1)=max(eps(1:ifull,1),tff/maxl,mineps)


  ! Calculate sources and sinks for TKE and eps
  ! prepare arrays for calculating buoyancy of saturated air
  ! (i.e., related to the saturated adiabatic lapse rate)
  qsatc=max(qsat,qvg(1:ifull,:))                                             ! assume qvg is saturated inside cloud
  ff=qfg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                              ! inside cloud value
  dd=qlg(1:ifull,:)/max(cfrac(1:ifull,:),1.E-8)                            &
    +qrg(1:ifull,:)/max(cfrac(1:ifull,:),cfrain(1:ifull,:),1.E-8)            ! inside cloud value assuming max overlap
  do k=1,kl
    tbb=max(1.-cfrac(1:ifull,k),1.E-8)
    qgnc=(qvg(1:ifull,k)-(1.-tbb)*qsatc(:,k))/tbb                            ! outside cloud value
    qgnc=min(max(qgnc,qgmin),qsatc(:,k))
    thetac(:,k)=thetal(:,k)+sigkap(k)*(lv*dd(:,k)+ls*ff(:,k))/cp             ! inside cloud value
    tempc(:,k)=thetac(:,k)/sigkap(k)                                         ! inside cloud value
    thetavnc(:,k)=thetal(:,k)*(1.+0.61*qgnc)                                 ! outside cloud value
  end do
  call updatekmo(thetahl,thetac,fzzh)                                        ! inside cloud value
  call updatekmo(thetavhl,thetavnc,fzzh)                                     ! outside cloud value
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

  ! Update TKE and eps terms

  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tke(1:ifull,kl)=mintke
  eps(1:ifull,kl)=mineps

  do k=2,kl-1
    ! Calculate buoyancy term
    ! saturated from Durran and Klemp JAS 1982 (see also WRF)
    tqq=(1.+lv*qsatc(:,k)/(rd*tempc(:,k)))/(1.+lv*lv*qsatc(:,k)/(cp*rv*tempc(:,k)*tempc(:,k)))
    tbb=-grav*km(:,k)*(tqq*((thetahl(:,k)-thetahl(:,k-1))/thetac(:,k)                              &
           +lv*(qshl(:,k)-qshl(:,k-1))/(cp*tempc(:,k)))-qshl(:,k)-qlhl(:,k)-qfhl(:,k)              &
           +qshl(:,k-1)+qlhl(:,k-1)+qfhl(:,k-1))/dz_fl(:,k)
    tbb=tbb+grav*(tqq*(gamtl(:,k)/thetac(:,k)+lv*gamqv(:,k)/(cp*tempc(:,k)))                       &
           -gamqv(:,k)-gamql(:,k)-gamqf(:,k)-gamqr(:,k))
    ! unsaturated
    tcc=-grav*km(:,k)*(thetavhl(:,k)-thetavhl(:,k-1))/(thetavnc(:,k)*dz_fl(:,k))
    tcc=tcc+grav*gamtv(:,k)/thetavnc(:,k)
    ppb(:,k)=(1.-cfrac(1:ifull,k))*tcc+cfrac(1:ifull,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)

    ! Calculate transport source term on full levels
    ppt(:,k)= kmo(:,k)*idzp(:,k)*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k)                                  &
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
  call updatekmo(gamhl,gamtl,fzzh)
  cc(:,1)     =rr(:,1)
  bb(:,1)     =1.-rr(:,1)
  aa(:,2:kl-1)=qq(:,2:kl-1)
  cc(:,2:kl-1)=rr(:,2:kl-1)
  bb(:,2:kl-1)=1.-qq(:,2:kl-1)-rr(:,2:kl-1)
  aa(:,kl)    =qq(:,kl)
  bb(:,kl)    =1.-qq(:,kl)
  dd(:,1)     =thetal(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)+ddts*rhos*wt0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=thetal(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =thetal(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wthl(:,1:kl-1)=-kmo(:,1:kl-1)*(thetal(1:ifull,2:kl)-thetal(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)+gamhl(:,1:kl-1)
#endif

  call updatekmo(gamhl,gamqv,fzzh)
  dd(:,1)     =qvg(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)+ddts*rhos*wq0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=qvg(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =qvg(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wqv(:,1:kl-1)=-kmo(:,1:kl-1)*(qvg(1:ifull,2:kl)-qvg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)+gamhl(:,1:kl-1)
#endif

  call updatekmo(gamhl,gamql,fzzh)
  dd(:,1)     =qlg(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qlg(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =qlg(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wql(:,1:kl-1)=-kmo(:,1:kl-1)*(qlg(1:ifull,2:kl)-qlg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)+gamhl(:,1:kl-1)
#endif

  call updatekmo(gamhl,gamqf,fzzh)
  dd(:,1)     =qfg(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qfg(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =qfg(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
#ifdef offline
  wqf(:,1:kl-1)=-kmo(:,1:kl-1)*(qfg(1:ifull,2:kl)-qfg(1:ifull,1:kl-1))/dz_hl(:,1:kl-1)+gamhl(:,1:kl-1)
#endif

  call updatekmo(gamhl,gamqr,fzzh)
  dd(:,1)     =qrg(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qrg(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =qrg(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qrg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  
  ! account for phase transitions
  do k=1,kl
    tgg=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k)+qrg(1:ifull,k),qgmin) ! qtot before phase transition
    qlg(1:ifull,k)=max(qlg(1:ifull,k),0.)
    qfg(1:ifull,k)=max(qfg(1:ifull,k),0.)
    qrg(1:ifull,k)=max(qrg(1:ifull,k),0.)
    qvg(1:ifull,k)=max(qvg(1:ifull,k),0.)
    tff=max(qvg(1:ifull,k)+qlg(1:ifull,k)+qfg(1:ifull,k)+qrg(1:ifull,k),qgmin)
    tgg=tgg/tff ! scale factor for conservation
    qvg(1:ifull,k)=qvg(1:ifull,k)*tgg
    qlg(1:ifull,k)=qlg(1:ifull,k)*tgg
    qfg(1:ifull,k)=qfg(1:ifull,k)*tgg
    qrg(1:ifull,k)=qrg(1:ifull,k)*tgg
    ! update theta for output or next time step
    theta(1:ifull,k)=thetal(:,k)+sigkap(k)*(lv*(qlg(1:ifull,k)+qrg(1:ifull,k))+ls*qfg(1:ifull,k))/cp
  end do

  ! update cloud fraction terms
  call updatekmo(gamhl,gamcf,fzzh)
  dd(:,1)     =cfrac(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=cfrac(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =cfrac(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(cfrac,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfrac(1:ifull,:)=min(max(cfrac(1:ifull,:),0.),1.)
  where (qlg(1:ifull,:)+qfg(1:ifull,:)>1.E-12)
    cfrac(1:ifull,:)=max(cfrac(1:ifull,:),1.E-8)
  end where

  call updatekmo(gamhl,gamcr,fzzh)
  dd(:,1)     =cfrain(1:ifull,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=cfrain(1:ifull,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)    =cfrain(1:ifull,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(cfrain,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  cfrain(1:ifull,:)=min(max(cfrain(1:ifull,:),0.),1.)
  where (qrg(1:ifull,:)>1.E-12)
    cfrain(1:ifull,:)=max(cfrain(1:ifull,:),1.E-8)
  end where
  
  ! Aerosols
  do j=1,naero
    call updatekmo(gamhl,gamar(:,:,j),fzzh)
    dd(:,1)     =aero(1:ifull,1,j)-ddts*gamhl(:,1)*idzp(:,1)
    dd(:,2:kl-1)=aero(1:ifull,2:kl-1,j)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
    dd(:,kl)    =aero(1:ifull,kl,j)+ddts*gamhl(:,kl-1)*idzm(:,kl)
    call thomas(aero(:,:,j),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  end do

  ! Winds
  aa(:,2:kl)  =qq(:,2:kl)
  cc(:,1:kl-1)=rr(:,1:kl-1)
  bb(:,1)=1.-cc(:,1)+ddts*rhos*cdrag*umag/(rhoa(:,1)*dz_fl(:,1)) 
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  bb(:,kl)=1.-aa(:,kl)
  dd(:,1:kl)=uo(1:ifull,1:kl)
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))
  dd(:,1:kl)=vo(1:ifull,1:kl)
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl))

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

real, dimension(ifull,kl), intent(in) :: km
real, dimension(ifull,kl-1), intent(in) :: fzhl
real, dimension(ifull,kl), intent(out) :: kmo

kmo(:,1:kl-1)=km(:,1:kl-1)+fzhl(:,1:kl-1)*(km(:,2:kl)-km(:,1:kl-1))
! These terms are never used
kmo(:,kl)=0.

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

real function entfn(zht,zi,zmin)

implicit none

real, intent(in) :: zht,zi,zmin

!entfn=0.002                                               ! Angevine (2005)
!entfn=2./max(100.,zi)                                     ! Angevine et al (2010)
!entfn=1./zht                                              ! Siebesma et al (2003)
!entfn=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))       ! Soares et al (2004)
entfn=ent0*(1./max(zht,1.)+1./max(zi-zht,100.))

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral detrainment

real function dtrfn(zht,zi,zmin,rat)

implicit none

real, intent(in) :: zht,zi,zmin,rat

!dtrfn=ent+0.05/max(zi-zht,zmin)   ! Angevine et al (2010)
dtrfn=rat/max(zi-zht,10.)+ent0/max(zi-zht,100.)

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

