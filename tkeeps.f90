
! This module calculates the turblent kinetic energy and mixing for the boundary layer based on Hurley 2009
! (eddy dissipation) and Angevine et al 2010 (mass flux).  Specifically, this version is adapted for CCAM
! with an implicit numerical formulation (i.e., a forward-backward scheme).

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
public tkeinit,tkemix,tkeend,tke,eps,tkesav,epssav,shear

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:,:), allocatable, save :: tkesav,epssav

! model constants
real, parameter :: b1      = 1.     ! Soares et al (2004) 1., Siebesma et al (2003) 2.
real, parameter :: b2      = 2.     ! Soares et al (2004) 2., Siebesma et al (2003) 1./3.
real, parameter :: be      = 1.     ! Hurley (2007) 1., Soares et al (2004) 0.3
real, parameter :: cm      = 0.09
real, parameter :: ce0     = 0.69
real, parameter :: ce1     = 1.46
real, parameter :: ce2     = 1.83
real, parameter :: cq      = 2.5


! physical constants
real, parameter :: grav  = 9.80616
real, parameter :: lv    = 2.5104e6
real, parameter :: lf    = 3.36e5
real, parameter :: ls    = lv+lf
real, parameter :: rd    = 287.04
real, parameter :: rv    = 461.5
real, parameter :: epsl  = rd/rv
real, parameter :: delta = 1./(epsl-1.)
real, parameter :: cp    = 1004.64
real, parameter :: vkar  = 0.4
real, parameter :: pi    = 3.1415927

! stability constants
real, parameter :: a_1   = 1.
real, parameter :: b_1   = 2./3.
real, parameter :: c_1   = 5.
real, parameter :: d_1   = 0.35
!real, parameter :: aa1 = 3.8 ! Luhar low wind
!real, parameter :: bb1 = 0.5 ! Luhar low wind
!real, parameter :: cc1 = 0.3 ! Luhar low wind

integer, parameter :: buoymeth   = 3 ! 0=Hurley (dry), 2=Smith, 3=Durran et al (cld wgt) method for calculating TKE buoyancy source
integer, parameter :: kmlimit    = 0 ! 0=No adjustment, 1=Limit decay in pbl top
integer, parameter :: icm  = 30      ! iterations for calculating pblh
integer, parameter :: icm2 = 10      ! iterations for calculating tke-eps
integer, parameter :: icm3 = 2       ! iterations for calculating kmo
real, parameter :: tol    = 1.E-9    ! tolarance for sectant method
real, parameter :: dz_ref = 500.     ! dz reference for kmlimit=1
real, parameter :: alpha  = 0.3      ! weight for updating pblh
real, parameter :: alpha2 = 0.7      ! weight for updating tke-eps non-linear terms
real, parameter :: beta1  = 1.0      ! weight for km(t+1) (beta=0.) and km(t) (beta=1.) when calculating tke non-linear terms
real, parameter :: beta2  = 0.0      ! weight for km(t+1) (beta=0.) and km(t) (beta=1.) when calculating eps non-linear terms
real, parameter :: mintke = 1.5E-8
real, parameter :: mineps = 1.0E-10

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifullin,iextrain,klin,diag)

implicit none

integer, intent(in) :: ifullin,iextrain,klin,diag

if (diag.gt.0) write(6,*) "Initialise TKE scheme"

ifull=ifullin
iextra=iextrain
kl=klin

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(tkesav(ifull,kl),epssav(ifull,kl))
allocate(shear(ifull,kl))

tke=mintke
eps=mineps
tkesav=tke(1:ifull,:)
epssav=eps(1:ifull,:)
shear=0.

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qg,qlg,qfg,cfrac,zi,wt0,wq0,ps,ustar,zz,sig,sigkap,dt,mode,diag)

implicit none

integer, intent(in) :: diag,mode
integer k,kk,i,klcl
integer icount,jcount,mc
real, intent(in) :: dt
real, dimension(ifull,kl), intent(inout) :: theta,qg
real, dimension(ifull,kl), intent(in) :: zz,cfrac,qlg,qfg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar
real, dimension(kl), intent(in) :: sigkap,sig
real, dimension(ifull,kl) :: km,ff,gg,templ,thetav,thetal,temp,qtot
real, dimension(ifull,kl) :: gamtt,gamtl,gamtv,gamth,gamqt,gamqv
real, dimension(ifull,kl) :: tkenew,epsnew,qsat
real, dimension(ifull,kl) :: bb,cc,dd
real, dimension(ifull,kl) :: mflx,ttup,tlup,tvup,thup,qtup,qvup
real, dimension(ifull,2:kl) :: aa,pps,ppt,ppb
real, dimension(ifull,kl) :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,kl-1) :: gam_hl
real, dimension(ifull) :: wstar,z_on_l,phim,hh,jj,dqsdt,ziold
real, dimension(ifull) :: dum,wtv0
real, dimension(kl-1) :: wpv_flux
real, dimension(kl) :: w2up
real, dimension(kl) :: qupsat,newsat,pres
real xp,ee,dtr,nn,as,bs,cs,cm34,qlup
real zht,dzht,zidry,zilcl
logical sconv

cm34=cm**0.75

if (diag.gt.0) write(6,*) "Update PBL mixing with TKE"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

! save current values of tke and eps
tkesav=tke(1:ifull,:)
epssav=eps(1:ifull,:)

! impose limits after host advection
tke(1:ifull,:)=max(tke(1:ifull,:),mintke)
tke(1:ifull,:)=min(tke(1:ifull,:),65.)
ff=cm34*(tke(1:ifull,:)**1.5)/5.
eps(1:ifull,:)=min(eps(1:ifull,:),ff)
ff=max(ff*5./500.,mineps)
eps(1:ifull,:)=max(eps(1:ifull,:),ff)

do k=1,kl
  ! calculate saturated mixing ratio
  temp(:,k)=theta(:,k)/sigkap(k)
  dum=ps*sig(k)
  call getqsat(ifull,qsat(:,k),temp(:,k),dum)
  ! prepare arrays for calculating buoyancy of saturated air
  ! (i.e., related to the saturated adiabatic lapse rate)
  gg(:,k)=(1.+lv*qsat(:,k)/(rd*temp(:,k)))/(1.+lv*lv*qsat(:,k)/(cp*rv*temp(:,k)*temp(:,k)))
end do

! Set-up thermodynamic variables theta_l,theta_v,qtot and surface fluxes
! Note that qfg is treated as passive
qtot=qg+qlg+qfg
thetav=theta*(1.+1.61*qg-qtot)
do k=1,kl
  thetal(:,k)=theta(:,k)-lv*sigkap(k)*qlg(:,k)/cp
end do
wtv0=wt0+0.61*theta(:,1)*wq0
!wtl0=wt0
!wqt0=wq0

! Calculate dz at full levels
dz_fl(:,1)=0.5*(zz(:,2)+zz(:,1))
do k=2,kl-1
  dz_fl(:,k)=0.5*(zz(:,k+1)-zz(:,k-1))
end do
dz_fl(:,kl)=zz(:,kl)-zz(:,kl-1)

! Calculate dz at half levels
do k=1,kl-1
  dz_hl(:,k)=zz(:,k+1)-zz(:,k)
end do

! Calculate first approximation to diffusion coeffs
km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-10)
call updatekmo(kmo,km,zz) ! interpolate diffusion coeffs to half levels


! Calculate non-local mass-flux terms for theta_l and qtot
mflx=0.
gamtl=0.
gamtv=0.
gamth=0.
gamqt=0.
gamqv=0.
tlup=thetal
tvup=thetav
thup=theta
ttup=temp
qtup=qtot
qvup=qg
ziold=zi
wstar=0.
where (wtv0.gt.0.)
  wstar=(grav*zi*wtv0/thetav(:,1))**(1./3.)
end where
if (mode.ne.1) then ! mass flux when mode is an even number
  do i=1,ifull
    if (wtv0(i).gt.0.) then ! unstable
      pres(:)=ps(i)*sig(:) ! pressure
      ! initial guess for plume state
      qupsat(:)=qsat(i,:)
      newsat(:)=qsat(i,:)
      zidry=zi(i)
      do icount=1,icm
        klcl=kl+1
        mflx(i,:)=0.
        ! first level -----------------
        ! dry convection case
        zht=zz(i,1)
        ! estimate entrainment and detrainment rates
        ee=2./max(100.,zidry)
        dtr=ee+0.05/max(zidry-zht,1.)
        ! initial thermodynamic state
        tlup(i,1)=thetal(i,1)+be*wt0(i)/sqrt(tke(i,1))
        qtup(i,1)=qtot(i,1)  +be*wq0(i)/sqrt(tke(i,1))
        ! diagnose thermodynamic variables assuming no condensation
        thup(i,1)=tlup(i,1)                                       ! theta,up
        tvup(i,1)=thup(i,1)+theta(i,1)*0.61*qtup(i,1)             ! thetav,up
        ttup(i,1)=thup(i,1)/sigkap(1)                             ! temp,up
        call getqsat(1,newsat(1),ttup(i,1),pres(1))               ! estimate of saturated mixing ratio in plume
        ! update vertical velocity
        w2up(1)=(0.5*wstar(i))**2
        ! update mass flux
        mflx(i,1)=0.03*wstar(i)
        ! check for lcl
        sconv=.false.
        if (qtup(i,1).ge.qupsat(1)) then ! LCL
          sconv=.true.
          zilcl=zz(i,1)
          klcl=1
        end if
        
        ! vertical advection ----------
        ! dry convection case
        do k=2,kl
          dzht=dz_hl(i,k-1)
          zht=zz(i,k)
          ! update detrainment rate
          dtr=ee+0.05/max(zidry-zht,1.)
          ! update thermodynamics of plume
          tlup(i,k)=(tlup(i,k-1)+dzht*ee*thetal(i,k))/(1.+dzht*ee)
          qtup(i,k)=(qtup(i,k-1)+dzht*ee*qtot(i,k)  )/(1.+dzht*ee)          
          ! diagnose thermodynamic variables assuming no condensation         
          thup(i,k)=tlup(i,k)                                        ! theta,up
          tvup(i,k)=thup(i,k)+theta(i,k)*0.61*qtup(i,k)              ! thetav,up
          ttup(i,k)=thup(i,k)/sigkap(k)                              ! temp,up
          call getqsat(1,newsat(k),ttup(i,k),pres(k))                ! estimate of saturated mixing ratio in plume
          ! calculate buoyancy
          nn=grav*(tvup(i,k)-thetav(i,k))/thetav(i,k)
          ! update updraft velocitiy
          w2up(k)=(w2up(k-1)+2.*dzht*b2*nn)/(1.+2.*dzht*b1*ee)
          ! update mass flux
          mflx(i,k)=mflx(i,k-1)/(1.+dzht*(dtr-ee))
          ! check for lcl
          if (.not.sconv) then
            if (qtup(i,k).ge.qupsat(k)) then
              xp=dzht*max(min((qtup(i,k-1)-qupsat(k-1))/(qupsat(k)-qupsat(k-1)-qtup(i,k)+qtup(i,k-1)),1.),0.)
              sconv=.true.
              zilcl=xp+zz(i,k-1)
              klcl=k
            end if
          end if
          ! test if maximum plume height is reached
          if (w2up(k).le.0.) then
            if (k.gt.2) then ! use quadratic interpolation
              bs=1./(zz(i,k)-zz(i,k-2))*                                              &
                  ((zz(i,k-1)-zz(i,k-2))*(w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))      &
                  +(zz(i,k)-zz(i,k-1))*(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))
              as=1./(zz(i,k)-zz(i,k-2))*                           &
                  ((w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))         &
                  -(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))            
              cs=w2up(k-1)
              xp=(-bs-sqrt(bs*bs-4.*as*cs))/(2.*as)
              zidry=xp+zz(i,k-1)
            else ! use linear interpolation
              xp=w2up(k-1)/(w2up(k-1)-w2up(k))
              zidry=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
            end if
            exit
          end if
        end do

        ! shallow convection case
        do k=klcl,kl
          dzht=dz_hl(i,k-1)
          zht=zz(i,k)
          ! update detrainment rate
          dtr=0.9*ee+0.006/pi*(atan((zht-(zilcl+(zi(i)-zilcl)/1.5))/((zi(i)-zilcl)/8.))+0.5*pi)
          ! update thermodynamics of plume
          tlup(i,k)=(tlup(i,k-1)+dzht*ee*thetal(i,k))/(1.+dzht*ee)
          qtup(i,k)=(qtup(i,k-1)+dzht*ee*qtot(i,k)  )/(1.+dzht*ee)          
          ! diagnose thermodynamic variables assuming condensation
          qvup(i,k)=min(qtup(i,k),qupsat(k))                         ! qv,up
          qlup=qtup(i,k)-qvup(i,k)                                   ! ql,up
          thup(i,k)=tlup(i,k)+lv*sigkap(k)*qlup/cp                   ! theta,up
          tvup(i,k)=thup(i,k)+theta(i,k)*(1.61*qvup(i,k)-qtup(i,k))  ! thetav,up
          ttup(i,k)=thup(i,k)/sigkap(k)                              ! temp,up
          call getqsat(1,newsat(k),ttup(i,k),pres(k))                ! estimate of saturated mixing ratio in plume
          ! calculate buoyancy
          nn=grav*(tvup(i,k)-thetav(i,k))/thetav(i,k)
          ! update updraft velocitiy
          w2up(k)=(w2up(k-1)+2.*dzht*b2*nn)/(1.+2.*dzht*b1*ee)
          ! update mass flux
          mflx(i,k)=mflx(i,k-1)*(1.-dzht*(dtr-ee))
          ! test if maximum plume height is reached
          if (w2up(k).le.0.) then
            if (k.gt.2) then ! use quadratic interpolation
              bs=1./(zz(i,k)-zz(i,k-2))*                                              &
                  ((zz(i,k-1)-zz(i,k-2))*(w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))      &
                  +(zz(i,k)-zz(i,k-1))*(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))
              as=1./(zz(i,k)-zz(i,k-2))*                           &
                  ((w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))         &
                  -(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))            
              cs=w2up(k-1)
              xp=(-bs-sqrt(bs*bs-4.*as*cs))/(2.*as)
              zi(i)=xp+zz(i,k-1)
            else ! use linear interpolation
              xp=w2up(k-1)/(w2up(k-1)-w2up(k))
              zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
            end if
            exit
          end if
        end do
        
        if (.not.sconv) zi(i)=zidry

        ! soft updating of saturated mixing ratio in plume
        qupsat=alpha*newsat+(1.-alpha)*qupsat  
        ! update boundary layer height
        zi(i)=alpha*zi(i)+(1.-alpha)*ziold(i)
        if (abs(zi(i)-ziold(i)).lt.0.1) exit
        ziold(i)=zi(i)
       
        ! update surface boundary conditions
        wstar(i)=(grav*zi(i)*wtv0(i)/thetav(i,1))**(1./3.)
        tke(i,1)=ustar(i)*ustar(i)/sqrt(cm)+0.45*wstar(i)*wstar(i)
        tke(i,1)=max(tke(i,1),mintke)
        tke(i,1)=min(tke(i,1),65.)
      end do
      ! update explicit counter gradient terms
      gamtl(i,:)=mflx(i,:)*(tlup(i,:)-thetal(i,:))
      gamqt(i,:)=mflx(i,:)*(qtup(i,:)-qtot(i,:))
      gamtl(i,:)=max(min(gamtl(i,:),5.E-1),-5.E-1)
      gamqt(i,:)=max(min(gamqt(i,:),5.E-4),-5.E-4)
      !gamql(i,:)=0.
      gamqv(i,:)=gamqt(i,:) !-gamql(i,:)
      gamth(i,:)=gamtl(i,:)+lv*sigkap(:)*(gamqt(i,:)-gamqv(i,:))/cp
      gamtv(i,:)=gamth(i,:)+theta(i,:)*(1.61*gamqv(i,:)-gamqt(i,:))
    else                   ! stable
      !wpv_flux is calculated at half levels
      wpv_flux(1)=-kmo(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1) !+gamt_hl(i,k)
      if (wpv_flux(1).le.0.) then
        zi(i)=zz(i,1)
      else
        do k=2,kl-1
          wpv_flux(k)=-kmo(i,k)*(thetav(i,k+1)-thetav(i,k))/dz_hl(i,k) !+gamt_hl(i,k)
          if (wpv_flux(k).le.0.05*wpv_flux(1)) then
            xp=(0.05*wpv_flux(1)-wpv_flux(k-1))/(wpv_flux(k)-wpv_flux(k-1))
            xp=min(max(xp,0.),1.)
            if (xp.lt.0.5) then
              xp=xp/0.5
              zi(i)=(1.-xp)*0.5*(zz(i,k-1)+zz(i,k))+xp*zz(i,k)
            else
              xp=(xp-0.5)/0.5
              zi(i)=(1.-xp)*zz(i,k)+xp*0.5*(zz(i,k)+zz(i,k+1))
            end if
            exit
          end if
        end do
      end if
    end if
  end do
end if

! calculate wstar and phim (from TAPM)
wstar=0.
where (wtv0.gt.0.)
  wstar=(grav*zi*wtv0/thetav(:,1))**(1./3.)
end where
z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*ustar*ustar*ustar)
z_on_l=min(z_on_l,10.)
where (z_on_l.lt.0.)
  phim=(1.-16.*z_on_l)**(-0.25)
elsewhere !(z_on_l.le.0.4)
  phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
!elsewhere
!  phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar (2007)
end where
tke(1:ifull,1)=ustar*ustar/sqrt(cm)+0.45*wstar*wstar
eps(1:ifull,1)=ustar*ustar*ustar*(phim-z_on_l)/(vkar*zz(:,1))
tke(1:ifull,1)=max(tke(1:ifull,1),mintke)
tke(1:ifull,1)=min(tke(1:ifull,1),65.)
aa(:,2)=cm34*(tke(1:ifull,1)**1.5)/5.
eps(1:ifull,1)=min(eps(1:ifull,1),aa(:,2))
aa(:,2)=max(aa(:,2)*5./500.,mineps)
eps(1:ifull,1)=max(eps(1:ifull,1),aa(:,2))


! Calculate shear term
pps(:,2:kl)=shear(:,2:kl)


! loop to update kmo
do jcount=1,icm3

  ! Update diffusion coeffs
  if (jcount.gt.1) then
    km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-10)
    call updatekmo(kmo,km,zz) ! interpolate diffusion coeffs to half levels
  end if

  ! Update transport term
  do k=2,kl-1
    ppt(:,k)=0.5*((kmo(:,k))*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k) &
                 -(kmo(:,k-1))*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1))/dz_fl(:,k)
  end do
  ppt(:,kl)=0.5*(-(kmo(:,kl-1))*(tke(1:ifull,kl)-tke(1:ifull,kl-1))/dz_hl(:,kl-1))/dz_fl(:,kl)
  ppt=ppt/km(:,2:kl)

  ! Calculate buoyancy term
  select case(buoymeth)
    case(0) ! Hurley 2007 (dry-PBL with counter gradient term)
      do k=2,kl
        ppb(:,k)=-grav*0.5*(thetav(:,k+1)-thetav(:,k-1))/(thetav(:,k)*dz_fl(:,k))
        ppb(:,k)=ppb(:,k)+grav*gamtv(:,k)/(km(:,k)*thetav(:,k))
      end do
     case(3) ! saturated conditions from Durran and Klemp JAS 1982 (see also WRF)
      do k=2,kl
      ! saturated
        bb(:,k)=-grav*0.5/dz_fl(:,k)*(gg(:,k)*((theta(:,k+1)-theta(:,k-1))/theta(:,k)              &
                   +lv*(qsat(:,k+1)-qsat(:,k-1))/(cp*temp(:,k)))-qsat(:,k+1)-qlg(:,k+1)-qfg(:,k+1) &
                   +qsat(:,k-1)+qlg(:,k-1)+qfg(:,k-1))
        bb(:,k)=bb(:,k)+grav/km(:,k)*(gg(:,k)*(gamth(:,k)/theta(:,k)+lv*gamqv(:,k)/(cp*temp(:,k)))-gamqt(:,k))
        ! unsaturated
        cc(:,k)=-grav*0.5*(thetav(:,k+1)-thetav(:,k-1))/(thetav(:,k)*dz_fl(:,k))
        cc(:,k)=cc(:,k)+grav*gamtv(:,k)/(km(:,k)*thetav(:,k))
      end do
      ppb=(1.-cfrac(:,2:kl))*cc(:,2:kl)+cfrac(:,2:kl)*bb(:,2:kl) ! cloud fraction weighted (e.g., Smith 1990)
    case(2) ! Smith (1990)
      do k=2,kl
        templ(:,k)=temp(:,k)-(lv/cp*qlg(:,k)+ls/cp*qfg(:,k))
        jj(:)=lv+lf*qfg(:,k)/max(qlg(:,k)+qfg(:,k),1.e-12)  ! L
        dqsdt(:)=epsl*jj*qsat(:,k)/(rv*temp(:,k)*temp(:,k))
        hh(:)=cfrac(:,k)*(jj/cp/temp(:,k)-delta/(1.-epsl)/(1.+delta*qg(:,k)-qlg(:,k)-qfg(:,k)))/(1.+jj/cp*dqsdt) ! betac
        ff(:,k)=1./temp(:,k)-dqsdt*hh(:)                         ! betatt
        gg(:,k)=delta/(1.+delta*qg(:,k)-qlg(:,k)-qfg(:,k))+hh(:) ! betaqt
      end do    
      do k=2,kl
        ppb(:,k)=-grav*ff(:,k)*0.5*((templ(:,k+1)-templ(:,k-1))/dz_fl(:,k)+grav/cp) &
                 -grav*gg(:,k)*0.5*(qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-qg(:,k-1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
      end do
    case DEFAULT
      write(6,*) "ERROR: Unsupported buoyancy option ",buoymeth
      stop
  end select

  ! update tke-eps non-linear terms
  ! implicit approach (split form - helps with numerical stability)
  tkenew(:,2:kl)=tke(1:ifull,2:kl) ! 1st guess
  do icount=1,icm2
    tkenew(:,2:kl)=max(tkenew(:,2:kl),mintke)
    tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)
    aa(:,2:kl)=dt*ce2/tkenew(:,2:kl)
    bb(:,2:kl)=1.-dt*ce1*beta2*km(:,2:kl)*(pps+max(ppb(:,2:kl),0.)+max(ppt,0.))/tkenew(:,2:kl)
    cc(:,2:kl)=-epssav(1:ifull,2:kl)-dt*ce1*(1.-beta2)*cm*tkenew(:,2:kl)*(pps+max(ppb(:,2:kl),0.)+max(ppt,0.))
    gg(:,2:kl)=sqrt(bb(:,2:kl)*bb(:,2:kl)-4.*aa(:,2:kl)*cc(:,2:kl))
    epsnew(:,2:kl)=(gg(:,2:kl)-bb(:,2:kl))/(2.*aa(:,2:kl))
    dd(:,2:kl)=-tkenew(:,2:kl)+tkesav(1:ifull,2:kl)+dt*(((1.-beta1)*cm*tkenew(:,2:kl)*tkenew(:,2:kl)/epsnew(:,2:kl) &
               +beta1*km(:,2:kl))*(pps+ppb(:,2:kl))-epsnew(:,2:kl)) ! error function
    ff(:,2:kl)=-1.+dt*2.*(1.-beta1)*cm*tkenew(:,2:kl)/epsnew(:,2:kl)*(pps+ppb(:,2:kl))                                &
               -dt*(1.+(1.-beta1)*cm*tkenew(:,2:kl)*tkenew(:,2:kl)*(pps+ppb(:,2:kl))/(epsnew(:,2:kl)*epsnew(:,2:kl))) &
                  *(epsnew(:,2:kl)/tkenew(:,2:kl)-0.5*(1.-bb(:,2:kl))/(aa(:,2:kl)*tkenew(:,2:kl))                      &
                   +(0.5*bb(:,2:kl)*(1.-bb(:,2:kl))/tkenew(:,2:kl)+aa(:,2:kl)*cc(:,2:kl)/tkenew(:,2:kl)                &
                   +aa(:,2:kl)*dt*ce1*(1.-beta2)*cm*(pps+max(ppb(:,2:kl),0.)+max(ppt,0.)))/(aa(:,2:kl)*gg(:,2:kl)))
    where (abs(ff(:,2:kl)).gt.tol) ! sectant method 
      tkenew(:,2:kl)=tkenew(:,2:kl)-alpha2*dd(:,2:kl)/ff(:,2:kl)
    end where
  end do
  
  tkenew(:,2:kl)=max(tkenew(:,2:kl),mintke)
  tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)
  aa(:,2:kl)=cm34*(tkenew(:,2:kl)**1.5)/5.
  epsnew(:,2:kl)=min(epsnew(:,2:kl),aa(:,2:kl))
  aa(:,2:kl)=max(aa(:,2:kl)*5./500.,mineps)
  epsnew(:,2:kl)=max(epsnew(:,2:kl),aa(:,2:kl))

  tke(1:ifull,2:kl)=tkenew(:,2:kl)
  eps(1:ifull,2:kl)=epsnew(:,2:kl)


  ! eps vertical mixing (done here as we skip level 1, instead of using trim)
  aa(:,2)=-dt*0.5*ce0*kmo(:,1)/(dz_fl(:,2)*dz_hl(:,1))
  cc(:,2)=-dt*0.5*ce0*kmo(:,2)/(dz_fl(:,2)*dz_hl(:,2))
  bb(:,2)=1.-aa(:,2)-cc(:,2)
  dd(:,2)=eps(1:ifull,2)-aa(:,2)*eps(1:ifull,1)
  do k=3,kl-1
    aa(:,k)=-dt*0.5*ce0*kmo(:,k-1)/(dz_fl(:,k)*dz_hl(:,k-1))
    cc(:,k)=-dt*0.5*ce0*kmo(:,k)/(dz_fl(:,k)*dz_hl(:,k))
    bb(:,k)=1.-aa(:,k)-cc(:,k)
    dd(:,k)=eps(1:ifull,k)
  end do
  aa(:,kl)=-dt*0.5*ce0*kmo(:,kl-1)/(dz_fl(:,kl)*dz_hl(:,kl-1))
  bb(:,kl)=1.-aa(:,kl)
  dd(:,kl)=eps(1:ifull,kl)
  call thomas(epsnew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

  ! TKE vertical mixing (done here as we skip level 1, instead of using trim)
  aa(:,2)=-dt*0.5*kmo(:,1)/(dz_fl(:,2)*dz_hl(:,1))
  cc(:,2)=-dt*0.5*kmo(:,2)/(dz_fl(:,2)*dz_hl(:,2))
  bb(:,2)=1.-aa(:,2)-cc(:,2)
  dd(:,2)=tke(1:ifull,2)-aa(:,2)*tke(1:ifull,1)
  do k=3,kl-1
    aa(:,k)=-dt*0.5*kmo(:,k-1)/(dz_fl(:,k)*dz_hl(:,k-1))
    cc(:,k)=-dt*0.5*kmo(:,k)/(dz_fl(:,k)*dz_hl(:,k))
    bb(:,k)=1.-aa(:,k)-cc(:,k)
    dd(:,k)=tke(1:ifull,k)
  end do
  aa(:,kl)=-dt*0.5*kmo(:,kl-1)/(dz_fl(:,kl)*dz_hl(:,kl-1))
  bb(:,kl)=1.-aa(:,kl)
  dd(:,kl)=tke(1:ifull,kl)
  call thomas(tkenew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

  ! Limits
  if (kmlimit.eq.1) then
    do k=2,kl
      where (zz(:,k).ge.0.55*zi.and.zz(:,k).le.0.95*zi.and.wstar.gt.0.5)
        tkenew(1:ifull,k)=max(tkenew(1:ifull,k),tkenew(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
        epsnew(1:ifull,k)=max(epsnew(1:ifull,k),epsnew(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
      end where
    end do
  end if
  tkenew(:,2:kl)=max(tkenew(:,2:kl),mintke)
  tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)
  aa(:,2:kl)=cm34*(tkenew(:,2:kl)**1.5)/5.
  epsnew(:,2:kl)=min(epsnew(:,2:kl),aa(:,2:kl))
  aa(:,2:kl)=max(aa(:,2:kl)*5./500.,mineps)
  epsnew(:,2:kl)=max(epsnew(:,2:kl),aa(:,2:kl))
  tke(1:ifull,2:kl)=tkenew(:,2:kl)
  eps(1:ifull,2:kl)=epsnew(:,2:kl)

end do

! Update theta and qg due to non-local term
if (mode.ne.1) then
  ! explicit
  call interphl(gamtl,gam_hl,zi,zz)
  thetal(:,1)=thetal(:,1)-dt*gam_hl(:,1)/dz_fl(:,1)
  do k=2,kl-1
    thetal(:,k)=thetal(:,k)-dt*(gam_hl(:,k)-gam_hl(:,k-1))/dz_fl(:,k)
  end do
end if
if (mode.eq.0) then
  ! explicit
  call interphl(gamqt,gam_hl,zi,zz)
  qtot(:,1)=qtot(:,1)-dt*gam_hl(:,1)/dz_fl(:,1)
  do k=2,kl-1
    qtot(:,k)=qtot(:,k)-dt*(gam_hl(:,k)-gam_hl(:,k-1))/dz_fl(:,k)
  end do
  qg=max(qtot-qlg-qfg,1.E-6)
end if
do k=1,kl
  theta(:,k)=thetal(:,k)+lv*sigkap(k)*qlg(:,k)/cp
end do

tkesav=tke(1:ifull,:) ! Not needed, but for consistancy when not using CCAM
epssav=eps(1:ifull,:) ! Not needed, but for consistancy when not using CCAM

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver

subroutine thomas(out,aa,bbi,cc,ddi)

implicit none

real, dimension(ifull,3:kl), intent(in) :: aa
real, dimension(ifull,2:kl), intent(in) :: bbi,ddi
real, dimension(ifull,2:kl-1), intent(in) :: cc
real, dimension(ifull,2:kl), intent(out) :: out
real, dimension(ifull,2:kl) :: bb,dd
real, dimension(ifull) :: n
integer k

bb=bbi
dd=ddi

do k=3,kl
  n=aa(:,k)/bb(:,k-1)
  bb(:,k)=bb(:,k)-n*cc(:,k-1)
  dd(:,k)=dd(:,k)-n*dd(:,k-1)
end do
out(:,kl)=dd(:,kl)/bb(:,kl)
do k=kl-1,2,-1
  out(:,k)=(dd(:,k)-cc(:,k)*out(:,k+1))/bb(:,k)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsat(ilen,qsat,temp,ps)

implicit none

integer, intent(in) :: ilen
real, dimension(ilen), intent(in) :: temp,ps
real, dimension(ilen), intent(out) :: qsat
real, dimension(0:220) :: table
real, dimension(ilen) :: esatf,tdiff,rx
integer, dimension(ilen) :: ix

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
table(219:220)=(/ 29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/(ps-esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

subroutine updatekmo(kmo,km,zz)

implicit none

integer k
real, dimension(ifull,kl), intent(in) :: km,zz
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull) :: zhl
real xp

zhl=0.5*sum(zz(:,1:2),2)
kmo(:,1)=km(:,2)+(zhl-zz(:,2))/(zz(:,3)-zz(:,1))*                    &
           ((zz(:,2)-zz(:,1))*(km(:,3)-km(:,2))/(zz(:,3)-zz(:,2))    &
           +(zz(:,3)-zz(:,2))*(km(:,2)-km(:,1))/(zz(:,2)-zz(:,1)))   &
         +(zhl-zz(:,2))**2/(zz(:,3)-zz(:,1))*                        &
           ((km(:,3)-km(:,2))/(zz(:,3)-zz(:,2))                      &
           -(km(:,2)-km(:,1))/(zz(:,2)-zz(:,1)))
do k=2,kl-1
  zhl=0.5*sum(zz(:,k:k+1),2)
  kmo(:,k)=km(:,k)+(zhl-zz(:,k))/(zz(:,k+1)-zz(:,k-1))*                      &
             ((zz(:,k)-zz(:,k-1))*(km(:,k+1)-km(:,k))/(zz(:,k+1)-zz(:,k))    &
             +(zz(:,k+1)-zz(:,k))*(km(:,k)-km(:,k-1))/(zz(:,k)-zz(:,k-1)))   &
           +(zhl-zz(:,k))**2/(zz(:,k+1)-zz(:,k-1))*                          &
             ((km(:,k+1)-km(:,k))/(zz(:,k+1)-zz(:,k))                        &
             -(km(:,k)-km(:,k-1))/(zz(:,k)-zz(:,k-1)))
end do
! These terms are never used
kmo(:,kl)=2.*kmo(:,kl-1)-kmo(:,kl-2)
kmo=max(kmo,1.E-10)

return
end subroutine updatekmo

subroutine interphl(gamfl,gamhl,zi,zz)

implicit none

integer k
real, dimension(ifull) :: zht,xp
real, dimension(ifull,kl), intent(in) :: gamfl,zz
real, dimension(ifull,kl-1), intent(out) :: gamhl
real, dimension(ifull), intent(in) :: zi
real, dimension(ifull) :: zhl

gamhl=0.

zhl=0.5*sum(zz(:,1:2),2)
where(zhl.gt.zi)
  gamhl(:,1)=0.
elsewhere (zz(:,2).gt.zi)
  gamhl(:,1)=gamfl(:,1) ! assume const in this case
elsewhere
  gamhl(:,1)=gamfl(:,2)+(zhl-zz(:,2))/(zz(:,3)-zz(:,1))*                     &
             ((zz(:,2)-zz(:,1))*(gamfl(:,3)-gamfl(:,2))/(zz(:,3)-zz(:,2))    &
             +(zz(:,3)-zz(:,2))*(gamfl(:,2)-gamfl(:,1))/(zz(:,2)-zz(:,1)))   &
           +(zhl-zz(:,2))**2/(zz(:,3)-zz(:,1))*                              &
             ((gamfl(:,3)-gamfl(:,2))/(zz(:,3)-zz(:,2))                      &
             -(gamfl(:,2)-gamfl(:,1))/(zz(:,2)-zz(:,1)))
end where

do k=2,kl-1
  zhl=0.5*sum(zz(:,k:k+1),2)
  where(zhl.gt.zi)
    gamhl(:,k)=0.
  elsewhere (zz(:,k+1).gt.zi)
    gamhl(:,k)=gamfl(:,k)+(zhl-zz(:,k))/(zi-zz(:,k-1))*                         &
                 ((zz(:,k)-zz(:,k-1))*(0.-gamfl(:,k))/(zi-zz(:,k))              &
                 +(zi-zz(:,k))*(gamfl(:,k)-gamfl(:,k-1))/(zz(:,k)-zz(:,k-1)))   &
               +(zhl-zz(:,k))**2/(zi-zz(:,k-1))*                                &
                 ((0.-gamfl(:,k))/(zi-zz(:,k-1))                                &
                 -(gamfl(:,k)-gamfl(:,k-1))/(zz(:,k)-zz(:,k-1)))
  elsewhere
    gamhl(:,k)=gamfl(:,k)+(zhl-zz(:,k))/(zz(:,k+1)-zz(:,k-1))*                         &
                 ((zz(:,k)-zz(:,k-1))*(gamfl(:,k+1)-gamfl(:,k))/(zz(:,k+1)-zz(:,k))    &
                 +(zz(:,k+1)-zz(:,k))*(gamfl(:,k)-gamfl(:,k-1))/(zz(:,k)-zz(:,k-1)))   &
               +(zhl-zz(:,k))**2/(zz(:,k+1)-zz(:,k-1))*                                &
                 ((gamfl(:,k+1)-gamfl(:,k))/(zz(:,k+1)-zz(:,k))                        &
                 -(gamfl(:,k)-gamfl(:,k-1))/(zz(:,k)-zz(:,k-1)))
  end where
end do

return
end subroutine interphl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE

subroutine tkeend(diag)

implicit none

integer, intent(in) :: diag

if (diag.gt.0) write(6,*) "Terminate TKE scheme"

deallocate(tke,eps)
deallocate(tkesav,epssav)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps
