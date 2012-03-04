
! This module calculates the turblent kinetic energy and mixing for the boundary layer based on Hurley 2009
! (eddy dissipation) and Angevine et al 2010 (mass flux).  Specifically, this version is adapted for CCAM
! with an implicit numerical formulation (i.e., a forward-backward scheme), as well as modifications for
! buoyancy within clouds.

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
real, parameter :: b1      = 2.     ! Soares et al (2004) 1., Siebesma et al (2003) 2.
real, parameter :: b2      = 1./3.  ! Soares et al (2004) 2., Siebesma et al (2003) 1./3.
real, parameter :: be      = 1.0    ! Hurley (2007) 1., Soares et al (2004) 0.3, Angevine (2005) 7.
real, parameter :: cm      = 0.03   ! Hurley (2007) 0.09 (0.03 recommended by Ashok Luhar)
real, parameter :: ce0     = 0.69
real, parameter :: ce1     = 1.46
real, parameter :: ce2     = 1.83
real, parameter :: ce3     = 0.45
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
integer, parameter :: icm1 = 40      ! iterations for calculating pblh
integer, parameter :: icm2 = 10      ! iterations for calculating tke-eps
integer, parameter :: icm3 = 5       ! iterations for calculating qupsat
real, parameter :: dz_ref = 500.     ! dz reference for kmlimit=1
real, parameter :: alpha  = 0.3      ! weight for updating pblh
real, parameter :: alpha3 = 0.7      ! weight for updating qupsat
real, parameter :: beta3  = 0.5      ! weight for explicit tke(t) (beta3=0.) or implicit tke(t+1) (beta3=1.) non-linear terms
real, parameter :: beta4  = 0.5      ! weight for explicit eps(t) (beta4=0.) or implicit eps(t+1) (beta4=1.) non-linear terms
real, parameter :: tol    = 2.E-4    ! tolerance for tke-eps non-linear solver
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

subroutine tkemix(kmo,theta,qg,qlg,qfg,cfrac,zi,wt0,wq0,ps,ustar,zz,zzh,sig,sigkap,dt,mode,diag)

implicit none

integer, intent(in) :: diag,mode
integer k,i,klcl,icount,jcount
real, intent(in) :: dt
real, dimension(ifull,kl), intent(inout) :: theta,qg
real, dimension(ifull,kl), intent(in) :: zz,cfrac,qlg,qfg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl-1), intent(in) :: zzh
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar
real, dimension(kl), intent(in) :: sigkap,sig
real, dimension(ifull,kl) :: km,ff,gg,templ,thetav,thetal,temp,qtot
real, dimension(ifull,kl) :: gamtt,gamtl,gamtv,gamth,gamqt,gamqv,gamhl
real, dimension(ifull,kl) :: tkenew,epsnew,qsat,bb,cc,dd
real, dimension(ifull,kl) :: mflx,tlup,qtup
real, dimension(ifull,2:kl) :: aa,pps,ppt,ppb
real, dimension(ifull,kl) :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,kl-1) :: gam_hl,zzm
real, dimension(ifull) :: wstar,z_on_l,phim,hh,jj,dqsdt
real, dimension(ifull) :: dum,wtv0
real, dimension(kl-1) :: wpv_flux
real, dimension(kl) :: w2up,ttup,tvup,thup,qvup
real, dimension(kl) :: qupsat,newsat,pres
real xp,ee,dtr,nn,as,bs,cs,cm34,qlup
real zht,dzht,zidry,zilcl,oldqupsat
real ziold
logical sconv

cm34=cm**0.75

if (diag.gt.0) write(6,*) "Update PBL mixing with TKE-eps turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

! impose limits after host advection
tke(1:ifull,:)=max(tke(1:ifull,:),mintke)
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
qtot=qg+qlg+qfg
thetav=theta*(1.+0.61*qg-qlg-qfg)
do k=1,kl
  thetal(:,k)=theta(:,k)-sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
end do
wtv0=wt0+theta(:,1)*0.61*wq0
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
  zzm(:,k)=0.5*(zz(:,k+1)+zz(:,k))
end do

! Calculate first approximation to diffusion coeffs
km=cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:)
call updatekmo(kmo,km,zz,zzm) ! interpolate diffusion coeffs to internal half levels (zzm)

! Calculate transport term on full levels
do k=2,kl-1
  ppt(:,k)=(kmo(:,k)*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k) &
           -kmo(:,k-1)*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1))/dz_fl(:,k)
end do
ppt(:,kl)=0.
ppt=ppt/km(:,2:kl)

! Calculate shear term on full levels (see hordifg.f for calculation of horizontal shear)
pps(:,2:kl-1)=shear(:,2:kl-1)
pps(:,kl)=0.



! Calculate non-local mass-flux terms for theta_l and qtot
mflx=0.
gamtl=0.
gamtv=0.
gamth=0.
gamqt=0.
gamqv=0.
tlup=thetal
qtup=qtot
wstar=(grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)
if (mode.ne.1) then ! mass flux
  do i=1,ifull
    if (wtv0(i).gt.0.) then ! unstable
      pres(:)=ps(i)*sig(:) ! pressure
      ! initial guess for plume state
      tvup=thetav(i,:)
      thup=theta(i,:)
      ttup=temp(i,:)
      qvup=qg(i,:)
      qupsat(:)=qsat(i,:)
      newsat(:)=qsat(i,:)
      zidry=zi(i)
      ziold=zi(i)
      do icount=1,icm1
        klcl=kl+1
        zilcl=zidry
        mflx(i,:)=0.
        ! entrainment rate from Angevine et al 2010
        ee=2./max(100.,zidry)
        ! first level -----------------
        ! initial thermodynamic state
        tlup(i,1)=thetal(i,1)+be*wt0(i)/sqrt(tke(i,1)) ! Hurley 2007
        qtup(i,1)=qtot(i,1)  +be*wq0(i)/sqrt(tke(i,1)) ! Hurley 2007
        ! diagnose thermodynamic variables assuming no condensation
        thup(1)=tlup(i,1)                                ! theta,up
        tvup(1)=thup(1)+theta(i,1)*0.61*qtup(i,1)        ! thetav,up
        ttup(1)=thup(1)/sigkap(1)                        ! temp,up
        call getqsat(1,qupsat(1),ttup(1),pres(1))        ! estimate of saturated mixing ratio in plume
        ! update updraft velocity and mass flux
        !nn=grav*wtv0(i)/thetav(i,1)
        !w2up(1)=2.*zz(i,1)*b2*nn/(1.+2.*zz(i,1)*b1*ee)  ! Hurley 2007
        !mflx(i,1)=0.1*sqrt(w2up(1))                     ! Hurley 2007
        w2up(1)=0.25*wstar(i)**2                         ! Angevine et al 2010
        mflx(i,1)=0.03*wstar(i)                          ! Angevine et al 2010
        ! check for lcl
        sconv=.false.
        if (qtup(i,1).ge.qupsat(1)) then
          sconv=.true.
          zilcl=zz(i,1)
          klcl=2
        end if
        
        ! dry convection case
        do k=2,kl-1
          dzht=dz_hl(i,k-1)
          zht=zz(i,k)
          ! update detrainment rate from Angevine et al 2010
          dtr=ee+0.05/max(zidry-zht,0.001)
          ! update thermodynamics of plume
          ! (use upwind as centred scheme requires vertical spacing less than 250m)
          tlup(i,k)=(tlup(i,k-1)+dzht*ee*thetal(i,k))/(1.+dzht*ee)
          qtup(i,k)=(qtup(i,k-1)+dzht*ee*qtot(i,k)  )/(1.+dzht*ee)
          ! diagnose thermodynamic variables assuming no condensation         
          thup(k)=tlup(i,k)                                      ! theta,up
          ttup(k)=thup(k)/sigkap(k)                              ! temp,up
          call getqsat(1,qupsat(k),ttup(k),pres(k))              ! estimate of saturated mixing ratio in plume
          tvup(k)=thup(k)+theta(i,k)*0.61*qtup(i,k)              ! thetav,up
          ! calculate buoyancy
          nn=grav*(tvup(k)-thetav(i,k))/thetav(i,k)
          ! update updraft velocity
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
              bs=1./(zz(i,k)-zz(i,k-2))*                                            &
                  ((zz(i,k-1)-zz(i,k-2))*(w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))    &
                  +(zz(i,k)-zz(i,k-1))*(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))
              as=1./(zz(i,k)-zz(i,k-2))*                        &
                  ((w2up(k)-w2up(k-1))/(zz(i,k)-zz(i,k-1))      &
                  -(w2up(k-1)-w2up(k-2))/(zz(i,k-1)-zz(i,k-2)))            
              cs=w2up(k-1)
              if (as.eq.0.) then
                xp=w2up(k-1)*dzht/(w2up(k-1)-w2up(k))
              else
                xp=(-bs-sqrt(bs*bs-4.*as*cs))/(2.*as)
              end if
            else ! use linear interpolation
              xp=w2up(k-1)*dzht/(w2up(k-1)-w2up(k))
            end if
            zidry=xp+zz(i,k-1)
            mflx(i,k)=0.
            if (sconv) then
              if (zidry.lt.zilcl) then
                sconv=.false.
                klcl=kl+1
              end if
            end if
            exit
          end if
        end do

        ! shallow convection case
        do k=klcl,kl-1
          qupsat(k)=newsat(k)
          dzht=dz_hl(i,k-1)
          zht=zz(i,k)
          ! update detrainment rate for cloudly air from Angevine et al 2010
          xp=max(zi(i)-zilcl,0.1)
          xp=8.*min(zht-zilcl,xp)/xp-16./3.
          dtr=max(0.9*ee+0.006/pi*(atan(xp)+0.5*pi),ee)
          ! update thermodynamics of plume
          tlup(i,k)=(tlup(i,k-1)+dzht*ee*thetal(i,k))/(1.+dzht*ee)
          qtup(i,k)=(qtup(i,k-1)+dzht*ee*qtot(i,k)  )/(1.+dzht*ee)
          ! diagnose thermodynamic variables assuming condensation
          do jcount=1,icm3
            ! qvup will always be saturated here.  The limit is just to help numerical solution
            qvup(k)=min(qtup(i,k),qupsat(k))                     ! qv,up
            qlup=qtup(i,k)-qvup(k)                               ! ql,up
            thup(k)=tlup(i,k)+lv*sigkap(k)*qlup/cp               ! theta,up
            ttup(k)=thup(k)/sigkap(k)                            ! temp,up
            call getqsat(1,newsat(k),ttup(k),pres(k))            ! estimate of saturated mixing ratio in plume
            qupsat(k)=alpha3*newsat(k)+(1.-alpha3)*qupsat(k)
          end do
          tvup(k)=thup(k)+theta(i,k)*(1.61*qvup(k)-qtup(i,k))    ! thetav,up
          ! calculate buoyancy
          nn=grav*(tvup(k)-thetav(i,k))/thetav(i,k)
          ! update updraft velocity
          w2up(k)=(w2up(k-1)+2.*dzht*b2*nn)/(1.+2.*dzht*b1*ee)
          ! update mass flux
          mflx(i,k)=mflx(i,k-1)/(1.+dzht*(dtr-ee))
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
              if (as.eq.0.) then
                xp=w2up(k-1)*dzht/(w2up(k-1)-w2up(k))
              else
                xp=(-bs-sqrt(bs*bs-4.*as*cs))/(2.*as)
              end if
            else ! use linear interpolation
              xp=w2up(k-1)*dzht/(w2up(k-1)-w2up(k))
            end if
            zi(i)=xp+zz(i,k-1)
            mflx(i,k)=0.
            exit
          end if
        end do
        
        if (.not.sconv) zi(i)=zidry

        ! update boundary layer height
        zi(i)=alpha*zi(i)+(1.-alpha)*ziold
        if (abs(zi(i)-ziold).lt.1.) exit
        ziold=zi(i)
       
        ! update surface boundary conditions
        wstar(i)=(grav*zi(i)*wtv0(i)/thetav(i,1))**(1./3.)
        tke(i,1)=ustar(i)*ustar(i)/sqrt(cm)+ce3*wstar(i)*wstar(i)
        tke(i,1)=max(tke(i,1),mintke)
      end do
      ! update explicit counter gradient terms
      gamtl(i,:)=mflx(i,:)*(tlup(i,:)-thetal(i,:))
      gamqt(i,:)=mflx(i,:)*(qtup(i,:)-qtot(i,:))
      !gamql(i,:)=0.
      gamqv(i,:)=gamqt(i,:) !-gamql(i,:)
      gamth(i,:)=gamtl(i,:) !+lv*sigkap(:)*gamql(i,:)/cp
      gamtv(i,:)=gamth(i,:)+theta(i,:)*0.61*gamqv(i,:) !-theta(i,:)*gamql(i,:))
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

! calculate tke and eps at 1st level
z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-3))
z_on_l=min(z_on_l,10.)
where (z_on_l.lt.0.)
  phim=(1.-16.*z_on_l)**(-0.25)
elsewhere !(z_on_l.le.0.4)
  phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
!elsewhere
!  phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar (2007)
end where
tke(1:ifull,1)=ustar*ustar/sqrt(cm)+ce3*wstar*wstar
eps(1:ifull,1)=ustar*ustar*ustar*(phim-z_on_l)/(vkar*zz(:,1))
tke(1:ifull,1)=max(tke(1:ifull,1),mintke)
aa(:,2)=cm34*(tke(1:ifull,1)**1.5)/5.
eps(1:ifull,1)=min(eps(1:ifull,1),aa(:,2))
aa(:,2)=max(aa(:,2)*5./500.,mineps)
eps(1:ifull,1)=max(eps(1:ifull,1),aa(:,2))



! Calculate buoyancy term
select case(buoymeth)
  case(0) ! Hurley 2007 (dry-PBL with counter gradient term)
    do k=2,kl-1
      ppb(:,k)=-grav*0.5*(thetav(:,k+1)-thetav(:,k-1))/(thetav(:,k)*dz_fl(:,k))
      ppb(:,k)=ppb(:,k)+grav*gamtv(:,k)/(km(:,k)*thetav(:,k))
    end do
    ppb(:,kl)=0.
  case(3) ! saturated conditions from Durran and Klemp JAS 1982 (see also WRF)
    do k=2,kl-1
      ! saturated
      bb(:,k)=-grav*0.5/dz_fl(:,k)*(gg(:,k)*((theta(:,k+1)-theta(:,k-1))/theta(:,k)              &
                 +lv*(qsat(:,k+1)-qsat(:,k-1))/(cp*temp(:,k)))-qsat(:,k+1)-qlg(:,k+1)-qfg(:,k+1) &
                 +qsat(:,k-1)+qlg(:,k-1)+qfg(:,k-1))
      bb(:,k)=bb(:,k)+grav/km(:,k)*(gg(:,k)*(gamth(:,k)/theta(:,k)+lv*gamqv(:,k)/(cp*temp(:,k)))-gamqt(:,k))
      ! unsaturated
      cc(:,k)=-grav*0.5*(thetav(:,k+1)-thetav(:,k-1))/(thetav(:,k)*dz_fl(:,k))
      cc(:,k)=cc(:,k)+grav*gamtv(:,k)/(km(:,k)*thetav(:,k))
    end do
    bb(:,kl)=0.
    cc(:,kl)=0.
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
    do k=2,kl-1
      ppb(:,k)=-grav*ff(:,k)*0.5*((templ(:,k+1)-templ(:,k-1))/dz_fl(:,k)+grav/cp) &
               -grav*gg(:,k)*0.5*(qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-qg(:,k-1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
    end do
    ppb(:,kl)=0.
  case DEFAULT
    write(6,*) "ERROR: Unsupported buoyancy option ",buoymeth
    stop
end select



! update tke-eps non-linear terms
! False position method (robust with highly non-linear error function for small tke,
! whereas the Newton-Raphson method was found to be unstable)
aa(:,2:kl)=-200.  ! lower tke limit
cc(:,2:kl)=200.   ! upper tke limit
call eenonlin(bb(:,2:kl),aa(:,2:kl),epsnew(:,2:kl),km(:,2:kl),pps,ppb,ppt,dt)
call eenonlin(ff(:,2:kl),cc(:,2:kl),epsnew(:,2:kl),km(:,2:kl),pps,ppb,ppt,dt)
!tkenew(:,2:kl)=0.5*(aa(:,2:kl)+cc(:,2:kl))
tkenew(:,2:kl)=(ff(:,2:kl)*aa(:,2:kl)-bb(:,2:kl)*cc(:,2:kl))/(ff(:,2:kl)-bb(:,2:kl))
do icount=1,icm2
  call eenonlin(dd(:,2:kl),tkenew(:,2:kl),epsnew(:,2:kl),km(:,2:kl),pps,ppb,ppt,dt)
  where (bb(:,2:kl)*dd(:,2:kl).le.0.)
    cc(:,2:kl)=tkenew(:,2:kl)
    ff(:,2:kl)=dd(:,2:kl)
  elsewhere
    aa(:,2:kl)=tkenew(:,2:kl)
    bb(:,2:kl)=dd(:,2:kl)
  end where
  !tkenew(:,2:kl)=0.5*(aa(:,2:kl)+cc(:,2:kl))
  tkenew(:,2:kl)=(ff(:,2:kl)*aa(:,2:kl)-bb(:,2:kl)*cc(:,2:kl))/(ff(:,2:kl)-bb(:,2:kl))
end do

tkenew(:,2:kl)=max(tkenew(:,2:kl),mintke)
aa(:,2:kl)=cm34*(tkenew(:,2:kl)**1.5)/5.
epsnew(:,2:kl)=min(epsnew(:,2:kl),aa(:,2:kl))
aa(:,2:kl)=max(aa(:,2:kl)/100.,mineps)
epsnew(:,2:kl)=max(epsnew(:,2:kl),aa(:,2:kl))

tke(1:ifull,2:kl)=tkenew(:,2:kl)
eps(1:ifull,2:kl)=epsnew(:,2:kl)



! eps vertical mixing (done here as we skip level 1, instead of using trim)
aa(:,2)=-dt*ce0*kmo(:,1)/(dz_fl(:,2)*dz_hl(:,1))
cc(:,2)=-dt*ce0*kmo(:,2)/(dz_fl(:,2)*dz_hl(:,2))
bb(:,2)=1.-aa(:,2)-cc(:,2)
dd(:,2)=eps(1:ifull,2)-aa(:,2)*eps(1:ifull,1)
do k=3,kl-1
  aa(:,k)=-dt*ce0*kmo(:,k-1)/(dz_fl(:,k)*dz_hl(:,k-1))
  cc(:,k)=-dt*ce0*kmo(:,k)/(dz_fl(:,k)*dz_hl(:,k))
  bb(:,k)=1.-aa(:,k)-cc(:,k)
  dd(:,k)=eps(1:ifull,k)
end do
aa(:,kl)=-dt*ce0*kmo(:,kl-1)/(dz_fl(:,kl)*dz_hl(:,kl-1))
bb(:,kl)=1.-aa(:,kl)
dd(:,kl)=eps(1:ifull,kl)
call thomas(epsnew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

! TKE vertical mixing (done here as we skip level 1, instead of using trim)
aa(:,2)=-dt*kmo(:,1)/(dz_fl(:,2)*dz_hl(:,1))
cc(:,2)=-dt*kmo(:,2)/(dz_fl(:,2)*dz_hl(:,2))
bb(:,2)=1.-aa(:,2)-cc(:,2)
dd(:,2)=tke(1:ifull,2)-aa(:,2)*tke(1:ifull,1)
do k=3,kl-1
  aa(:,k)=-dt*kmo(:,k-1)/(dz_fl(:,k)*dz_hl(:,k-1))
  cc(:,k)=-dt*kmo(:,k)/(dz_fl(:,k)*dz_hl(:,k))
  bb(:,k)=1.-aa(:,k)-cc(:,k)
  dd(:,k)=tke(1:ifull,k)
end do
aa(:,kl)=-dt*kmo(:,kl-1)/(dz_fl(:,kl)*dz_hl(:,kl-1))
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
aa(:,2:kl)=cm34*(tkenew(:,2:kl)**1.5)/5.
epsnew(:,2:kl)=min(epsnew(:,2:kl),aa(:,2:kl))
aa(:,2:kl)=max(aa(:,2:kl)*5./500.,mineps)
epsnew(:,2:kl)=max(epsnew(:,2:kl),aa(:,2:kl))
tke(1:ifull,2:kl)=tkenew(:,2:kl)
eps(1:ifull,2:kl)=epsnew(:,2:kl)



! Update diffusion coeffs
km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-10)
call updatekmo(kmo,km,zz,zzm) ! internal mid point levels

! Update thetal and qtot due to non-local and diffusion terms
! (Here we use quadratic interpolation to determine the explicit
! counter gradient terms)
call updategam(gamhl,gamtl,zz,zzm,zi)
cc(:,1)=-dt*kmo(:,1)/dz_hl(:,1)/dz_fl(:,1)
bb(:,1)=1.+dt*kmo(:,1)/dz_hl(:,1)/dz_fl(:,1)
dd(:,1)=thetal(:,1)-dt*gamhl(:,1)/dz_fl(:,1)+dt*wt0/dz_fl(:,1)
do k=2,kl-2
  aa(:,k)=-dt*kmo(:,k-1)/dz_hl(:,k-1)/dz_fl(:,k)
  cc(:,k)=-dt*kmo(:,k)/dz_hl(:,k)/dz_fl(:,k)
  bb(:,k)=1.+dt*(kmo(:,k-1)/dz_hl(:,k-1)+kmo(:,k)/dz_hl(:,k))/dz_fl(:,k)
  dd(:,k)=thetal(:,k)+dt*(gamhl(:,k-1)-gamtl(:,k))/dz_fl(:,k)
end do
aa(:,kl-1)=-dt*kmo(:,kl-2)/dz_hl(:,kl-2)/dz_fl(:,kl-1)
bb(:,kl-1)=1.+dt*kmo(:,kl-2)/dz_hl(:,kl-2)/dz_fl(:,kl-1)
dd(:,kl-1)=thetal(:,kl-1)+dt*gamhl(:,kl-2)/dz_fl(:,kl-1)
call thomas(thetal(:,1:kl-1),aa(:,2:kl-1),bb(:,1:kl-1),cc(:,1:kl-2),dd(:,1:kl-1))

call updategam(gamhl,gamqt,zz,zzm,zi)
dd(:,1)=qtot(:,1)-dt*gamhl(:,1)/dz_fl(:,1)+dt*wq0/dz_fl(:,1)
do k=2,kl-2
  dd(:,k)=qtot(:,k)+dt*(gamhl(:,k-1)-gamhl(:,k))/dz_fl(:,k)
end do
dd(:,kl-1)=qtot(:,kl-1)+dt*gamhl(:,kl-2)/dz_fl(:,kl-1)
call thomas(qtot(:,1:kl-1),aa(:,2:kl-1),bb(:,1:kl-1),cc(:,1:kl-2),dd(:,1:kl-1))

qg=max(qtot-qlg-qfg,1.E-6)
do k=1,kl
  theta(:,k)=thetal(:,k)+sigkap(k)*(lv*qlg(:,k)+ls*qfg(:,k))/cp
end do



! Update diffusion coeffs
km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-10)
call updatekmo(kmo,km,zz,zzh) ! output sigmh levels

tkesav=tke(1:ifull,:) ! Not needed, but for consistancy when not using CCAM
epssav=eps(1:ifull,:) ! Not needed, but for consistancy when not using CCAM

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Error function for tke-eps non-linear terms

subroutine eenonlin(dd,tkenew,epsnew,km,pps,ppb,ppt,dt)

implicit none

real, intent(in) :: dt
real, dimension(ifull,2:kl), intent(in) :: tkenew,km,pps,ppb,ppt
real, dimension(ifull,2:kl), intent(out) :: dd,epsnew
real, dimension(ifull,2:kl) :: aa,bb,cc,gg
real, dimension(ifull,2:kl) :: efftke,effeps

efftke=beta4*max(tkenew,mintke)+(1.-beta4)*tke(1:ifull,2:kl)
aa=(1.-beta3)**2*dt*ce2/efftke
bb=1.-(1.-beta3)*dt*ce1*km*(pps+max(ppb,0.)+max(ppt,0.))/efftke &
   +2.*(1.-beta3)*beta3*dt*ce2*eps(1:ifull,2:kl)/efftke
cc=-eps(1:ifull,2:kl)-eps(1:ifull,2:kl)*beta3*dt*ce1*km*(pps+max(ppb,0.)+max(ppt,0.))/efftke  &
   +beta3*beta3*dt*ce2*eps(1:ifull,2:kl)*eps(1:ifull,2:kl)/efftke
gg=sqrt(max(bb*bb-4.*aa*cc,0.))
epsnew=0.5*(gg-bb)/aa
effeps=beta3*epsnew+(1.-beta3)*eps(1:ifull,2:kl)
dd=-tkenew+tke(1:ifull,2:kl)+dt*(km*(pps+ppb)-effeps) ! error function

return
end subroutine eenonlin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tri-diagonal solver (single point version)
!subroutine thomas1(out,aa,bbi,cc,ddi)
!
!implicit none
!
!real, dimension(2:kl), intent(in) :: aa
!real, dimension(kl), intent(in) :: bbi,ddi
!real, dimension(kl-1), intent(in) :: cc
!real, dimension(kl), intent(out) :: out
!real, dimension(kl) :: bb,dd
!real n
!integer k
!
!bb=bbi
!dd=ddi
!
!do k=2,kl
!  n=aa(k)/bb(k-1)
!  bb(k)=bb(k)-n*cc(k-1)
!  dd(k)=dd(k)-n*dd(k-1)
!end do
!out(kl)=dd(kl)/bb(kl)
!do k=kl-1,1,-1
!  out(k)=(dd(k)-cc(k)*out(k+1))/bb(k)
!end do
!
!return
!end subroutine thomas1

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
qsat=0.622*esatf/max(ps-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

subroutine updatekmo(kmo,km,zz,zhl)

implicit none

integer k
real, dimension(ifull,kl), intent(in) :: km,zz
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl-1), intent(in) :: zhl
real xp

kmo(:,1)=km(:,2)+(zhl(:,1)-zz(:,2))/(zz(:,3)-zz(:,1))*               &
           ((zz(:,2)-zz(:,1))*(km(:,3)-km(:,2))/(zz(:,3)-zz(:,2))    &
           +(zz(:,3)-zz(:,2))*(km(:,2)-km(:,1))/(zz(:,2)-zz(:,1)))   &
         +(zhl(:,1)-zz(:,2))**2/(zz(:,3)-zz(:,1))*                   &
           ((km(:,3)-km(:,2))/(zz(:,3)-zz(:,2))                      &
           -(km(:,2)-km(:,1))/(zz(:,2)-zz(:,1)))
where (kmo(:,1).lt.2.E-7)
  kmo(:,1)=0.5*(km(:,1)+km(:,2))
end where
do k=2,kl-1
  kmo(:,k)=km(:,k)+(zhl(:,k)-zz(:,k))/(zz(:,k+1)-zz(:,k-1))*                 &
             ((zz(:,k)-zz(:,k-1))*(km(:,k+1)-km(:,k))/(zz(:,k+1)-zz(:,k))    &
             +(zz(:,k+1)-zz(:,k))*(km(:,k)-km(:,k-1))/(zz(:,k)-zz(:,k-1)))   &
           +(zhl(:,k)-zz(:,k))**2/(zz(:,k+1)-zz(:,k-1))*                     &
             ((km(:,k+1)-km(:,k))/(zz(:,k+1)-zz(:,k))                        &
             -(km(:,k)-km(:,k-1))/(zz(:,k)-zz(:,k-1)))
  where (kmo(:,k).lt.2.E-7)
    kmo(:,k)=0.5*(km(:,k)+km(:,k+1))
  end where
end do
! These terms are never used
kmo(:,kl)=2.E-7

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update mass flux coeffs at half levels

subroutine updategam(gamhl,gamin,zz,zhl,zi)

implicit none

integer k
real, dimension(ifull), intent(in) :: zi
real, dimension(ifull,kl), intent(in) :: gamin,zz
real, dimension(ifull,kl), intent(out) :: gamhl
real, dimension(ifull,kl-1), intent(in) :: zhl
real xp

gamhl=0.
gamhl(:,1)=gamin(:,2)+(zhl(:,1)-zz(:,2))/(zz(:,3)-zz(:,1))*                &
           ((zz(:,2)-zz(:,1))*(gamin(:,3)-gamin(:,2))/(zz(:,3)-zz(:,2))    &
           +(zz(:,3)-zz(:,2))*(gamin(:,2)-gamin(:,1))/(zz(:,2)-zz(:,1)))   &
         +(zhl(:,1)-zz(:,2))**2/(zz(:,3)-zz(:,1))*                         &
           ((gamin(:,3)-gamin(:,2))/(zz(:,3)-zz(:,2))                      &
           -(gamin(:,2)-gamin(:,1))/(zz(:,2)-zz(:,1)))
where (zhl(:,1).gt.zi)
  gamhl(:,1)=0.
end where
do k=2,kl-1
  gamhl(:,k)=gamin(:,k)+(zhl(:,k)-zz(:,k))/(zz(:,k+1)-zz(:,k-1))*                  &
             ((zz(:,k)-zz(:,k-1))*(gamin(:,k+1)-gamin(:,k))/(zz(:,k+1)-zz(:,k))    &
             +(zz(:,k+1)-zz(:,k))*(gamin(:,k)-gamin(:,k-1))/(zz(:,k)-zz(:,k-1)))   &
           +(zhl(:,k)-zz(:,k))**2/(zz(:,k+1)-zz(:,k-1))*                           &
             ((gamin(:,k+1)-gamin(:,k))/(zz(:,k+1)-zz(:,k))                        &
             -(gamin(:,k)-gamin(:,k-1))/(zz(:,k)-zz(:,k-1)))
  where (zhl(:,k).gt.zi)
    gamhl(:,k)=0.
  end where
end do
gamhl=max(gamhl,0.)

return
end subroutine updategam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE-eps

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
