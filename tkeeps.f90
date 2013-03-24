
! This module calculates the turblent kinetic energy and mixing for the boundary layer based on Hurley 2009
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

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:), allocatable, save :: zidry

! model constants
real, parameter :: b1      = 2.     ! Soares et al (2004) 1., Siebesma et al (2003) 2.
real, parameter :: b2      = 1./3.  ! Soares et al (2004) 2., Siebesma et al (2003) 1./3.
real, parameter :: be      = 7.     ! Angevine (2005) 7., Hurley (2007) 1., Soares et al (2004) 0.3
real, parameter :: cm0     = 0.09   ! Hurley (2007) 0.09, Duynkerke 1988 0.03
real, parameter :: ce0     = 0.69   ! Hurley (2007) 0.69, Duynkerke 1988 0.42
real, parameter :: ce1     = 1.46
real, parameter :: ce2     = 1.83
real, parameter :: ce3     = 0.45   ! Hurley (2007) 0.45, Dynkerke et al 1987 0.35
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
real, parameter :: pi    = 3.14159265

! stability constants
real, parameter :: a_1   = 1.
real, parameter :: b_1   = 2./3.
real, parameter :: c_1   = 5.
real, parameter :: d_1   = 0.35
!real, parameter :: aa1 = 3.8 ! Luhar low wind
!real, parameter :: bb1 = 0.5 ! Luhar low wind
!real, parameter :: cc1 = 0.3 ! Luhar low wind

integer, parameter :: icm1   = 5        ! max iterations for calculating pblh
real, parameter :: maxdts    = 300.     ! max timestep for split
real, parameter :: maxdtt    = 100.     ! max timestep for tke-eps
real, parameter :: mintke    = 1.E-8    ! min value for tke
real, parameter :: mineps    = 1.E-10   ! min value for eps
real, parameter :: minl      = 1.       ! min value for L (constraint on eps)
real, parameter :: maxl      = 1000.    ! max value for L (constraint on eps)

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
zidry=1000.

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qg,qlg,qfg,qrg,cfrac,cfrain,zi,wt0,wq0,ps,ustar,rhoas, &
                  zz,zzh,sig,sigkap,dt,qgmin,mode,diag,naero,aero)

implicit none

integer, intent(in) :: diag,mode,naero
integer k,i,j,klcl,ktopmax
integer icount,kcount,ncount,mcount
real, intent(in) :: dt,qgmin
real, dimension(:,:,:), intent(inout) :: aero
real, dimension(ifull,kl), intent(inout) :: theta,qg,qlg,qfg,qrg,cfrac,cfrain
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull,kl), intent(in) :: zz,zzh
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar,rhoas
real, dimension(kl), intent(in) :: sigkap,sig
real, dimension(ifull,kl,size(aero,3)) :: gamar
real, dimension(ifull,kl) :: km,thetav,thetal,temp,qsat
real, dimension(ifull,kl) :: gamtv,gamth,gamqv,gamql,gamqf
real, dimension(ifull,kl) :: gamqr,gamcf,gamcr,gamhl
real, dimension(ifull,kl) :: thetavnc,qsatc
real, dimension(ifull,kl) :: tkenew,epsnew,bb,cc,dd,ff,rr
real, dimension(ifull,kl) :: rhoa,rhoahl,thetavhl,thetahl
real, dimension(ifull,kl) :: qshl,qlhl,qfhl,qrhl
real, dimension(ifull,kl) :: pres,rmask
real, dimension(ifull,2:kl) :: idzm
real, dimension(ifull,1:kl-1) :: idzp
real, dimension(ifull,2:kl) :: aa,qq,pps,ppt,ppb
real, dimension(ifull,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,kl-1) :: fzzh
real, dimension(ifull) :: wstar,z_on_l,phim,wtv0,dum
real, dimension(ifull) :: tkeold,epsold
real, dimension(ifull) :: tbb,tcc,tff,tgg,tqq
real, dimension(ifull) :: qgnc
real, dimension(kl) :: w2up,nn
real, dimension(kl) :: qtup,qvup,qlup,qfup,qrup,qupsat
real, dimension(kl) :: mflx,thup,ttup,tvup,tlup,arup
real, dimension(1) :: tdum
real xp,dtr,as,bs,cs,cm12,cm34,qcup
real zht,dzht,nnc,zidryold,zilcl,ee
real w2c,mfc,ddtt,ddts
real lx,templ,fice,qxup,txup,dqsdt,al
logical, dimension(ifull,kl) :: lta
logical sconv

cm12=1./sqrt(cm0)
cm34=sqrt(sqrt(cm0**3))

if (diag>0) write(6,*) "Update PBL mixing with TKE-eps turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

do k=1,kl
  ! Impose limits after host advection
  tke(1:ifull,k)=max(tke(1:ifull,k),mintke)
  tff=cm34*tke(1:ifull,k)*sqrt(tke(1:ifull,k))/minl
  eps(1:ifull,k)=min(eps(1:ifull,k),tff)
  tff=max(tff*minl/maxl,mineps)
  eps(1:ifull,k)=max(eps(1:ifull,k),tff)

  ! Calculate air density - must use same theta for calculating dz
  pres(:,k)=ps(:)*sig(k) ! pressure
  temp(:,k)=theta(:,k)/sigkap(k)
  rhoa(:,k)=pres(:,k)/(rd*temp(:,k))

  ! Calculate first approximation to diffusion coeffs
  km(:,k)=cm0*tke(1:ifull,k)*tke(1:ifull,k)/eps(1:ifull,k)
end do

do k=1,kl-1
  ! Fraction for interpolation
  fzzh(:,k)=(zzh(:,k)-zz(:,k))/(zz(:,k+1)-zz(:,k))

  ! Calculate dz at half levels
  dz_hl(:,k)=zz(:,k+1)-zz(:,k)
end do

! Calculate dz at full levels
dz_fl(:,1)=zzh(:,1)
dz_fl(:,2:kl)=zzh(:,2:kl)-zzh(:,1:kl-1)

! Calculate shear term on full levels (see hordifg.f for calculation of horizontal shear)
pps(:,2:kl-1)=km(:,2:kl-1)*shear(:,2:kl-1)
pps(:,kl)=0.
ppb(:,kl)=0.
ppt(:,kl)=0.

! interpolate to half levels
call updatekmo(rhoahl,rhoa,fzzh)
call updatekmo(kmo,km,fzzh)

idzm(:,2:kl)=rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1)=rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

! Main loop to prevent time splitting errors
mcount=int(dt/(maxdts+0.01))+1
ddts=dt/real(mcount)
do kcount=1,mcount

  ! Set-up thermodynamic variables temp, theta_l, theta_v and surface fluxes
  thetav=theta*(1.+0.61*qg-qlg-qfg-qrg)
  wtv0=wt0+theta(:,1)*0.61*wq0
  !wtl0=wt0
  !wqt0=wq0
  tkeold=tke(1:ifull,1)
  epsold=eps(1:ifull,1)
  do k=1,kl
    temp(:,k)=theta(:,k)/sigkap(k)
    thetal(:,k)=theta(:,k)-sigkap(k)*(lv*(qlg(:,k)+qrg(:,k))+ls*qfg(:,k))/cp
    ! calculate saturated mixing ratio
    call getqsat(qsat(:,k),temp(:,k),pres(:,k))
  end do

  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is approximately constant in the plume (i.e., volume conserving)
  gamtv=0.
  gamth=0.
  gamqv=0.
  gamql=0.
  gamqf=0.
  gamqr=0.
  gamcf=0.
  gamcr=0.
  if (naero>0) then
    gamar=0.
  end if

  wstar=(grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)
  if (mode/=1) then ! mass flux
    do i=1,ifull
      if (wtv0(i)>0.) then ! unstable
        zidryold=zidry(i)
        do icount=1,icm1
          tke(i,1)=cm12*ustar(i)*ustar(i)+ce3*wstar(i)*wstar(i)
          tke(i,1)=max(tke(i,1),mintke)
          klcl=kl+1
          ktopmax=0
          w2up=0.
          zht=zz(i,1)
          dzht=zht
          ! Entrainment and detrainment rates (based on Angevine et al (2010), but scaled for CCAM)
          !ee=0.002                                     ! Angevine (2005)
          !ee=2./max(100.,zidry(i))                     ! Angevine et al (2010)
          !ee=1./max(zht,1.)                            ! Siebesma et al (2003)
          !ee=0.5*(1./max(zht,1.)+1./max(zi(i)-zht,1.)) ! Soares et al (2004)
          ee=1./max(100.,zidry(i))                      ! MJT suggestion for CCAM
          ! first level -----------------
          ! initial thermodynamic state
          ! split thetal and qtot into components (conservation is maintained)
          thup(1)=theta(i,1)+be*wt0(i)/sqrt(max(tke(i,1),1.E-4))   ! Hurley 2007
          qvup(1)=qg(i,1)+be*wq0(i)/sqrt(max(tke(i,1),1.E-4))      ! Hurley 2007
          qlup(1)=qlg(i,1)
          qfup(1)=qfg(i,1)
          qrup(1)=qrg(i,1)
          ! diagnose thermodynamic variables assuming no condensation
          tlup(1)=thup(1)-(lv*(qlup(1)+qrup(1))+ls*qfup(1))/cp ! thetal,up
          qtup(1)=qvup(1)+qlup(1)+qfup(1)+qrup(1)              ! qtot,up
          txup=tlup(1)                                                 ! theta,up after evaporation of ql,up and qf,up
          ttup(1)=txup/sigkap(1)                                       ! temp,up
          qxup=qtup(1)                                                 ! qv,up after evaporation of ql,up and qf,up
          call getqsat(qupsat(1:1),ttup(1:1),pres(i:i,1))              ! estimate of saturated mixing ratio in plume
          tvup(1)=txup+theta(i,1)*0.61*qxup                            ! thetav,up after evaporation of ql,up and qf,up
          ! update updraft velocity and mass flux
          nn(1)=grav*be*wtv0(i)/(thetav(i,1)*sqrt(max(tke(i,1),1.E-4)))          ! Hurley 2007
          !w2up(1)=2.*dzht*b2*nn(1)/(1.+2.*dzht*b1*ee(i))                        ! Hurley 2007
          w2up(1)=(0.25*wstar(i)*wstar(i)+2.*dzht*b2*nn(1))/(1.+2.*dzht*b1*ee)   ! Angevine et al (2010)
          mflx(1)=0.1*sqrt(w2up(1))                                              ! Hurley 2007
          ! check for lcl
          sconv=.false.
          if (qxup>=qupsat(1)) then
            sconv=.true.
            zilcl=zz(i,1)
            klcl=1
            w2c=w2up(1)
            mfc=mflx(1)
          end if
        
          ! dry convection case
          do k=2,kl
            dzht=dz_hl(i,k-1)
            zht=zz(i,k)
            ! update detrainment
            dtr=ee+0.025/max(zidry(i)-zht,1.)
            ! update thermodynamics of plume
            ! split thetal and qtot into components (conservation is maintained)
            ! (use upwind as centred scheme requires vertical spacing less than 250m)
            thup(k)=(thup(k-1)+dzht*ee*theta(i,k))/(1.+dzht*ee)
            qvup(k)=(qvup(k-1)+dzht*ee*qg(i,k)   )/(1.+dzht*ee)
            qlup(k)=(qlup(k-1)+dzht*ee*qlg(i,k)  )/(1.+dzht*ee)
            qfup(k)=(qfup(k-1)+dzht*ee*qfg(i,k)  )/(1.+dzht*ee)
            qrup(k)=(qrup(k-1)+dzht*ee*qrg(i,k)  )/(1.+dzht*ee)
            ! diagnose thermodynamic variables assuming no condensation
            tlup(k)=thup(k)-(lv*(qlup(k)+qrup(k))+ls*qfup(k))/cp  ! thetal,up
            qtup(k)=qvup(k)+qlup(k)+qfup(k)+qrup(k)               ! qtot,up
            txup=tlup(k)                                          ! theta,up after evaporation of ql,up and qf,up
            ttup(k)=txup/sigkap(k)                                ! temp,up
            qxup=qtup(k)                                          ! qv,up after evaporation of ql,up and qf,up
            tvup(k)=txup+theta(i,k)*0.61*qxup                     ! thetav,up after evaporation of ql,up and qf,up
            ! calculate buoyancy
            nn(k)=grav*(tvup(k)-thetav(i,k))/thetav(i,k)
            ! update updraft velocity
            w2up(k)=(w2up(k-1)+2.*dzht*b2*nn(k))/(1.+2.*dzht*b1*ee)
            ! update mass flux
            mflx(k)=mflx(k-1)/(1.+dzht*(dtr-ee))
            ! check for lcl
            if (.not.sconv) then
              call getqsat(qupsat(k:k),ttup(k:k),pres(i:i,k))     ! estimate of saturated mixing ratio in plume
              if (qxup>=qupsat(k)) then
                ! estimate LCL when saturation occurs
                as=ee*(qupsat(k)-qupsat(k-1))/dzht
                bs=(qupsat(k)-qupsat(k-1))/dzht+ee*(qupsat(k-1)-qg(i,k)-qlg(i,k)-qlg(i,k)-qrg(i,k))
                cs=qupsat(k-1)-qtup(k-1)
                xp=0.5*(-bs-sqrt(bs*bs-4.*as*cs))/as
                xp=min(max(xp,0.),dzht)
                sconv=.true.
                zilcl=xp+zz(i,k-1)
                klcl=k
                ! use dry convection to advect to lcl
                dtr=ee+0.05/max(zidry(i)-zilcl,1.)
                nnc=nn(k-1)+xp*(nn(k)-nn(k-1))/dzht
                w2c=(w2up(k-1)+2.*xp*b2*nnc)/(1.+2.*xp*b1*ee)
                mfc=mflx(k-1)/(1.+xp*(dtr-ee))
              end if
            end if
            ! test if maximum plume height is reached
            if (w2up(k)<=0.) then
              as=2.*b2*(nn(k)-nn(k-1))/dzht
              bs=2.*b2*nn(k-1)
              cs=w2up(k-1)
              xp=0.5*(-bs-sqrt(bs*bs-4.*as*cs))/as
              xp=min(max(xp,0.),dzht)
              zidry(i)=xp+zz(i,k-1)
              mflx(k)=0.
              ktopmax=max(ktopmax,k-1)
              exit
            end if
          end do

          if (sconv) then
            if (zidry(i)<zilcl) then
              sconv=.false.
              klcl=kl+1
            end if
          end if
          
          ! shallow convection case
          if (sconv) then
            ! advect from LCL to next model level
            dzht=zz(i,klcl)-zilcl
            zht=zz(i,klcl)
            xp=max(zi(i)-zilcl,1.)
            xp=8.*min(max(zht-zilcl,1.),xp)/xp-16./3.
            dtr=max(0.9*ee+0.003/pi*(atan(xp)+0.5*pi),ee)
            ! calculate conserved variables
            tlup(klcl)=thup(klcl)-(lv*(qlup(klcl)+qrup(klcl))+ls*qfup(klcl))/cp  ! thetal,up
            qtup(klcl)=qvup(klcl)+qlup(klcl)+qfup(klcl)+qrup(klcl)               ! qtot,up
            templ=tlup(klcl)/sigkap(klcl)                                        ! templ,up
            ! estimate saturated air temperature
            tdum=templ
            call getqsat(qupsat(klcl:klcl),tdum(1:1),pres(i:i,klcl))
            qxup=qupsat(klcl)
            fice=min(max(273.16-templ,0.),40.)/40.
            lx=lv+lf*fice
            dqsdt=qupsat(klcl)*lx/(rv*templ*templ)
            al=cp/(cp+lx*dqsdt)
            qcup=qtup(klcl)-qupsat(klcl)
            qcup=qcup*al
            ttup(klcl)=templ+lx*qcup/cp                                         ! temp,up
            txup=ttup(klcl)*sigkap(klcl)                                        ! theta,up after redistribution
            tvup(klcl)=txup+theta(i,klcl)*(1.61*qxup-qtup(klcl))                ! thetav,up after redistribution
            nn(klcl)=grav*(tvup(klcl)-thetav(i,klcl))/thetav(i,klcl)            ! buoyancy
            w2up(klcl)=(w2c+2.*dzht*b2*nn(klcl))/(1.+2.*dzht*b1*ee)
            mflx(klcl)=mfc/(1.+dzht*(dtr-ee))
            ! check for plume top
            if (w2up(klcl)<=0.) then
              as=2.*b2*(nn(klcl)-nnc)/dzht
              bs=2.*b2*nnc
              cs=w2c
              xp=0.5*(-bs-sqrt(bs*bs-4.*as*cs))/as
              xp=min(max(xp,0.),dzht)
              zi(i)=xp+zilcl
              mflx(klcl)=0.
              ktopmax=max(ktopmax,klcl-1)
            else
              do k=klcl+1,kl
                ! full level advection
                dzht=dz_hl(i,k-1)
                zht=zz(i,k)
                ! update detrainment rate for cloudly air from Angevine et al 2010
                xp=max(zi(i)-zilcl,1.)
                xp=8.*min(max(zht-zilcl,1.),xp)/xp-16./3.
                dtr=max(0.9*ee+0.003/pi*(atan(xp)+0.5*pi),ee)
                ! update thermodynamics of plume
                thup(k)=(thup(k-1)+dzht*ee*theta(i,k) )/(1.+dzht*ee)
                qvup(k)=(qvup(k-1)+dzht*ee*qg(i,k)    )/(1.+dzht*ee)
                qlup(k)=(qlup(k-1)+dzht*ee*qlg(i,k)   )/(1.+dzht*ee)
                qfup(k)=(qfup(k-1)+dzht*ee*qfg(i,k)   )/(1.+dzht*ee)
                qrup(k)=(qrup(k-1)+dzht*ee*qrg(i,k)   )/(1.+dzht*ee)
                ! calculate conserved variables
                tlup(k)=thup(k)-(lv*(qlup(k)+qrup(k))+ls*qfup(k))/cp  ! thetal,up
                qtup(k)=qvup(k)+qlup(k)+qfup(k)+qrup(k)               ! qtot,up
                templ=tlup(k)/sigkap(k)                               ! templ,up
                ! estimate saturated air temperature
                tdum=templ
                call getqsat(qupsat(k:k),tdum(1:1),pres(i:i,k))
                qxup=qupsat(k)
                fice=min(max(273.16-templ,0.),40.)/40.
                lx=lv+lf*fice
                dqsdt=qupsat(k)*lx/(rv*templ*templ)
                al=cp/(cp+lx*dqsdt)
                qcup=qtup(k)-qupsat(k)
                qcup=qcup*al
                ttup(k)=templ+lx*qcup/cp                               ! temp,up
                txup=ttup(k)*sigkap(k)                                 ! theta,up after redistribution
                tvup(k)=txup+theta(i,k)*(1.61*qxup-qtup(k))            ! thetav,up after redistribution
                ! calculate buoyancy
                nn(k)=grav*(tvup(k)-thetav(i,k))/thetav(i,k)
                ! update updraft velocity
                w2up(k)=(w2up(k-1)+2.*dzht*b2*nn(k))/(1.+2.*dzht*b1*ee)
                ! update mass flux
                mflx(k)=mflx(k-1)/(1.+dzht*(dtr-ee))
                ! test if maximum plume height is reached
                if (w2up(k)<=0.) then
                  as=2.*b2*(nn(k)-nn(k-1))/dzht
                  bs=2.*b2*nn(k-1)
                  cs=w2up(k-1)
                  xp=0.5*(-bs-sqrt(bs*bs-4.*as*cs))/as
                  xp=min(max(xp,0.),dzht)
                  zi(i)=xp+zz(i,k-1)
                  mflx(k)=0.
                  ktopmax=max(ktopmax,k-1)
                  exit
                end if
              end do
            end if
          else
            zi(i)=zidry(i)
          end if

          ! update surface boundary conditions
          wstar(i)=(grav*zi(i)*wtv0(i)/thetav(i,1))**(1./3.)

          ! check for convergence
          ! use zidry instead of zi to avoid oscillations
          ! for borderline shallow convection
          if (abs(zidry(i)-zidryold)<10.) exit
          zidryold=zidry(i)
        end do

        ! update explicit counter gradient terms
        do k=1,ktopmax
          gamth(i,k)=mflx(k)*(thup(k)-theta(i,k))
          gamqv(i,k)=mflx(k)*(qvup(k)-qg(i,k))
          gamql(i,k)=mflx(k)*(qlup(k)-qlg(i,k))
          gamqf(i,k)=mflx(k)*(qfup(k)-qfg(i,k))
          gamqr(i,k)=mflx(k)*(qrup(k)-qrg(i,k))
          gamtv(i,k)=gamth(i,k)+theta(i,k)*(0.61*gamqv(i,k)-gamql(i,k)-gamqf(i,k)-gamqr(i,k))
        end do

        ! update reamining scalars which are not used in the iterative loop
        arup(1)=cfrac(i,1)
        gamcf(i,1)=mflx(1)*(arup(1)-cfrac(i,1))
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          arup(k)=(arup(k-1)+dzht*ee*cfrac(i,k) )/(1.+dzht*ee)
          gamcf(i,k)=mflx(k)*(arup(k)-cfrac(i,k))
        end do
        arup(1)=cfrain(i,1)
        gamcr(i,1)=mflx(1)*(arup(1)-cfrain(i,1))
        do k=2,ktopmax
          dzht=dz_hl(i,k-1)
          arup(k)=(arup(k-1)+dzht*ee*cfrain(i,k))/(1.+dzht*ee)
          gamcr(i,k)=mflx(k)*(arup(k)-cfrain(i,k))
        end do
        do j=1,naero
          arup(1)=aero(i,1,j)
          gamar(i,1,j)=mflx(1)*(arup(1)-aero(i,1,j))
          do k=2,ktopmax
            dzht=dz_hl(i,k-1)
            arup(k)=(arup(k-1)+dzht*ee*aero(i,k,j))/(1.+dzht*ee)
            gamar(i,k,j)=mflx(k)*(arup(k)-aero(i,k,j))
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
        zidry(i)=zi(i)
      end if
    end do
       
  else ! other PBL scheme
    zidry=zi
  end if

  ! calculate tke and eps at 1st level
  z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-20))
  z_on_l=min(z_on_l,10.)
  where (z_on_l<0.)
    phim=(1.-16.*z_on_l)**(-0.25)
  elsewhere !(z_on_l<=0.4)
    phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
  !elsewhere
  !  phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar (2007)
  end where
  tke(1:ifull,1)=cm12*ustar*ustar+ce3*wstar*wstar
  eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
  tke(1:ifull,1)=max(tke(1:ifull,1),mintke)
  tff=cm34*tke(1:ifull,1)*sqrt(tke(1:ifull,1))/minl
  eps(1:ifull,1)=min(eps(1:ifull,1),tff)
  tff=max(tff*minl/maxl,mineps)
  eps(1:ifull,1)=max(eps(1:ifull,1),tff)

  ! Calculate sources and sinks for TKE and eps
  ! prepare arrays for calculating buoyancy of saturated air
  ! (i.e., related to the saturated adiabatic lapse rate)
  qsatc=max(qsat,qg)                           ! assume qg is saturated inside cloud
  dd=qlg/max(cfrac,1.E-6)                      ! inside cloud value
  ff=qfg/max(cfrac,1.E-6)                      ! inside cloud value
  rr=qrg/max(cfrac,cfrain,1.E-6)               ! inside cloud value assuming max overlap
  do k=1,kl
    tbb=max(1.-cfrac(:,k),1.E-6)
    qgnc=(qg(:,k)-(1.-tbb)*qsatc(:,k))/tbb     ! outside cloud value
    qgnc=min(max(qgnc,qgmin),qsatc(:,k))
    thetavnc(:,k)=theta(:,k)*(1.+0.61*qgnc)    ! outside cloud value
  end do
  call updatekmo(thetahl,theta,fzzh)
  call updatekmo(thetavhl,thetavnc,fzzh)       ! outside cloud value
  call updatekmo(qshl,qsatc,fzzh)              ! inside cloud value
  call updatekmo(qlhl,dd,fzzh)                 ! inside cloud value
  call updatekmo(qfhl,ff,fzzh)                 ! inside cloud value
  call updatekmo(qrhl,rr,fzzh)                 ! inside cloud value
  ! fixes for clear/cloudy interface
  lta(:,2:kl)=cfrac(:,2:kl)<=1.E-6
  do k=2,kl-1
    do i=1,ifull
      if (lta(i,k).and..not.lta(i,k+1)) then
        if (lta(i,k)) then
          qlhl(i,k)=dd(i,k+1)
          qfhl(i,k)=ff(i,k+1)
          qrhl(i,k)=rr(i,k+1)
        else
          qlhl(i,k)=dd(i,k)
          qfhl(i,k)=ff(i,k)
          qrhl(i,k)=rr(i,k)
        end if
      end if
    end do
  end do

  ! Calculate buoyancy term
  ! saturated conditions from Durran and Klemp JAS 1982 (see also WRF)
  do k=2,kl-1
    tqq=(1.+lv*qsatc(:,k)/(rd*temp(:,k)))/(1.+lv*lv*qsatc(:,k)/(cp*rv*temp(:,k)*temp(:,k)))
    tbb=-grav*km(:,k)*(tqq*((thetahl(:,k)-thetahl(:,k-1))/theta(:,k)                              &
           +lv*(qshl(:,k)-qshl(:,k-1))/(cp*temp(:,k)))-qshl(:,k)-qlhl(:,k)-qfhl(:,k)-qrhl(:,k)    &
           +qshl(:,k-1)+qlhl(:,k-1)+qfhl(:,k-1)+qrhl(:,k-1))/dz_fl(:,k)
    tbb=tbb+grav*(tqq*(gamth(:,k)/theta(:,k)+lv*gamqv(:,k)/(cp*temp(:,k)))         &
           -gamqv(:,k)-gamql(:,k)-gamqf(:,k)-gamqr(:,k))
    ! unsaturated
    tcc=-grav*km(:,k)*(thetavhl(:,k)-thetavhl(:,k-1))/(thetavnc(:,k)*dz_fl(:,k))
    tcc=tcc+grav*gamtv(:,k)/thetav(:,k)
    ppb(:,k)=(1.-cfrac(:,k))*tcc+cfrac(:,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)

    ! Calculate transport term on full levels
    ppt(:,k)=   kmo(:,k)*idzp(:,k)*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k)   &
             -kmo(:,k-1)*idzm(:,k)*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1)
  end do

  ! Update TKE and eps terms
  ncount=int(ddts/(maxdtt+0.01))+1
  ddtt=ddts/real(ncount)
  qq(:,2:kl-1)=-ddtt*idzm(:,2:kl-1)/dz_hl(:,1:kl-2)
  rr(:,2:kl-1)=-ddtt*idzp(:,2:kl-1)/dz_hl(:,2:kl-1)
  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tke(1:ifull,kl)=mintke
  eps(1:ifull,kl)=mineps
  do icount=1,ncount
    xp=real(icount)/real(ncount)

    ! eps vertical mixing (done here as we skip level 1, instead of using trim)
    aa(:,2:kl-1)=ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)
    cc(:,2:kl-1)=ce0*kmo(:,2:kl-1)*rr(:,2:kl-1)
    bb(:,2:kl-1)=ddtt*ce2*eps(1:ifull,2:kl-1)/tke(1:ifull,2:kl-1)
    bb(:,2)=bb(:,2)-aa(:,2)
    bb(:,kl-1)=bb(:,kl-1)-cc(:,kl-1)
    dd(:,2:kl-1)=eps(1:ifull,2:kl-1)+ddtt*ce1*eps(1:ifull,2:kl-1)/tke(1:ifull,2:kl-1) &
                *(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
    dd(:,2)=dd(:,2)-aa(:,2)*(epsold+xp*(eps(1:ifull,1)-epsold))
    dd(:,kl-1)=dd(:,kl-1)-cc(:,kl-1)*mineps
    call thomas(epsnew(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),kl-2)

    ! TKE vertical mixing (done here as we skip level 1, instead of using trim)
    aa(:,2:kl-1)=kmo(:,1:kl-2)*qq(:,2:kl-1)
    cc(:,2:kl-1)=kmo(:,2:kl-1)*rr(:,2:kl-1)
    bb(:,2)=-aa(:,2)
    bb(:,3:kl-2)=0.
    bb(:,kl-1)=-cc(:,kl-1)
    dd(:,2:kl-1)=tke(1:ifull,2:kl-1)+ddtt*(pps(:,2:kl-1)+ppb(:,2:kl-1)-epsnew(:,2:kl-1))
    dd(:,2)=dd(:,2)-aa(:,2)*(tkeold+xp*(tke(1:ifull,1)-tkeold))
    dd(:,kl-1)=dd(:,kl-1)-cc(:,kl-1)*mintke
    call thomas(tkenew(:,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),kl-2)

    do k=2,kl-1
      tke(1:ifull,k)=max(tkenew(:,k),mintke)
      tff=cm34*tke(1:ifull,k)*sqrt(tke(1:ifull,k))/minl
      eps(1:ifull,k)=min(epsnew(:,k),tff)
      tff=max(tff*minl/maxl,mineps)
      eps(1:ifull,k)=max(eps(1:ifull,k),tff)
    end do
    
    km=cm0*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:)
    call updatekmo(kmo,km,fzzh) ! interpolate diffusion coeffs to half levels

  end do

  ! map grid boxes under the boundary layer height
  do k=1,kl-1
    where (zi>=zzh(:,k))
      rmask(:,k)=1.
    elsewhere
      rmask(:,k)=0.
    end where
  end do
  rmask(:,kl)=0.

  ! updating diffusion and non-local terms for qtot and thetal
  aa(:,2:kl)=-ddts*kmo(:,1:kl-1)*idzm(:,2:kl)/dz_hl(:,1:kl-1)
  cc(:,1:kl-1)=-ddts*kmo(:,1:kl-1)*idzp(:,1:kl-1)/dz_hl(:,1:kl-1)
  bb(:,1:kl)=0.
  
  ! update non-local terms
  call updatekmo(gamhl,gamth,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=theta(:,1)-ddts*gamhl(:,1)*idzp(:,1)+ddts*wt0*rhoas/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=theta(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=theta(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(theta,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)
  
  call updatekmo(gamhl,gamqv,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=qg(:,1)-ddts*gamhl(:,1)*idzp(:,1)+ddts*wq0*rhoas/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=qg(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qg(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)

  call updatekmo(gamhl,gamql,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=qlg(:,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qlg(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qlg(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)

  call updatekmo(gamhl,gamqf,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=qfg(:,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qfg(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qfg(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)

  call updatekmo(gamhl,gamqr,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=qrg(:,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=qrg(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qrg(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(qrg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)

  ! account for phase transitions
  do k=1,kl
    rr(:,k)=theta(:,k)-sigkap(k)*(lv*(qlg(:,k)+qrg(:,k))+ls*qfg(:,k))/cp ! thetal
    tgg=max(qg(:,k)+qlg(:,k)+qfg(:,k)+qrg(:,k),qgmin) ! qtot before phase transition
    qlg(:,k)=max(qlg(:,k),0.)
    qfg(:,k)=max(qfg(:,k),0.)
    qrg(:,k)=max(qrg(:,k),0.)
    qg(:,k)=max(qg(:,k),0.)
    tff=max(qg(:,k)+qlg(:,k)+qfg(:,k)+qrg(:,k),qgmin) ! qtot after phase transition
    tgg=tgg/tff                                       ! scale factor for conservation
    qg(:,k)=qg(:,k)*tgg
    qlg(:,k)=qlg(:,k)*tgg
    qfg(:,k)=qfg(:,k)*tgg
    qrg(:,k)=qrg(:,k)*tgg
    theta(:,k)=rr(:,k)+sigkap(k)*(lv*(qlg(:,k)+qrg(:,k))+ls*qfg(:,k))/cp
  end do

  ! update cloud fraction terms
  call updatekmo(gamhl,gamcf,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=cfrac(:,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=cfrac(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=cfrac(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(cfrac,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)
  cfrac=min(max(cfrac,0.),1.)
  where (qlg+qfg>1.E-12)
    cfrac=max(cfrac,1.E-6)
  end where

  call updatekmo(gamhl,gamcr,fzzh)
  gamhl=gamhl*rmask
  dd(:,1)=cfrain(:,1)-ddts*gamhl(:,1)*idzp(:,1)
  dd(:,2:kl-1)=cfrain(:,2:kl-1)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=cfrain(:,kl)+ddts*gamhl(:,kl-1)*idzm(:,kl)
  call thomas(cfrain,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)
  cfrain=min(max(cfrain,0.),1.)
  where (qrg>1.E-12)
    cfrain=max(cfrain,1.E-6)
  end where
  
  ! Aerosols
  do j=1,naero
    call updatekmo(gamhl,gamar(:,:,j),fzzh)
    gamhl=gamhl*rmask
    dd(:,1)=aero(:,1,j)-ddts*gamhl(:,1)*idzp(:,1)
    dd(:,2:kl-1)=aero(:,2:kl-1,j)+ddts*(gamhl(:,1:kl-2)*idzm(:,2:kl-1)-gamhl(:,2:kl-1)*idzp(:,2:kl-1))
    dd(:,kl)=aero(:,kl,j)+ddts*gamhl(:,kl-1)*idzm(:,kl)
    call thomas(aero(:,:,j),aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),kl)
  end do

end do

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

subroutine thomas(out,aa,xtr,cc,ddi,klin)

implicit none

integer, intent(in) :: klin
real, dimension(ifull,2:klin), intent(in) :: aa
real, dimension(ifull,1:klin), intent(in) :: xtr,ddi
real, dimension(ifull,1:klin-1), intent(in) :: cc
real, dimension(ifull,1:klin), intent(out) :: out
real, dimension(ifull,1:klin) :: bb,dd
real, dimension(ifull) :: n
integer k

bb=xtr
bb(:,1)=1.-cc(:,1)+bb(:,1)
bb(:,2:klin-1)=1.-aa(:,2:klin-1)-cc(:,2:klin-1)+bb(:,2:klin-1)
bb(:,klin)=1.-aa(:,klin)+bb(:,klin)
dd=ddi

do k=2,klin
  n=aa(:,k)/bb(:,k-1)
  bb(:,k)=bb(:,k)-n*cc(:,k-1)
  dd(:,k)=dd(:,k)-n*dd(:,k-1)
end do
out(:,klin)=dd(:,klin)/bb(:,klin)
do k=klin-1,1,-1
  out(:,k)=(dd(:,k)-cc(:,k)*out(:,k+1))/bb(:,k)
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

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

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
