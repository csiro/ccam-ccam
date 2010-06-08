
! This module calculates the turblent kinetic energy and mixing for the boundary layer based on
! Hurley 2009.  Specifically, this version is adapted for the conformal cubic coordinate system.

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
public ww,dwdx,dwdy

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:,:), allocatable, save :: tkesav,epssav
real, dimension(:,:), allocatable, save :: ww,dwdx,dwdy

! model constants
real, parameter :: cm      = 0.09
real, parameter :: ce0     = 0.69
real, parameter :: ce1     = 1.46
real, parameter :: ce2     = 1.83
real, parameter :: aup     = 0.1
real, parameter :: b1      = 1.
real, parameter :: b2      = 2.
real, parameter :: entr    = 2.E-3
real, parameter :: detr    = 3.E-3
!real, parameter :: cq      = 2.5
real, parameter :: rcrit_l = 0.75
real, parameter :: rcrit_s = 0.85

! physical constants
real, parameter :: grav = 9.80616
real, parameter :: lv = 2.5104e6
real, parameter :: lf = 3.36e5
real, parameter :: ls = lv+lf
real, parameter :: epsl = 1./1.61
real, parameter :: delta=1./(epsl-1.)
real, parameter :: rd = 287.04
real, parameter :: rv = 461.5
real, parameter :: cp = 1004.64
real, parameter :: vkar = 0.4

! stability constants
real, parameter :: a_1   = 1.0
real, parameter :: b_1   = 2.0/3.0
real, parameter :: c_1   = 5.0
real, parameter :: d_1   = 0.35
!real, parameter :: aa1 = 3.8 ! Luhar low wind
!real, parameter :: bb1 = 0.5 ! Luhar low wind
!real, parameter :: cc1 = 0.3 ! Luhar low wind

integer, parameter :: buoymeth = 1 ! 0=Dry air, 1=Duynkerke, 2=Smith

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
allocate(shear(ifull,kl),ww(ifull+iextra,kl))
allocate(dwdx(ifull,kl),dwdy(ifull,kl))

tke=1.5E-4
eps=1.0E-6
tkesav=tke(1:ifull,:)
epssav=eps(1:ifull,:)
shear=0.
ww=0.
dwdx=0.
dwdy=0.

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux
! mode=2 mass flux without moist convection (no counter-gradient for moisture)

subroutine tkemix(kmo,theta,qg,qlg,qfg,u,v,cfrac,zi,land,wt0,wq0,ps,ustar,zz,sig,sigkap,dt,mode,diag)

implicit none

integer, intent(in) :: diag,mode
integer k,i,ktop,kbot
real, intent(in) :: dt
real, dimension(ifull,kl), intent(inout) :: theta,qg
real, dimension(ifull,kl), intent(in) :: u,v,zz,qlg,qfg,cfrac
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar
real, dimension(kl), intent(in) :: sigkap,sig
real, dimension(ifull,kl) :: km,gamt,ff,gg,templ,thetav,temp,gamq
real, dimension(ifull,kl) :: tkenew,epsnew,qsat,ppb
real, dimension(ifull,2:kl) :: pps,ppt
real, dimension(ifull,kl) :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,2:kl) :: aa,bb,cc,dd
real, dimension(ifull) :: wstar,z_on_l,phim,wtv0,hh,jj,dqsdt
real, dimension(kl-1) :: wpv_flux,tup,qup,w2up
real xp,wup,mflx,qupsat(1),ee
real cf,qc,rcrit,delq
real dz_ref,cm34
logical, dimension(ifull), intent(in) :: land
logical sconv
integer icount
integer, parameter :: icm = 10
real, parameter :: tol=0.001
!real, parameter :: cm34=cm**0.75 ! Crashes the SX6 cross-compiler

cm34=cm**0.75

if (diag.gt.0) write(6,*) "Update PBL mixing with TKE"

thetav=theta*(1.+0.61*qg)
wtv0=wt0+0.61*theta(:,1)*wq0

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

! impose limits after host advection
tke(1:ifull,:)=max(tke(1:ifull,:),1.5E-4)
tke(1:ifull,:)=min(tke(1:ifull,:),65.)
ff=cm34*(tke(1:ifull,:)**1.5)/5.
eps(1:ifull,:)=min(eps(1:ifull,:),ff)
ff=max(ff*5./500.,1.E-6)
eps(1:ifull,:)=max(eps(1:ifull,:),ff)

! Calculate diffusion coeffs
km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-3)

! calculate saturated mixing ratio
do k=1,kl
  temp(:,k)=theta(:,k)/sigkap(k)
  call getqsat(ifull,qsat(:,k),temp(:,k),ps*sig(k))
end do

! Calculate buoyancy terms (gamt included later)
select case(buoymeth)

  case(0) ! Hurley
    ppb(:,1)=-grav/thetav(:,1)*(thetav(:,2)-thetav(:,1))/dz_hl(:,1)
    do k=2,kl-1
      ppb(:,k)=-grav/thetav(:,k)*0.5*(thetav(:,k+1)-thetav(:,k-1))/dz_fl(:,k)
    end do
    ppb(:,kl)=-grav/thetav(:,kl)*(thetav(:,kl)-thetav(:,kl-1))/dz_hl(:,kl-1)

  case(1) ! Saturated conditions based on WRF formulation
    do k=1,kl
      ff(:,k)=theta(:,k)*(1.+lv*qsat(:,k)/(cp*temp(:,k)))
      gg(:,k)=(1.+1.61*epsl*lv*qg(:,k)/(rd*temp(:,k))) &
             /(1.+epsl*lv*lv*qg(:,k)/(cp*1.61*rd*temp(:,k)*temp(:,k)))/theta(:,k)
    end do
    !where (qg(:,1).lt.qsat(:,1).and.cfrac(:,1).lt.0.5)
    !  ppb(:,1)=-grav/thetav(:,1)*(thetav(:,2)-thetav(:,1))/dz_hl(:,1) &
    !            +grav*(qlg(:,2)+qfg(:,2)-qlg(:,1)-qfg(:,1))/dz_hl(:,1) ! same as nvmix=5
    !elsewhere
    !  ppb(:,1)=-grav*gg(:,1)*(ff(:,2)-ff(:,1))/dz_hl(:,1) &
    !            +grav*(qg(:,2)+qlg(:,2)+qfg(:,2)-qg(:,1)-qlg(:,1)-qfg(:,1))/dz_hl(:,1)
    !end where
    ppb(:,1)=-grav/dz_hl(:,1)*((1.-cfrac(:,1))*(thetav(:,2)-thetav(:,1))/thetav(:,1) &
                              +cfrac(:,1)*(ff(:,2)-ff(:,1))*gg(:,1)) &
             +grav/dz_hl(:,1)*cfrac(:,1)*(qg(:,2)-qg(:,1)) &
             +grav/dz_hl(:,1)*(qlg(:,2)+qfg(:,2)-qlg(:,1)-qfg(:,1))
    do k=2,kl-1
      !where (qg(:,k).lt.qsat(:,k).and.cfrac(:,k).lt.0.5)
      !  ppb(:,k)=-grav/thetav(:,k)*0.5*(thetav(:,k+1)-thetav(:,k-1))/dz_fl(:,k) &
      !           +grav*0.5*(qlg(:,k+1)+qfg(:,k+1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k) ! same as nvmix=5
      !elsewhere
      !  ppb(:,k)=-grav*gg(:,k)*0.5*(ff(:,k+1)-ff(:,k-1))/dz_fl(:,k) &
      !           +grav*0.5*(qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-qg(:,k-1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
      !end where
      ppb(:,k)=-grav*0.5/dz_fl(:,k)*((1.-cfrac(:,k))*(thetav(:,k+1)-thetav(:,k-1))/thetav(:,k) &
                                    +cfrac(:,k)*(ff(:,k+1)-ff(:,k-1))*gg(:,k)) &
               +grav*0.5/dz_fl(:,k)*cfrac(:,k)*(qg(:,k+1)-qg(:,k-1)) &
               +grav*0.5/dz_fl(:,k)*(qlg(:,k+1)+qfg(:,k+1)-qlg(:,k-1)-qfg(:,k-1))
    end do
    !where (qg(:,kl).lt.qsat(:,kl).and.cfrac(:,kl).lt.0.5)
    !  ppb(:,kl)=-grav/thetav(:,kl)*(thetav(:,kl)-thetav(:,kl-1))/dz_hl(:,kl-1) &
    !            +grav*(qlg(:,kl)+qfg(:,kl)-qlg(:,kl-1)-qfg(:,kl-1))/dz_hl(:,kl-1) ! same as nvmix=5
    !elsewhere
    !  ppb(:,kl)=-grav*gg(:,kl)*(ff(:,kl)-ff(:,kl-1))/dz_hl(:,kl-1) &
    !            +grav*(qg(:,kl)+qlg(:,kl)+qfg(:,kl)-qg(:,kl-1)-qlg(:,kl-1)-qfg(:,kl-1))/dz_hl(:,kl-1)
    !end where
    ppb(:,kl)=-grav/dz_hl(:,kl-1)*((1.-cfrac(:,kl))*(thetav(:,kl)-thetav(:,kl-1))/thetav(:,kl) &
                                  +cfrac(:,kl)*(ff(:,kl)-ff(:,kl-1))*gg(:,kl)) &
              +grav/dz_hl(:,kl-1)*cfrac(:,kl)*(qg(:,kl)-qg(:,kl-1)) &
              +grav/dz_hl(:,kl-1)*(qlg(:,kl)+qfg(:,kl)-qlg(:,kl-1)-qfg(:,kl-1))
              
    
  case(2) ! Smith (1990)
    do k=1,kl
      templ(:,k)=temp(:,k)-(lv/cp*qlg(:,k)+ls/cp*qfg(:,k))
      jj(:)=lv+lf*qfg(:,k)/max(qlg(:,k)+qfg(:,k),1.e-12) ! L
      dqsdt(:)=epsl*jj(:)*qsat(:,k)/(rd*temp(:,k)**2)
      hh(:)=cfrac(:,k)*(jj(:)/cp/temp(:,k) &
                        -delta/(1.-epsl)/(1.+delta*qg(:,k)-qlg(:,k)-qfg(:,k)))/(1.+jj(:)/cp*dqsdt) ! betac
      ff(:,k)=1./temp(:,k)-dqsdt*hh(:)                         ! betatt
      gg(:,k)=delta/(1.+delta*qg(:,k)-qlg(:,k)-qfg(:,k))+hh(:) ! betaqt
    end do    
    ppb(:,1)=-grav*ff(:,1)*((templ(:,2)-templ(:,1))/dz_hl(:,1)+grav/cp) &
             -grav*gg(:,1)*(qg(:,2)+qlg(:,2)+qfg(:,2)-qg(:,1)-qlg(:,1)-qfg(:,1))/dz_hl(:,1)
    do k=2,kl-1
      ppb(:,k)=-grav*ff(:,k)*(0.5*(templ(:,k+1)-templ(:,k-1))/dz_fl(:,k)+grav/cp) &
               -grav*gg(:,k)*0.5*(qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-qg(:,k-1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
    end do
    ppb(:,kl)=-grav*ff(:,kl)*((templ(:,kl)-templ(:,kl-1))/dz_hl(:,kl-1)+grav/cp) &
              -grav*gg(:,kl)*(qg(:,kl)+qlg(:,kl)+qfg(:,kl)-qg(:,kl-1)-qlg(:,kl-1)-qfg(:,kl-1))/dz_hl(:,kl-1)
      
  case DEFAULT
    write(6,*) "ERROR: Unsupported buoyancy option"
    stop
    
end select

! Calculate non-local terms for theta_v
gamt=0.
gamq=0.

if (mode.ne.1) then ! mass flux when mode is an even number
  zi=zz(:,kl)
  do i=1,ifull
    if (wtv0(i).gt.0.) then ! unstable
      sconv=.false.
      ee=0.5*(1./zz(i,1)+1./(zi(i)+zz(i,1)))    ! dz_hl(:,0)=zz(:,1)
      tup(1)=thetav(i,1)+wtv0(i)/sqrt(tke(i,1)) !*0.3 in Soares et al (2004)
      qup(1)=qg(i,1)+wq0(i)/sqrt(tke(i,1))      !*0.3 in Soares et al (2004)
      w2up(1)=(2.*zz(i,1)*b2*grav*(tup(1)-thetav(i,1))/thetav(i,1))/(1.+zz(i,1)*b1*ee)
      wup=sqrt(max(w2up(1),0.))
      aa(1,2)=ps(i)*sig(1)
      bb(1,2)=tup(1)/((1.+0.61*qup(1))*sigkap(1))
      call getqsat(1,qupsat(1),bb(1,2),aa(1,2))
      cf=0.
      if (qup(1).lt.qupsat(1).or.mode.eq.2) then ! dry convection
        mflx=aup*wup 
      else                                       ! moist convection (boundary condition)
        sconv=.true.
        qc=qup(1)-qupsat(1) ! cf calculation from LDR (1996)
        if (land(i)) then
          rcrit=max(rcrit_l,sig(1)**3)
        else
          rcrit=max(rcrit_s,sig(1)**3)
        end if
        delq=(1.-rcrit)*qupsat(1)
        if (qc.ge.delq) then
          cf=1.
        else if (qc.gt.0.) then
          cf=max(1.E-6,1.-0.5*((qc-delq)/delq)**2)
        else if (qc.gt.-delq) then
          cf=max(1.E-6,0.5*((qc+delq)/delq)**2)
        else
          cf=0.
        end if
        mflx=cf*aup*wup
      end if
      gamt(i,1)=mflx*(tup(1)-thetav(i,1))
      gamq(i,1)=mflx*(qup(1)-qg(i,1))
      do k=2,kl-1
        ee=0.5*(1./(zz(i,k-1)+dz_hl(i,k-1))+1./(max(zi(i)-zz(i,k-1),0.)+dz_hl(i,k-1)))
        tup(k)=(tup(k-1)*(1.-0.5*dz_hl(i,k-1)*ee)+0.5*dz_hl(i,k-1)*ee*(thetav(i,k)+thetav(i,k-1))) &
             /(1.+0.5*dz_hl(i,k-1)*ee)
        qup(k)=(qup(k-1)*(1.-0.5*dz_hl(i,k-1)*ee)+0.5*dz_hl(i,k-1)*ee*(qg(i,k)+qg(i,k-1))) &
             /(1.+0.5*dz_hl(i,k-1)*ee)
        w2up(k)=(w2up(k-1)*(1.-dz_hl(i,k-1)*b1*ee)+2.*dz_hl(i,k-1)*b2*grav*(tup(k)+tup(k-1)-thetav(i,k)-thetav(i,k-1)) &
              /(thetav(i,k)+thetav(i,k-1)))/(1.+dz_hl(i,k-1)*b1*ee)
        wup=sqrt(max(w2up(k),0.))
        if (.not.sconv) then
          aa(1,2)=ps(i)*sig(k)
          bb(1,2)=tup(k)/((1.+0.61*qup(k))*sigkap(k))
          call getqsat(1,qupsat(1),bb(1,2),aa(1,2))
          cf=0.
          if (qup(k).lt.qupsat(1).or.mode.eq.2) then ! dry convection
            mflx=aup*wup
          else                                       ! moist convection (boundary condition)
            sconv=.true.
            !sigup=max(1.E-6,-1.6*tke(i,k)/eps(i,k)*cq*km(i,k)*((qup(k)-qup(k-1))/dz_hl(i,k-1))**2
            !cf=0.5+0.36*arctan(1.55*(qup(k)-qupsat(1))/sqrt(sigup)) ! Cuijpers and Bechtold (1995)
            qc=qup(k)-qupsat(1) ! cf calculation from LDR (1996)
            if (land(i)) then
              rcrit=max(rcrit_l,sig(k)**3)
            else
              rcrit=max(rcrit_s,sig(k)**3)
            end if
            delq=(1.-rcrit)*qupsat(1)
            if (qc.ge.delq) then
              cf=1.
            else if (qc.gt.0.) then
              cf=max(1.E-6,1.-0.5*((qc-delq)/delq)**2)
            else if (qc.gt.-delq) then
              cf=max(1.E-6,0.5*((qc+delq)/delq)**2)
            else
              cf=0.
            end if
            mflx=cf*aup*wup
          end if
        else
          if (w2up(k).gt.0.) then       ! moist convection
            mflx=mflx*exp((entr-detr)*dz_hl(i,k-1))
          else
            mflx=0.
          end if
        end if
        gamt(i,k)=mflx*(tup(k)-thetav(i,k))
        gamq(i,k)=mflx*(qup(k)-qg(i,k))
        if (w2up(k).le.0.) then
          xp=w2up(k-1)/(w2up(k-1)-w2up(k))
          xp=min(max(xp,0.),1.)
          zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
          exit
        end if
      end do
    else                   ! stable
      wpv_flux(1)=km(i,1)*ppb(i,1)*thetav(i,1)/grav !+gamt(i,1)
      if (wpv_flux(1).lt.0.) then
        zi(i)=zz(i,1)
      else
        do k=2,kl-1
          wpv_flux(k)=km(i,k)*ppb(i,k)*thetav(i,k)/grav !+gamt(i,k)
          if (abs(wpv_flux(k)).lt.0.05*abs(wpv_flux(1))) then
            xp=(0.05*abs(wpv_flux(1))-abs(wpv_flux(k-1)))/(abs(wpv_flux(k))-abs(wpv_flux(k-1)))
            xp=min(max(xp,0.),1.)
            zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
            exit
          end if
          if (wpv_flux(k).le.0.) then
            xp=-wpv_flux(k-1)/(wpv_flux(k)-wpv_flux(k-1))
            xp=min(max(xp,0.),1.)
            zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
            exit
          end if
        end do
      end if
    end if
  end do

  gamt=km*min(max(gamt/km,0.),2.E-3)
  gamq=km*min(max(gamq/km,0.),2.E-5)
  
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

! boundary condition
tke(1:ifull,1)=ustar*ustar/sqrt(cm)+0.45*wstar*wstar
eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
tke(1:ifull,1)=max(tke(1:ifull,1),1.5E-4)
tke(1:ifull,1)=min(tke(1:ifull,1),65.)
aa(:,2)=cm34*(tke(1:ifull,1)**1.5)/5.
eps(1:ifull,1)=min(eps(1:ifull,1),aa(:,2))
aa(:,2)=max(aa(:,2)*5./500.,1.E-6)
eps(1:ifull,1)=max(eps(1:ifull,1),aa(:,2))

! Calculate buoyancy terms (include gamt)
if (mode.ne.1) then
  ppb(:,2:kl)=ppb(:,2:kl)+(grav/thetav(:,2:kl))*gamt(:,2:kl)/km(:,2:kl)
end if

! Calculate shear and transport terms
do k=2,kl-1
  pps(:,k)=((u(:,k+1)-u(:,k-1))/dz_fl(:,k)+dwdx(:,k))**2+((v(:,k+1)-v(:,k-1))/dz_fl(:,k)+dwdy(:,k))**2 &
          +(2.*(ww(1:ifull,k+1)-ww(1:ifull,k-1))/dz_fl(:,k))**2+shear(:,k)/km(:,k)
  where (wt0.le.0.)
    ppt(:,k)=0.5*((km(:,k+1)+km(:,k))*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k) &
                 -(km(:,k)+km(:,k-1))*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1))/dz_fl(:,k)
  elsewhere
    ppt(:,k)=0.
  end where  
end do
pps(:,kl)=((u(:,kl)-u(:,kl-1))/dz_hl(:,kl-1)+dwdx(:,kl))**2+((v(:,kl)-v(:,kl-1))/dz_hl(:,kl-1)+dwdy(:,kl))**2 &
         +(2.*(ww(1:ifull,kl)-ww(1:ifull,kl-1))/dz_hl(:,kl-1))**2+shear(:,kl)/km(:,kl)
where (wt0.le.0.)
  ppt(:,kl)=0.5*(-(km(:,kl)+km(:,kl-1))*(tke(1:ifull,kl)-tke(1:ifull,kl-1))/dz_hl(:,kl-1))/dz_fl(:,kl)
elsewhere
  ppt(:,kl)=0.
end where
ppt=ppt/km(:,2:kl)

! implicit approach (split form - helps with numerical stability)
tkenew(:,2:kl)=tke(1:ifull,2:kl) ! 1st guess
ee=1./dt
do icount=1,icm
  tkenew(:,2:kl)=max(tkenew(:,2:kl),1.5E-4)
  tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)  
  aa=ce2/tkenew(:,2:kl)
  cc=-eps(1:ifull,2:kl)/dt-ce1*cm*tkenew(:,2:kl)*(pps+max(ppb(:,2:kl),0.)+max(ppt,0.))
  epsnew(:,2:kl)=(-ee+sqrt(max(ee*ee-4.*aa*cc,0.)))/(2.*aa)
  !dd=-tkenew(:,2:kl)+tke(1:ifull,2:kl)+dt*(max(cm*tkenew(:,2:kl)*tkenew(:,2:kl)/epsnew(:,2:kl),1.E-3)*(pps+max(ppb(:,2:kl),0.))-epsnew(:,2:kl)) ! error function
  dd=-tkenew(:,2:kl)+tke(1:ifull,2:kl)+dt*(km(:,2:kl)*(pps+ppb(:,2:kl))-epsnew(:,2:kl)) ! error function
  ff(:,2:kl)=-1.-dt*(epsnew(:,2:kl)/tkenew(:,2:kl) &
                 +(cc/tkenew(:,2:kl)+ce1*cm*(pps+max(ppb(:,2:kl),0.)+max(ppt,0.)))/sqrt(max(ee*ee-4.*aa*cc,0.)))
  where (abs(ff(:,2:kl)).gt.tol) ! sectant method 
    tkenew(:,2:kl)=tkenew(:,2:kl)-0.7*dd/ff(:,2:kl)
  end where
end do

tkenew(:,2:kl)=max(tkenew(:,2:kl),1.5E-4)
tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)
aa=cm34*(tkenew(:,2:kl)**1.5)/5.
epsnew(:,2:kl)=min(epsnew(:,2:kl),aa)
aa=max(aa*5./500.,1.E-6)
epsnew(:,2:kl)=max(epsnew(:,2:kl),aa)

tke(1:ifull,2:kl)=tkenew(:,2:kl)
eps(1:ifull,2:kl)=epsnew(:,2:kl)

! eps vertical mixing (done here as we skip level 1, instead of using trim)
aa(:,2)=-0.5*ce0*(km(:,2)+km(:,1))/(dz_fl(:,2)*dz_hl(:,1))
cc(:,2)=-0.5*ce0*(km(:,3)+km(:,2))/(dz_fl(:,2)*dz_hl(:,2))
bb(:,2)=1./dt-cc(:,2)-aa(:,2)
dd(:,2)=eps(1:ifull,2)/dt-aa(:,2)*eps(1:ifull,1)
do k=3,kl-1
  aa(:,k)=-0.5*ce0*(km(:,k)+km(:,k-1))/(dz_fl(:,k)*dz_hl(:,k-1))
  cc(:,k)=-0.5*ce0*(km(:,k+1)+km(:,k))/(dz_fl(:,k)*dz_hl(:,k))
  bb(:,k)=1./dt-aa(:,k)-cc(:,k)
  dd(:,k)=eps(1:ifull,k)/dt
end do
aa(:,kl)=-0.5*ce0*(km(:,kl)+km(:,kl-1))/(dz_fl(:,kl)*dz_hl(:,kl-1))
bb(:,kl)=1./dt-aa(:,kl)
dd(:,kl)=eps(1:ifull,kl)/dt
call thomas_min(epsnew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

! TKE vertical mixing (done here as we skip level 1, instead of using trim)
aa(:,2)=-0.5*(km(:,2)+km(:,1))/(dz_fl(:,2)*dz_hl(:,1))
cc(:,2)=-0.5*(km(:,3)+km(:,2))/(dz_fl(:,2)*dz_hl(:,2))
bb(:,2)=1./dt-cc(:,2)-aa(:,2)
dd(:,2)=tke(1:ifull,2)/dt-aa(:,2)*tke(1:ifull,1)
do k=3,kl-1
  aa(:,k)=-0.5*(km(:,k)+km(:,k-1))/(dz_fl(:,k)*dz_hl(:,k-1))
  cc(:,k)=-0.5*(km(:,k+1)+km(:,k))/(dz_fl(:,k)*dz_hl(:,k))
  bb(:,k)=1./dt-aa(:,k)-cc(:,k)
  dd(:,k)=tke(1:ifull,k)/dt
end do
aa(:,kl)=-0.5*(km(:,kl)+km(:,kl-1))/(dz_fl(:,kl)*dz_hl(:,kl-1))
bb(:,kl)=1./dt-aa(:,kl)
dd(:,kl)=tke(1:ifull,kl)/dt
call thomas_min(tkenew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

! Limits
!dz_ref=500.
!do k=2,kl
!  where (zz(:,k).ge.0.55*zi.and.zz(:,k).le.0.95*zi.and.wstar.gt.0.5)
!    tke(1:ifull,k)=max(tke(1:ifull,k),tke(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
!    eps(1:ifull,k)=max(eps(1:ifull,k),eps(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
!  end where
!end do
tkenew(:,2:kl)=max(tkenew(:,2:kl),1.5E-4)
tkenew(:,2:kl)=min(tkenew(:,2:kl),65.)
aa=cm34*(tkenew(:,2:kl)**1.5)/5.
epsnew(:,2:kl)=min(epsnew(:,2:kl),aa)
aa=max(aa*5./500.,1.E-6)
epsnew(:,2:kl)=max(epsnew(:,2:kl),aa)
tke(1:ifull,2:kl)=tkenew(:,2:kl)
eps(1:ifull,2:kl)=epsnew(:,2:kl)

!! Update thetav and qg for mass-flux terms (implicit split form)
!cc(:,1)=-0.25*(mflx(:,2)+mflx(:,1))/dz_fl(:,1)
!bb(:,1)=1./dt-0.25*(mflx(:,2)+mflx(:,1))/dz_fl(:,1)
!dd(:,1)=thetav(:,1)/dt-0.25*(mflx(:,2)+mflx(:,1))*(tup(:,2)+tup(:,1))/dz_fl(:,1)
!do k=2,kl-1
!  aa(:,k)=0.25*(mflx(:,k)+mflx(:,k-1))/dz_fl(:,k)
!  cc(:,k)=-0.25*(mflx(:,k+1)+mflx(:,k))/dz_fl(:,k)
!  bb(:,k)=1./dt-0.25*(mflx(:,k+1)-mflx(:,k-1))/dz_fl(:,k)
!  dd(:,k)=thetav(:,k)/dt-0.25*(mflx(:,k+1)+mflx(:,k))*(tup(:,k+1)+tup(:,k))/dz_fl(:,k) &
!                        +0.25*(mflx(:,k)+mflx(:,k-1))*(tup(:,k)+tup(:,k-1))/dz_fl(:,k)
!end do
!aa(:,kl)=0.25*(mflx(:,kl)+mflx(:,kl-1))/dz_fl(:,kl)
!bb(:,kl)=1./dt+0.25*(mflx(:,kl)+mflx(:,kl-1))/dz_fl(:,kl)
!dd(:,kl)=thetav(:,kl)/dt+0.25*(mflx(:,kl)+mflx(:,kl-1))*(tup(:,kl)+tup(:,kl-1))/dz_fl(:,kl)
!call thomas(kl,thetav,aa(:,2:kl),bb,cc(:,1:kl-1),dd)
!
!dd(:,1)=qg(:,1)/dt-0.25*(mflx(:,2)+mflx(:,1))*(qup(:,2)+qup(:,1))/dz_fl(:,1)
!do k=2,kl-1
!  dd(:,k)=qg(:,k)/dt-0.25*(mflx(:,k+1)+mflx(:,k))*(qup(:,k+1)+qup(:,k))/dz_fl(:,k) &
!                    +0.25*(mflx(:,k)+mflx(:,k-1))*(qup(:,k)+qup(:,k-1))/dz_fl(:,k)
!end do
!dd(:,kl)=qg(:,kl)/dt+0.25*(mflx(:,kl)+mflx(:,kl-1))*(tup(:,kl)+tup(:,kl-1))/dz_fl(:,kl)
!call thomas(kl,qg,aa(:,2:kl),bb,cc(:,1:kl-1),dd)

! Update diffusion coeffs at half levels
do k=1,kl-1
  kmo(:,k)=0.5*(km(:,k+1)+km(:,k))
end do
! These terms are never used
kmo(:,kl)=2.*kmo(:,kl-1)-kmo(:,kl-2) 

! Update thetav and qg due to non-local term (explicit split form)
if (mode.ne.1) then
  thetav(:,1)=thetav(:,1)-dt*0.5*(gamt(:,2)+gamt(:,1))/dz_fl(:,1)
  do k=2,kl-1
    thetav(:,k)=thetav(:,k)-dt*0.5*(gamt(:,k+1)-gamt(:,k-1))/dz_fl(:,k)
  end do
  thetav(:,kl)=thetav(:,kl)+dt*0.5*(gamt(:,kl)+gamt(:,kl-1))/dz_fl(:,kl)
end if
if (mode.ne.2) then
  qg(:,1)=qg(:,1)-dt*0.5*(gamq(:,2)+gamq(:,1))/dz_fl(:,1)
  do k=2,kl-1
    qg(:,k)=qg(:,k)-dt*0.5*(gamq(:,k+1)-gamq(:,k-1))/dz_fl(:,k)
  end do
  qg(:,kl)=qg(:,kl)+dt*0.5*(gamq(:,kl)+gamq(:,kl-1))/dz_fl(:,kl)
end if
qg=max(qg,1.E-6)
theta=thetav/(1.+0.61*qg)

tkesav=tke(1:ifull,:) ! Not needed, but for consistancy when not using CCAM
epssav=eps(1:ifull,:) ! Not needed, but for consistancy when not using CCAM

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (skips level 1 for TKE and eps)

subroutine thomas_min(out,aa,bbi,cc,ddi)

implicit none

integer k
real, dimension(ifull,3:kl), intent(in) :: aa
real, dimension(ifull,2:kl), intent(in) :: bbi,ddi
real, dimension(ifull,2:kl-1), intent(in) :: cc
real, dimension(ifull,2:kl), intent(out) :: out
real, dimension(ifull,2:kl) :: bb,dd
real, dimension(ifull) :: n

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
end subroutine thomas_min

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from TAPM)

subroutine getqsat(ilen,qsat,temp,ps)

implicit none

integer, intent(in) :: ilen
real, dimension(ilen), intent(in) :: temp,ps
real, dimension(ilen), intent(out) :: qsat
real, dimension(ilen) :: esatf

where (temp.ge.273.15)
  esatf = 610.*exp(lv/rv*(1./273.15-1./min(temp,343.)))
elsewhere
  esatf = 610.*exp(ls/rv*(1./273.15-1./max(temp,123.)))
end where
qsat = 0.622*esatf/(ps-0.378*esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE

subroutine tkeend(diag)

implicit none

integer, intent(in) :: diag

if (diag.gt.0) write(6,*) "Terminate TKE scheme"

deallocate(tke,eps)
deallocate(tkesav,epssav)
deallocate(shear,ww)
deallocate(dwdx,dwdy)

return
end subroutine tkeend

end module tkeeps
