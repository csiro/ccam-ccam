
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
public tkeinit,tkemix,tkeend,tke,eps,tkesav,epssav,shear,tkedwn,epsdwn,pblhdwn,tkeotf,epsotf,pblhotf
public ww,dwdx,dwdy

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:,:), allocatable, save :: tkesav,epssav
real, dimension(:,:), allocatable, save :: ww,dwdx,dwdy
real, dimension(:,:), allocatable, save :: tkedwn,epsdwn ! These variables are for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: tkeotf,epsotf ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: pblhdwn         ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: pblhotf         ! These variables are for CCAM onthefly.f

! model constants
real, parameter :: cm  = 0.09
real, parameter :: ce0 = 0.69
real, parameter :: ce1 = 1.46
real, parameter :: ce2 = 1.83
real, parameter :: aup = 0.1
real, parameter :: b1  = 1.
real, parameter :: b2  = 2.
real, parameter :: entr = 2.E-3
real, parameter :: detr = 3.E-3
real, parameter :: cq   = 2.5

! physical constants
real, parameter :: grav = 9.80616
real, parameter :: lv = 2.5104e6
real, parameter :: epsl = 0.622
real, parameter :: rd = 287.04
real, parameter :: cp = 1004.64
real, parameter :: vkar = 0.4

! stability constants
real, parameter :: a_1   = 1.0
real, parameter :: b_1   = 2.0/3.0
real, parameter :: c_1   = 5.0
real, parameter :: d_1   = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

integer, parameter :: shallmeth = 1 ! 0 = Dry air, 1=Duynkerke, 2=Geleyn

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
allocate(shear(ifull,kl),ww(ifull,kl))
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

subroutine tkemix(kmo,theta,qg,qlg,qfg,u,v,zi,wt0,wq0,ps,ustar,zz,sig,sigkap,dt,diag)

implicit none

integer, intent(in) :: diag
integer k,i,ktop,kbot
real, intent(in) :: dt
real, dimension(ifull,kl), intent(inout) :: theta,qg
real, dimension(ifull,kl), intent(in) :: u,v,zz,qlg,qfg
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar
real, dimension(kl), intent(in) :: sigkap,sig
real, dimension(ifull,kl) :: km,gam,gamhl,ff,gg,thetav,temp
real, dimension(ifull,kl) :: tkenew,epsnew,qsat
real, dimension(ifull,2:kl) :: pps,ppb,ppt
real, dimension(ifull,kl) :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,2:kl) :: aa,bb,cc,dd
real, dimension(ifull) :: wstar,z_on_l,phim,wtv0
real, dimension(kl-1) :: wpv_flux,tup,qup,w2up
real xp,wup,mflx,qupsat(1),acldf,qt2,ee
real dz_ref,cm34
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

! Calculate diffusion coeffs
km=max(cm*tke(1:ifull,:)*tke(1:ifull,:)/eps(1:ifull,:),1.E-3)

! calculate saturated mixing ratio
do k=1,kl
  temp(:,k)=theta(:,k)/sigkap(k)
  call getqsat(ifull,qsat(:,k),temp(:,k),ps*sig(k))
end do

! Calculate non-local terms for theta_v
zi=zz(:,kl)
gam=0.

do i=1,ifull
  if (wtv0(i).gt.0.) then ! unstable
    sconv=.false.
    ee=0.5*(1./zz(i,1)+1./(zi(i)+zz(i,1))) ! dz_hl(:,0)=zz(:,1)
    tup(1)=thetav(i,1)+wtv0(i)/sqrt(tke(i,1))
    qup(1)=qg(i,1)+wq0(i)/sqrt(tke(i,1))
    w2up(1)=(2.*zz(i,1)*b2*grav*(tup(1)-thetav(i,1))/thetav(i,1))/(1.+zz(i,1)*b1*ee)
    wup=sqrt(max(w2up(1),0.))
    aa(1,2)=ps(i)*sig(1)
    call getqsat(1,qupsat(1),temp(i,1),aa(1,2))
    if (qup(1).lt.qupsat(1)) then ! dry convection
      mflx=aup*wup 
    else                       ! moist convection (boundary condition)
      sconv=.true.
      qt2=max(1.E-6,1.6*(tke(i,1)/eps(i,1))*cq*km(i,1)*((qg(i,2)-qg(i,1))/dz_hl(i,1))**2)      
      acldf=0.5+0.36*atan(1.55*(qg(i,1)-qsat(i,k))/sqrt(qt2))    
      mflx=acldf*aup*wup
    end if
    gam(i,1)=mflx*(tup(1)-thetav(i,1))
    !wpv_flux(1)=-km(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1)+gam(i,1)
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
        call getqsat(1,qupsat(1),temp(i,1),aa(1,2))      
        if (qup(k).lt.qupsat(1)) then ! dry convection
          mflx=aup*wup
        else                       ! moist convection (boundary condition)
          sconv=.true.
          qt2=max(1.E-6,1.6*(tke(i,k)/eps(i,k))*cq*km(i,k)*0.25*((qg(i,k+1)-qg(i,k-1))/dz_fl(i,k))**2)        
          acldf=0.5+0.36*atan(1.55*(qg(i,k)-qsat(i,k))/sqrt(qt2))
          mflx=acldf*aup*wup
        end if
      else
        if (w2up(k).gt.0.) then     ! moist convection
          mflx=mflx*exp((entr-detr)*dz_hl(i,k-1))
        else
          mflx=0.
        end if
      end if
      gam(i,k)=mflx*(tup(k)-thetav(i,k))
      !wpv_flux(k)=-km(i,k)*0.5*(thetav(i,k+1)-thetav(i,k-1))/dz_fl(i,k)+gam(i,k)
      if (w2up(k).le.0.) then
        xp=w2up(k-1)/(w2up(k-1)-w2up(k))
        xp=min(max(xp,0.),1.)
        zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
        exit
      end if
    end do
  else                   ! stable
    wpv_flux(1)=-km(i,1)*(thetav(i,2)-thetav(i,1))/dz_hl(i,1) !+gam(i,1)
    if (wpv_flux(1).lt.0.) then
      zi(i)=zz(i,1)
    else
      do k=2,kl-1
        wpv_flux(k)=-km(i,k)*0.5*(thetav(i,k+1)-thetav(i,k-1))/dz_fl(i,k) !+gam(i,k)
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

gam=max(gam,0.)

! calculate wstar and phim (from TAPM)
wstar=0.
where (wtv0.gt.0.)
  wstar=(grav*zi*wtv0/thetav(:,1))**(1./3.)
end where
z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*ustar*ustar*ustar)
z_on_l=min(z_on_l,10.)
where (z_on_l.lt.0.)
  phim=(1.-16.*z_on_l)**(-0.25)
elsewhere (z_on_l.le.0.4)
  phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
elsewhere
  phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar (2007)
end where

! boundary condition
tke(1:ifull,1)=ustar*ustar/sqrt(cm)+0.5*wstar*wstar
eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
tke(1:ifull,1)=max(tke(1:ifull,1),1.5E-4)
tke(1:ifull,1)=min(tke(1:ifull,1),65.)
aa(:,2)=cm34*(tke(1:ifull,1)**1.5)/5.
eps(1:ifull,1)=min(eps(1:ifull,1),aa(:,2))
aa(:,2)=max(aa(:,2)*5./500.,1.E-6)
eps(1:ifull,1)=max(eps(1:ifull,1),aa(:,2))

! Calculate buoyancy terms
select case(shallmeth)

  case(1) ! Saturated conditions based on WRF formulation
    do k=1,kl
      ff(:,k)=theta(:,k)*(1.+lv*qsat(:,k)/(cp*temp(:,k)))
      gg(:,k)=(1.+1.61*epsl*lv*qg(:,k)/(rd*temp(:,k))) &
             /(1.+epsl*lv*lv*qg(:,k)/(cp*1.61*rd*temp(:,k)*temp(:,k)))/theta(:,k)  
    end do
    do k=2,kl-1
      where (qg(:,k).lt.qsat(:,k))
        ppb(:,k)=-grav*(0.5*(thetav(:,k+1)-thetav(:,k-1))/dz_fl(:,k)-min(max(gam(:,k)/km(:,k),0.),0.002))/thetav(:,k) &
                 +grav*0.5*(qlg(:,k+1)+qfg(:,k+1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
      elsewhere
        ppb(:,k)=-grav*(gg(:,k)*0.5*(ff(:,k+1)-ff(:,k-1))/dz_fl(:,k)-min(max(gam(:,k)/km(:,k),0.),0.002)/thetav(:,k)) &
                 +grav*0.5*(qg(:,k+1)+qlg(:,k+1)+qfg(:,k+1)-qg(:,k-1)-qlg(:,k-1)-qfg(:,k-1))/dz_fl(:,k)
      end where
    end do
    where (qg(:,kl).lt.qsat(:,kl))
      ppb(:,kl)=-grav*((thetav(:,kl)-thetav(:,kl-1))/dz_hl(:,kl-1)-min(max(gam(:,kl)/km(:,kl),0.),0.002))/thetav(:,kl) &
                +grav*(qlg(:,kl)+qfg(:,kl)-qlg(:,kl-1)-qfg(:,kl-1))/dz_hl(:,kl-1)
    elsewhere
      ppb(:,kl)=-grav*(gg(:,kl)*(ff(:,kl)-ff(:,kl-1))/dz_hl(:,kl-1)-min(max(gam(:,kl)/km(:,kl),0.),0.002)/thetav(:,kl)) &
                +grav*(qg(:,kl)+qlg(:,kl)+qfg(:,kl)-qg(:,kl-1)-qlg(:,kl-1)-qfg(:,kl-1))/dz_hl(:,kl-1)
    end where
      
  case(2) ! Geleyn
    ktop=2
    do while(sig(ktop+1).gt.0.75)
      ktop=ktop+1
    end do
    kbot=2
    do while(sig(kbot).gt.0.99)
      kbot=kbot+1
    end do
    do k=kbot,ktop
      ppb(:,k)=-grav*(0.5*(theta(:,k+1)-theta(:,k-1))/(dz_fl(:,k)*theta(:,k))-min(max(gam(:,k)/km(:,k),0.),0.002)/thetav(:,k)) &
               -grav*(lv/cp)*min(0.,0.5*(qg(:,k+1)-qsat(:,k+1)-qg(:,k-1)+qsat(:,k-1))/dz_fl(:,k))
    end do
  
  case DEFAULT ! Hurley
    do k=2,kl-1
      ppb(:,k)=-grav*(0.5*(thetav(:,k+1)-thetav(:,k-1))/dz_fl(:,k)-min(max(gam(:,k)/km(:,k),0.),0.002))/thetav(:,k)
    end do
    ppb(:,kl)=-grav*((thetav(:,kl)-thetav(:,kl-1))/dz_hl(:,kl-1)-min(max(gam(:,kl)/km(:,kl),0.),0.002))/thetav(:,kl)
    
end select

! Calculate shear and transport terms
do k=2,kl-1
  pps(:,k)=((u(:,k+1)-u(:,k-1))/dz_fl(:,k)+dwdx(:,k))**2+((v(:,k+1)-v(:,k-1))/dz_fl(:,k)+dwdy(:,k))**2 &
          +(2.*(ww(:,k+1)-ww(:,k-1))/dz_fl(:,k))**2+shear(:,k)/km(:,k)
  where (wt0.le.0.)
    ppt(:,k)=0.5*((km(:,k+1)+km(:,k))*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k) &
                 -(km(:,k)+km(:,k-1))*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1))/dz_fl(:,k)
  elsewhere
    ppt(:,k)=0.
  end where  
end do
pps(:,kl)=((u(:,kl)-u(:,kl-1))/dz_hl(:,kl-1)+dwdx(:,kl))**2+((v(:,kl)-v(:,kl-1))/dz_hl(:,kl-1)+dwdy(:,kl))**2 &
         +(2.*(ww(:,kl)-ww(:,kl-1))/dz_hl(:,kl-1))**2+shear(:,kl)/km(:,kl)
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
  cc=-eps(1:ifull,2:kl)/dt-ce1*cm*tkenew(:,2:kl)*(pps+max(ppb,0.)+max(ppt,0.))
  epsnew(:,2:kl)=(-ee+sqrt(max(ee*ee-4.*aa*cc,0.)))/(2.*aa)
  !dd=-tkenew(:,2:kl)+tke(1:ifull,2:kl)+dt*(max(cm*tkenew(:,2:kl)*tkenew(:,2:kl)/epsnew(:,2:kl),1.E-3)*(pps+ppb)-epsnew(:,2:kl)) ! error function
  dd=-tkenew(:,2:kl)+tke(1:ifull,2:kl)+dt*(km(:,2:kl)*(pps+ppb)-epsnew(:,2:kl)) ! error function
  ff(:,2:kl)=-1.-dt*(epsnew(:,2:kl)/tkenew(:,2:kl) &
                 +(cc/tkenew(:,2:kl)+ce1*cm*(pps+max(ppb,0.)+max(ppt,0.)))/sqrt(max(ee*ee-4.*aa*cc,0.)))
  where (abs(ff(:,2:kl)).gt.tol)  
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

! Update diffusion coeffs at half levels
do k=1,kl-1
  kmo(:,k)=0.5*(km(:,k+1)+km(:,k))
  gamhl(:,k)=0.5*(gam(:,k)+gam(:,k+1))/kmo(:,k)
end do
! These terms are never used
kmo(:,kl)=2.*kmo(:,kl-1)-kmo(:,kl-2) 
gamhl(:,kl)=2.*gamhl(:,kl-1)-gamhl(:,kl-2)

gamhl=min(max(gamhl,0.),0.002)

! Update thetav due to non-local term (split form)
thetav(:,1)=thetav(:,1)-dt*(kmo(:,1)*gamhl(:,1))/dz_fl(:,1)
do k=2,kl-1
  thetav(:,k)=thetav(:,k)-dt*(kmo(:,k)*gamhl(:,k)-kmo(:,k-1)*gamhl(:,k-1))/dz_fl(:,k)
end do
thetav(:,kl)=thetav(:,kl)-dt*(-kmo(:,kl-1)*gamhl(:,kl-1))/dz_fl(:,kl)

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
real, parameter :: latent= 2.50E6
real, parameter :: latsub= 2.83E6
real, parameter :: rv    = 461.5

where (temp.ge.273.15)
  esatf = 610.*exp(latent/rv*(1./273.15-1./min(max(temp,123.),343.)))
elsewhere
  esatf = 610.*exp(latsub/rv*(1./273.15-1./min(max(temp,123.),343.)))
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
