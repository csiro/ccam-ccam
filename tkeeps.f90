
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

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: tke,eps
real, dimension(:,:), allocatable, save :: tkesav,epssav
real, dimension(:,:), allocatable, save :: shear
real, parameter :: cm  = 0.09
real, parameter :: ce0 = 0.69
real, parameter :: ce1 = 1.46
real, parameter :: ce2 = 1.83
real, parameter :: aup = 0.1
real, parameter :: b1  = 1.
real, parameter :: b2  = 2.
real, parameter :: grav = 9.80616
real, parameter :: entr = 2.E-3
real, parameter :: detr = 3.E-3
real, parameter :: cq   = 2.5
real, parameter :: vkar = 0.4
real, parameter :: a_1   = 1.0
real, parameter :: b_1   = 2.0/3.0
real, parameter :: c_1   = 5.0
real, parameter :: d_1   = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

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

tke=1.5E-4
eps=1.0E-6
tkesav=tke(1:ifull,:)
epssav=eps(1:ifull,:)
shear=0.

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

subroutine tkemix(kmo,theta,qg,u,v,zi,wt0,wq0,ps,ustar,zz,sigkap,dt,av_vmod,diag)

implicit none

integer, intent(in) :: diag
integer k,i,maxits
real, intent(in) :: dt,av_vmod
real, dimension(ifull,kl), intent(inout) :: theta,qg
real, dimension(ifull,kl), intent(in) :: u,v,zz
real, dimension(ifull,kl), intent(out) :: kmo
real, dimension(ifull), intent(inout) :: zi
real, dimension(ifull), intent(in) :: wt0,wq0,ps,ustar
real, dimension(kl), intent(in) :: sigkap
real, dimension(ifull,kl) :: km,kmsav,gam,gamhl
real, dimension(ifull,2:kl) :: jacobi11,jacobi21,jacobi12,jacobi22,det,fn,gn
real, dimension(ifull,kl) :: tkenew,epsnew
real, dimension(ifull,2:kl) :: pps,ppb,ppt
real, dimension(ifull,kl) :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(ifull,kl-1) :: dz_hl ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(ifull,2:kl) :: aa
real, dimension(ifull,2:kl) :: bb,dd
real, dimension(ifull,2:kl-1) :: cc
real, dimension(ifull) :: wstar,z_on_l,phim
real, dimension(kl-1) :: wpv_flux,tup,qup,w2up
real xp,wup,mflx,qupsat,qsat,temp,acldf,qt2,ee
real dt_s,dz_ref
logical sconv
logical, dimension(ifull,kl) :: teok
integer, parameter :: method = 1 ! 0 = explicit, 1 = implicit
real, parameter :: alpha=0.6
real, parameter :: dt_t = 100.01

if (diag.gt.0) write(6,*) "Update PBL mixing with TKE"

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
kmsav=max(cm*tkesav*tkesav/epssav,1.E-3)
km=av_vmod*km+(1.-av_vmod)*kmsav

! calculate wstar and phim (from TAPM)
wstar=0.
where (wt0.gt.0.)
  wstar=(grav*zi*wt0/theta(:,1))**(1./3.)
end where
z_on_l=-vkar*zz(:,1)*grav*wt0/(theta(:,1)*ustar*ustar*ustar)
z_on_l=min(z_on_l,10.)
where (z_on_l.lt.0.)
  phim=(1.-16.*z_on_l)**(-0.25)
else where (z_on_l.le.0.4)
  phim=1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
else where
  phim=aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar (2007)
end where

! Calculate non-local terms for theta (gamt) and qg (gamq)
zi=zz(:,kl)
gam=0.

do i=1,ifull
  if (wt0(i).gt.0.) then ! unstable
    sconv=.false.
    ee=0.5*(1./zz(i,1)+1./(zi(i)+zz(i,1))) ! dz_hl(:,0)=zz(:,1)
    tup(1)=theta(i,1)+wt0(i)/sqrt(tke(i,1))
    qup(1)=qg(i,1)+wq0(i)/sqrt(tke(i,1))
    w2up(1)=(2.*zz(i,1)*b2*grav*(tup(1)-theta(i,1))/theta(i,1))/(1.+zz(i,1)*b1*ee)
    wup=sqrt(max(w2up(1),0.))
    temp=tup(1)/sigkap(1)
    call getqsat(qupsat,temp,ps(i))
    if (qup(1).lt.qupsat) then ! dry convection
      mflx=aup*wup 
    else                       ! moist convection (boundary condition)
      sconv=.true.
      temp=theta(i,1)/sigkap(1)
      call getqsat(qsat,temp,ps(i))
      qt2=max(1.E-6,1.6*(tke(i,1)/eps(i,1))*cq*km(i,1)*((qg(i,2)-qg(i,1))/dz_hl(i,1))**2)      
      acldf=0.5+0.36*atan(1.55*(qg(i,1)-qsat)/sqrt(qt2))    
      mflx=acldf*aup*wup
    end if
    gam(i,1)=mflx*(tup(1)-theta(i,1))
    !wpv_flux(1)=-km(i,1)*(theta(i,2)-theta(i,1))/dz_hl(i,1)+km(i,1)*gam(i,1)
    do k=2,kl-1
      ee=0.5*(1./(zz(i,k-1)+dz_hl(i,k-1))+1./(max(zi(i)-zz(i,k-1),0.)+dz_hl(i,k-1)))
      tup(k)=(tup(k-1)*(1.-0.5*dz_hl(i,k-1)*ee)+0.5*dz_hl(i,k-1)*ee*(theta(i,k)+theta(i,k-1))) &
           /(1.+0.5*dz_hl(i,k-1)*ee)
      qup(k)=(qup(k-1)*(1.-0.5*dz_hl(i,k-1)*ee)+0.5*dz_hl(i,k-1)*ee*(qg(i,k)+qg(i,k-1))) &
           /(1.+0.5*dz_hl(i,k-1)*ee)
      w2up(k)=(w2up(k-1)*(1.-dz_hl(i,k-1)*b1*ee)+2.*dz_hl(i,k-1)*b2*grav*(tup(k)+tup(k-1)-theta(i,k)-theta(i,k-1)) &
            /(theta(i,k)+theta(i,k-1)))/(1.+dz_hl(i,k-1)*b1*ee)
      wup=sqrt(max(w2up(k),0.))
      if (.not.sconv) then
        temp=tup(k)/sigkap(k)
        call getqsat(qupsat,temp,ps(i))      
        if (qup(k).lt.qupsat) then ! dry convection
          mflx=aup*wup 
        else                       ! moist convection (boundary condition)
          sconv=.true.
          temp=theta(i,k)/sigkap(k)
          call getqsat(qsat,temp,ps(i))
          qt2=max(1.E-6,1.6*(tke(i,k)/eps(i,k))*cq*km(i,k)*0.25*((qg(i,k+1)-qg(i,k-1))/dz_fl(i,k))**2)        
          acldf=0.5+0.36*atan(1.55*(qg(i,k)-qsat)/sqrt(qt2))
          mflx=acldf*aup*wup
        end if
      else
        if (w2up(k).gt.0.) then     ! moist convection
          mflx=mflx*exp((entr-detr)*dz_hl(i,k-1))
        else
          mflx=0.
        end if
      end if
      gam(i,k)=mflx*(tup(k)-theta(i,k))
      !wpv_flux(k)=-km(i,k)*0.5*(theta(i,k+1)-theta(i,k-1))/dz_fl(i,k)+km(i,k)*gam(i,k)
      if (w2up(k).le.0.) then
        xp=w2up(k-1)/(w2up(k-1)-w2up(k))
        xp=min(max(xp,0.),1.)
        zi(i)=(1.-xp)*zz(i,k-1)+xp*zz(i,k)
        exit
      end if
    end do
  else                   ! stable
    wpv_flux(1)=-km(i,1)*(theta(i,2)-theta(i,1))/dz_hl(i,1) !+km(i,1)*gam(i,1)
    if (wpv_flux(1).lt.0.) then
      zi(i)=zz(i,1)
    else
      do k=2,kl-1
        wpv_flux(k)=-km(i,k)*0.5*(theta(i,k+1)-theta(i,k-1))/dz_fl(i,k) !+km(i,k)*gam(i,k)
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

! boundary condition
tke(1:ifull,1)=ustar*ustar/sqrt(cm)+0.5*wstar*wstar
eps(1:ifull,1)=ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wt0/theta(:,1)
tke(1:ifull,1)=max(tke(1:ifull,1),1.5E-4)
eps(1:ifull,1)=max(eps(1:ifull,1),1.E-6)

! Calculate source and sink terms for TKE and eps (vertical only)
do k=2,kl-1
  pps(:,k)=0.25*((u(:,k+1)-u(:,k-1))**2+(v(:,k+1)-v(:,k-1))**2)/dz_fl(:,k)**2+shear(:,k)/km(:,k)
  ppb(:,k)=-grav*(0.5*(theta(:,k+1)-theta(:,k-1))/dz_fl(:,k)-min(max(gam(:,k)/km(:,k),0.),0.002))/theta(:,k)
  where (wt0.le.0.)
    ppt(:,k)=0.5*((km(:,k+1)+km(:,k))*(tke(1:ifull,k+1)-tke(1:ifull,k))/dz_hl(:,k) &
                 -(km(:,k)+km(:,k-1))*(tke(1:ifull,k)-tke(1:ifull,k-1))/dz_hl(:,k-1))/dz_fl(:,k)
  else where
    ppt(:,k)=0.
  end where  
end do
pps(:,kl)=((u(:,kl)-u(:,kl-1))**2+(v(:,kl)-v(:,kl-1))**2)/dz_hl(:,kl-1)**2+shear(:,kl)/km(:,kl)
ppb(:,kl)=-grav*((theta(:,kl)-theta(:,kl-1))/dz_hl(:,kl-1)-min(max(gam(:,kl)/km(:,kl),0.),0.002))/theta(:,kl)
where (wt0.le.0.)
  ppt(:,kl)=0.5*(-(km(:,kl)+km(:,kl-1))*(tke(1:ifull,kl)-tke(1:ifull,kl-1))/dz_hl(:,kl-1))/dz_fl(:,kl)
else where
  ppt(:,kl)=0.
end where


select case(method)
  !******************************************************************
  case(0) ! explicit (yields correct results)

    maxits=1+int(dt/dt_t)
    dt_s=dt/real(maxits)

    do i=1,maxits

      ! TKE vertical mixing (done here as we skip level 1, instead of using trim)
      aa(:,2)=-0.5*(km(:,2)+km(:,1))/(dz_fl(:,2)*dz_hl(:,1))
      cc(:,2)=-0.5*(km(:,3)+km(:,2))/(dz_fl(:,2)*dz_hl(:,2))
      bb(:,2)=1./dt_s-cc(:,2)-aa(:,2)
      dd(:,2)=tke(1:ifull,2)/dt_s-aa(:,2)*tke(1:ifull,1)+km(:,2)*(pps(:,2)+ppb(:,2))-eps(1:ifull,2)
      do k=3,kl-1
        aa(:,k)=-0.5*(km(:,k)+km(:,k-1))/(dz_fl(:,k)*dz_hl(:,k-1))
        cc(:,k)=-0.5*(km(:,k+1)+km(:,k))/(dz_fl(:,k)*dz_hl(:,k))
        bb(:,k)=1./dt_s-aa(:,k)-cc(:,k)
        dd(:,k)=tke(1:ifull,k)/dt_s+km(:,k)*(pps(:,k)+ppb(:,k))-eps(1:ifull,k)
      end do
      aa(:,kl)=-0.5*(km(:,kl)+km(:,kl-1))/(dz_fl(:,kl)*dz_hl(:,kl-1))
      bb(:,kl)=1./dt_s-aa(:,kl)
      dd(:,kl)=tke(1:ifull,kl)/dt_s+km(:,kl)*(pps(:,kl)+ppb(:,kl))-eps(1:ifull,kl)
      call thomas_min(tkenew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

      ! eps vertical mixing (done here as we skip level 1, instead of using trim)
      aa(:,2)=-0.5*ce0*(km(:,2)+km(:,1))/(dz_fl(:,2)*dz_hl(:,1))
      cc(:,2)=-0.5*ce0*(km(:,3)+km(:,2))/(dz_fl(:,2)*dz_hl(:,2))
      bb(:,2)=1./dt_s-cc(:,2)-aa(:,2)+ce2*eps(1:ifull,2)/tke(1:ifull,2)
      dd(:,2)=eps(1:ifull,2)/dt_s-aa(:,2)*eps(1:ifull,1)+ce1*(eps(1:ifull,2)/tke(1:ifull,2)) &
                                  *(km(:,2)*(pps(:,2)+max(ppb(:,2),0.))+max(ppt(:,2),0.))
      do k=3,kl-1
        aa(:,k)=-0.5*ce0*(km(:,k)+km(:,k-1))/(dz_fl(:,k)*dz_hl(:,k-1))
        cc(:,k)=-0.5*ce0*(km(:,k+1)+km(:,k))/(dz_fl(:,k)*dz_hl(:,k))
        bb(:,k)=1./dt_s-aa(:,k)-cc(:,k)+ce2*eps(1:ifull,k)/tke(1:ifull,k)
        dd(:,k)=eps(1:ifull,k)/dt_s+ce1*(eps(1:ifull,k)/tke(1:ifull,k)) &
                                   *(km(:,k)*(pps(:,k)+max(ppb(:,k),0.))+max(ppt(:,k),0.))
      end do
      aa(:,kl)=-0.5*ce0*(km(:,kl)+km(:,kl-1))/(dz_fl(:,kl)*dz_hl(:,kl-1))
      bb(:,kl)=1./dt_s-aa(:,kl)+ce2*eps(1:ifull,kl)/tke(1:ifull,kl)
      dd(:,kl)=eps(1:ifull,kl)/dt_s+ce1*(eps(1:ifull,kl)/tke(1:ifull,kl)) &
                                   *(km(:,kl)*(pps(:,kl)+max(ppb(:,kl),0.))+max(ppt(:,kl),0.))
      call thomas_min(epsnew(:,2:kl),aa(:,3:kl),bb(:,2:kl),cc(:,2:kl-1),dd(:,2:kl))

      tke(1:ifull,2:kl)=tkenew(:,2:kl)
      eps(1:ifull,2:kl)=epsnew(:,2:kl)

      ! Limits
      !dz_ref=500.
      !do k=2,kl
      !  where (zz(:,k).ge.0.55*zi.and.zz(:,k).le.0.95*zi.and.wstar.gt.0.5)
      !    tke(1:ifull,k)=max(tke(1:ifull,k),tke(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
      !    eps(1:ifull,k)=max(eps(1:ifull,k),eps(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
      !  end where
      !end do
      tke=max(tke,1.5E-4)
      eps(:,2:kl)=min(eps(:,2:kl),(cm**0.75)*(tke(:,2:kl)**1.5)/5.)
      eps(:,2:kl)=max(eps(:,2:kl),(cm**0.75)*(tke(:,2:kl)**1.5)/500.)
      eps=max(eps,1.E-6)

    end do

  !******************************************************************
  case(1) ! implicit (currently yields incorrect results due to the split formulation)

    ! implicit approach
    maxits=10
    teok=.true.
    tkenew=tke(1:ifull,:)
    epsnew=eps(1:ifull,:)

    do i=1,maxits
 
      ! non-linear part
      fn(:,2:kl)=tke(1:ifull,2:kl)-tkenew(:,2:kl)+dt*(km(:,2:kl)*(pps(:,2:kl)+ppb(:,2:kl))-epsnew(:,2:kl))
      gn(:,2:kl)=eps(1:ifull,2:kl)-epsnew(:,2:kl)+dt*(epsnew(:,2:kl)/tkenew(:,2:kl)) &
                                                  *(ce1*km(:,2:kl)*(pps(:,2:kl)+max(ppb(:,2:kl),0.)) &
                                                   +ce1*max(ppt(:,2:kl),0.)-ce2*epsnew(:,2:kl))
                                         
      jacobi11(:,2:kl)=-1.
      jacobi21(:,2:kl)=-dt
      jacobi12(:,2:kl)=-dt*(epsnew(:,2:kl)/tkenew(:,2:kl)**2) &
                          *(ce1*km(:,2:kl)*(pps(:,2:kl)+max(ppb(:,2:kl),0.)) &
                           +ce1*max(ppt(:,2:kl),0.)-ce2*epsnew(:,2:kl))
      jacobi22(:,2:kl)=-1.+dt*(1./tkenew(:,2:kl)) &
                          *(ce1*km(:,2:kl)*(pps(:,2:kl)+max(ppb(:,2:kl),0.)) &
                           +ce1*max(ppt(:,2:kl),0.)-ce2*epsnew(:,2:kl)) &
                           +dt*(epsnew(:,2:kl)/tkenew(:,2:kl))*(-ce2)
      where (.not.teok(:,2:kl)) ! MJT suggestion - decouple equations when a bound is reached
        jacobi12(:,2:kl)=0.
        jacobi21(:,2:kl)=0.
      end where
      det(:,2:kl)=jacobi11(:,2:kl)*jacobi22(:,2:kl)-jacobi12(:,2:kl)*jacobi21(:,2:kl)
      where (abs(det(:,2:kl)).gt.0.001)
        tkenew(:,2:kl)=tkenew(:,2:kl)-alpha*(jacobi22(:,2:kl)*fn(:,2:kl)-jacobi21(:,2:kl)*gn(:,2:kl))/det(:,2:kl)
        epsnew(:,2:kl)=epsnew(:,2:kl)-alpha*(-jacobi12(:,2:kl)*fn(:,2:kl)+jacobi11(:,2:kl)*gn(:,2:kl))/det(:,2:kl)
      end where
      teok=.true.
      where (tkenew(:,2:kl).le.1.5E-4)
        tkenew(:,2:kl)=1.5E-4
        teok(:,2:kl)=.false.
      end where
      where (tkenew(:,2:kl).ge.65.)
        tkenew(:,2:kl)=65.
        teok(:,2:kl)=.false.
      end where      
      aa(:,2:kl)=(cm**0.75)*(tkenew(:,2:kl)**1.5)/5.
      where (epsnew(:,2:kl).ge.aa(:,2:kl))
        epsnew(:,2:kl)=aa(:,2:kl)
        teok(:,2:kl)=.false.
      end where
      aa(:,2:kl)=max(aa(:,2:kl)*5./500.,1.E-6)
      where (epsnew(:,2:kl).le.aa(:,2:kl))
        epsnew(:,2:kl)=aa(:,2:kl)
        teok(:,2:kl)=.false.
      end where
    end do

    tke(1:ifull,2:kl)=tkenew(:,2:kl)
    eps(1:ifull,2:kl)=epsnew(:,2:kl)

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

    ! Limits
    !dz_ref=500.
    !do k=2,kl
    !  where (zz(:,k).ge.0.55*zi.and.zz(:,k).le.0.95*zi.and.wstar.gt.0.5)
    !    tke(1:ifull,k)=max(tke(1:ifull,k),tke(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
    !    eps(1:ifull,k)=max(eps(1:ifull,k),eps(1:ifull,k-1)*(1.-0.05*dz_hl(:,k-1)/dz_ref))
    !  end where
    !end do
    tke(1:ifull,:)=max(tkenew,1.5E-4)
    eps(1:ifull,2:kl)=min(epsnew(:,2:kl),(cm**0.75)*(tkenew(:,2:kl)**1.5)/5.)
    eps(1:ifull,2:kl)=max(epsnew(:,2:kl),(cm**0.75)*(tkenew(:,2:kl)**1.5)/500.)
    eps(1:ifull,:)=max(epsnew,1.E-6)

  !******************************************************************
end select

! Update diffusion coeffs at half levels
do k=1,kl-1
  kmo(:,k)=0.5*(km(:,k+1)+km(:,k))
end do
kmo(:,kl)=2.*kmo(:,kl-1)-kmo(:,kl-2)

do k=1,kl-1
  gamhl(:,k)=min(max(0.5*(gam(:,k)+gam(:,k+1))/kmo(:,k),0.),0.002)
end do
gamhl(:,kl)=min(max(2.*gamhl(:,kl-1)-gamhl(:,kl-2),0.),0.002)

! Update theta due to non-local term (split form)
theta(:,1)=theta(:,1)-dt*(kmo(:,1)*gamhl(:,1))/dz_fl(:,1)
do k=2,kl-1
  theta(:,k)=theta(:,k)-dt*(kmo(:,k)*gamhl(:,k)-kmo(:,k-1)*gamhl(:,k-1))/dz_fl(:,k)
end do
theta(:,kl)=theta(:,kl)-dt*(-kmo(:,kl-1)*gamhl(:,kl-1))/dz_fl(:,kl)

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

subroutine getqsat(qsat,temp,ps)

implicit none

real, intent(in) :: temp,ps
real, intent(out) :: qsat
real esatf
real, parameter :: latent= 2.50E6
real, parameter :: latsub= 2.83E6
real, parameter :: rv    = 461.5

if (temp.ge.273.15) then
  esatf = 610.*exp(latent/rv*(1./273.15-1./min(max(temp,123.),343.)))
else
  esatf = 610.*exp(latsub/rv*(1./273.15-1./min(max(temp,123.),343.)))
end if
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
deallocate(shear)

return
end subroutine tkeend

end module tkeeps
