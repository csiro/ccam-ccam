
! This is a 1D, mixed layer ocean model for regioal climate simulations based on Large, et al (1994)
! (i.e., adapted from the GFDL code)

! Specifically, this version can assimilate SSTs from GCMs, using a convolution based digital filter
! which avoids problems with complex land-sea boundary conditions

module mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,wlev,depth

type tdata
  real temp,sal,u,v
end type tdata
type tprog2
  real mixdepth,bf
  integer mixind
end type tprog2
type tatm
  real sg,rg,fg,eg,rnd,taux,tauy,f
end type tatm
type tout
  real sst
end type tout
type tdiag2
  real b0,ustar,wu0,wv0,wt0,ws0
end type tdiag2
type tdiag3
  real rho,nsq,rad,alpha,beta
end type tdiag3

! parameters
integer, parameter :: wlev = 10
integer, parameter :: iqx = 500
! model arrays
integer, save :: wfull
integer, dimension(:), allocatable, save :: wgrid
real, dimension(:,:), allocatable, save :: depth,dz
real, dimension(:,:), allocatable, save :: depth_hl
real, dimension(:,:), allocatable, save :: dz_hl
type(tdata), dimension(:,:), allocatable, save :: water
type(tprog2), dimension(:), allocatable, save :: pg

! mode
integer, parameter :: incradbf  = 1 ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 0 ! include shortwave in non-local term

! max depth
real, parameter :: mxd    = 1000.   ! Max depth
real, parameter :: mindep = 1.      ! Thickness of first layer

! model parameters
real, parameter :: ric     = 0.3    ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1

! radiation parameters
real, parameter :: r_1 = 0.58
real, parameter :: r_2 = 1.-r_1
real, parameter :: mu_1 = 0.35
real, parameter :: mu_2 = 23.

! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: grav=9.80              ! graviational constant (m/s^2)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialise MLO


subroutine mloinit(ifull,depin,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iq,iqw,ii
real, dimension(ifull), intent(in) :: depin
real smxd,smnd

if (diag.ge.1) write(6,*) "Initialising MLO"

wfull=count(depin.gt.mindep)
if (wfull.eq.0) return

allocate(water(wfull,wlev),wgrid(wfull),pg(wfull))
allocate(depth(wfull,wlev),dz(wfull,wlev))
allocate(depth_hl(wfull,wlev+1))
allocate(dz_hl(wfull,2:wlev))

wgrid=0
iqw=0
do iq=1,ifull
  if (depin(iq).gt.mindep) then
    iqw=iqw+1
    wgrid(iqw)=iq
  end if
end do

if (iqw.ne.wfull) then
  write(6,*) "ERROR: Internal memory error in MLO"
  stop
end if

water%temp=288. ! K
water%sal=35.   ! PSU
water%u=0.      ! m/s
water%v=0.      ! m/s

pg%mixdepth=100. ! m
pg%mixind=wlev-1

!depth = (/ 0.5, 4.5, 13., 25., 41., 61., 85., 113., 145., 181., 221., 265., 313.,
!           365., 421., 481., 545., 613., 685., 761., 841, 925., 1013.  /)  ! 2.*(x-0.5)^2

smxd=maxval(depin(wgrid))
smnd=minval(depin(wgrid))

do iqw=1,wfull
  call vgrid(depin(wgrid(iqw)),depth(iqw,:),depth_hl(iqw,:))
  if (smxd.eq.depin(wgrid(iqw))) then
    print *,"MLO max depth ",depth(iqw,:)
    smxd=smxd+10.
  end if
  if (smnd.eq.depin(wgrid(iqw))) then
    print *,"MLO min depth ",depth(iqw,:)
    smnd=smnd-10.
  end if
end do
do ii=1,wlev
  dz(:,ii)=depth_hl(:,ii+1)-depth_hl(:,ii)
end do
do ii=2,wlev
  dz_hl(:,ii)=depth(:,ii)-depth(:,ii-1)
end do

return
end subroutine mloinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define vertical grid
!

subroutine vgrid(depin,depthout,depth_hlout)

implicit none

integer ii
real, intent(in) :: depin
real, dimension(wlev), intent(out) :: depthout
real, dimension(wlev+1), intent(out) :: depth_hlout
real dd,x,al,bt

dd=min(mxd,max(mindep*real(wlev),depin))
x=real(wlev)
al=(mindep*x-dd)/(x-x*x*x)
bt=(mindep*x*x*x-dd)/(x*x*x-x)
do ii=1,wlev+1
  x=real(ii)-1.
  depth_hlout(ii)=al*x*x*x+bt*x ! ii is for half level ii-0.5
end do
do ii=1,wlev
  depthout(ii)=0.5*(depth_hlout(ii)+depth_hlout(ii+1))
end do

return
end subroutine vgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate MLO arrays

subroutine mloend

implicit none

if (wfull.eq.0) return

deallocate(water,wgrid,pg)
deallocate(depth,dz,depth_hl,dz_hl)

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(ifull,datain,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,wlev,4), intent(in) :: datain

if (wfull.eq.0) return

do ii=1,wlev
  water(:,ii)%temp=datain(wgrid,ii,1)
  water(:,ii)%sal=datain(wgrid,ii,2)
  water(:,ii)%u=datain(wgrid,ii,3)
  water(:,ii)%v=datain(wgrid,ii,4)
end do

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(ifull,dataout,depout,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,wlev,4), intent(out) :: dataout
real, dimension(ifull), intent(out) :: depout

if (wfull.eq.0) return

do ii=1,wlev
  dataout(wgrid,ii,1)=water(:,ii)%temp
  dataout(wgrid,ii,2)=water(:,ii)%sal
  dataout(wgrid,ii,3)=water(:,ii)%u
  dataout(wgrid,ii,4)=water(:,ii)%v
end do
depout=0.
depout(wgrid)=depth_hl(:,wlev+1)

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mloimport(ifull,sst,ilev,diag)

implicit none

integer, intent(in) :: ifull,ilev,diag
real, dimension(ifull), intent(in) :: sst

if (wfull.eq.0) return
water(:,ilev)%temp=sst(wgrid)

return
end subroutine mloimport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mloexport(ifull,sst,ilev,diag)

implicit none

integer, intent(in) :: ifull,ilev,diag
real, dimension(ifull), intent(inout) :: sst

if (wfull.eq.0) return
sst(wgrid)=water(:,ilev)%temp

return
end subroutine mloexport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data
!

subroutine mloregrid(ifull,depin,mlodat)

implicit none

integer, intent(in) :: ifull
integer iqw,ii,pos(1)
real, dimension(ifull), intent(in) :: depin
real, dimension(ifull,wlev,4), intent(inout) :: mlodat
real, dimension(ifull,wlev,4) :: newdat
real, dimension(wlev) :: dpin
real, dimension(wlev+1) :: dp_hlin
real x

if (wfull.eq.0) return
if (all(depin(wgrid).eq.depth(:,wlev))) return

newdat=mlodat

do iqw=1,wfull
  call vgrid(depin(wgrid(iqw)),dpin,dp_hlin)
  do ii=1,wlev
    if (depth(iqw,ii).ge.dpin(wlev)) then
      newdat(wgrid(iqw),ii,:)=mlodat(wgrid(iqw),wlev,:)
    else if (depth(iqw,ii).le.dpin(1)) then
      newdat(wgrid(iqw),ii,:)=mlodat(wgrid(iqw),1,:)
    else
      pos=maxloc(dpin,dpin.lt.depth(iqw,ii))
      pos(1)=max(1,min(wlev-1,pos(1)))
      x=(depth(iqw,ii)-dpin(pos(1)))/(dpin(pos(1)+1)-dpin(pos(1)))
      x=max(0.,min(1.,x))
      newdat(wgrid(iqw),ii,:)=mlodat(wgrid(iqw),pos(1)+1,:)*x+mlodat(wgrid(iqw),pos(1),:)*(1.-x)
    end if
  end do  
end do

mlodat=newdat

return
end subroutine mloregrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

subroutine mloeval(ifull,sst,dt,fg,eg,sg,rg,precp,taux,tauy,f,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: dt
real, dimension(ifull), intent(in) :: sg,rg,fg,eg,taux,tauy,precp,f
real, dimension(ifull), intent(inout) :: sst
type(tatm), dimension(wfull) :: atm
type(tout), dimension(wfull) :: uo

if (wfull.eq.0) return

atm%sg=sg(wgrid)
atm%rg=rg(wgrid)
atm%fg=fg(wgrid)
atm%eg=eg(wgrid)
atm%rnd=precp(wgrid)
atm%taux=taux(wgrid)
atm%tauy=tauy(wgrid)
atm%f=f(wgrid)

call mlocalc(uo,dt,atm,diag)

sst(wgrid)=uo%sst

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs

subroutine mlocalc(uo,dt,atm,diag)

implicit none

integer, intent(in) :: diag
integer ii,iqw,pos(1)
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
real, dimension(wfull,2:wlev) :: aa
real, dimension(wfull,wlev) :: bb,dd
real, dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull) :: xp,xm,dumt0
type(tout), dimension(wfull), intent(out) :: uo
type(tatm), dimension(wfull), intent(in) :: atm
type(tdiag2), dimension(wfull) :: dg2
type(tdiag3), dimension(wfull,wlev) :: dg3
type(tdata), dimension(wfull,wlev) :: new

call getrho(dg2,dg3,atm)           ! calculate rho and bf.  Also calculate boundary conditions.
call getmixdepth(dg2,dg3,atm)      ! solve for mixed layer depth
call getstab(km,ks,gammas,dg2,dg3) ! solve for stability functions and non-local term

! TEMPERATURE
! use +ve for gamma terms as depth is down
if (incradgam.gt.0) then
  do iqw=1,wfull
    dumt0(iqw)=dg2(iqw)%wt0-sum(dg3(iqw,1:pg(iqw)%mixind)%rad)
  end do
else
  dumt0=dg2%wt0
end if
rhs(:,1)=(0.25/dz(:,1))*(ks(:,2)+ks(:,1))*(gammas(:,2)+gammas(:,1))*dumt0
do ii=2,wlev-1
  where (ii.le.pg%mixind)
    rhs(:,ii)=(0.25/dz(:,ii))*((ks(:,ii+1)+ks(:,ii))*(gammas(:,ii+1)+gammas(:,ii)) &
                          -(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dumt0 ! non-local
  elsewhere (ii.eq.pg%mixind+1)
    rhs(:,ii)=(0.25/dz(:,ii))*(-(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dumt0
  elsewhere
    rhs(:,ii)=0.
  end where
end do
where (wlev.eq.pg%mixind+1)
  rhs(:,wlev)=(0.25/dz(:,wlev))*(-(ks(:,wlev)+ks(:,wlev-1))*(gammas(:,wlev)+gammas(:,wlev-1)))*dumt0
elsewhere
  rhs(:,wlev)=0.
end where
rhs(:,:)=rhs(:,:)+dg3%rad/dz ! shortwave
bb(:,1)=1./dt+0.5*(ks(:,1)+ks(:,2))/(dz_hl(:,2)*dz(:,1))
cc(:,1)=-0.5*(ks(:,1)+ks(:,2))/(dz_hl(:,2)*dz(:,1))
! use -ve for BC as depth is down
dd(:,1)=water(:,1)%temp/dt+rhs(:,1)-dg2%wt0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-0.5*(ks(:,ii)+ks(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  bb(:,ii)=1./dt+0.5*((ks(:,ii+1)+ks(:,ii))/dz_hl(:,ii+1)+(ks(:,ii)+ks(:,ii-1))/dz_hl(:,ii))/dz(:,ii)
  cc(:,ii)=-0.5*(ks(:,ii+1)+ks(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  dd(:,ii)=water(:,ii)%temp/dt+rhs(:,ii)
end do
aa(:,wlev)=-0.5*(ks(:,wlev)+ks(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1./dt+0.5*(ks(:,wlev)+ks(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
dd(:,wlev)=water(:,wlev)%temp/dt+rhs(:,wlev)
call thomas(new%temp,aa,bb,cc,dd)

! SALINITY
rhs(:,1)=(0.25/dz(:,1))*(ks(:,2)+ks(:,1))*(gammas(:,2)+gammas(:,1))*dg2%ws0
do ii=2,wlev-1
  where(ii.le.pg%mixind)
    rhs(:,ii)=(0.25/dz(:,ii))*((ks(:,ii+1)+ks(:,ii))*(gammas(:,ii+1)+gammas(:,ii)) &
                          -(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dg2%ws0 ! non-local
  elsewhere(ii.eq.pg%mixind+1)
    rhs(:,ii)=(0.25/dz(:,ii))*(-(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dg2%ws0  
  elsewhere
    rhs(:,ii)=0.
  end where
end do
where(wlev.eq.pg%mixind+1)
  rhs(:,wlev)=(0.25/dz(:,wlev))*(-(ks(:,wlev)+ks(:,wlev-1))*(gammas(:,wlev)+gammas(:,wlev-1)))*dg2%ws0
elsewhere
  rhs(:,wlev)=0.
end where
dd(:,1)=water(:,1)%sal/dt+rhs(:,1)-dg2%ws0/dz(:,1)
do ii=2,wlev
  dd(:,ii)=water(:,ii)%sal/dt+rhs(:,ii)
end do
call thomas(new%sal,aa,bb,cc,dd)

! split U diffusion term
bb(:,1)=1./dt+0.5*(km(:,1)+km(:,2))/(dz_hl(:,2)*dz(:,1))
cc(:,1)=-0.5*(km(:,1)+km(:,2))/(dz_hl(:,2)*dz(:,1))
dd(:,1)=water(:,1)%u/dt-dg2%wu0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-0.5*(km(:,ii)+km(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  bb(:,ii)=1./dt+0.5*((km(:,ii+1)+km(:,ii))/dz_hl(:,ii+1)+(km(:,ii)+km(:,ii-1))/dz_hl(:,ii))/dz(:,ii)
  cc(:,ii)=-0.5*(km(:,ii+1)+km(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  dd(:,ii)=water(:,ii)%u/dt
end do
aa(:,wlev)=-0.5*(km(:,wlev)+km(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1./dt+0.5*(km(:,wlev)+km(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
dd(:,wlev)=water(:,wlev)%u/dt
call thomas(new%u,aa,bb,cc,dd)

! split V diffusion term
dd(:,1)=water(:,1)%v/dt-dg2%wv0/dz(:,1)
do ii=2,wlev
  dd(:,ii)=water(:,ii)%v/dt
end do
call thomas(new%v,aa,bb,cc,dd)

water%u=new%u
water%v=new%v

! Split U and V coriolis terms
xp=1.+(0.5*dt*atm%f)**2
xm=1.-(0.5*dt*atm%f)**2
do ii=1,wlev
  new(:,ii)%u=(water(:,ii)%u*xm+water(:,ii)%v*dt*atm%f)/xp
  new(:,ii)%v=(water(:,ii)%v*xm-water(:,ii)%u*dt*atm%f)/xp
end do

water=new

!--------------------------------------------------------------------
uo%sst=water(:,1)%temp

return
end subroutine mlocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tri-diagonal solver

subroutine thomas(out,aa,bbi,cc,ddi)

implicit none

integer ii
real, dimension(wfull,2:wlev), intent(in) :: aa
real, dimension(wfull,wlev), intent(in) :: bbi,ddi
real, dimension(wfull,wlev-1), intent(in) :: cc
real, dimension(wfull,wlev), intent(out) :: out
real, dimension(wfull,wlev) :: bb,dd
real, dimension(wfull) :: n

bb=bbi
dd=ddi

do ii=2,wlev
  n=aa(:,ii)/bb(:,ii-1)
  bb(:,ii)=bb(:,ii)-n*cc(:,ii-1)
  dd(:,ii)=dd(:,ii)-n*dd(:,ii-1)
end do
out(:,wlev)=dd(:,wlev)/bb(:,wlev)
do ii=wlev-1,1,-1
  out(:,ii)=(dd(:,ii)-cc(:,ii)*out(:,ii+1))/bb(:,ii)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for stability functions
! GFDL use a look-up table to speed-up the code...

subroutine getstab(km,ks,gammas,dg2,dg3)

implicit none

integer ii,iqw
real, dimension(wfull,wlev), intent(out) :: km,ks,gammas
real, dimension(wfull,wlev) :: num,nus,wm,ws,ri
real, dimension(wfull) :: sigma
real, dimension(wfull) :: a2m,a3m,a2s,a3s
real numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real xp,cg
type(tdiag2), dimension(wfull), intent(in) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(in) :: dg3
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.E-4
real, parameter :: nusw = 0.1E-4

do ii=2,wlev-1
  ri(:,ii)=dg3(:,ii)%nsq*(2.*dz(:,ii))**2 &
           /max((water(:,ii-1)%u-water(:,ii+1)%u)**2 &
               +(water(:,ii-1)%v-water(:,ii+1)%v)**2,1.E-10)
end do
ri(:,1)=2.*ri(:,2)-ri(:,3)
ri(:,wlev)=2.*ri(:,wlev-1)-ri(:,wlev-2)

! diffusion ---------------------------------------------------------
do ii=1,wlev

  ! s - shear
  where (ri(:,ii).lt.0.)
    num(:,ii)=nu0
  elsewhere (ri(:,ii).lt.ri0)
    num(:,ii)=nu0*(1.-(ri(:,ii)/ri0)**2)**3
  elsewhere
    num(:,ii)=0.
  endwhere
  nus(:,ii)=num(:,ii)
  
  ! d - double-diffusive
  ! double-diffusive mixing is neglected for now
  
end do
! w - internal wave
num=num+numw
nus=nus+nusw
!--------------------------------------------------------------------

! stability ---------------------------------------------------------
do ii=1,wlev
  call getwx(wm(:,ii),ws(:,ii),depth(:,ii),pg%bf,dg2%ustar,pg%mixdepth)
end do
!--------------------------------------------------------------------

! calculate G profile -----------------------------------------------
! -ve as z is down
do iqw=1,wfull
  xp=(pg(iqw)%mixdepth-depth(iqw,pg(iqw)%mixind))/(depth(iqw,pg(iqw)%mixind+1)-depth(iqw,pg(iqw)%mixind))
  xp=max(0.,min(1.,xp))
  numh=(1.-xp)*num(iqw,pg(iqw)%mixind)+xp*num(iqw,pg(iqw)%mixind+1)
  wm1=(1.-xp)*wm(iqw,pg(iqw)%mixind)+xp*wm(iqw,pg(iqw)%mixind+1)
  dnumhdz=-(num(iqw,pg(iqw)%mixind+1)-num(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dwm1ds=-pg(iqw)%mixdepth*(wm(iqw,pg(iqw)%mixind+1)-wm(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  nush=(1.-xp)*nus(iqw,pg(iqw)%mixind)+xp*nus(iqw,pg(iqw)%mixind+1)
  ws1=(1.-xp)*ws(iqw,pg(iqw)%mixind)+xp*ws(iqw,pg(iqw)%mixind+1)
  dnushdz=-(nus(iqw,pg(iqw)%mixind+1)-nus(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dws1ds=-pg(iqw)%mixdepth*(ws(iqw,pg(iqw)%mixind+1)-ws(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)

  wm1=max(wm1,1.E-10)
  ws1=max(ws1,1.E-10)
  
  g1m=numh/(pg(iqw)%mixdepth*wm1)
  dg1mds=-dnumhdz/wm1-numh*dwm1ds/(pg(iqw)%mixdepth*wm1*wm1)
  g1s=nush/(pg(iqw)%mixdepth*ws1)
  dg1sds=-dnushdz/ws1-nush*dws1ds/(pg(iqw)%mixdepth*ws1*ws1)

  a2m(iqw)=-2.+3.*g1m-dg1mds
  a3m(iqw)=1.-2.*g1m+dg1mds
  a2s(iqw)=-2.+3.*g1s-dg1sds
  a3s(iqw)=1.-2.*g1s+dg1sds
end do
!--------------------------------------------------------------------
! combine
do ii=1,wlev
  where (ii.le.pg%mixind)
    sigma=depth(:,ii)/pg%mixdepth
    km(:,ii)=pg%mixdepth*wm(:,ii)*(sigma+a2m*sigma**2+a3m*sigma**3)
    ks(:,ii)=pg%mixdepth*ws(:,ii)*(sigma+a2s*sigma**2+a3s*sigma**3)
  elsewhere
    km(:,ii)=num(:,ii)
    ks(:,ii)=nus(:,ii)
  end where
end do

! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.)
do ii=1,wlev
  where (pg%bf.lt.0.) ! unstable
    gammas(:,ii)=cg/(ws(:,ii)*pg%mixdepth)
  elsewhere            ! stable
    gammas(:,ii)=0.
  end where
end do

! enhancement is neglected for now

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixing layer depth

subroutine getmixdepth(dg2,dg3,atm)

implicit none

integer ii,iqw
real vtc,dvsq,vtsq,xp
real, dimension(wfull,wlev) :: ws,wm
real, dimension(wfull) :: dumbf,l,he
real, dimension(wlev) :: rib
type(tdiag2), dimension(wfull), intent(in) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(in) :: dg3
type(tatm), dimension(wfull), intent(in) :: atm

vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar**2*ric)

if (incradbf.gt.0) then
  do ii=1,wlev
    dumbf=dg2%b0+grav*sum(dg3(:,1:ii)%alpha*dg3(:,1:ii)%rad,2)
    call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dumbf,dg2%ustar,depth(:,ii))
  end do
else
  dumbf=dg2%b0
  do ii=1,wlev
    call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dumbf,dg2%ustar,depth(:,ii))
  end do
end if

pg%mixind=wlev-1
pg%mixdepth=depth(:,wlev)
do iqw=1,wfull
  ! should be averaged over 0 < sigma < epsilon instead of 1st level
  rib(1)=0.
  do ii=2,wlev
    vtsq=depth(iqw,ii)*ws(iqw,ii)*sqrt(abs(dg3(iqw,ii)%nsq))*vtc
    dvsq=(water(iqw,1)%u-water(iqw,ii)%u)**2+(water(iqw,1)%v-water(iqw,ii)%v)**2
    rib(ii)=(depth(iqw,ii)-depth(iqw,1))*(1.-dg3(iqw,1)%rho/dg3(iqw,ii)%rho)/max(dvsq+vtsq,1.E-10)
  
    if (rib(ii).gt.ric) then
      pg(iqw)%mixind=ii-1
      xp=min(max((ric-rib(ii-1))/max(rib(ii)-rib(ii-1),1.E-10),0.),1.)
      pg(iqw)%mixdepth=(1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii)
      exit
    end if
  
  end do 

end do

! calculate buoyancy forcing
call getbf(dg2,dg3)

! impose limits for stable conditions
l=dg2%ustar*dg2%ustar*dg2%ustar/(vkar*pg%bf)
he=0.7*dg2%ustar/abs(atm%f)
where(pg%bf.gt.0.)
  pg%mixdepth=min(pg%mixdepth,l)
  pg%mixdepth=min(pg%mixdepth,he)
end where
pg%mixdepth=max(pg%mixdepth,depth(:,1))

! recalculate index for mixdepth
do iqw=1,wfull
  do ii=2,wlev
    if (depth(iqw,ii).gt.pg(iqw)%mixdepth) then
      pg(iqw)%mixind=ii-1
      exit
    end if
  end do
end do

! recalculate buoyancy forcing
call getbf(dg2,dg3)

return
end subroutine getmixdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate bouyancy forcing
!

subroutine getbf(dg2,dg3)

implicit none

integer iqw
type(tdiag2), dimension(wfull), intent(in) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(in) :: dg3

if (incradbf.gt.0) then
  do iqw=1,wfull
    pg(iqw)%bf=dg2(iqw)%b0+grav*sum(dg3(iqw,1:pg(iqw)%mixind)%alpha*dg3(iqw,1:pg(iqw)%mixind)%rad)
  end do
else
  pg%bf=dg2%b0
end if

return
end subroutine getbf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the stability functions
!

subroutine getwx(wm,ws,dep,bf,ustar,mixdp)

implicit none

real, dimension(wfull), intent(out) :: wm,ws
real, dimension(wfull), intent(in) :: bf,ustar,mixdp,dep
real, dimension(wfull) :: zeta,sig,invl,uuu
real, parameter :: zetam=-0.2
real, parameter :: zetas=-1.0
real, parameter :: am=1.26
real, parameter :: cm=8.38
real, parameter :: as=-28.86
real, parameter :: cs=98.96

sig=dep/mixdp                       ! stable
where (bf.le.0..and.sig.gt.epsilon) ! unstable
  sig=epsilon
end where
uuu=ustar**3
invl=vkar*bf ! invl = ustar**3/L or L=ustar**3/(vkar*bf)
invl=max(invl,-10.*uuu/mixdp) ! MJT suggestion
invl=min(invl,1.*uuu/mixdp)   ! MJT suggestion
zeta=sig*mixdp*invl

where (zeta.gt.0.)
  wm=vkar*ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetam*uuu)
  wm=vkar*ustar*(1.-16.*zeta/uuu)**(1./4.)
elsewhere
  wm=vkar*(am*uuu-cm*zeta)**(1./3.)
end where

where (zeta.gt.0.)
  ws=vkar*ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetas*uuu)
  ws=vkar*ustar*(1.-16.*zeta/uuu)**(1./2.)
elsewhere
  ws=vkar*(as*uuu-cs*zeta)**(1./3.)
end where

return
end subroutine getwx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rho from equation of state
! From GFDL (MOM3)

subroutine getrho(dg2,dg3,atm)

implicit none

integer ii,i
integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(wfull) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32
real, dimension(wfull) :: drho0dt,drho0ds,dskdt,dskds,sk,sks
real, dimension(wfull) :: drhodt,drhods,rs0,rho0,hlra,hlrb
real, dimension(wfull,wlev) :: rs
real, parameter :: density = 1035.
type(tdiag2), dimension(wfull), intent(inout) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(inout) :: dg3
type(tatm), dimension(wfull), intent(in) :: atm

do ii=1,wlev
  dg3(:,ii)%rho=density
end do

do i=1,nits
  t = max(water(:,1)%temp-273.16,-2.)
  s = max(water(:,1)%sal,0.)
  t2 = t**2
  t3 = t**3
  t4 = t**4
  t5 = t**5
  s2 = s**2
  s3 = s**3
  s32 = sqrt(s3)

  rs0 = 999.842594 + 6.793952e-2*t(:) &
         - 9.095290e-3*t2(:) + 1.001685e-4*t3(:) &
         - 1.120083e-6*t4(:) + 6.536332e-9*t5(:) ! density for sal=0.
  rho0 = rs0+ s(:)*(0.824493 - 4.0899e-3*t(:) &
         + 7.6438e-5*t2(:) &
         - 8.2467e-7*t3(:) + 5.3875e-9*t4(:)) &
         + s32(:)*(-5.72466e-3 + 1.0227e-4*t(:) &
         - 1.6546e-6*t2(:)) + 4.8314e-4*s2(:)     ! + sal terms    
  drho0dt=6.793952e-2 &
         - 2.*9.095290e-3*t(:) + 3.*1.001685e-4*t2(:) &
         - 4.*1.120083e-6*t3(:) + 5.*6.536332e-9*t4(:) &
         + s(:)*( - 4.0899e-3 + 2.*7.6438e-5*t(:) &
         - 3.*8.2467e-7*t2(:) + 4.*5.3875e-9*t3(:)) &
         + s32(:)*(1.0227e-4 - 2.*1.6546e-6*t(:))
  drho0ds= (0.824493 - 4.0899e-3*t(:) &
         + 7.6438e-5*t2(:) &
         - 8.2467e-7*t3(:) + 5.3875e-9*t4(:)) &
         + 1.5*sqrt(s(:))*(-5.72466e-3 + 1.0227e-4*t(:) &
         - 1.6546e-6*t2(:)) + 2.*4.8314e-4*s(:)

  do ii=1,wlev
    t = max(water(:,ii)%temp-273.16,-2.)
    s = max(water(:,ii)%sal,0.)
    p1 = grav*depth(:,ii)*dg3(:,ii)%rho/100000.
    t2 = t**2
    t3 = t**3
    t4 = t**4
    t5 = t**5
    s2 = s**2
    s3 = s**3
    p2 = p1**2
    s32 = sqrt(s3)
    
    sks = 1.965933e4 + 1.444304e2*t(:) - 1.706103*t2(:) &
                + 9.648704e-3*t3(:)  - 4.190253e-5*t4(:) &
                + p1(:)*(3.186519 + 2.212276e-2*t(:) &
                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:))  &
                + p2(:)*(2.102898e-4 - 1.202016e-5*t(:) &
                + 1.394680e-7*t2(:)) ! sal=0.
    sk=sks+ s(:)*(52.84855 - 3.101089e-1*t(:) &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:)) &
                + s32(:)*(3.886640e-1 + 9.085835e-3*t(:) &
                - 4.619924e-4*t2(:)) &
                + p1(:)*s(:)*(6.704388e-3  -1.847318e-4*t(:) &
                + 2.059331e-7*t2(:)) + 1.480266e-4*p1(:)*s32(:) &
                +p2(:)*s(:)*(-2.040237e-6 &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:)) ! + sal terms             
    dskdt= 1.444304e2 - 2.*1.706103*t(:) &
                + 3.*9.648704e-3*t2(:)  - 4.*4.190253e-5*t3(:) &
                + s(:)*( - 3.101089e-1 &
                + 2.*6.283263e-3*t(:) -3.*5.084188e-5*t2(:)) &
                + s32(:)*(9.085835e-3 &
                - 2.*4.619924e-4*t(:)) &
                + p1(:)*(2.212276e-2 &
                - 2.*2.984642e-4*t(:) + 3.*1.956415e-6*t2(:))  &
                + p1(:)*s(:)*(-1.847318e-4 &
                + 2.*2.059331e-7*t(:)) &
                + p2(:)*(- 1.202016e-5 &
                + 2.*1.394680e-7*t(:)) +p2(:)*s(:)*( &
                + 6.128773e-8 + 2.*6.207323e-10*t(:))
    dskds=(52.84855 - 3.101089e-1*t(:) &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:)) &
                + 1.5*sqrt(s(:))*(3.886640e-1 + 9.085835e-3*t(:) &
                - 4.619924e-4*t2(:)) &
                + p1(:)*(6.704388e-3  -1.847318e-4*t(:) &
                + 2.059331e-7*t2(:)) + 1.5*1.480266e-4*p1(:)*sqrt(s(:)) &
                +p2(:)*(-2.040237e-6 &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
       
    dg3(:,ii)%rho=rho0/(1.-p1/sk)
    rs(:,ii)=rs0/(1.-p1/sks) ! sal=0.
  
    drhodt=drho0dt/(1.-p1/sk)-rho0*p1*dskdt/((sk-p1)**2)
    drhods=drho0ds/(1.-p1/sk)-rho0*p1*dskds/((sk-p1)**2)
    dg3(:,ii)%alpha=-drhodt/dg3(:,ii)%rho ! MOM convention
    dg3(:,ii)%beta=drhods/dg3(:,ii)%rho   ! MOM convention
  end do

end do

! buoyancy frequency
do ii=2,wlev-1
  dg3(:,ii)%nsq=-(grav/dg3(:,ii)%rho)*0.5*(dg3(:,ii-1)%rho-dg3(:,ii+1)%rho)/dz(:,ii)
end do
dg3(:,1)%nsq=2.*dg3(:,2)%nsq-dg3(:,3)%nsq
dg3(:,wlev)%nsq=2.*dg3(:,wlev-1)%nsq-dg3(:,wlev-2)%nsq

! shortwave
! use -ve as depth is down
hlrb=0.5*(dg3(:,1)%rho+dg3(:,2)%rho)
dg3(:,1)%rad=-atm%sg/cp0*((r_1*exp(-depth_hl(:,2)/mu_1)+r_2*exp(-depth_hl(:,2)/mu_2))/hlrb &
                         -(r_1*exp(-depth_hl(:,1)/mu_1)+r_2*exp(-depth_hl(:,1)/mu_2))/rho0)
do ii=2,wlev-1 
  hlra=0.5*(dg3(:,ii-1)%rho+dg3(:,ii)%rho)
  hlrb=0.5*(dg3(:,ii)%rho+dg3(:,ii+1)%rho)
  dg3(:,ii)%rad=-atm%sg/cp0*((r_1*exp(-depth_hl(:,ii+1)/mu_1)+r_2*exp(-depth_hl(:,ii+1)/mu_2))/hlrb &
                            -(r_1*exp(-depth_hl(:,ii)/mu_1)+r_2*exp(-depth_hl(:,ii)/mu_2))/hlra)
end do
hlra=0.5*(dg3(:,wlev-1)%rho+dg3(:,wlev)%rho)
dg3(:,wlev)%rad=atm%sg/cp0*(r_1*exp(-depth_hl(:,wlev)/mu_1)+r_2*exp(-depth_hl(:,wlev)/mu_2))/hlra ! remainder

! Boundary conditions (use rho at level 1 for consistancy with shortwave radiation)
dg2%wu0=-atm%taux/dg3(:,1)%rho                             ! BC
dg2%wv0=-atm%tauy/dg3(:,1)%rho                             ! BC
dg2%wt0=-(-atm%fg-atm%eg+atm%rg)/(dg3(:,1)%rho*cp0)        ! BC
dg2%ws0=(atm%rnd-atm%eg/lv)*water(:,1)%sal/rs(:,1)         ! BC ! negect ice for now

dg2%ustar=max(sqrt(sqrt(dg2%wu0*dg2%wu0+dg2%wv0*dg2%wv0)),1.E-10)
dg2%b0=-grav*(dg3(:,1)%alpha*dg2%wt0-dg3(:,1)%beta*dg2%ws0)
!dg2%b0=-grav*(-drho0dt*dg2%wt0-drho0ds*dg2%ws0)

return
end subroutine getrho

end module mlo
