
! This is a 1D, mixed layer ocean model for regioal climate simulations based on Large, et al (1994)
! (i.e., adapted from the GFDL code)

! Specifically, this version can assimilate SSTs from GCMs, using a convolution based digital filter
! which avoids problems with complex land-sea boundary conditions

module mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,wlev

type tdata
  real temp,sal,u,v
end type tdata
type tprog2
  real mixdepth
  integer mixind
end type tprog2
type tatm
  real sg,rg,fg,eg,rnd,taux,tauy,f
end type tatm
type tout
  real sst
end type tout
type tdiag2
  real bf,ustar,wu0,wv0,wt0,ws0,rho0
end type tdiag2
type tdiag3
  real rho,nsq
end type tdiag3

! parameters
integer, parameter :: wlev = 11
! model arrays
integer, save :: wfull
integer, dimension(:), allocatable, save :: wgrid
real, dimension(:,:), allocatable, save :: depth,dz
real, dimension(:,:), allocatable, save :: depth_hl
real, dimension(:,:), allocatable, save :: dz_hl
type(tdata), dimension(:,:), allocatable, save :: water
type(tprog2), dimension(:), allocatable, save :: pg

! max depth
real, parameter :: mxd =221.

! model parameters
real, parameter :: ric     = 0.3
real, parameter :: epsilon = 0.1

! radiation parameters
real, parameter :: r_lw = 0.58
real, parameter :: r_sw = 0.42
real, parameter :: mu_lw = 0.35
real, parameter :: mu_sw = 23.

! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: grav=9.80              ! graviational constant (m/s^2)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialise MLO


subroutine mloinit(ifull,sigw,depin,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iq,iqw,ii
real, dimension(ifull), intent(in) :: sigw,depin

if (diag.ge.1) write(6,*) "Initialising MLO"

wfull=count(sigw.gt.0)
if (wfull.eq.0) return

allocate(water(wfull,wlev),wgrid(wfull),pg(wfull))
allocate(depth(wfull,wlev),dz(wfull,wlev))
allocate(depth_hl(wfull,wlev+1))
allocate(dz_hl(wfull,2:wlev))

wgrid=0
iqw=0
do iq=1,ifull
  if (sigw(iq).gt.0.5) then
    iqw=iqw+1
    wgrid(iqw)=iq
  end if
end do

water%temp=288. ! K
water%sal=35.   ! PSU
water%u=0.      ! m/s
water%v=0.      ! m/s

pg%mixdepth=100. ! m
!pg%inflow=0.

!depth = (/ 0.5, 4.5, 13., 25., 41., 61., 85., 113., 145., 181., 221., 265., 313.,
!           365., 421., 481., 545., 613., 685., 761., 841, 925., 1013.  /)  ! 2.*(x-0.5)^2

do iqw=1,wfull
  call vgrid(depin(wgrid(iqw)),depth(iqw,:),depth_hl(iqw,:))
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
real dd,x,y,al,bt

dd=min(mxd,max(6.,depin))
x=real(wlev)-0.5
al=(x-dd)/(0.5*x-x*x)
bt=(x*x-0.5*dd)/(x*x-0.5*x)
do ii=1,wlev
  x=real(ii)-0.5
  depthout(ii)=al*x*x+bt*x
end do
do ii=1,wlev+1
  y=real(ii)-1.
  depth_hlout(ii)=al*y*y+bt*y ! ii is for half level ii-0.5
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
integer i
real, dimension(ifull,wlev,4), intent(in) :: datain

if (wfull.eq.0) return

do i=1,wlev
  water(:,i)%temp=datain(wgrid,i,1)
  water(:,i)%sal=datain(wgrid,i,2)
  water(:,i)%u=datain(wgrid,i,3)
  water(:,i)%v=datain(wgrid,i,4)
end do

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(ifull,dataout,depout,diag)

implicit none

integer, intent(in) :: ifull,diag
integer i
real, dimension(ifull,wlev,4), intent(out) :: dataout
real, dimension(ifull), intent(out) :: depout

if (wfull.eq.0) return

do i=1,wlev
  dataout(wgrid,i,1)=water(:,i)%temp
  dataout(wgrid,i,2)=water(:,i)%sal
  dataout(wgrid,i,3)=water(:,i)%u
  dataout(wgrid,i,4)=water(:,i)%v
end do
depout=0.
depout(wgrid)=depth(:,wlev)

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mloimport(ifull,sst,ilev,diag)

implicit none

integer, intent(in) :: ifull,ilev,diag
real, dimension(ifull), intent(in) :: sst

water(:,ilev)%temp=sst(wgrid)

return
end subroutine mloimport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mloexport(ifull,sst,ilev,diag)

implicit none

integer, intent(in) :: ifull,ilev,diag
real, dimension(ifull), intent(inout) :: sst

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
if (all(depin.eq.depth(:,wlev))) return

do iqw=1,wfull
  call vgrid(depin(wgrid(iqw)),dpin,dp_hlin)
  ii=1
  do while (ii.le.wlev)
    if (depth(iqw,ii).gt.dpin(wlev)) then
      newdat(iqw,ii,:)=mlodat(iqw,wlev,:)
    else if (depth(iqw,ii).lt.dpin(1)) then
      newdat(iqw,ii,:)=mlodat(iqw,1,:)
    else
      pos=maxloc(dpin,dpin.le.depth(iqw,ii))
      pos(1)=max(1,min(wlev-1,pos(1)))
      x=(depth(iqw,ii)-dpin(pos(1)))/(dpin(pos(1)+1)-dpin(pos(1)))
      x=max(0.,min(1.,x))
      newdat(iqw,ii,:)=mlodat(iqw,pos(1)+1,:)*x+mlodat(iqw,pos(1),:)*(1.-x)
    end if
    ii=ii+1
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
integer ii
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
real, dimension(wfull,2:wlev) :: aa
real, dimension(wfull,wlev) :: bb,dd
real, dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull) :: xp,xm
type(tout), dimension(wfull), intent(out) :: uo
type(tatm), dimension(wfull), intent(in) :: atm
type(tdiag2), dimension(wfull) :: dg2
type(tdiag3), dimension(wfull,wlev) :: dg3
type(tdata), dimension(wfull,wlev) :: new

call getrho(dg2,dg3,atm)           ! calculate rho0,rho and bf.  Also calculate boundary conditions.
call getmixdepth(dg2,dg3)          ! solve for mixed layer depth
call getstab(km,ks,gammas,dg2,dg3) ! solve for stability functions and non-local term

! TEMPERATURE
! radiative contribution to non-local term is negected for now (i.e., d2%wt0 only)
! use +ve for gamma terms as depth is down
rhs(:,1)=(0.25/dz(:,1))*(ks(:,2)+ks(:,1))*(gammas(:,2)+gammas(:,1))*dg2%wt0
do ii=2,wlev-1
  rhs(:,ii)=(0.25/dz(:,ii))*((ks(:,ii+1)+ks(:,ii))*(gammas(:,ii+1)+gammas(:,ii)) &
                          -(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dg2%wt0 ! non-local
end do
rhs(:,wlev)=0.
do ii=1,wlev-1
  ! use -ve as depth is down
  rhs(:,ii)=rhs(:,ii)-(atm%rg*exp(-depth_hl(:,ii+1)/mu_lw)+atm%sg*exp(-depth_hl(:,ii+1)/mu_sw) &
                      -atm%rg*exp(-depth_hl(:,ii)/mu_lw)-atm%sg*exp(-depth_hl(:,ii)/mu_sw))/(dg3(:,ii)%rho*cp0*dz(:,ii)) ! radiation
end do
rhs(:,wlev)=rhs(:,wlev)+(atm%rg*exp(-depth_hl(:,wlev)/mu_lw)+atm%sg*exp(-depth_hl(:,wlev)/mu_sw))/(dg3(:,wlev)%rho*cp0*dz(:,wlev)) ! remainder
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
  rhs(:,ii)=(0.25/dz(:,ii))*((ks(:,ii+1)+ks(:,ii))*(gammas(:,ii+1)+gammas(:,ii)) &
                          -(ks(:,ii)+ks(:,ii-1))*(gammas(:,ii)+gammas(:,ii-1)))*dg2%ws0 ! non-local
end do
rhs(:,wlev)=0.
dd(:,1)=water(:,1)%sal/dt+rhs(:,1)-dg2%ws0/dz(:,1)
do ii=2,wlev-1
  dd(:,ii)=water(:,ii)%sal/dt+rhs(:,ii)
end do
dd(:,wlev)=water(:,wlev)%sal/dt+rhs(:,wlev)
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
do ii=2,wlev-1
  dd(:,ii)=water(:,ii)%v/dt
end do
dd(:,wlev)=water(:,wlev)%v/dt
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
! tri-diaganol solver

subroutine thomas(out,aa,bbi,cc,ddi)

implicit none

integer ii
real, dimension(wfull,2:wlev), intent(in) :: aa
real, dimension(wfull,wlev), intent(in) :: bbi,ddi
real, dimension(wfull,1:wlev-1), intent(in) :: cc
real, dimension(wfull,wlev), intent(out) :: out
real, dimension(wfull,wlev) :: bb,dd
real, dimension(wfull) :: n

bb=bbi
dd=ddi

! solve for temperature
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for stability functions
! GFDL use a look-up table to speed-up the code...

subroutine getstab(km,ks,gammas,dg2,dg3)

implicit none

integer ii,iqw,iip1,iim1
real, dimension(wfull,wlev), intent(out) :: km,ks,gammas
real, dimension(wfull,wlev) :: num,nus,wm,ws
real, dimension(wfull) :: sigma,ri
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

! diffusion ---------------------------------------------------------
do ii=1,wlev

  iim1=max(ii-1,1)
  iip1=min(ii+1,wlev)
  ri=dg3(:,ii)%nsq*(depth(:,iim1)-depth(:,iip1))**2 &
           /max((water(:,iim1)%u-water(:,iip1)%u)**2 &
               +(water(:,iim1)%v-water(:,iip1)%v)**2,0.01)
  
  ! s - shear
  where (ri.lt.0.)
    num(:,ii)=nu0
  elsewhere (ri.lt.ri0)
    num(:,ii)=nu0*(1.-(ri/ri0)**2)**3
  elsewhere
    num(:,ii)=0.
  endwhere
  nus(:,ii)=num(:,ii)
  
  ! w - internal wave
  num(:,ii)=num(:,ii)+numw
  nus(:,ii)=nus(:,ii)+nusw
  
  ! d - double-diffusive
  ! double-diffusive mixing is neglected for now
  
end do
!--------------------------------------------------------------------

! stability ---------------------------------------------------------
do ii=1,wlev
  call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dg2%bf,dg2%ustar,pg%mixdepth)
end do
!--------------------------------------------------------------------

! calculate G profile -----------------------------------------------
do iqw=1,wfull
  xp=(pg(iqw)%mixdepth-depth(iqw,pg(iqw)%mixind))/(depth(iqw,pg(iqw)%mixind+1)-depth(iqw,pg(iqw)%mixind))
  numh=xp*num(iqw,pg(iqw)%mixind)+(1.-xp)*num(iqw,pg(iqw)%mixind+1)
  wm1=xp*wm(iqw,pg(iqw)%mixind)+(1.-xp)*wm(iqw,pg(iqw)%mixind+1)
  dnumhdz=(num(iqw,pg(iqw)%mixind+1)-num(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dwm1ds=pg(iqw)%mixdepth*(wm(iqw,pg(iqw)%mixind+1)-wm(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  nush=xp*nus(iqw,pg(iqw)%mixind)+(1.-xp)*nus(iqw,pg(iqw)%mixind+1)
  ws1=xp*ws(iqw,pg(iqw)%mixind)+(1.-xp)*ws(iqw,pg(iqw)%mixind+1)
  dnushdz=(nus(iqw,pg(iqw)%mixind+1)-nus(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dws1ds=pg(iqw)%mixdepth*(ws(iqw,pg(iqw)%mixind+1)-ws(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  
  g1m=numh/max(pg(iqw)%mixdepth*wm1,0.000001)
  dg1mds=-dnumhdz/max(wm1,0.000001)-numh*dwm1ds/max(pg(iqw)%mixdepth*wm1*wm1,0.000001)
  g1s=nush/max(pg(iqw)%mixdepth*ws1,0.000001)
  dg1sds=-dnushdz/max(ws1,0.000001)-nush*dws1ds/max(pg(iqw)%mixdepth*ws1*ws1,0.000001)


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
cg      = 10. * vkar * (98.96 * vkar * epsilon)**(1./3.)
do ii=1,wlev
  where (dg2%bf.lt.0.)
    gammas(:,ii)=cg/max(ws(:,ii)*pg%mixdepth,0.000001)
  elsewhere
    gammas(:,ii)=0.
  end where
end do

! enhancement is neglected for now

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixing layer depth

subroutine getmixdepth(dg2,dg3)

implicit none

integer ii,iqw,ie,pos(1)
real vtc,dvsq,vtsq,rsfc,xp,usfc,vsfc
real, dimension(wfull,wlev) :: ws,wm
real, dimension(wfull) :: dum
real, dimension(wlev) :: rib
type(tdiag2), dimension(wfull), intent(in) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(inout) :: dg3

vtc=1.8 * sqrt(0.2/(98.96*epsilon)) / (vkar**2*ric)

do ii=1,wlev
  dum=depth(:,ii) ! use current depth for boundary layer depth
  call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dg2%bf,dg2%ustar,dum)
end do

do iqw=1,wfull
  pos=maxloc(depth(iqw,:),depth(iqw,:).le.5.)
  ie=pos(1)

  pg(iqw)%mixind=wlev-1
  pg(iqw)%mixdepth=depth(iqw,wlev)
  ! should be averaged over 0 < sigma < epsilon instead of 1st level
  rsfc=sum(dg3(iqw,1:ie)%rho*dz(iqw,1:ie))/sum(dz(iqw,1:ie))
  usfc=sum(water(iqw,1:ie)%u*dz(iqw,1:ie))/sum(dz(iqw,1:ie))
  vsfc=sum(water(iqw,1:ie)%v*dz(iqw,1:ie))/sum(dz(iqw,1:ie))
  do ii=ie,wlev
    vtsq=depth(iqw,ii)*ws(iqw,ii)*sqrt(abs(dg3(iqw,ii)%nsq))*vtc
    dvsq=(usfc-water(iqw,ii)%u)**2+(vsfc-water(iqw,ii)%v)**2
    rib(ii)=(depth(iqw,ii)-depth(iqw,ie))*(1.-rsfc/dg3(iqw,ii)%rho)/max(dvsq+vtsq,0.000001)
  
    if (rib(ii).gt.ric.and.ii.gt.ie) then
      pg(iqw)%mixind=ii-1
      xp=min(max((rib(ii-1)-ric)/(rib(ii-1)-rib(ii)),0.),1.)
      pg(iqw)%mixdepth=(1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii)
      exit
    end if
  
  end do 

end do

return
end subroutine getmixdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the stability functions
!

subroutine getwx(wm,ws,dep,bf,ustar,dp)

implicit none

real, dimension(wfull), intent(out) :: wm,ws
real, dimension(wfull), intent(in) :: bf,ustar,dp,dep
real, dimension(wfull) :: zeta,sig,l
real, parameter :: zetam=-0.2
real, parameter :: zetas=-1.0
real, parameter :: am=1.26
real, parameter :: cm=8.38
real, parameter :: as=-28.86
real, parameter :: cs=98.96

where (bf.lt.0.)
  sig=epsilon
elsewhere
  sig=dep/dp
end where
l=vkar*bf/max(ustar**3,0.000001)
zeta=sig*dp/l

where (zeta.gt.0.)
  wm=vkar*ustar/max(1.+5.*zeta,0.000001)
elsewhere (zeta.gt.zetam)
  wm=vkar*ustar/max((1.-16.*zeta)**(-1./4.),0.000001)
elsewhere
  wm=vkar*ustar/max((am-cm*zeta)**(-1./3.),0.000001)
end where

where (zeta.gt.0.)
  ws=vkar*ustar/max(1.+5.*zeta,0.000001)
elsewhere (zeta.gt.zetas)
  ws=vkar*ustar/max((1.-16.*zeta)**(-1./2.),0.000001)
elsewhere
  ws=vkar*ustar/max((as-cs*zeta)**(-1./3.),0.000001)
end where

return
end subroutine getwx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rho from equation of state
! From GFDL (MOM3)

subroutine getrho(dg2,dg3,atm)

implicit none

integer ii,iim1,iip1,i
integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(wfull) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32
real, dimension(wfull) :: drho0dt,drho0ds,dskdt,dskds,sk
real, dimension(wfull) :: drhodt,drhods
real, dimension(wfull,wlev) :: alpha,beta
real, parameter :: density = 1035.
type(tdiag2), dimension(wfull), intent(inout) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(inout) :: dg3
type(tatm), dimension(wfull), intent(in) :: atm

do ii=1,wlev
  dg3(:,ii)%rho=density
end do

do i=1,nits
  t = max(water(:,1)%temp-273.16,-2.)
  s = water(:,1)%sal
  t2 = t**2
  t3 = t**3
  t4 = t**4
  t5 = t**5
  s2 = s**2
  s3 = s**3
  s32 = sqrt(s3)

  dg2%rho0 = 999.842594 + 6.793952e-2*t(:) &
         - 9.095290e-3*t2(:) + 1.001685e-4*t3(:) &
         - 1.120083e-6*t4(:) + 6.536332e-9*t5(:) &
         + s(:)*(0.824493 - 4.0899e-3*t(:) &
         + 7.6438e-5*t2(:) &
         - 8.2467e-7*t3(:) + 5.3875e-9*t4(:)) &
         + s32(:)*(-5.72466e-3 + 1.0227e-4*t(:) &
         - 1.6546e-6*t2(:)) + 4.8314e-4*s2(:)
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
    s = water(:,ii)%sal
    p1 = grav*depth(:,ii)*dg3(:,ii)%rho/100000.
    t2 = t**2
    t3 = t**3
    t4 = t**4
    t5 = t**5
    s2 = s**2
    s3 = s**3
    p2 = p1**2
    s32 = sqrt(s3)

    sk = 1.965933e4 + 1.444304e2*t(:) - 1.706103*t2(:) &
                + 9.648704e-3*t3(:)  - 4.190253e-5*t4(:) &
                + s(:)*(52.84855 - 3.101089e-1*t(:) &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:)) &
                + s32(:)*(3.886640e-1 + 9.085835e-3*t(:) &
                - 4.619924e-4*t2(:)) &
                + p1(:)*(3.186519 + 2.212276e-2*t(:) &
                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:))  &
                + p1(:)*s(:)*(6.704388e-3  -1.847318e-4*t(:) &
                + 2.059331e-7*t2(:)) + 1.480266e-4*p1(:)*s32(:) &
                + p2(:)*(2.102898e-4 - 1.202016e-5*t(:) &
                + 1.394680e-7*t2(:)) +p2(:)*s(:)*(-2.040237e-6 &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
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
       
    dg3(:,ii)%rho=dg2%rho0/(1.-p1/sk)
  
    drhodt=drho0dt/(1.-p1/sk)+dg2%rho0*p1*dskdt/((sk-p1)**2)
    drhods=drho0ds/(1.-p1/sk)+dg2%rho0*p1*dskds/((sk-p1)**2)
    alpha(:,ii)=-drhodt
    beta(:,ii)=drhods
  end do

end do

do ii=1,wlev
  iim1=max(ii-1,1)
  iip1=min(ii+1,wlev)
  ! +ve as depth is down
  dg3(:,ii)%nsq=(grav/dg3(:,ii)%rho)*(dg3(:,iim1)%rho-dg3(:,iip1)%rho)/(depth(:,iim1)-depth(:,iip1))
end do

! Boundary conditions
dg2%wu0=-atm%taux/dg2%rho0                             ! BC
dg2%wv0=-atm%tauy/dg2%rho0                             ! BC
dg2%wt0=-(-atm%fg-atm%eg)/(dg2%rho0*cp0)               ! BC
dg2%ws0=(atm%rnd-atm%eg/lv)*water(:,1)%sal/dg2%rho0    ! BC ! negect ice for now

dg2%ustar=sqrt(sqrt(dg2%wu0*dg2%wu0+dg2%wv0*dg2%wv0))
!dg2%bf=-grav*(alpha(:,1)*dg2%wt0-beta(:,1)*dg2%ws0)
dg2%bf=-grav*(-drho0dt*dg2%wt0-drho0ds*dg2%ws0)

return
end subroutine getrho

end module mlo